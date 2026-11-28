# Real dense gates and Kraus operations

## Results

PR #36 already sped up real density-matrix gates. Comparing the pre-PR build to the baseline, two-target gates improved 1.19× on low targets and 1.41× on spread targets. Three-target gates improved 1.23× and 1.78×. Both density-matrix gate passes use the statevector gate kernels.

Times below are microseconds per call with one thread. Speedup is baseline time divided by candidate time. For rows marked with an asterisk, the historical API lacked a three-target Kraus entry point and the baseline uses a composition described below.

| Dense real operation | Low targets, baseline → candidate | Speedup | Spread targets, baseline → candidate | Speedup |
|---|---:|---:|---:|---:|
| Statevector Kraus trajectory, 2 targets | 1,056 → 465 | 2.27× | 513 → 299 | 1.72× |
| Statevector Kraus trajectory, 3 targets* | 1,266 → 722 | 1.75× | 608 → 473 | 1.28× |
| Statevector local map with output norm, 2 targets | 361 → 253 | 1.42× | 201 → 140 | 1.44× |
| Statevector local map with output norm, 3 targets | 755 → 400 | 1.89× | 402 → 226 | 1.78× |
| Density gate, 2 targets | 359 → 308 | 1.16× | 302 → 297 | 1.01× |
| Density gate, 3 targets | 652 → 428 | 1.52× | 438 → 449 | 0.98× |
| Density Kraus channel, 1 target | 566 → 176 | 3.21× | 293 → 133 | 2.20× |
| Density Kraus channel, 2 targets | 1,881 → 468 | 4.02× | 959 → 534 | 1.80× |
| Density Kraus channel, 3 targets* | 1,718 → 1,203 | 1.43× | 1,334 → 1,185 | 1.13× |

At four threads the real two-target Kraus speedups were 1.68–2.01× for statevectors and 1.65–3.50× for density matrices. Three-target statevector speedups were 1.22–1.62×. The three-target density channel took 1,445 µs versus 4,727 µs on low targets and 1,555 µs versus 4,326 µs on spread targets, giving 3.27× and 2.78× speedups. It accumulates branches in tiles of up to 4,096 complex entries using two local buffers per worker.

For complex coefficients, two-target statevector trajectories sped up 1.35–1.82× and two-target density channels sped up 1.53–2.77×. Three-target complex density channels ran 1.07–1.13× faster than the reference. Three-target complex statevector trajectories took 9–10% longer than their reference. Sparse controls were within 5% of baseline. Across the 52 cases, the largest difference in final-state checksums between builds was 1.6 × 10⁻¹⁵.

A separate check isolated the two-target complex statevector map-with-norm operation on low targets. Five runs with seven batches of 32 calls gave 385 µs for the baseline, 429 µs for an intermediate implementation, and 373 µs after moving the real arithmetic into its own helper. The shorter batches in the full benchmark gave 369 µs for the baseline and 375 µs for the final implementation. The CSV labels the longer runs `focused_norm`.

## Measurement

The benchmark uses the public simulator APIs on an Intel Core i9-13900KS running Windows x64, with MSVC 19.52, AVX2, OpenMP and Eigen 5.0.0. All three builds use the same harness and compiler flags. The pre-PR build uses `4e5705b4`. The baseline uses `049bff4a`, which includes PR #36. The candidate uses the accompanying source changes.

Each gate is a fixed dense real orthogonal matrix. The complex controls add different phases to its rows. Each dense channel has two nonunitary Kraus operators with different singular values and a complete sum of effects. The sparse control mixes identity and a simultaneous bit flip with probabilities 0.6 and 0.4. All cases start from seeded complex amplitudes. Statevector trajectories alternate fixed draws of 0.13 and 0.87. The density channel applies both operators on every call.

The low layout uses targets (0, 1) or (0, 1, 2). In an *n*-qubit register, the spread layout uses (`n - 2`, 2) or (`n - 2`, 2, `floor(n / 2) + 1`). The first listed target is the low bit of the local matrix index. The timings include branch-weight computation, channel validation and application where the API performs those steps. Initial-state and input-matrix construction are outside the timer.

Serial measurements use 9 density qubits or 18 statevector qubits, each with 4 MiB of complex entries. Larger measurements use four OpenMP threads and 10 density qubits or 20 statevector qubits, each with 16 MiB. The benchmark pins the serial process to logical processor 4 and the four-thread process to processors 4, 6, 8 and 10. Each case gets one warmup operation. The serial result is the median of five process medians, each computed from seven batches of eight calls. The parallel result uses three processes with five batches of eight calls each. The revision order varies between runs.

The baseline already supports two-target Kraus trajectories and density channels. For three targets the baseline runs equivalent compositions of public operations. The trajectory reference computes the reduced density matrix, evaluates both branch weights, draws the branch and applies its normalized gate. The density reference copies the input to reusable storage, applies one gate per branch, adds the two matrices and installs the sum. Comparisons for three-target channels therefore measure the cost relative to these particular compositions.

The [sample CSV](real_local_operations.csv) records each process median in microseconds per call, its individual batch measurements and a final-state checksum. The [standalone harness](real_local_operations.cpp) defines all operators, target sets and seeds.

## Reproduction

Use a Visual Studio x64 developer command prompt. Point `QCSIM_HEADERS` to the chosen revision's `QCSim` directory and `EIGEN5_INCLUDE_DIR` to Eigen's include directory. The following builds and runs the candidate from the repository root.

```bat
cl /nologo /std:c++17 /O2 /EHsc /bigobj /arch:AVX2 /fp:precise /openmp /DNDEBUG /DPIN_BENCH /DNATIVE_THREE_KRAUS /I"%QCSIM_HEADERS%" /I"%EIGEN5_INCLUDE_DIR%" benchmarks\real_local_operations.cpp /Fe:real_local_operations.exe
real_local_operations.exe all 7 8 9 18 1
real_local_operations.exe real 5 8 10 20 4
```

Build the same file against each historical revision without `/DNATIVE_THREE_KRAUS` to enable the three-target reference compositions. Run `density-gate 7 8 9 18 1` for the pre-PR comparison. The `PIN_BENCH` flag selects the processor mask for this host and can be omitted on a different system.
