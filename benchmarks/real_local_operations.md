# Real dense gates and Kraus operations

## Results

PR #36 already sped up real density-matrix gates. Comparing the pre-PR build to the baseline, two-target gates improved 1.19× on low targets and 1.39× on spread targets. Three-target gates improved 1.20× and 1.78×. Both density-matrix gate passes use the statevector gate kernels.

Times below are microseconds per call with one thread. Speedup is baseline time divided by candidate time. For rows marked with an asterisk, the historical API lacked a three-target Kraus entry point and the baseline uses a composition described below.

| Dense real operation | Low targets, baseline → candidate | Speedup | Spread targets, baseline → candidate | Speedup |
|---|---:|---:|---:|---:|
| Statevector Kraus trajectory, 2 targets | 1,042 → 457 | 2.28× | 531 → 289 | 1.84× |
| Statevector Kraus trajectory, 3 targets* | 1,272 → 743 | 1.71× | 625 → 493 | 1.27× |
| Statevector local map with output norm, 2 targets | 375 → 255 | 1.47× | 205 → 140 | 1.46× |
| Statevector local map with output norm, 3 targets | 782 → 405 | 1.93× | 401 → 238 | 1.69× |
| Density gate, 2 targets | 375 → 305 | 1.23× | 305 → 300 | 1.02× |
| Density gate, 3 targets | 707 → 459 | 1.54× | 459 → 465 | 0.99× |
| Density Kraus channel, 1 target | 634 → 173 | 3.67× | 305 → 137 | 2.23× |
| Density Kraus channel, 2 targets | 1,908 → 568 | 3.36× | 983 → 573 | 1.72× |
| Density Kraus channel, 3 targets* | 1,784 → 1,186 | 1.50× | 1,484 → 1,216 | 1.22× |

At four threads the real two-target Kraus speedups were 1.91–2.15× for statevectors and 1.63–3.19× for density matrices. Three-target statevector speedups were 1.36–1.70×. The three-target density channel took 1,426 µs versus 4,553 µs on low targets and 1,467 µs versus 4,277 µs on spread targets, giving 3.19× and 2.91× speedups. It accumulates branches in tiles of up to 4,096 complex entries using two local buffers per worker.

For complex coefficients, two-target statevector trajectories sped up 1.36–1.81× and two-target density channels sped up 1.51–2.84×. Three-target complex density channels ran 1.14–1.21× faster than the reference. Three-target complex statevector trajectories took 1–6% longer than their reference. Sparse controls were within 8% of baseline. Across the 52 cases, the largest difference in final-state checksums between builds was 1.8 × 10⁻¹⁵.

The corrected SIMD kernel adds input columns in matrix order for both output rows. This preserves small residuals and avoids premature overflow in cancellation cases. Compared to the optimized implementation in `734263db`, the low-target real two-qubit density channel took 568 µs versus 456 µs serially, and 688 µs versus 585 µs on four threads. The original serial baseline took 3.36× as long as the corrected implementation. The other serial real Kraus trajectories and density channels were within 6% of `734263db`.

An earlier isolated check on the two-target complex statevector map-with-norm operation gave 385 µs for the baseline, 429 µs for an intermediate implementation, and 373 µs for `734263db` after separating the real loop. The `focused_norm` rows preserve those historical results. The current mixed-order benchmark gave 368 µs for the baseline and 366 µs for the corrected implementation.

## Measurement

The benchmark uses the public simulator APIs on an Intel Core i9-13900KS running Windows x64, with MSVC 19.52, AVX2, OpenMP and Eigen 5.0.0. All four builds use the same harness and compiler flags. The pre-PR build uses `4e5705b4`. The baseline (`head` in the CSV) uses `049bff4a`, which includes PR #36. The other two builds use `734263db` and the corrected implementation at `3d7b0703` (`candidate` in the CSV).

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

Omit `/DNATIVE_THREE_KRAUS` on `049bff4a` and `4e5705b4` to enable the three-target reference compositions. Use `/DNATIVE_THREE_KRAUS` with `734263db` and the corrected implementation. Run `density-gate 7 8 9 18 1` for the pre-PR comparison. The `PIN_BENCH` flag selects the processor mask for this host and can be omitted on a different system.
