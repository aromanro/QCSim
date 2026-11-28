# Real dense gates and Kraus operations

## Three-target complex statevector trajectories

The new kernel forms the eight-by-eight reduced density matrix from tiles of 64 amplitude groups. In an instrumented 18-qubit build of `3d7b0703`, the reduction took about 818–857 µs on low targets and 400–417 µs on spread targets. Applying the selected dense branch took about 387–416 µs. Computing the two weights and drawing the branch took about 1 µs. The largest opportunity was the reduction.

The focused timings below use the same dense complex channel and target layouts for all four implementations. Values are medians in microseconds per complete trajectory step. The final column divides the original composed time by the new native time.

| Register and targets | Original composition | Previous native | New native | New composition | Speedup |
|---|---:|---:|---:|---:|---:|
| 18 qubits, one thread, low | 1,297 | 1,227 | 693 | 702 | 1.87× |
| 18 qubits, one thread, spread | 804 | 824 | 722 | 718 | 1.11× |
| 20 qubits, four threads, low | 1,645 | 1,539 | 888 | 904 | 1.85× |
| 20 qubits, four threads, spread | 1,029 | 1,036 | 934 | 913 | 1.10× |

The new native trajectory beat the original composition in all eight matched rounds on each layout. Against the previous native implementation it was 1.77× and 1.14× faster at 18 qubits, and 1.73× and 1.11× faster at 20 qubits. The new composition uses the same optimized public reduction. Its medians and the new native timings differed by at most 2.4%.

The earlier full sweep measured the previous native trajectory 1–6% slower than the original composition at 18 qubits. In the longer paired runs, only the spread layout reproduced a slowdown (824 versus 804 µs). On low targets the previous native trajectory was faster (1,227 versus 1,297 µs).

At eight qubits there are 32 groups, below the 64-group dispatch threshold. Old and new native calls differed by about 0.01 µs at that size. At 9, 10 and 12 qubits, low-target medians fell from 3.63 to 2.82 µs, 5.85 to 4.02 µs and 19.63 to 11.63 µs. Spread-target medians fell from 2.83 to 2.80 µs, 4.28 to 4.12 µs and 13.32 to 12.35 µs.

## Other operations

PR #36 already sped up real density-matrix gates. Comparing the pre-PR build to the baseline, two-target gates improved 1.18× on low targets and 1.41× on spread targets. Three-target gates improved 1.24× and 1.71×. Both density-matrix gate passes use the statevector gate kernels.

Times below are microseconds per call with one thread. Speedup is baseline time divided by candidate time. For rows marked with an asterisk, the historical API lacked a three-target Kraus entry point and the baseline uses a composition described below.

| Dense real operation | Low targets, baseline → candidate | Speedup | Spread targets, baseline → candidate | Speedup |
|---|---:|---:|---:|---:|
| Statevector Kraus trajectory, 2 targets | 1,012 → 454 | 2.23× | 536 → 289 | 1.85× |
| Statevector Kraus trajectory, 3 targets* | 1,243 → 740 | 1.68× | 611 → 493 | 1.24× |
| Statevector local map with output norm, 2 targets | 370 → 247 | 1.50× | 202 → 137 | 1.47× |
| Statevector local map with output norm, 3 targets | 793 → 389 | 2.04× | 413 → 223 | 1.85× |
| Density gate, 2 targets | 363 → 303 | 1.20× | 305 → 302 | 1.01× |
| Density gate, 3 targets | 651 → 445 | 1.46× | 455 → 453 | 1.00× |
| Density Kraus channel, 1 target | 597 → 167 | 3.57× | 293 → 138 | 2.13× |
| Density Kraus channel, 2 targets | 1,884 → 569 | 3.31× | 980 → 538 | 1.82× |
| Density Kraus channel, 3 targets* | 1,841 → 1,196 | 1.54× | 1,465 → 1,233 | 1.19× |

At four threads the real two-target Kraus speedups were 1.92–2.10× for statevectors and 1.60–3.10× for density matrices. Three-target statevector speedups were 1.39–1.78×. The three-target density channel took 1,526 µs versus 5,074 µs on low targets and 1,510 µs versus 4,485 µs on spread targets, giving 3.33× and 2.97× speedups.

For complex coefficients, two-target statevector trajectories sped up 1.42–1.90× and two-target density channels sped up 1.57–2.67× over the baseline. Three-target complex density channels ran 1.17–1.20× faster than the reference. Sparse controls were within 7% of baseline. Across the full, focused and small-register comparisons, the largest difference in final-state checksums between implementations was 2.0 × 10⁻¹⁵.

The ordered dense-map kernel introduced in `3d7b0703` preserves small residuals and avoids premature overflow in cancellation cases. On low targets the real two-qubit density channel took 569 µs versus 481 µs for `734263db` serially, and 707 µs versus 597 µs on four threads. The original serial baseline took 3.31× as long as the current implementation.

An earlier isolated check on the two-target complex statevector map-with-norm operation gave 385 µs for the baseline, 429 µs for an intermediate implementation, and 373 µs for `734263db` after separating the real loop. The `focused_norm` rows preserve those historical results. The current mixed-order benchmark gave 367 µs for the baseline and 362 µs for `687b8e0`.

## Measurement

The benchmark uses the public simulator APIs on an Intel Core i9-13900KS running Windows x64, with MSVC 19.52, AVX2, OpenMP and Eigen 5.0.0. All five builds use the same harness and compiler flags. `4e5705b4` is the pre-PR build. `049bff4a` is the baseline (`head` in the CSV) and includes PR #36. `734263db` added the native three-target APIs. `3d7b0703` is the previous native build. `687b8e0` is the new build (`candidate` in the CSV). `candidate_composed` runs the reference operations in the same executable as `candidate`.

Each gate is a fixed dense real orthogonal matrix. The complex controls add different phases to its rows. Each dense channel has two nonunitary Kraus operators with different singular values and a complete sum of effects. The sparse control mixes identity and a simultaneous bit flip with probabilities 0.6 and 0.4. All cases start from seeded complex amplitudes. Statevector trajectories alternate fixed draws of 0.13 and 0.87. The density channel applies both operators on every call.

The low layout uses targets (0, 1) or (0, 1, 2). In an *n*-qubit register, the spread layout uses (`n - 2`, 2) or (`n - 2`, 2, `floor(n / 2) + 1`). The first listed target is the low bit of the local matrix index. The timings include branch-weight computation, channel validation and application where the API performs those steps. Initial-state and input-matrix construction are outside the timer.

Serial measurements use 9 density qubits or 18 statevector qubits, each with 4 MiB of complex entries. Larger measurements use four OpenMP threads and 10 density qubits or 20 statevector qubits, each with 16 MiB. The benchmark pins the serial process to logical processor 4 and the four-thread process to processors 4, 6, 8 and 10. Each case gets one warmup operation. The full serial result is the median of five process medians, each computed from seven batches of eight calls. The full parallel result uses three processes with five batches of eight calls each. The revision order varies between runs.

The focused `complex_trajectory` series uses eight rounds per configuration, 11 batches of 24 calls per process and a balanced order for the four implementations. The `complex_threshold` series alternates old and new native calls over four rounds for each register size. Each process measures seven batches of 512 calls at 8–10 qubits and 128 calls at 12 qubits.

The baseline already supports two-target Kraus trajectories and density channels. For three targets the baseline runs equivalent compositions of public operations. The trajectory reference computes the reduced density matrix, evaluates both branch weights, draws the branch and applies its normalized gate. The density reference copies the input to reusable storage, applies one gate per branch, adds the two matrices and installs the sum. Comparisons for three-target channels therefore measure the cost relative to these particular compositions.

The [sample CSV](real_local_operations.csv) records each process median in microseconds per call, its individual batch measurements and a final-state checksum. The [standalone harness](real_local_operations.cpp) defines all operators, target sets and seeds.

## Reproduction

Use a Visual Studio x64 developer command prompt. Point `QCSIM_HEADERS` to the chosen revision's `QCSim` directory and `EIGEN5_INCLUDE_DIR` to Eigen's include directory. The following builds and runs the candidate from the repository root.

```bat
cl /nologo /std:c++17 /O2 /EHsc /bigobj /arch:AVX2 /fp:precise /openmp /DNDEBUG /DPIN_BENCH /DNATIVE_THREE_KRAUS /I"%QCSIM_HEADERS%" /I"%EIGEN5_INCLUDE_DIR%" benchmarks\real_local_operations.cpp /Fe:real_local_operations.exe
real_local_operations.exe all 7 8 9 18 1
real_local_operations.exe real 5 8 10 20 4
real_local_operations.exe complex-trajectory 11 24 9 18 1
real_local_operations.exe composed-trajectory 11 24 9 18 1
real_local_operations.exe complex-trajectory 11 24 10 20 4
real_local_operations.exe composed-trajectory 11 24 10 20 4
```

Omit `/DNATIVE_THREE_KRAUS` on `049bff4a` and `4e5705b4` to enable the three-target reference compositions. Use `/DNATIVE_THREE_KRAUS` with `734263db`, `3d7b0703` and `687b8e0`. The `composed-trajectory` filter selects the composed reference even in a native build. Run `density-gate 7 8 9 18 1` for the pre-PR comparison. The `PIN_BENCH` flag selects the processor mask for this host and can be omitted on a different system.
