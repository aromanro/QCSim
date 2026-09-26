#pragma once


#ifdef _WIN32
#include <windows.h>
#undef min
#undef max
#endif // _WIN32

#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Eigen>

#include <random>
#include <complex>
#include <chrono>
#include <algorithm>

#include <iostream>
#include <iterator>
#include <iomanip>
#include <fstream>

#include <vector>
#include <array>
#include <atomic>
#include <cstdint>
#include <thread>
#include <type_traits>

#if !defined(QCSIM_DISABLE_SIMD) && (defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64))
#include <emmintrin.h>
#define QCSIM_SSE2 1
#endif

#if !defined(QCSIM_DISABLE_SIMD) && !defined(QCSIM_DISABLE_AVX2) && (defined(__AVX2__) || defined(_M_AVX2))
#include <immintrin.h>
#define QCSIM_AVX2 1
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "QuantumGate.h"


namespace QC {

	// Mixed-circuit measurements after SIMD favor serial execution below
	// 14 qubits on both tested platforms. This default remains tunable.
	inline constexpr size_t DefaultParallelMinBasisStates = 16384;

	// Shared by all calculator instantiations (statevector, density matrix rows and columns)
	inline std::atomic<size_t>& ParallelMinBasisStatesSetting()
	{
		static std::atomic<size_t> value{ DefaultParallelMinBasisStates };
		return value;
	}

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegisterCalculator {
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;

		QubitRegisterCalculator() = default;
		virtual ~QubitRegisterCalculator() = default;

		// No output register is required. The descriptor remains valid under
		// conjugation: all arithmetic coefficients come from the supplied matrix.
		static inline void ApplyGateInPlace(VectorClass& state, const MatrixClass& matrix,
			const Gates::GateStructure& structure, const std::array<size_t, 3>& bits,
			unsigned gateQubits, size_t count, bool parallel)
		{
			assert(gateQubits >= 1 && gateQubits <= 3 && structure.controls < gateQubits);
			if (parallel)
			{
				if (gateQubits == 1) DispatchGate<true, 1>(structure, state, matrix, bits, count);
				else if (gateQubits == 2) DispatchGate<true, 2>(structure, state, matrix, bits, count);
				else DispatchGate<true, 3>(structure, state, matrix, bits, count);
			}
			else
			{
				if (gateQubits == 1) DispatchGate<false, 1>(structure, state, matrix, bits, count);
				else if (gateQubits == 2) DispatchGate<false, 2>(structure, state, matrix, bits, count);
				else DispatchGate<false, 3>(structure, state, matrix, bits, count);
			}
		}

		// Compatibility entry points. The state is always updated in place.

		static inline void ApplyOneQubitGate(const GateClass& gate, VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t NrBasisStates)
		{
			DispatchGate<false, 1>(gate.getStructure(), registerStorage, gateMatrix, {qubitBit, 0, 0}, NrBasisStates);
		}

		static inline void ApplyOneQubitGateOmp(const GateClass& gate, VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t NrBasisStates)
		{
			DispatchGate<true, 1>(gate.getStructure(), registerStorage, gateMatrix, {qubitBit, 0, 0}, NrBasisStates);
		}

		static inline void ApplyTwoQubitsGate(const GateClass& gate, VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			DispatchGate<false, 2>(gate.getStructure(), registerStorage, gateMatrix, {qubitBit, ctrlQubitBit, 0}, NrBasisStates);
		}

		static inline void ApplyTwoQubitsGateOmp(const GateClass& gate, VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			DispatchGate<true, 2>(gate.getStructure(), registerStorage, gateMatrix, {qubitBit, ctrlQubitBit, 0}, NrBasisStates);
		}

		static inline void ApplyThreeQubitsGate(const GateClass& gate, VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			DispatchGate<false, 3>(gate.getStructure(), registerStorage, gateMatrix, {qubitBit, qubitBit2, ctrlQubitBit}, NrBasisStates);
		}

		static inline void ApplyThreeQubitsGateOmp(const GateClass& gate, VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			DispatchGate<true, 3>(gate.getStructure(), registerStorage, gateMatrix, {qubitBit, qubitBit2, ctrlQubitBit}, NrBasisStates);
		}

	private:
		// Deposit the compact index into the free bit positions. A tile covers the lowest run of free bits, so its
		// inner loop has indices with a constant stride: 1, or 2^k if the k lowest bits (qubit 0 upwards) are fixed.
		struct IndexPlan
		{
			size_t lowMasks[3]{};
			unsigned fixedBits = 0;
			unsigned strideShift = 0;
			size_t tileSize = 0;
			IndexPlan(size_t mask)
			{
				while ((mask >> strideShift) & 1) ++strideShift;
				const size_t higherFixed = mask >> strideShift;
				tileSize = higherFixed ? std::min<size_t>(higherFixed & (~higherFixed + 1), 1024) : 1024;
				while (mask)
				{
					const size_t bit = mask & (~mask + 1);
					assert(fixedBits < 3);
					lowMasks[fixedBits++] = bit - 1;
					mask ^= bit;
				}
			}
			template<unsigned Bits> size_t Expand(size_t index) const
			{
				if constexpr (Bits > 0) index = (index & lowMasks[0]) | ((index & ~lowMasks[0]) << 1);
				if constexpr (Bits > 1) index = (index & lowMasks[1]) | ((index & ~lowMasks[1]) << 1);
				if constexpr (Bits > 2) index = (index & lowMasks[2]) | ((index & ~lowMasks[2]) << 1);
				return index;
			}
		};

		// Direct access to the amplitudes; for a contiguous vector the stride is known at compile time. With GCC,
		// auto-vectorized loops that went through Eigen's operator() were measured to be up to 12 times slower.
		struct Amplitudes
		{
			using Scalar = typename VectorClass::Scalar;
			static constexpr bool contiguous = VectorClass::InnerStrideAtCompileTime == 1;

			explicit Amplitudes(VectorClass& vector) : data(vector.data()), stride(vector.innerStride()) {}

			Scalar& operator()(size_t i) const
			{
				if constexpr (contiguous) return data[i];
				else return data[static_cast<Eigen::Index>(i) * stride];
			}

			Scalar* data;
			Eigen::Index stride;
		};

		// The same rule for every compiler and OpenMP runtime: a parallel kernel always uses the whole team.
		// Whether to go parallel at all is decided per register, by its size (see GetParallelMinBasisStates), so
		// a circuit never alternates between serial and parallel kernels, nor changes the team size. Both were
		// measured to be very costly: changing the team size between regions made mixed circuits up to hundreds
		// of times slower with GCC's runtime, and serial kernels run while the team spin-waits are slowed down.
		// Nested callers (already in a parallel region) run serially, to avoid oversubscription.
		static int TeamThreads()
		{
#ifdef _OPENMP
			return omp_in_parallel() ? 1 : omp_get_max_threads();
#else
			return 1;
#endif
		}

		template<bool Parallel, class Function>
		static inline void ForEachTile(size_t tiles, const Function& function)
		{
			if constexpr (Parallel)
			{
				const int threads = TeamThreads();
				if (threads > 1)
				{
#pragma omp parallel for schedule(static) num_threads(threads)
					for (long long tile = 0; tile < static_cast<long long>(tiles); ++tile)
						function(static_cast<size_t>(tile));
					return;
				}
			}
			for (size_t tile = 0; tile < tiles; ++tile) function(tile);
		}

		template<bool Parallel, unsigned Bits, class Function>
		static inline void RunSelected(const IndexPlan& plan, size_t count, size_t oneMask,
			const Function& function)
		{
			const size_t active = count >> Bits;
			size_t tileSize = std::min(plan.tileSize, active);
			if constexpr (Parallel)
			{
				const int threads = TeamThreads();
				if (threads > 1)
					while (tileSize > 1 && active / tileSize < static_cast<size_t>(threads) * 4) tileSize >>= 1;
			}
			if (tileSize == 1)
				ForEachTile<Parallel>(active, [&](size_t i) { function(plan.template Expand<Bits>(i) | oneMask); });
			else if (plan.strideShift == 0)
				ForEachTile<Parallel>(active / tileSize, [&](size_t tile) {
					const size_t base = plan.template Expand<Bits>(tile * tileSize) | oneMask;
					for (size_t i = 0; i < tileSize; ++i) function(base + i);
				});
			else
				ForEachTile<Parallel>(active / tileSize, [&](size_t tile) {
					const size_t base = plan.template Expand<Bits>(tile * tileSize) | oneMask;
					const unsigned shift = plan.strideShift;
					for (size_t i = 0; i < tileSize; ++i) function(base + (i << shift));
				});
		}

		template<bool Parallel, class Function>
		static inline void ForEachSelected(size_t count, size_t zeroMask, size_t oneMask,
			const Function& function)
		{
			const IndexPlan plan(zeroMask | oneMask);
			switch (plan.fixedBits)
			{
			case 0: RunSelected<Parallel, 0>(plan, count, oneMask, function); break;
			case 1: RunSelected<Parallel, 1>(plan, count, oneMask, function); break;
			case 2: RunSelected<Parallel, 2>(plan, count, oneMask, function); break;
			case 3: RunSelected<Parallel, 3>(plan, count, oneMask, function); break;
			}
		}

		template<bool Parallel>
		static inline void ScaleSelected(const Amplitudes& state, size_t count, size_t zeroMask, size_t oneMask, const std::complex<double>& value)
		{
			if (value == std::complex<double>(1., 0.)) return;
			if (value == std::complex<double>(-1., 0.))
				ForEachSelected<Parallel>(count, zeroMask, oneMask, [&](size_t i) { state(i) = -state(i); });
			else if (value == std::complex<double>(0., 1.))
				ForEachSelected<Parallel>(count, zeroMask, oneMask, [&](size_t i) {
					const auto x = state(i); state(i) = std::complex<double>(-x.imag(), x.real());
				});
			else if (value == std::complex<double>(0., -1.))
				ForEachSelected<Parallel>(count, zeroMask, oneMask, [&](size_t i) {
					const auto x = state(i); state(i) = std::complex<double>(x.imag(), -x.real());
				});
			else
				ForEachSelected<Parallel>(count, zeroMask, oneMask, [&](size_t i) { state(i) *= value; });
		}

#ifdef QCSIM_SSE2
		// Complex coefficient times amplitude v = [re, im], with the coefficient prepared once as
		// re = [cr, cr] and im = [-ci, ci]: [cr * vr - ci * vi, cr * vi + ci * vr].
		static inline __m128d MultiplyByCoefficient(__m128d re, __m128d im, __m128d v)
		{
			return _mm_add_pd(_mm_mul_pd(re, v), _mm_mul_pd(im, _mm_shuffle_pd(v, v, 1)));
		}
#endif

		struct ControlMasks { size_t zeros = 0, ones = 0; };

#ifdef QCSIM_AVX2
		struct WideCoefficient
		{
			__m256d re, im;
			WideCoefficient(const std::complex<double>& first, const std::complex<double>& second)
				: re(_mm256_set_pd(second.real(), second.real(), first.real(), first.real())),
				  im(_mm256_set_pd(second.imag(), -second.imag(), first.imag(), -first.imag())) {}
			__m256d Multiply(__m256d v) const
			{
				return _mm256_add_pd(_mm256_mul_pd(re, v), _mm256_mul_pd(im, _mm256_permute_pd(v, 5)));
			}
		};

		template<bool Parallel>
		static inline bool TryWideDenseGate(const Amplitudes& state, size_t count, size_t target, ControlMasks controls,
			const std::complex<double>& a, const std::complex<double>& b,
			const std::complex<double>& c, const std::complex<double>& d, bool hadamard)
		{
			const size_t mask = controls.zeros | controls.ones;
			// For a higher target, two consecutive pairs can share a vector only
			// when bit zero is free. Limit to one control so the index plan still
			// needs at most three fixed bits (target, control, and batching bit).
			if (target != 1 && ((mask & 1) || (mask & (mask - 1)))) return false;
			if (hadamard)
			{
				const __m256d scale = _mm256_set1_pd(a.real());
				if (target == 1)
				{
					const __m256d upperSign = _mm256_set_pd(-0., -0., 0., 0.);
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						double* const p = reinterpret_cast<double*>(&state(i));
						const __m256d v = _mm256_loadu_pd(p), swap = _mm256_permute2f128_pd(v, v, 1);
						_mm256_storeu_pd(p, _mm256_mul_pd(_mm256_add_pd(_mm256_xor_pd(v, upperSign), swap), scale));
					});
				}
				else
					ForEachSelected<Parallel>(count, target | controls.zeros | 1, controls.ones, [&](size_t i) {
						double* const p = reinterpret_cast<double*>(&state(i)), * const q = reinterpret_cast<double*>(&state(i | target));
						const __m256d x = _mm256_loadu_pd(p), y = _mm256_loadu_pd(q);
						_mm256_storeu_pd(p, _mm256_mul_pd(_mm256_add_pd(x, y), scale));
						_mm256_storeu_pd(q, _mm256_mul_pd(_mm256_sub_pd(x, y), scale));
					});
			}
			else if (target == 1)
			{
				const WideCoefficient direct(a, d), swapped(b, c);
				ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
					double* const p = reinterpret_cast<double*>(&state(i));
					const __m256d v = _mm256_loadu_pd(p);
					_mm256_storeu_pd(p, _mm256_add_pd(direct.Multiply(v), swapped.Multiply(_mm256_permute2f128_pd(v, v, 1))));
				});
			}
			else
			{
				const WideCoefficient av(a, a), bv(b, b), cv(c, c), dv(d, d);
				ForEachSelected<Parallel>(count, target | controls.zeros | 1, controls.ones, [&](size_t i) {
					double* const p = reinterpret_cast<double*>(&state(i)), * const q = reinterpret_cast<double*>(&state(i | target));
					const __m256d x = _mm256_loadu_pd(p), y = _mm256_loadu_pd(q);
					_mm256_storeu_pd(p, _mm256_add_pd(av.Multiply(x), bv.Multiply(y)));
					_mm256_storeu_pd(q, _mm256_add_pd(cv.Multiply(x), dv.Multiply(y)));
				});
			}
			return true;
		}
#endif

		template<bool Parallel, class Operator>
		static inline void ApplySingleTargetGate(const Amplitudes& state, const Operator& matrix,
			size_t target, ControlMasks controls, size_t count, size_t offset, bool diagonal, bool antidiagonal)
		{
			if (diagonal)
			{
				const auto a = matrix(offset, offset), d = matrix(offset + 1, offset + 1);
				if (a == std::complex<double>(1., 0.))
					ScaleSelected<Parallel>(state, count, controls.zeros, controls.ones | target, d);
				else if (d == std::complex<double>(1., 0.))
					ScaleSelected<Parallel>(state, count, target | controls.zeros, controls.ones, a);
				else
					ForEachSelected<Parallel>(count, controls.zeros, controls.ones, [&](size_t i) { state(i) *= (i & target) ? d : a; });
			}
			else if (antidiagonal)
			{
				const auto b = matrix(offset, offset + 1), c = matrix(offset + 1, offset);
				if (b == std::complex<double>(1., 0.) && c == b)
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) { std::swap(state(i), state(i | target)); });
				else if (b == std::complex<double>(0., -1.) && c == -b)
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						const auto x = state(i), y = state(i | target);
						state(i) = std::complex<double>(y.imag(), -y.real());
						state(i | target) = std::complex<double>(-x.imag(), x.real());
					});
				else if (b == std::complex<double>(0., 1.) && c == -b)
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						const auto x = state(i), y = state(i | target);
						state(i) = std::complex<double>(-y.imag(), y.real());
						state(i | target) = std::complex<double>(x.imag(), -x.real());
					});
				else
				{
#ifdef QCSIM_SSE2
					// Whole 16 byte loads and stores, as in the Hadamard kernel below: GCC's code for this loop
					// stalled on reassembling the complex values, 7 times slower than with MSVC.
					const __m128d bRe = _mm_set1_pd(b.real()), bIm = _mm_set_pd(b.imag(), -b.imag());
					const __m128d cRe = _mm_set1_pd(c.real()), cIm = _mm_set_pd(c.imag(), -c.imag());
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						double* const p = reinterpret_cast<double*>(&state(i));
						double* const q = reinterpret_cast<double*>(&state(i | target));
						const __m128d x = _mm_loadu_pd(p), y = _mm_loadu_pd(q);
						_mm_storeu_pd(p, MultiplyByCoefficient(bRe, bIm, y));
						_mm_storeu_pd(q, MultiplyByCoefficient(cRe, cIm, x));
					});
#else
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						const auto x = state(i), y = state(i | target);
						state(i) = b * y;
						state(i | target) = c * x;
					});
#endif
				}
			}
			else
			{
				const auto a = matrix(offset, offset), b = matrix(offset, offset + 1);
				const auto c = matrix(offset + 1, offset), d = matrix(offset + 1, offset + 1);
				const bool hadamard = a.imag() == 0. && a == b && a == c && d == -a;
#ifdef QCSIM_AVX2
				if constexpr (Amplitudes::contiguous && std::is_same_v<typename Amplitudes::Scalar, std::complex<double>>)
					if (count >= 256 && TryWideDenseGate<Parallel>(state, count, target, controls, a, b, c, d, hadamard)) return;
#endif
				if (hadamard)
				{
					const double scale = a.real();
#ifdef QCSIM_SSE2
					// Each amplitude is loaded and stored as one 16 byte unit. Otherwise GCC assembles the values from
					// two 8 byte stack stores and a 16 byte reload, which cannot be forwarded: 2.5 times slower.
					const __m128d s = _mm_set1_pd(scale);
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						double* const p = reinterpret_cast<double*>(&state(i));
						double* const q = reinterpret_cast<double*>(&state(i | target));
						const __m128d x = _mm_loadu_pd(p), y = _mm_loadu_pd(q);
						_mm_storeu_pd(p, _mm_mul_pd(_mm_add_pd(x, y), s));
						_mm_storeu_pd(q, _mm_mul_pd(_mm_sub_pd(x, y), s));
					});
#else
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						const auto x = state(i), y = state(i | target);
						state(i) = (x + y) * scale;
						state(i | target) = (x - y) * scale;
					});
#endif
				}
				else
				{
#ifdef QCSIM_SSE2
					const __m128d aRe = _mm_set1_pd(a.real()), aIm = _mm_set_pd(a.imag(), -a.imag());
					const __m128d bRe = _mm_set1_pd(b.real()), bIm = _mm_set_pd(b.imag(), -b.imag());
					const __m128d cRe = _mm_set1_pd(c.real()), cIm = _mm_set_pd(c.imag(), -c.imag());
					const __m128d dRe = _mm_set1_pd(d.real()), dIm = _mm_set_pd(d.imag(), -d.imag());
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						double* const p = reinterpret_cast<double*>(&state(i));
						double* const q = reinterpret_cast<double*>(&state(i | target));
						const __m128d x = _mm_loadu_pd(p), y = _mm_loadu_pd(q);
						_mm_storeu_pd(p, _mm_add_pd(MultiplyByCoefficient(aRe, aIm, x), MultiplyByCoefficient(bRe, bIm, y)));
						_mm_storeu_pd(q, _mm_add_pd(MultiplyByCoefficient(cRe, cIm, x), MultiplyByCoefficient(dRe, dIm, y)));
					});
#else
					ForEachSelected<Parallel>(count, target | controls.zeros, controls.ones, [&](size_t i) {
						const auto x = state(i), y = state(i | target);
						state(i) = a * x + b * y;
						state(i | target) = c * x + d * y;
					});
#endif
				}
			}
		}

		template<bool Parallel>
		static inline void ApplySwapPairs(const Amplitudes& state, size_t first, size_t second, ControlMasks controls,
			size_t count, const std::complex<double>& phase)
		{
			const size_t flip = first | second;
			if (phase == std::complex<double>(1., 0.))
				ForEachSelected<Parallel>(count, first | controls.zeros, second | controls.ones, [&](size_t i) { std::swap(state(i), state(i ^ flip)); });
			else if (phase == std::complex<double>(0., 1.))
				ForEachSelected<Parallel>(count, first | controls.zeros, second | controls.ones, [&](size_t i) {
					const auto x = state(i), y = state(i ^ flip);
					state(i) = std::complex<double>(-y.imag(), y.real());
					state(i ^ flip) = std::complex<double>(-x.imag(), x.real());
				});
			else
				ForEachSelected<Parallel>(count, first | controls.zeros, second | controls.ones, [&](size_t i) {
					const auto x = state(i), y = state(i ^ flip);
					state(i) = std::complex<double>(y.imag(), -y.real());
					state(i ^ flip) = std::complex<double>(x.imag(), -x.real());
				});
		}

		// Save every input in a coupled block before any output is written.
		// Matrix entries are copied once per gate, outside the amplitude loop.
		template<unsigned Bits, bool Parallel, class Operator>
		static inline void ApplyDenseGate(const Amplitudes& state, const Operator& matrix,
			const std::array<size_t, 3>& bits, ControlMasks controls, size_t count, size_t offset)
		{
			constexpr size_t dimension = size_t{1} << Bits;
			using Scalar = typename MatrixClass::Scalar;
			std::array<size_t, dimension> offsets{};
			std::array<Scalar, dimension * dimension> coefficients;
			size_t targets = 0;
			for (unsigned b = 0; b < Bits; ++b) targets |= bits[b];
			for (size_t r = 0; r < dimension; ++r)
			{
				for (unsigned b = 0; b < Bits; ++b)
					if (r & (size_t{1} << b)) offsets[r] |= bits[b];
				for (size_t c = 0; c < dimension; ++c)
					coefficients[r * dimension + c] = matrix(offset + r, offset + c);
			}
			ForEachSelected<Parallel>(count, targets | controls.zeros, controls.ones, [&](size_t base) {
				constexpr size_t blockSize = size_t{1} << Bits;
				std::array<Scalar, blockSize> original;
				for (size_t c = 0; c < blockSize; ++c) original[c] = state(base | offsets[c]);
				for (size_t r = 0; r < blockSize; ++r)
				{
					Scalar value = coefficients[r * blockSize] * original[0];
					for (size_t c = 1; c < blockSize; ++c) value += coefficients[r * blockSize + c] * original[c];
					state(base | offsets[r]) = value;
				}
			});
		}

		// Multiplies each control-selected amplitude by its diagonal entry in a single pass, unless scaling only the
		// entries that differ from one, each in its own pass, touches fewer cache lines. A 64 byte cache line holds
		// 4 amplitudes, so a target on qubit 0 or 1 does not reduce the lines a pass touches. With all targets on
		// qubit 2 or above the passes never touch more lines than a single pass, and their contiguous loops are faster.
		template<unsigned Bits, bool Parallel, class Operator>
		static inline void ApplyDiagonalGate(const Amplitudes& state, const Operator& matrix,
			const std::array<size_t, 3>& bits, ControlMasks controls, size_t count, size_t offset)
		{
			constexpr size_t dimension = size_t{1} << Bits;
			using Scalar = typename MatrixClass::Scalar;
			std::array<Scalar, dimension> diagonal;
			size_t targets = 0;
			size_t notOne = 0;
			unsigned lineBits = 0;
			for (unsigned b = 0; b < Bits; ++b)
			{
				targets |= bits[b];
				if (bits[b] >= 4) ++lineBits;
			}
			for (size_t row = 0; row < dimension; ++row)
			{
				diagonal[row] = matrix(offset + row, offset + row);
				if (diagonal[row] != Scalar(1.)) ++notOne;
			}

			if (lineBits == Bits || notOne < (size_t{1} << lineBits))
			{
				for (size_t row = 0; row < dimension; ++row)
				{
					size_t selected = 0;
					for (unsigned b = 0; b < Bits; ++b)
						if (row & (size_t{1} << b)) selected |= bits[b];
					ScaleSelected<Parallel>(state, count, (targets ^ selected) | controls.zeros, controls.ones | selected, diagonal[row]);
				}
				return;
			}

			ForEachSelected<Parallel>(count, controls.zeros, controls.ones, [&](size_t i) {
				size_t local = 0;
				for (unsigned b = 0; b < Bits; ++b)
					if (i & bits[b]) local |= size_t{1} << b;
				state(i) *= diagonal[local];
			});
		}

		// Only the active block is mapped; coefficients are still loaded outside
		// the amplitude loops, from the supplied (possibly conjugated) matrix.
		struct ActiveMatrix
		{
			const MatrixClass& matrix;
			std::array<size_t, 8> indices{};
			auto operator()(size_t row, size_t col) const { return matrix(indices[row], indices[col]); }
		};

		template<unsigned Bits, bool Parallel>
		static inline void ApplyPermutationGate(const Amplitudes& state, const ActiveMatrix& matrix,
			const Gates::GateStructure& structure, const std::array<size_t, 3>& bits, ControlMasks controls, size_t count)
		{
			constexpr size_t dimension = size_t{1} << Bits;
			using Scalar = typename MatrixClass::Scalar;
			std::array<size_t, dimension> offsets{}, outputs{}, sources{};
			std::array<Scalar, dimension> coefficients;
			size_t targets = 0, changed = 0;
			bool allOne = true;
			for (unsigned b = 0; b < Bits; ++b) targets |= bits[b];
			for (size_t row = 0; row < dimension; ++row)
				for (unsigned b = 0; b < Bits; ++b)
					if (row & (size_t{1} << b)) offsets[row] |= bits[b];
			for (unsigned row = 0; row < dimension; ++row)
			{
				const unsigned source = structure.kind == Gates::GateStructure::Kind::Antidiagonal
					? unsigned(dimension - 1 - row) : structure.Source(row);
				const Scalar coefficient = matrix(row, source);
				if (row == source && coefficient == Scalar(1.)) continue;
				outputs[changed] = offsets[row]; sources[changed] = offsets[source];
				coefficients[changed++] = coefficient;
				allOne = allOne && coefficient == Scalar(1.);
			}
			if (changed == 0) return;
			if (allOne)
			{
				ForEachSelected<Parallel>(count, targets | controls.zeros, controls.ones, [&](size_t base) {
					std::array<Scalar, size_t{1} << Bits> original;
					for (size_t r = 0; r < changed; ++r) original[r] = state(base | sources[r]);
					for (size_t r = 0; r < changed; ++r) state(base | outputs[r]) = original[r];
				});
				return;
			}
#ifdef QCSIM_SSE2
			__m128d re[dimension], im[dimension];
			for (size_t r = 0; r < changed; ++r)
			{
				re[r] = _mm_set1_pd(coefficients[r].real());
				im[r] = _mm_set_pd(coefficients[r].imag(), -coefficients[r].imag());
			}
			ForEachSelected<Parallel>(count, targets | controls.zeros, controls.ones, [&](size_t base) {
				__m128d original[size_t{1} << Bits];
				for (size_t r = 0; r < changed; ++r) original[r] = _mm_loadu_pd(reinterpret_cast<const double*>(&state(base | sources[r])));
				for (size_t r = 0; r < changed; ++r)
					_mm_storeu_pd(reinterpret_cast<double*>(&state(base | outputs[r])), MultiplyByCoefficient(re[r], im[r], original[r]));
			});
#else
			ForEachSelected<Parallel>(count, targets | controls.zeros, controls.ones, [&](size_t base) {
				std::array<Scalar, size_t{1} << Bits> original;
				for (size_t r = 0; r < changed; ++r) original[r] = state(base | sources[r]);
				for (size_t r = 0; r < changed; ++r) state(base | outputs[r]) = coefficients[r] * original[r];
			});
#endif
		}

		template<bool Parallel, unsigned GateQubits>
		static inline void DispatchGate(const Gates::GateStructure& structure, VectorClass& vector, const MatrixClass& matrix,
			const std::array<size_t, 3>& bits, size_t count)
		{
			using Kind = Gates::GateStructure::Kind;
			const Amplitudes state(vector);
			if constexpr (GateQubits == 1)
			{
				ApplySingleTargetGate<Parallel>(state, matrix, bits[0], {}, count, 0,
					structure.kind == Kind::Diagonal, structure.kind == Kind::Antidiagonal);
				return;
			}
			const unsigned activeBits = GateQubits - structure.controls;
			ControlMasks controls;
			std::array<size_t, 3> targets{};
			unsigned targetIndex = 0;
			for (unsigned b = 0; b < GateQubits; ++b)
				if (structure.controlMask & (1u << b))
				{
					if (structure.controlValue & (1u << b)) controls.ones |= bits[b];
					else controls.zeros |= bits[b];
				}
				else targets[targetIndex++] = bits[b];
			ActiveMatrix active{matrix};
			if (activeBits == 1)
			{
				active.indices[0] = structure.controlValue;
				active.indices[1] = structure.controlValue | (((1u << GateQubits) - 1) ^ structure.controlMask);
				ApplySingleTargetGate<Parallel>(state, active, targets[0], controls, count, 0,
					structure.kind == Kind::Diagonal, structure.kind == Kind::Antidiagonal);
				return;
			}
			unsigned index = 0;
			for (unsigned i = 0; i < (1u << GateQubits); ++i)
				if ((i & structure.controlMask) == structure.controlValue) active.indices[index++] = i;
			if (structure.kind == Kind::Swap || structure.kind == Kind::ISwap || structure.kind == Kind::ISwapDag)
				ApplySwapPairs<Parallel>(state, targets[0], targets[1], controls, count, active(1, 2));
			else if (structure.kind == Kind::Diagonal)
			{
				if (activeBits == 2) ApplyDiagonalGate<2, Parallel>(state, active, targets, controls, count, 0);
				else ApplyDiagonalGate<3, Parallel>(state, active, targets, controls, count, 0);
			}
			else if (structure.kind == Kind::Antidiagonal || structure.kind == Kind::Permutation)
			{
				if (activeBits == 2) ApplyPermutationGate<2, Parallel>(state, active, structure, targets, controls, count);
				else ApplyPermutationGate<3, Parallel>(state, active, structure, targets, controls, count);
			}
			else if (activeBits == 2) ApplyDenseGate<2, Parallel>(state, active, targets, controls, count, 0);
			else ApplyDenseGate<3, Parallel>(state, active, targets, controls, count, 0);
		}

		// Measurements. The measured qubits firstQubit..secondQubit are contiguous, so the states with a given outcome
		// form runs of 2^firstQubit contiguous amplitudes, one run in each block of 2^(secondQubit + 1) amplitudes.

		template<bool Parallel, class Function>
		static inline double Sum(size_t count, const Function& term)
		{
			double sum = 0;
			if constexpr (Parallel)
			{
				const int threads = TeamThreads();
				if (threads > 1)
				{
#pragma omp parallel for reduction(+:sum) schedule(static) num_threads(threads)
					for (long long k = 0; k < static_cast<long long>(count); ++k)
						sum += term(static_cast<size_t>(k));
					return sum;
				}
			}
			for (size_t k = 0; k < count; ++k) sum += term(k);
			return sum;
		}

		// Returns the first state for which the cumulative probability reaches prob, as a sequential scan does.
		// The parallel version sums fixed size chunks first (independent of the number of threads), then scans
		// only the chunk where prob is reached. Their results differ only if prob is within rounding of a boundary.
		template<bool Parallel>
		static inline size_t FindSampledState(size_t count, const VectorClass& state, double prob, size_t notFound)
		{
			constexpr size_t chunkSize = 4096;
			double accum = 0;

			if constexpr (Parallel)
			{
				const size_t chunks = count / chunkSize;
				if (chunks > 1 && count % chunkSize == 0 && TeamThreads() > 1)
				{
					// Concentrated states (in particular |0>) should not pay for a
					// complete scan or a team launch. Continue the first chunk after
					// this prefix, without counting its probability twice.
					constexpr size_t prefixSize = 128;
					for (size_t i = 0; i < prefixSize; ++i)
					{
						accum += std::norm(state(i));
						if (prob <= accum) return i;
					}
					std::vector<double> chunkSums(chunks);
					ForEachTile<true>(chunks, [&](size_t chunk) {
						double sum = 0;
						for (size_t i = chunk ? chunk * chunkSize : prefixSize; i < (chunk + 1) * chunkSize; ++i)
							sum += std::norm(state(i));
						chunkSums[chunk] = sum;
					});

					for (size_t chunk = 0; chunk < chunks; ++chunk)
					{
						if (prob <= accum + chunkSums[chunk])
						{
							for (size_t i = chunk ? chunk * chunkSize : prefixSize; i < (chunk + 1) * chunkSize; ++i)
							{
								accum += std::norm(state(i));
								if (prob <= accum) return i;
							}
						}
						else
							accum += chunkSums[chunk];
					}

					return notFound;
				}
			}

			for (size_t i = 0; i < count; ++i)
			{
				accum += std::norm(state(i));
				if (prob <= accum) return i;
			}

			return notFound;
		}

		template<bool Parallel>
		static inline double OutcomeProbability(size_t count, const VectorClass& state, size_t firstQubit, size_t secondQubit, size_t outcome)
		{
			const size_t lowMask = (size_t{1} << firstQubit) - 1;
			const size_t selected = outcome << firstQubit;
			const size_t shift = secondQubit + 1 - firstQubit;

			return Sum<Parallel>(count >> shift, [&](size_t k) {
				return std::norm(state(((k & ~lowMask) << shift) | selected | (k & lowMask)));
			});
		}

		template<bool Parallel>
		static inline void Collapse(size_t count, VectorClass& state, size_t measuredMask, size_t selected, double norm)
		{
			const size_t tileSize = std::min<size_t>(count, 1024);
			ForEachTile<Parallel>(count / tileSize, [&](size_t tile) {
				for (size_t i = tile * tileSize; i < (tile + 1) * tileSize; ++i)
					state(i) = (i & measuredMask) == selected ? state(i) * norm : std::complex<double>(0., 0.);
			});
		}

		static inline size_t MeasuredMask(size_t firstQubit, size_t secondQubit)
		{
			return ((size_t{1} << (secondQubit + 1)) - 1) & ~((size_t{1} << firstQubit) - 1);
		}

		template<bool Parallel>
		static inline size_t MeasureRange(size_t count, VectorClass& state, size_t firstQubit, size_t secondQubit, double prob)
		{
			const size_t measuredMask = MeasuredMask(firstQubit, secondQubit);
			const size_t outcome = (FindSampledState<Parallel>(count, state, prob, 0) & measuredMask) >> firstQubit;
			const double norm = 1. / sqrt(OutcomeProbability<Parallel>(count, state, firstQubit, secondQubit, outcome));
			Collapse<Parallel>(count, state, measuredMask, outcome << firstQubit, norm);

			return outcome;
		}

	public:
		//*****************************************************************************************************************************************************************************************
		// 
		//  Measurements
		// 
		//*****************************************************************************************************************************************************************************************

		static inline size_t MeasureQubit(size_t NrBasisStates, VectorClass& registerStorage, size_t qubit, const double prob)
		{
			return MeasureRange<false>(NrBasisStates, registerStorage, qubit, qubit, prob);
		}

		static inline size_t MeasureQubitOmp(size_t NrBasisStates, VectorClass& registerStorage, size_t qubit, const double prob)
		{
			return MeasureRange<true>(NrBasisStates, registerStorage, qubit, qubit, prob);
		}

		static inline size_t MeasureQubitNoCollapse(size_t NrBasisStates, VectorClass& registerStorage, size_t qubit, const double prob)
		{
			return (FindSampledState<false>(NrBasisStates, registerStorage, prob, 0) >> qubit) & 1;
		}

		static inline double GetQubitProbability(size_t NrBasisStates, const VectorClass& registerStorage, size_t qubit)
		{
			return OutcomeProbability<false>(NrBasisStates, registerStorage, qubit, qubit, 1);
		}

		static inline double GetQubitProbabilityOmp(size_t NrBasisStates, const VectorClass& registerStorage, size_t qubit)
		{
			return OutcomeProbability<true>(NrBasisStates, registerStorage, qubit, qubit, 1);
		}

		static inline size_t Measure(size_t NrBasisStates, VectorClass& registerStorage, size_t firstQubit, size_t secondQubit, const double prob)
		{
			return MeasureRange<false>(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}

		static inline size_t MeasureOmp(size_t NrBasisStates, VectorClass& registerStorage, size_t firstQubit, size_t secondQubit, const double prob)
		{
			return MeasureRange<true>(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}

		static inline size_t MeasureNoCollapse(size_t NrBasisStates, VectorClass& registerStorage, size_t firstQubit, size_t secondQubit, const double prob)
		{
			return (FindSampledState<false>(NrBasisStates, registerStorage, prob, 0) & MeasuredMask(firstQubit, secondQubit)) >> firstQubit;
		}

		// The basis state selected by prob (the first one where the cumulative probability reaches it), without collapsing.
		static inline size_t SampleBasisState(size_t NrBasisStates, const VectorClass& registerStorage, const double prob, size_t notFound, bool parallel)
		{
			if (parallel)
				return FindSampledState<true>(NrBasisStates, registerStorage, prob, notFound);

			return FindSampledState<false>(NrBasisStates, registerStorage, prob, notFound);
		}

		// <psi|P|psi> for a Pauli string P in a single read-only pass, without copying the state.
		// X and Y qubits are set in xMask, Z and Y qubits in zMask. Since Y = -i * Z * X:
		// (P psi)(i) = (-i)^nrY * (-1)^parity(i & zMask) * psi(i ^ xMask)
		static inline std::complex<double> PauliExpectationValue(size_t NrBasisStates, const VectorClass& registerStorage, size_t xMask, size_t zMask, size_t nrY, bool parallel)
		{
			static const std::complex<double> minusIPowers[4] = { {1., 0.}, {0., -1.}, {-1., 0.}, {0., 1.} };
			const auto phase = minusIPowers[nrY & 3];
			const auto sum = [&](size_t count, const auto& term) {
				return parallel ? Sum<true>(count, term) : Sum<false>(count, term);
			};
			if (xMask == 0)
			{
				const double value = sum(NrBasisStates, [&](size_t i) {
					const double norm = std::norm(registerStorage(i));
					return OddParity(i & zMask) ? -norm : norm;
				});
				return phase * value;
			}

			// Each pair i,j=i^xMask contributes t+conj(t), or t-conj(t).
			// Enumerate just the member with one selected X/Y bit clear. This
			// halves amplitude loads and needs only one real-valued reduction.
			const size_t lowMask = (xMask & (~xMask + 1)) - 1;
			const auto pairs = [&](auto imaginary) {
				return sum(NrBasisStates / 2, [&](size_t k) {
					const size_t i = (k & lowMask) | ((k & ~lowMask) << 1);
					const auto x = registerStorage(i), y = registerStorage(i ^ xMask);
					double value;
					if constexpr (decltype(imaginary)::value) value = 2. * (x.real() * y.imag() - x.imag() * y.real());
					else value = 2. * (x.real() * y.real() + x.imag() * y.imag());
					return OddParity(i & zMask) ? -value : value;
				});
			};
			if (OddParity(xMask & zMask)) return phase * std::complex<double>(0., pairs(std::true_type{}));
			return phase * pairs(std::false_type{});
		}

		static int GetNumberOfThreads()
		{
			const size_t threads = std::thread::hardware_concurrency();
			return static_cast<int>(threads ? threads : GetCpuInfoNrThreads());
		}


		void SetMultithreading(bool enable = true)
		{
			enableMultithreading = enable;
		}

		bool GetMultithreading() const
		{
			return enableMultithreading;
		}

		// With multithreading enabled, registers with at least this many basis states run all their gates and
		// measurements with the whole OpenMP team, smaller ones run serially. The density matrix applies it to the
		// length of its rows and columns. It is a process-wide setting, to be tuned for the machine if needed.
		// Simulators running concurrently on several threads should have multithreading disabled instead.
		static void SetParallelMinBasisStates(size_t nrBasisStates)
		{
			ParallelMinBasisStatesSetting().store(nrBasisStates, std::memory_order_relaxed);
		}

		static size_t GetParallelMinBasisStates()
		{
			return ParallelMinBasisStatesSetting().load(std::memory_order_relaxed);
		}

	private:
		static inline bool OddParity(uint64_t v)
		{
			v ^= v >> 32;
			v ^= v >> 16;
			v ^= v >> 8;
			v ^= v >> 4;
			v ^= v >> 2;
			v ^= v >> 1;

			return (v & 1) != 0;
		}

		static size_t GetCpuInfoNrThreads()
		{
#ifdef _WIN32
			SYSTEM_INFO sysinfo;
			GetSystemInfo(&sysinfo);
			return sysinfo.dwNumberOfProcessors;
#else
			std::ifstream cpuinfo("/proc/cpuinfo");

			return std::count(std::istream_iterator<std::string>(cpuinfo), std::istream_iterator<std::string>(), std::string("processor"));
#endif
		}

		bool enableMultithreading = true;
	};

}
