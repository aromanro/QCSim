#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>

#include "BitOps.h"

namespace QC { namespace detail {

	// True when angle is a multiple of pi/2, which quarterTurns receives. The
	// one rule for every simulator: angle must equal the rounded product
	// turns * (pi/2), so k * M_PI_2 and similar expressions count as exact.
	inline bool TryGetQuarterTurns(double angle, long long& quarterTurns)
	{
		constexpr double quarterTurnAngle = 1.57079632679489661923132169163975144;
		const double turns = angle / quarterTurnAngle;
		if (!(std::abs(turns) < static_cast<double>(std::numeric_limits<long long>::max())))
			return false;
		const auto roundedTurns = static_cast<long long>(std::llround(turns));
		if (angle != static_cast<double>(roundedTurns) * quarterTurnAngle)
			return false;
		quarterTurns = roundedTurns;
		return true;
	}

	// cos and sin of angle, exact at multiples of pi/2. Returns whether angle is
	// such a multiple, by the rule of TryGetQuarterTurns.
	inline bool QuarterTurnCosSin(double angle, double& cosine, double& sine)
	{
		long long quarterTurns = 0;
		if (TryGetQuarterTurns(angle, quarterTurns))
		{
			static constexpr double values[4][2] = { { 1., 0. }, { 0., 1. }, { -1., 0. }, { 0., -1. } };
			const auto& value = values[((quarterTurns % 4) + 4) % 4];
			cosine = value[0];
			sine = value[1];
			return true;
		}
		cosine = std::cos(angle);
		sine = std::sin(angle);
		return false;
	}

	// A Pauli string on the at most three qubits of a gate. Local qubit j sets
	// bit j of x and z; both bits mean Y, with the convention
	// P(x, z) = i^|x&z| X^x Z^z of the packed tableaux.
	struct LocalPauli {
		static constexpr size_t MaxQubits = 3;

		uint8_t x = 0;
		uint8_t z = 0;

		size_t Key() const noexcept { return x | (size_t(z) << MaxQubits); }

		static LocalPauli FromKey(size_t key) noexcept
		{
			return { uint8_t(key & ((1U << MaxQubits) - 1)), uint8_t(key >> MaxQubits) };
		}

		static LocalPauli Single(size_t qubit, char pauli) noexcept
		{
			const uint8_t bit = uint8_t(1U << qubit);
			return { uint8_t(pauli != 'Z' ? bit : 0), uint8_t(pauli != 'X' ? bit : 0) };
		}

		// 'I', 'X', 'Y' or 'Z' on local qubit j.
		char Factor(size_t qubit) const noexcept
		{
			const bool hasX = (x >> qubit) & 1U, hasZ = (z >> qubit) & 1U;
			return hasX ? (hasZ ? 'Y' : 'X') : (hasZ ? 'Z' : 'I');
		}
	};

	// An operator on at most three qubits as a combination of Pauli strings.
	// Products and sums of exactly representable coefficients stay exact.
	class LocalPauliSum {
	public:
		using Complex = std::complex<double>;
		static constexpr size_t Keys = size_t(1) << (2 * LocalPauli::MaxQubits);

		// The zero operator.
		LocalPauliSum() = default;

		static LocalPauliSum Identity(Complex coefficient = 1.0)
		{
			LocalPauliSum sum;
			sum.Add(LocalPauli{}, coefficient);
			return sum;
		}

		static LocalPauliSum Pauli(size_t qubit, char pauli, Complex coefficient = 1.0)
		{
			LocalPauliSum sum;
			sum.Add(LocalPauli::Single(qubit, pauli), coefficient);
			return sum;
		}

		// (I + Z)/2 or (I - Z)/2: the projector onto |0> or |1> of a local qubit.
		static LocalPauliSum Projector(size_t qubit, bool one)
		{
			LocalPauliSum sum = Identity(0.5);
			sum.Add(LocalPauli::Single(qubit, 'Z'), one ? -0.5 : 0.5);
			return sum;
		}

		void Add(LocalPauli pauli, Complex coefficient) noexcept
		{
			const size_t key = pauli.Key();
			coefficients[key] += coefficient;
			if (coefficients[key] == Complex(0.0)) nonZero &= ~(uint64_t(1) << key);
			else nonZero |= uint64_t(1) << key;
		}

		LocalPauliSum operator+(const LocalPauliSum& other) const noexcept
		{
			LocalPauliSum sum(*this);
			other.ForEachTerm([&sum](LocalPauli pauli, Complex coefficient) { sum.Add(pauli, coefficient); });
			return sum;
		}

		LocalPauliSum operator*(const LocalPauliSum& right) const noexcept
		{
			LocalPauliSum product;
			ForEachTerm([&](LocalPauli a, Complex left) {
				right.ForEachTerm([&](LocalPauli b, Complex coefficient) {
					product.Add({ uint8_t(a.x ^ b.x), uint8_t(a.z ^ b.z) },
						left * coefficient * IPower(ProductPhase(a, b)));
				});
			});
			return product;
		}

		// Only the identity term is nonzero: a scalar multiple of I.
		bool IsScalar() const noexcept { return (nonZero & ~uint64_t(1)) == 0; }

		size_t Terms() const noexcept { return Clifford::detail::Popcount(nonZero); }

		template<class Function> void ForEachTerm(Function function) const
		{
			for (uint64_t keys = nonZero; keys; keys &= keys - 1)
			{
				const size_t key = Clifford::detail::TrailingZero(keys);
				function(LocalPauli::FromKey(key), coefficients[key]);
			}
		}

		// i^k for k modulo 4, exactly.
		static Complex IPower(unsigned k) noexcept
		{
			switch (k & 3U)
			{
			case 0: return { 1.0, 0.0 };
			case 1: return { 0.0, 1.0 };
			case 2: return { -1.0, 0.0 };
			default: return { 0.0, -1.0 };
			}
		}

		// P(a) P(b) = i^k P(a xor b); returns k modulo 4.
		static unsigned ProductPhase(LocalPauli a, LocalPauli b) noexcept
		{
			const auto count = [](unsigned value) { return Clifford::detail::Popcount(value); };
			return (count(a.x & a.z) + count(b.x & b.z) + 2 * count(a.z & b.x)
				- count((a.x ^ b.x) & (a.z ^ b.z))) & 3U;
		}

	private:
		static_assert(Keys <= 64, "Nonzero terms are tracked in one 64-bit mask");

		std::array<Complex, Keys> coefficients{};
		uint64_t nonZero = 0;
	};

}}
