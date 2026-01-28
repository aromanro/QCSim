#pragma once

#include "PauliStringXZ.h"

namespace QC
{

	class PauliStringXZWithCoefficient : public PauliStringXZ {
	public:
		PauliStringXZWithCoefficient() : Coefficient(1.0) {}
		PauliStringXZWithCoefficient(size_t nQubits) : PauliStringXZ(nQubits), Coefficient(1.0) {}
		PauliStringXZWithCoefficient(const PauliStringXZWithCoefficient& other) : PauliStringXZ(other), Coefficient(other.Coefficient) {}
		PauliStringXZWithCoefficient(PauliStringXZWithCoefficient&& other) noexcept : PauliStringXZ(std::move(other)), Coefficient(other.Coefficient) {}
		PauliStringXZWithCoefficient& operator=(const PauliStringXZWithCoefficient& other)
		{
			if (this != &other)
			{
				PauliStringXZ::operator=(other);
				Coefficient = other.Coefficient;
			}
			return *this;
		}
		PauliStringXZWithCoefficient& operator=(PauliStringXZWithCoefficient&& other) noexcept
		{
			if (this != &other)
			{
				Coefficient = other.Coefficient;
				PauliStringXZ::operator=(std::move(other));
			}
			return *this;
		}
		void Clear() override
		{
			PauliStringXZ::Clear();
			Coefficient = 1.0;
		}

		bool operator==(const PauliStringXZWithCoefficient& other) const
		{
			return (X == other.X) && (Z == other.Z);
		}

		double ExpectationValue() const override
		{
			const double expval = PauliStringXZ::ExpectationValue();

			return Coefficient * expval;
		}

		inline void ApplyH(size_t qubit) override
		{
			// how does this work?
			// looking at it might not reveal immediately what it does
			// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

			// if we have both X and Z, then it's a Y (with some global phase, given by the sign, Y = iXZ)
			// a Y is transformed to a -Y, so a sign change is needed

			if (X[qubit] && Z[qubit])
				Coefficient = -Coefficient;

			// swap X and Z
			PauliStringXZ::ApplyH(qubit);
		}

		inline void ApplyS(size_t qubit) override
		{
			if (X[qubit] && Z[qubit])
				Coefficient = -Coefficient;

			PauliStringXZ::ApplyS(qubit);
		}

		inline void ApplyX(size_t qubit) override
		{
			// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
			if (Z[qubit])
				Coefficient = -Coefficient;
		}

		inline void ApplyY(size_t qubit) override
		{
			// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
			// if both X and Z are present, the sign is flipped twice, so it remains unchanged

			// can be done with ifs, can be done with XORs
			if (XOR(Z[qubit], X[qubit]))
				Coefficient = -Coefficient;
		}

		inline void ApplyZ(size_t qubit) override
		{
			// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as ZXZ^t = -X
			if (X[qubit])
				Coefficient = -Coefficient;
		}

		inline void ApplyCX(size_t target, size_t control) override
		{
			if (X[control] && Z[target] && XOR(X[target], !Z[control]))
				Coefficient = -Coefficient;

			PauliStringXZ::ApplyCX(target, control);
		}

		std::string ToString() const override
		{
			std::string result = PauliStringXZ::ToString();
			return std::to_string(Coefficient) + " * " + result;
		}

		mutable double Coefficient;
	};

}
