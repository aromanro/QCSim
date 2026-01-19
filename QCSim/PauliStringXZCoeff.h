#pragma once

#include "PauliStringXZ.h"

namespace QC
{

	class PauliStringXZWithCoefficient : public PauliStringXZ {
	public:
		PauliStringXZWithCoefficient() : Coefficient(1.0) {}
		PauliStringXZWithCoefficient(size_t nQubits) : PauliStringXZ(nQubits), Coefficient(1.0) {}
		PauliStringXZWithCoefficient(const PauliStringXZWithCoefficient& other) : PauliStringXZ(other), Coefficient(other.Coefficient) {}
		PauliStringXZWithCoefficient(PauliStringXZWithCoefficient&& other) noexcept : PauliStringXZ(other), Coefficient(other.Coefficient) {}
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

		double ExpectationValue() const override
		{
			const double expval = PauliStringXZ::ExpectationValue();

			return Coefficient * expval;
		}

		inline void ApplyH(size_t qubit)
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

		inline void ApplyS(size_t qubit)
		{
			if (X[qubit] && Z[qubit])
				Coefficient = -Coefficient;

			PauliStringXZ::ApplyS(qubit);
		}

		inline void ApplyX(size_t qubit)
		{
			// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
			if (Z[qubit])
				Coefficient = -Coefficient;
		}

		inline void ApplyY(size_t qubit)
		{
			// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
			// if both X and Z are present, the sign is flipped twice, so it remains unchanged

			// can be done with ifs, can be done with XORs
			if (XOR(Z[qubit], X[qubit]))
				Coefficient = -Coefficient;
		}

		inline void ApplyZ(size_t qubit)
		{
			// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as ZXZ^t = -X
			if (X[qubit])
				Coefficient = -Coefficient;
		}

		inline void ApplyCX(size_t target, size_t control)
		{
			if (X[control] && Z[target] && XOR(X[target], !Z[control]))
				Coefficient = -Coefficient;

			PauliStringXZ::ApplyCX(target, control);
		}

		inline void ApplyK(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		inline void ApplySdg(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
		}

		inline void ApplySx(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyZ(qubit);
			ApplyS(qubit);
		}

		inline void ApplySxDag(size_t qubit)
		{
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		inline void ApplyCY(size_t target, size_t control)
		{
			ApplyZ(target);
			ApplyS(target);
			ApplyCX(target, control);
			ApplyS(target);
		}

		inline void ApplyCZ(size_t target, size_t control)
		{
			ApplyH(target);
			ApplyCX(target, control);
			ApplyH(target);
		}

		inline void ApplySwap(size_t qubit1, size_t qubit2)
		{
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
		}

		inline void ApplyISwap(size_t qubit1, size_t qubit2)
		{
			ApplyS(qubit1);
			ApplyH(qubit1);
			ApplyS(qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
			ApplyH(qubit2);
		}

		inline void ApplyISwapDag(size_t qubit1, size_t qubit2)
		{
			ApplyH(qubit2);
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyZ(qubit2);
			ApplyS(qubit2);
			ApplyH(qubit1);
			ApplyZ(qubit1);
			ApplyS(qubit1);
		}

		double Coefficient;
	};

}
