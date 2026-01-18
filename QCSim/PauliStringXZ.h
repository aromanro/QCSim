#pragma once

#include <vector>
#include <algorithm>

namespace QC
{

	class PauliStringXZ {
	public:
		PauliStringXZ() {}

		PauliStringXZ(size_t nQubits) : X(nQubits, false), Z(nQubits, false) {}

		PauliStringXZ(const PauliStringXZ& other) : X(other.X), Z(other.Z) {}

		PauliStringXZ(PauliStringXZ&& other) noexcept : X(std::move(other.X)), Z(std::move(other.Z)) {}

		virtual ~PauliStringXZ() = default;

		PauliStringXZ& operator=(const PauliStringXZ& other)
		{
			if (this != &other)
			{
				X = other.X;
				Z = other.Z;
			}

			return *this;
		}

		PauliStringXZ& operator=(PauliStringXZ&& other) noexcept
		{
			if (this != &other)
			{
				X.swap(other.X);
				Z.swap(other.Z);
			}

			return *this;
		}

		void Resize(size_t nQubits)
		{
			X.resize(nQubits, false);
			Z.resize(nQubits, false);
		}

		virtual void Clear()
		{
			std::fill(X.begin(), X.end(), false);
			std::fill(Z.begin(), Z.end(), false);
		}

		size_t GetNrQubits() const
		{
			return X.size();
		}

		inline static bool XOR(bool a, bool b)
		{
			return a != b;
		}

		virtual double ExpectationValue() const
		{
			for (size_t i = 0; i < X.size(); ++i)
				if (X[i]) // it's either X or Y
					return 0.0;

			return 1.0;
		}

		inline void ApplyH(size_t qubit)
		{
			// how does this work?
			// looking at it might not reveal immediately what it does
			// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

			// if we have both X and Z, then it's a Y (with some global phase, given by the sign, Y = iXZ)
			// a Y is transformed to a -Y, so a sign change is needed

			// swap X and Z
			const bool t = X[qubit];
			X[qubit] = Z[qubit];
			Z[qubit] = t;
		}

		inline void ApplyS(size_t qubit)
		{
			Z[qubit] = XOR(Z[qubit], X[qubit]);
		}

		inline void ApplyX(size_t qubit)
		{
			// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
		}

		inline void ApplyY(size_t qubit)
		{
			// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
			// if both X and Z are present, the sign is flipped twice, so it remains unchanged

			// can be done with ifs, can be done with XORs
		}

		inline void ApplyZ(size_t qubit)
		{
			// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as ZXZ^t = -X
		}

		inline void ApplyCX(size_t target, size_t control)
		{
			X[target] = XOR(X[target], X[control]);
			Z[control] = XOR(Z[control], Z[target]);
		}

		std::vector<bool> X;
		std::vector<bool> Z;
	};

	class PauliStringXZWithSign : public PauliStringXZ {
	public:
		PauliStringXZWithSign() : PhaseSign(false) {}

		PauliStringXZWithSign(size_t nQubits) : PauliStringXZ(nQubits), PhaseSign(false) {}

		PauliStringXZWithSign(const PauliStringXZWithSign& other) : PauliStringXZ(other), PhaseSign(other.PhaseSign) {}

		PauliStringXZWithSign(PauliStringXZWithSign&& other) noexcept : PauliStringXZ(other), PhaseSign(other.PhaseSign) {}

		PauliStringXZWithSign& operator=(const PauliStringXZWithSign& other)
		{
			if (this != &other)
			{
				PauliStringXZ::operator=(other);
				PhaseSign = other.PhaseSign;
			}

			return *this;
		}

		PauliStringXZWithSign& operator=(PauliStringXZWithSign&& other) noexcept
		{
			if (this != &other)
			{
				PhaseSign = other.PhaseSign;
				PauliStringXZ::operator=(std::move(other));
			}

			return *this;
		}

		void Clear() override
		{
			PauliStringXZ::Clear();
			PhaseSign = false;
		}

		double ExpectationValue() const override
		{
			const double expval = PauliStringXZ::ExpectationValue();

			return PhaseSign ? -expval : expval;
		}

		inline void ApplyH(size_t qubit)
		{
			// how does this work?
			// looking at it might not reveal immediately what it does
			// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

			// if we have both X and Z, then it's a Y (with some global phase, given by the sign, Y = iXZ)
			// a Y is transformed to a -Y, so a sign change is needed

			if (X[qubit] && Z[qubit])
				PhaseSign = !PhaseSign;

			// swap X and Z
			PauliStringXZ::ApplyH(qubit);
		}

		inline void ApplyS(size_t qubit)
		{
			if (X[qubit] && Z[qubit])
				PhaseSign = !PhaseSign;

			PauliStringXZ::ApplyS(qubit);
		}

		inline void ApplyX(size_t qubit)
		{
			// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
			if (Z[qubit])
				PhaseSign = !PhaseSign;
		}

		inline void ApplyY(size_t qubit)
		{
			// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
			// if both X and Z are present, the sign is flipped twice, so it remains unchanged

			// can be done with ifs, can be done with XORs
			if (XOR(Z[qubit], X[qubit]))
				PhaseSign = !PhaseSign;
		}

		inline void ApplyZ(size_t qubit)
		{
			// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as ZXZ^t = -X
			if (X[qubit])
				PhaseSign = !PhaseSign;
		}

		inline void ApplyCX(size_t target, size_t control)
		{
			if (X[control] && Z[target] && XOR(X[target], !Z[control]))
				PhaseSign = !PhaseSign;

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
			ApplyS(qubit2);
			ApplyH(qubit1);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
			ApplyH(qubit2);
		}

		inline void ApplyISwapDag(size_t qubit1, size_t qubit2)
		{
			ApplyH(qubit2);
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyH(qubit1);
			ApplySdg(qubit2);
			ApplySdg(qubit1);
		}

		bool PhaseSign; // true means negative phase
	};

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
			ApplyS(qubit2);
			ApplyH(qubit1);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
			ApplyH(qubit2);
		}

		inline void ApplyISwapDag(size_t qubit1, size_t qubit2)
		{
			ApplyH(qubit2);
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyH(qubit1);
			ApplySdg(qubit2);
			ApplySdg(qubit1);
		}
		
		double Coefficient;
	};

} // namespace QC