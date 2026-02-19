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

		PauliStringXZ(PauliStringXZ&& other) noexcept
		{
			X.swap(other.X);
			Z.swap(other.Z);
		}

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

		bool operator==(const PauliStringXZ& other) const
		{
			return (X == other.X) && (Z == other.Z);
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

		size_t PauliWeight() const
		{
			size_t weight = 0;
			for (size_t i = 0; i < X.size(); ++i)
				if (X[i] || Z[i])
					++weight;

			return weight;
		}

		bool IsX(size_t qubit) const
		{
			return X[qubit] && !Z[qubit];
		}

		bool IsY(size_t qubit) const
		{
			return X[qubit] && Z[qubit];
		}

		bool IsZ(size_t qubit) const
		{
			return !X[qubit] && Z[qubit];
		}

		bool IsI(size_t qubit) const
		{
			return !X[qubit] && !Z[qubit];
		}

		void SetI(size_t qubit)
		{
			X[qubit] = false;
			Z[qubit] = false;
		}

		void SetX(size_t qubit)
		{
			X[qubit] = true;
			Z[qubit] = false;
		}

		void SetY(size_t qubit)
		{
			X[qubit] = true;
			Z[qubit] = true;
		}

		void SetZ(size_t qubit)
		{
			X[qubit] = false;
			Z[qubit] = true;
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

		inline virtual void ApplyH(size_t qubit)
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

		inline virtual void ApplyS(size_t qubit)
		{
			Z[qubit] = XOR(Z[qubit], X[qubit]);
		}

		inline virtual void ApplyX(size_t /*qubit*/)
		{
			// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
		}

		inline virtual void ApplyY(size_t /*qubit*/)
		{
			// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
			// if both X and Z are present, the sign is flipped twice, so it remains unchanged

			// can be done with ifs, can be done with XORs
		}

		inline virtual void ApplyZ(size_t /*qubit*/)
		{
			// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as ZXZ^t = -X
		}

		inline virtual void ApplyCX(size_t target, size_t control)
		{
			X[target] = XOR(X[target], X[control]);
			Z[control] = XOR(Z[control], Z[target]);
		}



		inline void ApplyK(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		inline void ApplySdag(size_t qubit)
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

		virtual std::string ToString() const
		{
			std::string result;
			result.resize(X.size());
			for (size_t i = 0; i < X.size(); ++i)
			{
				if (X[i] && Z[i])
					result[i] = 'Y';
				else if (X[i])
					result[i] = 'X';
				else if (Z[i])
					result[i] = 'Z';
				else
					result[i] = 'I';
			}
			return result;
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

		inline void ApplyH(size_t qubit) override
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

		inline void ApplyS(size_t qubit) override
		{
			if (X[qubit] && Z[qubit])
				PhaseSign = !PhaseSign;

			PauliStringXZ::ApplyS(qubit);
		}

		inline void ApplyX(size_t qubit) override
		{
			// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
			if (Z[qubit])
				PhaseSign = !PhaseSign;
		}

		inline void ApplyY(size_t qubit) override
		{
			// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
			// if both X and Z are present, the sign is flipped twice, so it remains unchanged

			// can be done with ifs, can be done with XORs
			if (XOR(Z[qubit], X[qubit]))
				PhaseSign = !PhaseSign;
		}

		inline void ApplyZ(size_t qubit) override
		{
			// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as ZXZ^t = -X
			if (X[qubit])
				PhaseSign = !PhaseSign;
		}

		inline void ApplyCX(size_t target, size_t control) override
		{
			if (X[control] && Z[target] && XOR(X[target], !Z[control]))
				PhaseSign = !PhaseSign;

			PauliStringXZ::ApplyCX(target, control);
		}

		// the following are for the extended stabilizer simulator
		inline void ApplyXleft(size_t qubit)
		{
			X[qubit] = !X[qubit];
		}

		inline void ApplyYleft(size_t qubit)
		{
			if (X[qubit] && !Z[qubit])
				PhaseSign = !PhaseSign;

			X[qubit] = !X[qubit];
			Z[qubit] = !Z[qubit];
		}

		inline void ApplyZleft(size_t qubit)
		{
			if (X[qubit] && Z[qubit])
				PhaseSign = !PhaseSign;
			Z[qubit] = !Z[qubit];
		}

		std::string ToString() const override
		{
			std::string result = PauliStringXZ::ToString();
			if (PhaseSign)
				result = "-" + result;
			else
				result = "+" + result;

			return result;
		}

		// returns the exponent of the i that multiplies the product of the two corresponding Pauli matrices
		// 0, 1 or -1
		// for example for x1 = 1, z1 = 0, x2 = 0, z2 = 1
		// we have X * Z = -i Y so the result should be -1
		// for x1 = 1, z1 = 0, x2 = 1, z2 = 0
		// we have X * X = I so the result should be 0
		static inline int g(int x1, int z1, int x2, int z2)
		{
			if (0 == x1 && 0 == z1) return 0; // I for the 1st generator, 0 exponent no matter what the second generator is
			else if (1 == x1)
			{
				if (1 == z1) return z2 - x2;

				return z2 * (2 * x2 - 1);
			}

			return x2 * (1 - 2 * z2);
		}

		// multiplies the two generators and stores the result in the current one
		inline void Multiply(const PauliStringXZWithSign& j, bool enableMultithreading)
		{
			const size_t nrQubits = X.size();
			// phase sign is negative when 'PhaseSign' is true
			// 2 because i^2 = -1, 0 because i^0 = 1
			long long int m = (PhaseSign ? 2 : 0) + (j.PhaseSign ? 2 : 0); // this gets the sign from the signs of the two generators
			
			// we still need to add the contribution from the Pauli strings:

			if (!enableMultithreading || nrQubits < 1024)
			{
				for (size_t q = 0; q < nrQubits; ++q)
				{
					const int x1 = BoolToInt(j.X[q]);
					const int z1 = BoolToInt(j.Z[q]);
					const int x2 = BoolToInt(X[q]);
					const int z2 = BoolToInt(Z[q]);

					// add up all the exponents of i that contribute to the sign of the product
					m += g(x1, z1, x2, z2);

					// X * X = I, Z * Z = I, so the value is set when there is only one of them
					X[q] = (x1 ^ x2) == 1;
					Z[q] = (z1 ^ z2) == 1;
				}
			}
			else
			{
				//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();
				long long int mloc = 0;

#pragma omp parallel for reduction(+:mloc) 
				//num_threads(processor_count) schedule(static, 256)
				for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
				{
					const int x1 = BoolToInt(j.X[q]);
					const int z1 = BoolToInt(j.Z[q]);
					const int x2 = BoolToInt(X[q]);
					const int z2 = BoolToInt(Z[q]);

					// add up all the exponents of i that contribute to the sign of the product
					mloc += g(x1, z1, x2, z2);

					// X * X = I, Z * Z = I, so the value is set when there is only one of them
					X[q] = (x1 ^ x2) == 1;
					Z[q] = (z1 ^ z2) == 1;
				}

				m += mloc;
			}

			// the mod 4 that appears here is because the values for the powers of i keep repeating
			const int mod = m % 4;

			assert(mod == 0 || mod == 2 || mod == -2);

			PhaseSign = mod != 0;
		}

		static inline int BoolToInt(bool b)
		{
			return b ? 1 : 0;
		}

		bool PhaseSign; // true means negative phase
	};

} // namespace QC