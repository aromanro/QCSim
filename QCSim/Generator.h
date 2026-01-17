#pragma once

#include <vector>
#include <algorithm>

namespace QC {
	namespace Clifford {

		class Generator {
		public:
			Generator() : PhaseSign(false) {}

			Generator(size_t nQubits) : X(nQubits, false), Z(nQubits, false), PhaseSign(false) {}

			Generator(const Generator& other) : X(other.X), Z(other.Z), PhaseSign(other.PhaseSign) {}

			Generator(Generator&& other) noexcept : X(std::move(other.X)), Z(std::move(other.Z)), PhaseSign(other.PhaseSign) {}

			Generator& operator=(const Generator& other)
			{
				if (this != &other)
				{
					X = other.X;
					Z = other.Z;
					PhaseSign = other.PhaseSign;
				}

				return *this;
			}

			Generator& operator=(Generator&& other) noexcept
			{
				if (this != &other)
				{
					X.swap(other.X);
					Z.swap(other.Z);
					PhaseSign = other.PhaseSign;
				}

				return *this;
			}

			void Resize(size_t nQubits)
			{
				X.resize(nQubits, false);
				Z.resize(nQubits, false);
			}

			void Clear()
			{
				std::fill(X.begin(), X.end(), false);
				std::fill(Z.begin(), Z.end(), false);
				PhaseSign = false;
			}

			size_t GetNrQubits() const
			{
				return X.size();
			}

			inline static bool XOR(bool a, bool b)
			{
				return a != b;
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
				const bool t = X[qubit];
				X[qubit] = Z[qubit];
				Z[qubit] = t;
			}

			inline void ApplyS(size_t qubit)
			{
				if (X[qubit] && Z[qubit])
					PhaseSign = !PhaseSign;

				Z[qubit] = XOR(Z[qubit], X[qubit]);
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

				X[target] = XOR(X[target], X[control]);
				Z[control] = XOR(Z[control], Z[target]);
			}

			std::vector<bool> X;
			std::vector<bool> Z;
			bool PhaseSign; // true means negative phase
		};

	}
}
