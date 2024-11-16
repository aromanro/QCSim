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

			std::vector<bool> X;
			std::vector<bool> Z;
			bool PhaseSign; // true means negative phase
		};

	}
}
