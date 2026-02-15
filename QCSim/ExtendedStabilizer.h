#pragma once

#include "PauliStringXZ.h"

namespace QC {

	class Frame {
	public:
		using Signs = std::vector<bool>;

		Frame(size_t nrQubits) : amplitude(1., 0.), signs(nrQubits, std::vector<bool>(nrQubits, false)) {}

		size_t GetNrQubits() const
		{
			return signs.size();
		}

		void SetZDiagonal()
		{
			for (size_t s = 0; s < GetNrQubits(); ++s)
				stabilizers[s].Z[s] = true;
		}

		std::complex<double> amplitude;
		std::vector<Signs> signs;
		std::vector<PauliStringXZ> stabilizers;
	};

	class ExtendedStabilizer {
	public:
		ExtendedStabilizer() = delete;

		ExtendedStabilizer(size_t nrQubits)
		{
			frames.emplace_back(nrQubits);
			frames.back().SetZDiagonal();
		}

		std::vector<Frame> frames;
	};

}

