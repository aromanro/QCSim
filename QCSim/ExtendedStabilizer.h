#pragma once
		
#include <complex>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "PauliStringXZ.h"

namespace QC {

	class Frame {
	public:
		using Signs = std::vector<bool>;

		Frame(size_t nrQubits) : amplitudes(1, std::complex<double>(1., 0.)), signs(1, std::vector<bool>(nrQubits, false))
		{
			stabilizers.resize(nrQubits);
			for (size_t s = 0; s < nrQubits; ++s)
				stabilizers[s] = PauliStringXZ(nrQubits);
		}

		size_t GetNrQubits() const
		{
			if (signs.empty()) return 0;

			return signs.front().size();
		}

		void SetZDiagonal()
		{
			for (size_t s = 0; s < GetNrQubits(); ++s)
				stabilizers[s].Z[s] = true;
		}

		size_t GetFrameSize() const
		{
			return amplitudes.size();
		}

		std::vector<std::complex<double>> amplitudes;
		std::vector<Signs> signs;
		std::vector<PauliStringXZ> stabilizers;
	};

	class ExtendedStabilizer {
	public:
		ExtendedStabilizer() = delete;

		explicit ExtendedStabilizer(size_t nrQubits)
			: gen(std::random_device{}()), rnd(0.5)
		{
			frames.emplace_back(nrQubits);
			frames.back().SetZDiagonal();
		}

		size_t GetNrQubits() const
		{
			if (frames.empty()) return 0;
			return frames.front().GetNrQubits();
		}

		void Reset(size_t nrQubits)
		{
			frames.clear();
			frames.emplace_back(nrQubits);
			frames.back().SetZDiagonal();
		}

		void ApplyH(size_t /*qubit*/)
		{
			throw std::logic_error("ExtendedStabilizer::ApplyH not implemented yet");
		}

		void ApplyS(size_t /*qubit*/)
		{
			throw std::logic_error("ExtendedStabilizer::ApplyS not implemented yet");
		}

		void ApplyCX(size_t /*target*/, size_t /*control*/)
		{
			throw std::logic_error("ExtendedStabilizer::ApplyCX not implemented yet");
		}

		bool Measure(size_t /*qubit*/)
		{
			throw std::logic_error("ExtendedStabilizer::Measure not implemented yet");
		}

		double GetQubitProbability(size_t /*qubit*/)
		{
			throw std::logic_error("ExtendedStabilizer::GetQubitProbability not implemented yet");
		}

		std::vector<Frame> frames;

	private:
		std::mt19937 gen;
		std::bernoulli_distribution rnd;
	};

}

