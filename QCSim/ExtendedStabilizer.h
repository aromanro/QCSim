#pragma once

#include "Frame.h"

namespace QC {

	// for now it will be a single frame implementation, might be extended to multiple frames later, if needed
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

		void ApplyH(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (auto& frame : frames)
			{
				for (size_t l = 0; l < nrQubits; ++l)
				{
					if (frame.stabilizers[l].X[qubit] && frame.stabilizers[l].Z[qubit])
						for (size_t k = 0; k < frame.GetFrameSize(); ++k)
							frame.signs[k][l] = !frame.signs[k][l];

					frame.stabilizers[l].ApplyH(qubit);
				}
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplyS(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (auto& frame : frames)
			{
				for (size_t l = 0; l < nrQubits; ++l)
				{
					if (frame.stabilizers[l].X[qubit] && frame.stabilizers[l].Z[qubit])
						for (size_t k = 0; k < frame.GetFrameSize(); ++k)
							frame.signs[k][l] = !frame.signs[k][l];

					frame.stabilizers[l].ApplyS(qubit);
				}
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplyX(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (auto& frame : frames)
			{
				for (size_t l = 0; l < nrQubits; ++l)
				{
					if (frame.stabilizers[l].Z[qubit])
						for (size_t k = 0; k < frame.GetFrameSize(); ++k)
							frame.signs[k][l] = !frame.signs[k][l];
				}
			}
		}

		void ApplyY(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (auto& frame : frames)
			{
				for (size_t l = 0; l < nrQubits; ++l)
				{
					if (frame.stabilizers[l].X[qubit] != frame.stabilizers[l].Z[qubit])
						for (size_t k = 0; k < frame.GetFrameSize(); ++k)
							frame.signs[k][l] = !frame.signs[k][l];
				}
			}
		}

		void ApplyZ(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (auto& frame : frames)
			{
				for (size_t l = 0; l < nrQubits; ++l)
				{
					if (frame.stabilizers[l].X[qubit])
						for (size_t k = 0; k < frame.GetFrameSize(); ++k)
							frame.signs[k][l] = !frame.signs[k][l];
				}
			}
		}

		void ApplySdg(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
		}

		void ApplyK(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		void ApplySx(size_t qubit)
		{
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyZ(qubit);
			ApplyS(qubit);
		}

		void ApplySxDag(size_t qubit)
		{
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		void ApplyCY(size_t target, size_t control)
		{
			ApplyZ(target);
			ApplyS(target);
			ApplyCX(target, control);
			ApplyS(target);
		}

		void ApplyCZ(size_t target, size_t control)
		{
			ApplyH(target);
			ApplyCX(target, control);
			ApplyH(target);
		}

		void ApplySwap(size_t qubit1, size_t qubit2)
		{
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
		}

		void ApplyISwap(size_t qubit1, size_t qubit2)
		{
			ApplyS(qubit1);
			ApplyH(qubit1);
			ApplyS(qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
			ApplyH(qubit2);
		}

		void ApplyISwapDag(size_t qubit1, size_t qubit2)
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

		void ApplyCX(size_t target, size_t control)
		{
			const size_t nrQubits = GetNrQubits();
			for (auto& frame : frames)
			{
				for (size_t l = 0; l < nrQubits; ++l)
				{
					if (frame.stabilizers[l].X[control] && frame.stabilizers[l].Z[target]
						&& PauliStringXZ::XOR(frame.stabilizers[l].X[target], !frame.stabilizers[l].Z[control]))
						for (size_t k = 0; k < frame.GetFrameSize(); ++k)
							frame.signs[k][l] = !frame.signs[k][l];

					frame.stabilizers[l].ApplyCX(target, control);
				}
				frame.ReduceToRowEchelonFormForColumns(target, control);
			}
		}

		bool Measure(size_t qubit)
		{
			auto& frame = frames.front();
			const size_t nrQubits = GetNrQubits();

			size_t p = nrQubits;
			for (size_t s = 0; s < nrQubits; ++s)
				if (frame.stabilizers[s].X[qubit])
				{
					p = s;
					break;
				}

			if (p < nrQubits)
			{
				for (size_t s = 0; s < nrQubits; ++s)
					if (s != p && frame.stabilizers[s].X[qubit])
						frame.MultiplyGeneratorInto(p, s);

				frame.stabilizers[p].Clear();
				frame.stabilizers[p].Z[qubit] = true;

				const bool outcome = rnd(gen);

				for (size_t k = 0; k < frame.GetFrameSize(); ++k)
					frame.signs[k][p] = outcome;

				frame.ReduceToRowEchelonForm();

				return outcome;
			}

			return frame.GetDeterministicOutcome(qubit, 0);
		}

		double GetQubitProbability(size_t qubit)
		{
			const auto& frame = frames.front();
			const size_t nrQubits = GetNrQubits();

			for (size_t s = 0; s < nrQubits; ++s)
				if (frame.stabilizers[s].X[qubit])
					return 0.5;

			double prob = 0.0;
			for (size_t k = 0; k < frame.GetFrameSize(); ++k)
				if (frame.GetDeterministicOutcome(qubit, k))
					prob += std::norm(frame.amplitudes[k]);
			
			return prob;
		}

		double ExpectationValue(const std::string& pauliString) const
		{
			if (pauliString.empty()) return 1.0;

			const auto& frame = frames.front();
			const size_t nrQubits = GetNrQubits();

			PauliStringXZ op(nrQubits);
			std::vector<size_t> pos;
			size_t phase = 0;
			SetPauliString(pauliString, op, pos, phase);
			if (pos.empty()) return 1.0;

			if (frame.CheckStabilizersAnticommutation(op, pos))
				return 0.0;

			double result = 0.0;
			for (size_t k = 0; k < frame.GetFrameSize(); ++k)
				result += std::norm(frame.amplitudes[k]) * frame.GetPauliExpectation(op, k);

			return result;
		}

		void SaveState()
		{
			savedFrames = frames;
		}

		void RestoreState()
		{
			frames = savedFrames;
		}

		std::vector<Frame> frames;

	private:
		static void SetPauliString(const std::string& pauliString, PauliStringXZ& g, std::vector<size_t>& pos, size_t& phase)
		{
			pos.reserve(pauliString.size());
			for (size_t i = 0; i < pauliString.size(); ++i)
			{
				const char c = toupper(pauliString[i]);
				switch (c)
				{
				case 'I':
					break;
				case 'X':
					g.X[i] = true;
					pos.push_back(i);
					break;
				case 'Y':
					g.X[i] = true;
					g.Z[i] = true;
					++phase;
					pos.push_back(i);
					break;
				case 'Z':
					g.Z[i] = true;
					pos.push_back(i);
					break;
				default:
					throw std::runtime_error("Invalid operator in the Pauli string");
				}
			}
		}

		std::vector<Frame> savedFrames;

		std::mt19937 gen;
		std::bernoulli_distribution rnd;
	};

}

