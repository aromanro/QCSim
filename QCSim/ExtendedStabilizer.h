#pragma once

#include "Frame.h"

namespace QC {

	// for now it will be a single frame implementation, might be extended to multiple frames later, if needed
	class ExtendedStabilizer {
	public:
		ExtendedStabilizer() = delete;

		explicit ExtendedStabilizer(size_t nrQubits)
			: gen(std::random_device{}()), rnd(0.5), dist(0.0, 1.0)
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
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyHNoReduction(qubit);
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplyS(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplySNoReduction(qubit);
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplyX(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
				frame.ApplyX(qubit);
		}

		void ApplyY(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
				frame.ApplyY(qubit);
		}

		void ApplyZ(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
				frame.ApplyZ(qubit);
		}

		void ApplySdg(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyZ(qubit);
				frame.ApplySNoReduction(qubit);
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplyK(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyZ(qubit);
				frame.ApplySNoReduction(qubit);
				frame.ApplyHNoReduction(qubit);
				frame.ApplySNoReduction(qubit);
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplySx(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyZ(qubit);
				frame.ApplySNoReduction(qubit);
				frame.ApplyHNoReduction(qubit);
				frame.ApplyZ(qubit);
				frame.ApplySNoReduction(qubit);
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplySxDag(size_t qubit)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplySNoReduction(qubit);
				frame.ApplyHNoReduction(qubit);
				frame.ApplySNoReduction(qubit);
				frame.ReduceToRowEchelonFormForColumn(qubit);
			}
		}

		void ApplyCY(size_t target, size_t control)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyZ(target);
				frame.ApplySNoReduction(target);
				frame.ApplyCXNoReduction(target, control);
				frame.ApplySNoReduction(target);
				frame.ReduceToRowEchelonFormForColumns(target, control);
			}
		}

		void ApplyCZ(size_t target, size_t control)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyCZNoReduction(target, control);
				frame.ReduceToRowEchelonFormForColumns(target, control);
			}
		}

		void ApplySwap(size_t qubit1, size_t qubit2)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyCXNoReduction(qubit1, qubit2);
				frame.ApplyCXNoReduction(qubit2, qubit1);
				frame.ApplyCXNoReduction(qubit1, qubit2);
				frame.ReduceToRowEchelonFormForColumns(qubit1, qubit2);
			}
		}

		void ApplyISwap(size_t qubit1, size_t qubit2)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplySNoReduction(qubit1);
				frame.ApplyHNoReduction(qubit1);
				frame.ApplySNoReduction(qubit2);
				frame.ApplyCXNoReduction(qubit2, qubit1);
				frame.ApplyCXNoReduction(qubit1, qubit2);
				frame.ApplyHNoReduction(qubit2);
				frame.ReduceToRowEchelonFormForColumns(qubit1, qubit2);
			}
		}

		void ApplyISwapDag(size_t qubit1, size_t qubit2)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyHNoReduction(qubit2);
				frame.ApplyCXNoReduction(qubit1, qubit2);
				frame.ApplyCXNoReduction(qubit2, qubit1);
				frame.ApplyZ(qubit2);
				frame.ApplySNoReduction(qubit2);
				frame.ApplyHNoReduction(qubit1);
				frame.ApplyZ(qubit1);
				frame.ApplySNoReduction(qubit1);
				frame.ReduceToRowEchelonFormForColumns(qubit1, qubit2);
			}
		}

		void ApplyCX(size_t target, size_t control)
		{
			// TODO: also need to update the amplitudes
			for (auto& frame : frames)
			{
				frame.ApplyCXNoReduction(target, control);
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
				// nondeterministic
				for (size_t s = 0; s < nrQubits; ++s)
					if (s != p && frame.stabilizers[s].X[qubit])
						frame.MultiplyGeneratorInto(p, s);

				frame.stabilizers[p].Clear();
				frame.stabilizers[p].Z[qubit] = true;

				const bool outcome = rnd(gen);
				if (outcome)
				{
					for (size_t k = 0; k < frame.GetFrameSize(); ++k)
						if (frame.signs[k][p])
							frame.amplitudes[k] = -frame.amplitudes[k];
				}

				for (size_t k = 0; k < frame.GetFrameSize(); ++k)
					frame.signs[k][p] = outcome;

				frame.ReduceToRowEchelonForm();

				return outcome;
			}

			// deterministic
			const double prob1 = frame.GetDeterministicQubitProbability(qubit);
			const bool outcome = dist(gen) < prob1;

			double normSq = 0.0;
			size_t writeIdx = 0;
			for (size_t k = 0; k < frame.GetFrameSize(); ++k)
			{
				if (frame.GetDeterministicOutcome(qubit, k) == outcome)
				{
					if (writeIdx != k)
					{
						frame.amplitudes[writeIdx] = frame.amplitudes[k];
						frame.signs[writeIdx] = std::move(frame.signs[k]);
					}
					normSq += std::norm(frame.amplitudes[writeIdx]);
					++writeIdx;
				}
			}

			frame.amplitudes.resize(writeIdx);
			frame.signs.resize(writeIdx);

			if (normSq > 0.0)
			{
				const double invNorm = 1.0 / std::sqrt(normSq);
				for (auto& a : frame.amplitudes)
					a *= invNorm;
			}

			return outcome;
		}

		double GetQubitProbability(size_t qubit) const
		{
			const auto& frame = frames.front();

			if (!frame.IsDeterministic(qubit))
				return 0.5;

			return frame.GetDeterministicQubitProbability(qubit);
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

			return frame.GetPauliExpectation(op);
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
		std::uniform_real_distribution<double> dist;
	};

}

