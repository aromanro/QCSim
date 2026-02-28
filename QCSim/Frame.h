#pragma once

#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <cassert>
#include <vector>

#include "PauliStringXZ.h"

namespace QC {

	enum class NormalizationGateType: unsigned char{
		H = 0,
		S,
		CX,
		CZ
	};

	struct NormalizationGate {
		NormalizationGateType type = NormalizationGateType::H;
		size_t qubit1 = 0;
		size_t qubit2 = 0;
	};

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

		void ApplyHNoReduction(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t l = 0; l < nrQubits; ++l)
			{
				if (stabilizers[l].X[qubit] && stabilizers[l].Z[qubit])
					for (size_t k = 0; k < GetFrameSize(); ++k)
						signs[k][l] = !signs[k][l];

				stabilizers[l].ApplyH(qubit);
			}
		}

		void ApplySNoReduction(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t l = 0; l < nrQubits; ++l)
			{
				if (stabilizers[l].X[qubit] && stabilizers[l].Z[qubit])
					for (size_t k = 0; k < GetFrameSize(); ++k)
						signs[k][l] = !signs[k][l];

				stabilizers[l].ApplyS(qubit);
			}
		}

		void ApplyX(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t l = 0; l < nrQubits; ++l)
			{
				if (stabilizers[l].Z[qubit])
					for (size_t k = 0; k < GetFrameSize(); ++k)
						signs[k][l] = !signs[k][l];
			}
		}

		void ApplyY(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t l = 0; l < nrQubits; ++l)
			{
				if (stabilizers[l].X[qubit] != stabilizers[l].Z[qubit])
					for (size_t k = 0; k < GetFrameSize(); ++k)
						signs[k][l] = !signs[k][l];
			}
		}

		void ApplyZ(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t l = 0; l < nrQubits; ++l)
			{
				if (stabilizers[l].X[qubit])
					for (size_t k = 0; k < GetFrameSize(); ++k)
						signs[k][l] = !signs[k][l];
			}
		}

		void ApplyCXNoReduction(size_t target, size_t control)
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t l = 0; l < nrQubits; ++l)
			{
				if (stabilizers[l].X[control] && stabilizers[l].Z[target]
					&& PauliStringXZ::XOR(stabilizers[l].X[target], !stabilizers[l].Z[control]))
					for (size_t k = 0; k < GetFrameSize(); ++k)
						signs[k][l] = !signs[k][l];
				stabilizers[l].ApplyCX(target, control);
			}
		}

		void ApplyCZNoReduction(size_t target, size_t control)
		{
			ApplyHNoReduction(target);
			ApplyCXNoReduction(target, control);
			ApplyHNoReduction(target);
		}

		static inline bool CommutesWithGenerator(const PauliStringXZ& g, const PauliStringXZ& generator, const std::vector<size_t>& pos)
		{
			bool commutes = true;
			for (size_t j = 0; j < pos.size(); ++j)
			{
				const size_t pauliOpQubit = pos[j];
				if (g.X[pauliOpQubit] && generator.Z[pauliOpQubit])
					commutes = !commutes;
				if (g.Z[pauliOpQubit] && generator.X[pauliOpQubit])
					commutes = !commutes;
			}

			return commutes;
		}

		bool CheckStabilizersAnticommutation(const PauliStringXZ& g, const std::vector<size_t>& pos) const
		{
			for (size_t i = 0; i < GetNrQubits(); ++i)
				if (!CommutesWithGenerator(g, stabilizers[i], pos))
					return true;

			return false;
		}

		void MultiplyGeneratorInto(size_t src, size_t dst)
		{
			const size_t nrQubits = GetNrQubits();
			int phaseExponent = 0;

			for (size_t q = 0; q < nrQubits; ++q)
			{
				phaseExponent += PauliStringXZWithSign::g(
					PauliStringXZWithSign::BoolToInt(stabilizers[src].X[q]),
					PauliStringXZWithSign::BoolToInt(stabilizers[src].Z[q]),
					PauliStringXZWithSign::BoolToInt(stabilizers[dst].X[q]),
					PauliStringXZWithSign::BoolToInt(stabilizers[dst].Z[q])
				);

				stabilizers[dst].X[q] = stabilizers[dst].X[q] != stabilizers[src].X[q];
				stabilizers[dst].Z[q] = stabilizers[dst].Z[q] != stabilizers[src].Z[q];
			}

			for (size_t k = 0; k < GetFrameSize(); ++k)
			{
				const int sgnCnt = phaseExponent + (signs[k][src] ? 2 : 0) + (signs[k][dst] ? 2 : 0);
				signs[k][dst] = (sgnCnt % 4) != 0;
			}
		}

		bool GetDeterministicOutcome(size_t qubit, size_t componentIndex) const
		{
			const size_t nrQubits = GetNrQubits();

			PauliStringXZ target(nrQubits);
			target.Z[qubit] = true;

			bool sign = false;

			for (size_t s = 0; s < nrQubits; ++s)
			{
				size_t pivotCol = nrQubits;
				for (size_t c = 0; c < nrQubits; ++c)
					if (target.X[c])
					{
						pivotCol = c;
						break;
					}
				if (pivotCol >= nrQubits)
				{
					for (size_t c = 0; c < nrQubits; ++c)
						if (target.Z[c])
						{
							pivotCol = c;
							break;
						}
				}

				if (pivotCol >= nrQubits || (target.X[pivotCol] && !stabilizers[s].X[pivotCol]) || (!target.X[pivotCol] && target.Z[pivotCol] && !stabilizers[s].Z[pivotCol]))
					continue;

				int phaseExponent = 0;
				for (size_t q = 0; q < nrQubits; ++q)
				{
					phaseExponent += PauliStringXZWithSign::g(
						PauliStringXZWithSign::BoolToInt(stabilizers[s].X[q]),
						PauliStringXZWithSign::BoolToInt(stabilizers[s].Z[q]),
						PauliStringXZWithSign::BoolToInt(target.X[q]),
						PauliStringXZWithSign::BoolToInt(target.Z[q])
					);
					target.X[q] = target.X[q] != stabilizers[s].X[q];
					target.Z[q] = target.Z[q] != stabilizers[s].Z[q];
				}

				phaseExponent += signs[componentIndex][s] ? 2 : 0;

				sign = sign != (phaseExponent % 4 != 0);
			}

			return sign;
		}

		double GetPauliExpectation(const PauliStringXZ& op) const
		{
			const size_t nrQubits = GetNrQubits();

			PauliStringXZ target(op);
			int phase = 0;

			std::vector<bool> addPhase(nrQubits, false);
			for (size_t s = 0; s < nrQubits; ++s)
			{
				size_t pivotCol = nrQubits;
				for (size_t c = 0; c < nrQubits; ++c)
					if (target.X[c])
					{
						pivotCol = c;
						break;
					}
				if (pivotCol >= nrQubits)
				{
					for (size_t c = 0; c < nrQubits; ++c)
						if (target.Z[c])
						{
							pivotCol = c;
							break;
						}
				}

				if (pivotCol >= nrQubits || (target.X[pivotCol] && !stabilizers[s].X[pivotCol]) || (!target.X[pivotCol] && target.Z[pivotCol] && !stabilizers[s].Z[pivotCol]))
					continue;

				addPhase[s] = true;

				for (size_t q = 0; q < nrQubits; ++q)
				{
					phase += PauliStringXZWithSign::g(
						PauliStringXZWithSign::BoolToInt(stabilizers[s].X[q]),
						PauliStringXZWithSign::BoolToInt(stabilizers[s].Z[q]),
						PauliStringXZWithSign::BoolToInt(target.X[q]),
						PauliStringXZWithSign::BoolToInt(target.Z[q])
					);
					target.X[q] = target.X[q] != stabilizers[s].X[q];
					target.Z[q] = target.Z[q] != stabilizers[s].Z[q];
				}
			}

			int phaseSave = phase;
			double result = 0.0;
			for (size_t componentIndex = 0; componentIndex < GetFrameSize(); ++componentIndex)
			{
				for (size_t s = 0; s < nrQubits; ++s)
					if (addPhase[s] && signs[componentIndex][s])
						phase += 2;

				result += std::norm(amplitudes[componentIndex]) * (phase % 4 == 0 ? 1.0 : -1.0);
				phase = phaseSave;
			}

			return result;
		}

		bool IsDeterministic(size_t qubit) const
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t s = 0; s < nrQubits; ++s)
				if (stabilizers[s].X[qubit])
					return false;
			return true;
		}

		double GetQubitProbability(size_t qubit) const
		{
			const size_t nrQubits = GetNrQubits();

			for (size_t s = 0; s < nrQubits; ++s)
				if (stabilizers[s].X[qubit])
					return 0.5;

			double prob = 0.0;
			for (size_t k = 0; k < GetFrameSize(); ++k)
				if (GetDeterministicOutcome(qubit, k))
					prob += std::norm(amplitudes[k]);

			return prob;
		}

		double GetDeterministicQubitProbability(size_t qubit) const
		{
			const size_t nrQubits = GetNrQubits();

			double prob = 0.0;
			for (size_t k = 0; k < GetFrameSize(); ++k)
				if (GetDeterministicOutcome(qubit, k))
					prob += std::norm(amplitudes[k]);

			return prob;
		}

		// this is not the most efficient way to do it
		// it can be done considering that only one or two columns are modified by gates
		// it's here for completion and tests
		// see 'Efficient Inner-Product Algorithm for Stabilizer States' by Garcia, Markov and Cross for details
		void ReduceToRowEchelonForm()
		{
			const size_t nrQubits = GetNrQubits();

			size_t currentRow = 0;
			// X block
			for (size_t col = 0; col < nrQubits; ++col)
			{
				size_t pivotRow = nrQubits;
				for (size_t r = currentRow; r < nrQubits; ++r)
					if (stabilizers[r].X[col])
					{
						pivotRow = r;
						break;
					}
				if (pivotRow >= nrQubits) continue;
				if (pivotRow != currentRow)
					Swap(pivotRow, currentRow);
				for (size_t r = 0; r < nrQubits; ++r)
				{
					if (r == currentRow) continue;
					if (stabilizers[r].X[col])
						MultiplyGeneratorInto(currentRow, r);
				}
				++currentRow;
			}
			
			// Z block
			for (size_t col = 0; col < nrQubits && currentRow < nrQubits; ++col)
			{
				size_t pivotRow = nrQubits;
				for (size_t r = currentRow; r < nrQubits; ++r)
					if (stabilizers[r].Z[col])
					{
						pivotRow = r;
						break;
					}
				if (pivotRow >= nrQubits) continue;
				if (pivotRow != currentRow)
					Swap(pivotRow, currentRow);
				for (size_t r = 0; r < nrQubits; ++r)
				{
					if (r == currentRow) continue;
					if (stabilizers[r].Z[col])
						MultiplyGeneratorInto(currentRow, r);
				}
				++currentRow;
			}
		}

		void ReduceToRowEchelonFormForColumn(size_t affectedCol)
		{
			const size_t nrQubits = GetNrQubits();

			size_t currentRow = 0;
			for (size_t col = 0; col < affectedCol; ++col)
				if (stabilizers[currentRow].X[col])
					++currentRow;

			for (size_t col = affectedCol; col < nrQubits; ++col)
			{
				size_t pivotRow = nrQubits;
				for (size_t r = currentRow; r < nrQubits; ++r)
					if (stabilizers[r].X[col])
					{
						pivotRow = r;
						break;
					}
				if (pivotRow >= nrQubits) continue;
				if (pivotRow != currentRow)
					Swap(pivotRow, currentRow);
				for (size_t r = 0; r < nrQubits; ++r)
				{
					if (r == currentRow) continue;
					if (stabilizers[r].X[col])
						MultiplyGeneratorInto(currentRow, r);
				}
				++currentRow;
			}

			for (size_t col = 0; col < nrQubits && currentRow < nrQubits; ++col)
			{
				size_t pivotRow = nrQubits;
				for (size_t r = currentRow; r < nrQubits; ++r)
					if (stabilizers[r].Z[col])
					{
						pivotRow = r;
						break;
					}
				if (pivotRow >= nrQubits) continue;
				if (pivotRow != currentRow)
					Swap(pivotRow, currentRow);
				for (size_t r = 0; r < nrQubits; ++r)
				{
					if (r == currentRow) continue;
					if (stabilizers[r].Z[col])
						MultiplyGeneratorInto(currentRow, r);
				}
				++currentRow;
			}
		}

		void ReduceToRowEchelonFormForColumns(size_t col1, size_t col2)
		{
			ReduceToRowEchelonFormForColumn(std::min(col1, col2));
		}

		void Swap(size_t i, size_t j)
		{
			if (i == j) return;
			std::swap(stabilizers[i], stabilizers[j]);
			for (size_t k = 0; k < GetFrameSize(); ++k)
			{
				const bool t = signs[k][i];
				signs[k][i] = signs[k][j];
				signs[k][j] = t;
			}
		}

		void Cofactor(size_t qubit)
		{
			const size_t nrQubits = GetNrQubits();

			// check the qubit to see if it's deterministic or not
			size_t p = nrQubits;
			for (size_t s = 0; s < nrQubits; ++s)
				if (stabilizers[s].X[qubit])
				{
					p = s;
					break;
				}

			// deterministic, duplication not needed
			if (p >= nrQubits) return;

			// remove X from all the other stabilizers
			for (size_t s = 0; s < nrQubits; ++s)
				if (s != p && stabilizers[s].X[qubit])
					MultiplyGeneratorInto(p, s);

			// collapse the stabilizer to Z on the qubit
			stabilizers[p].Clear();
			stabilizers[p].Z[qubit] = true;

			// duplicate and update the signs and amplitudes
			const size_t oldSize = GetFrameSize();
			const double invSqrt2 = 1.0 / std::sqrt(2.0);

			const size_t newSize = 2 * oldSize;
			amplitudes.resize(newSize);
			signs.resize(newSize);

			for (size_t k = 0; k < oldSize; ++k)
			{
				amplitudes[k] *= invSqrt2;

				amplitudes[oldSize + k] = amplitudes[k];
				signs[oldSize + k] = signs[k];
				signs[oldSize + k][p] = !signs[k][p];
			}

			ReduceToRowEchelonForm();
		}

		// assumes the frame is already in row echelon form
		std::vector<NormalizationGate> BasisNormalize()
		{
			std::vector<NormalizationGate> gates;
			const size_t nrQubits = GetNrQubits();

			// apply block of Hadamard gates
			for (size_t j = 0; j < nrQubits; ++j)
			{
				size_t k = nrQubits;
				for (size_t r = j; r < nrQubits; ++r)
					if (stabilizers[r].X[j]) // X or Y
					{
						k = r;
						break;
					}
				if (k < nrQubits)
				{
					Swap(j, k);
					for (size_t r = j + 1; r < nrQubits; ++r)
						if (stabilizers[r].X[j])
							MultiplyGeneratorInto(j, r);
				}
				else {
					for (long long int r = nrQubits - 1; r >= j; --r)
						if (stabilizers[r].Z[j]) // Z
						{
							k = r;
							break;
						}
					if (k < nrQubits)
					{
						Swap(j, k);

						for (size_t r = j + 1; r < nrQubits; ++r)
							if (stabilizers[r].Z[j])
								MultiplyGeneratorInto(j, r);

						for (size_t c = j + 1; c < nrQubits; ++c)
						{
							if (stabilizers[j].X[c] || stabilizers[j].Z[c]) // has X, Y or Z
							{
								ApplyHNoReduction(j);
								gates.push_back({ NormalizationGateType::H, j });

								break;
							}
						}
					}
				}
			}

			// check: below the diagonal are only I and Z
			/*
			for (size_t j = 0; j < nrQubits; ++j)
				for (size_t k = 0; k < j; ++k)
					if (stabilizers[j].X[k]) // X or Y
					{
						std::cerr << "Below diagonal there aren't only I and Z, can't be normalized to basis form" << std::endl;
						exit(1);
					}
			*/

			// above the diagonal
			// apply block of CX gates
			for (size_t j = 0; j < nrQubits; ++j)
				for (size_t k = j + 1; k < nrQubits; ++k)
					if (stabilizers[j].X[k]) // X or Y
					{
						ApplyCXNoReduction(k, j);
						gates.push_back({ NormalizationGateType::CX, k, j });
					}

			// above the diagonal
			// apply block of CZ gates
			for (size_t j = 0; j < nrQubits; ++j)
				for (size_t k = j + 1; k < nrQubits; ++k)
					if (stabilizers[j].Z[k]) // Z
					{
						ApplyCZNoReduction(k, j);
						gates.push_back({ NormalizationGateType::CZ, k, j });
					}

			// diagonal
			// apply block of S gates
			for (size_t j = 0; j < nrQubits; ++j)
				if (stabilizers[j].X[j] && stabilizers[j].Z[j]) // Y
				{
					ApplySNoReduction(j);
					gates.push_back({ NormalizationGateType::S, j });
				}

			// diagonal
			// apply block of H gates
			for (size_t j = 0; j < nrQubits; ++j)
				if (stabilizers[j].X[j]) // X
				{
					ApplyHNoReduction(j);
					gates.push_back({ NormalizationGateType::H, j });
				}

			// eliminate trailing Zs below the diagonal to ensure basis form
			for (size_t j = 0; j < nrQubits; ++j)
				for (size_t k = j + 1; k < nrQubits; ++k)
					if (stabilizers[k].Z[j]) // Z
						MultiplyGeneratorInto(j, k);

			return gates;
		}

		void ApplyCircuit(const std::vector<NormalizationGate>& gates)
		{
			for (const auto& gate : gates)
			{
				switch (gate.type)
				{
				case NormalizationGateType::H:
					ApplyHNoReduction(gate.qubit1);
					break;
				case NormalizationGateType::S:
					ApplySNoReduction(gate.qubit1);
					break;
				case NormalizationGateType::CX:
					ApplyCXNoReduction(gate.qubit1, gate.qubit2);
					break;
				case NormalizationGateType::CZ:
					ApplyCZNoReduction(gate.qubit1, gate.qubit2);
					break;
				default:
					std::cerr << "Unknown gate type in ApplyCircuit" << std::endl;
					break;
				}
			}
		}

		void Print() const
		{
			const size_t nrQubits = GetNrQubits();
			for (size_t s = 0; s < nrQubits; ++s)
			{
				std::string line;
				for (size_t q = 0; q < nrQubits; ++q)
				{
					if (stabilizers[s].X[q] && stabilizers[s].Z[q])
						line += "Y";
					else if (stabilizers[s].X[q])
						line += "X";
					else if (stabilizers[s].Z[q])
						line += "Z";
					else
						line += "I";
				}
				std::cout << line << std::endl;
			}
		}

		std::vector<std::complex<double>> amplitudes;
		std::vector<Signs> signs;
		std::vector<PauliStringXZ> stabilizers;
	};

}
