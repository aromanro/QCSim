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

		double GetPauliExpectation(const PauliStringXZ& op, size_t componentIndex) const
		{
			const size_t nrQubits = GetNrQubits();

			PauliStringXZ target(op);
			int phase = 0;

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

				if (signs[componentIndex][s])
					phase += 2;
			}

			return phase % 4 == 0 ? 1.0 : -1.0;
		}

		// this is not the most efficient way to do it
		// it can be done considering that only one or two columns are modified by gates
		// it's here for completion and tests
		void ReduceToRowEchelonForm()
		{
			const size_t nrQubits = GetNrQubits();

			size_t currentRow = 0;
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
			std::swap(stabilizers[i], stabilizers[j]);
			for (size_t k = 0; k < GetFrameSize(); ++k)
			{
				const bool t = signs[k][i];
				signs[k][i] = signs[k][j];
				signs[k][j] = t;
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
