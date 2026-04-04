#pragma once

#include <vector>

#include "QuantumGate.h"

namespace QC {
	namespace PathIntegral {
		class PathIntegralSimulator {
		public:
			// the start state is |0>
			static std::complex<double> Propagate(const std::vector<QC::Gates::AppliedGate<>>& circuit, const std::vector<bool>& endState)
			{
				const std::vector<bool> startState(endState.size(), false);

				return Propagate(circuit, startState, endState);
			}

			static std::complex<double> Propagate(const std::vector<QC::Gates::AppliedGate<>>& circuit, const std::vector<bool>& startState, const std::vector<bool>& endState)
			{
				assert(startState.size() == endState.size());

				const size_t possibleQubitsChanges = CountPossibleQubitsChanges(circuit);

				return Propagate(circuit, startState, endState, 0, possibleQubitsChanges);
			}

		private:
			static std::complex<double> Propagate(const std::vector<QC::Gates::AppliedGate<>>& circuit, const std::vector<bool>& currentState, const std::vector<bool>& endState, size_t gateIndex, size_t possibleQubitsChanges)
			{
				if (gateIndex >= circuit.size())
					return currentState == endState ? std::complex<double>(1., 0.) : std::complex<double>(0., 0.);

				const double epsilon = 1e-15;
	
				const auto& gate = circuit[gateIndex];
				const auto& U = gate.getRawOperatorMatrix();
				const size_t gateQubits = gate.getQubitsNumber();

				assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

				possibleQubitsChanges -= gateQubits;

				std::complex<double> amplitude = 0.;
				if (gateQubits == 1)
				{
					const size_t qubit = gate.getQubit1();
					assert(qubit < currentState.size());

					const Eigen::Index col = (currentState[qubit] ? 1 : 0);

					for (Eigen::Index row = 0; row < 2; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::abs(std::real(val)) > epsilon || std::abs(std::imag(val)) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit] = row == 1;

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							amplitude += val * Propagate(circuit, nextState, endState, gateIndex + 1, possibleQubitsChanges);
						}
					}
				}
				else if (gateQubits == 2)
				{
					const size_t qubit1 = gate.getQubit1();
					const size_t qubit2 = gate.getQubit2();
					assert(qubit1 < currentState.size() && qubit2 < currentState.size());

					const Eigen::Index col = ((currentState[qubit2] ? 2 : 0) | (currentState[qubit1] ? 1 : 0));

					for (Eigen::Index row = 0; row < 4; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::abs(std::real(val)) > epsilon || std::abs(std::imag(val)) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							amplitude += val * Propagate(circuit, nextState, endState, gateIndex + 1, possibleQubitsChanges);
						}
					}
				}
				else // only up to three qubits gates are supported
				{
					const size_t qubit1 = gate.getQubit1();
					const size_t qubit2 = gate.getQubit2();
					const size_t qubit3 = gate.getQubit3();
					assert(qubit1 < currentState.size() && qubit2 < currentState.size() && qubit3 < currentState.size());

					const Eigen::Index col = ((currentState[qubit3] ? 4 : 0) | (currentState[qubit2] ? 2 : 0) | (currentState[qubit1] ? 1 : 0));

					for (Eigen::Index row = 0; row < 8; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::abs(std::real(val)) > epsilon || std::abs(std::imag(val)) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;
							nextState[qubit3] = (row & 4) == 4;

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							amplitude += val * Propagate(circuit, nextState, endState, gateIndex + 1, possibleQubitsChanges);
						}
					}
				}

				return amplitude;
			}

			static size_t CountPossibleQubitsChanges(const std::vector<QC::Gates::AppliedGate<>>& circuit)
			{
				size_t count = 0;

				for (const auto& gate : circuit)
					count += gate.getQubitsNumber();

				return count;
			}

			static size_t CountDifferentQubits(const std::vector<bool>& curState, const std::vector<bool>& endState)
			{
				assert(curState.size() == endState.size());
				size_t count = 0;
				for (size_t i = 0; i < curState.size(); ++i)
					if (curState[i] != endState[i])
						++count;
				return count;
			}
		};
	}
}