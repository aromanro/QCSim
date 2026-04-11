#pragma once

#include <vector>

#include "QuantumGate.h"

namespace QC {
	namespace PathIntegral {
		class PathIntegralSimulator {
		public:
			void SetCircuit(const std::vector<QC::Gates::AppliedGate<>>& circuit)
			{
				intermediateAmplitudes.clear();
				circuitBack.clear();

				if (doublingsLimit > 0)
				{
					size_t doublings = 0;
					for (const auto& gate : circuit)
					{
						if (!gate.isBranching())
							continue;

						doublings += gate.getMaxBranching();
					}

					const size_t halfDoublings = doublings / 2;

					// iterate backwards 
					doublings = 0;
					size_t gatesAdded = 0;
					for (auto it = circuit.rbegin(); it != circuit.rend(); ++it)
					{
						const auto& gate = *it;
						
						if (!gate.isBranching())
						{
							circuitBack.emplace_back(gate.getRawOperatorMatrix().transpose(), gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
							++gatesAdded;
							continue;
						}

						doublings += gate.getMaxBranching();
						if (doublings > doublingsLimit || doublings > halfDoublings)
							break;

						circuitBack.emplace_back(gate.getRawOperatorMatrix().transpose(), gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
						++gatesAdded;
					}

					this->circuit = std::vector<QC::Gates::AppliedGate<>>(circuit.begin(), circuit.end() - gatesAdded);
				}
				else
					this->circuit = circuit;
			}

			// the start state is |0>
			std::complex<double> Propagate(const std::vector<bool>& endState)
			{
				const std::vector<bool> startState(endState.size(), false);

				return Propagate(startState, endState);
			}

			std::complex<double> Propagate(const std::vector<bool>& startState, const std::vector<bool>& endState)
			{
				assert(startState.size() == endState.size());

				const size_t possibleQubitsChanges = CountPossibleQubitsChanges(circuit);

				if (doublingsLimit > 0 && !circuitBack.empty())
				{
					intermediateAmplitudes.clear();

					PropagateBackward(endState, std::complex<double>(1., 0.), 0);

					return PropagateForward(startState, 0, possibleQubitsChanges);
				}

				return Propagate(startState, endState, 0, possibleQubitsChanges);
			}

			void SetTrimValue(double val)
			{
				epsilon = val;
			}

			double GetTrimValue() const
			{
				return epsilon;
			}

			void SetMaxDoublingsForBackwardPaths(size_t doublings)
			{
				doublingsLimit = doublings;
			}

			size_t GetMaxDoublingsForBackwardPaths() const
			{
				return doublingsLimit;
			}

		private:
			// TODO: this can be parallelized!
			std::complex<double> Propagate(const std::vector<bool>& currentState, const std::vector<bool>& endState, size_t gateIndex, size_t possibleQubitsChanges)
			{
				if (gateIndex >= circuit.size())
					return currentState == endState ? std::complex<double>(1., 0.) : std::complex<double>(0., 0.);
	
				const auto& gate = circuit[gateIndex];
				const auto& U = gate.getRawOperatorMatrix();
				const size_t gateQubits = gate.getQubitsNumber();

				assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

				possibleQubitsChanges -= gateQubits;

				std::complex<double> amplitude(0., 0.);

				if (gateQubits == 1)
				{
					const size_t qubit = gate.getQubit1();
					assert(qubit < currentState.size());

					const Eigen::Index col = (currentState[qubit] ? 1 : 0);

					for (Eigen::Index row = 0; row < 2; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit] = row == 1;

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
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
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
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
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;
							nextState[qubit3] = (row & 4) == 4;

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
						}
					}
				}

				return amplitude;
			}

			void PropagateBackward(const std::vector<bool>& currentState, std::complex<double> amplitude, size_t gateIndex)
			{
				if (gateIndex >= circuitBack.size())
				{
					intermediateAmplitudes[currentState] += amplitude;
					return;
				}

				const auto& gate = circuitBack[gateIndex];
				const auto& U = gate.getRawOperatorMatrix();
				const size_t gateQubits = gate.getQubitsNumber();

				assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

				if (gateQubits == 1)
				{
					const size_t qubit = gate.getQubit1();
					assert(qubit < currentState.size());

					const Eigen::Index col = (currentState[qubit] ? 1 : 0);

					for (Eigen::Index row = 0; row < 2; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit] = row == 1;

							PropagateBackward(nextState, val * amplitude, gateIndex + 1);
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
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;

							PropagateBackward(nextState, val * amplitude, gateIndex + 1);
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
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;
							nextState[qubit3] = (row & 4) == 4;

							PropagateBackward(nextState, val * amplitude, gateIndex + 1);
						}
					}
				}
			}

			std::complex<double> PropagateForward(const std::vector<bool>& currentState, size_t gateIndex, size_t possibleQubitsChanges)
			{
				if (gateIndex >= circuit.size())
				{
					if (intermediateAmplitudes.find(currentState) != intermediateAmplitudes.end())
						return intermediateAmplitudes[currentState];
					
					return std::complex<double>(0., 0.);
				}

				const auto& gate = circuit[gateIndex];
				const auto& U = gate.getRawOperatorMatrix();
				const size_t gateQubits = gate.getQubitsNumber();

				assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

				possibleQubitsChanges -= gateQubits;

				std::complex<double> amplitude(0., 0.);

				if (gateQubits == 1)
				{
					const size_t qubit = gate.getQubit1();
					assert(qubit < currentState.size());

					const Eigen::Index col = (currentState[qubit] ? 1 : 0);

					for (Eigen::Index row = 0; row < 2; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit] = row == 1;

							const auto localAmplitude = val * PropagateForward(nextState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
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
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;

							const auto localAmplitude = val * PropagateForward(nextState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
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
						if (std::norm(val) > epsilon)
						{
							std::vector<bool> nextState{ currentState };
							nextState[qubit1] = (row & 1) == 1;
							nextState[qubit2] = (row & 2) == 2;
							nextState[qubit3] = (row & 4) == 4;

							const auto localAmplitude = val * PropagateForward(nextState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
						}
					}
				}

				return amplitude;
			}

			// TODO: this could be improved, this sums the max possible considering only the number of qubits that the gates act on, not the action of the particular gates
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

			std::vector<QC::Gates::AppliedGate<>> circuit;
			std::vector<QC::Gates::AppliedGate<>> circuitBack;
			std::unordered_map<std::vector<bool>, std::complex<double>> intermediateAmplitudes;

			size_t doublingsLimit = std::numeric_limits<size_t>::max();
			double epsilon = 1e-15;
		};
	}
}