#pragma once

#include <array>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include "QuantumGate.h"

namespace QC {
	namespace PathIntegral {

		struct FastVectorBool {
			static constexpr size_t MaxWords = 16;

			FastVectorBool() = default;

			explicit FastVectorBool(const std::vector<bool>& v) 
				: nBits(v.size())
			{
				for (size_t i = 0; i < v.size(); ++i)
					if (v[i]) words[i / 64] |= (1ULL << (i % 64));
				for (size_t i = (v.size() + 63) / 64; i < MaxWords; ++i)
					words[i] = 0;
			}

			bool get(size_t i) const { return (words[i / 64] >> (i % 64)) & 1; }

			void set(size_t i, bool val)
			{
				if (val) words[i / 64] |= (1ULL << (i % 64));
				else words[i / 64] &= ~(1ULL << (i % 64));
			}

			size_t size() const { return nBits; }
			size_t nWords() const { return (nBits + 63) / 64; }

			bool operator==(const FastVectorBool& o) const
			{
				const size_t w = nWords();
				for (size_t i = 0; i < w; ++i)
					if (words[i] != o.words[i]) return false;
				return true;
			}

			bool operator!=(const FastVectorBool& o) const { return !(*this == o); }

			const std::array<uint64_t, MaxWords>& getWords() const { return words; }

		private:
			std::array<uint64_t, MaxWords> words{};
			size_t nBits = 0;
		};

		struct FastVectorBoolHash {
			size_t operator()(const FastVectorBool& s) const
			{
				const std::array<uint64_t, FastVectorBool::MaxWords>& words = s.getWords();
				size_t seed = 0;
				for (int i = static_cast<int>((s.size() + 63) / 64) - 1; i >= 0; --i)
					seed ^= std::hash<uint64_t>{}(words[i]) + (seed << 6);
				return seed;
			}
		};

		class PathIntegralSimulator {
		public:
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

			std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash>& GetAmplitudes()
			{
				return intermediateAmplitudes;
			}

			void SaveAmplitudes()
			{
				savedAmplitudes = intermediateAmplitudes;
			}

			void RestoreAmplitudes()
			{
				intermediateAmplitudes = savedAmplitudes;
			}

			void SwapAmplitudes	()
			{
				intermediateAmplitudes.swap(savedAmplitudes);
			}

			void ClearSavedAmplitudes()
			{
				savedAmplitudes.clear();
			}

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

						++doublings;
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

						++doublings;

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

				const FastVectorBool startBits(startState);
				const FastVectorBool endBits(endState);

				if (doublingsLimit > 0 && !circuitBack.empty())
				{
					intermediateAmplitudes.clear();

					PropagateBackward(endBits, std::complex<double>(1., 0.));

					return PropagateForward(startBits);
				}

				const size_t possibleQubitsChanges = CountPossibleQubitsChanges(circuit);

				return Propagate(startBits, endBits, 0, possibleQubitsChanges);
			}

			// the following methods are for the case one needs all non-zero amplitudes
			// of course, except the trimmed out ones
			void PropagateAll(const std::vector<QC::Gates::AppliedGate<>>& circuit)
			{
				size_t nQubits = 0;
				for (const auto& gate : circuit)
				{
					nQubits = std::max(nQubits, gate.getQubit1() + 1);
					if (gate.getQubitsNumber() > 1)
						nQubits = std::max(nQubits, gate.getQubit2() + 1);
					if (gate.getQubitsNumber() > 2)
						nQubits = std::max(nQubits, gate.getQubit3() + 1);
				}

				const std::vector<bool> startState(nQubits, false);

				PropagateAll(circuit, startState);
			}

			void PropagateAll(const std::vector<QC::Gates::AppliedGate<>>& circuit, const std::vector<bool>& startState)
			{
				assert(startState.size() > 0);
				const FastVectorBool startBits(startState);
				intermediateAmplitudes.clear();

				circuitBack.clear();
				this->circuit = circuit;

				intermediateAmplitudes = PropagateAll(startBits);
			}

			void PropagateStep(const QC::Gates::AppliedGate<>& gate, std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash>& currentAmplitudes)
			{
				const auto& U = gate.getRawOperatorMatrix();
				const size_t gateQubits = gate.getQubitsNumber();

				assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

				std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash> nextAmplitudes;
				nextAmplitudes.reserve(currentAmplitudes.size() * 2);

				if (gateQubits == 1)
				{
					const size_t qubit = gate.getQubit1();

					for (const auto& [state, amp] : currentAmplitudes)
					{
						if (std::norm(amp) < epsilon) continue;

						assert(qubit < state.size());
						const Eigen::Index col = (state.get(qubit) ? 1 : 0);

						for (Eigen::Index row = 0; row < 2; ++row)
						{
							const std::complex<double> val = U(row, col);
							if (std::norm(val) > epsilon)
							{
								auto nextState = state;
								nextState.set(qubit, row == 1);

								nextAmplitudes[nextState] += val * amp;
							}
						}
					}
				}
				else if (gateQubits == 2)
				{
					const size_t qubit1 = gate.getQubit1();
					const size_t qubit2 = gate.getQubit2();

					for (const auto& [state, amp] : currentAmplitudes)
					{
						if (std::norm(amp) < epsilon) continue;

						assert(qubit1 < state.size() && qubit2 < state.size());
						const Eigen::Index col = ((state.get(qubit2) ? 2 : 0) | (state.get(qubit1) ? 1 : 0));

						for (Eigen::Index row = 0; row < 4; ++row)
						{
							const std::complex<double> val = U(row, col);
							if (std::norm(val) > epsilon)
							{
								auto nextState = state;
								nextState.set(qubit1, (row & 1) == 1);
								nextState.set(qubit2, (row & 2) == 2);

								nextAmplitudes[nextState] += val * amp;
							}
						}
					}
				}
				else // only up to three qubits gates are supported
				{
					const size_t qubit1 = gate.getQubit1();
					const size_t qubit2 = gate.getQubit2();
					const size_t qubit3 = gate.getQubit3();

					for (const auto& [state, amp] : currentAmplitudes)
					{
						if (std::norm(amp) < epsilon) continue;

						assert(qubit1 < state.size() && qubit2 < state.size() && qubit3 < state.size());
						const Eigen::Index col = ((state.get(qubit3) ? 4 : 0) | (state.get(qubit2) ? 2 : 0) | (state.get(qubit1) ? 1 : 0));

						for (Eigen::Index row = 0; row < 8; ++row)
						{
							const std::complex<double> val = U(row, col);
							if (std::norm(val) > epsilon)
							{
								auto nextState = state;
								nextState.set(qubit1, (row & 1) == 1);
								nextState.set(qubit2, (row & 2) == 2);
								nextState.set(qubit3, (row & 4) == 4);

								nextAmplitudes[nextState] += val * amp;
							}
						}
					}
				}

				currentAmplitudes.swap(nextAmplitudes);
			}

		private:
			// TODO: this can be parallelized!
			std::complex<double> Propagate(const FastVectorBool& currentState, const FastVectorBool& endState, size_t gateIndex, size_t possibleQubitsChanges)
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

					const Eigen::Index col = (currentState.get(qubit) ? 1 : 0);

					for (Eigen::Index row = 0; row < 2; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::norm(val) > epsilon)
						{
							auto nextState = currentState;
							nextState.set(qubit, row == 1);

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

					const Eigen::Index col = ((currentState.get(qubit2) ? 2 : 0) | (currentState.get(qubit1) ? 1 : 0));

					for (Eigen::Index row = 0; row < 4; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::norm(val) > epsilon)
						{
							auto nextState = currentState;
							nextState.set(qubit1, (row & 1) == 1);
							nextState.set(qubit2, (row & 2) == 2);

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

					const Eigen::Index col = ((currentState.get(qubit3) ? 4 : 0) | (currentState.get(qubit2) ? 2 : 0) | (currentState.get(qubit1) ? 1 : 0));

					for (Eigen::Index row = 0; row < 8; ++row)
					{
						const std::complex<double> val = U(row, col);
						if (std::norm(val) > epsilon)
						{
							auto nextState = currentState;
							nextState.set(qubit1, (row & 1) == 1);
							nextState.set(qubit2, (row & 2) == 2);
							nextState.set(qubit3, (row & 4) == 4);

							if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
								continue;

							const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
							amplitude += localAmplitude;
						}
					}
				}

				return amplitude;
			}

			void PropagateBackward(const FastVectorBool& endState, std::complex<double> amplitude)
			{
				intermediateAmplitudes.clear();
				intermediateAmplitudes[endState] = amplitude;

				for (const auto& gate : circuitBack)
					PropagateStep(gate, intermediateAmplitudes);
			}

			std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash> PropagateAll(const FastVectorBool& startState)
			{
				std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash> currentAmplitudes;
				currentAmplitudes[startState] = std::complex<double>(1., 0.);

				for (const auto& gate : circuit)
					PropagateStep(gate, currentAmplitudes);

				return currentAmplitudes;
			}

			std::complex<double> PropagateForward(const FastVectorBool& startState)
			{
				const auto currentAmplitudes = PropagateAll(startState);

				std::complex<double> result(0., 0.);
				for (const auto& [state, amp] : currentAmplitudes)
				{
					auto it = intermediateAmplitudes.find(state);
					if (it != intermediateAmplitudes.end())
						result += amp * it->second;
				}

				return result;
			}

			// TODO: this could be improved, this sums the max possible considering only the number of qubits that the gates act on, not the action of the particular gates
			static size_t CountPossibleQubitsChanges(const std::vector<QC::Gates::AppliedGate<>>& circuit)
			{
				size_t count = 0;

				for (const auto& gate : circuit)
					count += gate.getQubitsNumber();

				return count;
			}

			static size_t CountDifferentQubits(const FastVectorBool& curState, const FastVectorBool& endState)
			{
				size_t count = 0;
				const size_t words = curState.nWords();
				for (size_t i = 0; i < words; ++i)
				{
					const uint64_t diff = curState.getWords()[i] ^ endState.getWords()[i];
#ifdef _MSC_VER
					count += __popcnt64(diff);
#else
					count += __builtin_popcountll(diff);
#endif
				}
				return count;
			}

			std::vector<QC::Gates::AppliedGate<>> circuit;
			std::vector<QC::Gates::AppliedGate<>> circuitBack;
			std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash> intermediateAmplitudes;

			std::unordered_map<FastVectorBool, std::complex<double>, FastVectorBoolHash> savedAmplitudes;

			size_t doublingsLimit = std::numeric_limits<size_t>::max();
			double epsilon = 1e-15;
		};
	}
}
