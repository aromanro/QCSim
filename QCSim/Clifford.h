#pragma once

#include <random>
#include <cassert>

#include "QubitRegisterCalculator.h"
#include "Generator.h"

namespace QC {
	namespace Clifford {

		class StabilizerSimulator {
		public:
			StabilizerSimulator() = delete;

			explicit StabilizerSimulator(size_t nQubits)
				: destabilizerGenerators(nQubits), stabilizerGenerators(nQubits), gen(std::random_device{}()), rnd(0.5)
			{
				// this puts it in the |0> state
				for (size_t q = 0; q < nQubits; ++q)
				{
					destabilizerGenerators[q].Resize(nQubits);
					stabilizerGenerators[q].Resize(nQubits);

					destabilizerGenerators[q].X[q] = true;
					stabilizerGenerators[q].Z[q] = true;
				}
			}

			void ApplyH(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyH(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyH(qubit, q);
				}
			}

			void ApplyK(size_t qubit)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplyS(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyS(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyS(qubit, q);
				}
			}

			void ApplySdg(size_t qubit)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplySx(size_t qubit)
			{
				if (stabilizerGenerators.size() < 512)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 128)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplySxDag(size_t qubit)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplyX(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyX(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyX(qubit, q);
				}
			}

			void ApplyY(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyY(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyY(qubit, q);
				}
			}

			void ApplyZ(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyZ(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyZ(qubit, q);
				}
			}

			void ApplyCX(size_t target, size_t control)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyCX(target, control, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyCX(target, control, q);
				}
			}

			void ApplyCY(size_t target, size_t control)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyZ(target, q);
						ApplyS(target, q);
						ApplyCX(target, control, q);
						ApplyS(target, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyZ(target, q);
						ApplyS(target, q);
						ApplyCX(target, control, q);
						ApplyS(target, q);
					}
				}
			}

			void ApplyCZ(size_t target, size_t control)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyH(target, q);
						ApplyCX(target, control, q);
						ApplyH(target, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyH(target, q);
						ApplyCX(target, control, q);
						ApplyH(target, q);
					}
				}
			}

			void ApplySwap(size_t qubit1, size_t qubit2)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
					}
				}
			}

			void ApplyISwap(size_t qubit1, size_t qubit2)
			{
				if (stabilizerGenerators.size() < 512)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyS(qubit1, q);
						ApplyS(qubit2, q);
						ApplyH(qubit1, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
						ApplyH(qubit2, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 128)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyS(qubit1, q);
						ApplyS(qubit2, q);
						ApplyH(qubit1, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
						ApplyH(qubit2, q);
					}
				}
			}

			bool MeasureQubit(size_t qubit)
			{
				size_t p;
				if (IsRandomResult(qubit, p))
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						if (destabilizerGenerators[q].X[qubit])
							rowsum(destabilizerGenerators[q], stabilizerGenerators[p], true);

						if (p != q && stabilizerGenerators[q].X[qubit])
							rowsum(stabilizerGenerators[q], stabilizerGenerators[p], false);
					}

					destabilizerGenerators[p] = stabilizerGenerators[p];

					stabilizerGenerators[p].Clear();
					stabilizerGenerators[p].Z[qubit] = true;
					stabilizerGenerators[p].PhaseSign = rnd(gen);

					return stabilizerGenerators[p].PhaseSign;
				}

				// case 2 - Z (on measured qubit) commutes with all generators
				// no change to generators, just need to compute the sign in order to get the measurement result
				return GetTheDeterministicOutcome(qubit);
			}

			double GetQubitProbability(size_t qubit)
			{
				size_t p;
				if (IsRandomResult(qubit, p))
					return 0.5;

				return GetTheDeterministicOutcome(qubit) ? 1. : 0.;
			}

			double getBasisStateProbability(size_t State)
			{
				std::vector<bool> state(getNrQubits());

				for (size_t i = 0; i < state.size(); ++i)
				{
					state[i] = (State & 1) == 1;
					State >>= 1;
				}

				return getBasisStateProbability(state);
			}

			double getBasisStateProbability(const std::vector<bool>& state)
			{
				const size_t nrQubits = getNrQubits();
				size_t p;

				std::vector<bool> handledQubits(nrQubits, false);

				size_t firstRandomQubit = 0;
				size_t firstP = 0;

				// first deal with the deterministic qubits, it might turn out that the probability is 0
				// in that case, we can return immediately
				size_t countRandomQubits = 0;
				for (size_t qubit = 0; qubit < nrQubits; ++qubit)
				{
					if (IsRandomResult(qubit, p))
					{
						if (0 == countRandomQubits)
						{
							firstRandomQubit = qubit;
							firstP = p;
						}
						++countRandomQubits;
					}
					else
					{
						if (GetTheDeterministicOutcome(qubit) != state[qubit])
							return 0.;

						handledQubits[qubit] = true;
					}
				}

				if (countRandomQubits == 0)
					return 1.;
				else if (countRandomQubits == 1)
					return 0.5;

				// if there is no random result and we reached here, the probability will stay 1
				double prob = 1.0;

				// we're going to modify the generators, so let's save the current state, to be restored at the end
				auto saveDest = destabilizerGenerators;
				auto saveStab = stabilizerGenerators;

				do {
					prob *= 0.5; // a random qubit has the 0.5 probability

					for (size_t q = 0; q < nrQubits; ++q)
					{
						if (destabilizerGenerators[q].X[firstRandomQubit])
							rowsum(destabilizerGenerators[q], stabilizerGenerators[firstP], true);

						if (firstP != q && stabilizerGenerators[q].X[firstRandomQubit])
							rowsum(stabilizerGenerators[q], stabilizerGenerators[firstP], false);
					}

					destabilizerGenerators[firstP] = stabilizerGenerators[firstP];

					stabilizerGenerators[firstP].Clear();
					stabilizerGenerators[firstP].Z[firstRandomQubit] = true;

					// set the measured outcome to the expected value for this state
					const bool expectedQubitOutcome = state[firstRandomQubit];
					stabilizerGenerators[firstP].PhaseSign = expectedQubitOutcome;

					handledQubits[firstRandomQubit] = true; // not really needed, we won't look back
					countRandomQubits = 0;

					// this measurement might have been turned some not measured yet qubits to deterministic, so let's check them
					// this also checks if we still have some non deterministic qubits left			
					for (size_t qubit = firstRandomQubit + 1; qubit < nrQubits; ++qubit)
					{
						if (!handledQubits[qubit])
						{
							if (IsRandomResult(qubit, p))
							{
								// there is still at least one more random qubit
								if (0 == countRandomQubits)
								{
									firstRandomQubit = qubit;
									firstP = p;
								}
								++countRandomQubits;
							}
							else
							{
								if (GetTheDeterministicOutcome(qubit) != state[qubit])
								{
									// restore the state before returning
									destabilizerGenerators.swap(saveDest);
									stabilizerGenerators.swap(saveStab);

									return 0;
								}

								handledQubits[qubit] = true;
							}
						}
					}

					if (countRandomQubits == 1)
						prob *= 0.5;
				} while (countRandomQubits > 1);

				// we're done, restore the state
				destabilizerGenerators.swap(saveDest);
				stabilizerGenerators.swap(saveStab);

				return prob;
			}

			std::vector<double> AllProbabilities()
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits > 32) throw std::runtime_error("The simulator has too many qubits for computing all probabilities");

				const size_t nrStates = 1ULL << nrQubits;
				std::vector<double> probs(nrStates, 0);

				for (size_t state = 0; state < nrStates; ++state)
					probs[state] = getBasisStateProbability(state);

				return probs;
			}

			size_t getNrQubits() const { return stabilizerGenerators.size(); }

			void SaveState()
			{
				savedDestabilizerGenerators = destabilizerGenerators;
				savedStabilizerGenerators = stabilizerGenerators;
			}

			void RestoreSavedState()
			{
				destabilizerGenerators = savedDestabilizerGenerators;
				stabilizerGenerators = savedStabilizerGenerators;
			}

			void ClearSavedState()
			{
				savedDestabilizerGenerators.clear();
				savedStabilizerGenerators.clear();
			}

		private:
			inline bool IsRandomResult(size_t qubit, size_t& p) const
			{
				for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					if (stabilizerGenerators[q].X[qubit])
					{
						// Z anticommutes with X
						p = q;
						return true;
					}

				return false;
			}

			// call it only in the case 2, Z (on measured qubit) commutes with all generators
			inline bool GetTheDeterministicOutcome(size_t qubit)
			{
				// no change to generators, just need to compute the sign in order to get the measurement result
				Generator h(destabilizerGenerators.size());

				for (size_t q = 0; q < destabilizerGenerators.size(); ++q)
					if (destabilizerGenerators[q].X[qubit])
						rowsum(h, stabilizerGenerators[q], true);

				return h.PhaseSign;
			}

			inline static bool XOR(bool a, bool b)
			{
				//return ((a ? 1 : 0) ^ (b ? 1 : 0)) == 1;
				return a != b;
			}

			inline void ApplyH(size_t qubit, size_t q)
			{
				// how does this work?
				// looking at it might not reveal immediately what it does
				// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

				// if we have both X and Z, then it's an Y (with some global phase, given by the sign, Y = iXZ)
				// an Y is transformed to a -Y, so a sign change is needed

				if (destabilizerGenerators[q].X[qubit] && destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				// swap X and Z
				bool t = destabilizerGenerators[q].X[qubit];
				destabilizerGenerators[q].X[qubit] = destabilizerGenerators[q].Z[qubit];
				destabilizerGenerators[q].Z[qubit] = t;

				if (stabilizerGenerators[q].X[qubit] && stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;

				// swap X and Z
				t = stabilizerGenerators[q].X[qubit];
				stabilizerGenerators[q].X[qubit] = stabilizerGenerators[q].Z[qubit];
				stabilizerGenerators[q].Z[qubit] = t;
			}

			inline void ApplyS(size_t qubit, size_t q)
			{
				if (destabilizerGenerators[q].X[qubit] && destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				destabilizerGenerators[q].Z[qubit] = XOR(destabilizerGenerators[q].Z[qubit], destabilizerGenerators[q].X[qubit]);

				if (stabilizerGenerators[q].X[qubit] && stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;

				stabilizerGenerators[q].Z[qubit] = XOR(stabilizerGenerators[q].Z[qubit], stabilizerGenerators[q].X[qubit]);
			}

			inline void ApplyX(size_t qubit, size_t q)
			{
				if (destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				if (stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
			}

			inline void ApplyY(size_t qubit, size_t q)
			{
				// can be done with ifs, can be done with XORs
				destabilizerGenerators[q].PhaseSign = XOR(destabilizerGenerators[q].PhaseSign, XOR(destabilizerGenerators[q].Z[qubit], destabilizerGenerators[q].X[qubit]));
				stabilizerGenerators[q].PhaseSign = XOR(stabilizerGenerators[q].PhaseSign, XOR(stabilizerGenerators[q].Z[qubit], stabilizerGenerators[q].X[qubit]));
			}

			inline void ApplyZ(size_t qubit, size_t q)
			{
				if (destabilizerGenerators[q].X[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				if (stabilizerGenerators[q].X[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
			}

			inline void ApplyCX(size_t target, size_t control, size_t q)
			{
				destabilizerGenerators[q].PhaseSign = XOR(destabilizerGenerators[q].PhaseSign, destabilizerGenerators[q].X[control] && destabilizerGenerators[q].Z[target] &&
					XOR(destabilizerGenerators[q].X[target], XOR(destabilizerGenerators[q].Z[control], true)));

				destabilizerGenerators[q].X[target] = XOR(destabilizerGenerators[q].X[target], destabilizerGenerators[q].X[control]);
				destabilizerGenerators[q].Z[control] = XOR(destabilizerGenerators[q].Z[control], destabilizerGenerators[q].Z[target]);

				stabilizerGenerators[q].PhaseSign = XOR(stabilizerGenerators[q].PhaseSign, stabilizerGenerators[q].X[control] && stabilizerGenerators[q].Z[target] &&
					XOR(stabilizerGenerators[q].X[target], XOR(stabilizerGenerators[q].Z[control], true)));

				stabilizerGenerators[q].X[target] = XOR(stabilizerGenerators[q].X[target], stabilizerGenerators[q].X[control]);
				stabilizerGenerators[q].Z[control] = XOR(stabilizerGenerators[q].Z[control], stabilizerGenerators[q].Z[target]);
			}

			static inline int g(int x1, int z1, int x2, int z2)
			{
				if (0 == x1 && 0 == z1) return 0;
				else if (1 == x1)
				{
					if (1 == z1) return z2 - x2;

					return z2 * (2 * x2 - 1);
				}

				return x2 * (1 - 2 * z2);
			}

			inline void rowsumDestabilizers(Generator& h, Generator& j, long long int& m)
			{
				m += j.PhaseSign ? 2 : 0;

				if (destabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < destabilizerGenerators.size(); ++q)
					{
						const int x1 = j.X[q] ? 1 : 0;
						const int z1 = j.Z[q] ? 1 : 0;
						const int x2 = h.X[q] ? 1 : 0;
						const int z2 = h.Z[q] ? 1 : 0;
						m += g(x1, z1, x2, z2);

						h.X[q] = XOR(h.X[q], j.X[q]);
						h.Z[q] = XOR(h.Z[q], j.Z[q]);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

					long long int mloc = 0;

#pragma omp parallel for reduction(+:mloc) num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(destabilizerGenerators.size()); ++q)
					{
						const int x1 = j.X[q] ? 1 : 0;
						const int z1 = j.Z[q] ? 1 : 0;
						const int x2 = h.X[q] ? 1 : 0;
						const int z2 = h.Z[q] ? 1 : 0;
						mloc += g(x1, z1, x2, z2);

						h.X[q] = XOR(h.X[q], j.X[q]);
						h.Z[q] = XOR(h.Z[q], j.Z[q]);
					}

					m += mloc;
				}
			}

			inline void rowsumStabilizers(Generator& h, Generator& j, long long int& m)
			{
				m += j.PhaseSign ? 2 : 0;

				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						const int x1 = j.X[q] ? 1 : 0;
						const int z1 = j.Z[q] ? 1 : 0;
						const int x2 = h.X[q] ? 1 : 0;
						const int z2 = h.Z[q] ? 1 : 0;
						m += g(x1, z1, x2, z2);

						h.X[q] = XOR(h.X[q], j.X[q]);
						h.Z[q] = XOR(h.Z[q], j.Z[q]);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

					long long int mloc = 0;

#pragma omp parallel for reduction(+:mloc) num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < stabilizerGenerators.size(); ++q)
					{
						const int x1 = j.X[q] ? 1 : 0;
						const int z1 = j.Z[q] ? 1 : 0;
						const int x2 = h.X[q] ? 1 : 0;
						const int z2 = h.Z[q] ? 1 : 0;
						mloc += g(x1, z1, x2, z2);

						h.X[q] = XOR(h.X[q], j.X[q]);
						h.Z[q] = XOR(h.Z[q], j.Z[q]);
					}

					m += mloc;
				}
			}

			inline void rowsum(Generator& h, Generator& j, bool destabilizers = false)
			{
				long long int m = h.PhaseSign ? 2 : 0;

				if (destabilizers)
					rowsumDestabilizers(h, j, m);
				else
					rowsumStabilizers(h, j, m);

				m %= 4;
				assert(m == 0 || m == 2);

				h.PhaseSign = m != 0;
			}

			std::vector<Generator> destabilizerGenerators;
			std::vector<Generator> stabilizerGenerators;

			std::vector<Generator> savedDestabilizerGenerators;
			std::vector<Generator> savedStabilizerGenerators;

			std::default_random_engine gen;
			std::bernoulli_distribution rnd;
		};
	}
}
