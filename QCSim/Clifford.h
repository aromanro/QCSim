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
				// for each stabilizer generator there is a corresponding destabilizer generator
				// a stabilizer generator anticommutes with the corresponding destabilizer generator
				// but commutes with all the other destabilizer generators
				// this is preserved during the simulation
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
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyH(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyH(qubit, q);
				}
			}

			void ApplyK(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
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
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyS(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyS(qubit, q);
				}
			}

			void ApplySdg(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplySx(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 512)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
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
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplyX(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyX(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyX(qubit, q);
				}
			}

			void ApplyY(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyY(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyY(qubit, q);
				}
			}

			void ApplyZ(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyZ(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyZ(qubit, q);
				}
			}

			void ApplyCX(size_t target, size_t control)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyCX(target, control, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyCX(target, control, q);
				}
			}

			void ApplyCY(size_t target, size_t control)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
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
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyH(target, q);
						ApplyCX(target, control, q);
						ApplyH(target, q);
					}
				}
			}

			void ApplySwap(size_t qubit1, size_t qubit2)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
					}
				}
			}

			void ApplyISwap(size_t qubit1, size_t qubit2)
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits < 512)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
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
					// grab the first anticommuting generator
					Generator& h = stabilizerGenerators[p];
					// then multiply each other anticommuting one (stabilizer or destabilizer) with it 
					for (size_t q = 0; q < getNrQubits(); ++q)
					{
						if (p == q) continue;

						if (destabilizerGenerators[q].X[qubit])
							rowsum(destabilizerGenerators[q], h);

						if (stabilizerGenerators[q].X[qubit])
							rowsum(stabilizerGenerators[q], h);
					}

					destabilizerGenerators[p] = h;

					h.Clear();
					h.Z[qubit] = true;
					h.PhaseSign = rnd(gen);

					return h.PhaseSign;
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
				const size_t nrQubits = getNrQubits();
				std::vector<bool> state(nrQubits);

				for (size_t i = 0; i < nrQubits; ++i)
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

				if (countRandomQubits == 0) // if there is no random result and we reached here, the probability will stay 1
					return 1.;
				else if (countRandomQubits == 1)
					return 0.5;

				double prob = 1.0;

				// we're going to modify the generators, so let's save the current state, to be restored at the end
				auto saveDest = destabilizerGenerators;
				auto saveStab = stabilizerGenerators;

				do {
					prob *= 0.5; // a random qubit has the 0.5 probability

					Generator& h = stabilizerGenerators[firstP];

					for (size_t q = 0; q < nrQubits; ++q)
					{
						if (firstP == q) continue;

						if (destabilizerGenerators[q].X[firstRandomQubit])
							rowsum(destabilizerGenerators[q], h);

						if (stabilizerGenerators[q].X[firstRandomQubit])
							rowsum(stabilizerGenerators[q], h);
					}

					destabilizerGenerators[firstP] = h;

					h.Clear();
					h.Z[firstRandomQubit] = true;

					// set the measured outcome to the expected value for this state
					h.PhaseSign = state[firstRandomQubit];

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

			void RestoreSavedStateDestructive()
			{
				destabilizerGenerators.swap(savedDestabilizerGenerators);
				stabilizerGenerators.swap(savedStabilizerGenerators);
				ClearSavedState();
			}

			void ClearSavedState()
			{
				savedDestabilizerGenerators.clear();
				savedStabilizerGenerators.clear();
			}

		private:
			inline bool IsRandomResult(size_t qubit, size_t& p) const
			{
				for (size_t q = 0; q < getNrQubits(); ++q)
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
				const size_t nrQubits = getNrQubits();

				// no change to generators, just need to compute the sign in order to get the measurement result
				Generator h(nrQubits);

				// all the stabilizer generators for which the corresponding destabilizer anticommuntes with Z are multiplied together
				// if this is called, all stabilizer generators commute with Z, by the way
				for (size_t q = 0; q < nrQubits; ++q)
					if (destabilizerGenerators[q].X[qubit])
						rowsum(h, stabilizerGenerators[q]);

				return h.PhaseSign;
			}

			inline static bool XOR(bool a, bool b)
			{
				return a != b;
			}

			inline void ApplyH(size_t qubit, size_t q)
			{
				// how does this work?
				// looking at it might not reveal immediately what it does
				// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

				// if we have both X and Z, then it's an Y (with some global phase, given by the sign, Y = iXZ)
				// a Y is transformed to a -Y, so a sign change is needed

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

			// returns the exponent of the i that multiplies the product of the two corresponding Pauli matrices
			// 0, 1 or -1
			// for example for x1 = 1, z1 = 0, x2 = 0, z2 = 1
			// we have X * Z = -i Y so the result should be -1
			// for x1 = 1, z1 = 0, x2 = 1, z2 = 0
			// we have X * X = I so the result should be 0
			static inline int g(int x1, int z1, int x2, int z2)
			{
				if (0 == x1 && 0 == z1) return 0; // I for the 1st generator, 0 exponent no matter what the second generator is
				else if (1 == x1)
				{
					if (1 == z1) return z2 - x2;

					return z2 * (2 * x2 - 1);
				}

				return x2 * (1 - 2 * z2);
			}

			// multiplies the two generators and stores the result in the first one
			static inline void rowsum(Generator& h, Generator& j)
			{
				const size_t nrQubits = h.X.size();
				// phase sign is negative when 'PhaseSign' is true
				// 2 because i^2 = -1
				long long int m = (h.PhaseSign ? 2 : 0) + (j.PhaseSign ? 2 : 0);

				if (nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						const int x1 = j.X[q] ? 1 : 0;
						const int z1 = j.Z[q] ? 1 : 0;
						const int x2 = h.X[q] ? 1 : 0;
						const int z2 = h.Z[q] ? 1 : 0;

						// add up all the exponents of i that contribute to the sign of the product
						m += g(x1, z1, x2, z2);

						// X * X = I, Z * Z = I, so the value is set when there is only one of them
						h.X[q] = (x1 ^ x2) == 1;
						h.Z[q] = (z1 ^ z2) == 1;
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

					long long int mloc = 0;

#pragma omp parallel for reduction(+:mloc) num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						const int x1 = j.X[q] ? 1 : 0;
						const int z1 = j.Z[q] ? 1 : 0;
						const int x2 = h.X[q] ? 1 : 0;
						const int z2 = h.Z[q] ? 1 : 0;

						// add up all the exponents of i that contribute to the sign of the product
						mloc += g(x1, z1, x2, z2);

						// X * X = I, Z * Z = I, so the value is set when there is only one of them
						h.X[q] = (x1 ^ x2) == 1;
						h.Z[q] = (z1 ^ z2) == 1;
					}

					m += mloc;
				}
				
				assert(m % 4 == 0 || m % 4 == 2 || m % 4 == -2);
;
				h.PhaseSign = m % 4 != 0;
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
