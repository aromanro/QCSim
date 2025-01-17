#pragma once

#include <random>
#include <cassert>

#include "QubitRegisterCalculator.h"
#include "Generator.h"

namespace QC {
	namespace Clifford {

		class StabilizerState {
		public:
			StabilizerState() = delete;

			explicit StabilizerState(size_t nQubits)
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

			void Reset()
			{
				for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
				{
					for (size_t p = 0; p < stabilizerGenerators.size(); ++p)
					{
						if (p == q)
						{
							destabilizerGenerators[q].X[q] = true;
							stabilizerGenerators[q].Z[q] = true;
						}
						else
						{
							destabilizerGenerators[q].X[p] = false;
							stabilizerGenerators[q].Z[p] = false;
						}

						destabilizerGenerators[q].PhaseSign = false;
						stabilizerGenerators[q].PhaseSign = false;
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
				std::vector<bool> handledQubits(nrQubits, false);

				size_t firstRandomQubit = 0;
				size_t firstP = 0;

				double prob = 1.0;
				size_t countRandomQubits = DealWithDeterministicQubits(state, handledQubits, firstRandomQubit, firstP, prob);
				if (countRandomQubits == 0)
					return prob;

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

					countRandomQubits = DealWithDeterministicQubits(state, handledQubits, firstRandomQubit, firstP, prob);
				} while (countRandomQubits > 0);

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

			void RestoreState()
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

			void SetMultithreading(bool enable = true)
			{
				enableMultithreading = enable;
			}

			bool GetMultithreading() const
			{
				return enableMultithreading;
			}

		protected:
			inline size_t DealWithDeterministicQubits(const std::vector<bool>& state, std::vector<bool>& handledQubits, size_t& firstRandomQubit, size_t& firstP, double& prob)
			{
				const size_t nrQubits = getNrQubits();
				size_t p;
				size_t countRandomQubits = 0;

				// first deal with the deterministic qubits, it might turn out that the probability is 0
				// in that case, we can return immediately
				for (size_t qubit = firstRandomQubit; qubit < nrQubits; ++qubit)
				{
					if (handledQubits[qubit]) continue;
					
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
						{
							prob = 0.;
							return 0;
						}

						handledQubits[qubit] = true;
					}
				}

				if (countRandomQubits == 1)
				{
					prob *= 0.5;
					countRandomQubits = 0;
				}

				return countRandomQubits;
			}

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
			inline void rowsum(Generator& h, Generator& j)
			{
				const size_t nrQubits = h.X.size();
				// phase sign is negative when 'PhaseSign' is true
				// 2 because i^2 = -1
				long long int m = (h.PhaseSign ? 2 : 0) + (j.PhaseSign ? 2 : 0);

				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						const int x1 = BoolToInt(j.X[q]);
						const int z1 = BoolToInt(j.Z[q]);
						const int x2 = BoolToInt(h.X[q]);
						const int z2 = BoolToInt(h.Z[q]);

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
						const int x1 = BoolToInt(j.X[q]);
						const int z1 = BoolToInt(j.Z[q]);
						const int x2 = BoolToInt(h.X[q]);
						const int z2 = BoolToInt(h.Z[q]);

						// add up all the exponents of i that contribute to the sign of the product
						mloc += g(x1, z1, x2, z2);

						// X * X = I, Z * Z = I, so the value is set when there is only one of them
						h.X[q] = (x1 ^ x2) == 1;
						h.Z[q] = (z1 ^ z2) == 1;
					}

					m += mloc;
				}

				// the mod 4 that appears here is because the values for the powers of i keep repeating
				assert(m % 4 == 0 || m % 4 == 2 || m % 4 == -2);

				h.PhaseSign = m % 4 != 0;
			}

			static inline int BoolToInt(bool b)
			{
				return b ? 1 : 0;
			}

			std::vector<Generator> destabilizerGenerators;
			std::vector<Generator> stabilizerGenerators;

			std::vector<Generator> savedDestabilizerGenerators;
			std::vector<Generator> savedStabilizerGenerators;

			std::default_random_engine gen;
			std::bernoulli_distribution rnd;

			bool enableMultithreading = true;
		};
	}
}
