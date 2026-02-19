#pragma once

#include <random>
#include <cassert>

#include "QubitRegisterCalculator.h"
#include "PauliStringXZ.h"

namespace QC {
	namespace Clifford {

		class StabilizerState {
		public:
			using Generator = PauliStringXZWithSign;

			StabilizerState()
				: gen(std::random_device{}()), rnd(0.5)
			{
			}

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

			// copy and move ctors and assignment operators
			StabilizerState(const StabilizerState& other)
				: gen(std::random_device{}()), rnd(0.5)
			{
				destabilizerGenerators = other.destabilizerGenerators;
				stabilizerGenerators = other.stabilizerGenerators;
				savedDestabilizerGenerators = other.savedDestabilizerGenerators;
				savedStabilizerGenerators = other.savedStabilizerGenerators;
				enableMultithreading = other.enableMultithreading;
			}

			StabilizerState(StabilizerState&& other) noexcept
				: gen(std::random_device{}()), rnd(0.5)
			{
				destabilizerGenerators.swap(other.destabilizerGenerators);
				stabilizerGenerators.swap(other.stabilizerGenerators);
				savedDestabilizerGenerators.swap(other.savedDestabilizerGenerators);
				savedStabilizerGenerators.swap(other.savedStabilizerGenerators);
				enableMultithreading = other.enableMultithreading;
			}

			StabilizerState& operator=(const StabilizerState& other)
			{
				if (this != &other)
				{
					destabilizerGenerators = other.destabilizerGenerators;
					stabilizerGenerators = other.stabilizerGenerators;
					savedDestabilizerGenerators = other.savedDestabilizerGenerators;
					savedStabilizerGenerators = other.savedStabilizerGenerators;
					enableMultithreading = other.enableMultithreading;
				}
				return *this;
			}

			StabilizerState& operator=(StabilizerState&& other) noexcept
			{
				if (this != &other)
				{
					destabilizerGenerators.swap(other.destabilizerGenerators);
					stabilizerGenerators.swap(other.stabilizerGenerators);
					savedDestabilizerGenerators.swap(other.savedDestabilizerGenerators);
					savedStabilizerGenerators.swap(other.savedStabilizerGenerators);
					enableMultithreading = other.enableMultithreading;
				}
				return *this;
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

					// X present means it's either X or Y, both of which anticommute with Z
					// otherwise it's either I or Z, which commute with Z
					for (size_t q = p + 1; q < getNrQubits(); ++q)
						if (stabilizerGenerators[q].X[qubit])
							stabilizerGenerators[q].Multiply(h, enableMultithreading);

					for (size_t q = 0; q < getNrQubits(); ++q)
						if (destabilizerGenerators[q].X[qubit])
							destabilizerGenerators[q].Multiply(h, enableMultithreading);

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
				if (countRandomQubits == 0 || prob == 0.0)
					return prob;

				// we're going to modify the generators, so let's save the current state, to be restored at the end
				auto saveDest = destabilizerGenerators;
				auto saveStab = stabilizerGenerators;

				do {
					prob *= 0.5; // a random qubit has the 0.5 probability

					Generator& h = stabilizerGenerators[firstP];

					for (size_t q = firstP + 1; q < nrQubits; ++q)
						if (stabilizerGenerators[q].X[firstRandomQubit])
							stabilizerGenerators[q].Multiply(h, enableMultithreading);

					for (size_t q = 0; q < nrQubits; ++q)
						if (destabilizerGenerators[q].X[firstRandomQubit])
							destabilizerGenerators[q].Multiply(h, enableMultithreading);

					destabilizerGenerators[firstP] = h;

					h.Clear();
					h.Z[firstRandomQubit] = true;

					// set the measured outcome to the expected value for this state
					h.PhaseSign = state[firstRandomQubit];
					handledQubits[firstRandomQubit] = true; // not really needed, we won't look back

					++firstRandomQubit;

					countRandomQubits = DealWithDeterministicQubits(state, handledQubits, firstRandomQubit, firstP, prob);
					if (prob == 0.0)
						break;
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
				if (savedDestabilizerGenerators.empty() || savedStabilizerGenerators.empty()) return;
				
				destabilizerGenerators = savedDestabilizerGenerators;
				stabilizerGenerators = savedStabilizerGenerators;
			}

			void RestoreSavedStateDestructive()
			{
				if (savedDestabilizerGenerators.empty() || savedStabilizerGenerators.empty()) return;

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

				// all the stabilizer generators for which the corresponding destabilizer anticommutes with Z are multiplied together
				// if this is called, all stabilizer generators commute with Z, by the way
				for (size_t q = 0; q < nrQubits; ++q)
					if (destabilizerGenerators[q].X[qubit])
						h.Multiply(stabilizerGenerators[q], enableMultithreading);

				return h.PhaseSign;
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
