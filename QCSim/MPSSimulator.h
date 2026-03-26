#pragma once

#include "MPSSimulatorImpl.h"

#include <functional>
#include <vector>
#include <unordered_set>

namespace QC
{

	namespace TensorNetworks
	{
		class MPSSimulatorState : public MPSSimulatorBaseState
		{
		public:
			MPSSimulatorState() = default;

			MPSSimulatorState(const MPSSimulatorState&) = default;
			MPSSimulatorState(MPSSimulatorState&&) = default;

			MPSSimulatorState& operator=(const MPSSimulatorState&) = default;
			MPSSimulatorState& operator=(MPSSimulatorState&&) = default;

			std::vector<MPSSimulatorInterface::IndexType> qubitsMap;
			std::vector<MPSSimulatorInterface::IndexType> qubitsMapInv;
		};


		// this is going to allow two qubit gates to be applied on qubits that are not adjacent
		// it's a sort of decorator pattern, contains the simulator inside and exposes the same interface
		// adding on top of the implementation the qubits mapping
		class MPSSimulator : public MPSSimulatorInterface
		{
		public:
			// Callback signature: (qubitsMap, bondDims) -> meetPosition
			using MeetingPositionCallback = std::function<IndexType(
				const std::vector<IndexType>&,
				const std::vector<IndexType>&)>;

			MPSSimulator() = delete;

			MPSSimulator(size_t N, unsigned int addseed = 0)
				: impl(N, addseed)
			{
				InitQubitsMap();
			}

			size_t getNrQubits() const final
			{
				return impl.getNrQubits();
			}

			void Clear() override
			{
				impl.Clear();
				InitQubitsMap();
			}

			void InitOnesState() override
			{
				impl.InitOnesState();
				InitQubitsMap();
			}

			void setToQubitState(IndexType q) override
			{
				impl.setToQubitState(q);
				InitQubitsMap();
			}

			void setToBasisState(size_t State) override
			{
				impl.setToBasisState(State);
				InitQubitsMap();
			}

			void setToBasisState(const std::vector<bool>& State) override
			{
				impl.setToBasisState(State);
				InitQubitsMap();
			}

			void SetInitialQubitsMap(const std::vector<long long int>& initialMap) {
				assert(initialMap.size() == impl.getNrQubits());

				for (size_t i = 0; i < initialMap.size(); ++i)
				{
					qubitsMap[i] = initialMap[i];
					qubitsMapInv[initialMap[i]] = i;
				}
			}

			void setLimitBondDimension(IndexType chival) override
			{
				impl.setLimitBondDimension(chival);
			}

			void setLimitEntanglement(double svdThreshold) override
			{
				impl.setLimitEntanglement(svdThreshold);
			}

			void dontLimitBondDimension() override
			{
				impl.dontLimitBondDimension();
			}

			void dontLimitEntanglement() override
			{
				impl.dontLimitEntanglement();
			}

			VectorClass getRegisterStorage() const override
			{
				const auto statev = impl.getRegisterStorage();

				VectorClass res(statev.size());

				// the real/internal qubits are actually in some other positions than 'seen' from outside
				// the correspondence is in the qubitsMap
				// they need to be translated from the internal state to the external state
				
				for (IndexType state = 0; state < res.size(); ++state)
				{
					size_t tmp = state;

					size_t mappedState = 0;
					for (IndexType i = 0; i < static_cast<IndexType>(getNrQubits()); ++i)
					{
						if (tmp & 1ULL)
						{
							const IndexType mappedQubit = qubitsMap.at(i);
							mappedState |= (1ULL << mappedQubit);
						}

						tmp >>= 1;
					}

					res[state] = statev[mappedState];
				}

				return res;
			}

			void print() const override
			{
				impl.print();

				std::cout << "Qubits map: ";
				for (int q = 0; q < static_cast<int>(qubitsMap.size()); ++q)
					std::cout << q << "->" << qubitsMap[q] << " ";
				std::cout << std::endl;
			}

			void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) override
			{
				ApplyGate(gate, gate.getQubit1(), gate.getQubit2());
			}

			void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");
				else if (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = qubitsMap[controllingQubit1];

				// for two qubit gates:
				// if the qubits are not adjacent, apply swap gates until they are
				// don't forget to update the qubitsMap
				if (gate.getQubitsNumber() > 1 && abs(qubit1 - qubit2) > 1)
				{
					if (meetingPositionCallback)
					{
						// External lookahead optimizer supplies the meeting position
						const IndexType meetPos = meetingPositionCallback(
							qubitsMap, impl.getBondDimensions());
						SwapQubitsToPosition(qubit, controllingQubit1, meetPos);
					}
					else if (useOptimalMeetingPosition)
					{
						// Use actual bond dimensions for immediate cost optimization
						const IndexType meetPos = FindBestMeetingPositionLocal(
							qubit, controllingQubit1);
						SwapQubitsToPosition(qubit, controllingQubit1, meetPos);
					}
					else
					{
						// Existing heuristic
						SwapQubits(qubit, controllingQubit1);
					}
					qubit1 = qubitsMap[qubit];
					qubit2 = qubitsMap[controllingQubit1];
					assert(abs(qubit1 - qubit2) == 1);
				}

				impl.ApplyGate(gate, qubit1, qubit2);
			}

			void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) override
			{
				for (const auto& gate : gates)
					ApplyGate(gate);
			}

			bool MeasureQubit(IndexType qubit) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");

				return impl.MeasureQubit(qubitsMap[qubit]);
			}

			std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) override
			{
				std::set<IndexType> mappedQubits;
				for (const auto qubit : qubits)
					mappedQubits.insert(qubitsMap[qubit]);
				
				const auto measuredQubits = impl.MeasureQubits(mappedQubits);

				std::unordered_map<IndexType, bool> res;
				for (const auto& [qubit, val] : measuredQubits)
					res[qubitsMapInv[qubit]] = val;

				return res;
			}

			std::unordered_map<IndexType, bool> MeasureNoCollapse() override
			{
				const auto measuredQubits = impl.MeasureNoCollapse();
				const size_t nrQubits = measuredQubits.size();

				std::unordered_map<IndexType, bool> res;
				for (const auto& [qubit, val] : measuredQubits)
					res[qubitsMapInv[qubit]] = val;
				
				return res;
			}

			double GetProbability(IndexType qubit, bool zeroVal = true) const override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");

				return impl.GetProbability(qubitsMap.at(qubit), zeroVal);
			}

			std::complex<double> getBasisStateAmplitude(size_t State) const override
			{
				std::vector<bool> state(getNrQubits());

				for (size_t i = 0; i < state.size(); ++i)
				{
					state[i] = (State & 1) == 1;
					State >>= 1;
				}

				return getBasisStateAmplitude(state);
			}

			std::complex<double> getBasisStateAmplitude(std::vector<bool>& State) const override
			{
				std::vector<bool> state(getNrQubits(), false);

				for (size_t q = 0; q < State.size(); ++q)
					state[qubitsMap.at(q)] = State[q];

				return impl.getBasisStateAmplitude(state);
			}

			double getBasisStateProbability(size_t State) const override
			{
				return std::norm(getBasisStateAmplitude(State));
			}

			double getBasisStateProbability(std::vector<bool>& State) const override
			{
				return std::norm(getBasisStateAmplitude(State));
			}

			std::shared_ptr<MPSSimulatorStateInterface> getState() const override
			{
				auto baseState = std::static_pointer_cast<MPSSimulatorBaseState>(impl.getState());
				if (!baseState) return nullptr;

				auto state = std::make_shared<MPSSimulatorState>();
				state->gammas.swap(baseState->gammas);
				state->lambdas.swap(baseState->lambdas);
				state->qubitsMap = qubitsMap;
				state->qubitsMapInv = qubitsMapInv;
				return state;
			}

			void setState(const std::shared_ptr<MPSSimulatorStateInterface>& state) override
			{
				if (!state) return;

				auto baseState = std::static_pointer_cast<MPSSimulatorBaseState>(state);
				impl.setState(baseState);

				auto simState = std::static_pointer_cast<MPSSimulatorState>(state);
				qubitsMap = simState->qubitsMap;
				qubitsMapInv = simState->qubitsMapInv;
			}

			void setStateDestructive(std::shared_ptr<MPSSimulatorStateInterface>& state) override
			{
				if (!state) return;
				impl.setStateDestructive(state);

				auto simState = std::static_pointer_cast<MPSSimulatorState>(state);
				qubitsMap.swap(simState->qubitsMap);
				qubitsMapInv.swap(simState->qubitsMapInv);
				
				state.reset();
			}

			void SaveState()
			{
				savedState = getState();
			}

			void RestoreState()
			{
				setState(savedState);
			}

			void RestoreStateDestructive()
			{
				if (!savedState) return;
				impl.setStateDestructive(savedState);

				auto simState = std::static_pointer_cast<MPSSimulatorState>(savedState);
				qubitsMap = simState->qubitsMap;
				qubitsMapInv = simState->qubitsMapInv;
				savedState.reset();
			}

			std::unique_ptr<MPSSimulator> Clone() const
			{
				auto sim = std::make_unique<MPSSimulator>(1);

				sim->qubitsMap = qubitsMap;
				sim->qubitsMapInv = qubitsMapInv;
				sim->impl.limitSize = impl.limitSize;
				sim->impl.limitEntanglement = impl.limitEntanglement;
				sim->impl.chi = impl.chi;
				sim->impl.singularValueThreshold = impl.singularValueThreshold;

				sim->impl.lambdas = impl.lambdas;
				sim->impl.gammas = impl.gammas;

				sim->useOptimalMeetingPosition = useOptimalMeetingPosition;
				sim->meetingPositionCallback = meetingPositionCallback;

				if (savedState) {
					auto simState = std::static_pointer_cast<MPSSimulatorState>(savedState);
					auto clonedSavedState = std::make_shared<MPSSimulatorState>();
					
					clonedSavedState->gammas = simState->gammas;
					clonedSavedState->lambdas = simState->lambdas;
					clonedSavedState->qubitsMap = simState->qubitsMap;
					clonedSavedState->qubitsMapInv = simState->qubitsMapInv;
					
					sim->savedState = clonedSavedState;
				}

				return sim;
			}

			// Enable/disable immediate meeting position optimization
			// using actual bond dimensions (no lookahead, just cheapest path)
			void SetUseOptimalMeetingPosition(bool enable)
			{
				useOptimalMeetingPosition = enable;
			}

			// Set a callback for external lookahead-based meeting position.
			// When set, this takes priority over both the heuristic and the
			// local optimizer.  Pass nullptr to clear.
			void SetMeetingPositionCallback(MeetingPositionCallback callback)
			{
				meetingPositionCallback = std::move(callback);
			}

			// Get actual bond dimensions from the underlying simulator
			std::vector<IndexType> getBondDimensions() const
			{
				return impl.getBondDimensions();
			}

			void MoveAtBeginningOfChain(const std::set<IndexType>& qubits) override
			{
				std::unordered_set<IndexType> handledQubits;

				IndexType currentQubitPos = 0;
				// skip over all qubits that are already in the right position, that is, at the beginning of the chain
				for (; currentQubitPos < static_cast<IndexType>(qubitsMap.size()); ++currentQubitPos)
				{
					const IndexType logicalQubit = qubitsMapInv[currentQubitPos];
					if (qubits.find(logicalQubit) == qubits.end())
						break;
					handledQubits.insert(logicalQubit);
				}

				if (handledQubits.size() == qubits.size())
					return;

				for (IndexType pos = currentQubitPos; pos < static_cast<IndexType>(qubitsMap.size()); ++pos)
				{
					const IndexType logicalQubit = qubitsMapInv[pos];
					// is this a qubit that doesn't need to be moved?
					if (qubits.find(logicalQubit) == qubits.end())
						continue;

					// needs to be moved if it's not already in the right position
					const IndexType currentLogicalPosQubit = qubitsMapInv[currentQubitPos];
					SwapQubits(currentLogicalPosQubit, logicalQubit, true);

					// they are brought together, now swap them
					const IndexType movingQubitReal = qubitsMap[logicalQubit];
					const IndexType toQubitReal = qubitsMap[currentLogicalPosQubit];

					assert(abs(movingQubitReal - toQubitReal) == 1);

					impl.ApplyGate(swapGate, movingQubitReal, toQubitReal);

					// update the maps
					qubitsMap[logicalQubit] = toQubitReal;
					qubitsMap[currentLogicalPosQubit] = movingQubitReal;
					qubitsMapInv[toQubitReal] = logicalQubit;
					qubitsMapInv[movingQubitReal] = currentLogicalPosQubit;

					handledQubits.insert(logicalQubit);
					if (handledQubits.size() == qubits.size())
						break;

					++currentQubitPos;
				}
			}

			// does not check for hermicity, that's why it returns a complex number
		    // the caller should ensure the hermicity and extract the real part
			// also (for now, at least) it supports only one qubit ops
			// the problem with two qubit gates is that they will swap qubits around
			// so instead of only saving the state to compute <psi|U|psi>, and then use it to restore the state,
			// it would need to save the state twice, once for restoring and one for computing the expectation value - the last one having the qubits swapped as the one on which the gates are applied
			// anyway, this would be probably used mostly on Pauli strings, so...
			std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) override
			{
				if (gates.empty()) return 1.;

				std::vector<Gates::AppliedGate<MatrixClass>> translatedOps;
				translatedOps.reserve(gates.size());

				for (const auto& gate : gates)
				{
					if (gate.getQubitsNumber() > 1)
						throw std::invalid_argument("Expectation value for ops applied on more than one qubit not supported yet");

					Gates::AppliedGate<MatrixClass> translated(gate);
					translated.setQubit1(qubitsMap[gate.getQubit1()]);

					translatedOps.push_back(std::move(translated));
				}

				return impl.ExpectationValue(translatedOps);
			}

			std::complex<double> ProjectOnZero() const override
			{
				return impl.ProjectOnZero();
			}

			// needs calling MoveAtBeginningOfChain before this (with the same qubits, of course), otherwise it will give wrong results
			std::unordered_map<IndexType, bool> MeasureNoCollapse(const std::set<IndexType>& qubits) override
			{
				if (qubits.empty()) return {};

				std::set<IndexType> mappedQubits;
				for (const auto qubit : qubits)
					mappedQubits.insert(qubitsMap[qubit]);

				const auto measuredQubits = impl.MeasureNoCollapse(mappedQubits);

				std::unordered_map<IndexType, bool> res;
				for (const auto& [qubit, val] : measuredQubits)
					res[qubitsMapInv[qubit]] = val;

				return res;
			}

		private:
			void InitQubitsMap()
			{
				qubitsMap.resize(getNrQubits());
				qubitsMapInv.resize(getNrQubits());

				for (IndexType i = 0; i < static_cast<IndexType>(getNrQubits()); ++i)
					qubitsMapInv[i] = qubitsMap[i] = i;
			}

			void SwapQubits(IndexType qubit1, IndexType qubit2, bool forceSwapDown = false)
			{
				// if the qubits are not adjacent, apply swap gates until they are
				// don't forget to update the qubitsMap
				IndexType realq1 = qubitsMap[qubit1];
				IndexType realq2 = qubitsMap[qubit2];
				if (realq1 > realq2)
				{
					std::swap(realq1, realq2);
					std::swap(qubit1, qubit2);
				}

				if (realq2 - realq1 <= 1) return;

				//const bool swapDown =
				//	realq1 + realq2 <= static_cast<IndexType>(qubitsMap.size()) - 1;

				if (!forceSwapDown)
				{
					const IndexType mid = (qubitsMap.size() - 1) >> 1;
					if (realq1 < mid && realq2 > mid) // is the middle between the two qubits?
					{
						const IndexType mappedMid = qubitsMapInv[mid];
						SwapQubits(qubit1, mappedMid); // this brings qubit1 near the middle
						realq1 = qubitsMap[qubit1];
						// the other qubit is above the middle, so it won't be affected by the swap
						// the code that follows will bring qubit2 in the middle
					} // otherwise the qubit that's near an end of the chain will be moved towards the other qubit
				}

				// this is just a heuristic, better solutions that minimize the number of swaps would be possible
				const bool swapDown = forceSwapDown ? true : static_cast<IndexType>(qubitsMap.size()) - realq2 <= realq1;

				const IndexType targetQubitReal = swapDown ? realq1 + 1 : realq2 - 1;
				IndexType movingQubitReal = swapDown ? realq2 : realq1;
				const IndexType movingQubitInv = swapDown ? qubit2 : qubit1;

				do
				{
					const IndexType toQubitReal = movingQubitReal + (swapDown ? -1 : 1);
					const IndexType toQubitInv = qubitsMapInv[toQubitReal];

					impl.ApplyGate(swapGate, movingQubitReal, toQubitReal);

					qubitsMap[toQubitInv] = movingQubitReal;
					qubitsMapInv[movingQubitReal] = toQubitInv;

					qubitsMap[movingQubitInv] = toQubitReal;
					qubitsMapInv[toQubitReal] = movingQubitInv;

					movingQubitReal = toQubitReal;
				} while (movingQubitReal != targetQubitReal);

				assert(abs(qubitsMap[qubit1] - qubitsMap[qubit2]) == 1);
				/*
				if (abs(qubitsMap[qubit1] - qubitsMap[qubit2]) != 1)
				{
					std::cerr << "Error: qubits not adjacent after SwapQubits" << std::endl;
					exit(1);
				}
				*/
			}

			// Swap two logical qubits so they meet at a specified bond position.
			// meetPosition is the bond index (in real/chain coordinates) where
			// the two qubits will end up adjacent: one at meetPosition, the other
			// at meetPosition+1.
			void SwapQubitsToPosition(IndexType qubit1, IndexType qubit2,
									  IndexType meetPosition)
			{
				IndexType realq1 = qubitsMap[qubit1];
				IndexType realq2 = qubitsMap[qubit2];
				if (realq1 > realq2)
				{
					std::swap(realq1, realq2);
					std::swap(qubit1, qubit2);
				}

				if (realq2 - realq1 <= 1) return;

				assert(meetPosition >= realq1 && meetPosition < realq2);

				// Move lower qubit (qubit1) rightward from realq1 to meetPosition
				{
					IndexType movingReal = realq1;
					while (movingReal < meetPosition)
					{
						const IndexType toReal = movingReal + 1;
						const IndexType toInv = qubitsMapInv[toReal];

						impl.ApplyGate(swapGate, movingReal, toReal);

						qubitsMap[toInv] = movingReal;
						qubitsMapInv[movingReal] = toInv;

						qubitsMap[qubit1] = toReal;
						qubitsMapInv[toReal] = qubit1;

						movingReal = toReal;
					}
				}

				// Move upper qubit (qubit2) leftward from realq2 to meetPosition+1
				{
					IndexType movingReal = realq2;
					while (movingReal > meetPosition + 1)
					{
						const IndexType toReal = movingReal - 1;
						const IndexType toInv = qubitsMapInv[toReal];

						impl.ApplyGate(swapGate, toReal, movingReal);

						qubitsMap[toInv] = movingReal;
						qubitsMapInv[movingReal] = toInv;

						qubitsMap[qubit2] = toReal;
						qubitsMapInv[toReal] = qubit2;

						movingReal = toReal;
					}
				}

				assert(abs(qubitsMap[qubit1] - qubitsMap[qubit2]) == 1);
				/*
				if (abs(qubitsMap[qubit1] - qubitsMap[qubit2]) != 1)
				{
					std::cerr << "Error: qubits not adjacent after SwapQubits" << std::endl;
					exit(1);
				}
				*/
			}


			// Heuristic: swap towards positions in chain with smaller bond dimensions, the cost of a generic gate is bigger than the cost of a swap
			IndexType FindBestMeetingPositionLocal(IndexType logicalQ1, IndexType logicalQ2) const
			{
				IndexType realq1 = qubitsMap[logicalQ1];
				IndexType realq2 = qubitsMap[logicalQ2];
				if (realq1 > realq2) std::swap(realq1, realq2);

				if (realq2 - realq1 <= 1) return realq1;

				const auto bondDims = impl.getBondDimensions();

				IndexType bestPos = realq1;
				IndexType bestBond = bondDims[realq1];
				
				for (IndexType m = realq1 + 1; m < realq2; ++m)
				{
					if (bondDims[m] < bestBond)
					{
						bestBond = bondDims[m];
						bestPos = m;
					}
				}

				return bestPos;
			}


			MPSSimulatorImpl impl;
			std::vector<IndexType> qubitsMap;
			std::vector<IndexType> qubitsMapInv;
			QC::Gates::SwapGate<MatrixClass> swapGate;

			bool useOptimalMeetingPosition = true;

			MeetingPositionCallback meetingPositionCallback;

			std::shared_ptr<MPSSimulatorStateInterface> savedState;
		};

	}

}