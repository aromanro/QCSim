#pragma once

#include "MPSSimulatorImpl.h"

#include <unordered_map>

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

			std::unordered_map<MPSSimulatorInterface::IndexType, MPSSimulatorInterface::IndexType> qubitsMap;
			std::unordered_map<MPSSimulatorInterface::IndexType, MPSSimulatorInterface::IndexType> qubitsMapInv;
		};


		// this is going to allow two qubit gates to be applied on qubits that are not adjacent
		class MPSSimulator : public MPSSimulatorInterface
		{
		public:
			MPSSimulator() = delete;

			MPSSimulator(size_t N, int addseed = 0)
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

			void setLimitBondDimension(IndexType chival) override
			{
				impl.setLimitBondDimension(chival);
				InitQubitsMap();
			}

			void setLimitEntanglement(double svdThreshold) override
			{
				impl.setLimitEntanglement(svdThreshold);
				InitQubitsMap();
			}

			void dontLimitBondDimension() override
			{
				impl.dontLimitBondDimension();
				InitQubitsMap();
			}

			void dontLimitEntanglement() override
			{
				impl.dontLimitEntanglement();
				InitQubitsMap();
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
							const IndexType mappedQubit = qubitsMapInv.at(i);
							mappedState |= (1ULL << mappedQubit);
						}

						tmp >>= 1;
					}

					res[mappedState] = statev[state];
				}

				return res;
			}

			void print() const override
			{
				impl.print();

				std::cout << "Qubits map: ";
				for (const auto [q1, q2] : qubitsMap)
					std::cout << q1 << "->" << q2 << " ";
				std::cout << std::endl;
			}

			void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) override
			{
				ApplyGate(gate, gate.getQubit1(), gate.getQubit2());
			}

			void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::runtime_error("Qubit index out of bounds");
				else if (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::runtime_error("Qubit index out of bounds");

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = qubitsMap[controllingQubit1];

				// for two qubit gates:
				// if the qubits are not adjacent, apply swap gates until they are
				// don't forget to update the qubitsMap
				if (gate.getQubitsNumber() > 1 && abs(qubit1 - qubit2) > 1)
				{
					SwapQubits(qubit, controllingQubit1);
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
					throw std::runtime_error("Qubit index out of bounds");

				return impl.MeasureQubit(qubitsMap[qubit]);
			}

			double GetProbability(IndexType qubit, bool zeroVal = true) const override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::runtime_error("Qubit index out of bounds");

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
				auto state = std::make_shared<MPSSimulatorState>();
				state->gammas.swap(baseState->gammas);
				state->lambdas.swap(baseState->lambdas);
				state->qubitsMap = qubitsMap;
				state->qubitsMapInv = qubitsMapInv;
				return state;
			}

			void setState(const std::shared_ptr<MPSSimulatorStateInterface>& state) override
			{
				auto baseState = std::static_pointer_cast<MPSSimulatorBaseState>(state);
				impl.setState(baseState);

				auto simState = std::static_pointer_cast<MPSSimulatorState>(state);
				qubitsMap = simState->qubitsMap;
				qubitsMapInv = simState->qubitsMapInv;
			}

		private:
			void InitQubitsMap()
			{
				qubitsMap.clear();
				qubitsMapInv.clear();

				for (IndexType i = 0; i < static_cast<IndexType>(getNrQubits()); ++i)
					qubitsMapInv[i] = qubitsMap[i] = i;
			}

			void SwapQubits(IndexType qubit1, IndexType qubit2)
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

				swapDown = !swapDown;
			}

			MPSSimulatorImpl impl;
			std::unordered_map<IndexType, IndexType> qubitsMap;
			std::unordered_map<IndexType, IndexType> qubitsMapInv;
			QC::Gates::SwapGate<MatrixClass> swapGate;
			bool swapDown = false;
		};

	}

}