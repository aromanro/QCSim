#pragma once

#include "MPOSimulatorImpl.h"

#include <algorithm>
#include <limits>
#include <vector>
#include <utility>

namespace QC
{

	namespace TensorNetworks
	{
		class MPOSimulatorState : public MPOSimulatorBaseState
		{
		public:
			MPOSimulatorState() = default;

			MPOSimulatorState(const MPOSimulatorState&) = default;
			MPOSimulatorState(MPOSimulatorState&&) = default;

			MPOSimulatorState& operator=(const MPOSimulatorState&) = default;
			MPOSimulatorState& operator=(MPOSimulatorState&&) = default;

			std::vector<MPOSimulatorInterface::IndexType> qubitsMap;
			std::vector<MPOSimulatorInterface::IndexType> qubitsMapInv;
		};

		// Decorator over MPOSimulatorImpl that allows two qubit gates to be applied on
		// non adjacent qubits: it keeps a logical -> physical qubit mapping and inserts
		// swap gates to bring the involved qubits together before applying the gate.
		// This mirrors the MPSSimulator decorator over MPSSimulatorImpl.
		class MPOSimulator : public MPOSimulatorInterface
		{
		public:
			MPOSimulator() = delete;

			MPOSimulator(size_t N, unsigned int addseed = 0)
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

			void setToMixtureOfBasisStates(const std::vector<std::pair<size_t, double>>& mixture) override
			{
				impl.setToMixtureOfBasisStates(mixture);
				InitQubitsMap();
			}

			void setToMixtureOfBasisStates(const std::vector<std::pair<std::vector<bool>, double>>& mixture) override
			{
				impl.setToMixtureOfBasisStates(mixture);
				InitQubitsMap();
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

			void Trim() override
			{
				impl.Trim();
			}

			void ReCanonicalize() override
			{
				impl.ReCanonicalize();
			}

			MatrixClass getDensityMatrix() const override
			{
				const MatrixClass rhoInternal = impl.getDensityMatrix();
				const IndexType n = rhoInternal.rows();

				// the logical qubits are actually in some other physical positions,
				// the correspondence is in qubitsMap; translate the basis state indices
				std::vector<IndexType> mapState(n);
				for (IndexType s = 0; s < n; ++s)
				{
					size_t tmp = static_cast<size_t>(s);
					size_t mapped = 0;
					for (IndexType i = 0; i < static_cast<IndexType>(getNrQubits()); ++i)
					{
						if (tmp & 1ULL)
							mapped |= (1ULL << qubitsMap[i]);
						tmp >>= 1;
					}
					mapState[s] = static_cast<IndexType>(mapped);
				}

				MatrixClass rho(n, n);
				for (IndexType r = 0; r < n; ++r)
					for (IndexType c = 0; c < n; ++c)
						rho(r, c) = rhoInternal(mapState[r], mapState[c]);

				return rho;
			}

			void print() const override
			{
				impl.print();

				std::cout << "Qubits map: ";
				for (int q = 0; q < static_cast<int>(qubitsMap.size()); ++q)
					std::cout << q << "->" << qubitsMap[q] << " ";
				std::cout << std::endl;
			}

			// remap the logical Pauli string into the physical (impl) qubit ordering, then delegate
			std::complex<double> ExpectationValue(const std::string& pauliString) const override
			{
				const size_t nrQubits = getNrQubits();
				if (pauliString.size() != nrQubits)
					throw std::invalid_argument("Pauli string length must match the number of qubits");

				std::string mapped(nrQubits, 'I');
				for (size_t i = 0; i < nrQubits; ++i)
					mapped[static_cast<size_t>(qubitsMap[i])] = pauliString[i];

				return impl.ExpectationValue(mapped);
			}

			void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) override
			{
				ApplyOperator(gate);
			}

			void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				ApplyOperator(gate, qubit, controllingQubit1);
			}

			void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) override
			{
				for (const auto& gate : gates)
					ApplyGate(gate);
			}

			void ApplyOperator(const Gates::AppliedGate<MatrixClass>& op) override
			{
				ApplyOperator(op, op.getQubit1(), op.getQubit2());
			}

			void ApplyOperator(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");
				else if (op.getQubitsNumber() > 1 && (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(impl.getNrQubits())))
					throw std::invalid_argument("Qubit index out of bounds");

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = op.getQubitsNumber() > 1 ? qubitsMap[controllingQubit1] : qubit1;

				// for two qubit operators: if the qubits are not adjacent, swap them until they are
				if (op.getQubitsNumber() > 1 && std::abs(qubit1 - qubit2) > 1)
				{
					SwapQubits(qubit, controllingQubit1);

					qubit1 = qubitsMap[qubit];
					qubit2 = qubitsMap[controllingQubit1];
					assert(std::abs(qubit1 - qubit2) == 1);
				}

				impl.ApplyOperator(op, qubit1, qubit2);
			}

			void ApplyOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) override
			{
				for (const auto& op : ops)
					ApplyOperator(op);
			}

			void ApplyOperatorAndNormalize(const Gates::AppliedGate<MatrixClass>& op) override
			{
				ApplyOperatorAndNormalize(op, op.getQubit1(), op.getQubit2());
			}

			void ApplyOperatorAndNormalize(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");
				else if (op.getQubitsNumber() > 1 && (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(impl.getNrQubits())))
					throw std::invalid_argument("Qubit index out of bounds");

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = op.getQubitsNumber() > 1 ? qubitsMap[controllingQubit1] : qubit1;

				if (op.getQubitsNumber() > 1 && std::abs(qubit1 - qubit2) > 1)
				{
					SwapQubits(qubit, controllingQubit1);

					qubit1 = qubitsMap[qubit];
					qubit2 = qubitsMap[controllingQubit1];
					assert(std::abs(qubit1 - qubit2) == 1);
				}

				impl.ApplyOperatorAndNormalize(op, qubit1, qubit2);
			}

			void ApplyKrausOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) override
			{
				if (ops.empty()) return;

				const size_t qubit = ops.front().getQubit1();
				const size_t controllingQubit1 = ops.front().getQubit2();
				for (const auto& op : ops)
					if (op.getQubit1() != qubit || op.getQubit2() != controllingQubit1)
						throw std::invalid_argument("All Kraus operators need to act on the same qubits");

				ApplyKrausOperatorsImpl(ops, qubit, controllingQubit1);
			}

			void ApplyKrausOperators(const std::vector<MatrixClass>& ops, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				ApplyKrausOperatorsImpl(ops, qubit, controllingQubit1);
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

			double GetProbability(IndexType qubit, bool zeroVal = true) const override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");

				return impl.GetProbability(qubitsMap.at(qubit), zeroVal);
			}

			std::complex<double> getBasisStateMatrixElement(size_t row, size_t col) const override
			{
				const size_t nrQubits = getNrQubits();

				std::vector<bool> rowState(nrQubits, false);
				std::vector<bool> colState(nrQubits, false);

				for (size_t i = 0; i < nrQubits; ++i)
				{
					rowState[i] = (row & 1) == 1;
					colState[i] = (col & 1) == 1;
					row >>= 1;
					col >>= 1;
				}

				return getBasisStateMatrixElement(rowState, colState);
			}

			std::complex<double> getBasisStateMatrixElement(const std::vector<bool>& row, const std::vector<bool>& col) const override
			{
				const size_t nrQubits = getNrQubits();

				std::vector<bool> rowState(nrQubits, false);
				std::vector<bool> colState(nrQubits, false);

				for (size_t q = 0; q < row.size(); ++q)
					rowState[qubitsMap.at(q)] = row[q];
				for (size_t q = 0; q < col.size(); ++q)
					colState[qubitsMap.at(q)] = col[q];

				return impl.getBasisStateMatrixElement(rowState, colState);
			}

			double getBasisStateProbability(size_t State) const override
			{
				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return 0.;

				return ClampProbability((getBasisStateMatrixElement(State, State) / tr).real());
			}

			double getBasisStateProbability(const std::vector<bool>& State) const override
			{
				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return 0.;

				return ClampProbability((getBasisStateMatrixElement(State, State) / tr).real());
			}

			std::complex<double> Trace() const override
			{
				return impl.Trace();
			}

			std::shared_ptr<MPOSimulatorStateInterface> getState() const override
			{
				auto baseState = std::static_pointer_cast<MPOSimulatorBaseState>(impl.getState());
				if (!baseState) return nullptr;

				auto state = std::make_shared<MPOSimulatorState>();
				state->gammas.swap(baseState->gammas);
				state->lambdas.swap(baseState->lambdas);
				state->qubitsMap = qubitsMap;
				state->qubitsMapInv = qubitsMapInv;
				return state;
			}

			void setState(const std::shared_ptr<MPOSimulatorStateInterface>& state) override
			{
				if (!state) return;

				auto baseState = std::static_pointer_cast<MPOSimulatorBaseState>(state);
				impl.setState(baseState);

				auto simState = std::static_pointer_cast<MPOSimulatorState>(state);
				qubitsMap = simState->qubitsMap;
				qubitsMapInv = simState->qubitsMapInv;
			}

			std::vector<IndexType> getBondDimensions() const
			{
				return impl.getBondDimensions();
			}

		private:
			template<class OperatorsContainer> void ApplyKrausOperatorsImpl(const OperatorsContainer& ops, IndexType qubit, IndexType controllingQubit1)
			{
				if (ops.empty()) return;

				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");
				else if (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");

				const size_t qubitsNumber = KrausOperatorQubitsNumber(ops.front());
				for (const auto& op : ops)
					if (KrausOperatorQubitsNumber(op) != qubitsNumber)
						throw std::invalid_argument("All Kraus operators need to have the same number of qubits");

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = qubitsMap[controllingQubit1];

				if (qubitsNumber > 1 && std::abs(qubit1 - qubit2) > 1)
				{
					SwapQubits(qubit, controllingQubit1);

					qubit1 = qubitsMap[qubit];
					qubit2 = qubitsMap[controllingQubit1];
					assert(std::abs(qubit1 - qubit2) == 1);
				}

				std::vector<MatrixClass> mappedOps;
				mappedOps.reserve(ops.size());
				for (const auto& op : ops)
					mappedOps.push_back(KrausOperatorMatrix(op));

				impl.ApplyKrausOperators(mappedOps, qubit1, qubit2);
			}

			static size_t KrausOperatorQubitsNumber(const Gates::AppliedGate<MatrixClass>& op)
			{
				return op.getQubitsNumber();
			}

			static size_t KrausOperatorQubitsNumber(const MatrixClass& op)
			{
				return static_cast<size_t>(std::log2(op.rows()));
			}

			static const MatrixClass& KrausOperatorMatrix(const Gates::AppliedGate<MatrixClass>& op)
			{
				return op.getRawOperatorMatrix();
			}

			static const MatrixClass& KrausOperatorMatrix(const MatrixClass& op)
			{
				return op;
			}

			static double ClampProbability(double probability)
			{
				constexpr double tolerance = 1E-12;
				if (probability < 0. && probability > -tolerance) return 0.;
				if (probability > 1. && probability < 1. + tolerance) return 1.;

				return probability;
			}

			void InitQubitsMap()
			{
				qubitsMap.resize(getNrQubits());
				qubitsMapInv.resize(getNrQubits());

				for (IndexType i = 0; i < static_cast<IndexType>(getNrQubits()); ++i)
					qubitsMapInv[i] = qubitsMap[i] = i;
			}

			// brings the two logical qubits to adjacent physical positions using nearest neighbour swaps
			void SwapQubits(IndexType logical1, IndexType logical2)
			{
				IndexType r1 = qubitsMap[logical1];
				IndexType r2 = qubitsMap[logical2];
				if (r1 > r2)
					std::swap(r1, r2);

				// move the qubit at physical position r2 down until it sits right next to r1
				while (r2 > r1 + 1)
				{
					const IndexType from = r2;
					const IndexType to = r2 - 1;

					impl.ApplyGate(swapGate, from, to);

					const IndexType lFrom = qubitsMapInv[from];
					const IndexType lTo = qubitsMapInv[to];

					std::swap(qubitsMap[lFrom], qubitsMap[lTo]);
					qubitsMapInv[from] = lTo;
					qubitsMapInv[to] = lFrom;

					r2 = to;
				}
			}

			MPOSimulatorImpl impl;
			std::vector<IndexType> qubitsMap;
			std::vector<IndexType> qubitsMapInv;
			QC::Gates::SwapGate<MatrixClass> swapGate;
		};

	}

}
