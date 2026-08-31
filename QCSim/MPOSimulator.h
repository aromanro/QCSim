#pragma once

#include "MPOSimulatorImpl.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>
#include <unordered_set>

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
			// Callback signature: (bondDims) -> meetPosition
			using MeetingPositionCallback = std::function<IndexType(
				const std::vector<IndexType>&)>;

			using BondDimensionCallback = std::function<void(const std::vector<IndexType>&)>;

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

			void SetInitialQubitsMap(const std::vector<long long int>& initialMap)
			{
				const size_t nrQubits = getNrQubits();
				if (initialMap.size() != nrQubits)
					throw std::invalid_argument("Initial qubit-map size must match the number of qubits");

				std::vector<IndexType> newQubitsMap(nrQubits);
				std::vector<IndexType> newQubitsMapInv(nrQubits);
				std::vector<bool> seenPhysical(nrQubits, false);
				for (size_t logical = 0; logical < nrQubits; ++logical)
				{
					const long long int physical = initialMap[logical];
					if (physical < 0 || static_cast<unsigned long long>(physical) >= nrQubits ||
						seenPhysical[static_cast<size_t>(physical)])
						throw std::invalid_argument("Initial qubit map must be a permutation");

					seenPhysical[static_cast<size_t>(physical)] = true;
					newQubitsMap[logical] = static_cast<IndexType>(physical);
					newQubitsMapInv[static_cast<size_t>(physical)] = static_cast<IndexType>(logical);
				}

				qubitsMap.swap(newQubitsMap);
				qubitsMapInv.swap(newQubitsMapInv);
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

			bool setTruncationMode(TruncationMode mode) override
			{
				return impl.setTruncationMode(mode);
			}

			TruncationMode getTruncationMode() const override
			{
				return impl.getTruncationMode();
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
				const size_t qubitsNumber = ValidateOperator(op);
				const IndexType qubit = CheckedAppliedQubit(op.getQubit1());
				const IndexType controllingQubit1 = qubitsNumber > 1 ? CheckedAppliedQubit(op.getQubit2()) : 0;

				ApplyValidatedOperator(op, qubit, controllingQubit1, qubitsNumber);
			}

			void ApplyOperator(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				const size_t qubitsNumber = ValidateOperator(op);
				ApplyValidatedOperator(op, qubit, controllingQubit1, qubitsNumber);
			}

			void ApplyOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) override
			{
				for (const auto& op : ops)
					ApplyOperator(op);
			}

			void ApplyOperatorAndNormalize(const Gates::AppliedGate<MatrixClass>& op) override
			{
				const size_t qubitsNumber = ValidateOperator(op);
				const IndexType qubit = CheckedAppliedQubit(op.getQubit1());
				const IndexType controllingQubit1 = qubitsNumber > 1 ? CheckedAppliedQubit(op.getQubit2()) : 0;

				ApplyValidatedOperatorAndNormalize(op, qubit, controllingQubit1, qubitsNumber);
			}

			void ApplyOperatorAndNormalize(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				const size_t qubitsNumber = ValidateOperator(op);
				ApplyValidatedOperatorAndNormalize(op, qubit, controllingQubit1, qubitsNumber);
			}

			void ApplyKrausOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) override
			{
				if (ops.empty()) return;

				const size_t qubitsNumber = ValidateOperator(ops.front());
				const size_t rawQubit = ops.front().getQubit1();
				const size_t rawControllingQubit1 = ops.front().getQubit2();
				for (const auto& op : ops)
				{
					if (ValidateOperator(op) != qubitsNumber)
						throw std::invalid_argument("All Kraus operators need to have the same number of qubits");
					if (op.getQubit1() != rawQubit ||
						(qubitsNumber > 1 && op.getQubit2() != rawControllingQubit1))
						throw std::invalid_argument("All Kraus operators need to act on the same qubits");
				}

				const IndexType qubit = CheckedAppliedQubit(rawQubit);
				const IndexType controllingQubit1 = qubitsNumber > 1 ? CheckedAppliedQubit(rawControllingQubit1) : 0;
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

				const bool result = impl.MeasureQubit(qubitsMap[qubit]);
				if (bondDimensionCallback)
					bondDimensionCallback(impl.getBondDimensions());
				return result;
			}

			std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) override
			{
				ValidateQubitSet(qubits);

				std::set<IndexType> mappedQubits;
				for (const auto qubit : qubits)
					mappedQubits.insert(qubitsMap[qubit]);

				const auto measuredQubits = impl.MeasureQubits(mappedQubits);

				std::unordered_map<IndexType, bool> res;
				for (const auto& [qubit, val] : measuredQubits)
					res[qubitsMapInv[qubit]] = val;

				if (bondDimensionCallback)
					bondDimensionCallback(impl.getBondDimensions());

				return res;
			}

			std::unordered_map<IndexType, bool> MeasureNoCollapse() override
			{
				const auto measuredQubits = impl.MeasureNoCollapse();

				std::unordered_map<IndexType, bool> res;
				for (const auto& [qubit, val] : measuredQubits)
					res[qubitsMapInv[qubit]] = val;

				return res;
			}

			std::unordered_map<IndexType, bool> MeasureNoCollapse(const std::set<IndexType>& qubits) override
			{
				if (qubits.empty()) return {};
				ValidateQubitSet(qubits);

				std::set<IndexType> mappedQubits;
				for (const auto qubit : qubits)
					mappedQubits.insert(qubitsMap[qubit]);

				const auto measuredQubits = impl.MeasureNoCollapse(mappedQubits);

				std::unordered_map<IndexType, bool> res;
				for (const auto& [qubit, val] : measuredQubits)
					res[qubitsMapInv[qubit]] = val;

				return res;
			}

			// moves the given qubits at the beginning of the chain (helps with sampling a subset of qubits faster)
			void MoveAtBeginningOfChain(const std::set<IndexType>& qubits) override
			{
				ValidateQubitSet(qubits);

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

					assert(std::abs(movingQubitReal - toQubitReal) == 1);

					impl.ApplyGate(swapGate, movingQubitReal, toQubitReal);

					// update the maps
					qubitsMap[logicalQubit] = toQubitReal;
					qubitsMap[currentLogicalPosQubit] = movingQubitReal;
					qubitsMapInv[toQubitReal] = logicalQubit;
					qubitsMapInv[movingQubitReal] = currentLogicalPosQubit;

					if (bondDimensionCallback)
						bondDimensionCallback(impl.getBondDimensions());

					handledQubits.insert(logicalQubit);
					if (handledQubits.size() == qubits.size())
						break;

					++currentQubitPos;
				}
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
					row >>= 1;
				}

				for (size_t i = 0; i < nrQubits; ++i)
				{
					colState[i] = (col & 1) == 1;
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
				auto baseState = std::dynamic_pointer_cast<MPOSimulatorBaseState>(impl.getState());
				if (!baseState)
					throw std::logic_error("MPO implementation returned an incompatible state type");

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

				const auto simState = std::dynamic_pointer_cast<const MPOSimulatorState>(state);
				if (!simState)
					throw std::invalid_argument("State is not compatible with MPOSimulator");
				ValidateStateStructure(*simState);

				// Complete every allocation before changing either the MPO or its logical mapping.
				auto implState = MakeImplStateCopy(*simState);
				auto newQubitsMap = simState->qubitsMap;
				auto newQubitsMapInv = simState->qubitsMapInv;
				std::shared_ptr<MPOSimulatorStateInterface> implStateInterface = implState;

				// The temporary has the exact base-state type expected by MPOSimulatorImpl.
				// Destructive restore avoids a second tensor copy while leaving the caller's
				// non-destructive state untouched.
				impl.setStateDestructive(implStateInterface);
				qubitsMap.swap(newQubitsMap);
				qubitsMapInv.swap(newQubitsMapInv);
			}

			void setStateDestructive(std::shared_ptr<MPOSimulatorStateInterface>& state) override
			{
				if (!state) return;

				auto simState = std::dynamic_pointer_cast<MPOSimulatorState>(state);
				if (!simState)
					throw std::invalid_argument("State is not compatible with MPOSimulator");
				ValidateStateStructure(*simState);

				// Allocate the exact base-state wrapper before moving anything out of the
				// caller's state. From this point onward only noexcept vector swaps occur.
				auto implState = std::make_shared<MPOSimulatorBaseState>();
				implState->gammas.swap(simState->gammas);
				implState->lambdas.swap(simState->lambdas);
				std::shared_ptr<MPOSimulatorStateInterface> implStateInterface = implState;

				try
				{
					impl.setStateDestructive(implStateInterface);
				}
				catch (...)
				{
					// A conforming implementation validates before mutation. Restore the source
					// as well so a rejected destructive restore has strong exception safety.
					simState->gammas.swap(implState->gammas);
					simState->lambdas.swap(implState->lambdas);
					throw;
				}

				qubitsMap.swap(simState->qubitsMap);
				qubitsMapInv.swap(simState->qubitsMapInv);
				// Preserve the previous destructive semantics for any outstanding aliases:
				// they observe the simulator state which was replaced.
				simState->gammas.swap(implState->gammas);
				simState->lambdas.swap(implState->lambdas);
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
				setStateDestructive(savedState);
			}

			std::unique_ptr<MPOSimulator> Clone() const
			{
				auto sim = std::make_unique<MPOSimulator>(getNrQubits());

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
				sim->bondDimensionCallback = bondDimensionCallback;

				if (savedState)
				{
					const auto state = std::dynamic_pointer_cast<const MPOSimulatorState>(savedState);
					if (!state)
						throw std::logic_error("Saved MPO state has an incompatible type");
					ValidateStateStructure(*state);

					auto stateClone = std::make_shared<MPOSimulatorState>();
					stateClone->gammas = state->gammas;
					stateClone->lambdas = state->lambdas;
					stateClone->qubitsMap = state->qubitsMap;
					stateClone->qubitsMapInv = state->qubitsMapInv;
					sim->savedState = stateClone;
				}

				return sim;
			}

			std::vector<IndexType> getBondDimensions() const
			{
				return impl.getBondDimensions();
			}

			// Enable/disable immediate meeting position optimization
			// using actual bond dimensions (no lookahead, just cheapest path)
			void SetUseOptimalMeetingPosition(bool enable)
			{
				useOptimalMeetingPosition = enable;
			}

			// Set a callback for external lookahead-based meeting position.
			// When set, this takes priority over both the heuristic and the
			// local optimizer. Pass nullptr to clear.
			void SetMeetingPositionCallback(MeetingPositionCallback callback)
			{
				meetingPositionCallback = std::move(callback);
			}

			// Set a callback that receives the current bond dimensions after an
			// operation that can change them. Pass nullptr to clear.
			void SetBondDimensionCallback(BondDimensionCallback callback)
			{
				bondDimensionCallback = std::move(callback);
			}

		private:
			static size_t ValidateOperatorMatrix(const MatrixClass& op)
			{
				if (op.rows() != op.cols() || (op.rows() != 2 && op.rows() != 4))
					throw std::invalid_argument("MPO operators must be finite 2x2 or 4x4 matrices");
				if (!op.allFinite())
					throw std::invalid_argument("MPO operator contains a non-finite value");

				return op.rows() == 2 ? 1 : 2;
			}

			static size_t ValidateOperator(const GateClass& op)
			{
				const size_t matrixQubits = ValidateOperatorMatrix(op.getRawOperatorMatrix());
				if (op.getQubitsNumber() != matrixQubits)
					throw std::invalid_argument("Operator matrix dimensions do not match its declared arity");
				return matrixQubits;
			}

			static IndexType CheckedAppliedQubit(size_t qubit)
			{
				if (qubit > static_cast<size_t>(std::numeric_limits<IndexType>::max()))
					throw std::invalid_argument("Qubit index out of bounds");
				return static_cast<IndexType>(qubit);
			}

			void ValidateLogicalQubit(IndexType qubit) const
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(impl.getNrQubits()))
					throw std::invalid_argument("Qubit index out of bounds");
			}

			void ValidateOperatorQubits(size_t qubitsNumber, IndexType qubit, IndexType controllingQubit1) const
			{
				ValidateLogicalQubit(qubit);
				if (qubitsNumber > 1)
				{
					ValidateLogicalQubit(controllingQubit1);
					if (qubit == controllingQubit1)
						throw std::invalid_argument("Two-qubit operators require distinct qubits");
				}
			}

			void ValidateQubitSet(const std::set<IndexType>& qubits) const
			{
				// Validate the complete set before the caller reads or changes either map.
				for (const IndexType qubit : qubits)
					ValidateLogicalQubit(qubit);
			}

			std::pair<IndexType, IndexType> RouteOperatorQubits(
				IndexType qubit, IndexType controllingQubit1, size_t qubitsNumber)
			{
				ValidateOperatorQubits(qubitsNumber, qubit, controllingQubit1);

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = qubitsNumber > 1 ? qubitsMap[controllingQubit1] : qubit1;

				if (qubitsNumber > 1 && std::abs(qubit1 - qubit2) > 1)
				{
					BringQubitsTogether(qubit, controllingQubit1);
					qubit1 = qubitsMap[qubit];
					qubit2 = qubitsMap[controllingQubit1];
					assert(std::abs(qubit1 - qubit2) == 1);
				}

				return { qubit1, qubit2 };
			}

			void ApplyValidatedOperator(const GateClass& op, IndexType qubit,
				IndexType controllingQubit1, size_t qubitsNumber)
			{
				const auto [qubit1, qubit2] = RouteOperatorQubits(qubit, controllingQubit1, qubitsNumber);
				impl.ApplyOperator(op, qubit1, qubit2);
				if (bondDimensionCallback && qubitsNumber > 1)
					bondDimensionCallback(impl.getBondDimensions());
			}

			void ApplyValidatedOperatorAndNormalize(const GateClass& op, IndexType qubit,
				IndexType controllingQubit1, size_t qubitsNumber)
			{
				const auto [qubit1, qubit2] = RouteOperatorQubits(qubit, controllingQubit1, qubitsNumber);
				impl.ApplyOperatorAndNormalize(op, qubit1, qubit2);
				if (bondDimensionCallback && qubitsNumber > 1)
					bondDimensionCallback(impl.getBondDimensions());
			}

			template<class OperatorsContainer> void ApplyKrausOperatorsImpl(const OperatorsContainer& ops, IndexType qubit, IndexType controllingQubit1)
			{
				if (ops.empty()) return;

				// Validate every matrix and every relevant target before routing performs swaps.
				const size_t qubitsNumber = ValidateOperatorMatrix(KrausOperatorMatrix(ops.front()));
				for (const auto& op : ops)
					if (ValidateOperatorMatrix(KrausOperatorMatrix(op)) != qubitsNumber)
						throw std::invalid_argument("All Kraus operators need to have the same number of qubits");
				ValidateOperatorQubits(qubitsNumber, qubit, controllingQubit1);

				IndexType qubit1 = qubitsMap[qubit];
				IndexType qubit2 = qubitsNumber > 1 ? qubitsMap[controllingQubit1] : qubit1;

				if (qubitsNumber > 1 && std::abs(qubit1 - qubit2) > 1)
				{
					BringQubitsTogether(qubit, controllingQubit1);

					qubit1 = qubitsMap[qubit];
					qubit2 = qubitsMap[controllingQubit1];
					assert(std::abs(qubit1 - qubit2) == 1);
				}

				std::vector<MatrixClass> mappedOps;
				mappedOps.reserve(ops.size());
				for (const auto& op : ops)
					mappedOps.push_back(KrausOperatorMatrix(op));

				impl.ApplyKrausOperators(mappedOps, qubit1, qubit2);
				if (bondDimensionCallback && qubitsNumber > 1)
					bondDimensionCallback(impl.getBondDimensions());
			}

			static const MatrixClass& KrausOperatorMatrix(const Gates::AppliedGate<MatrixClass>& op)
			{
				return op.getRawOperatorMatrix();
			}

			static const MatrixClass& KrausOperatorMatrix(const MatrixClass& op)
			{
				return op;
			}

			void ValidateStateStructure(const MPOSimulatorState& state) const
			{
				const size_t nrQubits = getNrQubits();
				const size_t nrBonds = nrQubits > 0 ? nrQubits - 1 : 0;

				if (state.gammas.size() != nrQubits || state.lambdas.size() != nrBonds)
					throw std::invalid_argument("MPO state qubit or bond count does not match the simulator");
				if (state.qubitsMap.size() != nrQubits || state.qubitsMapInv.size() != nrQubits)
					throw std::invalid_argument("MPO state qubit-map size does not match the simulator");

				for (size_t q = 0; q < nrQubits; ++q)
				{
					const auto& gamma = state.gammas[q];
					if (gamma.dimension(0) <= 0 || gamma.dimension(1) != 2 ||
						gamma.dimension(2) != 2 || gamma.dimension(3) <= 0)
						throw std::invalid_argument("MPO state contains an invalid site tensor");
					if ((q == 0 && gamma.dimension(0) != 1) ||
						(q + 1 == nrQubits && gamma.dimension(3) != 1))
						throw std::invalid_argument("MPO state has invalid boundary bond dimensions");

					for (IndexType i = 0; i < gamma.size(); ++i)
					{
						const auto value = gamma.data()[i];
						if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
							throw std::invalid_argument("MPO state contains a non-finite tensor value");
					}
				}

				for (size_t bond = 0; bond < nrBonds; ++bond)
				{
					const IndexType dimension = state.lambdas[bond].size();
					if (dimension <= 0 || state.gammas[bond].dimension(3) != dimension ||
						state.gammas[bond + 1].dimension(0) != dimension ||
						!state.lambdas[bond].allFinite() || (state.lambdas[bond].array() < 0.).any())
						throw std::invalid_argument("MPO state has inconsistent bond dimensions");
				}

				std::vector<bool> seenPhysical(nrQubits, false);
				std::vector<bool> seenLogical(nrQubits, false);
				for (size_t logical = 0; logical < nrQubits; ++logical)
				{
					const IndexType physical = state.qubitsMap[logical];
					if (physical < 0 || static_cast<size_t>(physical) >= nrQubits || seenPhysical[physical])
						throw std::invalid_argument("MPO state qubit map is not a permutation");
					seenPhysical[physical] = true;
				}

				for (size_t physical = 0; physical < nrQubits; ++physical)
				{
					const IndexType logical = state.qubitsMapInv[physical];
					if (logical < 0 || static_cast<size_t>(logical) >= nrQubits || seenLogical[logical])
						throw std::invalid_argument("MPO state inverse qubit map is not a permutation");
					seenLogical[logical] = true;
					if (state.qubitsMap[logical] != static_cast<IndexType>(physical))
						throw std::invalid_argument("MPO state qubit maps are not inverses");
				}
			}

			static std::shared_ptr<MPOSimulatorBaseState> MakeImplStateCopy(const MPOSimulatorState& state)
			{
				auto implState = std::make_shared<MPOSimulatorBaseState>();
				implState->gammas = state.gammas;
				implState->lambdas = state.lambdas;
				return implState;
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

			void BringQubitsTogether(IndexType logical1, IndexType logical2)
			{
				if (meetingPositionCallback)
				{
					const IndexType meetPosition = meetingPositionCallback(impl.getBondDimensions());
					SwapQubitsToPosition(logical1, logical2, meetPosition);
				}
				else if (useOptimalMeetingPosition)
				{
					const IndexType meetPosition = FindBestMeetingPositionLocal(logical1, logical2);
					SwapQubitsToPosition(logical1, logical2, meetPosition);
				}
				else
					SwapQubits(logical1, logical2);
			}

			// Brings the two logical qubits to adjacent physical positions using nearest
			// neighbour swaps. The default heuristic prefers the chain middle when the
			// qubits straddle it, otherwise it moves the qubit that is closer to an end.
			void SwapQubits(IndexType qubit1, IndexType qubit2, bool forceSwapDown = false)
			{
				IndexType realq1 = qubitsMap[qubit1];
				IndexType realq2 = qubitsMap[qubit2];
				if (realq1 > realq2)
				{
					std::swap(realq1, realq2);
					std::swap(qubit1, qubit2);
				}

				if (realq2 - realq1 <= 1) return;

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

					if (bondDimensionCallback)
						bondDimensionCallback(impl.getBondDimensions());

					movingQubitReal = toQubitReal;
				} while (movingQubitReal != targetQubitReal);

				assert(std::abs(qubitsMap[qubit1] - qubitsMap[qubit2]) == 1);
			}

			// Swap two logical qubits so they meet at a specified bond position.
			// meetPosition is expressed in physical chain coordinates; the qubits
			// finish at meetPosition and meetPosition + 1.
			void SwapQubitsToPosition(IndexType logical1, IndexType logical2, IndexType meetPosition)
			{
				IndexType r1 = qubitsMap[logical1];
				IndexType r2 = qubitsMap[logical2];
				if (r1 > r2)
				{
					std::swap(r1, r2);
					std::swap(logical1, logical2);
				}

				if (r2 - r1 <= 1) return;

				// Fall back to the local optimizer or the existing heuristic if the
				// callback returns a bond outside the interval between the two qubits.
				if (meetPosition < r1 || meetPosition >= r2)
				{
					if (useOptimalMeetingPosition)
						meetPosition = FindBestMeetingPositionLocal(logical1, logical2);
					else
					{
						SwapQubits(logical1, logical2);
						return;
					}
				}

				while (r1 < meetPosition)
				{
					const IndexType to = r1 + 1;
					const IndexType displacedLogical = qubitsMapInv[to];

					impl.ApplyGate(swapGate, r1, to);

					qubitsMap[displacedLogical] = r1;
					qubitsMapInv[r1] = displacedLogical;
					qubitsMap[logical1] = to;
					qubitsMapInv[to] = logical1;

					if (bondDimensionCallback)
						bondDimensionCallback(impl.getBondDimensions());
					r1 = to;
				}

				while (r2 > meetPosition + 1)
				{
					const IndexType to = r2 - 1;
					const IndexType displacedLogical = qubitsMapInv[to];

					impl.ApplyGate(swapGate, to, r2);

					qubitsMap[displacedLogical] = r2;
					qubitsMapInv[r2] = displacedLogical;
					qubitsMap[logical2] = to;
					qubitsMapInv[to] = logical2;

					if (bondDimensionCallback)
						bondDimensionCallback(impl.getBondDimensions());
					r2 = to;
				}

				assert(std::abs(qubitsMap[logical1] - qubitsMap[logical2]) == 1);
			}

			// Heuristic: swap towards positions in chain with smaller bond dimensions,
			// the cost of a generic gate is bigger than the cost of a swap
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

			MPOSimulatorImpl impl;
			std::vector<IndexType> qubitsMap;
			std::vector<IndexType> qubitsMapInv;
			QC::Gates::SwapGate<MatrixClass> swapGate;

			bool useOptimalMeetingPosition = true;

			MeetingPositionCallback meetingPositionCallback;
			BondDimensionCallback bondDimensionCallback;

			std::shared_ptr<MPOSimulatorStateInterface> savedState;
		};

	}

}
