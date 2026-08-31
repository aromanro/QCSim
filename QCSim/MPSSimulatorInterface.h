#pragma once

#include <set>

namespace QC {

	namespace TensorNetworks {

		class MPSSimulatorStateInterface
		{
		public:
			virtual ~MPSSimulatorStateInterface() = default;
		};

		class MPSSimulatorInterface
		{
		public:
			using LambdaType = Eigen::VectorXd;
			using GammaType = Eigen::Tensor<std::complex<double>, 3>;
			using MatrixTensorType = Eigen::Tensor<std::complex<double>, 2>;
			using MatrixClass = Eigen::MatrixXcd;
			using VectorClass = Eigen::VectorXcd;
			using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
			using IndexType = Eigen::Index;
			using IntIndexPair = Eigen::IndexPair<int>;
			using Indexes = Eigen::array<IntIndexPair, 1>;
			using OneQubitGateTensor = Eigen::TensorFixedSize<std::complex<double>, Eigen::Sizes<2, 2>>;
			using TwoQubitsGateTensor = Eigen::TensorFixedSize<std::complex<double>, Eigen::Sizes<2, 2, 2, 2>>;

			// Selects how setLimitEntanglement's threshold is interpreted when truncating a bond's
			// singular values:
			//  - RelativeToMax: keep sigma_i iff it is strictly greater than threshold * sigma_max
			//    (this simulator's original behavior, matching Eigen::SVDBase::rank()'s own
			//    relative-cutoff convention).
			//  - DiscardedWeight: discard the smallest singular values while the cumulative sum of
			//    their squares, normalized by the sum of squares of all of them, stays below
			//    threshold. This is Qiskit Aer's and ITensor's convention (Aer's reduce_zeros in
			//    matrix_product_state/svd.cpp; ITensor's default `cutoff`; NVIDIA cuTensorNet's
			//    CUTENSORNET_TENSOR_SVD_CONFIG_DISCARDED_WEIGHT_CUTOFF).
			//
			// DiscardedWeight is the default (see MPSSimulatorBase) - this is a deliberate default
			// change from earlier versions of this simulator, which only ever implemented
			// RelativeToMax. Pass RelativeToMax explicitly to reproduce that old behavior exactly.
			enum class TruncationMode
			{
				RelativeToMax,
				DiscardedWeight
			};

			MPSSimulatorInterface() = default;
			virtual ~MPSSimulatorInterface() = default;

			virtual size_t getNrQubits() const = 0;
			virtual void Clear() = 0;
			virtual void InitOnesState() = 0;
			virtual void setToQubitState(IndexType q) = 0;
			virtual void setToBasisState(size_t State) = 0;
			// Bit i selects qubit i. This overload also represents basis states wider
			// than size_t; shorter vectors are zero-extended and oversized vectors are rejected.
			virtual void setToBasisState(const std::vector<bool>& State) = 0;
			virtual void setLimitBondDimension(IndexType chival) = 0;
			virtual void setLimitEntanglement(double svdThreshold) = 0;
			virtual void dontLimitBondDimension() = 0;
			virtual void dontLimitEntanglement() = 0;
			// Returns true once the mode is applied. Both TruncationMode values are always
			// supported here (unlike e.g. the Aer backend, which only supports DiscardedWeight);
			// an out-of-range enum value throws std::invalid_argument instead of returning false.
			// The bool return is kept for API symmetry with backends where "well-formed but
			// unsupported" is a real, reachable outcome.
			virtual bool setTruncationMode(TruncationMode mode) = 0;
			virtual TruncationMode getTruncationMode() const = 0;
			virtual void Trim() = 0;
			virtual void ReCanonicalize() = 0;
			virtual VectorClass getRegisterStorage() const = 0;
			virtual void print() const = 0;
			virtual void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) = 0;
			virtual void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) = 0;
			virtual void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) = 0;
			virtual bool MeasureQubit(IndexType qubit) = 0;
			virtual std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) = 0;
			virtual std::unordered_map<IndexType, bool> MeasureNoCollapse() = 0;
			virtual double GetProbability(IndexType qubit, bool zeroVal = true) const = 0;
			virtual std::complex<double> getBasisStateAmplitude(size_t State) const = 0;
			virtual std::complex<double> getBasisStateAmplitude(std::vector<bool>& State) const = 0;
			virtual double getBasisStateProbability(size_t State) const = 0;
			virtual double getBasisStateProbability(std::vector<bool>& State) const = 0;
			virtual std::shared_ptr<MPSSimulatorStateInterface> getState() const = 0;
			virtual void setState(const std::shared_ptr<MPSSimulatorStateInterface>& state) = 0;
			virtual void setStateDestructive(std::shared_ptr<MPSSimulatorStateInterface>& state) = 0;
			virtual void MoveAtBeginningOfChain(const std::set<IndexType>& qubits) = 0;
			virtual std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) = 0;
			virtual std::complex<double> ProjectOnZero() const = 0;
			virtual std::unordered_map<IndexType, bool> MeasureNoCollapse(const std::set<IndexType>& qubits) = 0;
		};

	}

}

