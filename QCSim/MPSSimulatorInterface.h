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

			MPSSimulatorInterface() = default;
			virtual ~MPSSimulatorInterface() = default;

			virtual size_t getNrQubits() const = 0;
			virtual void Clear() = 0;
			virtual void InitOnesState() = 0;
			virtual void setToQubitState(IndexType q) = 0;
			virtual void setToBasisState(size_t State) = 0;
			virtual void setToBasisState(const std::vector<bool>& State) = 0;
			virtual void setLimitBondDimension(IndexType chival) = 0;
			virtual void setLimitEntanglement(double svdThreshold) = 0;
			virtual void dontLimitBondDimension() = 0;
			virtual void dontLimitEntanglement() = 0;
			virtual VectorClass getRegisterStorage() const = 0;
			virtual void print() const = 0;
			virtual void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) = 0;
			virtual void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) = 0;
			virtual void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) = 0;
			virtual bool MeasureQubit(IndexType qubit) = 0;
			virtual std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) = 0;
			virtual std::vector<bool> MeasureNoCollapse() = 0;
			virtual double GetProbability(IndexType qubit, bool zeroVal = true) const = 0;
			virtual std::complex<double> getBasisStateAmplitude(size_t State) const = 0;
			virtual std::complex<double> getBasisStateAmplitude(std::vector<bool>& State) const = 0;
			virtual double getBasisStateProbability(size_t State) const = 0;
			virtual double getBasisStateProbability(std::vector<bool>& State) const = 0;
			virtual std::shared_ptr<MPSSimulatorStateInterface> getState() const = 0;
			virtual void setState(const std::shared_ptr<MPSSimulatorStateInterface>& state) = 0;
			virtual void MoveAtBeginningOfChain(const std::set<IndexType>& qubits) = 0;
		};

	}

}

