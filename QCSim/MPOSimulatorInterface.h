#pragma once

#include <set>
#include <cmath>
#include <vector>
#include <memory>
#include <complex>
#include <utility>
#include <unordered_map>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include "SimpleGates.h"

namespace QC {

	namespace TensorNetworks {

		// A Matrix Product Operator (MPO) simulator.
		//
		// While the MPS simulator (see MPSSimulatorInterface.h) is the tensor-network
		// analogue of a state vector |psi>, this one is the analogue of a density
		// matrix rho.
		//
		// In the MPS each qubit is described by a rank-3 tensor with two virtual
		// (bond) legs and one physical leg. Here, each qubit is described by a rank-4
		// tensor with two virtual (bond) legs and two physical legs: a 'ket' (row) leg
		// and a 'bra' (column) leg. The leg order used everywhere is:
		//
		//        (leftBond, ket, bra, rightBond)
		//
		// Applying a unitary gate U evolves the density matrix as rho -> U rho U^dagger,
		// so the gate tensor is contracted with the ket legs and its adjoint (the
		// conjugate gate tensor) is contracted with the bra legs.
		//
		// Compression caveat: with no SVD truncation this is an exact MPO representation
		// of the density matrix. When bond dimension or singular-value truncation is
		// enabled, the simulator becomes an approximate operator-space MPO simulator.
		// The compressed operator is not guaranteed to remain a physical density matrix:
		// trace may drift, Hermiticity may be approximate and positivity is not preserved
		// by ordinary MPO SVD truncation.

		class MPOSimulatorStateInterface
		{
		public:
			virtual ~MPOSimulatorStateInterface() = default;
		};

		class MPOSimulatorInterface
		{
		public:
			using LambdaType = Eigen::VectorXd;
			// rank-4 site tensor: (leftBond, ket, bra, rightBond)
			using TensorType = Eigen::Tensor<std::complex<double>, 4>;
			using MatrixTensorType = Eigen::Tensor<std::complex<double>, 2>;
			using MatrixClass = Eigen::MatrixXcd;
			using VectorClass = Eigen::VectorXcd;
			using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
			using IndexType = Eigen::Index;
			using IntIndexPair = Eigen::IndexPair<int>;
			using Indexes = Eigen::array<IntIndexPair, 1>;
			using Indexes2 = Eigen::array<IntIndexPair, 2>;
			using OneQubitGateTensor = Eigen::TensorFixedSize<std::complex<double>, Eigen::Sizes<2, 2>>;
			using TwoQubitsGateTensor = Eigen::TensorFixedSize<std::complex<double>, Eigen::Sizes<2, 2, 2, 2>>;

			MPOSimulatorInterface() = default;
			virtual ~MPOSimulatorInterface() = default;

			virtual size_t getNrQubits() const = 0;
			virtual void Clear() = 0;
			virtual void InitOnesState() = 0;
			virtual void setToQubitState(IndexType q) = 0;
			virtual void setToBasisState(size_t State) = 0;
			virtual void setToBasisState(const std::vector<bool>& State) = 0;

			// Unlike a state vector, a density matrix can represent a statistical
			// (classical) mixture of pure states. These set rho to a mixture of basis
			// states: rho = sum_i prob_i |state_i><state_i|. The probabilities should
			// be non negative; they are normalized so that Tr(rho) = 1.
			virtual void setToMixtureOfBasisStates(const std::vector<std::pair<size_t, double>>& mixture) = 0;
			virtual void setToMixtureOfBasisStates(const std::vector<std::pair<std::vector<bool>, double>>& mixture) = 0;

			virtual void setLimitBondDimension(IndexType chival) = 0;
			virtual void setLimitEntanglement(double svdThreshold) = 0;
			virtual void dontLimitBondDimension() = 0;
			virtual void dontLimitEntanglement() = 0;
			virtual void Trim() = 0;
			virtual void ReCanonicalize() = 0;

			// the analogue of MPS getRegisterStorage(): the full 2^N x 2^N density matrix
			virtual MatrixClass getDensityMatrix() const = 0;
			virtual void print() const = 0;

			virtual void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) = 0;
			virtual void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) = 0;
			virtual void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) = 0;

			// Explicit non-unitary-capable local operator API. Applies rho -> A rho A^dagger.
			// If A is not unitary the trace is generally not preserved; the normalized variants
			// condition on the operation succeeding by dividing the resulting state by its trace.
			virtual void ApplyOperator(const Gates::AppliedGate<MatrixClass>& op) = 0;
			virtual void ApplyOperator(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) = 0;
			virtual void ApplyOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) = 0;
			virtual void ApplyOperatorAndNormalize(const Gates::AppliedGate<MatrixClass>& op) = 0;
			virtual void ApplyOperatorAndNormalize(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) = 0;

			// Applies a multi-Kraus channel: rho -> sum_i K_i rho K_i^dagger.
			// The Kraus completeness relation is not checked; non trace-preserving maps are allowed.
			virtual void ApplyKrausOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) = 0;
			virtual void ApplyKrausOperators(const std::vector<MatrixClass>& ops, IndexType qubit, IndexType controllingQubit1 = 0) = 0;

			// ---- predefined single qubit noise channels ----
			//
			// Convenience wrappers over ApplyKrausOperators, expressed with the standard Kraus
			// operators of each channel. They are trace preserving (up to MPO truncation error).

			// bit flip: rho -> (1 - p) rho + p X rho X
			void ApplyBitFlipNoise(IndexType qubit, double p)
			{
				ApplyKrausOperators({ std::sqrt(1. - p) * NoisePauliI(), std::sqrt(p) * NoisePauliX() }, qubit);
			}

			// phase flip: rho -> (1 - p) rho + p Z rho Z
			void ApplyPhaseFlipNoise(IndexType qubit, double p)
			{
				ApplyKrausOperators({ std::sqrt(1. - p) * NoisePauliI(), std::sqrt(p) * NoisePauliZ() }, qubit);
			}

			// depolarizing: rho -> (1 - p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
			void ApplyDepolarizingNoise(IndexType qubit, double p)
			{
				const double s = std::sqrt(p / 3.);
				ApplyKrausOperators({ std::sqrt(1. - p) * NoisePauliI(), s * NoisePauliX(), s * NoisePauliY(), s * NoisePauliZ() }, qubit);
			}

			// amplitude damping (|1> -> |0> relaxation with probability gamma)
			void ApplyAmplitudeDamping(IndexType qubit, double gamma)
			{
				MatrixClass E0 = MatrixClass::Zero(2, 2);
				E0(0, 0) = 1.;
				E0(1, 1) = std::sqrt(1. - gamma);

				MatrixClass E1 = MatrixClass::Zero(2, 2);
				E1(0, 1) = std::sqrt(gamma);

				ApplyKrausOperators({ E0, E1 }, qubit);
			}

			// phase damping / dephasing, suppresses the off diagonal coherences by lambda = sqrt(1 - gamma)
			void ApplyPhaseDamping(IndexType qubit, double gamma)
			{
				MatrixClass E0 = MatrixClass::Zero(2, 2);
				E0(0, 0) = 1.;
				E0(1, 1) = std::sqrt(1. - gamma);

				MatrixClass E1 = MatrixClass::Zero(2, 2);
				E1(1, 1) = std::sqrt(gamma);

				ApplyKrausOperators({ E0, E1 }, qubit);
			}

			// reset a qubit to |0>: E0 = |0><0|, E1 = |0><1|
			void ApplyReset(IndexType qubit)
			{
				MatrixClass E0 = MatrixClass::Zero(2, 2);
				E0(0, 0) = 1.;

				MatrixClass E1 = MatrixClass::Zero(2, 2);
				E1(0, 1) = 1.;

				ApplyKrausOperators({ E0, E1 }, qubit);
			}

			virtual bool MeasureQubit(IndexType qubit) = 0;
			virtual std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) = 0;
			virtual double GetProbability(IndexType qubit, bool zeroVal = true) const = 0;

			// rho is a density matrix, so the basis-state 'amplitude' is a matrix element <row|rho|col>
			virtual std::complex<double> getBasisStateMatrixElement(size_t row, size_t col) const = 0;
			virtual std::complex<double> getBasisStateMatrixElement(const std::vector<bool>& row, const std::vector<bool>& col) const = 0;
			virtual double getBasisStateProbability(size_t State) const = 0;
			virtual double getBasisStateProbability(const std::vector<bool>& State) const = 0;

			// trace of the density matrix, should be 1 for a properly normalized state
			virtual std::complex<double> Trace() const = 0;

			virtual std::shared_ptr<MPOSimulatorStateInterface> getState() const = 0;
			virtual void setState(const std::shared_ptr<MPOSimulatorStateInterface>& state) = 0;

		protected:
			// single-qubit Pauli matrices used to build the predefined noise channels
			static MatrixClass NoisePauliI()
			{
				MatrixClass m = MatrixClass::Identity(2, 2);
				return m;
			}

			static MatrixClass NoisePauliX()
			{
				MatrixClass m = MatrixClass::Zero(2, 2);
				m(0, 1) = 1.;
				m(1, 0) = 1.;
				return m;
			}

			static MatrixClass NoisePauliY()
			{
				MatrixClass m = MatrixClass::Zero(2, 2);
				m(0, 1) = std::complex<double>(0., -1.);
				m(1, 0) = std::complex<double>(0., 1.);
				return m;
			}

			static MatrixClass NoisePauliZ()
			{
				MatrixClass m = MatrixClass::Zero(2, 2);
				m(0, 0) = 1.;
				m(1, 1) = -1.;
				return m;
			}
		};

	}

}
