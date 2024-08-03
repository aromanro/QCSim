#pragma once

#include <vector>

#include <unsupported/Eigen/CXX11/Tensor>

#include <Eigen/Eigen>

#include "QuantumGate.h"

namespace QC {

	namespace TensorNetworks {

		class MPSSimulator
		{
		public:
			using LambdaType = Eigen::VectorXd;
			using GammaType = Eigen::Tensor<std::complex<double>, 3>;
			using MatrixClass = Eigen::MatrixXcd;
			using VectorClass = Eigen::VectorXcd;
			using GateClass = Gates::QuantumGateWithOp<MatrixClass>;

			MPSSimulator(size_t N)
				: lambdas(N - 1, LambdaType::Ones(1)), gammas(N, GammaType(1, 2, 1))
			{
				for (size_t i = 0; i < gammas.size(); ++i)
				{
					gammas[i](0, 0, 0) = 1.;
					gammas[i](0, 1, 0) = 0.;
				}
			}

			void Clear()
			{
				const size_t szm1 = gammas.size() - 1;
				for (size_t i = 0; i < szm1; ++i)
				{
					gammas[i].resize(1, 2, 1);
					gammas[i](0, 0, 0) = 1.;
					gammas[i](0, 1, 0) = 0.;

					lambdas[i].resize(1);
					lambdas[i](0) = 1.;
				}

				gammas[szm1].resize(1, 2, 1);
				gammas[szm1](0, 0, 0) = 1.;
				gammas[szm1](0, 1, 0) = 0.;
			}

			void InitOnesState()
			{
				const size_t szm1 = gammas.size() - 1;
				for (size_t i = 0; i < szm1; ++i)
				{
					gammas[i].resize(1, 2, 1);
					gammas[i](0, 0, 0) = 0.;
					gammas[i](0, 1, 0) = 1.;

					lambdas[i].resize(1);
					lambdas[i](0) = 1.;
				}

				gammas[szm1].resize(1, 2, 1);
				gammas[szm1](0, 0, 0) = 0.;
				gammas[szm1](0, 1, 0) = 1.;
			}

			void setToQubitState(size_t q)
			{
				Clear();
				if (q >= gammas.size())
					return;

				gammas[q](0, 0, 0) = 0.;
				gammas[q](0, 1, 0) = 1.;
			}

			void setToBasisState(size_t State)
			{
				const size_t NrBasisStates = gammas.size() > sizeof(size_t) * 8 ? 64 : (1ULL << gammas.size());
				if (State >= NrBasisStates) return;

				Clear();

				size_t pos = 0;
				while (State)
				{
					if (State & 1)
					{
						gammas[pos](0, 0, 0) = 0.;
						gammas[pos](0, 1, 0) = 1.;
					}
					State >>= 1;
					++pos;
				}
			}

			void setToBasisState(const std::vector<bool>& State)
			{
				if (State.size() > gammas.size()) return;

				Clear();

				for (size_t i = 0; i < State.size(); ++i)
				{
					if (State[i] == 1)
					{
						gammas[i](0, 0, 0) = 0.;
						gammas[i](0, 1, 0) = 1.;
					}
				}
			}

			void setLimitBondDimension(size_t chival)
			{
				limitSize = true;
				chi = chival;
			}

			void setLimitEntanglement(double svdThreshold)
			{
				limitEntanglement = true;
				singularValueThreshold = svdThreshold;
			}

			void dontLimitBondDimension()
			{
				limitSize = false;
			}

			void dontLimitEntanglement()
			{
				limitEntanglement = false;
			}

			void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate)
			{
				ApplyGate(gate, gate.getQubit1(), gate.getQubit2());
			}

			void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates)
			{
				for (const auto& gate : gates)
					ApplyGate(gate);
			}

			// three qubit gates not supported, convert them to two qubit gates
			// also two qubit gates need to act on adjacent qubits
			// don't try to apply a gate that doesn't satisfy these conditions
			// use swap gates to move qubits around
			// maybe wrap this up into a higher level simulator that swaps the qubits for you and maps them to minimize swaps
			void ApplyGate(const GateClass& gate, size_t qubit, size_t controllingQubit1 = 0)
			{
				if (gate.getQubitsNumber() > 2) throw std::runtime_error("Three qubit gates not supported");
				else if (gate.getQubitsNumber() == 2 && std::abs(static_cast<int>(qubit) - static_cast<int>(controllingQubit1)) != 1)
					throw std::runtime_error("Two qubit gates need to act on adjacent qubits");

				if (gate.getQubitsNumber() == 1)
					ApplySingleQubitGate(gate, qubit);
				else if (gate.getQubitsNumber() == 2)
					ApplyTwoQubitGate(gate, qubit, controllingQubit1);
			}

			// false if measured 0, true if measured 1
			bool Measure(size_t qubit)
			{
				// TODO: Implement it

				return false; // for now
			}

			// this is for 'compatibility' with the statevector simulator (QubitRegister)
			// it's not stored as this and it's costly to compute, it will throw an exception for more than 64 qubits
			// but don't call it for such a large number of qubits
			VectorClass getRegisterStorage() const
			{
				if (gammas.size() > sizeof(size_t) * 8) throw std::runtime_error("Too many qubits to compute the state vector");
				const size_t NrBasisStates = 1ULL << gammas.size();

				VectorClass res = VectorClass::Zero(NrBasisStates);

				// TODO: need to contract all the qubit tensors (gammas and lambdas corresponding to each qubit) into a single tensor
				// the physical indices of the resulting tensor remain

				return res;
			}

		private:
			void ApplySingleQubitGate(const GateClass& gate, size_t qubit)
			{
				// TODO: Implement it
				
				// easy: shape the gate into a tensor and contract it with the qubit tensor
			}

			void ApplyTwoQubitGate(const GateClass& gate, size_t qubit, size_t controllingQubit1)
			{
				// TODO: Implement it

				// it's more complex than the single qubit gate
				// very shortly:
				// contract tensors for the two qubits, along with the correspnding lambdas
				// shape the gate intro a tensor
				// contract the gate tensor with the two qubit tensor
				// apply SVD to separate out the resulting tensor into the two qubit tensors and the lambdas
			}

			bool limitSize = false;
			bool limitEntanglement = false;
			size_t chi = 0; // if limitSize is true
			double singularValueThreshold = 0.; // if limitEntanglement is true

			std::vector<LambdaType> lambdas;
			std::vector<GammaType> gammas;
		};

	}
}
