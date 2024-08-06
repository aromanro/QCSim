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
			using IntIndexPair = Eigen::IndexPair<int>;
			using Indexes = Eigen::array<IntIndexPair, 1>;

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
				else if (qubit >= gammas.size() || controllingQubit1 >= gammas.size())
					throw std::runtime_error("Qubit index out of bounds");


				if (gate.getQubitsNumber() == 1)
					ApplySingleQubitGate(gate, qubit);
				else
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

			void print() const
			{
				for (size_t i = 0; i < gammas.size() - 1; ++i)
				{
					std::cout << "Lambda " << i << ":\n" << lambdas[i] << std::endl;
					std::cout << "Gamma " << i << ":\n" << gammas[i] << std::endl;
				}

				std::cout << "Gamma " << gammas.size() - 1 << ":\n" << gammas[gammas.size() - 1] << std::endl;
			}

		private:
			void ApplySingleQubitGate(const GateClass& gate, size_t qubit)
			{
				// easy: shape the gate into a tensor and contract it with the qubit tensor
				Eigen::Tensor<std::complex<double>, 2> opTensor(2, 2);

				const MatrixClass& opMat = gate.getRawOperatorMatrix();

				opTensor(0, 0) = opMat(0, 0);
				opTensor(0, 1) = opMat(0, 1);
				opTensor(1, 0) = opMat(1, 0);
				opTensor(1, 1) = opMat(1, 1);

				// contract the gate tensor with the qubit tensor

				static const Indexes product_dims1{ IntIndexPair(1, 1) };
				static const std::array<int, 3> permute{ 0, 2, 1 };
				gammas[qubit] = gammas[qubit].contract(opTensor, product_dims1).shuffle(permute);
			}

			void ApplyTwoQubitGate(const GateClass& gate, size_t qubit, size_t controllingQubit1)
			{
				// it's more complex than the single qubit gate
				// very shortly:
				// contract tensors for the two qubits, along with the correspnding lambdas
				// shape the gate into a tensor
				// contract the gate tensor with the two qubit tensor
				// apply SVD to separate out the resulting tensor into the two qubit tensors and the lambdas

				size_t qubit1 = controllingQubit1;
				size_t qubit2 = qubit;
				bool reversed = false;

				if (qubit1 > qubit2)
				{
					std::swap(qubit1, qubit2);
					reversed = true;
				}

				const Eigen::Tensor<std::complex<double>, 4> U = GetTwoQubitsGateTensor(gate, reversed);
				const Eigen::Tensor<std::complex<double>, 4> thetabar = ConstructTheta(qubit1, U);

				const Eigen::MatrixXcd thetaMatrix = ReshapeTheta(thetabar);

				Eigen::JacobiSVD<Eigen::MatrixXcd> SVD; 

				if (limitEntanglement)
					SVD.setThreshold(singularValueThreshold);

				SVD.compute(thetaMatrix, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV);

				const Eigen::MatrixXcd& UmatrixFull = SVD.matrixU();

				const int sz = limitSize ? std::min<size_t>(chi, UmatrixFull.cols()) : UmatrixFull.cols();

				const int Dchi = 2 * sz;
				const Eigen::MatrixXcd& Umatrix = UmatrixFull.topLeftCorner(Dchi, sz);
				const Eigen::MatrixXcd& Vmatrix = SVD.matrixV().topLeftCorner(Dchi, sz).adjoint();

				const LambdaType Svalues = SVD.singularValues().head(sz);
				
				// now set back lambdas and gammas

				// this with lambdas is not optimal
				Eigen::Tensor<std::complex<double>, 3> Utensor(sz, 2, sz);
				Eigen::Tensor<std::complex<double>, 3> Vtensor(sz, 2, sz);

				for (Eigen::Index i = 0; i < sz; ++i)
					for (Eigen::Index j = 0; j < 2; ++j)
						for (Eigen::Index k = 0; k < sz; ++k)
						{
							const Eigen::Index jchi = j * sz;
							Utensor(i, j, k) = Umatrix(jchi + i, k);
							Vtensor(i, j, k) = Vmatrix(i, jchi + k);
						}

				if (qubit1 == 0)
					gammas[qubit1] = Utensor;
				else
				{
					Eigen::Tensor<std::complex<double>, 2> lambdaLeftInv = GetInverseLambdaTensor(qubit1 - 1);

					static const Eigen::array<Eigen::IndexPair<int>, 1> product_dims{ Eigen::IndexPair<int>(1, 0) };
					gammas[qubit1] = lambdaLeftInv.contract(Utensor, product_dims);
				}

				lambdas[qubit1] = Svalues;
				lambdas[qubit1].normalize();

				if (qubit2 == lambdas.size())
					gammas[qubit2] = Vtensor;
				else
				{
					Eigen::Tensor<std::complex<double>, 2> lambdaRightInv = GetInverseLambdaTensor(qubit2);

					static const Eigen::array<Eigen::IndexPair<int>, 1> product_dims{ Eigen::IndexPair<int>(2, 0) };
					gammas[qubit2] = Vtensor.contract(lambdaRightInv, product_dims);
				}
			}


			Eigen::Tensor<std::complex<double>, 2> GetLambdaTensor(size_t pos) const
			{
				assert(pos < lambdas.size());

				Eigen::Tensor<std::complex<double>, 2> res(lambdas[pos].size(), lambdas[pos].size());
				res.setZero();
				
				for (size_t i = 0; i < lambdas[pos].size(); ++i)
					res(i, i) = lambdas[pos](i);

				return res;
			}

			Eigen::Tensor<std::complex<double>, 2> GetInverseLambdaTensor(size_t pos) const
			{
				assert(pos < lambdas.size());

				Eigen::Tensor<std::complex<double>, 2> res(lambdas[pos].size(), lambdas[pos].size());
				res.setZero();

				for (size_t i = 0; i < lambdas[pos].size(); ++i)
				{
					if (abs(lambdas[pos](i)) > 1E-17)
						res(i, i) = 1. / lambdas[pos](i);
				}

				return res;
			}


			Eigen::Tensor<std::complex<double>, 4> GetTwoQubitsGateTensor(const GateClass& gate, bool reversed) const
			{
				Eigen::Tensor<std::complex<double>, 4> result(2, 2, 2, 2);

				if (reversed)
					for (int q1 = 0; q1 < 2; ++q1)
					{
						const int m1 = q1 << 1;
						for (int q2 = 0; q2 < 2; ++q2)
						{
							const int m2 = q2 << 1;
							for (int q3 = 0; q3 < 2; ++q3)
								for (int q4 = 0; q4 < 2; ++q4)
									result(q3, q1, q4, q2) = gate.getRawOperatorMatrix()(m1 | q3, m2 | q4);
						}
					}
				else
					for (int q1 = 0; q1 < 2; ++q1)
					{
						const int m1 = q1 << 1;
						for (int q2 = 0; q2 < 2; ++q2)
						{
							const int m2 = q2 << 1;
							for (int q3 = 0; q3 < 2; ++q3)
								for (int q4 = 0; q4 < 2; ++q4)
									result(q1, q3, q2, q4) = gate.getRawOperatorMatrix()(m1 | q3, m2 | q4);
						}
					}

				return result;
			}

			Eigen::Tensor<std::complex<double>, 4> ContractTwoQubits(size_t qubit1) const
			{
				const size_t qubit2 = qubit1 + 1;

				static const Indexes product_dims1{ IntIndexPair(1, 0) };
				static const Indexes product_dims_int{ IntIndexPair(2, 0) };
				static const Indexes product_dims4{ IntIndexPair(3, 0) };

				const Eigen::Tensor<std::complex<double>, 2> lambdaMiddle = GetLambdaTensor(qubit1);

				if (qubit1 == 0 && qubit2 == lambdas.size())
				{
					return gammas[qubit1].contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int);
				}
				else if (qubit1 == 0 && qubit2 < lambdas.size())
				{
					const Eigen::Tensor<std::complex<double>, 2> lambdaRight = GetLambdaTensor(qubit2);

					return gammas[qubit1].contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int).contract(lambdaRight, product_dims4);
				}
				else if (qubit1 > 0 && qubit2 == lambdas.size())
				{
					const Eigen::Tensor<std::complex<double>, 2> lambdaLeft = GetLambdaTensor(qubit1 - 1);

					return lambdaLeft.contract(gammas[qubit1], product_dims1).contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int);
				}

				const Eigen::Tensor<std::complex<double>, 2> lambdaLeft = GetLambdaTensor(qubit1 - 1);
				const Eigen::Tensor<std::complex<double>, 2> lambdaRight = GetLambdaTensor(qubit2);

				return lambdaLeft.contract(gammas[qubit1], product_dims1).contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int).contract(lambdaRight, product_dims4);
			}

			// for the two qubit gate - U is the gate tensor
			Eigen::Tensor<std::complex<double>, 4> ConstructTheta(size_t qubit1, const Eigen::Tensor<std::complex<double>, 4>& U) const
			{
				const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);

				// apply evolution/gate operator
				using DimPair = Eigen::Tensor<std::complex<double>, 4>::DimensionPair;

				// from theta the physical indexes are contracted out
				// the last two become the physical indexes

				//static const Eigen::array<DimPair, 2> product_dim{ DimPair(1, 0), DimPair(2, 1) };
				// hmmm... I'm unsure at this moment, I'm too tired, I'll check it tomorrow after I'll finish all of the 2-qubit gate apply code
				static const Eigen::array<DimPair, 2> product_dim{ DimPair(1, 2), DimPair(2, 3) };

				// this applies the time evolution operator U
				return theta.contract(U, product_dim);
			}

			Eigen::MatrixXcd ReshapeTheta(const Eigen::Tensor<std::complex<double>, 4>& theta)
			{
				// get it into a matrix for SVD - use JacobiSVD

				int sz = static_cast<int>(theta.dimension(0));
				assert(sz == theta.dimension(1));
				assert(2 == theta.dimension(2));
				assert(2 == theta.dimension(3));

				const int Dchi = 2 * sz;
				Eigen::MatrixXcd thetaMatrix(Dchi, Dchi);

				for (Eigen::Index i = 0; i < sz; ++i)
					for (Eigen::Index j = 0; j < sz; ++j)
						for (Eigen::Index k = 0; k < 2; ++k)
							for (Eigen::Index l = 0; l < 2; ++l)
								// 2, 0, 3, 1 - k & l from theta are the physical indexes
								thetaMatrix(k * sz + i, l * sz + j) = theta(i, j, k, l);

				return thetaMatrix;
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
