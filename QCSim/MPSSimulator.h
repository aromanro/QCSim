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
			using IndexType = Eigen::Index;
			using IntIndexPair = Eigen::IndexPair<int>;
			using Indexes = Eigen::array<IntIndexPair, 1>;

			template<class MatrixClass = Eigen::MatrixXcd> class ZeroProjection : public Gates::SingleQubitGate<MatrixClass>
			{
			public:
				using BaseClass = Gates::SingleQubitGate<MatrixClass>;
				using OpClass = typename BaseClass::BaseClass;

				ZeroProjection()
				{
					OpClass::operatorMat(0, 0) = 1;
				}
			};

			template<class MatrixClass = Eigen::MatrixXcd> class OneProjection : public Gates::SingleQubitGate<MatrixClass>
			{
			public:
				using BaseClass = Gates::SingleQubitGate<MatrixClass>;
				using OpClass = typename BaseClass::BaseClass;

				OneProjection()
				{
					OpClass::operatorMat(1, 1) = 1;
				}
			};

			MPSSimulator(size_t N, int addseed = 0)
				: lambdas(N - 1, LambdaType::Ones(1)), gammas(N, GammaType(1, 2, 1)),
				uniformZeroOne(0, 1)
			{
				for (size_t i = 0; i < gammas.size(); ++i)
				{
					gammas[i](0, 0, 0) = 1.;
					gammas[i](0, 1, 0) = 0.;
				}

				const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
				std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
				rng.seed(seed);
			}

			void Clear()
			{
				const size_t szm1 = lambdas.size();
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
				const size_t szm1 = lambdas.size();
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
				else if (qubit >= gammas.size() || (gate.getQubitsNumber() == 2 && controllingQubit1 >= gammas.size()))
					throw std::runtime_error("Qubit index out of bounds");


				if (gate.getQubitsNumber() == 1)
					ApplySingleQubitGate(gate, qubit);
				else
					ApplyTwoQubitGate(gate, qubit, controllingQubit1);
			}

			// false if measured 0, true if measured 1
			bool Measure(size_t qubit)
			{
				const double rndVal = 1. - uniformZeroOne(rng);
				
				const double prob0 = GetProbability0(qubit);

				const bool zeroMeasured = rndVal < prob0;
				MatrixClass projMat;
				
				if (zeroMeasured)
				{
					projMat = zeroProjection.getRawOperatorMatrix();
					projMat *= 1. / sqrt(prob0);
				}
				else
				{
					projMat = oneProjection.getRawOperatorMatrix();
					projMat *= 1. / sqrt(1. - prob0);
				}

				const QC::Gates::SingleQubitGate projOp(projMat);

				ApplySingleQubitGate(projOp, qubit);

				// propagate to the other qubits to the left and right until the end or there is no entanglement with the next one

				for (size_t q = qubit; q > 0; --q)
				{
					const size_t q1 = q - 1;
					if (lambdas[q1].size() == 1) break;
					ApplyTwoQubitGate(projOp, q1, q, true);
				}

				for (size_t q = qubit; q < lambdas.size(); ++q)
				{
					if (lambdas[q].size() == 1) break;
					ApplyTwoQubitGate(projOp, q, q + 1, true);
				}

				return !zeroMeasured;
			}

			// this is for 'compatibility' with the statevector simulator (QubitRegister)
			// it's not stored as this and it's costly to compute, it will throw an exception for more than 64 qubits
			// but don't call it for such a large number of qubits
			VectorClass getRegisterStorage() const
			{
				if (gammas.size() > sizeof(size_t) * 8) throw std::runtime_error("Too many qubits to compute the state vector");

				// TODO: implement this with variadic templates, perhaps
				if (gammas.size() == 1)
					return GenerateStatevector<1>(gammas[0]);
				else if (gammas.size() == 2)
					return GenerateStatevector<2>(ContractNQubits<2>(gammas[0], lambdas[0], gammas[1]));
				else if (gammas.size() == 3)
					return GenerateStatevector<3>(ContractNQubits<3>(ContractNQubits<2>(gammas[0], lambdas[0], gammas[1]), lambdas[1], gammas[2]));
				else if (gammas.size() == 4)
					return GenerateStatevector<4>(ContractNQubits<4>(ContractNQubits<3>(ContractNQubits<2>(gammas[0], lambdas[0], gammas[1]), lambdas[1], gammas[2]), lambdas[2], gammas[3]));
				else if (gammas.size() == 5)
					return GenerateStatevector<5>(ContractNQubits<5>(ContractNQubits<4>(ContractNQubits<3>(ContractNQubits<2>(gammas[0], lambdas[0], gammas[1]), lambdas[1], gammas[2]), lambdas[2], gammas[3]), lambdas[3], gammas[4]));
				else if (gammas.size() == 6)
					return GenerateStatevector<6>(ContractNQubits<6>(ContractNQubits<5>(ContractNQubits<4>(ContractNQubits<3>(ContractNQubits<2>(gammas[0], lambdas[0], gammas[1]), lambdas[1], gammas[2]), lambdas[2], gammas[3]), lambdas[3], gammas[4]), lambdas[4], gammas[5]));
				else
					throw std::runtime_error("Not implemented yet");

				return {};
			}

			

			double GetProbability0(size_t qubit) const
			{
				MatrixClass qubitMatrix = GetQubitMatrix(qubit, 0);

				if (qubit > 0)
				{
					const size_t qbit1 = qubit - 1;
					for (IndexType col = 0; col < qubitMatrix.cols(); col++)
						for (IndexType row = 0; row < qubitMatrix.rows(); row++)
							qubitMatrix(row, col) *= (row < lambdas[qbit1].size()) ? lambdas[qbit1][row] : 0;
				}

				if (qubit < lambdas.size())
				{
					for (IndexType col = 0; col < qubitMatrix.cols(); col++)
						for (IndexType row = 0; row < qubitMatrix.rows(); row++)
							qubitMatrix(row, col) *= (row < lambdas[qubit].size()) ? lambdas[qubit][col] : 0;
				}

				return qubitMatrix.cwiseProduct(qubitMatrix.conjugate()).sum().real();
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

			void ApplyTwoQubitGate(const GateClass& gate, size_t qubit, size_t controllingQubit1, bool dontApplyGate = false)
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

				MatrixClass thetaMatrix;

				if (dontApplyGate)
				{
					const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);

					thetaMatrix = ReshapeTheta(theta);
				}
				else
				{
					const Eigen::Tensor<std::complex<double>, 4> U = GetTwoQubitsGateTensor(gate, reversed);
					const Eigen::Tensor<std::complex<double>, 4> thetabar = ConstructThetaBar(qubit1, U);

					thetaMatrix = ReshapeThetaBar(thetabar);
				}

				Eigen::JacobiSVD<MatrixClass> SVD;

				if (limitEntanglement)
					SVD.setThreshold(singularValueThreshold);

				SVD.compute(thetaMatrix, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV);


				const MatrixClass& UmatrixFull = SVD.matrixU();
				const MatrixClass& VmatrixFull = SVD.matrixV();

				const Eigen::VectorXd& SvaluesFull = SVD.singularValues();

				long long szm = SVD.nonzeroSingularValues();

				if (szm == 0) szm = 1; // Shouldn't happen (unless some big limit was put on 'zero')!

				//szm = std::min<long long>(szm, UmatrixFull.cols());
				const long long sz = limitSize ? std::min<long long>(chi, szm) : szm;

				long long Dchi = 2 * sz;
				const MatrixClass& Umatrix = UmatrixFull.topLeftCorner(std::min<long long>(Dchi, UmatrixFull.rows()), sz);
				const MatrixClass& Vmatrix = VmatrixFull.topLeftCorner(std::min<long long>(Dchi, VmatrixFull.rows()), sz).adjoint();

				const LambdaType Svalues = SvaluesFull.head(szm);
				
				// now set back lambdas and gammas

				const long long szl = (qubit1 == 0) ? 1 : sz;
				const long long szr = (qubit2 == lambdas.size()) ? 1 : sz;

				Eigen::Tensor<std::complex<double>, 3> Utensor(szl, 2, sz);
				Eigen::Tensor<std::complex<double>, 3> Vtensor(sz, 2, szr);

				if (sz != szl || sz != szr)
				{
					for (IndexType k = 0; k < sz; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < szl; ++i)
							{
								const IndexType jchi = j * szl;
								const IndexType jind = jchi + i;
								Utensor(i, j, k) = jind < Umatrix.rows() ? Umatrix(jind, k) : 0;
							}

					for (IndexType k = 0; k < szr; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < sz; ++i)
							{
								const IndexType jchi = j * szr;
								const IndexType jind = jchi + k;
								Vtensor(i, j, k) = jind < Vmatrix.cols() ? Vmatrix(i, jind) : 0;
							}
				}
				else
				{
					for (IndexType k = 0; k < sz; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < sz; ++i)
							{
								const IndexType jchi = j * sz;
								IndexType jind = jchi + i;
								Utensor(i, j, k) = jind < Umatrix.rows() ? Umatrix(jind, k) : 0;

								jind = jchi + k;
								Vtensor(i, j, k) = jind < Vmatrix.cols() ? Vmatrix(i, jind) : 0;
							}
				}

				if (qubit1 == 0)
					gammas[qubit1] = Utensor;
				else
				{
					const Eigen::Tensor<std::complex<double>, 2> lambdaLeftInv = GetInverseLambdaTensor(qubit1 - 1, Utensor.dimension(0));

					static const Eigen::array<Eigen::IndexPair<int>, 1> product_dims{ Eigen::IndexPair<int>(1, 0) };
					gammas[qubit1] = lambdaLeftInv.contract(Utensor, product_dims);
				}

				lambdas[qubit1] = Svalues;
				lambdas[qubit1].normalize();
				//if (lambdas[qubit1][0] == 0.) lambdas[qubit1][0] = 1;

				if (qubit2 == lambdas.size())
					gammas[qubit2] = Vtensor;
				else
				{
					const Eigen::Tensor<std::complex<double>, 2> lambdaRightInv = GetInverseLambdaTensor(qubit2, Vtensor.dimension(2));

					static const Eigen::array<Eigen::IndexPair<int>, 1> product_dims{ Eigen::IndexPair<int>(2, 0) };
					gammas[qubit2] = Vtensor.contract(lambdaRightInv, product_dims);
				}
			}


			Eigen::Tensor<std::complex<double>, 2> GetLambdaTensor(size_t pos, size_t dim1, size_t dim2) const
			{
				assert(pos < lambdas.size());

				Eigen::Tensor<std::complex<double>, 2> res(dim1, dim2);
				res.setZero();
				
				for (size_t i = 0; i < std::min<size_t>(lambdas[pos].size(), std::min(dim1, dim2)); ++i)
					res(i, i) = lambdas[pos](i);

				return res;
			}

			Eigen::Tensor<std::complex<double>, 2> GetInverseLambdaTensor(size_t pos, size_t dim) const
			{
				assert(pos < lambdas.size());

				Eigen::Tensor<std::complex<double>, 2> res(dim, dim);
				res.setZero();

				for (size_t i = 0; i < std::min<size_t>(lambdas[pos].size(), dim); ++i)
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
					for (int q0l = 0; q0l < 2; ++q0l)
					{
						const int l0 = q0l << 1;
						for (int q0c = 0; q0c < 2; ++q0c)
						{
							const int c0 = q0c << 1;
							for (int q1l = 0; q1l < 2; ++q1l)
								for (int q1c = 0; q1c < 2; ++q1c)
									result(q1l, q0l, q1c, q0c) = gate.getRawOperatorMatrix()(l0 | q1l, c0 | q1c);
						}
					}
				else
					for (int q0l = 0; q0l < 2; ++q0l) // ctrl qubit
					{
						const int l0 = q0l << 1;
						for (int q0c = 0; q0c < 2; ++q0c) // ctrl qubit
						{
							const int c0 = q0c << 1;
							for (int q1l = 0; q1l < 2; ++q1l)
								for (int q1c = 0; q1c < 2; ++q1c)
									result(q0l, q1l, q0c, q1c) = gate.getRawOperatorMatrix()(l0 | q1l, c0 | q1c);
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

				const Eigen::Tensor<std::complex<double>, 2> lambdaMiddle = GetLambdaTensor(qubit1, gammas[qubit1].dimension(2), gammas[qubit2].dimension(0) );

				if (qubit1 == 0 && qubit2 == lambdas.size())
				{
					return gammas[qubit1].contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int);
				}
				else if (qubit1 == 0 && qubit2 < lambdas.size())
				{
					const Eigen::Tensor<std::complex<double>, 2> lambdaRight = GetLambdaTensor(qubit2, gammas[qubit2].dimension(2), lambdas[qubit2].size());

					return gammas[qubit1].contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int).contract(lambdaRight, product_dims4);
				}
				else if (qubit1 > 0 && qubit2 == lambdas.size())
				{
					const Eigen::Tensor<std::complex<double>, 2> lambdaLeft = GetLambdaTensor(qubit1 - 1, lambdas[qubit1].size(), gammas[qubit1].dimension(0));

					return lambdaLeft.contract(gammas[qubit1], product_dims1).contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int);
				}

				const Eigen::Tensor<std::complex<double>, 2> lambdaLeft = GetLambdaTensor(qubit1 - 1, lambdas[qubit1 - 1].size(), gammas[qubit1].dimension(0));
				const Eigen::Tensor<std::complex<double>, 2> lambdaRight = GetLambdaTensor(qubit2, gammas[qubit2].dimension(2), lambdas[qubit2].size());

				return lambdaLeft.contract(gammas[qubit1], product_dims1).contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int).contract(lambdaRight, product_dims4);
			}

			// for the two qubit gate - U is the gate tensor
			Eigen::Tensor<std::complex<double>, 4> ConstructThetaBar(size_t qubit1, const Eigen::Tensor<std::complex<double>, 4>& U) const
			{
				const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);

				// apply evolution/gate operator
				using DimPair = Eigen::Tensor<std::complex<double>, 4>::DimensionPair;

				// from theta the physical indexes are contracted out
				// the last two become the physical indexes

				static const Eigen::array<DimPair, 2> product_dim{ DimPair(1, 2), DimPair(2, 3) };

				// this applies the time evolution operator U
				return theta.contract(U, product_dim);
			}

			MatrixClass ReshapeThetaBar(const Eigen::Tensor<std::complex<double>, 4>& theta)
			{
				// get it into a matrix for SVD - use JacobiSVD

				int sz0 = static_cast<int>(theta.dimension(0));
				int sz1 = static_cast<int>(theta.dimension(1));

				assert(2 == theta.dimension(2));
				assert(2 == theta.dimension(3));

				MatrixClass thetaMatrix(2 * sz0, 2 * sz1);


				for (IndexType l = 0; l < 2; ++l)
					for (IndexType j = 0; j < sz1; ++j)
						for (IndexType k = 0; k < 2; ++k)
							for (IndexType i = 0; i < sz0; ++i)
								// 2, 0, 3, 1 - k & l from theta are the physical indexes
								thetaMatrix(k * sz0 + i, l * sz1 + j) = theta(i, j, k, l);

				return thetaMatrix;
			}

			MatrixClass ReshapeTheta(const Eigen::Tensor<std::complex<double>, 4>& theta)
			{
				// get it into a matrix for SVD - use JacobiSVD

				int sz0 = static_cast<int>(theta.dimension(0));
				assert(2 == theta.dimension(1));
				assert(2 == theta.dimension(2));
				int sz3 = theta.dimension(3);

				MatrixClass thetaMatrix(2 * sz0, 2 * sz3);

				for (IndexType k = 0; k < 2; ++k)
					for (IndexType l = 0; l < sz3; ++l)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < sz0; ++i)
								// j & k from theta are the physical indexes
								thetaMatrix(j * sz0 + i, k * sz3 + l) = theta(i, j, k, l);

				return thetaMatrix;
			}

			MatrixClass GetQubitMatrix(size_t qubit, unsigned short outcome) const
			{
				assert(qubit < gammas.size());
				assert(outcome < 2);
				
				const GammaType& qubitTensor = gammas[qubit];
				MatrixClass res(qubitTensor.dimension(0), qubitTensor.dimension(2));

				for (IndexType j = 0; j < res.cols(); ++j)
					for (IndexType i = 0; i < res.rows(); ++i)
						res(i, j) = qubitTensor(i, outcome, j);

				return res;
			}

			template<int N> static Eigen::Tensor<std::complex<double>, N + 2> ContractNQubits(const Eigen::Tensor<std::complex<double>, N + 1>& left, const LambdaType& lambdaVal, const GammaType& nextQubit)
			{
				const IndexType dim1 = left.dimension(N);
				const IndexType dim2 = nextQubit.dimension(0);

				Eigen::Tensor<std::complex<double>, 2> lambdaTensor(dim1, dim2);
				lambdaTensor.setZero();

				for (size_t i = 0; i < std::min<size_t>(lambdaVal.size(), std::min(dim1, dim2)); ++i)
					lambdaTensor(i, i) = lambdaVal(i);

				static const Indexes productDim{ IntIndexPair(N, 0) };

				return left.contract(lambdaTensor, productDim).contract(nextQubit, productDim);
			}

			template<int N> static VectorClass GenerateStatevector(const Eigen::Tensor<std::complex<double>, N + 2>& tensor)
			{
				const size_t NrBasisStates = 1ULL << N;
				VectorClass res(NrBasisStates);

				for (size_t state = 0; state < NrBasisStates; ++state)
				{
					// index for tensor
					std::array<IndexType, N + 2> indices;

					size_t tmp = state;

					indices[0] = 0;
					indices[N + 1] = 0;
					for (size_t q = 1; q < N + 1; ++q)
					{
						indices[q] = tmp & 1;
						tmp >>= 1;
					}

					res(state) = tensor(indices);
				}

				return res;
			}
				
			bool limitSize = false;
			bool limitEntanglement = false;
			size_t chi = 10; // if limitSize is true
			double singularValueThreshold = 0.; // if limitEntanglement is true

			std::vector<LambdaType> lambdas;
			std::vector<GammaType> gammas;

			std::mt19937_64 rng;
			std::uniform_real_distribution<double> uniformZeroOne;

			const ZeroProjection<MatrixClass> zeroProjection;
			const OneProjection<MatrixClass> oneProjection;
		};

	}
}
