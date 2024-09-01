#pragma once

#include "MPSSimulatorBase.h"

namespace QC {

	namespace TensorNetworks {

		class MPSSimulatorImpl : public MPSSimulatorBase
		{
		public:
			MPSSimulatorImpl(size_t N, int addseed = 0)
				: MPSSimulatorBase(N, addseed)
			{
			}

			void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) override
			{
				ApplyGate(gate, gate.getQubit1(), gate.getQubit2());
			}

			// three qubit gates not supported, convert them to two qubit gates
			// also two qubit gates need to act on adjacent qubits
			// don't try to apply a gate that doesn't satisfy these conditions
			// use swap gates to move qubits around
			// maybe wrap this up into a higher level simulator that swaps the qubits for you and maps them to minimize swaps
			void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (gate.getQubitsNumber() > 2) throw std::runtime_error("Three qubit gates not supported");
				else if ((qubit < 0 || qubit >= static_cast<IndexType>(gammas.size())) || (gate.getQubitsNumber() == 2 && (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(gammas.size()))))
					throw std::runtime_error("Qubit index out of bounds");
				else if (gate.getQubitsNumber() == 2 && std::abs(static_cast<int>(qubit) - static_cast<int>(controllingQubit1)) != 1)
					throw std::runtime_error("Two qubit gates need to act on adjacent qubits");


				if (gate.getQubitsNumber() == 1)
					ApplySingleQubitGate(gate, qubit);
				else
					ApplyTwoQubitGate(gate, qubit, controllingQubit1);
			}

			void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) override
			{
				for (const auto& gate : gates)
					ApplyGate(gate);
			}

			// false if measured 0, true if measured 1
			bool MeasureQubit(IndexType qubit) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
					throw std::runtime_error("Qubit index out of bounds");

				const double rndVal = 1. - uniformZeroOne(rng);

				const double prob0 = GetProbability(qubit);

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

				for (IndexType q = qubit; q > 0; --q)
				{
					const IndexType q1 = q - 1;
					if (lambdas[q1].size() == 1) break;
					ApplyTwoQubitGate(projOp, q1, q, true);
				}

				for (IndexType q = qubit; q < static_cast<IndexType>(lambdas.size()); ++q)
				{
					if (lambdas[q].size() == 1) break;
					ApplyTwoQubitGate(projOp, q, q + 1, true);
				}

				return !zeroMeasured;
			}

			double GetProbability(IndexType qubit, bool zeroVal = true) const
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
					throw std::runtime_error("Qubit index out of bounds");

				MatrixClass qubitMatrix = GetQubitMatrix(qubit, zeroVal ? 0 : 1);

				if (qubit > 0)
				{
					const IndexType qbit1 = qubit - 1;
					for (IndexType col = 0; col < qubitMatrix.cols(); col++)
						for (IndexType row = 0; row < qubitMatrix.rows(); row++)
							qubitMatrix(row, col) *= (row < lambdas[qbit1].size()) ? lambdas[qbit1][row] : 0;
				}

				if (qubit < static_cast<IndexType>(lambdas.size()))
				{
					for (IndexType col = 0; col < qubitMatrix.cols(); col++)
						for (IndexType row = 0; row < qubitMatrix.rows(); row++)
							qubitMatrix(row, col) *= (col < lambdas[qubit].size()) ? lambdas[qubit][col] : 0;
				}

				return qubitMatrix.cwiseProduct(qubitMatrix.conjugate()).sum().real();
			}

		private:
			void ApplySingleQubitGate(const GateClass& gate, IndexType qubit)
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

			void ApplyTwoQubitGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1, bool dontApplyGate = false)
			{
				// it's more complex than the single qubit gate
				// very shortly:
				// contract tensors for the two qubits, along with the correspnding lambdas
				// shape the gate into a tensor
				// contract the gate tensor with the two qubit tensor
				// apply SVD to separate out the resulting tensor into the two qubit tensors and the lambdas

				IndexType qubit1 = controllingQubit1;
				IndexType qubit2 = qubit;
				bool reversed = false;

				if (qubit1 > qubit2)
				{
					std::swap(qubit1, qubit2);
					reversed = true;
				}

				MatrixClass thetaMatrix;

				if (dontApplyGate)
				{
					// this case is for dealing with entanglement with the qubits entangled with a measured one
					// see 'Measure' function for details
					const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);

					thetaMatrix = ReshapeTheta(theta);
				}
				else
				{
					const Eigen::Tensor<std::complex<double>, 4> U = GetTwoQubitsGateTensor(gate, reversed);
					const Eigen::Tensor<std::complex<double>, 4> thetabar = ConstructThetaBar(qubit1, U);

					thetaMatrix = ReshapeThetaBar(thetabar);
				}

				Eigen::JacobiSVD<MatrixClass/*, Eigen::FullPivHouseholderQRPreconditioner*/> SVD;
				//Eigen::BDCSVD<MatrixClass> SVD;

				if (limitEntanglement)
					SVD.setThreshold(singularValueThreshold);

#ifdef _DEBUG
				std::cout << "Gamma1 dims: " << gammas[qubit1].dimension(0) << " x " << gammas[qubit1].dimension(1) << " x " << gammas[qubit1].dimension(2) << std::endl;
				std::cout << "Gamma2 dims: " << gammas[qubit2].dimension(0) << " x " << gammas[qubit2].dimension(1) << " x " << gammas[qubit2].dimension(2) << std::endl;
				std::cout << "Lambda size: " << lambdas[qubit1].size() << std::endl;
				std::cout << "Theta matrix size: " << thetaMatrix.rows() << " x " << thetaMatrix.cols() << std::endl;
#endif

				SVD.compute(thetaMatrix, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV);


				const MatrixClass& UmatrixFull = SVD.matrixU();
				const MatrixClass& VmatrixFull = SVD.matrixV();
				const LambdaType& SvaluesFull = SVD.singularValues();

#ifdef _DEBUG
				std::cout << "U matrix size: " << UmatrixFull.rows() << " x " << UmatrixFull.cols() << std::endl;
				std::cout << UmatrixFull << std::endl;
				std::cout << "V matrix size: " << VmatrixFull.rows() << " x " << VmatrixFull.cols() << std::endl;
				std::cout << VmatrixFull << std::endl;
				std::cout << "Singular values: " << SvaluesFull << std::endl;
#endif

				//IndexType szm = SvaluesFull.size();
				IndexType szm = SVD.nonzeroSingularValues(); // or SvaluesFull.size() for tests

				//if (szm < SvaluesFull.size())
				//	std::cout << "Singular values truncated" << std::endl;

				//assert(szm > 0);

				if (szm == 0) szm = 1; // Shouldn't happen (unless some big limit was put on 'zero')!

				const IndexType sz = limitSize ? std::min<IndexType>(chi, szm) : szm;


				const IndexType szl = qubit1 == 0 ? 1 : lambdas[qubit1 - 1].size();
				const IndexType szr = qubit2 == static_cast<IndexType>(lambdas.size()) ? 1 : lambdas[qubit2].size();

				assert(UmatrixFull.cols() == VmatrixFull.cols()); // for 'thin'
				assert(sz <= UmatrixFull.cols());

				//const MatrixClass Umatrix = UmatrixFull;
				//const MatrixClass Vmatrix = VmatrixFull.adjoint();

				const MatrixClass Umatrix = UmatrixFull.topLeftCorner(std::min<IndexType>(2 * szl, UmatrixFull.rows()), sz);
				const MatrixClass Vmatrix = VmatrixFull.topLeftCorner(std::min<IndexType>(2 * szr, VmatrixFull.rows()), sz).adjoint();

#ifdef _DEBUG
				std::cout << "Trunc U matrix size: " << Umatrix.rows() << " x " << Umatrix.cols() << std::endl;
				std::cout << "Trunc V matrix size: " << Vmatrix.rows() << " x " << Vmatrix.cols() << std::endl;
#endif

				// now set back lambdas and gammas

				lambdas[qubit1] = SvaluesFull.head(sz);
				assert(lambdas[qubit1][0] != 0.);
				//if (lambdas[qubit1][0] == 0.) lambdas[qubit1][0] = 1; // this should not happen
				lambdas[qubit1].normalize();

#ifdef _DEBUG
				std::cout << "Normalized: " << lambdas[qubit1] << std::endl;
#endif

				SetNewGammas(Umatrix, Vmatrix, qubit1, qubit2, szl, sz, szr);
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


			Eigen::Tensor<std::complex<double>, 4> GetTwoQubitsGateTensor(const GateClass& gate, bool reversed) const
			{
				Eigen::Tensor<std::complex<double>, 4> result(2, 2, 2, 2);

				const auto& gateMat = gate.getRawOperatorMatrix();

				if (reversed)
					for (int q0l = 0; q0l < 2; ++q0l)
					{
						const int l0 = q0l << 1;
						for (int q0c = 0; q0c < 2; ++q0c)
						{
							const int c0 = q0c << 1;
							for (int q1l = 0; q1l < 2; ++q1l)
								for (int q1c = 0; q1c < 2; ++q1c)
									result(q1l, q0l, q1c, q0c) = gateMat(l0 | q1l, c0 | q1c);
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
									result(q0l, q1l, q0c, q1c) = gateMat(l0 | q1l, c0 | q1c);
						}
					}


				return result;
			}

			inline void SetNewGammas(const MatrixClass& Umatrix, const MatrixClass& Vmatrix, IndexType qubit1, IndexType qubit2, IndexType szl, IndexType sz, IndexType szr)
			{
				if (sz != szl || sz != szr)
					SetNewGammasDif(Umatrix, Vmatrix, qubit1, qubit2, szl, sz, szr);
				else
					SetNewGammasSame(Umatrix, Vmatrix, qubit1, qubit2, sz);

				DivideGammasWithLambdas(qubit1, qubit2, szl, sz, szr);
			}

			inline void SetNewGammasDif(const MatrixClass& Umatrix, const MatrixClass& Vmatrix, IndexType qubit1, IndexType qubit2, IndexType szl, IndexType sz, IndexType szr)
			{
				Eigen::Tensor<std::complex<double>, 3> Utensor(szl, 2, sz);
				Eigen::Tensor<std::complex<double>, 3> Vtensor(sz, 2, szr);

#ifdef _DEBUG
				std::cout << "Different sizes, szl=" << szl << " sz=" << sz << " szr=" << szr << std::endl;
#endif

				for (IndexType k = 0; k < sz; ++k)
					for (IndexType j = 0; j < 2; ++j)
						for (IndexType i = 0; i < szl; ++i)
						{
							const IndexType jind = j * szl + i;
							Utensor(i, j, k) = (jind < Umatrix.rows()) ? Umatrix(jind, k) : 0;
						}

				for (IndexType k = 0; k < szr; ++k)
					for (IndexType j = 0; j < 2; ++j)
						for (IndexType i = 0; i < sz; ++i)
						{
							const IndexType jind = j * szr + k;
							Vtensor(i, j, k) = (jind < Vmatrix.cols()) ? Vmatrix(i, jind) : 0;
						}

				gammas[qubit1] = Utensor;
				gammas[qubit2] = Vtensor;
			}

			inline void SetNewGammasSame(const MatrixClass& Umatrix, const MatrixClass& Vmatrix, IndexType qubit1, IndexType qubit2, IndexType sz)
			{
				Eigen::Tensor<std::complex<double>, 3> Utensor(sz, 2, sz);
				Eigen::Tensor<std::complex<double>, 3> Vtensor(sz, 2, sz);

#ifdef _DEBUG
				std::cout << "Same sizes, sz=" << sz << std::endl;
#endif

				for (IndexType k = 0; k < sz; ++k)
					for (IndexType j = 0; j < 2; ++j)
						for (IndexType i = 0; i < sz; ++i)
						{
							const IndexType jchi = j * sz;
							IndexType jind = jchi + i;
							Utensor(i, j, k) = (jind < Umatrix.rows()) ? Umatrix(jind, k) : 0;

							jind = jchi + k;
							Vtensor(i, j, k) = (jind < Vmatrix.cols()) ? Vmatrix(i, jind) : 0;
						}

				gammas[qubit1] = Utensor;
				gammas[qubit2] = Vtensor;
			}

			inline void DivideGammasWithLambdas(IndexType qubit1, IndexType qubit2, IndexType szl, IndexType sz, IndexType szr)
			{
				assert(gammas[qubit1].dimension(0) == szl);
				assert(gammas[qubit1].dimension(2) == sz);
				assert(gammas[qubit2].dimension(0) == sz);
				assert(gammas[qubit2].dimension(2) == szr);

				if (qubit1 != 0)
				{
					const IndexType prev = qubit1 - 1;
					for (IndexType k = 0; k < sz; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < szl; ++i)
								if (lambdas[prev][i] > std::numeric_limits<double>::epsilon() * 1E-10) gammas[qubit1](i, j, k) /= lambdas[prev][i];
								else gammas[qubit1](i, j, k) = 0;
				}

				if (qubit2 != static_cast<IndexType>(lambdas.size()))
				{
					for (IndexType j = 0; j < 2; ++j)
						for (IndexType i = 0; i < sz; ++i)
							for (IndexType k = 0; k < szr; ++k)
								if (lambdas[qubit2][k] > std::numeric_limits<double>::epsilon() * 1E-10) gammas[qubit2](i, j, k) /= lambdas[qubit2][k];
								else gammas[qubit2](i, j, k) = 0;
				}
			}

			Eigen::Tensor<std::complex<double>, 4> ContractTwoQubits(IndexType qubit1)
			{
				const IndexType qubit2 = qubit1 + 1;

				static const Indexes product_dims_int{ IntIndexPair(2, 0) };

				assert(gammas[qubit1].dimension(2) == gammas[qubit2].dimension(0));

				const Eigen::Tensor<std::complex<double>, 2> lambdaMiddle = GetLambdaTensor(qubit1, gammas[qubit1].dimension(2), gammas[qubit2].dimension(0));

				const IndexType szl = gammas[qubit1].dimension(0);
				const IndexType sz = gammas[qubit1].dimension(2);
				const IndexType szr = gammas[qubit2].dimension(2);

				if (qubit1 != 0)
				{
					const IndexType prev = qubit1 - 1;
					for (IndexType k = 0; k < sz; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < szl; ++i)
								gammas[qubit1](i, j, k) *= lambdas[prev][i];
				}

				if (qubit2 != static_cast<IndexType>(lambdas.size()))
				{
					for (IndexType j = 0; j < 2; ++j)
						for (IndexType i = 0; i < sz; ++i)
							for (IndexType k = 0; k < szr; ++k)
								gammas[qubit2](i, j, k) *= lambdas[qubit2][k];
				}

				// contract first gamma with the lambda in the middle
				// the resulting tensor has three legs, 1 is the physical one

				// then

				// contract the result with the next gamma
				// the resulting tensor has four legs, 1 and 2 are the physical ones

				return gammas[qubit1].contract(lambdaMiddle, product_dims_int).contract(gammas[qubit2], product_dims_int);
			}

			// for the two qubit gate - U is the gate tensor
			Eigen::Tensor<std::complex<double>, 4> ConstructThetaBar(IndexType qubit1, const Eigen::Tensor<std::complex<double>, 4>& U)
			{
				const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);

				// apply evolution/gate operator
				using DimPair = Eigen::Tensor<std::complex<double>, 4>::DimensionPair;

				// from theta the physical indexes are contracted out
				// the last two become the physical indexes

				static const Eigen::array<DimPair, 2> product_dim{ DimPair(1, 2), DimPair(2, 3) };

				// this applies the time evolution/gate operator U
				return theta.contract(U, product_dim);
			}

			MatrixClass ReshapeThetaBar(const Eigen::Tensor<std::complex<double>, 4>& theta)
			{
				// get it into a matrix for SVD - use JacobiSVD
				IndexType sz0 = theta.dimension(0);
				IndexType sz1 = theta.dimension(1);

				assert(sz0 > 0);
				assert(sz1 > 0);

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

				IndexType sz0 = theta.dimension(0);
				assert(sz0 > 0);

				assert(2 == theta.dimension(1));
				assert(2 == theta.dimension(2));

				IndexType sz3 = theta.dimension(3);
				assert(sz3 > 0);

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
		};

	}
}
