#pragma once

#include "MPSSimulatorBase.h"

#include <future>

namespace QC {

	namespace TensorNetworks {

		class MPSSimulatorImpl : public MPSSimulatorBase
		{
		public:
			friend class MPSSimulator;

			MPSSimulatorImpl(size_t N, unsigned int addseed = 0)
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
			// this is wrapped up into a higher level simulator that swaps the qubits for you and maps them to minimize swaps
			// so this should not be used directly
			void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (gate.getQubitsNumber() > 2) throw std::invalid_argument("Three qubit gates not supported");
				else if ((qubit < 0 || qubit >= static_cast<IndexType>(gammas.size())) || (gate.getQubitsNumber() == 2 && (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(gammas.size()))))
					throw std::invalid_argument("Qubit index out of bounds");
				else if (gate.getQubitsNumber() == 2 && std::abs(static_cast<int>(qubit) - static_cast<int>(controllingQubit1)) != 1)
					throw std::invalid_argument("Two qubit gates need to act on adjacent qubits");


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
					throw std::invalid_argument("Qubit index out of bounds");

				const double rndVal = 1. - uniformZeroOne(rng);
				const double prob0 = GetProbability(qubit);
				const bool zeroMeasured = rndVal < prob0;
				
				MatrixClass projMat = zeroMeasured ? zeroProjection.getRawOperatorMatrix() * 1. / sqrt(prob0) : oneProjection.getRawOperatorMatrix() * 1. / sqrt(1. - prob0);
				const QC::Gates::SingleQubitGate projOp(std::move(projMat));

				ApplySingleQubitGate(projOp, qubit);

				// propagate to the other qubits to the right and left until the end or there is no entanglement with the next one

				for (IndexType q = qubit; q < static_cast<IndexType>(lambdas.size()); ++q)
				{
					if (lambdas[q].size() == 1) break;
					ApplyTwoQubitGate(projOp, q, q + 1, true);
				}

				for (IndexType q = qubit; q > 0; --q)
				{
					const IndexType q1 = q - 1;
					if (lambdas[q1].size() == 1) break;
					ApplyTwoQubitGate(projOp, q1, q, true);
				}

				return !zeroMeasured;
			}

			std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) override
			{
				if (qubits.empty()) return {};

				std::unordered_map<IndexType, bool> res;
				if (qubits.size() == 1)
				{
					const IndexType qubit = *qubits.begin();
					res[qubit] = MeasureQubit(qubit);
					return res;
				}

				const auto lastIt = --qubits.end();

				for (auto it = qubits.begin(); it != lastIt; )
				{
					const auto qubit = *it;
					++it;
					const auto nextQubit = *it;

					const double rndVal = 1. - uniformZeroOne(rng);
					const double prob0 = GetProbability(qubit);
					const bool zeroMeasured = rndVal < prob0;

					res[qubit] = !zeroMeasured;
					
					MatrixClass projMat = zeroMeasured ? zeroProjection.getRawOperatorMatrix() * 1. / sqrt(prob0) : oneProjection.getRawOperatorMatrix() * 1. / sqrt(1. - prob0);
					const QC::Gates::SingleQubitGate projOp(std::move(projMat));

					ApplySingleQubitGate(projOp, qubit);
					
					for (IndexType q = qubit; q < nextQubit; ++q)
					{
						if (lambdas[q].size() == 1) break;
						ApplyTwoQubitGate(projOp, q, q + 1, true);
					}

					for (IndexType q = qubit; q > 0; --q)
					{
						const IndexType q1 = q - 1;
						if (lambdas[q1].size() == 1) break;
						ApplyTwoQubitGate(projOp, q1, q, true);
					}
				}

				const auto lastQubit = *lastIt;
				res[lastQubit] = MeasureQubit(lastQubit);

				return res;
			}

			// does not check for hermicity, that's why it returns a complex number
			// the caller should ensure the hermicity and extract the real part
			// also (for now, at least) it supports only one qubit ops
			// the problem with two qubit gates is that they will swap qubits around
			// so instead of only saving the state to compute <psi|U|psi>, and then use it to restore the state,
			// it would need to save the state twice, once for restoring and one for computing the expectation value - the last one having the qubits swapped as the one on which the gates are applied
			// anyway, this would be probably used mostly on Pauli strings, so...
			std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) override
			{
				if (gates.empty()) return 1.;

				// at this point it doesn't check if the gates are one qubit, that's done at the higher level
				// 	
				const IndexType lastQubit = static_cast<IndexType>(lambdas.size());
				
				// no need to zip up the whole chain, contracting the ends of the chain that are not touched by the operators should give 1.
				IndexType minQubit = lastQubit;
				IndexType maxQubit = 0;

				for (const auto& gate : gates)
				{
					const IndexType qubit = gate.getQubit1();
					minQubit = std::min(minQubit, qubit);
					maxQubit = std::max(maxQubit, qubit);
				}

				const IndexType nrSites = maxQubit - minQubit + 1;

				// lambdas are not modified by the single qubit gates

				std::vector<GammaType> modGammas(nrSites);
				for (IndexType s = 0; s < nrSites; ++s)
					modGammas[s] = gammas[minQubit + s];

				// the lambdas multiplication goes for both dagger and non-dagger gammas,
				// so the dagger are computed after the lambdas are applied
				MultiplyModGammasWithLambdas(modGammas, minQubit, nrSites);

				std::vector<GammaType> daggerGammas(nrSites);
				for (IndexType s = 0; s < nrSites; ++s)
					daggerGammas[s] = modGammas[s].conjugate();

				// apply the gates
				for (const auto& gate : gates)
					ApplySingleQubitGate(modGammas[gate.getQubit1() - minQubit], gate);

				// contract the saved dagger chain with the one with the gates applied
				// we need to do that only for the qubits in the range [minQubit, maxQubit]

				static const Eigen::array<IntIndexPair, 2> contract_dim{ IntIndexPair(0, 0), IntIndexPair(1, 1) };
				static const Indexes contract_dim1{ IntIndexPair(0, 0) };

				// start by contracting the first gamma with the first dagger gamma
				//  -O-         |
				// | |     ==>  O
				//  -O-         |
				Eigen::Tensor<std::complex<double>, 2> resTensor = modGammas[0].contract(daggerGammas[0], contract_dim);

				Eigen::Tensor<std::complex<double>, 2> resTensor2 = resTensor;

				// now contract the rest of the gammas with the dagger gammas
				for (IndexType s = 1; s < nrSites; ++s)
				{
					//  /-O-       -O-        |
					// O  |   ==> | |    ==>  O
					//  \-O-       -O-        |
					
					// bug in eigen, does not work as expected without using an intermediate result?
					// fails even with eval() after the first contract, so I use a separate variable for the first contract
					Eigen::Tensor<std::complex<double>, 3> intermediateResult = resTensor.contract(modGammas[s], contract_dim1);
					resTensor = intermediateResult.contract(daggerGammas[s], contract_dim);
				}

				const Eigen::Tensor<std::complex<double>, 0> t = resTensor.trace();
				std::complex<double> res = t(0);

				return res;
			}

			std::vector<bool> MeasureNoCollapse() override
			{
				const size_t nrQubits = gammas.size();
				return MeasureNoCollapse(static_cast<IndexType>(nrQubits - 1));
			}

			std::vector<bool> MeasureNoCollapse(const std::set<IndexType>& qubits) override
			{
				return MeasureNoCollapse(*qubits.crbegin());
			}

		private:
			void MultiplyModGammasWithLambdas(std::vector<GammaType>& modGammas, IndexType minQubit, IndexType nrSites)
			{
				const IndexType lastQubit = static_cast<IndexType>(lambdas.size());

				if (minQubit > 0)
				{
					// multiply with the left lambda as well
					const IndexType prev = minQubit - 1;
					const IndexType szl = modGammas[0].dimension(0);
					const IndexType szr = modGammas[0].dimension(2);

					for (IndexType r = 0; r < szr; ++r)
						for (IndexType p = 0; p < 2; ++p)
							for (IndexType l = 0; l < szl; ++l)
								modGammas[0](l, p, r) *= lambdas[prev][l];
				}

				// multiply with the right lambdas
				for (IndexType s = 0; s < nrSites; ++s)
				{
					const IndexType q = minQubit + s;
					if (q >= lastQubit) break; // no right lambdas for the last qubit

					const IndexType szl = modGammas[s].dimension(0);
					const IndexType szr = modGammas[s].dimension(2);

					for (IndexType r = 0; r < szr; ++r)
						for (IndexType p = 0; p < 2; ++p)
							for (IndexType l = 0; l < szl; ++l)
								modGammas[s](l, p, r) *= lambdas[q][r];
				}
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
					// TODO: optimize for some gates? The most important would be the swap gate, as it's applied a lot to move qubits near each other before applying other two qubit gates
					// something like:
					if (gate.isSwapGate())
					{
						Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);
						// apply the swap gate, optimized
						
						// NOTE: with proper data structures this could be way more optimized
						// instead of those generic eigen tensors for example, a more optimized tensor for this would hold Eigen matrices inside and
						// the physical indices would select among them
						// that way swapping will reduce to swapping pointers, it would be very fast
						// but I'm lazy so I use the generic implementation
						
						// qiskit aer does it, though

						// it's less important than it appears, because usually SVD will take the most time (unless severely limited or having some particularly easy circuits)

						for (IndexType j = 0; j < theta.dimension(3); ++j)
							for (IndexType i = 0; i < theta.dimension(0); ++i)
								std::swap(theta(i, 0, 1, j), theta(i, 1, 0, j));

						thetaMatrix = ReshapeTheta(theta);
					}
					else
					{
						const Eigen::TensorFixedSize<std::complex<double>, Eigen::Sizes<2, 2, 2, 2>> U = GetTwoQubitsGateTensor(gate, reversed);
						const Eigen::Tensor<std::complex<double>, 4> thetabar = ConstructThetaBar(qubit1, U);

						thetaMatrix = ReshapeThetaBar(thetabar);
					}
				}

				Eigen::JacobiSVD<MatrixClass> SVD;
				
				// This one is supposed to be faster on big matrices but from my tests it seems that it's not so accurate, so I'll use Jacobi for now 
				//Eigen::BDCSVD<MatrixClass> SVD(thetaMatrix.rows(), thetaMatrix.cols());

				if (limitEntanglement)
					SVD.setThreshold(singularValueThreshold);

				SVD.compute(thetaMatrix, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV);

				const MatrixClass& UmatrixFull = SVD.matrixU();
				const MatrixClass& VmatrixFull = SVD.matrixV();
				const LambdaType& SvaluesFull = SVD.singularValues();

				IndexType szm = SVD.nonzeroSingularValues(); // or SvaluesFull.size() for tests
				if (szm == 0) szm = 1; // Shouldn't happen (unless some big limit was put on 'zero')!

				const IndexType sz = limitSize ? std::min<IndexType>(chi, szm) : szm;


				const IndexType szl = qubit1 == 0 ? 1 : lambdas[qubit1 - 1].size();
				const IndexType szr = qubit2 == static_cast<IndexType>(lambdas.size()) ? 1 : lambdas[qubit2].size();

				assert(UmatrixFull.cols() == VmatrixFull.cols()); // for 'thin'
				assert(sz <= UmatrixFull.cols());


				const MatrixClass Umatrix = UmatrixFull.topLeftCorner(std::min<IndexType>(2 * szl, UmatrixFull.rows()), sz);
				const MatrixClass Vmatrix = VmatrixFull.topLeftCorner(std::min<IndexType>(2 * szr, VmatrixFull.rows()), sz).adjoint();

				// now set back lambdas and gammas

				lambdas[qubit1] = SvaluesFull.head(sz);
				assert(lambdas[qubit1][0] != 0.);
				lambdas[qubit1].normalize();

				SetNewGammas(Umatrix, Vmatrix, qubit1, qubit2, szl, sz, szr);
			}


			TwoQubitsGateTensor GetTwoQubitsGateTensor(const GateClass& gate, bool reversed) const
			{
				TwoQubitsGateTensor result;

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

				for (IndexType k = 0; k < sz; ++k)
					for (IndexType j = 0; j < 2; ++j)
					{
						const IndexType jszl = j * szl;
						for (IndexType i = 0; i < szl; ++i)
						{
							const IndexType jind = jszl + i;
							Utensor(i, j, k) = (jind < Umatrix.rows()) ? Umatrix(jind, k) : 0;
						}
					}

				for (IndexType k = 0; k < szr; ++k)
					for (IndexType j = 0; j < 2; ++j)
					{
						const IndexType jind = j * szr + k;
						for (IndexType i = 0; i < sz; ++i)
						{	
							Vtensor(i, j, k) = (jind < Vmatrix.cols()) ? Vmatrix(i, jind) : 0;
						}
					}

				gammas[qubit1] = Utensor;
				gammas[qubit2] = Vtensor;
			}

			inline void SetNewGammasSame(const MatrixClass& Umatrix, const MatrixClass& Vmatrix, IndexType qubit1, IndexType qubit2, IndexType sz)
			{
				Eigen::Tensor<std::complex<double>, 3> Utensor(sz, 2, sz);
				Eigen::Tensor<std::complex<double>, 3> Vtensor(sz, 2, sz);

				for (IndexType k = 0; k < sz; ++k)
					for (IndexType j = 0; j < 2; ++j)
					{
						const IndexType jchi = j * sz;
						const IndexType jchik = jchi + k;
						for (IndexType i = 0; i < sz; ++i)
						{
							const IndexType jind = jchi + i;
							Utensor(i, j, k) = (jind < Umatrix.rows()) ? Umatrix(jind, k) : 0;

							Vtensor(i, j, k) = (jchik < Vmatrix.cols()) ? Vmatrix(i, jchik) : 0;
						}
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
					for (IndexType k = 0; k < szr; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < sz; ++i)
								if (lambdas[qubit2][k] > std::numeric_limits<double>::epsilon() * 1E-10) gammas[qubit2](i, j, k) /= lambdas[qubit2][k];
								else gammas[qubit2](i, j, k) = 0;
				}
			}

			Eigen::Tensor<std::complex<double>, 4> ContractTwoQubits(IndexType qubit1)
			{
				const IndexType qubit2 = qubit1 + 1;

				static const Indexes product_dims_int{ IntIndexPair(2, 0) };

				assert(gammas[qubit1].dimension(2) == gammas[qubit2].dimension(0));

				const IndexType szl = gammas[qubit1].dimension(0);
				const IndexType sz = gammas[qubit1].dimension(2);
				const IndexType szr = gammas[qubit2].dimension(2);

				if (qubit1 != 0)
				{
					const IndexType prev = qubit1 - 1;
					for (IndexType k = 0; k < sz; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < szl; ++i)
								gammas[qubit1](i, j, k) *= lambdas[prev][i] * lambdas[qubit1][k];
				}
				else
				{
					for (IndexType k = 0; k < sz; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < szl; ++i)
								gammas[qubit1](i, j, k) *= lambdas[qubit1][k];
				}

				if (qubit2 != static_cast<IndexType>(lambdas.size()))
				{
					for (IndexType k = 0; k < szr; ++k)
						for (IndexType j = 0; j < 2; ++j)
							for (IndexType i = 0; i < sz; ++i)
								gammas[qubit2](i, j, k) *= lambdas[qubit2][k];
				}

				// contract first gamma with the lambda in the middle
				// the resulting tensor has three legs, 1 is the physical one

				// then

				// contract the result with the next gamma
				// the resulting tensor has four legs, 1 and 2 are the physical ones

				return gammas[qubit1].contract(gammas[qubit2], product_dims_int);
			}

			// for the two qubit gate - U is the gate tensor
			Eigen::Tensor<std::complex<double>, 4> ConstructThetaBar(IndexType qubit1, const TwoQubitsGateTensor& U)
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
				{
					const IndexType lsz1 = l * sz1;
					for (IndexType k = 0; k < 2; ++k)
					{
						const IndexType ksz0 = k * sz0;
						for (IndexType j = 0; j < sz1; ++j)
						{
							const IndexType lsz1j = lsz1 + j;
							for (IndexType i = 0; i < sz0; ++i)
								// 2, 0, 3, 1 - k & l from theta are the physical indexes
								thetaMatrix(ksz0 + i, lsz1j) = theta(i, j, k, l);
						}
					}
				}

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

				for (IndexType l = 0; l < sz3; ++l)
					for (IndexType k = 0; k < 2; ++k)
					{
						const IndexType ksz3l = k * sz3 + l;
						for (IndexType j = 0; j < 2; ++j)
						{
							const IndexType jsz0 = j * sz0;
							for (IndexType i = 0; i < sz0; ++i)
								// j & k from theta are the physical indexes
								thetaMatrix(jsz0 + i, ksz3l) = theta(i, j, k, l);
						}
					}

				return thetaMatrix;
			}

			std::vector<bool> MeasureNoCollapse(IndexType limit)
			{
				if (gammas.empty()) return {};

				const IndexType limit1 = limit + 1;
				std::vector<bool> res(limit1);

				Eigen::MatrixXcd mat;

				double totalProb = 1.;

				for (IndexType qubit = 0; qubit < limit1; ++qubit)
				{
					// 1. First, compute probability for measuring 0
					// zero matrix
					MatrixTensorType qubitMat = gammas[qubit].chip(0, 1);
					MatrixClass mq = Eigen::Map<const MatrixClass>(qubitMat.data(), qubitMat.dimension(0), qubitMat.dimension(1));
					MultiplyMatrixWithLambda(qubit, mq);
					if (qubit != 0)
						mq = mat * mq;

					// this is the probability of measuring all the qubits with the picked up values using the random number generator, up to this one
					// including a measured zero value for the current qubit 
					const double allProbability = mq.cwiseProduct(mq.conjugate()).sum().real();
					// to get the probability for the current qubit to be 0, we need to divide by the probability of measuring all the previous qubits
					const double prob0 = allProbability / totalProb;

					// 2. use that probability to measure the qubit
					const double rndVal = 1. - uniformZeroOne(rng);
					const bool zeroMeasured = rndVal < prob0;
					res[qubit] = !zeroMeasured;

					// accumulate the probability for measuring the current qubit to whatever was picked by using the random number generator
					totalProb *= zeroMeasured ? prob0 : 1. - prob0;

					// now update the matrix
					if (zeroMeasured) // no need to compute it again if 0 was measured, it was already computed above
						mat = std::move(mq);
					else
					{
						qubitMat = gammas[qubit].chip(1, 1);
						mq = Eigen::Map<const MatrixClass>(qubitMat.data(), qubitMat.dimension(0), qubitMat.dimension(1));
						MultiplyMatrixWithLambda(qubit, mq);
						if (qubit == 0)
							mat = std::move(mq);
						else
							mat = mat * mq;
					}
				}

				return res;
			}
		};

	}
}
