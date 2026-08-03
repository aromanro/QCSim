#pragma once

#include "MPOSimulatorBase.h"

#define USE_FAST_SVD 1

namespace QC {

	namespace TensorNetworks {

		// The actual MPO simulator working on adjacent qubits.
		// Two qubit gates need to act on adjacent qubits; the MPOSimulator decorator
		// (see MPOSimulator.h) takes care of swapping non adjacent qubits together.
		// If compression is enabled, the SVD split below truncates the operator-space
		// MPO, not a positivity-preserving purification. Therefore the result remains
		// an approximation to the density matrix, but it is not guaranteed to be a valid
		// density matrix after truncation.
		class MPOSimulatorImpl : public MPOSimulatorBase
		{
		public:
			friend class MPOSimulator;

			MPOSimulatorImpl(size_t N, unsigned int addseed = 0)
				: MPOSimulatorBase(N, addseed)
			{
#ifdef USE_FAST_SVD
				SVD.setSwitchSize(blockSizeLimit); // lower sizes will use Jacobi
#endif
			}

			void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate) override
			{
				ApplyOperator(gate);
			}

			// three qubit gates not supported, two qubit gates need to act on adjacent qubits
			// (the higher level MPOSimulator swaps qubits around to satisfy this)
			void ApplyGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				ApplyOperator(gate, qubit, controllingQubit1);
			}

			void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates) override
			{
				for (const auto& gate : gates)
					ApplyGate(gate);
			}

			// Walks over the chain and, wherever a bond dimension is bigger than the currently set
			// bond dimension limit, contracts the two neighbour sites and re-splits them with a
			// truncated SVD to bring the bond dimension down to the limit. No gate is applied.
			// Useful after lowering the bond dimension limit with setLimitBondDimension once the
			// bonds have already grown beyond the new limit.
			void Trim() override
			{
				if (!limitSize) return; // nothing to trim against

				for (IndexType qubit1 = 0; qubit1 < static_cast<IndexType>(lambdas.size()); ++qubit1)
				{
					if (lambdas[qubit1].size() <= chi) continue;

					// contract the two neighbour sites (no gate is applied), then re-split with a truncated SVD
					const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
					const MatrixClass thetaMatrix = ReshapeTheta(theta);

					DecomposeAndSetGammas(thetaMatrix, qubit1, qubit1 + 1);
				}
			}

			void ApplyOperator(const Gates::AppliedGate<MatrixClass>& op) override
			{
				ApplyOperator(op, op.getQubit1(), op.getQubit2());
			}

			// three qubit operators not supported, two qubit operators need to act on adjacent qubits
			// (the higher level MPOSimulator swaps qubits around to satisfy this)
			void ApplyOperator(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (op.getQubitsNumber() > 2) throw std::invalid_argument("Three qubit operators not supported");
				else if ((qubit < 0 || qubit >= static_cast<IndexType>(gammas.size())) || (op.getQubitsNumber() == 2 && (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(gammas.size()))))
					throw std::invalid_argument("Qubit index out of bounds");
				else if (op.getQubitsNumber() == 2 && std::abs(static_cast<int>(qubit) - static_cast<int>(controllingQubit1)) != 1)
					throw std::invalid_argument("Two qubit operators need to act on adjacent qubits");

				if (op.getQubitsNumber() == 1)
					ApplySingleQubitGate(gammas[qubit], op);
				else
					ApplyTwoQubitGate(op, qubit, controllingQubit1);
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
				ApplyOperator(op, qubit, controllingQubit1);
				NormalizeByTrace();
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
				std::vector<Gates::AppliedGate<MatrixClass>> appliedOps;
				appliedOps.reserve(ops.size());
				for (const auto& op : ops)
					appliedOps.emplace_back(op, qubit, controllingQubit1);

				ApplyKrausOperatorsImpl(appliedOps, qubit, controllingQubit1);
			}

			// false if measured 0, true if measured 1
			bool MeasureQubit(IndexType qubit) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
					throw std::invalid_argument("Qubit index out of bounds");

				const double rndVal = 1. - uniformZeroOne(rng);
				const double prob0 = ValidMeasurementProbability(GetProbability(qubit, true));
				const bool zeroMeasured = rndVal < prob0;

				// for a density matrix, the post measurement state is P rho P (suitably normalized)
				// P acts only on the measured site, so this is a purely local operation - no need to
				// propagate the collapse along the chain like the MPS simulator does
				const MatrixClass& projMat = zeroMeasured ? zeroProjection.getRawOperatorMatrix() : oneProjection.getRawOperatorMatrix();
				ApplySingleQubitGate(gammas[qubit], projMat);

				// Tr(P rho P) is the outcome probability; rescaling a single site tensor rescales the whole trace
				const double prob = zeroMeasured ? prob0 : 1. - prob0;
				if (prob > std::numeric_limits<double>::epsilon())
					ScaleSite(qubit, 1. / prob);

				return !zeroMeasured;
			}

			std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) override
			{
				std::unordered_map<IndexType, bool> res;
				for (const auto qubit : qubits)
					res[qubit] = MeasureQubit(qubit);

				return res;
			}

		private:
			template<class OperatorsContainer> void ApplyKrausOperatorsImpl(const OperatorsContainer& ops, IndexType qubit, IndexType controllingQubit1)
			{
				if (ops.empty()) return;

				const auto originalLambdas = lambdas;
				const auto originalGammas = gammas;
				std::vector<LambdaType> resultLambdas;
				std::vector<TensorType> resultGammas;
				bool hasResult = false;

				for (const auto& op : ops)
				{
					lambdas = originalLambdas;
					gammas = originalGammas;

					ApplyOperator(op, qubit, controllingQubit1);

					if (!hasResult)
					{
						resultLambdas = lambdas;
						resultGammas = gammas;
						hasResult = true;
					}
					else
					{
						AddState(resultLambdas, resultGammas, lambdas, gammas);
					}
				}

				lambdas = std::move(resultLambdas);
				gammas = std::move(resultGammas);
			}

			static double ValidMeasurementProbability(double probability)
			{
				constexpr double tolerance = 1E-9;
				if (probability < -tolerance || probability > 1. + tolerance)
					std::cerr << "Invalid measurement probability produced by the MPO state" << std::endl;

				return ClampProbability(std::clamp(probability, 0., 1.));
			}

			void ApplyTwoQubitGate(const GateClass& gate, IndexType qubit, IndexType controllingQubit1)
			{
				// contract the tensors for the two qubits (folding in the corresponding lambdas)
				// apply the gate on the kets and its conjugate on the bras (rho -> U rho U^dagger)
				// SVD the result to separate the two qubit tensors and the new lambda in between.
				// If the SVD is truncated, this is an operator-space approximation only; it does
				// not preserve positivity of the represented density matrix.

				IndexType qubit1 = controllingQubit1;
				IndexType qubit2 = qubit;
				bool reversed = false;

				if (qubit1 > qubit2)
				{
					std::swap(qubit1, qubit2);
					reversed = true;
				}

				const TwoQubitsGateTensor U = GetTwoQubitsGateTensor(gate, reversed);
				const Eigen::Tensor<std::complex<double>, 6> thetaBar = ConstructThetaBar(qubit1, U);

				// (4 * leftBond) x (4 * rightBond) matrix, the physical dimension per site is 4 = 2 (ket) x 2 (bra)
				const MatrixClass thetaMatrix = ReshapeThetaBar(thetaBar);

				DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2);
			}

			// SVD the (already built) theta matrix, truncate it according to the bond dimension / entanglement
			// limits and write back the two new site tensors and the lambda in between.
			// Shared by the two qubit gate application and by Trim.
			void DecomposeAndSetGammas(const MatrixClass& thetaMatrix, IndexType qubit1, IndexType qubit2)
			{
#ifdef USE_FAST_SVD
				const bool computeWithJacobi = thetaMatrix.rows() < blockSizeLimit && thetaMatrix.cols() < blockSizeLimit;
#endif

				if (limitEntanglement)
				{
#ifdef USE_FAST_SVD
					if (computeWithJacobi)
#endif
						jacobiSVD.setThreshold(singularValueThreshold);
#ifdef USE_FAST_SVD
					else
						SVD.setThreshold(singularValueThreshold);
#endif
				}

#ifdef USE_FAST_SVD
				if (computeWithJacobi)
#endif
					jacobiSVD.compute(thetaMatrix);
#ifdef USE_FAST_SVD
				else
					SVD.compute(thetaMatrix);

				const MatrixClass& UmatrixFull = computeWithJacobi ? jacobiSVD.matrixU() : SVD.matrixU();
				const MatrixClass& VmatrixFull = computeWithJacobi ? jacobiSVD.matrixV() : SVD.matrixV();
				const LambdaType& SvaluesFull = computeWithJacobi ? jacobiSVD.singularValues() : SVD.singularValues();

				// or SvaluesFull.size() for tests
				IndexType szm = limitEntanglement ? (computeWithJacobi ? jacobiSVD.rank() : SVD.rank()) : (computeWithJacobi ? jacobiSVD.nonzeroSingularValues() : SVD.nonzeroSingularValues());
#else
				const MatrixClass& UmatrixFull = jacobiSVD.matrixU();
				const MatrixClass& VmatrixFull = jacobiSVD.matrixV();
				const LambdaType& SvaluesFull = jacobiSVD.singularValues();

				IndexType szm = limitEntanglement ? jacobiSVD.rank() : jacobiSVD.nonzeroSingularValues(); // or SvaluesFull.size() for tests
#endif

				if (szm == 0) szm = 1;

				const IndexType sz = limitSize ? std::min<IndexType>(chi, szm) : szm;

				const IndexType L = qubit1 == 0 ? 1 : lambdas[qubit1 - 1].size();
				const IndexType R = qubit2 == static_cast<IndexType>(lambdas.size()) ? 1 : lambdas[qubit2].size();

				assert(UmatrixFull.cols() == VmatrixFull.cols()); // for 'thin'
				assert(sz <= UmatrixFull.cols());
				assert(4 * L == UmatrixFull.rows());
				assert(4 * R == VmatrixFull.rows());

				const MatrixClass Umatrix = UmatrixFull.topLeftCorner(UmatrixFull.rows(), sz);
				const MatrixClass Vmatrix = VmatrixFull.topLeftCorner(VmatrixFull.rows(), sz).adjoint();

				// the new lambda is the raw singular values: NOT normalized, to keep Tr(rho) preserved when untruncated
				lambdas[qubit1] = SvaluesFull.head(sz);

				SetNewGammas(Umatrix, Vmatrix, qubit1, qubit2, L, sz, R);
			}

			static TwoQubitsGateTensor GetTwoQubitsGateTensor(const GateClass& gate, bool reversed)
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

			static IndexType CountNonZeroSingularValues(const LambdaType& singularValues)
			{
				IndexType count = 0;
				while (count < singularValues.size() && singularValues[count] > 0.)
					++count;

				return count;
			}

			// contracts the two adjacent sites (folding in the lambdas) into a rank-6 tensor
			// with leg order (leftBond, ket1, bra1, ket2, bra2, rightBond)
			Eigen::Tensor<std::complex<double>, 6> ContractTwoQubits(IndexType qubit1)
			{
				const IndexType qubit2 = qubit1 + 1;

				static const Indexes contractMid{ IntIndexPair(3, 0) };

				assert(gammas[qubit1].dimension(3) == gammas[qubit2].dimension(0));

				const IndexType L = gammas[qubit1].dimension(0);
				const IndexType mid = gammas[qubit1].dimension(3);
				const IndexType R = gammas[qubit2].dimension(3);

				if (qubit1 != 0)
				{
					const IndexType prev = qubit1 - 1;
					for (IndexType r = 0; r < mid; ++r)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType l = 0; l < L; ++l)
									gammas[qubit1](l, ket, bra, r) *= lambdas[prev][l] * lambdas[qubit1][r];
				}
				else
				{
					for (IndexType r = 0; r < mid; ++r)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType l = 0; l < L; ++l)
									gammas[qubit1](l, ket, bra, r) *= lambdas[qubit1][r];
				}

				if (qubit2 != static_cast<IndexType>(lambdas.size()))
				{
					for (IndexType r = 0; r < R; ++r)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType l = 0; l < mid; ++l)
									gammas[qubit2](l, ket, bra, r) *= lambdas[qubit2][r];
				}

				return gammas[qubit1].contract(gammas[qubit2], contractMid);
			}

			// applies the two qubit gate U on the kets and its conjugate on the bras
			// returns a rank-6 tensor with leg order (leftBond, rightBond, ket1', ket2', bra1', bra2')
			Eigen::Tensor<std::complex<double>, 6> ConstructThetaBar(IndexType qubit1, const TwoQubitsGateTensor& U)
			{
				const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);

				const TwoQubitsGateTensor Uconj = U.conjugate();

				// contract the kets: theta(ket1=1, ket2=3) with U(in1=2, in2=3)
				// (leftBond, ket1, bra1, ket2, bra2, rightBond) -> (leftBond, bra1, bra2, rightBond, ket1', ket2')
				static const Indexes2 ketDims{ IntIndexPair(1, 2), IntIndexPair(3, 3) };
				const Eigen::Tensor<std::complex<double>, 6> theta2 = theta.contract(U, ketDims);

				// contract the bras: theta2(bra1=1, bra2=2) with U*(in1=2, in2=3)
				// (leftBond, bra1, bra2, rightBond, ket1', ket2') -> (leftBond, rightBond, ket1', ket2', bra1', bra2')
				static const Indexes2 braDims{ IntIndexPair(1, 2), IntIndexPair(2, 3) };
				return theta2.contract(Uconj, braDims);
			}

			// reshapes the rank-6 theta into a (4 * leftBond) x (4 * rightBond) matrix for the SVD
			// the combined physical index per site is p = ket * 2 + bra
			static MatrixClass ReshapeThetaBar(const Eigen::Tensor<std::complex<double>, 6>& theta)
			{
				// theta dims: 0=leftBond, 1=rightBond, 2=ket1, 3=ket2, 4=bra1, 5=bra2
				const IndexType L = theta.dimension(0);
				const IndexType R = theta.dimension(1);

				assert(L > 0);
				assert(R > 0);

				MatrixClass thetaMatrix(4 * L, 4 * R);

				for (IndexType ket1 = 0; ket1 < 2; ++ket1)
					for (IndexType bra1 = 0; bra1 < 2; ++bra1)
					{
						const IndexType p1L = (ket1 * 2 + bra1) * L;
						for (IndexType ket2 = 0; ket2 < 2; ++ket2)
							for (IndexType bra2 = 0; bra2 < 2; ++bra2)
							{
								const IndexType p2R = (ket2 * 2 + bra2) * R;
								for (IndexType r = 0; r < R; ++r)
									for (IndexType l = 0; l < L; ++l)
										thetaMatrix(p1L + l, p2R + r) = theta(l, r, ket1, ket2, bra1, bra2);
							}
					}

				return thetaMatrix;
			}

			// reshapes the raw two-site contraction (leg order leftBond, ket1, bra1, ket2, bra2, rightBond,
			// as produced by ContractTwoQubits) into a (4 * leftBond) x (4 * rightBond) matrix for the SVD.
			// Same target layout as ReshapeThetaBar (physical index per site p = ket * 2 + bra), but for the
			// no-gate leg order used by Trim.
			static MatrixClass ReshapeTheta(const Eigen::Tensor<std::complex<double>, 6>& theta)
			{
				// theta dims: 0=leftBond, 1=ket1, 2=bra1, 3=ket2, 4=bra2, 5=rightBond
				const IndexType L = theta.dimension(0);
				const IndexType R = theta.dimension(5);

				assert(L > 0);
				assert(R > 0);

				MatrixClass thetaMatrix(4 * L, 4 * R);

				for (IndexType ket1 = 0; ket1 < 2; ++ket1)
					for (IndexType bra1 = 0; bra1 < 2; ++bra1)
					{
						const IndexType p1L = (ket1 * 2 + bra1) * L;
						for (IndexType ket2 = 0; ket2 < 2; ++ket2)
							for (IndexType bra2 = 0; bra2 < 2; ++bra2)
							{
								const IndexType p2R = (ket2 * 2 + bra2) * R;
								for (IndexType r = 0; r < R; ++r)
									for (IndexType l = 0; l < L; ++l)
										thetaMatrix(p1L + l, p2R + r) = theta(l, ket1, bra1, ket2, bra2, r);
							}
					}

				return thetaMatrix;
			}

			void SetNewGammas(const MatrixClass& Umatrix, const MatrixClass& Vmatrix, IndexType qubit1, IndexType qubit2, IndexType L, IndexType sz, IndexType R)
			{
				// Umatrix is (4L x sz), Vmatrix is the adjoint of V, that is (sz x 4R)
				TensorType Utensor(L, 2, 2, sz);
				TensorType Vtensor(sz, 2, 2, R);

				for (IndexType ket = 0; ket < 2; ++ket)
					for (IndexType bra = 0; bra < 2; ++bra)
					{
						const IndexType pL = (ket * 2 + bra) * L;
						for (IndexType m = 0; m < sz; ++m)
							for (IndexType l = 0; l < L; ++l)
							{
								const IndexType idx = pL + l;
								Utensor(l, ket, bra, m) = (idx < Umatrix.rows()) ? Umatrix(idx, m) : std::complex<double>(0.);
							}
					}

				for (IndexType ket = 0; ket < 2; ++ket)
					for (IndexType bra = 0; bra < 2; ++bra)
					{
						const IndexType pR = (ket * 2 + bra) * R;
						for (IndexType r = 0; r < R; ++r)
						{
							const IndexType idx = pR + r;
							for (IndexType m = 0; m < sz; ++m)
								Vtensor(m, ket, bra, r) = (idx < Vmatrix.cols()) ? Vmatrix(m, idx) : std::complex<double>(0.);
						}
					}

				gammas[qubit1] = Utensor;
				gammas[qubit2] = Vtensor;

				DivideGammasWithLambdas(qubit1, qubit2, L, sz, R);
			}

			void DivideGammasWithLambdas(IndexType qubit1, IndexType qubit2, IndexType L, IndexType sz, IndexType R)
			{
				assert(gammas[qubit1].dimension(0) == L);
				assert(gammas[qubit1].dimension(3) == sz);
				assert(gammas[qubit2].dimension(0) == sz);
				assert(gammas[qubit2].dimension(3) == R);

				constexpr double thr = std::numeric_limits<double>::epsilon() * 1E-10;

				if (qubit1 != 0)
				{
					const IndexType prev = qubit1 - 1;
					for (IndexType m = 0; m < sz; ++m)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType l = 0; l < L; ++l)
									if (lambdas[prev][l] > thr) gammas[qubit1](l, ket, bra, m) /= lambdas[prev][l];
									else gammas[qubit1](l, ket, bra, m) = 0;
				}

				if (qubit2 != static_cast<IndexType>(lambdas.size()))
				{
					for (IndexType r = 0; r < R; ++r)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType m = 0; m < sz; ++m)
									if (lambdas[qubit2][r] > thr) gammas[qubit2](m, ket, bra, r) /= lambdas[qubit2][r];
									else gammas[qubit2](m, ket, bra, r) = 0;
				}
			}

#ifdef USE_FAST_SVD
			constexpr static IndexType blockSizeLimit = 64;
			Eigen::BDCSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV> SVD;
#endif
			Eigen::JacobiSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV> jacobiSVD;
		};

	}
}
