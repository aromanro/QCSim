#pragma once

#include "MPOSimulatorBase.h"

#include <iostream>

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

				bool truncated = false;
				for (IndexType qubit1 = 0; qubit1 < static_cast<IndexType>(lambdas.size()); ++qubit1)
				{
					if (lambdas[qubit1].size() <= chi) continue;

					// contract the two neighbour sites (no gate is applied), then re-split with a truncated SVD
					const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
					const MatrixClass thetaMatrix = ReshapeTheta(theta);

					truncated = DecomposeAndSetGammas(thetaMatrix, qubit1, qubit1 + 1) || truncated;
				}

				ApplyPostTruncationPatches(truncated);
			}

			void ReCanonicalize() override
			{
				// Gauge only: do not apply user-requested chi / singular-value cuts.
				for (IndexType qubit1 = 0; qubit1 < static_cast<IndexType>(lambdas.size()); ++qubit1)
				{
					const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
					const MatrixClass thetaMatrix = ReshapeTheta(theta);
					DecomposeAndSetGammas(thetaMatrix, qubit1, qubit1 + 1, false);
				}
				for (IndexType qubit1 = static_cast<IndexType>(lambdas.size()) - 1; qubit1 >= 0; --qubit1)
				{
					const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
					const MatrixClass thetaMatrix = ReshapeTheta(theta);
					DecomposeAndSetGammas(thetaMatrix, qubit1, qubit1 + 1, false);
				}
			}

			void Hermitize() override
			{
				const bool wasApplying = applyingPostTruncationPatches;
				applyingPostTruncationPatches = true;

				std::vector<LambdaType> adjointLambdas = lambdas;
				std::vector<TensorType> adjointGammas;
				adjointGammas.reserve(gammas.size());
				for (const auto& gamma : gammas)
					adjointGammas.emplace_back(AdjointSite(gamma));

				AddState(lambdas, gammas, adjointLambdas, adjointGammas);
				ScaleSite(0, 0.5);

				ReCanonicalize();
				if (limitSize || limitEntanglement)
				{
					for (IndexType qubit1 = 0; qubit1 < static_cast<IndexType>(lambdas.size()); ++qubit1)
					{
						if (limitSize && !limitEntanglement && lambdas[qubit1].size() <= chi)
							continue;

						const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
						DecomposeAndSetGammas(ReshapeTheta(theta), qubit1, qubit1 + 1);
					}
				}

				applyingPostTruncationPatches = wasApplying;
				if (restoreTraceAfterTruncation && !wasApplying)
					RestoreTraceIfSafe();
			}

			void ApplyOperator(const Gates::AppliedGate<MatrixClass>& op) override
			{
				const size_t operatorQubits = ValidateOperator(op);
				ValidateQubits(operatorQubits, op.getQubit1(), op.getQubit2());

				ApplyValidatedOperator(op.getRawOperatorMatrix(), operatorQubits,
					static_cast<IndexType>(op.getQubit1()),
					operatorQubits > 1 ? static_cast<IndexType>(op.getQubit2()) : 0);
			}

			// three qubit operators not supported, two qubit operators need to act on adjacent qubits
			// (the higher level MPOSimulator swaps qubits around to satisfy this)
			void ApplyOperator(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				const size_t operatorQubits = ValidateOperator(op);
				ValidateQubits(operatorQubits, qubit, controllingQubit1);

				ApplyValidatedOperator(op.getRawOperatorMatrix(), operatorQubits, qubit, controllingQubit1);
			}

			void ApplyOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) override
			{
				for (const auto& op : ops)
					ApplyOperator(op);
			}

			void ApplyOperatorAndNormalize(const Gates::AppliedGate<MatrixClass>& op) override
			{
				const size_t operatorQubits = ValidateOperator(op);
				ValidateQubits(operatorQubits, op.getQubit1(), op.getQubit2());

				ApplyValidatedOperatorAndNormalize(op.getRawOperatorMatrix(), operatorQubits,
					static_cast<IndexType>(op.getQubit1()),
					operatorQubits > 1 ? static_cast<IndexType>(op.getQubit2()) : 0);
			}

			void ApplyOperatorAndNormalize(const GateClass& op, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				const size_t operatorQubits = ValidateOperator(op);
				ValidateQubits(operatorQubits, qubit, controllingQubit1);

				ApplyValidatedOperatorAndNormalize(op.getRawOperatorMatrix(), operatorQubits, qubit, controllingQubit1);
			}

			void ApplyKrausOperators(const std::vector<Gates::AppliedGate<MatrixClass>>& ops) override
			{
				if (ops.empty()) return;

				const size_t qubit = ops.front().getQubit1();
				const size_t controllingQubit1 = ops.front().getQubit2();
				const size_t operatorQubits = ValidateOperator(ops.front());

				// Validate the complete channel before changing any tensor.
				for (const auto& op : ops)
				{
					if (ValidateOperator(op) != operatorQubits)
						throw std::invalid_argument("All Kraus operators need to have the same arity");
				}

				ValidateQubits(operatorQubits, qubit, controllingQubit1);
				for (const auto& op : ops)
					if (op.getQubit1() != qubit ||
						(operatorQubits > 1 && op.getQubit2() != controllingQubit1))
						throw std::invalid_argument("All Kraus operators need to act on the same qubits");

				ApplyKrausOperatorsImpl(ops, operatorQubits,
					static_cast<IndexType>(qubit),
					operatorQubits > 1 ? static_cast<IndexType>(controllingQubit1) : 0);
			}

			void ApplyKrausOperators(const std::vector<MatrixClass>& ops, IndexType qubit, IndexType controllingQubit1 = 0) override
			{
				if (ops.empty()) return;

				const size_t operatorQubits = ValidateOperatorMatrix(ops.front());
				for (const auto& op : ops)
					if (ValidateOperatorMatrix(op) != operatorQubits)
						throw std::invalid_argument("All Kraus operators need to have the same arity");

				ValidateQubits(operatorQubits, qubit, controllingQubit1);
				ApplyKrausOperatorsImpl(ops, operatorQubits, qubit, controllingQubit1);
			}

			// false if measured 0, true if measured 1
			bool MeasureQubit(IndexType qubit) override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
					throw std::invalid_argument("Qubit index out of bounds");

				const double prob0 = ValidMeasurementProbability(GetProbability(qubit, true));
				const bool zeroMeasured = uniformZeroOne(rng) < prob0;

				// for a density matrix, the post measurement state is P rho P (suitably normalized)
				// P acts only on the measured site, so this is a purely local operation - no need to
				// propagate the collapse along the chain like the MPS simulator does
				const MatrixClass& projMat = zeroMeasured ? zeroProjection.getRawOperatorMatrix() : oneProjection.getRawOperatorMatrix();
				TensorType collapsedGamma = gammas[qubit];
				ApplySingleQubitGate(collapsedGamma, projMat);

				// Normalize by the actual post-projection trace, rather than by the probability
				// normalized against the input trace. This also handles deliberately unnormalized MPOs.
				const std::complex<double> postTrace = ContractChain([this, qubit, &collapsedGamma](IndexType q) {
					return q == qubit ? TraceSiteTensor(collapsedGamma) : SiteTraceMatrix(q);
				});
				if (!IsFinite(postTrace) || std::abs(postTrace) <= std::numeric_limits<double>::epsilon())
					throw std::runtime_error("Cannot collapse an MPO state onto a zero-probability outcome");

				ScaleTensor(collapsedGamma, 1. / postTrace);
				gammas[qubit] = std::move(collapsedGamma);

				return !zeroMeasured;
			}

			std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType>& qubits) override
			{
				for (const IndexType qubit : qubits)
					if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
						throw std::invalid_argument("Qubit index out of bounds");

				std::unordered_map<IndexType, bool> res;
				for (const auto qubit : qubits)
					res[qubit] = MeasureQubit(qubit);

				return res;
			}

		private:
			static bool IsFinite(const std::complex<double>& value)
			{
				return std::isfinite(value.real()) && std::isfinite(value.imag());
			}

			static size_t ValidateOperatorMatrix(const MatrixClass& matrix)
			{
				if (matrix.rows() <= 0 || matrix.cols() <= 0 || matrix.rows() != matrix.cols())
					throw std::invalid_argument("Operator matrix must be a non-empty square matrix");
				if (matrix.rows() != 2 && matrix.rows() != 4)
					throw std::invalid_argument("MPO operators must have dimensions 2x2 or 4x4");

				for (IndexType col = 0; col < matrix.cols(); ++col)
					for (IndexType row = 0; row < matrix.rows(); ++row)
						if (!IsFinite(matrix(row, col)))
							throw std::invalid_argument("Operator matrix elements must be finite");

				return matrix.rows() == 2 ? 1 : 2;
			}

			static size_t ValidateOperator(const GateClass& op)
			{
				const size_t matrixQubits = ValidateOperatorMatrix(op.getRawOperatorMatrix());
				if (op.getQubitsNumber() != matrixQubits)
					throw std::invalid_argument("Operator matrix dimensions do not match its declared arity");
				return matrixQubits;
			}

			void ValidateQubits(size_t operatorQubits, IndexType qubit, IndexType controllingQubit1) const
			{
				const IndexType nrQubits = static_cast<IndexType>(gammas.size());
				if (qubit < 0 || qubit >= nrQubits ||
					(operatorQubits == 2 && (controllingQubit1 < 0 || controllingQubit1 >= nrQubits)))
					throw std::invalid_argument("Qubit index out of bounds");
				if (operatorQubits == 2 && qubit == controllingQubit1)
					throw std::invalid_argument("Two qubit operators require distinct qubits");
				if (operatorQubits == 2 && std::abs(qubit - controllingQubit1) != 1)
					throw std::invalid_argument("Two qubit operators need to act on adjacent qubits");
			}

			void ValidateQubits(size_t operatorQubits, size_t qubit, size_t controllingQubit1) const
			{
				const size_t nrQubits = gammas.size();
				if (qubit >= nrQubits || (operatorQubits == 2 && controllingQubit1 >= nrQubits))
					throw std::invalid_argument("Qubit index out of bounds");
				if (operatorQubits == 2 && qubit == controllingQubit1)
					throw std::invalid_argument("Two qubit operators require distinct qubits");
				if (operatorQubits == 2 && qubit + 1 != controllingQubit1 && controllingQubit1 + 1 != qubit)
					throw std::invalid_argument("Two qubit operators need to act on adjacent qubits");
			}

			static const MatrixClass& GetOperatorMatrix(const MatrixClass& op)
			{
				return op;
			}

			static const MatrixClass& GetOperatorMatrix(const GateClass& op)
			{
				return op.getRawOperatorMatrix();
			}

			static MatrixClass TraceSiteTensor(const TensorType& gamma)
			{
				const IndexType L = gamma.dimension(0);
				const IndexType R = gamma.dimension(3);
				MatrixClass result = MatrixClass::Zero(L, R);

				for (IndexType state = 0; state < 2; ++state)
					for (IndexType r = 0; r < R; ++r)
						for (IndexType l = 0; l < L; ++l)
							result(l, r) += gamma(l, state, state, r);

				return result;
			}

			static void ScaleTensor(TensorType& gamma, std::complex<double> factor)
			{
				for (IndexType r = 0; r < gamma.dimension(3); ++r)
					for (IndexType bra = 0; bra < 2; ++bra)
						for (IndexType ket = 0; ket < 2; ++ket)
							for (IndexType l = 0; l < gamma.dimension(0); ++l)
								gamma(l, ket, bra, r) *= factor;
			}

			void ApplyValidatedOperator(const MatrixClass& op, size_t operatorQubits, IndexType qubit, IndexType controllingQubit1)
			{
				if (operatorQubits == 1)
					ApplySingleQubitGate(gammas[qubit], op);
				else
					ApplyTwoQubitGate(op, qubit, controllingQubit1);
			}

			void ApplyValidatedOperatorAndNormalize(const MatrixClass& op, size_t operatorQubits, IndexType qubit, IndexType controllingQubit1)
			{
				const auto originalLambdas = lambdas;
				const auto originalGammas = gammas;

				try
				{
					ApplyValidatedOperator(op, operatorQubits, qubit, controllingQubit1);
					const std::complex<double> trace = Trace();
					if (!IsFinite(trace) || std::abs(trace) <= std::numeric_limits<double>::epsilon())
						throw std::runtime_error("Cannot normalize an MPO state with zero or non-finite trace");
					ScaleSite(0, 1. / trace);
				}
				catch (...)
				{
					lambdas = originalLambdas;
					gammas = originalGammas;
					throw;
				}
			}

			template<class OperatorsContainer> void ApplyKrausOperatorsImpl(const OperatorsContainer& ops, size_t operatorQubits, IndexType qubit, IndexType controllingQubit1)
			{
				CheckKrausCompleteness(ops);

				if (operatorQubits == 1)
					ApplySingleQubitChannel(ops, qubit);
				else
					ApplyTwoQubitChannel(ops, qubit, controllingQubit1);
			}

			template<class OperatorsContainer> void CheckKrausCompleteness(const OperatorsContainer& ops) const
			{
				if (krausCompletenessCheck == KrausCompletenessCheck::Ignore) return;

				const MatrixClass& first = GetOperatorMatrix(ops.front());
				MatrixClass sum = MatrixClass::Zero(first.rows(), first.cols());
				for (const auto& op : ops)
				{
					const MatrixClass& K = GetOperatorMatrix(op);
					sum.noalias() += K.adjoint() * K;
				}

				const double residual = (sum - MatrixClass::Identity(sum.rows(), sum.cols())).norm();
				if (residual <= krausCompletenessTolerance) return;

				if (krausCompletenessCheck == KrausCompletenessCheck::Warn)
				{
					std::cerr << "Kraus operators do not satisfy the completeness relation (residual "
						<< residual << ")" << std::endl;
					return;
				}

				throw std::invalid_argument("Kraus operators do not satisfy the completeness relation");
			}

			template<class OperatorsContainer> void ApplySingleQubitChannel(const OperatorsContainer& ops, IndexType qubit)
			{
				const TensorType& original = gammas[qubit];
				const IndexType L = original.dimension(0);
				const IndexType R = original.dimension(3);
				TensorType result(L, 2, 2, R);
				result.setZero();

				for (const auto& op : ops)
				{
					const MatrixClass& E = GetOperatorMatrix(op);
					for (IndexType ketOut = 0; ketOut < 2; ++ketOut)
						for (IndexType braOut = 0; braOut < 2; ++braOut)
							for (IndexType ketIn = 0; ketIn < 2; ++ketIn)
								for (IndexType braIn = 0; braIn < 2; ++braIn)
								{
									const std::complex<double> factor = E(ketOut, ketIn) * std::conj(E(braOut, braIn));
									if (factor == std::complex<double>(0., 0.)) continue;

									for (IndexType r = 0; r < R; ++r)
										for (IndexType l = 0; l < L; ++l)
											result(l, ketOut, braOut, r) += factor * original(l, ketIn, braIn, r);
								}
				}

				gammas[qubit] = std::move(result);
			}

			template<class OperatorsContainer> void ApplyTwoQubitChannel(const OperatorsContainer& ops, IndexType qubit, IndexType controllingQubit1)
			{
				IndexType qubit1 = controllingQubit1;
				IndexType qubit2 = qubit;
				bool reversed = false;
				if (qubit1 > qubit2)
				{
					std::swap(qubit1, qubit2);
					reversed = true;
				}

				const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
				const IndexType L = theta.dimension(0);
				const IndexType R = theta.dimension(5);
				Eigen::Tensor<std::complex<double>, 6> result(L, R, 2, 2, 2, 2);
				result.setZero();

				for (const auto& op : ops)
				{
					const TwoQubitsGateTensor E = GetTwoQubitsGateTensor(GetOperatorMatrix(op), reversed);
					for (IndexType ket1Out = 0; ket1Out < 2; ++ket1Out)
						for (IndexType ket2Out = 0; ket2Out < 2; ++ket2Out)
							for (IndexType bra1Out = 0; bra1Out < 2; ++bra1Out)
								for (IndexType bra2Out = 0; bra2Out < 2; ++bra2Out)
									for (IndexType ket1In = 0; ket1In < 2; ++ket1In)
										for (IndexType ket2In = 0; ket2In < 2; ++ket2In)
											for (IndexType bra1In = 0; bra1In < 2; ++bra1In)
												for (IndexType bra2In = 0; bra2In < 2; ++bra2In)
												{
													const std::complex<double> factor = E(ket1Out, ket2Out, ket1In, ket2In) *
														std::conj(E(bra1Out, bra2Out, bra1In, bra2In));
													if (factor == std::complex<double>(0., 0.)) continue;

													for (IndexType r = 0; r < R; ++r)
														for (IndexType l = 0; l < L; ++l)
															result(l, r, ket1Out, ket2Out, bra1Out, bra2Out) += factor *
																theta(l, ket1In, bra1In, ket2In, bra2In, r);
												}
				}

				const MatrixClass thetaMatrix = ReshapeThetaBar(result);
				ApplyPostTruncationPatches(DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2));
			}

			void ApplyTwoQubitGate(const MatrixClass& gate, IndexType qubit, IndexType controllingQubit1)
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

				ApplyPostTruncationPatches(DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2));
			}

			void ApplyPostTruncationPatches(bool truncated)
			{
				if (!truncated || applyingPostTruncationPatches) return;
				if (!hermitizeAfterTruncation && !restoreTraceAfterTruncation) return;

				applyingPostTruncationPatches = true;
				if (hermitizeAfterTruncation)
					Hermitize();
				if (restoreTraceAfterTruncation)
					RestoreTraceIfSafe();
				applyingPostTruncationPatches = false;
			}

			// SVD the (already built) theta matrix and write back the two new site tensors and the
			// lambda in between. Shared by two-qubit gates, Trim, and ReCanonicalize. User-requested
			// chi / entanglement cuts are applied only when applyUserCompression is true (the
			// default); ReCanonicalize passes false so it is a gauge restore.
			// Returns true when a user-requested cut actually dropped singular values.
			bool DecomposeAndSetGammas(const MatrixClass& thetaMatrix, IndexType qubit1, IndexType qubit2, bool applyUserCompression = true)
			{
#ifdef USE_FAST_SVD
				const bool computeWithJacobi = thetaMatrix.rows() < blockSizeLimit && thetaMatrix.cols() < blockSizeLimit;
#endif

				// Eigen's default rank threshold is close enough to machine epsilon that
				// roundoff-only singular values can survive and later make the Vidal
				// pseudoinverse ill-conditioned. Always filter those out first with the fixed
				// scale-relative floor, regardless of the requested truncation mode or whether
				// user-requested compression (limitEntanglement) is enabled at all - compression
				// selection, below, is a separate step applied on top of this floor.
#ifdef USE_FAST_SVD
				if (computeWithJacobi)
#endif
					jacobiSVD.setThreshold(numericalRankThreshold);
#ifdef USE_FAST_SVD
				else
					SVD.setThreshold(numericalRankThreshold);
#endif

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

				const IndexType floorRank = computeWithJacobi ? jacobiSVD.rank() : SVD.rank();
#else
				const MatrixClass& UmatrixFull = jacobiSVD.matrixU();
				const MatrixClass& VmatrixFull = jacobiSVD.matrixV();
				const LambdaType& SvaluesFull = jacobiSVD.singularValues();

				const IndexType floorRank = jacobiSVD.rank();
#endif
				// If user-requested compression is enabled, further reduce the rank according to
				// the configured truncation mode, applied on top of the already floor-filtered
				// (still descending-sorted) singular values.
				IndexType szm = (applyUserCompression && limitEntanglement) ?
					ComputeCompressedRank(SvaluesFull, floorRank, truncationMode, singularValueThreshold) : floorRank;

				if (szm == 0) szm = 1;

				const IndexType sz = (applyUserCompression && limitSize) ? std::min<IndexType>(chi, szm) : szm;
				const bool truncated = applyUserCompression && sz < floorRank;

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
				return truncated;
			}

			// See MPSSimulatorImpl::ComputeCompressedRank for the full explanation of both modes -
			// identical logic, duplicated here per this codebase's per-file convention (no shared
			// header between the MPS and MPO impls today). Always keeps at least the largest
			// singular value.
			static IndexType ComputeCompressedRank(const LambdaType& sortedDescendingSVs, IndexType rank, TruncationMode mode, double threshold)
			{
				if (rank <= 1) return rank;

				if (mode == TruncationMode::RelativeToMax)
				{
					const double effectiveThreshold = std::max(threshold, numericalRankThreshold);
					const double premultipliedThreshold = std::max(effectiveThreshold * sortedDescendingSVs[0], std::numeric_limits<double>::min());
					IndexType i = rank - 1;
					while (i >= 0 && sortedDescendingSVs[i] < premultipliedThreshold) --i;
					return i + 1;
				}

				const double effectiveThreshold = std::max(threshold, numericalRankThresholdDiscardedWeight);
				const double total = sortedDescendingSVs.head(rank).squaredNorm();
				if (total <= 0.) return rank;

				double discarded = 0.;
				IndexType keep = rank;
				for (IndexType i = rank - 1; i > 0; --i)
				{
					const double sq = sortedDescendingSVs[i] * sortedDescendingSVs[i];
					if ((discarded + sq) / total >= effectiveThreshold) break;

					discarded += sq;
					keep = i;
				}

				return keep;
			}

			static TwoQubitsGateTensor GetTwoQubitsGateTensor(const MatrixClass& gateMat, bool reversed)
			{
				TwoQubitsGateTensor result;

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

				if (qubit1 != 0)
				{
					const IndexType prev = qubit1 - 1;
					// Treat numerically null Schmidt sectors as zero when applying the
					// pseudoinverse. Dividing by SVD noise near machine epsilon makes the
					// Vidal tensors ill-conditioned and can amplify roundoff to visible
					// density-matrix errors after subsequent gates.
					const double threshold = numericalRankThreshold * lambdas[prev].maxCoeff();
					for (IndexType m = 0; m < sz; ++m)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType l = 0; l < L; ++l)
									if (lambdas[prev][l] > threshold) gammas[qubit1](l, ket, bra, m) /= lambdas[prev][l];
									else gammas[qubit1](l, ket, bra, m) = 0;
				}

				if (qubit2 != static_cast<IndexType>(lambdas.size()))
				{
					const double threshold = numericalRankThreshold * lambdas[qubit2].maxCoeff();
					for (IndexType r = 0; r < R; ++r)
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType m = 0; m < sz; ++m)
									if (lambdas[qubit2][r] > threshold) gammas[qubit2](m, ket, bra, r) /= lambdas[qubit2][r];
									else gammas[qubit2](m, ket, bra, r) = 0;
				}
			}

#ifdef USE_FAST_SVD
			constexpr static IndexType blockSizeLimit = 64;
			Eigen::BDCSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV> SVD;
#endif
			constexpr static double numericalRankThreshold = 1E-12;
			constexpr static double krausCompletenessTolerance = 1E-10;
			// See MPSSimulatorImpl::numericalRankThresholdDiscardedWeight for why this is squared
			// rather than reusing numericalRankThreshold as-is.
			constexpr static double numericalRankThresholdDiscardedWeight = numericalRankThreshold * numericalRankThreshold;
			Eigen::JacobiSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV> jacobiSVD;
			bool applyingPostTruncationPatches = false;
		};

	}
}
