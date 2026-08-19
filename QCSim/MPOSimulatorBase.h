#pragma once

#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <map>
#include <utility>
#include <limits>
#include <algorithm>

#include <unsupported/Eigen/CXX11/Tensor>

#include "MPOSimulatorInterface.h"
#include "Operators.h"

namespace QC {

	namespace TensorNetworks {

		class MPOSimulatorBaseState : public MPOSimulatorStateInterface
		{
		public:
			MPOSimulatorBaseState() = default;
			MPOSimulatorBaseState(const MPOSimulatorBaseState&) = default;
			MPOSimulatorBaseState(MPOSimulatorBaseState&&) = default;
			MPOSimulatorBaseState& operator=(const MPOSimulatorBaseState&) = default;
			MPOSimulatorBaseState& operator=(MPOSimulatorBaseState&&) = default;
			virtual ~MPOSimulatorBaseState() = default;

			std::vector<MPOSimulatorInterface::LambdaType> lambdas;
			std::vector<MPOSimulatorInterface::TensorType> gammas;
		};

		// As for the MPS simulator, this base class is separated from the actual
		// simulator to reduce class complexity. It holds the data structures, the
		// initialization functions, the single qubit gate (which is local, so it does
		// not need the SVD machinery) and the 'observables': the trace, qubit
		// probabilities, basis-state matrix elements and the (costly) reconstruction of
		// the full density matrix, used for comparing against other simulators.
		//
		// Each site is a rank-4 tensor with leg order (leftBond, ket, bra, rightBond).
		// The represented density matrix is the Vidal-like product
		//        rho = Gamma[0] Lambda[0] Gamma[1] Lambda[1] ... Gamma[N-1]
		// where the physical indices are pairs (ket, bra).
		//
		// IMPORTANT difference from the MPS simulator: the singular values (lambdas)
		// are NOT renormalized after the SVD. The trace is linear in rho, so keeping
		// the raw singular values keeps Tr(rho) = 1 exactly under unitary evolution
		// when the SVD is not truncated (L2-normalizing them, as the MPS does to keep
		// <psi|psi> = 1, would instead rescale the trace).
		//
		// Compression caveat: limiting the bond dimension or dropping singular values
		// is ordinary operator-space MPO truncation. It minimizes a local SVD error, but
		// it does not enforce the density-matrix constraints globally. After truncation
		// the represented operator can have trace drift, small non-Hermitian components
		// and negative eigenvalues. Query functions normalize by the current trace where
		// appropriate, but they cannot make a truncated MPO positive semidefinite.
		class MPOSimulatorBase : public MPOSimulatorInterface
		{
		public:
			MPOSimulatorBase() = delete;

			MPOSimulatorBase(size_t N, unsigned int addseed = 0)
				: lambdas(N > 0 ? N - 1 : 0, LambdaType::Ones(1)), gammas(N, TensorType(1, 2, 2, 1))
			{
				if (N == 0)
					throw std::invalid_argument("MPOSimulator requires at least one qubit");

				for (auto& gamma : gammas)
					SetSiteToBasis(gamma, 0); // |0><0|

				const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
				std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
				rng.seed(seed);
			}

			MPOSimulatorBase(const MPOSimulatorBase&) = default;
			MPOSimulatorBase(MPOSimulatorBase&&) = default;
			MPOSimulatorBase& operator=(const MPOSimulatorBase&) = default;
			MPOSimulatorBase& operator=(MPOSimulatorBase&&) = default;

			size_t getNrQubits() const override
			{
				return gammas.size();
			}

			void Clear() override
			{
				const size_t szm1 = lambdas.size();
				for (size_t i = 0; i < szm1; ++i)
				{
					gammas[i].resize(1, 2, 2, 1);
					SetSiteToBasis(gammas[i], 0);

					lambdas[i].resize(1);
					lambdas[i](0) = 1.;
				}

				gammas[szm1].resize(1, 2, 2, 1);
				SetSiteToBasis(gammas[szm1], 0);
			}

			void InitOnesState() override
			{
				const size_t szm1 = lambdas.size();
				for (size_t i = 0; i < szm1; ++i)
				{
					gammas[i].resize(1, 2, 2, 1);
					SetSiteToBasis(gammas[i], 1);

					lambdas[i].resize(1);
					lambdas[i](0) = 1.;
				}

				gammas[szm1].resize(1, 2, 2, 1);
				SetSiteToBasis(gammas[szm1], 1);
			}

			void setToQubitState(IndexType q) override
			{
				Clear();
				if (q < 0 || q >= static_cast<IndexType>(gammas.size()))
					return;

				SetSiteToBasis(gammas[q], 1);
			}

			void setToBasisState(size_t State) override
			{
				const size_t NrBasisStates = gammas.size() > sizeof(size_t) * 8 ? 64 : (1ULL << gammas.size());
				if (State >= NrBasisStates) return;

				Clear();

				size_t pos = 0;
				while (State)
				{
					if (State & 1)
						SetSiteToBasis(gammas[pos], 1);

					State >>= 1;
					++pos;
				}
			}

			void setToBasisState(const std::vector<bool>& State) override
			{
				if (State.size() > gammas.size()) return;

				Clear();

				for (size_t i = 0; i < State.size(); ++i)
					if (State[i])
						SetSiteToBasis(gammas[i], 1);
			}

			// rho = sum_k prob_k |state_k><state_k|, with the basis states given as bit
			// masks (bit i is qubit i). The probabilities are normalized so Tr(rho) = 1.
			void setToMixtureOfBasisStates(const std::vector<std::pair<size_t, double>>& mixture) override
			{
				const size_t nrQubits = gammas.size();

				std::vector<std::pair<std::vector<bool>, double>> bitMixture;
				bitMixture.reserve(mixture.size());

				for (const auto& [state, prob] : mixture)
				{
					std::vector<bool> bits(nrQubits, false);
					size_t s = state;
					for (size_t i = 0; i < nrQubits; ++i)
					{
						bits[i] = (s & 1) == 1;
						s >>= 1;
					}
					bitMixture.emplace_back(std::move(bits), prob);
				}

				setToMixtureOfBasisStates(bitMixture);
			}

			// rho = sum_k prob_k |state_k><state_k|, with the basis states given as bit
			// vectors. The probabilities are normalized so Tr(rho) = 1.
			//
			// A statistical mixture of basis states is a diagonal density matrix, which
			// is represented exactly by a 'diagonal' MPO: the bond index labels the
			// mixture term k, each site tensor is diagonal in that bond index (and has
			// ket == bra == the basis bit of qubit i for term k), and the leftmost site
			// carries the (normalized) probabilities. The bond lambdas are all ones.
			void setToMixtureOfBasisStates(const std::vector<std::pair<std::vector<bool>, double>>& mixture) override
			{
				const size_t nrQubits = gammas.size();

				// merge duplicate basis states and accumulate their probabilities,
				// dropping non positive weights
				std::map<std::vector<bool>, double> merged;
				double total = 0.;
				for (const auto& [state, prob] : mixture)
				{
					if (prob <= 0.) continue;

					std::vector<bool> bits(nrQubits, false);
					for (size_t i = 0; i < nrQubits && i < state.size(); ++i)
						bits[i] = state[i];

					merged[bits] += prob;
					total += prob;
				}

				// fall back to |0...0><0...0| if there is nothing usable
				if (merged.empty() || total <= std::numeric_limits<double>::epsilon())
				{
					Clear();
					return;
				}

				const IndexType nrTerms = static_cast<IndexType>(merged.size());

				std::vector<std::vector<bool>> states;
				std::vector<double> probs;
				states.reserve(merged.size());
				probs.reserve(merged.size());
				for (const auto& [state, prob] : merged)
				{
					states.push_back(state);
					probs.push_back(prob / total); // normalize so Tr(rho) = 1
				}

				// the bond dimension is the number of distinct mixture terms; each lambda is all ones
				for (size_t b = 0; b < lambdas.size(); ++b)
					lambdas[b] = LambdaType::Ones(nrTerms);

				// each site is diagonal in the bond index k: gamma(k, bit_i(k), bit_i(k), k).
				// the first site folds in the probabilities so the chain contraction reproduces
				// rho = sum_k prob_k |state_k><state_k|
				for (size_t i = 0; i < nrQubits; ++i)
				{
					const IndexType leftBond = (i == 0) ? 1 : nrTerms;
					const IndexType rightBond = (i + 1 == nrQubits) ? 1 : nrTerms;

					TensorType gamma(leftBond, 2, 2, rightBond);
					gamma.setZero();

					for (IndexType k = 0; k < nrTerms; ++k)
					{
						const int bit = states[k][i] ? 1 : 0;
						const IndexType l = (i == 0) ? 0 : k;
						const IndexType r = (i + 1 == nrQubits) ? 0 : k;
						const std::complex<double> val = (i == 0) ? std::complex<double>(probs[k], 0.) : std::complex<double>(1., 0.);
						gamma(l, bit, bit, r) = val;
					}

					gammas[i] = std::move(gamma);
				}
			}

			void setLimitBondDimension(IndexType chival) override
			{
				limitSize = true;
				chi = chival;
			}

			void setLimitEntanglement(double svdThreshold) override
			{
				limitEntanglement = true;
				singularValueThreshold = svdThreshold;
			}

			void dontLimitBondDimension() override
			{
				limitSize = false;
			}

			void dontLimitEntanglement() override
			{
				limitEntanglement = false;
			}

			std::complex<double> Trace() const override
			{
				return ContractChain([this](IndexType q) { return SiteTraceMatrix(q); });
			}

			std::complex<double> ExpectationValue(const std::string& pauliString) const override
			{
				const size_t nrQubits = getNrQubits();
				if (pauliString.size() != nrQubits)
					throw std::invalid_argument("Pauli string length must match the number of qubits");

				// per site single qubit Pauli matrices (index by qubit); identity where the character is 'I'
				std::vector<MatrixClass> siteOps(nrQubits);
				for (size_t i = 0; i < nrQubits; ++i)
					siteOps[i] = PauliMatrixFromChar(pauliString[i]);

				const std::complex<double> num = ContractChain([this, &siteOps](IndexType q) {
					return SitePauliMatrix(q, siteOps[static_cast<size_t>(q)]);
				});

				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return 0.;

				return num / tr;
			}

			double GetProbability(IndexType qubit, bool zeroVal = true) const override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
					throw std::invalid_argument("Qubit index out of bounds");

				const int physIndex = zeroVal ? 0 : 1;

				const std::complex<double> num = ContractChain([this, qubit, physIndex](IndexType q) {
					return (q == qubit) ? SiteSelectMatrix(q, physIndex, physIndex) : SiteTraceMatrix(q);
				});

				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return 0.;

				return ClampProbability((num / tr).real());
			}

			// samples a full computational basis outcome from the density matrix populations without
			// collapsing the state, by sampling qubit after qubit conditioned on the previous outcomes
			std::unordered_map<IndexType, bool> MeasureNoCollapse() override
			{
				const IndexType n = static_cast<IndexType>(gammas.size());
				if (n == 0) return {};

				return MeasureNoCollapseUpTo(n - 1);
			}

			// samples a subset of qubits without collapsing the state. As with the MPS simulator, the
			// sampling proceeds along the chain up to the largest requested qubit; only the requested
			// qubits are reported. The map keys are the qubit indices, the values the outcomes.
			std::unordered_map<IndexType, bool> MeasureNoCollapse(const std::set<IndexType>& qubits) override
			{
				if (qubits.empty()) return {};

				const auto sampled = MeasureNoCollapseUpTo(*qubits.crbegin());

				std::unordered_map<IndexType, bool> res;
				for (const IndexType qubit : qubits)
				{
					const auto it = sampled.find(qubit);
					if (it != sampled.end())
						res[qubit] = it->second;
				}

				return res;
			}

			std::complex<double> getBasisStateMatrixElement(size_t row, size_t col) const override
			{
				const size_t nrQubits = getNrQubits();

				std::vector<bool> rowState(nrQubits, false);
				std::vector<bool> colState(nrQubits, false);

				for (size_t i = 0; i < nrQubits; ++i)
				{
					rowState[i] = (row & 1) == 1;
					colState[i] = (col & 1) == 1;
					row >>= 1;
					col >>= 1;
				}

				return getBasisStateMatrixElement(rowState, colState);
			}

			std::complex<double> getBasisStateMatrixElement(const std::vector<bool>& row, const std::vector<bool>& col) const override
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits == 0) return 0.;

				return ContractChain([this, &row, &col](IndexType q) {
					const int ket = (static_cast<size_t>(q) < row.size() && row[q]) ? 1 : 0;
					const int bra = (static_cast<size_t>(q) < col.size() && col[q]) ? 1 : 0;
					return SiteSelectMatrix(q, ket, bra);
				});
			}

			double getBasisStateProbability(size_t State) const override
			{
				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return 0.;

				return ClampProbability((getBasisStateMatrixElement(State, State) / tr).real());
			}

			double getBasisStateProbability(const std::vector<bool>& State) const override
			{
				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return 0.;

				return ClampProbability((getBasisStateMatrixElement(State, State) / tr).real());
			}

			// this is costly (it builds the full 2^N x 2^N matrix) and it's meant only
			// for comparing the results against other simulators, not for simulation
			MatrixClass getDensityMatrix() const override
			{
				const size_t sz = gammas.size();
				if (sz == 0) return {};
				if (sz > 13) throw std::runtime_error("Too many qubits to build the full density matrix");

				const size_t NrBasisStates = 1ULL << sz;
				MatrixClass rho(NrBasisStates, NrBasisStates);

				for (size_t r = 0; r < NrBasisStates; ++r)
					for (size_t c = 0; c < NrBasisStates; ++c)
						rho(r, c) = getBasisStateMatrixElement(r, c);

				return rho;
			}

			std::shared_ptr<MPOSimulatorStateInterface> getState() const override
			{
				auto state = std::make_shared<MPOSimulatorBaseState>();
				state->lambdas = lambdas;
				state->gammas = gammas;

				return state;
			}

			void setState(const std::shared_ptr<MPOSimulatorStateInterface>& state) override
			{
				if (!state) return;

				auto stateRef = std::static_pointer_cast<MPOSimulatorBaseState>(state);
				lambdas = stateRef->lambdas;
				gammas = stateRef->gammas;
			}

			void setStateDestructive(std::shared_ptr<MPOSimulatorStateInterface>& state)
			{
				if (!state) return;

				auto stateRef = std::static_pointer_cast<MPOSimulatorBaseState>(state);
				lambdas.swap(stateRef->lambdas);
				gammas.swap(stateRef->gammas);
			}

			void print() const override
			{
				for (size_t i = 0; i < gammas.size(); ++i)
				{
					std::cout << std::endl << "Gamma " << i << " (leftBond, ket, bra, rightBond):" << std::endl;
					PrintGamma(i);
					if (i < lambdas.size())
						std::cout << "Lambda " << i << ":\n" << lambdas[i] << std::endl;
				}
			}

			std::vector<IndexType> getBondDimensions() const
			{
				std::vector<IndexType> dims(lambdas.size());
				for (size_t i = 0; i < lambdas.size(); ++i)
					dims[i] = lambdas[i].size();
				return dims;
			}

			void printBondDimensions() const
			{
				std::cout << "Bond dimensions: ";
				for (const auto& lambda : lambdas)
					std::cout << lambda.size() << " ";
				std::cout << std::endl;
			}

		protected:
			static double ClampProbability(double probability)
			{
				constexpr double tolerance = 1E-12;
				if (probability < 0. && probability > -tolerance) return 0.;
				if (probability > 1. && probability < 1. + tolerance) return 1.;

				return probability;
			}

			static double ValidMeasurementProbability(double probability)
			{
				constexpr double tolerance = 1E-9;
				if (probability < -tolerance || probability > 1. + tolerance)
					std::cerr << "Invalid measurement probability produced by the MPO state" << std::endl;

				return ClampProbability(std::clamp(probability, 0., 1.));
			}

			// samples qubits [0, limit] one after the other, each conditioned on the previous outcomes,
			// without collapsing the state; the remaining qubits are traced out. Shared by both the
			// full and the subset MeasureNoCollapse overloads.
			std::unordered_map<IndexType, bool> MeasureNoCollapseUpTo(IndexType limit)
			{
				std::unordered_map<IndexType, bool> res;

				const IndexType n = static_cast<IndexType>(gammas.size());
				if (n == 0) return res;

				if (limit < 0) limit = 0;
				if (limit >= n) limit = n - 1;

				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					return res;

				std::vector<int> outcomes(n, 0);

				// joint probability that qubits [0, upTo) have the outcomes stored in 'outcomes'
				// (the remaining qubits are traced out)
				auto jointProbability = [this, &outcomes, tr](IndexType upTo) -> double {
					const std::complex<double> num = ContractChain([this, &outcomes, upTo](IndexType q) {
						if (q < upTo)
							return SiteSelectMatrix(q, outcomes[q], outcomes[q]);
						return SiteTraceMatrix(q);
					});

					return ClampProbability((num / tr).real());
				};

				double priorProb = 1.;
				for (IndexType qubit = 0; qubit <= limit; ++qubit)
				{
					// conditional probability of measuring 0 on the current qubit given the previous outcomes
					outcomes[qubit] = 0;
					const double joint0 = jointProbability(qubit + 1);
					const double prob0 = ValidMeasurementProbability(priorProb > std::numeric_limits<double>::epsilon() ? joint0 / priorProb : 0.);

					const double rndVal = 1. - uniformZeroOne(rng);
					const bool zeroMeasured = rndVal < prob0;

					outcomes[qubit] = zeroMeasured ? 0 : 1;
					res[qubit] = !zeroMeasured;

					priorProb *= zeroMeasured ? prob0 : 1. - prob0;
				}

				return res;
			}

			static void SetSiteToBasis(TensorType& gamma, int basis)
			{
				gamma.setZero();
				gamma(0, basis, basis, 0) = 1.; // |basis><basis|
			}

			void PrintGamma(size_t i) const
			{
				assert(i < gammas.size());

				const auto& g = gammas[i];
				for (IndexType ket = 0; ket < 2; ++ket)
					for (IndexType bra = 0; bra < 2; ++bra)
					{
						std::cout << "ket " << ket << ", bra " << bra << " matrix:" << std::endl;
						for (IndexType l = 0; l < g.dimension(0); ++l)
						{
							for (IndexType r = 0; r < g.dimension(3); ++r)
								std::cout << g(l, ket, bra, r) << " ";
							std::cout << std::endl;
						}
					}
			}

			// matrix (leftBond x rightBond) obtained by tracing out the physical index of a site (ket == bra, summed)
			MatrixClass SiteTraceMatrix(IndexType q) const
			{
				const auto& g = gammas[q];
				const IndexType L = g.dimension(0);
				const IndexType R = g.dimension(3);

				MatrixClass m = MatrixClass::Zero(L, R);
				for (IndexType s = 0; s < 2; ++s)
					for (IndexType r = 0; r < R; ++r)
						for (IndexType l = 0; l < L; ++l)
							m(l, r) += g(l, s, s, r);

				return m;
			}

			// matrix (leftBond x rightBond) obtained by selecting fixed ket and bra physical indices of a site
			MatrixClass SiteSelectMatrix(IndexType q, int ket, int bra) const
			{
				const auto& g = gammas[q];
				const IndexType L = g.dimension(0);
				const IndexType R = g.dimension(3);

				MatrixClass m(L, R);
				for (IndexType r = 0; r < R; ++r)
					for (IndexType l = 0; l < L; ++l)
						m(l, r) = g(l, ket, bra, r);

				return m;
			}

			// single qubit Pauli matrix from a character in a Pauli string ('I', 'X', 'Y', 'Z')
			static MatrixClass PauliMatrixFromChar(char c)
			{
				MatrixClass m = MatrixClass::Zero(2, 2);
				switch (toupper(static_cast<unsigned char>(c)))
				{
				case 'I':
					m(0, 0) = 1.; m(1, 1) = 1.;
					break;
				case 'X':
					m(0, 1) = 1.; m(1, 0) = 1.;
					break;
				case 'Y':
					m(0, 1) = std::complex<double>(0., -1.); m(1, 0) = std::complex<double>(0., 1.);
					break;
				case 'Z':
					m(0, 0) = 1.; m(1, 1) = -1.;
					break;
				default:
					throw std::invalid_argument("Invalid operator in the Pauli string");
				}

				return m;
			}

			// matrix (leftBond x rightBond) obtained by contracting a single qubit operator P into
			// the physical legs of a site: m(l, r) = sum_{ket,bra} g(l, ket, bra, r) P(bra, ket).
			// With P = Identity this reduces to SiteTraceMatrix.
			MatrixClass SitePauliMatrix(IndexType q, const MatrixClass& P) const
			{
				const auto& g = gammas[q];
				const IndexType L = g.dimension(0);
				const IndexType R = g.dimension(3);

				MatrixClass m = MatrixClass::Zero(L, R);
				for (IndexType ket = 0; ket < 2; ++ket)
					for (IndexType bra = 0; bra < 2; ++bra)
					{
						const std::complex<double> p = P(bra, ket);
						if (p == std::complex<double>(0., 0.)) continue;

						for (IndexType r = 0; r < R; ++r)
							for (IndexType l = 0; l < L; ++l)
								m(l, r) += g(l, ket, bra, r) * p;
					}

				return m;
			}

			// contracts the whole chain, picking at each site a (leftBond x rightBond) matrix
			// supplied by 'siteMatrix', with the lambdas folded in on the bonds in between
			template<typename SiteMatrixFunc> std::complex<double> ContractChain(SiteMatrixFunc siteMatrix) const
			{
				const size_t n = gammas.size();
				if (n == 0) return 0.;

				MatrixClass res = siteMatrix(0);

				for (size_t q = 1; q < n; ++q)
				{
					const auto& lam = lambdas[q - 1];
					for (IndexType c = 0; c < res.cols(); ++c)
					{
						const double l = c < lam.size() ? lam[c] : 0.;
						for (IndexType r = 0; r < res.rows(); ++r)
							res(r, c) *= l;
					}

					const MatrixClass sm = siteMatrix(static_cast<IndexType>(q));
					res = (res * sm).eval();
				}

				return res(0, 0);
			}

			void NormalizeByTrace()
			{
				const std::complex<double> tr = Trace();
				if (std::abs(tr) < std::numeric_limits<double>::epsilon())
					throw std::runtime_error("Cannot normalize an MPO state with zero trace");

				ScaleSite(0, 1. / tr);
			}

			void ScaleSite(IndexType q, std::complex<double> factor)
			{
				auto& g = gammas[q];
				for (IndexType r = 0; r < g.dimension(3); ++r)
					for (IndexType bra = 0; bra < 2; ++bra)
						for (IndexType ket = 0; ket < 2; ++ket)
							for (IndexType l = 0; l < g.dimension(0); ++l)
								g(l, ket, bra, r) *= factor;
			}

			static void AbsorbLambdasIntoGammas(std::vector<LambdaType>& stateLambdas, std::vector<TensorType>& stateGammas)
			{
				for (size_t q = 0; q < stateLambdas.size(); ++q)
				{
					auto& gamma = stateGammas[q];
					const auto& lambda = stateLambdas[q];

					for (IndexType r = 0; r < gamma.dimension(3); ++r)
					{
						const double lval = r < lambda.size() ? lambda[r] : 0.;
						for (IndexType bra = 0; bra < 2; ++bra)
							for (IndexType ket = 0; ket < 2; ++ket)
								for (IndexType l = 0; l < gamma.dimension(0); ++l)
									gamma(l, ket, bra, r) *= lval;
					}

					stateLambdas[q] = LambdaType::Ones(gamma.dimension(3));
				}
			}

			static void AddState(std::vector<LambdaType>& lhsLambdas, std::vector<TensorType>& lhsGammas, const std::vector<LambdaType>& rhsLambdas, const std::vector<TensorType>& rhsGammas)
			{
				assert(lhsGammas.size() == rhsGammas.size());
				assert(lhsLambdas.size() == rhsLambdas.size());

				std::vector<LambdaType> rhsLambdasAbsorbed = rhsLambdas;
				std::vector<TensorType> rhsGammasAbsorbed = rhsGammas;

				AbsorbLambdasIntoGammas(lhsLambdas, lhsGammas);
				AbsorbLambdasIntoGammas(rhsLambdasAbsorbed, rhsGammasAbsorbed);

				const size_t n = lhsGammas.size();
				if (n == 1)
				{
					lhsGammas[0] = (lhsGammas[0] + rhsGammasAbsorbed[0]).eval();
					return;
				}

				std::vector<TensorType> sumGammas(n);
				std::vector<LambdaType> sumLambdas(n - 1);

				for (size_t q = 0; q < n; ++q)
				{
					const auto& lg = lhsGammas[q];
					const auto& rg = rhsGammasAbsorbed[q];

					const IndexType lL = lg.dimension(0);
					const IndexType lR = lg.dimension(3);
					const IndexType rL = rg.dimension(0);
					const IndexType rR = rg.dimension(3);
					const IndexType sumL = q == 0 ? 1 : lL + rL;
					const IndexType sumR = q + 1 == n ? 1 : lR + rR;

					TensorType gamma(sumL, 2, 2, sumR);
					gamma.setZero();

					CopyGammaBlock(lg, gamma, 0, 0);
					CopyGammaBlock(rg, gamma, q == 0 ? 0 : lL, q + 1 == n ? 0 : lR);

					sumGammas[q] = std::move(gamma);

					if (q + 1 < n)
						sumLambdas[q] = LambdaType::Ones(sumR);
				}

				lhsGammas = std::move(sumGammas);
				lhsLambdas = std::move(sumLambdas);
			}

			static void CopyGammaBlock(const TensorType& source, TensorType& target, IndexType leftOffset, IndexType rightOffset)
			{
				for (IndexType r = 0; r < source.dimension(3); ++r)
					for (IndexType bra = 0; bra < 2; ++bra)
						for (IndexType ket = 0; ket < 2; ++ket)
							for (IndexType l = 0; l < source.dimension(0); ++l)
								target(leftOffset + l, ket, bra, rightOffset + r) = source(l, ket, bra, r);
			}

			// applies a single qubit local operator to a site tensor: rho_site -> A rho_site A^dagger
			// the operator A is contracted with the ket leg, the conjugate operator A* with the bra leg
			static void ApplySingleQubitGate(TensorType& gamma, const GateClass& gate)
			{
				ApplySingleQubitGate(gamma, gate.getRawOperatorMatrix());
			}

			static void ApplySingleQubitGate(TensorType& gamma, const MatrixClass& opMat)
			{
				static const Indexes contractKet{ IntIndexPair(1, 1) }; // gamma ket (dim 1) with U column (dim 1)
				static const Indexes contractBra{ IntIndexPair(1, 1) }; // intermediate bra (dim 1) with U* column (dim 1)
				static const std::array<int, 4> permute{ 0, 2, 3, 1 };

				const Eigen::TensorMap<const OneQubitGateTensor> Utensor(opMat.data(), opMat.rows(), opMat.cols());
				const OneQubitGateTensor Uconj = Utensor.conjugate();

				// (leftBond, ket, bra, rightBond) x A(ket', ket) over ket -> (leftBond, bra, rightBond, ket')
				const TensorType tmp = gamma.contract(Utensor, contractKet);
				// (leftBond, bra, rightBond, ket') x A*(bra', bra) over bra -> (leftBond, rightBond, ket', bra')
				const TensorType res = tmp.contract(Uconj, contractBra);
				// shuffle to (leftBond, ket', bra', rightBond)
				gamma = res.shuffle(permute);
			}

		protected:
			bool limitSize = false;
			bool limitEntanglement = false;
			IndexType chi = 10; // if limitSize is true
			double singularValueThreshold = 0.; // if limitEntanglement is true

			std::vector<LambdaType> lambdas;
			std::vector<TensorType> gammas;

			std::mt19937_64 rng;
			std::uniform_real_distribution<double> uniformZeroOne{ 0, 1 };

			const Operators::ZeroProjection<MatrixClass> zeroProjection;
			const Operators::OneProjection<MatrixClass> oneProjection;
		};

	}

}
