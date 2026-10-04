#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <typeinfo>
#include <utility>
#include <vector>

#include <unsupported/Eigen/CXX11/Tensor>

#include "MPOSimulatorInterface.h"
#include "Operators.h"
#include "SingleThreaded.h"

namespace QC
{

namespace TensorNetworks
{

class MPOSimulatorBaseState : public MPOSimulatorStateInterface
{
  public:
    MPOSimulatorBaseState() = default;
    MPOSimulatorBaseState(const MPOSimulatorBaseState &) = default;
    MPOSimulatorBaseState(MPOSimulatorBaseState &&) = default;
    MPOSimulatorBaseState &operator=(const MPOSimulatorBaseState &) = default;
    MPOSimulatorBaseState &operator=(MPOSimulatorBaseState &&) = default;
    virtual ~MPOSimulatorBaseState() = default;

    std::vector<MPOSimulatorInterface::LambdaType> lambdas;
    std::vector<MPOSimulatorInterface::TensorType> gammas;
    // rho = 2^scaleExponent * B[0] B[1] ... B[N-1] (see MPOSimulatorBase)
    int64_t scaleExponent = 0;
};

// As for the MPS simulator, this base class is separated from the actual
// simulator to reduce class complexity. It holds the data structures, the
// initialization functions, the single qubit gate (which is local, so it does
// not need the SVD machinery) and the 'observables': the trace, qubit
// probabilities, basis-state matrix elements and the (costly) reconstruction of
// the full density matrix, used for comparing against other simulators.
//
// Each site is a rank-4 tensor with leg order (leftBond, ket, bra, rightBond).
// As for the MPS simulator, the 'gammas' are NOT the Vidal Gamma tensors, they hold
// B[i] = Gamma[i] Lambda[i] (Hastings' form, see M. B. Hastings, J. Math. Phys. 50, 095207 (2009)),
// the last site has no right lambda, so there B = Gamma. The represented density matrix is the product
//        rho = B[0] B[1] ... B[N-1]   ( = Gamma[0] Lambda[0] Gamma[1] Lambda[1] ... Gamma[N-1] )
// where the physical indices are pairs (ket, bra).
// The lambdas (the operator space singular values on the bonds) are needed only as the left
// environment in the two site SVD, so the two qubit gates never divide by them.
//
// IMPORTANT difference from the MPS simulator: the singular values (lambdas) and the B tensors
// are NOT renormalized after the SVD. The trace is linear in rho, so keeping
// the raw singular values keeps Tr(rho) = 1 to numerical precision under
// unitary evolution when user-requested compression is disabled
// (L2-normalizing the lambdas, as the MPS does to keep <psi|psi> = 1, would
// instead rescale the trace).
// Local nonunitary operations leave the lambdas stale. Queries still contract the exact
// represented operator, but SVD rank decisions require valid environments. The implementation
// uses local QR transport to supply orthonormal environments for a two-site update.
// Between full ReCanonicalize calls the tensors may use a mixed canonical gauge;
// then only the most recently split bond has current Schmidt weights.
//
// Scale: the represented operator is 2^scaleExponent * B[0] ... B[N-1]. The exponent stays 0
// while the magnitudes are ordinary; the QR transport and the two-site split move a power of
// two out of a tensor only when its largest component leaves [2^-64, 2^64]. Without that, the
// transport piles the norm of every visited site into one tensor (2^-n/2 for a maximally mixed
// chain), and Householder QR silently drops columns whose squared norm underflows: the trace
// of a 1030-qubit maximally mixed chain became 2^-9 after ReCanonicalize. The lambdas only
// weight SVD inputs, so they are rescaled freely and do not enter the exponent. Readers
// rescale their contractions the same way and restore the exponent in the result. Powers of
// two are exact, so nothing changes for ordinary magnitudes.
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
    struct CanonicalMetadata
    {
        bool valid = true;
        IndexType first = 0, last = 0;
    };

    MPOSimulatorBase() = delete;

    MPOSimulatorBase(size_t N, unsigned int addseed = 0)
        : lambdas(N > 0 ? N - 1 : 0, LambdaType::Ones(1)), gammas(N, TensorType(1, 2, 2, 1))
    {
        if (N == 0)
            throw std::invalid_argument("MPOSimulator requires at least one qubit");

        for (auto &gamma : gammas)
            SetSiteToBasis(gamma, 0); // |0><0|

        if (addseed == 0)
        {
            std::random_device rdl;
            addseed = rdl();
        }

        const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
        std::seed_seq seed{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(seed);
    }

    MPOSimulatorBase(const MPOSimulatorBase &) = default;
    MPOSimulatorBase(MPOSimulatorBase &&) = default;
    MPOSimulatorBase &operator=(const MPOSimulatorBase &) = default;
    MPOSimulatorBase &operator=(MPOSimulatorBase &&) = default;

    void SetSeed(uint64_t theSeed)
    {
        std::seed_seq seed{uint32_t(theSeed & 0xffffffff), uint32_t(theSeed >> 32)};
        rng.seed(seed);
    }

    size_t getNrQubits() const override
    {
        return gammas.size();
    }

    void Clear() override
    {
        InvalidateSamplingCache();
        scaleExponent = 0;
        canonicalFormValid = true;
        centerFirst = centerLast = 0;
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
        InvalidateSamplingCache();
        scaleExponent = 0;
        canonicalFormValid = true;
        centerFirst = centerLast = 0;
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
        if (q < 0 || q >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        Clear();

        SetSiteToBasis(gammas[q], 1);
    }

    void setToBasisState(size_t State) override
    {
        constexpr size_t stateBits = std::numeric_limits<size_t>::digits;
        if (gammas.size() < stateBits && State >= (size_t{1} << gammas.size()))
            throw std::invalid_argument("Basis state is outside the MPO register");

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

    void setToBasisState(const std::vector<bool> &State) override
    {
        if (State.size() > gammas.size())
            throw std::invalid_argument("Basis state has more bits than the MPO register");

        Clear();

        for (size_t i = 0; i < State.size(); ++i)
            if (State[i])
                SetSiteToBasis(gammas[i], 1);
    }

    // rho = sum_k prob_k |state_k><state_k|, with the basis states given as bit
    // masks (bit i is qubit i). The probabilities are normalized so Tr(rho) = 1.
    void setToMixtureOfBasisStates(const std::vector<std::pair<size_t, double>> &mixture) override
    {
        const size_t nrQubits = gammas.size();
        constexpr size_t stateBits = std::numeric_limits<size_t>::digits;

        std::vector<std::pair<std::vector<bool>, double>> bitMixture;
        bitMixture.reserve(mixture.size());

        for (const auto &[state, prob] : mixture)
        {
            if (!std::isfinite(prob))
                throw std::invalid_argument("Mixture weights must be finite");
            if (prob <= 0.)
                continue;
            if (nrQubits < stateBits && state >= (size_t{1} << nrQubits))
                continue;

            std::vector<bool> bits(nrQubits, false);
            size_t s = state;
            for (size_t i = 0; i < nrQubits && i < stateBits; ++i)
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
    void setToMixtureOfBasisStates(const std::vector<std::pair<std::vector<bool>, double>> &mixture) override
    {
        const size_t nrQubits = gammas.size();

        // Validate the complete input before building replacement storage. Oversized
        // bit vectors describe out-of-range states and, like out-of-range integer
        // states, are ignored. Shorter vectors are zero-extended.
        double maxWeight = 0.;
        for (const auto &[state, prob] : mixture)
        {
            if (!std::isfinite(prob))
                throw std::invalid_argument("Mixture weights must be finite");
            if (prob > 0. && state.size() <= nrQubits)
                maxWeight = std::max(maxWeight, prob);
        }

        if (maxWeight <= 0.)
            throw std::invalid_argument("Mixture must contain at least one valid positive weight");

        // Scale by the largest usable weight before accumulation. This keeps both
        // subnormal-only and near-DBL_MAX mixtures normalizable without underflowing
        // the total to zero or overflowing it to infinity.
        std::map<std::vector<bool>, double> merged;
        double total = 0.;
        for (const auto &[state, prob] : mixture)
        {
            if (prob <= 0. || state.size() > nrQubits)
                continue;

            std::vector<bool> bits(nrQubits, false);
            for (size_t i = 0; i < state.size(); ++i)
                bits[i] = state[i];

            const double scaledWeight = prob / maxWeight;
            merged[bits] += scaledWeight;
            total += scaledWeight;
        }

        if (merged.empty() || !std::isfinite(total) || total <= 0.)
            throw std::invalid_argument("Mixture must contain at least one valid positive weight");

        const IndexType nrTerms = static_cast<IndexType>(merged.size());

        std::vector<std::vector<bool>> states;
        std::vector<double> probs;
        states.reserve(merged.size());
        probs.reserve(merged.size());
        for (const auto &[state, prob] : merged)
        {
            states.push_back(state);
            probs.push_back(prob / total); // normalize so Tr(rho) = 1
        }

        // Build the complete replacement before mutating the simulator, so invalid
        // input and allocation failures leave the existing state untouched.
        std::vector<LambdaType> newLambdas(lambdas.size(), LambdaType::Ones(nrTerms));
        std::vector<TensorType> newGammas;
        newGammas.reserve(nrQubits);

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
                const std::complex<double> val =
                    (i == 0) ? std::complex<double>(probs[k], 0.) : std::complex<double>(1., 0.);
                gamma(l, bit, bit, r) = val;
            }

            newGammas.emplace_back(std::move(gamma));
        }

        lambdas.swap(newLambdas);
        gammas.swap(newGammas);
        scaleExponent = 0;
        InvalidateCanonicalForm();
    }

    void setLimitBondDimension(IndexType chival) override
    {
        if (chival <= 0)
            throw std::invalid_argument("Bond dimension limit must be positive");

        limitSize = true;
        chi = chival;
    }

    void setLimitEntanglement(double svdThreshold) override
    {
        if (!std::isfinite(svdThreshold) || svdThreshold < 0.)
            throw std::invalid_argument("Singular-value threshold must be finite and non-negative");

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

    bool setTruncationMode(TruncationMode mode) override
    {
        switch (mode)
        {
        case TruncationMode::RelativeToMax:
        case TruncationMode::DiscardedWeight:
            truncationMode = mode;
            return true;
        default:
            throw std::invalid_argument("Unrecognized truncation mode");
        }
    }

    TruncationMode getTruncationMode() const override
    {
        return truncationMode;
    }

    void SetMultithreading(bool enable = true) override
    {
        enableMultithreading = enable;
    }

    bool GetMultithreading() const override
    {
        return enableMultithreading;
    }

    bool setKrausCompletenessCheck(KrausCompletenessCheck mode) override
    {
        switch (mode)
        {
        case KrausCompletenessCheck::Ignore:
        case KrausCompletenessCheck::Warn:
        case KrausCompletenessCheck::Strict:
            krausCompletenessCheck = mode;
            return true;
        default:
            throw std::invalid_argument("Unrecognized Kraus completeness check mode");
        }
    }

    KrausCompletenessCheck getKrausCompletenessCheck() const override
    {
        return krausCompletenessCheck;
    }

    void setRestoreTraceAfterTruncation(bool enable) override
    {
        restoreTraceAfterTruncation = enable;
    }

    bool getRestoreTraceAfterTruncation() const override
    {
        return restoreTraceAfterTruncation;
    }

    void setHermitizeAfterTruncation(bool enable) override
    {
        hermitizeAfterTruncation = enable;
    }

    bool getHermitizeAfterTruncation() const override
    {
        return hermitizeAfterTruncation;
    }

    void RestoreTrace() override
    {
        const std::complex<double> tr = Trace();
        if (!HasSafelyPositiveTrace(tr))
            throw std::runtime_error("Cannot restore trace of an MPO operator whose trace is not safely positive");

        ScaleSite(canonicalFormValid ? 0 : centerFirst, 1. / tr);
    }

    std::complex<double> Trace() const override
    {
        return ContractChain([this](IndexType q) { return SiteTraceMatrix(q); });
    }

    std::complex<double> TraceOfSquare() const override
    {
        // the environment matrix products are parallelized by Eigen
        std::complex<double> result;
        RunMaybeSingleThreaded(enableMultithreading, [&]() { result = TraceOfSquareImpl(); });

        return result;
    }

    std::complex<double> TraceOfSquareImpl() const
    {
        const size_t n = gammas.size();
        if (n == 0)
            return 0.;

        MatrixClass env = MatrixClass::Ones(1, 1);
        int64_t exponent = 2 * scaleExponent;

        for (size_t q = 0; q < n; ++q)
        {
            const auto &g = gammas[q];
            const IndexType R = g.dimension(3);
            MatrixClass next = MatrixClass::Zero(R, R);

            for (IndexType ket = 0; ket < 2; ++ket)
                for (IndexType bra = 0; bra < 2; ++bra)
                {
                    const auto Gkb = MapSiteSlice(g, ket, bra);
                    const auto Gbk = MapSiteSlice(g, bra, ket);
                    next.noalias() += Gkb.transpose() * env * Gbk;
                }

            // the lambdas are already included in the B tensors
            env = std::move(next);
            RescaleIfOutOfRange(env, exponent);
        }

        return ScaleByPowerOfTwo(env(0, 0), exponent);
    }

    double Purity() const override
    {
        const std::complex<double> tr = Trace();
        if (!HasSafelyPositiveTrace(tr))
            throw std::runtime_error("Cannot compute purity of an MPO operator whose trace is not safely positive");

        return (TraceOfSquare() / (tr * tr)).real();
    }

    double HermiticityResidual() const override
    {
        if (gammas.empty())
            return 0.;
        // QR the difference MPO itself. Subtracting two squared norms would lose
        // sensitivity near Hermiticity and cannot support the default 1e-10 test.
        auto difference = gammas;
        auto bonds = lambdas;
        std::vector<TensorType> negativeAdjoint;
        negativeAdjoint.reserve(gammas.size());
        for (const auto &gamma : gammas)
            negativeAdjoint.emplace_back(AdjointSite(gamma));
        for (IndexType i = 0; i < negativeAdjoint[0].size(); ++i)
            negativeAdjoint[0].data()[i] *= -1.;
        AddState(bonds, difference, negativeAdjoint);
        double residual = 0.;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            MatrixClass transfer = MatrixClass::Ones(1, 1);
            int64_t exponent = scaleExponent;
            for (const auto &gamma : difference)
            {
                const IndexType R = gamma.dimension(3), rows = transfer.rows();
                const Eigen::Map<const MatrixClass> site(gamma.data(), gamma.dimension(0), 4 * R);
                MatrixClass expanded = transfer * site;
                // keep the QR input in range: the transfer accumulates the norm of the sweep
                RescaleIfOutOfRange(expanded, exponent);
                const Eigen::Map<const MatrixClass> matrix(expanded.data(), 4 * rows, R);
                Eigen::HouseholderQR<MatrixClass> qr(matrix);
                transfer = qr.matrixQR().topRows(std::min(4 * rows, R)).template triangularView<Eigen::Upper>();
            }
            residual = std::ldexp(transfer.stableNorm(), ClampExponent(exponent));
        });
        return residual;
    }

    bool IsHermitian(double eps = 1E-10) const override
    {
        return HermiticityResidual() < eps;
    }

    MatrixClass PartialTrace(const std::vector<IndexType> &keepQubits) const override
    {
        const size_t nrQubits = getNrQubits();
        const size_t numKeep = keepQubits.size();
        if (numKeep > nrQubits)
            throw std::invalid_argument("Keep qubits set size exceeds total qubits");

        std::vector<bool> isKept(nrQubits, false);
        for (IndexType q : keepQubits)
        {
            if (q < 0 || static_cast<size_t>(q) >= nrQubits)
                throw std::invalid_argument("Qubit index out of bounds");
            if (isKept[static_cast<size_t>(q)])
                throw std::invalid_argument("Duplicate qubit index in keepQubits");
            isKept[static_cast<size_t>(q)] = true;
        }

        const size_t dimA = CheckedDensityMatrixDimension(numKeep);
        const std::complex<double> tr = Trace();
        RequireNormalizableTrace(tr);
        // Absorb traced segments once, including gaps between retained sites.
        // Reuse prefix vectors along the output tree: temporary storage stays
        // polynomial in bond dimension, even for the largest dense result.
        MatrixClass gap = MatrixClass::Ones(1, 1);
        // each rescale of the gap scales exactly one factor of every output element
        int64_t exponent = scaleExponent;
        std::vector<std::array<MatrixClass, 4>> reduced;
        std::vector<size_t> orderedPositions;
        for (size_t q = 0; q < nrQubits; ++q)
        {
            if (!isKept[q])
            {
                gap = (gap * SiteTraceMatrix(q)).eval();
                RescaleIfOutOfRange(gap, exponent);
                continue;
            }
            reduced.emplace_back();
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType ket = 0; ket < 2; ++ket)
                    reduced.back()[ket + 2 * bra].noalias() = gap * SiteSelectMatrix(q, ket, bra);
            const IndexType R = gammas[q].dimension(3);
            gap = MatrixClass::Identity(R, R);
            orderedPositions.push_back(
                static_cast<size_t>(std::find(keepQubits.begin(), keepQubits.end(), q) - keepQubits.begin()));
        }
        if (reduced.empty())
            return MatrixClass::Ones(1, 1);
        for (auto &last : reduced.back())
            last = (last * gap).eval();
        std::vector<MatrixClass> prefix(reduced.size() + 1);
        prefix[0] = MatrixClass::Ones(1, 1);
        MatrixClass rhoA(dimA, dimA);
        auto expand = [&](auto &&self, size_t depth, size_t row, size_t col) -> void {
            if (depth == reduced.size())
            {
                rhoA(row, col) = ScaleByPowerOfTwo(prefix[depth](0, 0), exponent) / tr;
                return;
            }
            for (size_t physical = 0; physical < 4; ++physical)
            {
                prefix[depth + 1].noalias() = prefix[depth] * reduced[depth][physical];
                self(self, depth + 1, row | ((physical & 1) << orderedPositions[depth]),
                     col | ((physical >> 1) << orderedPositions[depth]));
            }
        };
        expand(expand, 0, 0, 0);
        return rhoA;
    }

    std::complex<double> HilbertSchmidtOverlap(const MPOSimulatorInterface &other) const override
    {
        const size_t n = getNrQubits();
        if (other.getNrQubits() != n)
            throw std::invalid_argument("MPO register sizes do not match");
        const auto tr1 = Trace(), tr2 = other.Trace();
        RequireNormalizableTrace(tr1);
        RequireNormalizableTrace(tr2);

        const auto *otherBase = dynamic_cast<const MPOSimulatorBase *>(&other);
        if (!otherBase)
        {
            // A decorator knows its logical mapping and can align a copy with this
            // physical chain. Conjugation reverses the inner-product arguments.
            return std::conj(other.OverlapWithPhysicalChain(*this));
        }

        if (n == 0)
            return 0.;

        MatrixClass env = MatrixClass::Ones(1, 1);
        int64_t exponent = scaleExponent + otherBase->scaleExponent;

        // the environment matrix products are parallelized by Eigen
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            for (size_t q = 0; q < n; ++q)
            {
                const auto &g1 = gammas[q];
                const auto &g2 = otherBase->gammas[q];

                const IndexType R1 = g1.dimension(3);
                const IndexType R2 = g2.dimension(3);

                MatrixClass next = MatrixClass::Zero(R1, R2);

                for (IndexType ket = 0; ket < 2; ++ket)
                    for (IndexType bra = 0; bra < 2; ++bra)
                    {
                        const auto G1 = MapSiteSlice(g1, ket, bra);
                        const auto G2 = MapSiteSlice(g2, ket, bra);
                        next.noalias() += G1.adjoint() * env * G2;
                    }

                // the lambdas are already included in the B tensors
                env = std::move(next);
                RescaleIfOutOfRange(env, exponent);
            }
        });

        return ScaleByPowerOfTwo(env(0, 0), exponent) / (std::conj(tr1) * tr2);
    }

    double FidelityWithStatevector(const VectorClass &psi) const override
    {
        const size_t nrQubits = getNrQubits();
        const size_t dim = CheckedStatevectorDimension(nrQubits);
        if (psi.size() < 0 || static_cast<size_t>(psi.size()) != dim)
            throw std::invalid_argument("Statevector dimension does not match the register");

        const std::complex<double> tr = Trace();
        RequireNormalizableTrace(tr);
        // Preserve the cheap path for basis states and other sparse vectors.
        // Dense vectors stop this scan early and use the contraction below.
        const auto nonzero = CollectSparseStatevector(psi);
        if (nonzero.size() <= sparseStatevectorLimit)
            return FidelityWithSparseStatevector(nonzero, tr);
        // Apply the MPO to the dense vector. At step q, low bits are output
        // ket bits and high bits are uncontracted input bra bits: 2^N entries
        // per open bond rather than 4^N separately contracted matrix elements.
        MatrixClass current = psi;
        for (size_t q = 0; q < nrQubits; ++q)
        {
            const auto &gamma = gammas[q];
            const IndexType L = gamma.dimension(0), R = gamma.dimension(3);
            const size_t bit = size_t{1} << q;
            if (dim > static_cast<size_t>(std::numeric_limits<IndexType>::max()) / sizeof(std::complex<double>) /
                          static_cast<size_t>(R))
                throw std::length_error("MPO-statevector contraction workspace is too large");
            MatrixClass next = MatrixClass::Zero(static_cast<IndexType>(dim), R);
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < L; ++l)
                    for (IndexType bra = 0; bra < 2; ++bra)
                        for (IndexType ket = 0; ket < 2; ++ket)
                        {
                            const auto value = gamma(l, ket, bra, r);
                            if (value == std::complex<double>{})
                                continue;
                            for (size_t high = 0; high < dim; high += 2 * bit)
                                next.col(r).segment(high + ket * bit, bit) +=
                                    value * current.col(l).segment(high + bra * bit, bit);
                        }
            current = std::move(next);
        }
        return std::clamp((ScaleByPowerOfTwo(psi.dot(current.col(0)), scaleExponent) / tr).real(), 0., 1.);
    }

    std::complex<double> UnnormalizedExpectationValue(const std::string &pauliString) const override
    {
        const size_t nrQubits = getNrQubits();
        if (pauliString.size() != nrQubits)
            throw std::invalid_argument("Pauli string length must match the number of qubits");

        std::vector<MatrixClass> siteOps(nrQubits);
        for (size_t i = 0; i < nrQubits; ++i)
            siteOps[i] = PauliMatrixFromChar(pauliString[i]);

        return ContractChain(
            [this, &siteOps](IndexType q) { return SitePauliMatrix(q, siteOps[static_cast<size_t>(q)]); });
    }

    std::complex<double> ExpectationValue(const std::string &pauliString) const override
    {
        const std::complex<double> num = UnnormalizedExpectationValue(pauliString);
        const std::complex<double> tr = Trace();
        RequireNormalizableTrace(tr);

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
        RequireNormalizableTrace(tr);

        return ClampProbability((num / tr).real());
    }

    // samples a full computational basis outcome from the density matrix populations without
    // collapsing the state, by sampling qubit after qubit conditioned on the previous outcomes
    std::unordered_map<IndexType, bool> MeasureNoCollapse() override
    {
        const IndexType n = static_cast<IndexType>(gammas.size());
        if (n == 0)
            return {};

        return MeasureNoCollapseUpTo(n - 1);
    }

    // samples a subset of qubits without collapsing the state. As with the MPS simulator, the
    // sampling proceeds along the chain up to the largest requested qubit; only the requested
    // qubits are reported. The map keys are the qubit indices, the values the outcomes.
    std::unordered_map<IndexType, bool> MeasureNoCollapse(const std::set<IndexType> &qubits) override
    {
        if (qubits.empty())
            return {};
        for (const IndexType qubit : qubits)
            if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
                throw std::invalid_argument("Qubit index out of bounds");

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

    // the impl works directly on the physical chain, so there's nothing to move here;
    // the qubit reordering is handled by the MPOSimulator decorator
    void MoveAtBeginningOfChain(const std::set<IndexType> &qubits) override
    {
        for (const IndexType qubit : qubits)
            if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
                throw std::invalid_argument("Qubit index out of bounds");

        // do nothing, it's here just to provide an implementation
    }

    std::complex<double> getBasisStateMatrixElement(size_t row, size_t col) const override
    {
        const size_t nrQubits = getNrQubits();

        std::vector<bool> rowState(nrQubits, false);
        std::vector<bool> colState(nrQubits, false);

        for (size_t i = 0; i < nrQubits; ++i)
        {
            rowState[i] = (row & 1) == 1;
            row >>= 1;
        }

        for (size_t i = 0; i < nrQubits; ++i)
        {
            colState[i] = (col & 1) == 1;
            col >>= 1;
        }

        return getBasisStateMatrixElement(rowState, colState);
    }

    std::complex<double> getBasisStateMatrixElement(const std::vector<bool> &row,
                                                    const std::vector<bool> &col) const override
    {
        const size_t nrQubits = getNrQubits();
        if (nrQubits == 0)
            return 0.;

        return ContractChain([this, &row, &col](IndexType q) {
            const int ket = (static_cast<size_t>(q) < row.size() && row[q]) ? 1 : 0;
            const int bra = (static_cast<size_t>(q) < col.size() && col[q]) ? 1 : 0;
            return SiteSelectMatrix(q, ket, bra);
        });
    }

    double getBasisStateProbability(size_t State) const override
    {
        const std::complex<double> tr = Trace();
        RequireNormalizableTrace(tr);

        return ClampProbability((getBasisStateMatrixElement(State, State) / tr).real());
    }

    double getBasisStateProbability(const std::vector<bool> &State) const override
    {
        const std::complex<double> tr = Trace();
        RequireNormalizableTrace(tr);

        return ClampProbability((getBasisStateMatrixElement(State, State) / tr).real());
    }

    // this is costly (it builds the full 2^N x 2^N matrix) and it's meant only
    // for comparing the results against other simulators, not for simulation.
    //
    // getDensityMatrix returns rho / Tr(rho), matching getBasisStateProbability() and
    // ExpectationValue(). That scale matches a state; it cannot restore positivity.
    // It throws if Re(Tr(rho)) is not safely positive — dividing by a vanished or
    // negative trace would produce a sign-flipped or non-finite "state".
    // getUnnormalizedDensityMatrix returns the raw MPO operator.
    MatrixClass getUnnormalizedDensityMatrix() const override
    {
        return ReconstructOperatorMatrix();
    }

    MatrixClass getDensityMatrix() const override
    {
        MatrixClass rho = ReconstructOperatorMatrix();
        if (rho.size() == 0)
            return rho;

        const std::complex<double> tr = Trace();
        if (!HasSafelyPositiveTrace(tr))
            throw std::runtime_error("Cannot normalize an MPO operator whose trace is not safely positive");

        rho /= tr;
        return rho;
    }

    std::shared_ptr<MPOSimulatorStateInterface> getState() const override
    {
        auto state = std::make_shared<MPOSimulatorBaseState>();
        state->lambdas = lambdas;
        state->gammas = gammas;
        state->scaleExponent = scaleExponent;

        return state;
    }

    void setState(const std::shared_ptr<MPOSimulatorStateInterface> &state) override
    {
        if (!state)
            return;

        const auto stateRef = CheckedBaseState(state);
        std::vector<LambdaType> newLambdas = stateRef->lambdas;
        std::vector<TensorType> newGammas = stateRef->gammas;
        lambdas.swap(newLambdas);
        gammas.swap(newGammas);
        scaleExponent = stateRef->scaleExponent;
        InvalidateCanonicalForm();
    }

    void setStateDestructive(std::shared_ptr<MPOSimulatorStateInterface> &state) override
    {
        if (!state)
            return;

        auto stateRef = CheckedBaseState(state);
        lambdas.swap(stateRef->lambdas);
        gammas.swap(stateRef->gammas);
        std::swap(scaleExponent, stateRef->scaleExponent);
        InvalidateCanonicalForm();
    }

    void print() const override
    {
        for (size_t i = 0; i < gammas.size(); ++i)
        {
            std::cout << std::endl << "B (Gamma * Lambda) " << i << " (leftBond, ket, bra, rightBond):" << std::endl;
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
        for (const auto &lambda : lambdas)
            std::cout << lambda.size() << " ";
        std::cout << std::endl;
    }

  protected:
    using SiteSliceMap = Eigen::Map<const MatrixClass, 0, Eigen::OuterStride<>>;
    static SiteSliceMap MapSiteSlice(const TensorType &gamma, IndexType ket, IndexType bra)
    {
        const IndexType L = gamma.dimension(0);
        // Tensor storage is (left, ket, bra, right); adjacent columns are 4L apart.
        return SiteSliceMap(gamma.data() + L * (ket + 2 * bra), L, gamma.dimension(3), Eigen::OuterStride<>(4 * L));
    }

    static constexpr size_t sparseStatevectorLimit = 32;
    using SparseStatevector = std::vector<std::pair<IndexType, std::complex<double>>>;
    static SparseStatevector CollectSparseStatevector(const VectorClass &psi)
    {
        SparseStatevector nonzero;
        nonzero.reserve(sparseStatevectorLimit + 1);
        for (IndexType i = 0; i < psi.size(); ++i)
            if (psi[i] != std::complex<double>{})
            {
                nonzero.emplace_back(i, psi[i]);
                if (nonzero.size() > sparseStatevectorLimit)
                    break;
            }
        return nonzero;
    }
    double FidelityWithSparseStatevector(const SparseStatevector &psi, const std::complex<double> &trace) const
    {
        std::complex<double> overlap = 0.;
        for (const auto &row : psi)
            for (const auto &col : psi)
                overlap += std::conj(row.second) * getBasisStateMatrixElement(row.first, col.first) * col.second;
        return std::clamp((overlap / trace).real(), 0., 1.);
    }

    // Outside this interval the left prefix / right suffix are orthonormal.
    // A valid Vidal form additionally has current Schmidt weights on every bond.
    CanonicalMetadata GetCanonicalMetadata() const
    {
        return {canonicalFormValid, centerFirst, centerLast};
    }
    void RestoreCanonicalMetadata(const CanonicalMetadata &metadata)
    {
        canonicalFormValid = metadata.valid;
        centerFirst = metadata.first;
        centerLast = metadata.last;
        InvalidateSamplingCache();
    }
    void InvalidateSamplingCache()
    {
        samplingRight.clear();
    }
    void InvalidateCanonicalForm(IndexType first, IndexType last)
    {
        if (canonicalFormValid)
            centerFirst = centerLast = 0;
        canonicalFormValid = false;
        centerFirst = std::min(centerFirst, first);
        centerLast = std::max(centerLast, last);
        InvalidateSamplingCache();
    }
    void InvalidateCanonicalForm()
    {
        InvalidateCanonicalForm(0, static_cast<IndexType>(gammas.size()) - 1);
    }

    std::shared_ptr<MPOSimulatorBaseState> CheckedBaseState(
        const std::shared_ptr<MPOSimulatorStateInterface> &state) const
    {
        if (typeid(*state) != typeid(MPOSimulatorBaseState))
            throw std::invalid_argument("State type is incompatible with MPOSimulatorImpl");

        auto baseState = std::dynamic_pointer_cast<MPOSimulatorBaseState>(state);
        if (!baseState || baseState->gammas.size() != gammas.size() || baseState->lambdas.size() != lambdas.size())
            throw std::invalid_argument("MPO state dimensions do not match the simulator");

        const size_t n = baseState->gammas.size();
        for (size_t q = 0; q < n; ++q)
        {
            const auto &gamma = baseState->gammas[q];
            const IndexType leftBond = gamma.dimension(0);
            const IndexType rightBond = gamma.dimension(3);
            if (leftBond <= 0 || rightBond <= 0 || gamma.dimension(1) != 2 || gamma.dimension(2) != 2 ||
                (q == 0 && leftBond != 1) || (q + 1 == n && rightBond != 1))
                throw std::invalid_argument("MPO state contains an invalid site tensor");

            for (IndexType i = 0; i < gamma.size(); ++i)
            {
                const auto value = gamma.data()[i];
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::invalid_argument("MPO state contains a non-finite tensor value");
            }

            if (q + 1 < n)
            {
                const auto &lambda = baseState->lambdas[q];
                if (lambda.size() <= 0 || rightBond != lambda.size() ||
                    baseState->gammas[q + 1].dimension(0) != rightBond || !lambda.allFinite() ||
                    (lambda.array() < 0.).any())
                    throw std::invalid_argument("MPO state contains an invalid bond");
            }
        }

        return baseState;
    }

    static double ClampProbability(double probability)
    {
        constexpr double tolerance = 1E-12;
        if (probability < 0. && probability > -tolerance)
            return 0.;
        if (probability > 1. && probability < 1. + tolerance)
            return 1.;

        return probability;
    }

    // ---- power-of-two scale handling (see the class comment) ----
    static constexpr int scaleWindowExponent = 64;

    // The exponent e with x / 2^e in [0.5, 1) when the magnitude x is outside the window,
    // otherwise 0. Zero and non-finite magnitudes are left alone.
    static int OutOfRangeExponent(double magnitude)
    {
        if (!(magnitude > 0.) || !std::isfinite(magnitude))
            return 0;
        int e = 0;
        std::frexp(magnitude, &e);
        return (e < -scaleWindowExponent || e > scaleWindowExponent) ? e : 0;
    }

    // The largest real or imaginary component: within sqrt(2) of the largest magnitude, which
    // is all the window needs, and without a hypot per element.
    template <class Derived> static double MaxAbs(const Eigen::MatrixBase<Derived> &m)
    {
        if (m.size() == 0)
            return 0.;
        return std::max(m.real().cwiseAbs().maxCoeff(), m.imag().cwiseAbs().maxCoeff());
    }

    // Divides m by 2^e, with e from its largest component when that is out of the window,
    // and adds e to exponent. Returns e.
    template <class Derived> static int RescaleIfOutOfRange(Eigen::MatrixBase<Derived> &m, int64_t &exponent)
    {
        const int e = OutOfRangeExponent(MaxAbs(m));
        if (e != 0)
        {
            m *= std::ldexp(1., -e);
            exponent += e;
        }
        return e;
    }

    static int ClampExponent(int64_t exponent)
    {
        constexpr int64_t limit = int64_t{1} << 20; // far beyond the double range; ldexp saturates
        return static_cast<int>(std::clamp<int64_t>(exponent, -limit, limit));
    }

    static std::complex<double> ScaleByPowerOfTwo(const std::complex<double> &value, int64_t exponent)
    {
        const int e = ClampExponent(exponent);
        return {std::ldexp(value.real(), e), std::ldexp(value.imag(), e)};
    }

    static bool HasSafelyPositiveTrace(const std::complex<double> &trace)
    {
        return std::isfinite(trace.real()) && std::isfinite(trace.imag()) &&
               trace.real() > std::numeric_limits<double>::epsilon();
    }
    static void RequireNormalizableTrace(const std::complex<double> &trace)
    {
        if (!HasSafelyPositiveTrace(trace))
            throw std::runtime_error("Cannot normalize an MPO operator whose trace is not safely positive");
    }

    MatrixClass ReconstructOperatorMatrix() const
    {
        const size_t sz = gammas.size();
        if (sz == 0)
            return {};
        const size_t NrBasisStates = CheckedDensityMatrixDimension(sz);
        MatrixClass rho(NrBasisStates, NrBasisStates);

        for (size_t r = 0; r < NrBasisStates; ++r)
            for (size_t c = 0; c < NrBasisStates; ++c)
                rho(r, c) = getBasisStateMatrixElement(r, c);

        return rho;
    }

    static double ValidMeasurementProbability(double probability)
    {
        constexpr double toleranceLow = 1E-9;
        constexpr double toleranceHigh = 1E-5;
        if (probability < -toleranceLow || probability > 1. + toleranceHigh)
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
        if (n == 0)
            return res;

        if (limit < 0)
            limit = 0;
        if (limit >= n)
            limit = n - 1;

        // Reuse traced suffixes until the tensor representation changes. Rescale locally:
        // their common scale cancels in each conditional probability, avoiding small
        // joint probabilities and repeated whole-chain contractions.
        using SiteMap = Eigen::Map<const MatrixClass, 0, Eigen::OuterStride<>>;
        if (samplingRight.empty())
        {
            std::vector<MatrixClass> right(static_cast<size_t>(n) + 1);
            right[n] = MatrixClass::Ones(1, 1);
            for (IndexType q = n - 1; q >= 0; --q)
            {
                const auto &g = gammas[q];
                const IndexType L = g.dimension(0), R = g.dimension(3);
                const SiteMap zero(g.data(), L, R, Eigen::OuterStride<>(4 * L));
                const SiteMap one(g.data() + 3 * L, L, R, Eigen::OuterStride<>(4 * L));
                right[q].noalias() = zero * right[q + 1];
                right[q].noalias() += one * right[q + 1];
                RescaleSamplingEnvironment(right[q]);
            }
            samplingRight.swap(right);
        }

        MatrixClass left = MatrixClass::Ones(1, 1);
        MatrixClass left0, left1;
        res.reserve(static_cast<size_t>(limit) + 1);
        for (IndexType qubit = 0; qubit <= limit; ++qubit)
        {
            const auto &g = gammas[qubit];
            const IndexType L = g.dimension(0), R = g.dimension(3);
            left0.noalias() = left * SiteMap(g.data(), L, R, Eigen::OuterStride<>(4 * L));
            left1.noalias() = left * SiteMap(g.data() + 3 * L, L, R, Eigen::OuterStride<>(4 * L));
            const std::complex<double> weight0 = (left0 * samplingRight[qubit + 1])(0, 0);
            const std::complex<double> weight1 = (left1 * samplingRight[qubit + 1])(0, 0);
            const std::complex<double> total = weight0 + weight1;
            if (!std::isfinite(total.real()) || !std::isfinite(total.imag()) || std::abs(total) == 0.)
                throw std::runtime_error("Cannot sample an MPO with zero or non-finite conditional weight");
            const double probability = (weight0 / total).real();
            if (!std::isfinite(probability))
                throw std::runtime_error("Cannot sample an MPO with a non-finite probability");
            const double prob0 = ValidMeasurementProbability(probability);
            const bool zeroMeasured = uniformZeroOne(rng) < prob0;
            res[qubit] = !zeroMeasured;
            left.swap(zeroMeasured ? left0 : left1);
            RescaleSamplingEnvironment(left);
        }

        return res;
    }

    static void RescaleSamplingEnvironment(MatrixClass &environment)
    {
        const double scale = environment.cwiseAbs().maxCoeff();
        if (!environment.allFinite() || !std::isfinite(scale) || scale == 0.)
            throw std::runtime_error("Cannot sample an MPO with zero or non-finite weight");
        environment /= scale;
    }

    static void SetSiteToBasis(TensorType &gamma, int basis)
    {
        gamma.setZero();
        gamma(0, basis, basis, 0) = 1.; // |basis><basis|
    }

    void PrintGamma(size_t i) const
    {
        assert(i < gammas.size());

        const auto &g = gammas[i];
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
        const auto &g = gammas[q];
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
        const auto &g = gammas[q];
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
            m(0, 0) = 1.;
            m(1, 1) = 1.;
            break;
        case 'X':
            m(0, 1) = 1.;
            m(1, 0) = 1.;
            break;
        case 'Y':
            m(0, 1) = std::complex<double>(0., -1.);
            m(1, 0) = std::complex<double>(0., 1.);
            break;
        case 'Z':
            m(0, 0) = 1.;
            m(1, 1) = -1.;
            break;
        default:
            throw std::invalid_argument("Invalid operator in the Pauli string");
        }

        return m;
    }

    // matrix (leftBond x rightBond) obtained by contracting a single qubit operator P into
    // the physical legs of a site: m(l, r) = sum_{ket,bra} g(l, ket, bra, r) P(bra, ket).
    // With P = Identity this reduces to SiteTraceMatrix.
    MatrixClass SitePauliMatrix(IndexType q, const MatrixClass &P) const
    {
        const auto &g = gammas[q];
        const IndexType L = g.dimension(0);
        const IndexType R = g.dimension(3);

        MatrixClass m = MatrixClass::Zero(L, R);
        for (IndexType ket = 0; ket < 2; ++ket)
            for (IndexType bra = 0; bra < 2; ++bra)
            {
                const std::complex<double> p = P(bra, ket);
                if (p == std::complex<double>(0., 0.))
                    continue;

                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                        m(l, r) += g(l, ket, bra, r) * p;
            }

        return m;
    }

    // contracts the whole chain, picking at each site a (leftBond x rightBond) matrix
    // supplied by 'siteMatrix' (the lambdas are already included in the B tensors)
    template <typename SiteMatrixFunc> std::complex<double> ContractChain(SiteMatrixFunc siteMatrix) const
    {
        const size_t n = gammas.size();
        if (n == 0)
            return 0.;

        MatrixClass res = siteMatrix(0);
        int64_t exponent = scaleExponent;
        RescaleIfOutOfRange(res, exponent);

        for (size_t q = 1; q < n; ++q)
        {
            const MatrixClass sm = siteMatrix(static_cast<IndexType>(q));
            res = (res * sm).eval();
            RescaleIfOutOfRange(res, exponent);
        }

        return ScaleByPowerOfTwo(res(0, 0), exponent);
    }

    void RestoreTraceIfSafe()
    {
        const std::complex<double> tr = Trace();
        if (!HasSafelyPositiveTrace(tr))
            return;

        ScaleSite(canonicalFormValid ? 0 : centerFirst, 1. / tr);
    }

    static TensorType AdjointSite(const TensorType &gamma)
    {
        TensorType adjoint(gamma.dimension(0), 2, 2, gamma.dimension(3));
        for (IndexType r = 0; r < gamma.dimension(3); ++r)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType ket = 0; ket < 2; ++ket)
                    for (IndexType l = 0; l < gamma.dimension(0); ++l)
                        adjoint(l, ket, bra, r) = std::conj(gamma(l, bra, ket, r));

        return adjoint;
    }

    void ScaleSite(IndexType q, std::complex<double> factor)
    {
        InvalidateSamplingCache();
        if (canonicalFormValid && q == 0 && std::isfinite(std::abs(factor)) && std::abs(factor) > 0.)
        {
            for (auto &lambda : lambdas)
                lambda *= std::abs(factor);
        }
        else
            InvalidateCanonicalForm(q, q);
        auto &g = gammas[q];
        for (IndexType r = 0; r < g.dimension(3); ++r)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType ket = 0; ket < 2; ++ket)
                    for (IndexType l = 0; l < g.dimension(0); ++l)
                        g(l, ket, bra, r) *= factor;
    }

    // lhs = lhs + rhs, the operators are the products of the site tensors (the lambdas are already included in them)
    // the resulting bonds are the direct sums of the bonds, the lambdas are set to ones, they are only
    // placeholders for the left environment until the next ReCanonicalize
    static void AddState(std::vector<LambdaType> &lhsLambdas, std::vector<TensorType> &lhsGammas,
                         const std::vector<TensorType> &rhsGammas)
    {
        assert(lhsGammas.size() == rhsGammas.size());

        const size_t n = lhsGammas.size();
        if (n == 1)
        {
            lhsGammas[0] = (lhsGammas[0] + rhsGammas[0]).eval();
            return;
        }

        std::vector<TensorType> sumGammas(n);
        std::vector<LambdaType> sumLambdas(n - 1);

        for (size_t q = 0; q < n; ++q)
        {
            const auto &lg = lhsGammas[q];
            const auto &rg = rhsGammas[q];

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

    static void CopyGammaBlock(const TensorType &source, TensorType &target, IndexType leftOffset,
                               IndexType rightOffset)
    {
        for (IndexType r = 0; r < source.dimension(3); ++r)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType ket = 0; ket < 2; ++ket)
                    for (IndexType l = 0; l < source.dimension(0); ++l)
                        target(leftOffset + l, ket, bra, rightOffset + r) = source(l, ket, bra, r);
    }

    // applies a single qubit local operator to a site tensor: rho_site -> A rho_site A^dagger
    // the operator A is contracted with the ket leg, the conjugate operator A* with the bra leg
    static void ApplySingleQubitGate(TensorType &gamma, const GateClass &gate)
    {
        ApplySingleQubitGate(gamma, gate.getRawOperatorMatrix());
    }

    static void ApplySingleQubitGate(TensorType &gamma, const MatrixClass &opMat)
    {
        static const Indexes contractKet{IntIndexPair(1, 1)}; // gamma ket (dim 1) with U column (dim 1)
        static const Indexes contractBra{IntIndexPair(1, 1)}; // intermediate bra (dim 1) with U* column (dim 1)
        static const std::array<int, 4> permute{0, 2, 3, 1};

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
    IndexType chi = 10;                 // if limitSize is true
    double singularValueThreshold = 0.; // if limitEntanglement is true

    // Default is DiscardedWeight, not RelativeToMax - see MPSSimulatorBase's identical
    // field for the full rationale. Applies equally to the MPO simulator: the selection
    // criterion normalizes against this bond's own pre-truncation Schmidt weight
    // regardless of the fact that MPO does not renormalize the *kept* lambdas afterward
    // (that is a separate, unrelated choice made to keep Tr(rho) stable - see the class
    // comment above and DecomposeAndSetGammas in MPOSimulatorImpl.h).
    TruncationMode truncationMode = TruncationMode::DiscardedWeight;
    KrausCompletenessCheck krausCompletenessCheck = KrausCompletenessCheck::Ignore;
    // if false, the SVDs and the matrix products are done single threaded, see SetMultithreading
    bool enableMultithreading = true;
    bool restoreTraceAfterTruncation = false;
    bool hermitizeAfterTruncation = false;
    // Public snapshots are untrusted; private SaveState snapshots retain this metadata.
    bool canonicalFormValid = true;
    IndexType centerFirst = 0, centerLast = 0;
    std::vector<MatrixClass> samplingRight;

    std::vector<LambdaType> lambdas;
    std::vector<TensorType> gammas;
    // rho = 2^scaleExponent * B[0] ... B[N-1], see the class comment
    int64_t scaleExponent = 0;

    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniformZeroOne{0, 1};

    const Operators::ZeroProjection<MatrixClass> zeroProjection;
    const Operators::OneProjection<MatrixClass> oneProjection;
};

} // namespace TensorNetworks

} // namespace QC
