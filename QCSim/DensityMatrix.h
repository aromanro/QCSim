#pragma once

#include "QubitRegister.h"
#include "QubitRegisterCalculator.h"
#include "SimpleGates.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// A density-matrix quantum computing simulator.
//
// The state is an owning 2^n x 2^n matrix. A zero-copy, contiguous vector view lets the
// statevector kernels act on all its entries at once: U on the ket bits and conjugate(U)
// on the bra bits. For column-major storage the ket bits are low and the bra bits high;
// row-major storage reverses those roles. No full-register operator is constructed.
//
// Qubits are numbered from right to left, starting with zero (same convention as QubitRegister).

namespace QC
{

inline constexpr size_t DefaultDensityParallelMinElements = 65536;

inline std::atomic<size_t> &DensityParallelMinElementsSetting()
{
    static std::atomic<size_t> value{DefaultDensityParallelMinElements};
    return value;
}

template <class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DensityMatrix
{
  public:
    using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
    using Scalar = typename MatrixClass::Scalar;
    using FlatVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using FlatView = Eigen::Map<FlatVector>;
    using FlatCalculator = QubitRegisterCalculator<FlatView, MatrixClass>;

    // An immutable probability snapshot. It stays valid after the simulator changes;
    // callers explicitly prepare a new snapshot when they want the new distribution.
    class PreparedSampler
    {
        friend class DensityMatrix;

      public:
        // Empty snapshots can be stored before preparation. Moving to a distinct snapshot
        // leaves its source empty, independently of vector move semantics; self-move is a no-op.
        PreparedSampler() noexcept = default;
        PreparedSampler(const PreparedSampler &) = default;

        PreparedSampler(PreparedSampler &&other) noexcept
        {
            Swap(other);
        }

        PreparedSampler &operator=(const PreparedSampler &other)
        {
            PreparedSampler copy(other);
            Swap(copy);
            return *this;
        }

        PreparedSampler &operator=(PreparedSampler &&other) noexcept
        {
            if (this != &other)
            {
                PreparedSampler moved(std::move(other));
                Swap(moved);
            }
            return *this;
        }

        bool isValid() const noexcept
        {
            return sourceBasisStates != 0;
        }

        bool isFullRegister() const noexcept
        {
            return isValid() && getNrOutcomes() == sourceBasisStates;
        }

        size_t getNrOutcomes() const noexcept
        {
            return cumulativeProbabilities.size();
        }

        // Returned outcome bit k corresponds to register qubit getFirstQubit() + k.
        size_t getFirstQubit() const noexcept
        {
            return firstQubit;
        }

        size_t getNrQubits() const noexcept
        {
            return nrQubits;
        }

      private:
        PreparedSampler(size_t basisStates, size_t first, size_t width, std::vector<double> &&cumulative, double mass)
            : sourceBasisStates(basisStates), firstQubit(first), nrQubits(width),
              cumulativeProbabilities(std::move(cumulative)), total(mass)
        {
        }

        void Swap(PreparedSampler &other) noexcept
        {
            std::swap(sourceBasisStates, other.sourceBasisStates);
            std::swap(firstQubit, other.firstQubit);
            std::swap(nrQubits, other.nrQubits);
            cumulativeProbabilities.swap(other.cumulativeProbabilities);
            std::swap(total, other.total);
        }

        size_t sourceBasisStates = 0;
        size_t firstQubit = 0;
        size_t nrQubits = 0;
        std::vector<double> cumulativeProbabilities;
        double total = 0.;
    };

    // Legacy public aliases retained for source compatibility; execution now uses FlatView.
    using ColXpr = decltype(std::declval<MatrixClass &>().col(0));
    using RowXpr = decltype(std::declval<MatrixClass &>().row(0));

    using ColCalculator = QubitRegisterCalculator<ColXpr, MatrixClass>;
    using RowCalculator = QubitRegisterCalculator<RowXpr, MatrixClass>;

    DensityMatrix(size_t N = 3, unsigned int addseed = 0)
        : NrQubits(N), NrBasisStates(CheckedBasisStateCount(N)), rho(MatrixClass::Zero(NrBasisStates, NrBasisStates)),
          uniformZeroOne(0, 1)
    {
        if (addseed == 0)
        {
            std::random_device rdl;
            addseed = rdl();
        }

        const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
        std::seed_seq seed{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(seed);

        rho(0, 0) = 1.; // |0...0><0...0|
    }

    void SetSeed(uint64_t theSeed)
    {
        std::seed_seq seed{uint32_t(theSeed & 0xffffffff), uint32_t(theSeed >> 32)};
        rng.seed(seed);
    }

    size_t getNrQubits() const
    {
        return NrQubits;
    }

    size_t getNrBasisStates() const
    {
        return NrBasisStates;
    }

    void SetMultithreading(bool enable = true)
    {
        enableMultithreading = enable;
    }

    bool GetMultithreading() const
    {
        return enableMultithreading;
    }

    // Independent of SetParallelMinBasisStates: density work is measured in total matrix elements.
    // The default is 65,536 elements (8 qubits), rather than the former 16,384-element row/column cutoff.
    static void SetParallelMinElements(size_t count)
    {
        DensityParallelMinElementsSetting().store(count, std::memory_order_relaxed);
    }

    static size_t GetParallelMinElements()
    {
        return DensityParallelMinElementsSetting().load(std::memory_order_relaxed);
    }

    // Allows simulations and statistical tests to be reproduced exactly.
    void SetRandomSeed(uint64_t seed)
    {
        rng.seed(seed);
    }

    const MatrixClass &getDensityMatrix() const
    {
        return rho;
    }

    void Clear()
    {
        rho.setZero();
    }

    void setToBasisState(size_t State)
    {
        if (State >= NrBasisStates)
            throw std::invalid_argument("Basis state is outside the register");

        rho.setZero();
        rho(State, State) = 1.;
    }

    void setToBasisState(const std::vector<bool> &State)
    {
        if (State.size() > NrQubits)
            throw std::invalid_argument("Basis state has more bits than the register");

        size_t stateIndex = 0;
        for (size_t i = 0; i < State.size(); ++i)
            if (State[i])
                stateIndex |= (1ULL << i);

        setToBasisState(stateIndex);
    }

    void Reset()
    {
        setToBasisState(0);
    }

    void SaveState()
    {
        savedStateStorage = rho;
    }

    void RestoreState()
    {
        if (savedStateStorage.size() == 0)
            return;

        rho = savedStateStorage;
    }

    void RestoreStateDestructive()
    {
        if (savedStateStorage.size() == 0)
            return;

        rho.swap(savedStateStorage);
        savedStateStorage.resize(0, 0);
    }

    std::unique_ptr<DensityMatrix<VectorClass, MatrixClass>> Clone() const
    {
        return std::make_unique<DensityMatrix<VectorClass, MatrixClass>>(*this);
    }

    // Initialize rho = |psi><psi| from a statevector, normalizing any finite non-zero input.
    // This is very convenient for comparing against the statevector simulator.
    void setFromStatevector(const VectorClass &psi)
    {
        if (psi.size() < 0 || static_cast<size_t>(psi.size()) != NrBasisStates)
            throw std::invalid_argument("Statevector dimension does not match the register");
        if (!psi.allFinite())
            throw std::invalid_argument("Statevector contains a non-finite value");

        const double normSquared = psi.squaredNorm();
        if (!std::isfinite(normSquared) || normSquared <= 1E-20)
            throw std::invalid_argument("Statevector must have a finite non-zero norm");

        rho = (psi * psi.adjoint()) / normSquared;
    }

    // set the whole density matrix directly (the caller is responsible for it being a valid state:
    // Hermitian, positive semidefinite and unit trace)
    void setDensityMatrix(const MatrixClass &newRho)
    {
        if (newRho.rows() < 0 || newRho.cols() < 0 || static_cast<size_t>(newRho.rows()) != NrBasisStates ||
            static_cast<size_t>(newRho.cols()) != NrBasisStates)
            throw std::invalid_argument("Density matrix dimensions do not match the register");

        rho = newRho;
    }

    // set the state to a classical mixture of computational basis states: rho = sum_k p_k |s_k><s_k|
    // the weights are normalized to sum to 1 (negative or zero weights are ignored)
    void setToMixtureOfBasisStates(const std::vector<std::pair<size_t, double>> &mixture)
    {
        double total = 0.;
        for (const auto &[state, weight] : mixture)
        {
            if (!std::isfinite(weight))
                throw std::invalid_argument("Mixture weights must be finite");
            if (weight > 0. && state < NrBasisStates)
                total += weight;
        }

        if (!std::isfinite(total) || total <= 0.)
            throw std::invalid_argument("Mixture must contain at least one valid positive weight");

        rho.setZero();

        for (const auto &[state, weight] : mixture)
            if (weight > 0. && state < NrBasisStates)
                rho(static_cast<Eigen::Index>(state), static_cast<Eigen::Index>(state)) += weight / total;
    }

    void setToMixtureOfBasisStates(std::initializer_list<std::pair<size_t, double>> mixture)
    {
        setToMixtureOfBasisStates(std::vector<std::pair<size_t, double>>(mixture));
    }

    void setToMixtureOfBasisStates(const std::vector<std::pair<std::vector<bool>, double>> &mixture)
    {
        std::vector<std::pair<size_t, double>> integerMixture;
        integerMixture.reserve(mixture.size());

        for (const auto &[bits, weight] : mixture)
        {
            if (!std::isfinite(weight))
                throw std::invalid_argument("Mixture weights must be finite");
            if (weight <= 0. || bits.size() > NrQubits)
                continue;

            size_t stateIndex = 0;
            for (size_t i = 0; i < bits.size(); ++i)
                if (bits[i])
                    stateIndex |= (1ULL << i);

            integerMixture.emplace_back(stateIndex, weight);
        }

        setToMixtureOfBasisStates(integerMixture);
    }

    // rho' = U rho U^dagger
    void ApplyGate(const GateClass &gate, size_t qubit, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0)
    {
        const size_t gateQubits = ValidateGateAndQubits(gate, qubit, controllingQubit1, controllingQubit2);
        const MatrixClass &U = gate.getRawOperatorMatrix();
        const auto structure = gate.getStructure();
        if (structure.kind == Gates::GateStructure::Kind::Diagonal)
        {
            const std::array<size_t, 3> qubits{qubit, controllingQubit1, controllingQubit2};
            if (gateQubits == 1)
                ApplyDensityDiagonal<1>(U, qubits);
            else if (gateQubits == 2)
                ApplyDensityDiagonal<2>(U, qubits);
            else
                ApplyDensityDiagonal<3>(U, qubits);
            return;
        }
        const MatrixClass Uconj = U.conjugate(); // only the small operator
        ApplyGateToBuffer(U, structure, gateQubits, qubit, controllingQubit1, controllingQubit2, true);
        ApplyGateToBuffer(Uconj, structure, gateQubits, qubit, controllingQubit1, controllingQubit2, false);
    }

    void ApplyGate(const Gates::AppliedGate<MatrixClass> &gate)
    {
        ApplyGate(gate, gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
    }

    void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>> &gates)
    {
        for (const auto &gate : gates)
            ApplyGate(gate);
    }

    // Generic completely positive trace preserving channel: rho' = sum_k E_k rho E_k^dagger.
    // The completeness relation sum_k E_k^dagger E_k = I is validated before changing the state.
    // The Kraus operators are small matrices (2x2 for a single qubit, 4x4 for two qubits) acting
    // on the given qubit(s). This is the fundamental non-unitary operation; a unitary gate is just
    // the special case of a single Kraus operator E_0 = U.
    void ApplyChannel(const std::vector<MatrixClass> &kraus, size_t qubit, size_t controllingQubit1 = 0)
    {
        if (kraus.empty())
            throw std::invalid_argument("A channel must contain at least one Kraus operator");

        const size_t gateQubits = GetOperatorQubits(kraus.front(), 2);
        ValidateQubits(gateQubits, qubit, controllingQubit1, 0);

        const size_t operatorDimension = static_cast<size_t>(kraus.front().rows());
        std::array<Scalar, 16> completeness{};
        for (const auto &E : kraus)
        {
            if (GetOperatorQubits(E, 2) != gateQubits)
                throw std::invalid_argument("All Kraus operators must have the same dimensions");

            if (!E.allFinite())
                throw std::invalid_argument("Kraus operators must contain only finite values");

            for (size_t r = 0; r < operatorDimension; ++r)
                for (size_t c = 0; c < operatorDimension; ++c)
                    for (size_t k = 0; k < operatorDimension; ++k)
                        completeness[r * operatorDimension + c] += std::conj(E(k, r)) * E(k, c);
        }

        double errorSquared = 0.;
        for (size_t r = 0; r < operatorDimension; ++r)
            for (size_t c = 0; c < operatorDimension; ++c)
                errorSquared += std::norm(completeness[r * operatorDimension + c] - Scalar(r == c ? 1. : 0.));
        const double completenessTolerance = 1E-10 * operatorDimension;
        if (!std::isfinite(errorSquared) || errorSquared > completenessTolerance * completenessTolerance)
            throw std::invalid_argument("Kraus operators do not define a trace-preserving channel");

        if (gateQubits == 1)
            ApplyLocalChannel<1>(kraus, {qubit, 0});
        else
            ApplyLocalChannel<2>(kraus, {qubit, controllingQubit1});
    }

    // ---- predefined single qubit noise channels ----

    // bit flip: rho' = (1 - p) rho + p X rho X
    void ApplyBitFlipNoise(size_t qubit, double p)
    {
        ValidateProbability(p, "Bit-flip probability");
        ValidateQubit(qubit);
        if (p == 0.)
            return;
        ForEachLocalBlock<1>({qubit, 0}, [&](Scalar *base, const auto &o) {
            if (p == 1.)
            {
                std::swap(base[o[0]], base[o[3]]);
                std::swap(base[o[1]], base[o[2]]);
                return;
            }
            const auto a = base[o[0]], b = base[o[1]], c = base[o[2]], d = base[o[3]];
            base[o[0]] = (1. - p) * a + p * d;
            base[o[3]] = (1. - p) * d + p * a;
            base[o[1]] = (1. - p) * b + p * c;
            base[o[2]] = (1. - p) * c + p * b;
        });
    }

    // phase flip: rho' = (1 - p) rho + p Z rho Z
    void ApplyPhaseFlipNoise(size_t qubit, double p)
    {
        ValidateProbability(p, "Phase-flip probability");
        ValidateQubit(qubit);
        ScaleCoherences(qubit, 1. - 2. * p);
    }

    // depolarizing: rho' = (1 - p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
    void ApplyDepolarizingNoise(size_t qubit, double p)
    {
        ValidateProbability(p, "Depolarizing probability");
        ValidateQubit(qubit);
        if (p == 0.)
            return;
        const double transfer = 2. * p / 3.;
        const double coherence = 1. - 4. * p / 3.;
        ForEachLocalBlock<1>({qubit, 0}, [&](Scalar *base, const auto &o) {
            const auto a = base[o[0]], d = base[o[3]];
            base[o[0]] = (1. - transfer) * a + transfer * d;
            base[o[3]] = (1. - transfer) * d + transfer * a;
            base[o[1]] *= coherence;
            base[o[2]] *= coherence;
        });
    }

    // amplitude damping (|1> -> |0> relaxation with probability gamma)
    void ApplyAmplitudeDamping(size_t qubit, double gamma)
    {
        ValidateProbability(gamma, "Amplitude-damping probability");
        ValidateQubit(qubit);
        if (gamma == 0.)
            return;
        if (gamma == 1.)
        {
            ApplyReset(qubit);
            return;
        }
        const double coherence = std::sqrt(1. - gamma);
        ForEachLocalBlock<1>({qubit, 0}, [&](Scalar *base, const auto &o) {
            base[o[0]] += gamma * base[o[3]];
            base[o[3]] *= 1. - gamma;
            base[o[1]] *= coherence;
            base[o[2]] *= coherence;
        });
    }

    // phase damping / dephasing, suppresses the off diagonal coherences by lambda = sqrt(1 - gamma)
    void ApplyPhaseDamping(size_t qubit, double gamma)
    {
        ValidateProbability(gamma, "Phase-damping probability");
        ValidateQubit(qubit);
        ScaleCoherences(qubit, std::sqrt(1. - gamma));
    }

    // reset a qubit to |0>: E0 = |0><0|, E1 = |0><1|
    void ApplyReset(size_t qubit)
    {
        ValidateQubit(qubit);
        ForEachLocalBlock<1>({qubit, 0}, [](Scalar *base, const auto &o) {
            base[o[0]] += base[o[3]];
            base[o[1]] = base[o[2]] = base[o[3]] = 0.;
        });
    }

    // ---- measurement ----

    double GetQubitProbability(size_t qubit) const
    {
        ValidateQubit(qubit);
        const size_t mask = 1ULL << qubit;

        double p1 = 0;
        for (size_t i = 0; i < NrBasisStates; ++i)
            if (i & mask)
                p1 += rho(i, i).real();

        return p1;
    }

    // sample a computational basis outcome for a single qubit and collapse
    size_t MeasureQubit(size_t qubit)
    {
        ValidateQubit(qubit);
        const size_t mask = 1ULL << qubit;

        double p0 = 0;
        double p1 = 0;
        for (size_t i = 0; i < NrBasisStates; ++i)
        {
            const double population = ValidatedPopulation(i);
            if ((i & mask) == 0)
                p0 += population;
            else
                p1 += population;
        }

        const double total = p0 + p1;
        if (!std::isfinite(total) || total <= 1E-20)
            throw std::domain_error("Cannot measure a state with no probability mass");

        const double r = std::min(uniformZeroOne(rng) * total, std::nextafter(total, 0.));
        const size_t result = (r < p0) ? 0 : 1;
        const double pm = (result == 0) ? p0 : p1;

        CollapseQubit(qubit, result, pm);

        return result;
    }

    // non selective measurement of a qubit: rho' = P0 rho P0 + P1 rho P1
    // destroys the coherence between the two measurement sectors but keeps the populations
    void DephaseMeasure(size_t qubit)
    {
        ValidateQubit(qubit);
        ScaleCoherences(qubit, 0.);
    }

    // sample a full computational basis outcome from the diagonal populations without collapsing
    // the state - useful for repeated sampling that avoids re-executing the circuit each time.
    // The returned value is the measured basis state (bit k corresponds to qubit k).
    size_t MeasureNoCollapse()
    {
        // Reuse capacity, but rebuild every time: even protected state mutations by
        // subclasses must be reflected without requiring cache-invalidation hooks.
        const double total = BuildCumulativeProbabilities(samplingWorkspace);
        const double probability = std::min(uniformZeroOne(rng) * total, std::nextafter(total, 0.));
        const auto it = std::upper_bound(samplingWorkspace.begin(), samplingWorkspace.end(), probability);
        return it == samplingWorkspace.end() ? NrBasisStates - 1 : static_cast<size_t>(it - samplingWorkspace.begin());
    }

    PreparedSampler PrepareSampler() const
    {
        std::vector<double> cumulative;
        const double total = BuildCumulativeProbabilities(cumulative);
        return PreparedSampler(NrBasisStates, 0, NrQubits, std::move(cumulative), total);
    }

    // Packed marginal outcomes: firstQubit becomes bit zero. This opt-in sampler
    // has the same distribution as sampling the full state and masking, but a different
    // mapping from random draws to outcomes. Existing measurement APIs keep their sequence.
    PreparedSampler PrepareSampler(size_t firstQubit, size_t secondQubit) const
    {
        ValidateMeasurementRange(firstQubit, secondQubit);
        const size_t outcomes = size_t{1} << (secondQubit - firstQubit + 1);
        std::vector<double> cumulative(outcomes, 0.);
        for (size_t state = 0; state < NrBasisStates; ++state)
            cumulative[(state >> firstQubit) & (outcomes - 1)] += ValidatedPopulation(state);
        double total = 0.;
        for (double &value : cumulative)
        {
            total += value;
            value = total;
        }
        ValidateProbabilityMass(total);
        return PreparedSampler(NrBasisStates, firstQubit, secondQubit - firstQubit + 1, std::move(cumulative), total);
    }

    // Draw from the snapshot using this simulator's RNG. Any valid snapshot from a
    // register of the same size is accepted; its qubit metadata describes the returned bits.
    size_t MeasureNoCollapse(const PreparedSampler &sampler)
    {
        if (!sampler.isValid() || sampler.sourceBasisStates != NrBasisStates)
            throw std::invalid_argument("Prepared sampler is empty or does not match the register");
        const double probability = std::min(uniformZeroOne(rng) * sampler.total, std::nextafter(sampler.total, 0.));
        const auto &cumulative = sampler.cumulativeProbabilities;
        const auto it = std::upper_bound(cumulative.begin(), cumulative.end(), probability);
        return it == cumulative.end() ? cumulative.size() - 1 : static_cast<size_t>(it - cumulative.begin());
    }

    // sample a subset of qubits from the diagonal populations without collapsing the state.
    // Internally a full basis outcome is drawn from the exact marginal distribution and only the
    // requested qubits are reported. The map keys are the qubit indices, the values the outcomes.
    std::unordered_map<size_t, bool> MeasureNoCollapse(const std::set<size_t> &qubits)
    {
        std::unordered_map<size_t, bool> res;
        if (qubits.empty())
            return res;

        for (const size_t qubit : qubits)
            ValidateQubit(qubit);

        const size_t state = MeasureNoCollapse();
        res.reserve(qubits.size());
        for (const size_t qubit : qubits)
            res[qubit] = (state & (1ULL << qubit)) != 0;

        return res;
    }

    // repeatedly sample the full computational basis distribution without collapsing the state.
    // The cumulative distribution is built once, then reused for every shot.
    std::map<size_t, size_t> RepeatedMeasure(size_t nrTimes = 1000)
    {
        return RepeatedMeasureImpl<std::map<size_t, size_t>>(0, NrBasisStates - 1, nrTimes);
    }

    std::unordered_map<size_t, size_t> RepeatedMeasureUnordered(size_t nrTimes = 1000)
    {
        return RepeatedMeasureImpl<std::unordered_map<size_t, size_t>>(0, NrBasisStates - 1, nrTimes);
    }

    // repeatedly sample a contiguous subregister. The returned outcomes are packed so firstQubit
    // becomes bit zero, matching the statevector simulator's RepeatedMeasure overloads.
    std::map<size_t, size_t> RepeatedMeasure(size_t firstQubit, size_t secondQubit, size_t nrTimes = 1000)
    {
        ValidateMeasurementRange(firstQubit, secondQubit);
        const size_t firstPartMask = (1ULL << firstQubit) - 1;
        const size_t measuredPartMask = (1ULL << (secondQubit + 1)) - 1 - firstPartMask;
        return RepeatedMeasureImpl<std::map<size_t, size_t>>(firstQubit, measuredPartMask, nrTimes);
    }

    std::unordered_map<size_t, size_t> RepeatedMeasureUnordered(size_t firstQubit, size_t secondQubit,
                                                                size_t nrTimes = 1000)
    {
        ValidateMeasurementRange(firstQubit, secondQubit);
        const size_t firstPartMask = (1ULL << firstQubit) - 1;
        const size_t measuredPartMask = (1ULL << (secondQubit + 1)) - 1 - firstPartMask;
        return RepeatedMeasureImpl<std::unordered_map<size_t, size_t>>(firstQubit, measuredPartMask, nrTimes);
    }

    // ---- diagnostics ----

    std::complex<double> Trace() const
    {
        return rho.trace();
    }

    // Tr(rho^2), equals 1 for a pure state, < 1 for a mixed one
    double Purity() const
    {
        // For a Hermitian density matrix Tr(rho^2) is its squared Frobenius norm.
        return rho.squaredNorm();
    }

    bool IsHermitian(double eps = 1E-10) const
    {
        return (rho - rho.adjoint()).norm() < eps;
    }

    // Computes reduced density matrix rho_A = Tr_B(rho) for qubits in keepQubits
    MatrixClass PartialTrace(const std::vector<size_t> &keepQubits) const
    {
        const size_t numKeep = keepQubits.size();
        if (numKeep > NrQubits)
            throw std::invalid_argument("Keep qubits set size exceeds total qubits");

        std::vector<bool> isKept(NrQubits, false);
        for (size_t q : keepQubits)
        {
            if (q >= NrQubits)
                throw std::invalid_argument("Qubit index out of bounds");
            if (isKept[q])
                throw std::invalid_argument("Duplicate qubit index in keepQubits");
            isKept[q] = true;
        }

        const size_t dimA = 1ULL << numKeep;
        const size_t dimB = NrBasisStates / dimA;
        std::vector<size_t> keptOffsets(dimA, 0), tracedOffsets(dimB, 0);
        for (size_t k = 0; k < numKeep; ++k)
            for (size_t i = 0; i < (size_t{1} << k); ++i)
                keptOffsets[i | (size_t{1} << k)] = keptOffsets[i] | (size_t{1} << keepQubits[k]);
        size_t filled = 1;
        for (size_t q = 0; q < NrQubits; ++q)
            if (!isKept[q])
            {
                for (size_t i = 0; i < filled; ++i)
                    tracedOffsets[i + filled] = tracedOffsets[i] | (size_t{1} << q);
                filled *= 2;
            }

        MatrixClass rhoA(dimA, dimA);
        // Each output is owned by one worker, with a fixed summation order and no atomics.
        // Visit exactly dimA^2 * dimB contributing entries, in the requested qubit order.
        ForEachRange(dimA * dimA, NrBasisStates * dimA, [&](size_t begin, size_t end) {
            for (size_t i = begin; i < end; ++i)
            {
                const size_t row = keptOffsets[MatrixClass::IsRowMajor ? i / dimA : i % dimA];
                const size_t col = keptOffsets[MatrixClass::IsRowMajor ? i % dimA : i / dimA];
                Scalar value = 0.;
                for (size_t offset : tracedOffsets)
                    value += rho(row | offset, col | offset);
                rhoA.data()[i] = value;
            }
        });

        return rhoA;
    }

    // Tr(rho_1^\dagger rho_2) = Tr(rho_1 rho_2) Hilbert-Schmidt inner product / state overlap
    std::complex<double> HilbertSchmidtOverlap(const DensityMatrix<VectorClass, MatrixClass> &other) const
    {
        if (other.NrQubits != NrQubits)
            throw std::invalid_argument("Register dimensions do not match");

        return rho.cwiseProduct(other.rho.conjugate()).sum();
    }

    // <psi|rho|psi> / Tr(rho) fidelity with pure statevector psi
    double FidelityWithStatevector(const VectorClass &psi) const
    {
        if (psi.size() < 0 || static_cast<size_t>(psi.size()) != NrBasisStates)
            throw std::invalid_argument("Statevector dimension does not match the register");

        const std::complex<double> tr = Trace();
        if (std::abs(tr) < std::numeric_limits<double>::epsilon())
            return 0.;

        const std::complex<double> val = (psi.adjoint() * (rho * psi))(0);
        return std::clamp((val / tr).real(), 0., 1.);
    }

    double getBasisStateProbability(size_t State) const
    {
        if (State >= NrBasisStates)
            return 0;

        return rho(State, State).real();
    }

    // <O> = Tr(rho O), the caller should ensure O is Hermitian and take the real part
    std::complex<double> ExpectationValue(const MatrixClass &O) const
    {
        if (O.rows() != rho.rows() || O.cols() != rho.cols())
            throw std::invalid_argument("Observable dimensions do not match the density matrix");

        return rho.cwiseProduct(O.transpose()).sum();
    }

    // <P> = Tr(rho P) for a Pauli string P = (x) P_i, where character i of the string is
    // the single qubit Pauli ('I', 'X', 'Y' or 'Z') acting on qubit i (qubit 0 is the
    // rightmost / least significant bit, character 0 in the string).
    //
    // This avoids the full 2^N x 2^N matrix product used by the operator overload: a single
    // qubit Pauli P_i maps a basis state |i> to a phase times a single other basis state, so
    // P|k> = phase(k) |k ^ flipMask> and Tr(rho P) = sum_k phase(k) rho(k ^ flipMask, k),
    // which only touches 2^N entries of rho.
    std::complex<double> ExpectationValue(const std::string &pauliString) const
    {
        if (pauliString.size() != NrQubits)
            throw std::invalid_argument("Pauli string length must match the number of qubits");

        size_t flipMask = 0; // qubits where X or Y act (the bit flipped by P)
        size_t signMask = 0; // qubits where Y or Z act (contribute a (-1)^bit sign)
        size_t yCount = 0;   // number of Y factors (each contributes a factor of i)

        for (size_t i = 0; i < NrQubits; ++i)
        {
            const size_t bit = 1ULL << i;
            switch (toupper(static_cast<unsigned char>(pauliString[i])))
            {
            case 'I':
                break;
            case 'X':
                flipMask |= bit;
                break;
            case 'Y':
                flipMask |= bit;
                signMask |= bit;
                ++yCount;
                break;
            case 'Z':
                signMask |= bit;
                break;
            default:
                throw std::invalid_argument("Invalid operator in the Pauli string");
            }
        }

        std::complex<double> result = 0.;
        for (size_t k = 0; k < NrBasisStates; ++k)
        {
            // parity of the sign-contributing bits set in k
            size_t s = signMask & k;
            bool negative = false;
            while (s)
            {
                negative = !negative;
                s &= s - 1;
            }

            const std::complex<double> term = rho(k, k ^ flipMask);
            result += negative ? -term : term;
        }

        // each Y contributed a factor of i, fold i^yCount in at the end
        static const std::complex<double> iPow[4] = {{1., 0.}, {0., 1.}, {-1., 0.}, {0., -1.}};
        result *= iPow[yCount & 3];

        return result;
    }

  protected:
    template <class Measurements>
    Measurements RepeatedMeasureImpl(size_t firstQubit, size_t measuredPartMask, size_t nrTimes)
    {
        Measurements measurements;
        if (nrTimes == 0)
            return measurements;

        if (nrTimes == 1)
        {
            const size_t state = MeasureNoCollapse();
            ++measurements[(state & measuredPartMask) >> firstQubit];
            return measurements;
        }

        auto &cumulativeProbabilities = samplingWorkspace;
        const double total = BuildCumulativeProbabilities(cumulativeProbabilities);
        const double upperEndpoint = std::nextafter(total, 0.);
        const size_t outcomes = (measuredPartMask >> firstQubit) + 1;
        // Bound auxiliary memory, and avoid clearing a large histogram for a few shots.
        std::vector<size_t> counts;
        if (outcomes <= 65536 && outcomes <= nrTimes)
            counts.resize(outcomes, 0);

        for (size_t shot = 0; shot < nrTimes; ++shot)
        {
            const double probability = std::min(uniformZeroOne(rng) * total, upperEndpoint);
            const auto it =
                std::upper_bound(cumulativeProbabilities.begin(), cumulativeProbabilities.end(), probability);
            const size_t state = it == cumulativeProbabilities.end()
                                     ? NrBasisStates - 1
                                     : static_cast<size_t>(it - cumulativeProbabilities.begin());
            const size_t outcome = (state & measuredPartMask) >> firstQubit;
            if (counts.empty())
                ++measurements[outcome];
            else
                ++counts[outcome];
        }
        for (size_t outcome = 0; outcome < counts.size(); ++outcome)
            if (counts[outcome])
                measurements.emplace(outcome, counts[outcome]);

        return measurements;
    }

    static size_t CheckedBasisStateCount(size_t nrQubits)
    {
        if (nrQubits == 0)
            throw std::invalid_argument("Qubit number must be positive");
        if (nrQubits >= std::numeric_limits<size_t>::digits)
            throw std::invalid_argument("Qubit number is too large for basis-state indexing");

        const size_t nrBasisStates = size_t{1} << nrQubits;
        if (nrBasisStates > static_cast<size_t>(std::numeric_limits<Eigen::Index>::max()))
            throw std::length_error("Register dimension exceeds Eigen's index range");
        const size_t maxElements = std::numeric_limits<size_t>::max() / sizeof(typename MatrixClass::Scalar);
        if (nrBasisStates > maxElements / nrBasisStates)
            throw std::length_error("Density-matrix storage size overflows size_t");
        if (nrBasisStates > static_cast<size_t>(std::numeric_limits<Eigen::Index>::max()) / nrBasisStates)
            throw std::length_error("Density-matrix element count exceeds Eigen's index range");

        return nrBasisStates;
    }

    static size_t GetOperatorQubits(const MatrixClass &op, size_t maxQubits = 3)
    {
        if (op.rows() <= 0 || op.cols() <= 0 || op.rows() != op.cols())
            throw std::invalid_argument("Operator must be a non-empty square matrix");

        const size_t dimension = static_cast<size_t>(op.rows());
        if ((dimension & (dimension - 1)) != 0)
            throw std::invalid_argument("Operator dimension must be a power of two");

        size_t qubits = 0;
        for (size_t value = dimension; value > 1; value >>= 1)
            ++qubits;
        if (qubits == 0 || qubits > maxQubits)
            throw std::invalid_argument("Operator acts on an unsupported number of qubits");

        return qubits;
    }

    size_t ValidateGateAndQubits(const GateClass &gate, size_t qubit, size_t controllingQubit1,
                                 size_t controllingQubit2) const
    {
        const MatrixClass &op = gate.getRawOperatorMatrix();
        const size_t gateQubits = GetOperatorQubits(op);
        if (gate.getQubitsNumber() != gateQubits)
            throw std::invalid_argument("Gate arity does not match its operator matrix");
        if (!op.allFinite())
            throw std::invalid_argument("Gate operator contains a non-finite value");

        ValidateQubits(gateQubits, qubit, controllingQubit1, controllingQubit2);
        return gateQubits;
    }

    void ValidateQubits(size_t gateQubits, size_t qubit, size_t controllingQubit1, size_t controllingQubit2) const
    {
        ValidateQubit(qubit);
        if (gateQubits >= 2)
        {
            ValidateQubit(controllingQubit1);
            if (qubit == controllingQubit1)
                throw std::invalid_argument("Gate qubits must be distinct");
        }
        if (gateQubits == 3)
        {
            ValidateQubit(controllingQubit2);
            if (qubit == controllingQubit2 || controllingQubit1 == controllingQubit2)
                throw std::invalid_argument("Gate qubits must be distinct");
        }
    }

    void ValidateQubit(size_t qubit) const
    {
        if (qubit >= NrQubits)
            throw std::invalid_argument("Qubit number is outside the register");
    }

    void ValidateMeasurementRange(size_t firstQubit, size_t secondQubit) const
    {
        if (firstQubit > secondQubit)
            throw std::invalid_argument("First measured qubit must not exceed the second");
        ValidateQubit(secondQubit);
    }

    static void ValidateProbability(double probability, const char *name)
    {
        if (!std::isfinite(probability) || probability < 0. || probability > 1.)
            throw std::invalid_argument(std::string(name) + " must be finite and in [0, 1]");
    }

    double ValidatedPopulation(size_t state) const
    {
        const auto diagonal = rho(state, state);
        const double population = diagonal.real();
        if (!std::isfinite(population) || !std::isfinite(diagonal.imag()) || std::abs(diagonal.imag()) > 1E-10)
            throw std::domain_error("Density-matrix populations must be finite and real");
        if (population < -1E-12)
            throw std::domain_error("Density-matrix populations must be non-negative");
        return std::max(0., population);
    }

    double BuildCumulativeProbabilities(std::vector<double> &cumulativeProbabilities) const
    {
        cumulativeProbabilities.resize(NrBasisStates);
        double cumulativeProbability = 0.;
        for (size_t state = 0; state < NrBasisStates; ++state)
        {
            cumulativeProbability += ValidatedPopulation(state);
            cumulativeProbabilities[state] = cumulativeProbability;
        }

        ValidateProbabilityMass(cumulativeProbability);
        return cumulativeProbability;
    }

    static void ValidateProbabilityMass(double mass)
    {
        if (!std::isfinite(mass) || mass <= 1E-20)
            throw std::domain_error("Cannot sample a state with no probability mass");
    }

    bool UseParallel(size_t elements) const
    {
        return GetMultithreading() && elements >= GetParallelMinElements();
    }

    // Disjoint tiles share one team. Calls from an existing OpenMP team remain serial.
    template <class Function> void ForEachRange(size_t count, size_t work, const Function &function) const
    {
#ifdef _OPENMP
        if (UseParallel(work) && !omp_in_parallel())
        {
            const int threads = omp_get_max_threads();
            if (threads > 1 && count > 1)
            {
                size_t tileSize = std::min<size_t>(1024, count);
                while (tileSize > 1 && count / tileSize < static_cast<size_t>(threads) * 4)
                    tileSize >>= 1;
                const size_t tiles = (count + tileSize - 1) / tileSize;
#pragma omp parallel for schedule(static) num_threads(threads)
                for (long long tile = 0; tile < static_cast<long long>(tiles); ++tile)
                {
                    const size_t begin = static_cast<size_t>(tile) * tileSize;
                    function(begin, std::min(begin + tileSize, count));
                }
                return;
            }
        }
#endif
        function(0, count);
    }

    // Enumerate independent d-by-d local density blocks. The small block uses
    // column-major local indices (row + d * column), for either physical storage order.
    template <unsigned Qubits, class Function>
    void ForEachLocalBlock(const std::array<size_t, 3> &qubits, const Function &function)
    {
        constexpr size_t dimension = size_t{1} << Qubits;
        std::array<size_t, 2 * Qubits> fixedBits{};
        std::array<size_t, dimension * dimension> offsets{};
        const size_t ketShift = MatrixClass::IsRowMajor ? NrQubits : 0;
        const size_t braShift = MatrixClass::IsRowMajor ? 0 : NrQubits;
        for (unsigned k = 0; k < Qubits; ++k)
        {
            fixedBits[k] = size_t{1} << (qubits[k] + ketShift);
            fixedBits[k + Qubits] = size_t{1} << (qubits[k] + braShift);
        }
        for (size_t i = 0; i < offsets.size(); ++i)
            for (unsigned k = 0; k < 2 * Qubits; ++k)
                if (i & (size_t{1} << k))
                    offsets[i] |= fixedBits[k];
        std::sort(fixedBits.begin(), fixedBits.end());
        for (auto &bit : fixedBits)
            --bit; // masks below each fixed bit
        const size_t elements = static_cast<size_t>(rho.size());
        Scalar *const data = rho.data();
        ForEachRange(elements >> (2 * Qubits), elements, [&](size_t begin, size_t end) {
            for (size_t i = begin; i < end; ++i)
            {
                size_t base = i;
                for (size_t mask : fixedBits)
                    base = (base & mask) | ((base & ~mask) << 1);
                function(data + base, offsets);
            }
        });
    }

    template <unsigned Qubits>
    void ApplyLocalChannel(const std::vector<MatrixClass> &kraus, const std::array<size_t, 2> &qubits)
    {
        constexpr size_t dimension = size_t{1} << Qubits;
        constexpr size_t blockSize = dimension * dimension;
        std::array<std::array<Scalar, blockSize>, blockSize> coefficients{};
        std::array<std::array<unsigned, blockSize>, blockSize> sources{};
        std::array<unsigned, blockSize> counts{};
        // Compile only the local 4x4 or 16x16 map. Exact zeros are omitted; small
        // nonzero terms are retained. Scratch never scales with the register size.
        for (size_t output = 0; output < blockSize; ++output)
            for (size_t input = 0; input < blockSize; ++input)
            {
                Scalar value = 0.;
                for (const auto &e : kraus)
                    value +=
                        e(output % dimension, input % dimension) * std::conj(e(output / dimension, input / dimension));
                if (value != Scalar(0.))
                {
                    const unsigned pos = counts[output]++;
                    coefficients[output][pos] = value;
                    sources[output][pos] = static_cast<unsigned>(input);
                }
            }
        ForEachLocalBlock<Qubits>({qubits[0], qubits[1], 0}, [&](Scalar *base, const auto &offsets) {
            std::array<Scalar, size_t{1} << (2 * Qubits)> input;
            for (size_t i = 0; i < blockSize; ++i)
                input[i] = base[offsets[i]];
            for (size_t output = 0; output < blockSize; ++output)
            {
                Scalar value = 0.;
                for (unsigned k = 0; k < counts[output]; ++k)
                    value += coefficients[output][k] * input[sources[output][k]];
                base[offsets[output]] = value;
            }
        });
    }

    void ScaleCoherences(size_t qubit, double scale)
    {
        if (scale == 1.)
            return;
        ForEachLocalBlock<1>({qubit, 0}, [&](Scalar *base, const auto &o) {
            if (scale == 0.)
                base[o[1]] = base[o[2]] = 0.;
            else
            {
                base[o[1]] *= scale;
                base[o[2]] *= scale;
            }
        });
    }

    // One operator, one traversal: multiply rho(r,c) by u(r)*conjugate(u(c)).
    // This also handles non-unitary diagonal operators without assuming unit-modulus entries.
    template <unsigned Qubits> void ApplyDensityDiagonal(const MatrixClass &matrix, const std::array<size_t, 3> &qubits)
    {
        constexpr size_t dimension = size_t{1} << Qubits;
        std::array<Scalar, dimension * dimension> factors;
        bool identity = true;
        for (size_t c = 0; c < dimension; ++c)
            for (size_t r = 0; r < dimension; ++r)
            {
                factors[r + dimension * c] = matrix(r, r) * std::conj(matrix(c, c));
                identity = identity && factors[r + dimension * c] == Scalar(1.);
            }
        if (identity)
            return;
        size_t lowQubit = qubits[0];
        for (unsigned k = 1; k < Qubits; ++k)
            lowQubit = std::min(lowQubit, qubits[k]);
        // Short runs otherwise repeat bit extraction for virtually every entry. Enumerate
        // local blocks instead, amortizing index expansion over their 4, 16, or 64 entries.
        if (lowQubit < 2)
        {
            ForEachLocalBlock<Qubits>(qubits, [&](Scalar *base, const auto &offsets) {
                for (size_t i = 0; i < offsets.size(); ++i)
                    if (factors[i] != Scalar(1.))
                        base[offsets[i]] *= factors[i];
            });
            return;
        }
        const size_t runMask = (size_t{1} << lowQubit) - 1;
        const size_t ketShift = MatrixClass::IsRowMajor ? NrQubits : 0;
        const size_t braShift = MatrixClass::IsRowMajor ? 0 : NrQubits;
        const size_t count = static_cast<size_t>(rho.size());
        Scalar *const data = rho.data();
        ForEachRange(count, count, [&](size_t begin, size_t end) {
            for (size_t i = begin; i < end;)
            {
                size_t r = 0, c = 0;
                for (unsigned k = 0; k < Qubits; ++k)
                {
                    r |= ((i >> (qubits[k] + ketShift)) & 1) << k;
                    c |= ((i >> (qubits[k] + braShift)) & 1) << k;
                }
                const Scalar factor = factors[r + dimension * c];
                const size_t stop = std::min((i | runMask) + 1, end);
                if (factor == Scalar(1.))
                    i = stop;
                else if (factor == Scalar(-1.))
                    for (; i < stop; ++i)
                        data[i] = -data[i];
                else
                    for (; i < stop; ++i)
                        data[i] *= factor;
            }
        });
    }

    void ApplyGateToBuffer(const MatrixClass &matrix, const Gates::GateStructure &structure, size_t gateQubits,
                           size_t qubit, size_t qubit2, size_t qubit3, bool ket)
    {
        const size_t shift = (ket == bool(MatrixClass::IsRowMajor)) ? NrQubits : 0;
        const std::array<size_t, 3> bits{size_t{1} << (qubit + shift),
                                         gateQubits > 1 ? size_t{1} << (qubit2 + shift) : 0,
                                         gateQubits > 2 ? size_t{1} << (qubit3 + shift) : 0};
        FlatView state(rho.data(), rho.size());
        FlatCalculator::ApplyGateInPlace(state, matrix, structure, bits, static_cast<unsigned>(gateQubits),
                                         static_cast<size_t>(rho.size()), UseParallel(static_cast<size_t>(rho.size())));
    }

    // Retain the protected entry points for subclasses; each is now one full-buffer pass.
    void ApplyGateToColumns(const GateClass &gate, const MatrixClass &gateMatrix, size_t gateQubits, size_t qubit,
                            size_t controllingQubit1, size_t controllingQubit2)
    {
        ApplyGateToBuffer(gateMatrix, gate.getStructure(), gateQubits, qubit, controllingQubit1, controllingQubit2,
                          true);
    }

    // Right multiplication by U^dagger applies conjugate(U) to each row.
    // Even specialized Y/phase/iSWAP paths read that supplied matrix; no
    // replacement gate object or hardcoded conjugation exception is needed.
    void ApplyGateToRows(const GateClass &gate, const MatrixClass &gateMatrix, size_t gateQubits, size_t qubit,
                         size_t controllingQubit1, size_t controllingQubit2)
    {
        ApplyGateToBuffer(gateMatrix, gate.getStructure(), gateQubits, qubit, controllingQubit1, controllingQubit2,
                          false);
    }

    void CollapseQubit(size_t qubit, size_t result, double pm)
    {
        const double invpm = (pm > 1E-20) ? 1. / pm : 0.;
        ForEachLocalBlock<1>({qubit, 0}, [&](Scalar *base, const auto &o) {
            base[o[1]] = base[o[2]] = 0.;
            if (result == 0)
            {
                base[o[0]] *= invpm;
                base[o[3]] = 0.;
            }
            else
            {
                base[o[3]] *= invpm;
                base[o[0]] = 0.;
            }
        });
    }

    size_t NrQubits;
    size_t NrBasisStates;

    MatrixClass rho;
    MatrixClass savedStateStorage;
    std::vector<double> samplingWorkspace; // capacity reuse only; values are never treated as a cached distribution

    bool enableMultithreading = true;

    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniformZeroOne;
};

} // namespace QC
