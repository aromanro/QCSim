#pragma once

#include "CliffordProbability.h"
#include "PauliStringXZ.h"
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace QC
{
namespace Clifford
{

// The rows represent U^dagger X_q U and U^dagger Z_q U for |psi> = U|0>.
// Measurements rebase this inverse Clifford map while retaining a zero logical
// input. No frame amplitudes or extended-stabilizer approximation are involved.
class StabilizerState
{
  public:
    using Generator = PauliStringXZWithSign;

    StabilizerState() : StabilizerState(0)
    {
    }

    explicit StabilizerState(size_t n)
        : inverseX(n), inverseZ(n), measurementScratch(1, n), gen(std::random_device{}()), rnd(0.5)
    {
        Reset();
    }

    // Copies retain the previous independent-RNG policy; moves transfer the RNG.
    StabilizerState(const StabilizerState &other)
        : inverseX(other.inverseX), inverseZ(other.inverseZ), savedX(other.savedX), savedZ(other.savedZ),
          measurementScratch(1, other.getNrQubits()), gen(std::random_device{}()), rnd(0.5),
          enableMultithreading(other.enableMultithreading)
    {
    }

    StabilizerState(StabilizerState &&other) noexcept : gen(std::move(other.gen)), rnd(std::move(other.rnd))
    {
        SwapQuantumState(other);
    }

    StabilizerState &operator=(const StabilizerState &other)
    {
        if (this != &other)
        {
            StabilizerState copy(other);
            SwapQuantumState(copy);
        }
        return *this;
    }

    StabilizerState &operator=(StabilizerState &&other) noexcept
    {
        if (this != &other)
        {
            SwapQuantumState(other);
            std::swap(gen, other.gen);
            std::swap(rnd, other.rnd);
        }
        return *this;
    }

    void SetSeed(uint64_t seed)
    {
        std::seed_seq sequence{uint32_t(seed), uint32_t(seed >> 32)};
        gen.seed(sequence);
    }

    void Reset() noexcept
    {
        InvalidateDistribution();
        Map().SetIdentity();
    }

    bool MeasureQubit(size_t qubit)
    {
        ValidateQubit(qubit);
        size_t pivot;
        if (!IsRandomResult(qubit, pivot))
            return bool(inverseZ[qubit].PhaseSign);
        const bool outcome = rnd(gen);
        CollapseRandomQubit(qubit, pivot, outcome);
        return outcome;
    }

    double GetQubitProbability(size_t qubit) const
    {
        ValidateQubit(qubit);
        return inverseZ[qubit].HasX() ? 0.5 : (inverseZ[qubit].PhaseSign ? 1.0 : 0.0);
    }

    double getBasisStateProbability(size_t state)
    {
        const size_t n = getNrQubits();
        if (n < std::numeric_limits<size_t>::digits && (state >> n) != 0)
            return 0.0;
        EnsureDistribution();
        return distribution.Probability(state);
    }

    double getBasisStateProbability(const std::vector<bool> &state)
    {
        if (state.size() != getNrQubits())
            throw std::invalid_argument("Basis state must contain one bit per qubit");
        EnsureDistribution();
        return distribution.Probability(state);
    }

    // These queries distinguish impossible outcomes from probabilities too
    // small for double. Unsupported outcomes have log2 probability -infinity.
    bool ContainsBasisState(size_t state)
    {
        const size_t n = getNrQubits();
        if (n < std::numeric_limits<size_t>::digits && (state >> n) != 0)
            return false;
        EnsureDistribution();
        return distribution.Contains(state);
    }

    bool ContainsBasisState(const std::vector<bool> &state)
    {
        if (state.size() != getNrQubits())
            throw std::invalid_argument("Basis state must contain one bit per qubit");
        EnsureDistribution();
        return distribution.Contains(state);
    }

    double Log2BasisStateProbability(size_t state)
    {
        const size_t n = getNrQubits();
        if (n < std::numeric_limits<size_t>::digits && (state >> n) != 0)
            return -std::numeric_limits<double>::infinity();
        EnsureDistribution();
        return distribution.Log2Probability(state);
    }

    double Log2BasisStateProbability(const std::vector<bool> &state)
    {
        if (state.size() != getNrQubits())
            throw std::invalid_argument("Basis state must contain one bit per qubit");
        EnsureDistribution();
        return distribution.Log2Probability(state);
    }

    // Terminal samples preserve the quantum state and advance this simulator's RNG.
    std::vector<bool> SampleBasisState()
    {
        EnsureDistribution();
        return distribution.Sample(gen, rnd);
    }

    std::vector<std::vector<bool>> SampleBasisStates(size_t shots)
    {
        std::vector<std::vector<bool>> samples;
        if (shots == 0)
            return samples;
        EnsureDistribution();
        samples.reserve(shots);
        for (size_t shot = 0; shot < shots; ++shot)
            samples.push_back(distribution.Sample(gen, rnd));
        return samples;
    }

    // Bit i in each key corresponds to qubits[i]. Repeated indices repeat the
    // same measured value. Sampling preserves both the live and saved state.
    std::unordered_map<size_t, size_t> SampleCounts(const std::vector<size_t> &qubits, size_t shots)
    {
        if (qubits.size() > std::numeric_limits<size_t>::digits)
            throw std::invalid_argument("Use SampleCountsMany for outcomes wider than size_t");
        ValidateSampleQubits(qubits);
        std::unordered_map<size_t, size_t> counts;
        if (shots == 0 || qubits.empty())
            return counts;
        ForEachSample(qubits, shots, [&](const auto &bits) { ++counts[static_cast<size_t>(bits[0])]; });
        return counts;
    }

    std::unordered_map<std::vector<bool>, size_t> SampleCountsMany(const std::vector<size_t> &qubits, size_t shots)
    {
        ValidateSampleQubits(qubits);
        std::unordered_map<std::vector<bool>, size_t> counts;
        if (shots == 0 || qubits.empty())
            return counts;
        std::vector<bool> result(qubits.size());
        ForEachSample(qubits, shots, [&](const auto &bits) {
            for (size_t q = 0; q < result.size(); ++q)
                result[q] = (bits[q / 64] >> (q % 64)) & 1;
            ++counts[result];
        });
        return counts;
    }

    std::vector<double> AllProbabilities()
    {
        const size_t nrQubits = getNrQubits();
        if (nrQubits > 32)
            throw std::runtime_error("The simulator has too many qubits for computing all probabilities");
        if (nrQubits >= std::numeric_limits<size_t>::digits)
            throw std::length_error("The simulator has too many qubits for computing all probabilities");

        const size_t nrStates = 1ULL << nrQubits;
        if (nrStates > std::vector<double>().max_size())
            throw std::length_error("Probability vector exceeds container capacity");
        std::vector<double> probs(nrStates, 0);

        EnsureDistribution();
        distribution.FillProbabilities(probs);

        return probs;
    }

    size_t getNrQubits() const noexcept
    {
        return inverseZ.size();
    }

    void SaveState()
    {
        if (savedX.HasSameShape(inverseX) && savedZ.HasSameShape(inverseZ))
        {
            savedX.CopyFrom(inverseX);
            savedZ.CopyFrom(inverseZ);
            return;
        }
        // First save (or save after ClearSavedState): publish both allocations
        // together. Repeated saves above cannot throw or partially fail.
        detail::PackedTableau x(inverseX), z(inverseZ);
        savedX.swap(x);
        savedZ.swap(z);
    }

    void RestoreState() noexcept
    {
        if (savedX.empty())
            return;
        inverseX.CopyFrom(savedX);
        inverseZ.CopyFrom(savedZ);
        InvalidateDistribution();
    }

    void RestoreSavedStateDestructive() noexcept
    {
        if (savedX.empty())
            return;
        inverseX.swap(savedX);
        inverseZ.swap(savedZ);
        ClearSavedState();
        InvalidateDistribution();
    }

    void ClearSavedState() noexcept
    {
        savedX.clear();
        savedZ.clear();
    }

    void SetMultithreading(bool enable = true) noexcept
    {
        enableMultithreading = enable;
    }

    bool GetMultithreading() const noexcept
    {
        return enableMultithreading;
    }

  protected:
    void FlipDistributionBit(size_t qubit) noexcept
    {
        if (validDistributions & FullDistribution)
            distribution.FlipBit(qubit);
        if (validDistributions & MarginalDistribution)
            for (size_t bit = 0; bit < marginalQubits.size(); ++bit)
                if (marginalQubits[bit] == qubit)
                    marginalDistribution.FlipBit(bit);
    }

    // Keep both validity flags in one byte: ordinary gates invalidate both
    // caches with one store, regardless of the amount of cached storage.
    void InvalidateDistribution() noexcept
    {
        validDistributions = 0;
    }

    void EnsureDistribution()
    {
        if (!(validDistributions & FullDistribution))
        {
            distribution.BuildFromInverse(inverseZ);
            validDistributions |= FullDistribution;
        }
    }

    void ValidateSampleQubits(const std::vector<size_t> &qubits) const
    {
        for (size_t q : qubits)
            ValidateQubit(q);
    }

    void EnsureMarginalDistribution(const std::vector<size_t> &qubits)
    {
        if ((validDistributions & MarginalDistribution) && marginalQubits == qubits)
            return;
        std::vector<size_t> nextQubits(qubits);
        marginalDistribution.BuildFromInverse(inverseZ, qubits);
        marginalQubits.swap(nextQubits);
        validDistributions |= MarginalDistribution;
    }

    template <class Consumer> bool TrySimpleSample(const std::vector<size_t> &qubits, Consumer consume)
    {
        size_t first = getNrQubits();
        for (size_t q : qubits)
        {
            const auto row = inverseZ[q];
            if (!row.HasX())
                continue;
            if (first == getNrQubits())
            {
                first = q;
                continue;
            }
            const auto previous = inverseZ[first];
            for (size_t w = 0; w < row.Words(); ++w)
                if (row.X.words[w] != previous.X.words[w])
                    return false;
        }
        // Rank zero/one needs neither a tableau copy nor elimination. Rows
        // with equal X parts differ by a deterministic diagonal observable.
        std::vector<detail::Word> bits(qubits.size() / 64 + (qubits.size() % 64 != 0));
        const bool random = first != getNrQubits() && rnd(gen);
        unsigned referenceY = 0;
        bool referenceSign = false;
        if (first != getNrQubits())
        {
            const auto reference = inverseZ[first];
            referenceSign = bool(reference.PhaseSign);
            for (size_t w = 0; w < reference.Words(); ++w)
                referenceY += detail::Popcount(reference.X.words[w] & reference.Z.words[w]);
        }
        for (size_t output = 0; output < qubits.size(); ++output)
        {
            const auto row = inverseZ[qubits[output]];
            bool value = bool(row.PhaseSign);
            if (row.HasX())
            {
                unsigned phase = referenceY + 2 * unsigned(value != referenceSign);
                for (size_t w = 0; w < row.Words(); ++w)
                    phase += 3 * detail::Popcount(row.X.words[w] & row.Z.words[w]);
                assert((phase & 1) == 0);
                value = random != bool(phase & 2);
            }
            if (value)
                bits[output / 64] |= detail::Word(1) << (output % 64);
        }
        consume(bits);
        return true;
    }

    template <class Consumer> void ForEachSample(const std::vector<size_t> &qubits, size_t shots, Consumer consume)
    {
        // A single wide cold shot need not pay for full Gaussian elimination.
        // The prepared path uses the same random free variables as measurement,
        // so warming the cache never changes seeded outcomes.
        if (shots == 1 && !((validDistributions & MarginalDistribution) && marginalQubits == qubits))
        {
            if (TrySimpleSample(qubits, consume))
                return;
            if (qubits.size() >= 128)
            {
                StabilizerState sample(*this, SamplingCopy{});
                std::vector<detail::Word> bits(qubits.size() / 64 + (qubits.size() % 64 != 0));
                for (size_t q = 0; q < qubits.size(); ++q)
                    if (sample.MeasureQubit(qubits[q]))
                        bits[q / 64] |= detail::Word(1) << (q % 64);
                consume(bits);
                gen = sample.gen;
                rnd = sample.rnd;
                return;
            }
        }
        EnsureMarginalDistribution(qubits);
        std::vector<detail::Word> bits(marginalDistribution.Words());
        for (size_t shot = 0; shot < shots; ++shot)
        {
            marginalDistribution.SampleInto(bits, gen, rnd);
            consume(bits);
        }
    }

    void ValidateQubit(size_t q) const
    {
        if (q >= getNrQubits())
            throw std::out_of_range("Qubit index out of range");
    }

    void ValidatePair(size_t a, size_t b, bool allowEqual = false) const
    {
        ValidateQubit(a);
        ValidateQubit(b);
        if (!allowEqual && a == b)
            throw std::invalid_argument("Two-qubit gate requires distinct qubits");
    }

    detail::InverseMap Map() noexcept
    {
        return {inverseX, inverseZ};
    }

    bool IsRandomResult(size_t qubit, size_t &pivot) const noexcept
    {
        return detail::FindPivot(inverseZ[qubit], pivot);
    }

    void CollapseRandomQubit(size_t qubit, size_t pivot, bool outcome) noexcept
    {
        InvalidateDistribution();
        Map().CollapseZ(qubit, pivot, outcome, measurementScratch[0],
                        detail::CollapseWorkers(getNrQubits(), enableMultithreading));
    }

    void SwapQuantumState(StabilizerState &other) noexcept
    {
        inverseX.swap(other.inverseX);
        inverseZ.swap(other.inverseZ);
        savedX.swap(other.savedX);
        savedZ.swap(other.savedZ);
        measurementScratch.swap(other.measurementScratch);
        std::swap(enableMultithreading, other.enableMultithreading);
        InvalidateDistribution();
        other.InvalidateDistribution();
    }

    detail::PackedTableau inverseX, inverseZ;
    detail::PackedTableau savedX, savedZ;
    detail::PackedTableau measurementScratch;
    std::mt19937_64 gen;
    std::bernoulli_distribution rnd{0.5};
    bool enableMultithreading = true;

    enum : uint8_t
    {
        FullDistribution = 1,
        MarginalDistribution = 2
    };

    uint8_t validDistributions = 0;
    detail::BasisDistribution distribution;
    detail::BasisDistribution marginalDistribution;
    std::vector<size_t> marginalQubits;

  private:
    struct SamplingCopy
    {
    };

    // Copy only the live tableau. Preserve the parent's RNG stream without
    // invoking random_device or copying its saved state and probability caches.
    StabilizerState(const StabilizerState &other, SamplingCopy)
        : inverseX(other.inverseX), inverseZ(other.inverseZ), measurementScratch(1, other.getNrQubits()),
          gen(other.gen), rnd(other.rnd), enableMultithreading(other.enableMultithreading)
    {
    }
};

} // namespace Clifford
} // namespace QC
