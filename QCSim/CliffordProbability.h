#pragma once

#include "CliffordTableau.h"
#include <cmath>
#include <random>

namespace QC
{
namespace Clifford
{
namespace detail
{

// Computational-basis support is an affine binary space. Pure-Z stabilizers
// supply its parity constraints; all supported outcomes have probability 2^-r.
// The state is U|label>, with U given by its inverse Z images; a null label is
// the zero input. U|b> = (U X^b)|0>, so a label only flips each image's sign
// by the parity of its Z part on b.
class BasisDistribution
{
  public:
    void BuildFromInverse(const PackedTableau &inverseZ, const Word *label = nullptr)
    {
        PackedTableau work(inverseZ);
        if (label)
            for (size_t row = 0; row < work.size(); ++row)
                ApplyLabel(work[row], label);
        Build(std::move(work));
    }

    void BuildFromInverse(const PackedTableau &inverseZ, const std::vector<size_t> &qubits, const Word *label = nullptr)
    {
        const size_t n = qubits.size();
        PackedTableau work(n, inverseZ.GetNrQubits()), expressions(n, n);
        std::vector<size_t> freeOutputs, logicalPivots;
        freeOutputs.reserve(n);
        logicalPivots.reserve(n);
        for (size_t output = 0; output < n; ++output)
        {
            auto row = work[output], expression = expressions[output];
            row.CopyFrom(inverseZ[qubits[output]]);
            if (label)
                ApplyLabel(row, label);
            for (size_t p = 0; p < freeOutputs.size(); ++p)
                if (row.X[logicalPivots[p]])
                {
                    row.Multiply(work[freeOutputs[p]]);
                    const auto previous = expressions[freeOutputs[p]];
                    for (size_t w = 0; w < expression.Words(); ++w)
                        expression.X.words[w] ^= previous.X.words[w];
                    expression.PhaseSign ^= bool(previous.PhaseSign);
                }
            size_t word = 0;
            while (word < row.Words() && !row.X.words[word])
                ++word;
            if (word < row.Words())
            {
                // This physical output is the next independent random bit.
                // Record the reduced observable's value as that bit XOR the
                // already eliminated observables' affine expressions.
                logicalPivots.push_back(word * 64 + TrailingZero(row.X.words[word]));
                freeOutputs.push_back(output);
                expression.X[output] = true;
            }
            else
                expression.PhaseSign ^= bool(row.PhaseSign);
        }

        PackedTableau nextConstraints(n - freeOutputs.size(), n), nextBasis(freeOutputs.size(), n);
        std::vector<Word> nextOffset(n / 64 + (n % 64 != 0));
        std::vector<size_t> nextPivots;
        nextPivots.reserve(nextConstraints.size());
        size_t free = 0;
        for (size_t output = 0; output < n; ++output)
        {
            if (free < freeOutputs.size() && freeOutputs[free] == output)
            {
                nextBasis[free++].X[output] = true;
                continue;
            }
            const auto expression = expressions[output];
            auto constraint = nextConstraints[nextPivots.size()];
            constraint.Z[output] = true;
            constraint.PhaseSign = bool(expression.PhaseSign);
            nextPivots.push_back(output);
            if (expression.PhaseSign)
                nextOffset[output / 64] |= Word(1) << (output % 64);
            for (size_t p = 0; p < freeOutputs.size(); ++p)
                if (expression.X[freeOutputs[p]])
                {
                    constraint.Z[freeOutputs[p]] = true;
                    nextBasis[p].X[output] = true;
                }
        }
        // No second elimination is needed: these expressions already solve
        // dependent outputs in terms of earlier random outputs, consuming RNG
        // draws in exactly the same order as sequential measurements.
        constraints.swap(nextConstraints);
        basis.swap(nextBasis);
        offset.swap(nextOffset);
        pivotColumns.swap(nextPivots);
        probability = freeOutputs.size() > 1074 ? 0.0 : std::ldexp(1.0, -static_cast<int>(freeOutputs.size()));
    }

  private:
    static void ApplyLabel(TableauRow<false> row, const Word *label) noexcept
    {
        Word parity = 0;
        for (size_t w = 0; w < row.Words(); ++w)
            parity ^= row.Z.words[w] & label[w];
        row.PhaseSign ^= (Popcount(parity) & 1) != 0;
    }

    void Build(PackedTableau work)
    {
        const size_t n = work.size(), logicalQubits = work.GetNrQubits();
        // Track which physical Z operators produced each inverse image.
        PackedTableau labels;
        const auto ensureLabels = [&] {
            if (!labels.empty())
                return;
            labels = PackedTableau(n, n);
            for (size_t q = 0; q < n; ++q)
                labels[q].Z[q] = true;
        };
        size_t rank = 0;
        for (size_t q = 0; q < logicalQubits && rank < n; ++q)
        {
            size_t pivot = rank;
            while (pivot < n && !work[pivot].X[q])
                ++pivot;
            if (pivot == n)
                continue;
            if (rank != pivot)
            {
                ensureLabels();
                work.SwapRows(rank, pivot);
                labels.SwapRows(rank, pivot);
            }
            for (size_t row = rank + 1; row < n; ++row)
                if (work[row].X[q])
                {
                    ensureLabels();
                    work[row].Multiply(work[rank]);
                    for (size_t w = 0; w < labels[row].Words(); ++w)
                        labels[row].Z.words[w] ^= labels[rank].Z.words[w];
                }
            ++rank;
        }
        PackedTableau constraints(n - rank, n);
        for (size_t row = 0; row < constraints.size(); ++row)
        {
            if (labels.empty())
                constraints[row].Z[rank + row] = true;
            else
                constraints[row].CopyFrom(labels[rank + row]);
            constraints[row].PhaseSign = bool(work[rank + row].PhaseSign);
        }
        Finish(std::move(constraints), rank);
    }

  public:
    bool Contains(size_t state) const noexcept
    {
        // Higher qubits are zero in the size_t overload.
        for (size_t row = 0; row < constraints.size(); ++row)
        {
            const auto r = constraints[row];
            const bool parity = r.Words() && (Popcount(r.Z.words[0] & static_cast<Word>(state)) & 1);
            if (parity != bool(r.PhaseSign))
                return false;
        }
        return true;
    }

    bool Contains(const std::vector<bool> &state) const
    {
        if (constraints.empty())
            return true;
        std::vector<Word> bits(offset.size());
        for (size_t q = 0; q < state.size(); ++q)
            if (state[q])
                bits[q / 64] |= Word(1) << (q % 64);
        for (size_t row = 0; row < constraints.size(); ++row)
        {
            const auto r = constraints[row];
            Word parity = 0;
            for (size_t w = 0; w < bits.size(); ++w)
                parity ^= r.Z.words[w] & bits[w];
            if (bool(Popcount(parity) & 1) != bool(r.PhaseSign))
                return false;
        }
        return true;
    }

    template <class State> double Probability(const State &state) const
    {
        return Contains(state) ? probability : 0.0;
    }

    template <class State> double Log2Probability(const State &state) const
    {
        return Contains(state) ? -static_cast<double>(basis.size()) : -std::numeric_limits<double>::infinity();
    }

    void FlipBit(size_t qubit) noexcept
    {
        // Keep the canonical offset's free bits zero. Only the right-hand
        // sides of constraints containing this output bit change.
        for (size_t row = 0; row < constraints.size(); ++row)
            if (constraints[row].Z[qubit])
            {
                constraints[row].PhaseSign ^= true;
                const size_t pivot = pivotColumns[row];
                offset[pivot / 64] ^= Word(1) << (pivot % 64);
            }
    }

    void FillProbabilities(std::vector<double> &probabilities) const noexcept
    {
        size_t state = offset.empty() ? 0 : static_cast<size_t>(offset[0]);
        probabilities[state] = probability;
        const size_t count = size_t(1) << basis.size();
        for (size_t i = 1; i < count; ++i)
        {
            // Consecutive Gray-code labels differ in one free variable.
            const size_t bit = TrailingZero(i);
            state ^= static_cast<size_t>(basis[bit].X.words[0]);
            probabilities[state] = probability;
        }
    }

    // random(engine) supplies one fair random bit per independent output.
    template <class Engine, class Random> std::vector<bool> Sample(Engine &engine, Random &random) const
    {
        std::vector<Word> bits(offset.size());
        std::vector<bool> result(basis.GetNrQubits());
        SampleInto(bits, engine, random);
        for (size_t q = 0; q < result.size(); ++q)
            result[q] = (bits[q / 64] >> (q % 64)) & 1;
        return result;
    }

    size_t Words() const noexcept
    {
        return offset.size();
    }

    template <class Engine, class Random> void SampleInto(std::vector<Word> &bits, Engine &engine, Random &random) const
    {
        assert(bits.size() == offset.size());
        std::copy(offset.begin(), offset.end(), bits.begin());
        for (size_t row = 0; row < basis.size(); ++row)
            if (random(engine))
                for (size_t w = 0; w < bits.size(); ++w)
                    bits[w] ^= basis[row].X.words[w];
    }

  private:
    void Finish(PackedTableau nextConstraints, size_t dimension)
    {
        const size_t n = nextConstraints.GetNrQubits();
        std::vector<size_t> pivots;
        pivots.reserve(nextConstraints.size());
        // Reduce the diagonal constraints to solve for pivot variables.
        for (size_t q = 0; q < n && pivots.size() < nextConstraints.size(); ++q)
        {
            const size_t next = pivots.size();
            size_t pivot = next;
            while (pivot < nextConstraints.size() && !nextConstraints[pivot].Z[q])
                ++pivot;
            if (pivot == nextConstraints.size())
                continue;
            nextConstraints.SwapRows(next, pivot);
            for (size_t row = 0; row < nextConstraints.size(); ++row)
                if (row != next && nextConstraints[row].Z[q])
                {
                    auto r = nextConstraints[row];
                    const auto p = nextConstraints[next];
                    for (size_t w = 0; w < r.Words(); ++w)
                        r.Z.words[w] ^= p.Z.words[w];
                    r.PhaseSign ^= bool(p.PhaseSign);
                }
            pivots.push_back(q);
        }
        assert(pivots.size() == nextConstraints.size());
        PackedTableau nextBasis(dimension, n);
        std::vector<Word> nextOffset(n / 64 + (n % 64 != 0));
        for (size_t row = 0; row < pivots.size(); ++row)
            if (nextConstraints[row].PhaseSign)
                nextOffset[pivots[row] / 64] |= Word(1) << (pivots[row] % 64);
        size_t free = 0, pivot = 0;
        for (size_t q = 0; q < n; ++q)
        {
            if (pivot < pivots.size() && pivots[pivot] == q)
            {
                ++pivot;
                continue;
            }
            auto r = nextBasis[free++];
            r.X[q] = true;
            for (size_t row = 0; row < pivots.size(); ++row)
                r.X[pivots[row]] = nextConstraints[row].Z[q];
        }
        assert(free == dimension);
        // Publish only after every allocation and elimination has succeeded.
        constraints.swap(nextConstraints);
        basis.swap(nextBasis);
        offset.swap(nextOffset);
        pivotColumns.swap(pivots);
        probability = dimension > 1074 ? 0.0 : std::ldexp(1.0, -static_cast<int>(dimension));
    }

    PackedTableau constraints, basis;
    std::vector<Word> offset;
    std::vector<size_t> pivotColumns;
    double probability = 1.0;
};

} // namespace detail
} // namespace Clifford
} // namespace QC
