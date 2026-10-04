#pragma once

#include "PathIntegralGate.h"
#include "PathIntegralStorage.h"
#include <array>
#include <cmath>
#include <random>
#include <unordered_map>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include "QuantumGate.h"

namespace QC
{
namespace PathIntegral
{

class PathIntegralSimulator
{
  public:
    PathIntegralSimulator()
        : doublingsLimit(std::numeric_limits<size_t>::max()), epsilon(1E-15), rng(std::random_device{}()),
          uniformZeroOne(0, 1)
    {
    }

    void SetSeed(uint64_t theSeed)
    {
        std::seed_seq seed{uint32_t(theSeed & 0xffffffff), uint32_t(theSeed >> 32)};
        rng.seed(seed);
    }

    void SetTrimValue(double val)
    {
        epsilon = val;
    }

    double GetTrimValue() const
    {
        return epsilon;
    }

    void SetMaxDoublingsForBackwardPaths(size_t doublings)
    {
        doublingsLimit = doublings;
    }

    size_t GetMaxDoublingsForBackwardPaths() const
    {
        return doublingsLimit;
    }

    AmplitudeMap &GetAmplitudes()
    {
        return intermediateAmplitudes;
    }

    void SaveAmplitudes()
    {
        savedAmplitudes = intermediateAmplitudes;
    }

    void RestoreAmplitudes()
    {
        intermediateAmplitudes = savedAmplitudes;
    }

    void SwapAmplitudes()
    {
        intermediateAmplitudes.swap(savedAmplitudes);
    }

    void ClearSavedAmplitudes()
    {
        savedAmplitudes.Release();
    }

    void Reset()
    {
        intermediateAmplitudes.Release();
        savedAmplitudes.Release();
        scratchAmplitudes.Release();
        circuit.clear();
        circuitBack.clear();
    }

    PathIntegralSimulator Clone() const
    {
        PathIntegralSimulator theClone;
        theClone.doublingsLimit = doublingsLimit;
        theClone.epsilon = epsilon;
        theClone.intermediateAmplitudes = intermediateAmplitudes;
        theClone.savedAmplitudes = savedAmplitudes;
        theClone.circuit = circuit;
        theClone.circuitBack = circuitBack;

        return theClone;
    }

    void SetCircuit(const std::vector<QC::Gates::AppliedGate<>> &circuit)
    {
        intermediateAmplitudes.clear();
        circuitBack.clear();

        if (doublingsLimit > 0)
        {
            size_t doublings = 0;
            for (const auto &gate : circuit)
            {
                if (!gate.isBranching())
                    continue;

                ++doublings;
            }

            const size_t halfDoublings = doublings / 2;

            // iterate backwards
            doublings = 0;
            size_t gatesAdded = 0;
            for (auto it = circuit.rbegin(); it != circuit.rend(); ++it)
            {
                const auto &gate = *it;

                if (!gate.isBranching())
                {
                    circuitBack.emplace_back(gate.getRawOperatorMatrix().transpose(), gate.getQubit1(),
                                             gate.getQubit2(), gate.getQubit3());
                    ++gatesAdded;
                    continue;
                }

                ++doublings;

                if (doublings > doublingsLimit || doublings > halfDoublings)
                    break;

                circuitBack.emplace_back(gate.getRawOperatorMatrix().transpose(), gate.getQubit1(), gate.getQubit2(),
                                         gate.getQubit3());
                ++gatesAdded;
            }

            this->circuit = std::vector<QC::Gates::AppliedGate<>>(circuit.begin(), circuit.end() - gatesAdded);
        }
        else
            this->circuit = circuit;
    }

    // the start state is |0>
    std::complex<double> Propagate(const std::vector<bool> &endState)
    {
        const std::vector<bool> startState(endState.size(), false);

        return Propagate(startState, endState);
    }

    std::complex<double> Propagate(const std::vector<bool> &startState, const std::vector<bool> &endState)
    {
        assert(startState.size() == endState.size());

        const FastVectorBool startBits(startState);
        const FastVectorBool endBits(endState);

        if (doublingsLimit > 0 && !circuitBack.empty())
        {
            intermediateAmplitudes.clear();

            PropagateBackward(endBits, std::complex<double>(1., 0.));

            return PropagateForward(startBits);
        }

        const size_t possibleQubitsChanges = CountPossibleQubitsChanges(circuit);

        return Propagate(startBits, endBits, 0, possibleQubitsChanges);
    }

    // the following methods are for the case one needs all non-zero amplitudes
    // of course, except the trimmed out ones
    void PropagateAll(const std::vector<QC::Gates::AppliedGate<>> &circuit)
    {
        size_t nQubits = 0;
        for (const auto &gate : circuit)
        {
            nQubits = std::max(nQubits, gate.getQubit1() + 1);
            if (gate.getQubitsNumber() > 1)
                nQubits = std::max(nQubits, gate.getQubit2() + 1);
            if (gate.getQubitsNumber() > 2)
                nQubits = std::max(nQubits, gate.getQubit3() + 1);
        }

        const std::vector<bool> startState(nQubits, false);

        PropagateAll(circuit, startState);
    }

    void PropagateAll(const std::vector<QC::Gates::AppliedGate<>> &circuit, const std::vector<bool> &startState)
    {
        assert(startState.size() > 0);
        const FastVectorBool startBits(startState);
        intermediateAmplitudes.clear();

        circuitBack.clear();
        this->circuit = circuit;

        intermediateAmplitudes = PropagateAll(startBits);
    }

    void PropagateStep(const QC::Gates::AppliedGate<> &gate, AmplitudeMap &currentAmplitudes)
    {
        if (currentAmplitudes.empty())
        {
            scratchAmplitudes.ReleaseIfOversized(0, currentAmplitudes.QubitCount());
            return;
        }
        currentAmplitudes.CompactIfSparse();
        const Detail::SparseGate transitions(gate, epsilon);
        for (unsigned i = 0; i < transitions.arity; ++i)
            if (transitions.qubits[i] >= currentAmplitudes.QubitCount())
                throw std::out_of_range("Path integral gate qubit is outside the register");
        switch (transitions.arity)
        {
        case 1:
            PropagateSparse<1>(transitions, currentAmplitudes);
            break;
        case 2:
            PropagateSparse<2>(transitions, currentAmplitudes);
            break;
        case 3:
            PropagateSparse<3>(transitions, currentAmplitudes);
            break;
        }
        scratchAmplitudes.ReleaseIfOversized(currentAmplitudes.size(), currentAmplitudes.QubitCount());
    }

  private:
    template <unsigned Arity>
    bool UseGroups(const Detail::SparseGate &gate, const AmplitudeMap &current, unsigned &expansion) const
    {
        if (Arity < 2 || gate.maxFanout < 4 || !gate.distinctQubits || current.size() < 128)
            return false;
        constexpr unsigned dimension = 1u << Arity;
        constexpr size_t samples = 8;
        size_t present = 0, outputs = 0;
        std::array<uint64_t, FastVectorBool::MaxWords> words;
        for (size_t sample = 0; sample < samples; ++sample)
        {
            std::copy_n(current.State(sample * (current.size() / samples)).getWords(), current.wordCount,
                        words.begin());
            Detail::MutableStateView state(words.data(), current.QubitCount());
            unsigned rows = 0;
            for (unsigned column = 0; column < dimension; ++column)
            {
                gate.SetRow<Arity>(state, column);
                if (current.FindIndex(state) == current.size())
                    continue;
                ++present;
                for (unsigned t = 0; t < gate.counts[column]; ++t)
                    rows |= 1u << gate.rows[column][t];
            }
            for (; rows; rows &= rows - 1)
                ++outputs;
        }
        // Sparse local groups really can expand four- or eightfold. Do
        // not force them through repeated allocation/rehash growth.
        expansion = static_cast<unsigned>(std::max<size_t>(1, (outputs + present - 1) / present));
        // The sample is only an estimate. Round up to leave room for
        // less-occupied groups without a mid-gate growth allocation.
        unsigned rounded = 1;
        while (rounded < expansion)
            rounded *= 2;
        expansion = rounded;
        return present * 2 >= samples * dimension;
    }

    template <unsigned Arity>
    void PropagateGroups(const Detail::SparseGate &gate, AmplitudeMap &current, unsigned expansion)
    {
        constexpr unsigned dimension = 1u << Arity;
        auto &next = scratchAmplitudes;
        next.PrepareOutput(current.QubitCount(),
                           gate.ReserveSize(current.size(), current.QubitCount(), next.max_size(), expansion));
        std::vector<unsigned char> visited(current.size(), 0);
        std::array<uint64_t, FastVectorBool::MaxWords> words;
        for (size_t i = 0; i < current.size(); ++i)
        {
            if (visited[i])
                continue;
            std::copy_n(current.State(i).getWords(), current.wordCount, words.begin());
            Detail::MutableStateView state(words.data(), current.QubitCount());
            const unsigned sourceColumn = gate.Column<Arity>(state);
            std::array<std::complex<double>, dimension> output{};
            unsigned retainedRows = 0;
            for (unsigned column = 0; column < dimension; ++column)
            {
                gate.SetRow<Arity>(state, column);
                const size_t index = column == sourceColumn ? i : current.FindIndex(state);
                if (index == current.size())
                    continue;
                visited[index] = 1;
                const auto amplitude = current.values[index];
                if (std::norm(amplitude) < epsilon)
                    continue;
                for (unsigned t = 0; t < gate.counts[column]; ++t)
                {
                    const unsigned row = gate.rows[column][t];
                    output[row] += gate.values[column][t] * amplitude;
                    retainedRows |= 1u << row;
                }
            }
            for (unsigned row = 0; row < dimension; ++row)
            {
                if (!(retainedRows & (1u << row)))
                    continue;
                gate.SetRow<Arity>(state, row);
                // Distinct untouched-bit groups cannot produce the same key.
                // Keep entries even when their contributions cancel to zero.
                next.InsertUnique(state, output[row]);
            }
        }
        current.swap(next);
        current.CompactIfSparse();
    }

    template <unsigned Arity> void PropagateSparse(const Detail::SparseGate &gate, AmplitudeMap &current)
    {
        if (gate.diagonal || gate.injective)
        {
            current.Transform(!gate.diagonal,
                              [&](Detail::MutableStateView &state, std::complex<double> &amplitude) noexcept {
                                  if (std::norm(amplitude) < epsilon)
                                      return false;
                                  const unsigned column = gate.Column<Arity>(state);
                                  if (!gate.counts[column])
                                      return false;
                                  if (!gate.diagonal)
                                      gate.SetRow<Arity>(state, gate.rows[column][0]);
                                  amplitude = gate.values[column][0] * amplitude;
                                  return true;
                              });
            return;
        }
        unsigned expansion = 0;
        if (UseGroups<Arity>(gate, current, expansion))
        {
            PropagateGroups<Arity>(gate, current, expansion);
            return;
        }
        // Retain capacities in both buffers across gates. No state nodes or
        // per-key allocations are needed, including for wide registers.
        auto &next = scratchAmplitudes;
        next.PrepareOutput(current.QubitCount(),
                           gate.ReserveSize(current.size(), current.QubitCount(), next.max_size(), expansion));
        std::array<uint64_t, FastVectorBool::MaxWords> words;
        const size_t width = current.QubitCount();
        const size_t wordCount = (width + 63) / 64;
        for (const auto &item : current)
        {
            if (std::norm(item.second) < epsilon)
                continue;
            const unsigned column = gate.Column<Arity>(item.first);
            std::copy_n(item.first.getWords(), wordCount, words.begin());
            Detail::MutableStateView state(words.data(), width);
            for (unsigned t = 0; t < gate.counts[column]; ++t)
            {
                gate.SetRow<Arity>(state, gate.rows[column][t]);
                next[state] += gate.values[column][t] * item.second;
            }
        }
        current.swap(next);
        current.CompactIfSparse();
    }

  public:
    double QubitProbability(size_t qubit, bool value = true) const
    {
        double prob = 0.;
        for (const auto &[state, amp] : intermediateAmplitudes)
        {
            if (state.get(qubit) == value)
                prob += std::norm(amp);
        }
        return prob;
    }

    FastVectorBool MeasureNoCollapse()
    {
        const auto &values = intermediateAmplitudes.values;
        if (values.empty())
            throw std::domain_error("Cannot measure an empty path integral state");
        // Rebuild from current values: callers can retain mutable amplitude
        // references, so a persistent cached distribution could be stale.
        // Independent sums scan all amplitudes without a per-entry CDF
        // dependency. Only the selected block needs a second, short scan.
        std::array<double, 64> cumulativeBlocks;
        if (values.size() <= cumulativeBlocks.size())
        {
            double total = 0.;
            for (size_t i = 0; i < values.size(); ++i)
            {
                total += std::norm(values[i]);
                cumulativeBlocks[i] = total;
            }
            const double target = MeasurementTarget(total);
            size_t index = 0;
            while (index + 1 < values.size() && target >= cumulativeBlocks[index])
                ++index;
            return intermediateAmplitudes.State(index);
        }
        const size_t blockSize = std::max<size_t>(64, 1 + (values.size() - 1) / cumulativeBlocks.size());
        const size_t blocks = 1 + (values.size() - 1) / blockSize;
        double total = 0.;
        for (size_t block = 0; block < blocks; ++block)
        {
            const size_t end = std::min(values.size(), (block + 1) * blockSize);
            size_t i = block * blockSize;
            double a = 0., b = 0., c = 0., d = 0.;
            for (; i + 3 < end; i += 4)
            {
                a += std::norm(values[i]);
                b += std::norm(values[i + 1]);
                c += std::norm(values[i + 2]);
                d += std::norm(values[i + 3]);
            }
            for (; i < end; ++i)
                a += std::norm(values[i]);
            total += (a + b) + (c + d);
            cumulativeBlocks[block] = total;
        }
        double target = MeasurementTarget(total);
        size_t block = 0;
        while (block + 1 < blocks && target >= cumulativeBlocks[block])
            ++block;
        if (block)
            target -= cumulativeBlocks[block - 1];
        double cumulative = 0.;
        size_t lastPositive = block * blockSize;
        const size_t end = std::min(values.size(), (block + 1) * blockSize);
        for (size_t i = block * blockSize; i < end; ++i)
        {
            const double weight = std::norm(values[i]);
            if (weight > 0.)
                lastPositive = i;
            cumulative += weight;
            if (target < cumulative)
                return intermediateAmplitudes.State(i);
        }
        return intermediateAmplitudes.State(lastPositive);
    }

    bool MeasureQubit(size_t qubit)
    {
        if (intermediateAmplitudes.empty())
            throw std::domain_error("Cannot measure an empty path integral state");
        if (qubit >= intermediateAmplitudes.QubitCount())
            throw std::out_of_range("Path integral measurement qubit is outside the register");
        // Draw directly from the two marginals instead of sampling a full
        // basis state, then compact survivors and rebuild their index once.
        double zeroMass = 0., oneMass = 0.;
        unsigned present = 0;
        const size_t word = qubit / 64;
        const uint64_t mask = uint64_t{1} << (qubit % 64);
        for (size_t i = 0; i < intermediateAmplitudes.size(); ++i)
        {
            const bool bit = (intermediateAmplitudes.keys[i * intermediateAmplitudes.wordCount + word] & mask) != 0;
            const double weight = std::norm(intermediateAmplitudes.values[i]);
            if (bit)
                oneMass += weight;
            else
                zeroMass += weight;
            present |= 1u << unsigned(bit);
        }
        const bool result = MeasurementTarget(zeroMass + oneMass) >= zeroMass;
        const double normFactor = 1. / std::sqrt(result ? oneMass : zeroMass);
        // A normalized state with this bit already fixed needs no writes.
        // Track keys as well as mass so opposite-sector zero entries still
        // get removed when that sector has zero probability.
        if (present == 3 || normFactor != 1.)
            intermediateAmplitudes.Transform(
                false,
                [&](Detail::MutableStateView &state, std::complex<double> &amplitude) noexcept {
                    if (state.get(qubit) != result)
                        return false;
                    amplitude *= normFactor;
                    return true;
                },
                true);
        intermediateAmplitudes.CompactIfSparse();
        scratchAmplitudes.ReleaseIfOversized(intermediateAmplitudes.size(), intermediateAmplitudes.QubitCount());

        return result;
    }

  private:
    // Measurement is conditional on the retained state, including after
    // pruning. Reject undefined distributions before consuming randomness.
    double MeasurementTarget(double total)
    {
        if (!(total > 0.) || !std::isfinite(total))
            throw std::domain_error("Path integral measurement requires finite positive probability mass");
        const double target = uniformZeroOne(rng) * total;
        return target < total ? target : std::nextafter(total, 0.);
    }

    // TODO: this can be parallelized!
    std::complex<double> Propagate(const FastVectorBool &currentState, const FastVectorBool &endState, size_t gateIndex,
                                   size_t possibleQubitsChanges)
    {
        if (gateIndex >= circuit.size())
            return currentState == endState ? std::complex<double>(1., 0.) : std::complex<double>(0., 0.);

        const auto &gate = circuit[gateIndex];
        const auto &U = gate.getRawOperatorMatrix();
        const size_t gateQubits = gate.getQubitsNumber();

        assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

        possibleQubitsChanges -= gateQubits;

        std::complex<double> amplitude(0., 0.);

        if (gateQubits == 1)
        {
            const size_t qubit = gate.getQubit1();
            assert(qubit < currentState.size());

            const Eigen::Index col = (currentState.get(qubit) ? 1 : 0);

            for (Eigen::Index row = 0; row < 2; ++row)
            {
                const std::complex<double> val = U(row, col);
                if (std::norm(val) > epsilon)
                {
                    auto nextState = currentState;
                    nextState.set(qubit, row == 1);

                    if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
                        continue;

                    const auto localAmplitude =
                        val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
                    amplitude += localAmplitude;
                }
            }
        }
        else if (gateQubits == 2)
        {
            const size_t qubit1 = gate.getQubit1();
            const size_t qubit2 = gate.getQubit2();
            assert(qubit1 < currentState.size() && qubit2 < currentState.size());

            const Eigen::Index col = ((currentState.get(qubit2) ? 2 : 0) | (currentState.get(qubit1) ? 1 : 0));

            for (Eigen::Index row = 0; row < 4; ++row)
            {
                const std::complex<double> val = U(row, col);
                if (std::norm(val) > epsilon)
                {
                    auto nextState = currentState;
                    nextState.set(qubit1, (row & 1) == 1);
                    nextState.set(qubit2, (row & 2) == 2);

                    if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
                        continue;

                    const auto localAmplitude =
                        val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
                    amplitude += localAmplitude;
                }
            }
        }
        else // only up to three qubits gates are supported
        {
            const size_t qubit1 = gate.getQubit1();
            const size_t qubit2 = gate.getQubit2();
            const size_t qubit3 = gate.getQubit3();
            assert(qubit1 < currentState.size() && qubit2 < currentState.size() && qubit3 < currentState.size());

            const Eigen::Index col = ((currentState.get(qubit3) ? 4 : 0) | (currentState.get(qubit2) ? 2 : 0) |
                                      (currentState.get(qubit1) ? 1 : 0));

            for (Eigen::Index row = 0; row < 8; ++row)
            {
                const std::complex<double> val = U(row, col);
                if (std::norm(val) > epsilon)
                {
                    auto nextState = currentState;
                    nextState.set(qubit1, (row & 1) == 1);
                    nextState.set(qubit2, (row & 2) == 2);
                    nextState.set(qubit3, (row & 4) == 4);

                    if (CountDifferentQubits(nextState, endState) > possibleQubitsChanges)
                        continue;

                    const auto localAmplitude =
                        val * Propagate(nextState, endState, gateIndex + 1, possibleQubitsChanges);
                    amplitude += localAmplitude;
                }
            }
        }

        return amplitude;
    }

    void PropagateBackward(const FastVectorBool &endState, std::complex<double> amplitude)
    {
        intermediateAmplitudes.clear();
        intermediateAmplitudes[endState] = amplitude;

        for (const auto &gate : circuitBack)
            PropagateStep(gate, intermediateAmplitudes);
    }

    AmplitudeMap PropagateAll(const FastVectorBool &startState)
    {
        AmplitudeMap currentAmplitudes;
        currentAmplitudes[startState] = std::complex<double>(1., 0.);

        for (const auto &gate : circuit)
            PropagateStep(gate, currentAmplitudes);

        return currentAmplitudes;
    }

    std::complex<double> PropagateForward(const FastVectorBool &startState)
    {
        const auto currentAmplitudes = PropagateAll(startState);

        std::complex<double> result(0., 0.);
        for (const auto &[state, amp] : currentAmplitudes)
        {
            auto it = intermediateAmplitudes.find(state);
            if (it != intermediateAmplitudes.end())
                result += amp * it->second;
        }

        return result;
    }

    // TODO: this could be improved, this sums the max possible considering only the number of qubits that the gates act
    // on, not the action of the particular gates
    static size_t CountPossibleQubitsChanges(const std::vector<QC::Gates::AppliedGate<>> &circuit)
    {
        size_t count = 0;

        for (const auto &gate : circuit)
            count += gate.getQubitsNumber();

        return count;
    }

    static size_t CountDifferentQubits(const FastVectorBool &curState, const FastVectorBool &endState)
    {
        size_t count = 0;
        const size_t words = curState.nWords();
        for (size_t i = 0; i < words; ++i)
        {
            const uint64_t diff = curState.getWords()[i] ^ endState.getWords()[i];
#ifdef _MSC_VER
            count += __popcnt64(diff);
#else
            count += __builtin_popcountll(diff);
#endif
        }
        return count;
    }

    std::vector<QC::Gates::AppliedGate<>> circuit;
    std::vector<QC::Gates::AppliedGate<>> circuitBack;
    AmplitudeMap intermediateAmplitudes;

    AmplitudeMap savedAmplitudes;
    AmplitudeMap scratchAmplitudes;

    size_t doublingsLimit;
    double epsilon;

    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniformZeroOne;
};
} // namespace PathIntegral
} // namespace QC
