#pragma once

#include "PathIntegralExpectation.h"
#include "PathIntegralGate.h"
#include "PathIntegralStorage.h"
#include <array>
#include <cmath>
#include <memory>
#include <new>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include "QuantumGate.h"

#define QC_PATH_INTEGRAL_COMPILED_TRANSITIONS 1
#define QC_PATH_INTEGRAL_GROWTH_SPLIT 1
#define QC_PATH_INTEGRAL_QUERY_CACHE 1

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

    PathIntegralSimulator(const PathIntegralSimulator &other)
        : circuit(other.circuit), circuitBack(other.circuitBack), forwardTouches(other.forwardTouches),
          backwardTouches(other.backwardTouches), intermediateAmplitudes(other.intermediateAmplitudes),
          savedAmplitudes(other.savedAmplitudes), compiledTables(other.compiledTables),
          compilationAttempted(other.compilationAttempted), recordedQueryCount(other.recordedQueryCount),
          compiledGateCacheMaxBytes(other.compiledGateCacheMaxBytes), useGrowthAwareSplit(other.useGrowthAwareSplit),
          forwardCache(other.forwardCache), forwardCacheMaxBytes(other.forwardCacheMaxBytes),
          doublingsLimit(other.doublingsLimit), epsilon(other.epsilon), rng(other.rng),
          uniformZeroOne(other.uniformZeroOne)
    {
        // Scratch holds a prior frontier after a swap. It is disposable work
        // storage; normal copies still preserve logical state and RNG position.
    }

    PathIntegralSimulator(PathIntegralSimulator &&) noexcept = default;
    PathIntegralSimulator &operator=(PathIntegralSimulator &&) noexcept = default;

    PathIntegralSimulator &operator=(const PathIntegralSimulator &other)
    {
        if (this != &other)
        {
            // Circuit, schedules, and caches form one coherent revision. A
            // failed allocation while copying must not mix two revisions.
            PathIntegralSimulator next(other);
            *this = std::move(next);
        }
        return *this;
    }

    void SetSeed(uint64_t theSeed)
    {
        std::seed_seq seed{uint32_t(theSeed & 0xffffffff), uint32_t(theSeed >> 32)};
        rng.seed(seed);
    }

    void SetTrimValue(double val)
    {
        if (epsilon == val)
            return;
        epsilon = val;
        InvalidateCompiled();
    }

    double GetTrimValue() const
    {
        return epsilon;
    }

    // Recorded gates are immutable until SetCircuit. Cap their compiled cache
    // independently of the state maps, especially for very long circuits.
    void SetCompiledGateCacheMaxBytes(size_t bytes)
    {
        if (compiledGateCacheMaxBytes != bytes)
        {
            compiledGateCacheMaxBytes = bytes;
            InvalidateCompiled();
        }
    }

    size_t GetCompiledGateCacheMaxBytes() const
    {
        return compiledGateCacheMaxBytes;
    }

    size_t GetCompiledGateCount() const
    {
        return compiledTables ? compiledTables->forward.size() + compiledTables->backward.size() : 0;
    }

    void SetForwardCacheMaxBytes(size_t bytes)
    {
        if (forwardCacheMaxBytes != bytes)
        {
            forwardCacheMaxBytes = bytes;
            forwardCache.reset();
        }
    }

    size_t GetForwardCacheMaxBytes() const
    {
        return forwardCacheMaxBytes;
    }

    size_t GetCachedForwardEntries() const
    {
        return forwardCache ? forwardCache->amplitudes.size() : 0;
    }

    void SetMaxDoublingsForBackwardPaths(size_t doublings)
    {
        if (doublings == doublingsLimit)
            return;
        // Changing the cut after SetCircuit must retain both recorded halves.
        // Only reconstruct on a setting change; do not store another copy of
        // every matrix merely to support this uncommon operation.
        if (!circuit.empty() || !circuitBack.empty())
        {
            auto fullCircuit = circuit;
            fullCircuit.reserve(circuit.size() + circuitBack.size());
            for (auto it = circuitBack.rbegin(); it != circuitBack.rend(); ++it)
                fullCircuit.emplace_back(it->getRawOperatorMatrix().transpose(), it->getQubit1(), it->getQubit2(),
                                         it->getQubit3());
            PrepareCircuit(fullCircuit, doublings, false, useGrowthAwareSplit);
        }
        doublingsLimit = doublings;
    }

    size_t GetMaxDoublingsForBackwardPaths() const
    {
        return doublingsLimit;
    }

    // Changing the cut changes threshold checkpoints in approximate runs.
    // Keep this heuristic explicit so the historical cut remains the default.
    void SetUseGrowthAwareSplit(bool enabled)
    {
        if (enabled == useGrowthAwareSplit)
            return;
        if (!circuit.empty() || !circuitBack.empty())
        {
            auto fullCircuit = circuit;
            fullCircuit.reserve(circuit.size() + circuitBack.size());
            for (auto it = circuitBack.rbegin(); it != circuitBack.rend(); ++it)
                fullCircuit.emplace_back(it->getRawOperatorMatrix().transpose(), it->getQubit1(), it->getQubit2(),
                                         it->getQubit3());
            PrepareCircuit(fullCircuit, doublingsLimit, false, enabled);
        }
        useGrowthAwareSplit = enabled;
    }

    bool GetUseGrowthAwareSplit() const
    {
        return useGrowthAwareSplit;
    }

    size_t GetForwardGateCount() const
    {
        return circuit.size();
    }

    AmplitudeMap &GetAmplitudes()
    {
        return intermediateAmplitudes;
    }

    double ExpectationValue(const std::string &pauli) const
    {
        return Detail::PauliExpectation(intermediateAmplitudes, epsilon, pauli);
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
        forwardTouches = TouchSchedule{};
        backwardTouches = TouchSchedule{};
        InvalidateCompiled();
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
        theClone.forwardTouches = forwardTouches;
        theClone.backwardTouches = backwardTouches;
        theClone.compiledGateCacheMaxBytes = compiledGateCacheMaxBytes;
        theClone.compiledTables = compiledTables;
        theClone.compilationAttempted = compilationAttempted;
        theClone.recordedQueryCount = recordedQueryCount;
        theClone.useGrowthAwareSplit = useGrowthAwareSplit;
        theClone.forwardCacheMaxBytes = forwardCacheMaxBytes;
        theClone.forwardCache = forwardCache;

        return theClone;
    }

    void SetCircuit(const std::vector<QC::Gates::AppliedGate<>> &circuit)
    {
        PrepareCircuit(circuit, doublingsLimit, true, useGrowthAwareSplit);
    }

    // the start state is |0>
    std::complex<double> Propagate(const std::vector<bool> &endState)
    {
        const std::vector<bool> startState(endState.size(), false);

        return Propagate(startState, endState);
    }

    std::complex<double> Propagate(const std::vector<bool> &startState, const std::vector<bool> &endState)
    {
        if (startState.size() != endState.size())
            throw std::invalid_argument("Path integral endpoints must have the same width");

        const FastVectorBool startBits(startState);
        const FastVectorBool endBits(endState);
        if (std::max(forwardTouches.width, backwardTouches.width) > startState.size())
            throw std::out_of_range("Path integral gate qubit is outside the register");
        for (size_t word = 0; word < startBits.nWords(); ++word)
            if ((startBits.getWords()[word] ^ endBits.getWords()[word]) &
                ~(forwardTouches.touched[word] | backwardTouches.touched[word]))
            {
                // Split queries replace the map with their backward half;
                // unsplit queries never touch it (e.g. a PropagateAll state).
                if (!circuitBack.empty())
                    intermediateAmplitudes.clear();
                return 0.;
            }

        // First use keeps the original stack-compiled gate path. Only repeat
        // queries amortize allocating and retaining a circuit-wide table cache.
        if (recordedQueryCount < 2)
            ++recordedQueryCount;
        if (!circuitBack.empty())
        {
            intermediateAmplitudes.clear();

            PropagateBackward(endBits, startBits);

            return PropagateForward(startBits, endBits);
        }

        if (const auto *compiled = Compiled())
            return PropagateCached(startBits, endBits, 0, compiled->forward);
        return Propagate(startBits, endBits, 0);
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
        PrepareCircuit(circuit, 0);
        recordedQueryCount = 1;

        intermediateAmplitudes = PropagateAll(startBits);
    }

    void PropagateStep(const QC::Gates::AppliedGate<> &gate, AmplitudeMap &currentAmplitudes)
    {
        PropagateStepConstrained(gate, currentAmplitudes, 0, nullptr);
    }

  private:
    struct CompiledTables
    {
        std::vector<Detail::SparseGate> forward, backward;
    };

    struct ForwardCache
    {
        FastVectorBool start;
        std::array<uint64_t, FastVectorBool::MaxWords> fixedEnd{};
        AmplitudeMap amplitudes;
    };

    void InvalidateCompiled()
    {
        compiledTables.reset();
        compilationAttempted = false;
        recordedQueryCount = 0;
        forwardCache.reset();
    }

    const CompiledTables *Compiled()
    {
        if (compilationAttempted)
            return compiledTables.get();
        if (recordedQueryCount < 2)
            return nullptr;
        return BuildCompiled();
    }

#if defined(_MSC_VER)
    __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    __attribute__((noinline))
#endif
    const CompiledTables *
    BuildCompiled()
    {
        const size_t count = circuit.size() + circuitBack.size();
        if (!count || count > compiledGateCacheMaxBytes / sizeof(Detail::SparseGate))
        {
            compilationAttempted = true;
            return nullptr;
        }
        try
        {
            auto next = std::make_shared<CompiledTables>();
            next->forward.reserve(circuit.size());
            next->backward.reserve(circuitBack.size());
            for (const auto &gate : circuit)
                next->forward.emplace_back(gate, epsilon);
            for (const auto &gate : circuitBack)
                next->backward.emplace_back(gate, epsilon);
            compiledTables = std::move(next);
        }
        catch (const std::bad_alloc &)
        {
            // Optional caching must not make an otherwise affordable query
            // fail. Validation errors still propagate; only allocation falls
            // back to the original per-gate path until settings change.
            compilationAttempted = true;
            return nullptr;
        }
        compilationAttempted = true;
        return compiledTables.get();
    }

    // One byte per gate records its final local writes. A qubit repeated in a
    // gate is represented only by its last overwrite, matching SetRow.
    struct TouchSchedule
    {
        std::vector<unsigned char> finishing;
        std::array<uint64_t, FastVectorBool::MaxWords> touched{};
        size_t width = 0;
    };

    static std::array<size_t, 3> GateQubits(const QC::Gates::AppliedGate<> &gate)
    {
        return {gate.getQubit1(), gate.getQubit2(), gate.getQubit3()};
    }

    static TouchSchedule BuildTouches(const std::vector<QC::Gates::AppliedGate<>> &gates)
    {
        TouchSchedule result;
        result.finishing.resize(gates.size());
        for (size_t index = gates.size(); index;)
        {
            const auto &gate = gates[--index];
            const size_t arity = gate.getQubitsNumber();
            if (!arity || arity > 3)
                throw std::invalid_argument("Path integral gates must act on one to three qubits");
            const auto qubits = GateQubits(gate);
            for (size_t local = arity; local;)
            {
                const size_t bit = --local, qubit = qubits[bit];
                if (qubit >= FastVectorBool::MaxWords * 64)
                    throw std::out_of_range("Path integral gate qubit exceeds the state capacity");
                const uint64_t mask = uint64_t{1} << (qubit % 64);
                if (!(result.touched[qubit / 64] & mask))
                    result.finishing[index] |= static_cast<unsigned char>(1u << bit);
                result.touched[qubit / 64] |= mask;
                result.width = std::max(result.width, qubit + 1);
            }
        }
        return result;
    }

    static void LimitFinishing(TouchSchedule &schedule, const std::vector<QC::Gates::AppliedGate<>> &gates,
                               const TouchSchedule &opposite)
    {
        for (size_t index = 0; index < gates.size(); ++index)
        {
            const auto qubits = GateQubits(gates[index]);
            for (size_t bit = 0; bit < gates[index].getQubitsNumber(); ++bit)
                if (opposite.touched[qubits[bit] / 64] & (uint64_t{1} << (qubits[bit] % 64)))
                    schedule.finishing[index] &= static_cast<unsigned char>(~(1u << bit));
        }
    }

    static size_t MaskWeight(const std::array<uint64_t, FastVectorBool::MaxWords> &mask)
    {
        size_t result = 0;
        for (uint64_t word : mask)
#if defined(_MSC_VER) && defined(_M_X64) && defined(__AVX2__)
            result += static_cast<size_t>(__popcnt64(word));
#elif defined(__GNUC__) || defined(__clang__)
            result += static_cast<size_t>(__builtin_popcountll(word));
#else
            for (; word; word &= word - 1)
                ++result;
#endif
        return result;
    }

    // Estimate work rather than merely count branching gates. A varying-bit
    // mask bounds support by active qubits, while fanout bounds correlations.
    // Last touches are removed just as endpoint pruning removes their rows.
    std::vector<double> EstimateWork(const std::vector<QC::Gates::AppliedGate<>> &source, bool backwards) const
    {
        std::array<size_t, FastVectorBool::MaxWords * 64> remaining{};
        for (const auto &gate : source)
        {
            const auto qubits = GateQubits(gate);
            for (size_t local = 0; local < gate.getQubitsNumber(); ++local)
                ++remaining[qubits[local]];
        }
        std::array<uint64_t, FastVectorBool::MaxWords> varying{};
        std::vector<double> work(source.size() + 1);
        double fanoutLog = 0., total = 0.;
        for (size_t step = 0; step < source.size(); ++step)
        {
            const size_t index = backwards ? source.size() - 1 - step : step;
            const auto &gate = source[index];
            const auto qubits = GateQubits(gate);
            const Detail::SparseGate transitions =
                backwards
                    ? Detail::SparseGate(QC::Gates::AppliedGate<>(gate.getRawOperatorMatrix().transpose(),
                                                                  gate.getQubit1(), gate.getQubit2(), gate.getQubit3()),
                                         epsilon)
                    : Detail::SparseGate(gate, epsilon);
            bool active = transitions.maxFanout > 1;
            for (size_t local = 0; local < transitions.arity; ++local)
                active = active || bool(varying[qubits[local] / 64] & (uint64_t{1} << (qubits[local] % 64)));
            if (transitions.maxFanout > 1)
                fanoutLog += std::log2(double(transitions.maxFanout));
            for (size_t local = 0; local < transitions.arity; ++local)
            {
                const size_t qubit = qubits[local];
                const uint64_t bit = uint64_t{1} << (qubit % 64);
                if (active)
                    varying[qubit / 64] |= bit;
                if (--remaining[qubit] == 0)
                    varying[qubit / 64] &= ~bit;
            }
            const double supportLog = std::min(fanoutLog, double(MaskWeight(varying)));
            total += std::exp2(std::min(60., supportLog)) * (1. + transitions.arity);
            work[backwards ? index : index + 1] = total;
        }
        return work;
    }

    size_t GrowthAwareCut(const std::vector<QC::Gates::AppliedGate<>> &source, size_t splitLimit,
                          size_t legacyCut) const
    {
        if (source.empty())
            return 0;
        const auto forwardWork = EstimateWork(source, false), backwardWork = EstimateWork(source, true);
        size_t best = legacyCut;
        double bestWork = forwardWork[best] + backwardWork[best];
        size_t backwardBranches = 0;
        for (size_t cut = source.size(); cut;)
        {
            --cut;
            backwardBranches += source[cut].isBranching() ? 1 : 0;
            if (backwardBranches > splitLimit)
                break;
            const double candidate = forwardWork[cut] + backwardWork[cut];
            if (candidate < bestWork)
            {
                best = cut;
                bestWork = candidate;
            }
        }
        return best;
    }

    void PrepareCircuit(const std::vector<QC::Gates::AppliedGate<>> &source, size_t splitLimit,
                        bool clearAmplitudes = true, bool growthAware = false)
    {
        // The transpose of a local matrix is a backwards step only when its
        // physical targets are distinct. Preserve the supported overwrite
        // behavior of repeated targets by evaluating such circuits forwards.
        for (const auto &gate : source)
        {
            const auto qubits = GateQubits(gate);
            const size_t arity = gate.getQubitsNumber();
            if (!arity || arity > 3)
                throw std::invalid_argument("Path integral gates must act on one to three qubits");
            const auto &matrix = gate.getRawOperatorMatrix();
            const size_t dimension = size_t{1} << arity;
            if (matrix.rows() != dimension || matrix.cols() != dimension)
                throw std::invalid_argument("Invalid path integral gate dimensions");
            for (size_t bit = 0; bit < arity; ++bit)
            {
                if (qubits[bit] >= FastVectorBool::MaxWords * 64)
                    throw std::out_of_range("Path integral gate qubit exceeds the state capacity");
                for (size_t earlier = 0; earlier < bit; ++earlier)
                    if (qubits[earlier] == qubits[bit])
                        splitLimit = 0;
            }
        }
        std::vector<QC::Gates::AppliedGate<>> nextBack;
        size_t gatesAdded = 0;
        if (splitLimit)
        {
            size_t doublings = 0;
            for (const auto &gate : source)
                doublings += gate.isBranching() ? 1 : 0;
            const size_t halfDoublings = doublings / 2;
            doublings = 0;
            for (auto it = source.rbegin(); it != source.rend(); ++it)
            {
                if (it->isBranching() && (++doublings > splitLimit || doublings > halfDoublings))
                    break;
                ++gatesAdded;
            }
            if (growthAware)
                gatesAdded = source.size() - GrowthAwareCut(source, splitLimit, source.size() - gatesAdded);
            for (size_t index = source.size(); index > source.size() - gatesAdded;)
            {
                const auto &gate = source[--index];
                nextBack.emplace_back(gate.getRawOperatorMatrix().transpose(), gate.getQubit1(), gate.getQubit2(),
                                      gate.getQubit3());
            }
        }
        std::vector<QC::Gates::AppliedGate<>> nextCircuit(source.begin(), source.end() - gatesAdded);
        auto nextForwardTouches = BuildTouches(nextCircuit);
        auto nextBackwardTouches = BuildTouches(nextBack);
        LimitFinishing(nextForwardTouches, nextCircuit, nextBackwardTouches);
        LimitFinishing(nextBackwardTouches, nextBack, nextForwardTouches);
        // Commit only after all allocation and validation has succeeded.
        circuit.swap(nextCircuit);
        circuitBack.swap(nextBack);
        forwardTouches = std::move(nextForwardTouches);
        backwardTouches = std::move(nextBackwardTouches);
        InvalidateCompiled();
        if (clearAmplitudes)
            intermediateAmplitudes.clear();
    }

    static void ConstrainRows(Detail::SparseGate &gate, unsigned mask, const FastVectorBool &endpoint)
    {
        unsigned value = 0;
        for (unsigned bit = 0; bit < gate.arity; ++bit)
            if ((mask & (1u << bit)) && endpoint.get(gate.qubits[bit]))
                value |= 1u << bit;
        gate.maxFanout = 0;
        for (unsigned column = 0; column < (1u << gate.arity); ++column)
        {
            unsigned retained = 0;
            for (unsigned transition = 0; transition < gate.counts[column]; ++transition)
                if ((gate.rows[column][transition] & mask) == value)
                {
                    gate.rows[column][retained] = gate.rows[column][transition];
                    gate.values[column][retained++] = gate.values[column][transition];
                }
            gate.counts[column] = retained;
            gate.maxFanout = std::max(gate.maxFanout, retained);
        }
    }

    void PropagateStepConstrained(const QC::Gates::AppliedGate<> &gate, AmplitudeMap &currentAmplitudes,
                                  unsigned fixedRows, const FastVectorBool *endpoint)
    {
        if (currentAmplitudes.empty())
        {
            scratchAmplitudes.ReleaseIfOversized(0, currentAmplitudes.QubitCount());
            return;
        }
        Detail::SparseGate transitions(gate, epsilon);
        currentAmplitudes.CompactIfSparse();
        for (unsigned i = 0; i < transitions.arity; ++i)
            if (transitions.qubits[i] >= currentAmplitudes.QubitCount())
                throw std::out_of_range("Path integral gate qubit is outside the register");
        // A freshly compiled table is already query-local: constrain it in
        // place rather than copying the immutable-cache representation.
        if (fixedRows)
            ConstrainRows(transitions, fixedRows, *endpoint);
        PropagateTransitions(transitions, currentAmplitudes);
    }

    void PropagateCompiled(const Detail::SparseGate &transitions, AmplitudeMap &currentAmplitudes, unsigned fixedRows,
                           const FastVectorBool *endpoint)
    {
        if (currentAmplitudes.empty())
            return;
        currentAmplitudes.CompactIfSparse();
        for (unsigned i = 0; i < transitions.arity; ++i)
            if (transitions.qubits[i] >= currentAmplitudes.QubitCount())
                throw std::out_of_range("Path integral gate qubit is outside the register");
        if (fixedRows)
        {
            auto constrained = transitions;
            ConstrainRows(constrained, fixedRows, *endpoint);
            PropagateTransitions(constrained, currentAmplitudes);
        }
        else
            PropagateTransitions(transitions, currentAmplitudes);
    }

    void PropagateTransitions(const Detail::SparseGate &transitions, AmplitudeMap &currentAmplitudes)
    {
        if (!transitions.maxFanout)
        {
            currentAmplitudes.Transform(
                false, [](Detail::MutableStateView &, std::complex<double> &) noexcept { return false; });
            scratchAmplitudes.ReleaseIfOversized(0, currentAmplitudes.QubitCount());
            return;
        }
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
        if (!intermediateAmplitudes.empty() && qubit >= intermediateAmplitudes.QubitCount())
            throw std::out_of_range("Path integral probability qubit is outside the register");
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
    template <unsigned Arity>
    std::complex<double> PropagateCachedRows(const FastVectorBool &currentState, const FastVectorBool &endState,
                                             size_t gateIndex, const std::vector<Detail::SparseGate> &gates,
                                             unsigned fixedRows, unsigned requiredRow)
    {
        const auto &gate = gates[gateIndex];
        const unsigned column = gate.Column<Arity>(currentState);
        std::complex<double> amplitude(0., 0.);
        for (unsigned transition = 0; transition < gate.counts[column]; ++transition)
        {
            const unsigned row = gate.rows[column][transition];
            if ((row & fixedRows) != requiredRow)
                continue;
            auto nextState = currentState;
            gate.SetRow<Arity>(nextState, row);
            amplitude += gate.values[column][transition] * PropagateCached(nextState, endState, gateIndex + 1, gates);
        }
        return amplitude;
    }

    std::complex<double> PropagateCached(const FastVectorBool &currentState, const FastVectorBool &endState,
                                         size_t gateIndex, const std::vector<Detail::SparseGate> &gates)
    {
        if (gateIndex == gates.size())
            return currentState == endState ? 1. : 0.;
        const auto &gate = gates[gateIndex];
        const unsigned fixedRows = forwardTouches.finishing[gateIndex];
        unsigned requiredRow = 0;
        for (unsigned bit = 0; bit < gate.arity; ++bit)
            if ((fixedRows & (1u << bit)) && endState.get(gate.qubits[bit]))
                requiredRow |= 1u << bit;
        switch (gate.arity)
        {
        case 1:
            return PropagateCachedRows<1>(currentState, endState, gateIndex, gates, fixedRows, requiredRow);
        case 2:
            return PropagateCachedRows<2>(currentState, endState, gateIndex, gates, fixedRows, requiredRow);
        default:
            return PropagateCachedRows<3>(currentState, endState, gateIndex, gates, fixedRows, requiredRow);
        }
    }

    std::complex<double> Propagate(const FastVectorBool &currentState, const FastVectorBool &endState, size_t gateIndex)
    {
        if (gateIndex >= circuit.size())
            return currentState == endState ? std::complex<double>(1., 0.) : std::complex<double>(0., 0.);

        const auto &gate = circuit[gateIndex];
        const auto &U = gate.getRawOperatorMatrix();
        const size_t gateQubits = gate.getQubitsNumber();

        assert(gateQubits > 0 && gateQubits <= 3); // only up to three qubits gates are supported

        // Previously this used a global mismatch count against the sum of
        // remaining gate arities. Only these final writes matter: after this
        // gate, no remaining gate can repair a mismatching finished bit.
        const unsigned fixedRows = forwardTouches.finishing[gateIndex];
        const auto qubits = GateQubits(gate);
        unsigned requiredRow = 0;
        for (size_t bit = 0; bit < gateQubits; ++bit)
            if ((fixedRows & (1u << bit)) && endState.get(qubits[bit]))
                requiredRow |= 1u << bit;

        std::complex<double> amplitude(0., 0.);

        if (gateQubits == 1)
        {
            const size_t qubit = gate.getQubit1();
            assert(qubit < currentState.size());

            const Eigen::Index col = (currentState.get(qubit) ? 1 : 0);

            for (Eigen::Index row = 0; row < 2; ++row)
            {
                if ((static_cast<unsigned>(row) & fixedRows) != requiredRow)
                    continue;
                const std::complex<double> val = U(row, col);
                if (std::norm(val) > epsilon)
                {
                    auto nextState = currentState;
                    nextState.set(qubit, row == 1);

                    const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1);
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
                if ((static_cast<unsigned>(row) & fixedRows) != requiredRow)
                    continue;
                const std::complex<double> val = U(row, col);
                if (std::norm(val) > epsilon)
                {
                    auto nextState = currentState;
                    nextState.set(qubit1, (row & 1) == 1);
                    nextState.set(qubit2, (row & 2) == 2);

                    const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1);
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
                if ((static_cast<unsigned>(row) & fixedRows) != requiredRow)
                    continue;
                const std::complex<double> val = U(row, col);
                if (std::norm(val) > epsilon)
                {
                    auto nextState = currentState;
                    nextState.set(qubit1, (row & 1) == 1);
                    nextState.set(qubit2, (row & 2) == 2);
                    nextState.set(qubit3, (row & 4) == 4);

                    const auto localAmplitude = val * Propagate(nextState, endState, gateIndex + 1);
                    amplitude += localAmplitude;
                }
            }
        }

        return amplitude;
    }

    void PropagateBackward(const FastVectorBool &endState, const FastVectorBool &startState)
    {
        intermediateAmplitudes.clear();
        intermediateAmplitudes[endState] = 1.;

        const auto *compiled = Compiled();
        for (size_t index = 0; index < circuitBack.size() && !intermediateAmplitudes.empty(); ++index)
            if (compiled)
                PropagateCompiled(compiled->backward[index], intermediateAmplitudes, backwardTouches.finishing[index],
                                  &startState);
            else
                PropagateStepConstrained(circuitBack[index], intermediateAmplitudes, backwardTouches.finishing[index],
                                         &startState);
    }

    AmplitudeMap PropagateAll(const FastVectorBool &startState, const FastVectorBool *endState = nullptr)
    {
        AmplitudeMap currentAmplitudes;
        currentAmplitudes[startState] = std::complex<double>(1., 0.);

        const auto *compiled = Compiled();
        for (size_t index = 0; index < circuit.size() && !currentAmplitudes.empty(); ++index)
            if (compiled)
                PropagateCompiled(compiled->forward[index], currentAmplitudes,
                                  endState ? forwardTouches.finishing[index] : 0, endState);
            else
                PropagateStepConstrained(circuit[index], currentAmplitudes,
                                         endState ? forwardTouches.finishing[index] : 0, endState);

        return currentAmplitudes;
    }

    std::complex<double> PropagateForward(const FastVectorBool &startState, const FastVectorBool &endState)
    {
        if (intermediateAmplitudes.empty())
            return 0.;
        std::array<uint64_t, FastVectorBool::MaxWords> fixedEnd{};
        if (forwardCacheMaxBytes && !circuit.empty())
        {
            // Only endpoint bits fixed in this half affect its pruning. The
            // shared qubits may vary across targets without changing this map.
            for (size_t word = 0; word < endState.nWords(); ++word)
                fixedEnd[word] =
                    endState.getWords()[word] & forwardTouches.touched[word] & ~backwardTouches.touched[word];
            if (forwardCache && forwardCache->start == startState && forwardCache->fixedEnd == fixedEnd)
                return JoinForward(forwardCache->amplitudes);
        }
        auto currentAmplitudes = PropagateAll(startState, &endState);
        // Tiny one-shot maps cost less to rebuild than to retain. Bound actual
        // retained capacities, including the hash index, rather than live size.
        if (forwardCacheMaxBytes && !circuit.empty() && (currentAmplitudes.size() >= 64 || circuit.size() >= 64) &&
            FitsForwardCache(currentAmplitudes))
        {
            std::shared_ptr<ForwardCache> next;
            try
            {
                next = std::make_shared<ForwardCache>();
            }
            catch (const std::bad_alloc &)
            {
                return JoinForward(currentAmplitudes);
            }
            next->start = startState;
            next->fixedEnd = fixedEnd;
            next->amplitudes.swap(currentAmplitudes);
            forwardCache = std::move(next);
            return JoinForward(forwardCache->amplitudes);
        }
        return JoinForward(currentAmplitudes);
    }

    bool FitsForwardCache(const AmplitudeMap &map) const
    {
        if (forwardCacheMaxBytes < sizeof(ForwardCache))
            return false;
        size_t remaining = forwardCacheMaxBytes - sizeof(ForwardCache);
        const auto consume = [&](size_t count, size_t bytes) {
            if (count > remaining / bytes)
                return false;
            remaining -= count * bytes;
            return true;
        };
        return consume(map.values.capacity(), sizeof(std::complex<double>)) &&
               consume(map.keys.capacity(), sizeof(uint64_t)) && consume(map.buckets.capacity(), sizeof(size_t));
    }

    std::complex<double> JoinForward(const AmplitudeMap &currentAmplitudes) const
    {
        std::complex<double> result(0., 0.);
        const auto &smaller =
            currentAmplitudes.size() <= intermediateAmplitudes.size() ? currentAmplitudes : intermediateAmplitudes;
        const auto &larger =
            currentAmplitudes.size() <= intermediateAmplitudes.size() ? intermediateAmplitudes : currentAmplitudes;
        for (const auto &[state, amp] : smaller)
        {
            auto it = larger.find(state);
            if (it != larger.end())
                result += amp * it->second;
        }

        return result;
    }

    std::vector<QC::Gates::AppliedGate<>> circuit;
    std::vector<QC::Gates::AppliedGate<>> circuitBack;
    TouchSchedule forwardTouches, backwardTouches;
    AmplitudeMap intermediateAmplitudes;

    AmplitudeMap savedAmplitudes;
    AmplitudeMap scratchAmplitudes;
    std::shared_ptr<const CompiledTables> compiledTables;
    bool compilationAttempted = false;
    unsigned char recordedQueryCount = 0;
    size_t compiledGateCacheMaxBytes = 16 * 1024 * 1024;
    bool useGrowthAwareSplit = false;
    std::shared_ptr<const ForwardCache> forwardCache;
    size_t forwardCacheMaxBytes = 16 * 1024 * 1024;

    size_t doublingsLimit;
    double epsilon;

    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniformZeroOne;
};
} // namespace PathIntegral
} // namespace QC
