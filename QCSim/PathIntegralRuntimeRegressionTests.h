#pragma once

#include "PathIntegralLightconeRegressionTests.h"

namespace QC
{
namespace PathIntegralRuntimeRegression
{
using namespace PathIntegralLightconeRegression;

inline void Transitions()
{
    size_t checks = 0;
    {
        Simulator cold;
        std::vector<Gate> circuit{{MatrixFor(0, 0.), 0}, {MatrixFor(0, 0.), 0}};
        cold.SetCircuit(circuit);
        Close(cold.Propagate(Bits(0, 1)), 1., checks);
        Require(cold.GetCompiledGateCount() == 0, "One-shot query allocated the transition cache");
        Close(cold.Propagate(Bits(0, 1)), 1., checks);
        Require(cold.GetCompiledGateCount() == 2, "Repeat query did not compile transition tables");
        auto copy = cold.Clone();
        Require(copy.GetCompiledGateCount() == 2, "Clone discarded a warm transition cache");
        cold.SetTrimValue(.125);
        Require(cold.GetCompiledGateCount() == 0, "Threshold change retained compiled transitions");
        cold.PropagateAll(circuit, Bits(0, 1));
        Require(cold.GetCompiledGateCount() == 0, "One-shot full state compiled a transition cache");
        Close(cold.Propagate(Bits(0, 1)), 1., checks);
        Require(cold.GetCompiledGateCount() == 2, "Full-state reuse did not compile transition tables");
    }
    std::mt19937 random(0x276147);
    for (size_t limit : {size_t(0), size_t(1), std::numeric_limits<size_t>::max()})
        for (double threshold : {0., 1e-15, .0625})
            for (size_t trial = 0; trial < 8; ++trial)
            {
                std::vector<Gate> circuit;
                for (size_t step = 0; step < 12; ++step)
                {
                    const auto first = random() % 4;
                    circuit.emplace_back(MatrixFor(random() % 6, .1 + (random() % 100) * .027), first, (first + 1) % 4,
                                         (first + 2) % 4);
                }
                Simulator cached, uncached;
                cached.SetTrimValue(threshold);
                uncached.SetTrimValue(threshold);
                cached.SetMaxDoublingsForBackwardPaths(limit);
                uncached.SetMaxDoublingsForBackwardPaths(limit);
                uncached.SetCompiledGateCacheMaxBytes(0);
                cached.SetCircuit(circuit);
                uncached.SetCircuit(circuit);
                for (size_t start : {size_t(0), size_t(9)})
                    for (size_t target = 0; target < 16; ++target)
                    {
                        const auto expected = DenseAmplitude(circuit, 4, start, target, limit, threshold);
                        Close(cached.Propagate(Bits(start, 4), Bits(target, 4)), expected, checks);
                        Close(uncached.Propagate(Bits(start, 4), Bits(target, 4)), expected, checks);
                    }
                auto copy = cached.Clone();
                cached.SetTrimValue(.125);
                cached.SetMaxDoublingsForBackwardPaths(0);
                for (size_t target = 0; target < 16; ++target)
                {
                    Close(cached.Propagate(Bits(0, 4), Bits(target, 4)), DenseAmplitude(circuit, 4, 0, target, 0, .125),
                          checks);
                    Close(copy.Propagate(Bits(0, 4), Bits(target, 4)),
                          DenseAmplitude(circuit, 4, 0, target, limit, threshold), checks);
                }
            }
    // Public incremental gates remain mutable and therefore do not use a
    // matrix-identity cache. Replacing the matrix must take effect immediately.
    Simulator incremental;
    std::vector<bool> zero(2);
    incremental.GetAmplitudes()[PathIntegral::FastVectorBool(zero)] = 1.;
    Gate gate(MatrixFor(0, 0.), 0);
    incremental.PropagateStep(gate, incremental.GetAmplitudes());
    gate.setOperator(Gates::PauliXGate<>().getRawOperatorMatrix());
    incremental.PropagateStep(gate, incremental.GetAmplitudes());
    const double scale = 1. / std::sqrt(2.);
    for (const auto &entry : incremental.GetAmplitudes())
        Close(entry.second, scale, checks);
    std::cout << "PASS " << checks << " compiled path transition and invalidation comparisons\n";
}

inline Complex DenseCut(const std::vector<Gate> &circuit, size_t cut, size_t width, size_t start, size_t target,
                        double epsilon)
{
    if (cut == circuit.size())
        return DenseAmplitude(circuit, width, start, target, 0, epsilon);
    std::vector<Complex> forward(size_t{1} << width), backward(forward.size());
    forward[start] = backward[target] = 1.;
    for (size_t gate = 0; gate < cut; ++gate)
        forward = DenseStep(forward, circuit[gate], epsilon, true);
    for (size_t gate = circuit.size(); gate > cut;)
    {
        const auto &source = circuit[--gate];
        Gate reverse(source.getRawOperatorMatrix().transpose(), source.getQubit1(), source.getQubit2(),
                     source.getQubit3());
        backward = DenseStep(backward, reverse, epsilon, true);
    }
    Complex result = 0.;
    for (size_t state = 0; state < forward.size(); ++state)
        result += forward[state] * backward[state];
    return result;
}

inline void UnusualThresholds()
{
    size_t checks = 0;
    Matrix contraction(2, 2);
    contraction << .25, .5, 0., .75;
    Matrix zero = Matrix::Zero(2, 2);
    for (bool killed : {false, true})
        for (size_t limit : {size_t(0), size_t(1), std::numeric_limits<size_t>::max()})
            for (double threshold : {-1., 0., .125, 1., 2., std::numeric_limits<double>::infinity(),
                                     std::numeric_limits<double>::quiet_NaN()})
            {
                std::vector<Gate> circuit{{MatrixFor(0, 0.), 0},
                                          {contraction, 1},
                                          {contraction, 0},
                                          {MatrixFor(3, 0.), 1, 0},
                                          {MatrixFor(0, 0.), 1}};
                if (killed)
                    circuit.emplace_back(zero, 0);
                Simulator cached, uncached;
                cached.SetTrimValue(threshold);
                uncached.SetTrimValue(threshold);
                cached.SetMaxDoublingsForBackwardPaths(limit);
                uncached.SetMaxDoublingsForBackwardPaths(limit);
                uncached.SetCompiledGateCacheMaxBytes(0);
                cached.SetCircuit(circuit);
                uncached.SetCircuit(circuit);
                for (size_t target = 0; target < 4; ++target)
                {
                    const auto expected = DenseAmplitude(circuit, 2, 0, target, limit, threshold);
                    // Repeat to exercise the warm transition cache as well as
                    // the intentionally uncached first query.
                    Close(cached.Propagate(Bits(0, 2), Bits(target, 2)), expected, checks);
                    Close(cached.Propagate(Bits(0, 2), Bits(target, 2)), expected, checks);
                    Close(uncached.Propagate(Bits(0, 2), Bits(target, 2)), expected, checks);
                }
            }
    std::cout << "PASS " << checks << " zero/nonunitary and unusual-threshold comparisons\n";
}

inline void Cuts()
{
    size_t checks = 0;
    std::mt19937 random(817723);
    for (size_t trial = 0; trial < 12; ++trial)
    {
        std::vector<Gate> circuit;
        for (size_t step = 0; step < 10; ++step)
        {
            const auto first = random() % 5;
            circuit.emplace_back(MatrixFor(random() % 6, .2 + .01 * (random() % 100)), first, (first + 1) % 5,
                                 (first + 2) % 5);
        }
        for (double epsilon : {0., .0625})
            for (size_t limit : {size_t(0), size_t(1), size_t(3), std::numeric_limits<size_t>::max()})
            {
                Simulator adaptive;
                adaptive.SetTrimValue(epsilon);
                adaptive.SetMaxDoublingsForBackwardPaths(limit);
                adaptive.SetUseGrowthAwareSplit(true);
                adaptive.SetCircuit(circuit);
                size_t backwardsBranches = 0;
                for (size_t index = adaptive.GetForwardGateCount(); index < circuit.size(); ++index)
                    backwardsBranches += circuit[index].isBranching() ? 1 : 0;
                Require(backwardsBranches <= limit, "Growth-aware cut exceeded the backward limit");
                for (size_t start : {size_t(0), size_t(19)})
                    for (size_t target = 0; target < 32; ++target)
                        Close(adaptive.Propagate(Bits(start, 5), Bits(target, 5)),
                              DenseCut(circuit, adaptive.GetForwardGateCount(), 5, start, target, epsilon), checks);
                auto copy = adaptive.Clone();
                Require(copy.GetUseGrowthAwareSplit(), "Growth-aware split setting was lost in clone");
                adaptive.SetUseGrowthAwareSplit(false);
                for (size_t target = 0; target < 32; ++target)
                {
                    Close(adaptive.Propagate(Bits(0, 5), Bits(target, 5)),
                          DenseAmplitude(circuit, 5, 0, target, limit, epsilon), checks);
                    Close(copy.Propagate(Bits(0, 5), Bits(target, 5)),
                          DenseCut(circuit, copy.GetForwardGateCount(), 5, 0, target, epsilon), checks);
                }
            }
    }
    std::cout << "PASS " << checks << " growth-aware cut, smaller join and lifecycle comparisons\n";
}

inline void QueryReuse()
{
    size_t checks = 0;
    std::vector<Gate> logical{{MatrixFor(1, .61), 6}};
    for (size_t qubit = 0; qubit < 6; ++qubit)
        logical.emplace_back(MatrixFor(0, 0.), qubit);
    for (size_t index = 0; index < 72; ++index)
        logical.emplace_back(MatrixFor(2, .07 + .003 * (index % 11)), index % 6);
    for (size_t qubit = 0; qubit < 6; ++qubit)
        logical.emplace_back(MatrixFor(0, 0.), qubit);
    for (size_t width : {size_t(7), size_t(129)})
    {
        const auto wideBits = [&](size_t value) {
            std::vector<bool> result(width);
            for (size_t bit = 0; bit < 6; ++bit)
                result[bit] = ((value >> bit) & 1) != 0;
            result[width - 1] = ((value >> 6) & 1) != 0;
            return result;
        };
        std::vector<Gate> circuit;
        for (const auto &gate : logical)
        {
            const auto map = [&](size_t qubit) { return qubit == 6 ? width - 1 : qubit; };
            circuit.emplace_back(gate.getRawOperatorMatrix(), map(gate.getQubit1()), map(gate.getQubit2()),
                                 map(gate.getQubit3()));
        }
        Simulator cached, uncached;
        cached.SetTrimValue(0.);
        uncached.SetTrimValue(0.);
        uncached.SetForwardCacheMaxBytes(0);
        cached.SetCircuit(circuit);
        uncached.SetCircuit(circuit);
        for (size_t start : {size_t(0), size_t(19)})
            for (size_t target = 0; target < 128; ++target)
            {
                auto startBits = wideBits(start), targetBits = wideBits(target);
                const auto expected = DenseAmplitude(logical, 7, start, target, std::numeric_limits<size_t>::max(), 0.);
                const auto actual = cached.Propagate(startBits, targetBits);
                const auto control = uncached.Propagate(startBits, targetBits);
                if (std::abs(actual - expected) >= 3e-11 || std::abs(control - expected) >= 3e-11)
                    std::cerr << "reuse width=" << width << " start=" << start << " target=" << target
                              << " actual=" << actual << " uncached=" << control << " expected=" << expected << '\n';
                Close(actual, expected, checks);
                Close(control, expected, checks);
            }
        Require(cached.GetCachedForwardEntries() == 64, "Forward cache regression never exercised the cache");
        auto clone = cached.Clone();
        auto copy = cached;
        Require(clone.GetCachedForwardEntries() == 64 && copy.GetCachedForwardEntries() == 64,
                "Forward map cache was not retained by copy/clone");
        cached.SetTrimValue(.0625);
        Require(cached.GetCachedForwardEntries() == 0, "Threshold change retained a stale forward cache");
        cached.SetForwardCacheMaxBytes(1);
        for (size_t target = 0; target < 16; ++target)
        {
            const auto startBits = wideBits(0), targetBits = wideBits(target);
            Close(cached.Propagate(startBits, targetBits),
                  DenseAmplitude(logical, 7, 0, target, std::numeric_limits<size_t>::max(), .0625), checks);
            Close(clone.Propagate(startBits, targetBits),
                  DenseAmplitude(logical, 7, 0, target, std::numeric_limits<size_t>::max(), 0.), checks);
        }
        Require(cached.GetCachedForwardEntries() == 0, "Forward cache exceeded its memory budget");
        cached.SetTrimValue(0.);
        cached.SetCircuit({Gate(MatrixFor(1, .37), width - 1)});
        auto end = wideBits(0);
        end[width - 1] = true;
        Close(cached.Propagate(wideBits(0), end), std::sin(.37 / 2), checks);
        clone.Reset();
        Require(clone.GetCachedForwardEntries() == 0, "Reset retained forward cache");
    }
    std::cout << "PASS " << checks << " repeated-target cache, endpoint keys, wide states and budget comparisons\n";
}

inline void Run()
{
    Transitions();
    UnusualThresholds();
    Cuts();
    QueryReuse();
}
} // namespace PathIntegralRuntimeRegression
} // namespace QC
