#pragma once

#include "PathIntegralRuntimeRegressionTests.h"

namespace QC
{
namespace PathIntegralColdRegression
{
inline void MalformedMatrices()
{
    using namespace PathIntegralLightconeRegression;
    size_t checks = 0;
    const std::vector<Gate> valid{{MatrixFor(0, 0.), 0}, {MatrixFor(2, .3), 0}};
    for (size_t split : {size_t(0), std::numeric_limits<size_t>::max()})
        for (size_t cache : {size_t(0), size_t(16 * 1024 * 1024)})
            for (const auto shape : {std::pair<int, int>{0, 0},
                                     {0, 2},
                                     {1, 1},
                                     {2, 0},
                                     {2, 1},
                                     {2, 3},
                                     {3, 3},
                                     {4, 3},
                                     {4, 5},
                                     {5, 5},
                                     {6, 6},
                                     {7, 8},
                                     {8, 7},
                                     {8, 9},
                                     {9, 9},
                                     {16, 16}})
            {
                Simulator sim;
                sim.SetMaxDoublingsForBackwardPaths(split);
                sim.SetCompiledGateCacheMaxBytes(cache);
                sim.PropagateAll(valid, std::vector<bool>(3));
                sim.SaveAmplitudes();
                sim.SetSeed(78321);
                sim.MeasureNoCollapse();
                const auto before = sim.GetAmplitudes();
                auto control = sim;
                const auto forward = sim.GetForwardGateCount();
                const Matrix matrix = Matrix::Identity(shape.first, shape.second);
                // Include repeated physical targets, which force forward-only evaluation.
                const std::vector<Gate> invalid{{matrix, 0, 0, 0}};
                const auto rejected = [&](auto action) {
                    bool caught = false;
                    try
                    {
                        action();
                    }
                    catch (const std::invalid_argument &)
                    {
                        caught = true;
                    }
                    Require(caught, "Malformed path matrix was not rejected upfront");
                    ++checks;
                };
                rejected([&] { sim.SetCircuit(invalid); });
                rejected([&] { sim.PropagateAll(invalid, std::vector<bool>(3)); });
                rejected([&] { sim.PropagateAll(invalid); });
                auto external = before;
                rejected([&] { sim.PropagateStep(invalid.front(), external); });
                Require(external.size() == before.size(), "Rejected incremental matrix changed the state");
                Require(sim.GetForwardGateCount() == forward, "Rejected matrix changed the recorded circuit");
                Require(sim.GetAmplitudes().size() == before.size(), "Rejected matrix changed state size");
                for (const auto &entry : before)
                    Close(sim.GetAmplitudes().at(entry.first), entry.second, checks);
                sim.RestoreAmplitudes();
                for (const auto &entry : before)
                    Close(sim.GetAmplitudes().at(entry.first), entry.second, checks);
                for (size_t draw = 0; draw < 4; ++draw)
                    Require(sim.MeasureNoCollapse() == control.MeasureNoCollapse(),
                            "Rejected matrix changed RNG state");
                for (size_t query = 0; query < 2; ++query)
                    Close(sim.Propagate(std::vector<bool>(3)), control.Propagate(std::vector<bool>(3)), checks);
            }
    // Exact supported dimensions remain valid, including repeated targets.
    for (int dimension : {2, 4, 8})
    {
        Simulator sim;
        sim.SetCircuit({Gate(Matrix::Identity(dimension, dimension), 0, 0, 0)});
        for (size_t query = 0; query < 3; ++query)
            Close(sim.Propagate(std::vector<bool>(3)), 1., checks);
    }
    std::cout << "PASS " << checks
              << " malformed dimensions, transactional rejection and valid repeated-target comparisons\n";
}

inline void ProbabilityBounds()
{
    using namespace PathIntegralLightconeRegression;
    size_t checks = 0;
    const std::vector<Gate> circuit{{MatrixFor(0, 0.), 0}};
    Simulator empty;
    Close(empty.QubitProbability(std::numeric_limits<size_t>::max()), 0., checks);
    for (size_t width : {size_t(1), size_t(3), size_t(63), size_t(64), size_t(65), size_t(129), size_t(1024)})
    {
        Simulator sim;
        sim.PropagateAll(circuit, std::vector<bool>(width));
        sim.SetSeed(397);
        auto control = sim;
        for (bool value : {false, true})
        {
            Close(sim.QubitProbability(0, value), .5, checks);
            if (width > 1)
                Close(sim.QubitProbability(width - 1, value), value ? 0. : 1., checks);
            for (size_t qubit : {width, std::numeric_limits<size_t>::max()})
            {
                bool caught = false;
                try
                {
                    sim.QubitProbability(qubit, value);
                }
                catch (const std::out_of_range &)
                {
                    caught = true;
                }
                Require(caught, "Out-of-range path probability was not rejected");
                ++checks;
            }
        }
        for (size_t draw = 0; draw < 8; ++draw)
            Require(sim.MeasureNoCollapse() == control.MeasureNoCollapse(), "Probability query changed RNG state");
    }
    Simulator inferred;
    inferred.PropagateAll(circuit);
    Require(inferred.GetAmplitudes().QubitCount() == 1, "One-argument propagation changed inferred width");
    bool caught = false;
    try
    {
        inferred.QubitProbability(1);
    }
    catch (const std::out_of_range &)
    {
        caught = true;
    }
    Require(caught, "Inferred-width probability was not bounds checked");
    std::cout << "PASS " << checks << " probability bounds and explicit/inferred register widths\n";
}

inline void Run()
{
    MalformedMatrices();
    ProbabilityBounds();
    using namespace PathIntegralLightconeRegression;
    size_t checks = 0;
    std::vector<Gate> circuit;
    for (size_t qubit = 0; qubit < 6; ++qubit)
        circuit.emplace_back(MatrixFor(0, 0.), qubit);
    for (size_t index = 0; index < 42; ++index)
        circuit.emplace_back(MatrixFor(2, .13 + .003 * (index % 7)), index % 6);
    for (size_t qubit = 0; qubit < 6; ++qubit)
        circuit.emplace_back(MatrixFor(0, 0.), qubit);
    Simulator sim;
    sim.SetTrimValue(0.);
    sim.SetCircuit(circuit);
    // The same recorded circuit can be queried at different register widths.
    // A matching low-word key must not reuse a map with a different width.
    for (size_t width : {size_t(6), size_t(65), size_t(129), size_t(1024), size_t(6)})
    {
        for (size_t start : {size_t(0), size_t(3)})
            for (size_t target = 0; target < 8; ++target)
            {
                std::vector<bool> initial(width), endpoint(width);
                for (size_t bit = 0; bit < 6; ++bit)
                {
                    initial[bit] = ((start >> bit) & 1) != 0;
                    endpoint[bit] = ((target >> bit) & 1) != 0;
                }
                Close(sim.Propagate(initial, endpoint),
                      DenseAmplitude(circuit, 6, start, target, std::numeric_limits<size_t>::max(), 0.), checks);
            }
        Require(sim.GetCachedForwardEntries() == 64, "Cold key test did not exercise forward cache");
        if (width > 6)
        {
            std::vector<bool> initial(width), endpoint(width);
            endpoint[width - 1] = true;
            Close(sim.Propagate(initial, endpoint), 0., checks);
            endpoint[width - 1] = false;
            Close(sim.Propagate(initial, endpoint),
                  DenseAmplitude(circuit, 6, 0, 0, std::numeric_limits<size_t>::max(), 0.), checks);
        }
    }

    // Endpoint rejection in the unsplit evaluator must not discard an existing
    // forward state, including a state exposed through GetAmplitudes().
    Simulator unsplit;
    const std::vector<bool> initial(3);
    const std::vector<Gate> oneGate{{MatrixFor(0, 0.), 0}};
    unsplit.PropagateAll(oneGate, initial);
    const auto stored = unsplit.GetAmplitudes();
    auto impossible = initial;
    impossible[2] = true;
    Close(unsplit.Propagate(initial, impossible), 0., checks);
    Require(unsplit.GetAmplitudes().size() == stored.size(), "Cold unsplit query discarded a forward state");
    for (const auto &entry : stored)
        Close(unsplit.GetAmplitudes().at(entry.first), entry.second, checks);
    std::cout << "PASS " << checks << " lazy cache keys, changing widths, idle bits and unsplit state comparisons\n";
}
} // namespace PathIntegralColdRegression
} // namespace QC
