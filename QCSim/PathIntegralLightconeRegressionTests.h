#pragma once

#include "PathIntegral.h"
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace QC
{
namespace PathIntegralLightconeRegression
{
using Complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Gate = Gates::AppliedGate<>;
using Simulator = PathIntegral::PathIntegralSimulator;

inline void Require(bool value, const char *message)
{
    if (!value)
        throw std::runtime_error(message);
}

inline void Close(Complex actual, Complex expected, size_t &checks)
{
    Require(std::abs(actual - expected) < 3e-11, "Path integral lightcone amplitude mismatch");
    ++checks;
}

inline std::vector<bool> Bits(size_t value, size_t width)
{
    std::vector<bool> result(width);
    for (size_t bit = 0; bit < width; ++bit)
        result[bit] = ((value >> bit) & 1) != 0;
    return result;
}

// Independent dense evolution: enumerate every physical basis state and every
// local matrix row, without using the simulator's sparse kernels or schedules.
inline std::vector<Complex> DenseStep(const std::vector<Complex> &input, const Gate &gate, double epsilon,
                                      bool trimInputs)
{
    std::vector<Complex> output(input.size());
    const std::array<size_t, 3> qubits{gate.getQubit1(), gate.getQubit2(), gate.getQubit3()};
    const size_t arity = gate.getQubitsNumber();
    const auto &matrix = gate.getRawOperatorMatrix();
    for (size_t state = 0; state < input.size(); ++state)
    {
        if (trimInputs && std::norm(input[state]) < epsilon)
            continue;
        unsigned column = 0;
        for (size_t bit = 0; bit < arity; ++bit)
            column |= unsigned((state >> qubits[bit]) & 1) << bit;
        for (unsigned row = 0; row < (1u << arity); ++row)
        {
            if (!(std::norm(matrix(row, column)) > epsilon))
                continue;
            size_t next = state;
            for (size_t bit = 0; bit < arity; ++bit)
            {
                next &= ~(size_t{1} << qubits[bit]);
                next |= size_t((row >> bit) & 1) << qubits[bit];
            }
            output[next] += matrix(row, column) * input[state];
        }
    }
    return output;
}

inline bool RepeatedTargets(const std::vector<Gate> &circuit)
{
    for (const auto &gate : circuit)
    {
        const std::array<size_t, 3> qubits{gate.getQubit1(), gate.getQubit2(), gate.getQubit3()};
        for (size_t bit = 0; bit < gate.getQubitsNumber(); ++bit)
            for (size_t other = 0; other < bit; ++other)
                if (qubits[bit] == qubits[other])
                    return true;
    }
    return false;
}

inline Complex DenseAmplitude(const std::vector<Gate> &circuit, size_t width, size_t start, size_t target, size_t limit,
                              double epsilon)
{
    std::vector<Complex> forward(size_t{1} << width), backward(forward.size());
    forward[start] = backward[target] = 1.;
    size_t cut = circuit.size();
    if (limit && !RepeatedTargets(circuit))
    {
        size_t branches = 0;
        for (const auto &gate : circuit)
            branches += gate.isBranching() ? 1 : 0;
        const size_t half = branches / 2;
        branches = 0;
        while (cut)
        {
            const auto &gate = circuit[cut - 1];
            if (gate.isBranching() && (++branches > limit || branches > half))
                break;
            --cut;
        }
    }
    const bool split = cut != circuit.size();
    for (size_t gate = 0; gate < cut; ++gate)
        forward = DenseStep(forward, circuit[gate], epsilon, split);
    if (!split)
        return forward[target];
    for (size_t gate = circuit.size(); gate > cut;)
    {
        const auto &original = circuit[--gate];
        Gate transpose(original.getRawOperatorMatrix().transpose(), original.getQubit1(), original.getQubit2(),
                       original.getQubit3());
        backward = DenseStep(backward, transpose, epsilon, true);
    }
    Complex amplitude = 0.;
    for (size_t state = 0; state < forward.size(); ++state)
        amplitude += forward[state] * backward[state];
    return amplitude;
}

inline Matrix Fourier(unsigned arity)
{
    const unsigned dimension = 1u << arity;
    Matrix result(dimension, dimension);
    const double pi = std::acos(-1.);
    for (unsigned row = 0; row < dimension; ++row)
        for (unsigned column = 0; column < dimension; ++column)
            result(row, column) = std::polar(1. / std::sqrt(double(dimension)), 2. * pi * row * column / dimension);
    return result;
}

inline Matrix MatrixFor(unsigned kind, double angle)
{
    if (kind == 0)
        return Gates::HadamardGate<>().getRawOperatorMatrix();
    if (kind == 1)
        return Gates::RyGate<>(angle).getRawOperatorMatrix();
    if (kind == 2)
        return Gates::RzGate<>(angle).getRawOperatorMatrix();
    if (kind == 3)
        return Gates::CNOTGate<>().getRawOperatorMatrix();
    if (kind == 4)
        return Gates::ToffoliGate<>().getRawOperatorMatrix();
    return Fourier(2);
}

inline void CheckCircuit(const std::vector<Gate> &circuit, size_t width, size_t &checks)
{
    for (size_t limit : {size_t(0), size_t(1), size_t(2), std::numeric_limits<size_t>::max()})
        for (double epsilon : {0., 1e-15, .0625})
        {
            Simulator sim;
            sim.SetTrimValue(epsilon);
            sim.SetMaxDoublingsForBackwardPaths(limit);
            sim.SetCircuit(circuit);
            for (size_t start : {size_t(0), (size_t{1} << width) - 3})
                for (size_t target = 0; target < (size_t{1} << width); ++target)
                    Close(sim.Propagate(Bits(start, width), Bits(target, width)),
                          DenseAmplitude(circuit, width, start, target, limit, epsilon), checks);
        }
}

inline void Run()
{
    size_t checks = 0;
    std::mt19937 rng(0x1c0e71);
    // Early and late local blocks surround qubits shared across the cut.
    std::vector<Gate> structured{{Fourier(1), 0},         {Fourier(1), 2},        {MatrixFor(3, 0.), 3, 2},
                                 {Fourier(1), 1},         {MatrixFor(1, .73), 2}, {Fourier(1), 4},
                                 {MatrixFor(2, -.41), 3}, {MatrixFor(1, .62), 4}};
    CheckCircuit(structured, 6, checks);
    for (size_t trial = 0; trial < 12; ++trial)
    {
        const size_t width = 3 + rng() % 3;
        std::vector<Gate> circuit;
        for (size_t step = 0; step < 8; ++step)
        {
            const unsigned kind = rng() % 6;
            const size_t first = rng() % width;
            circuit.emplace_back(MatrixFor(kind, .19 + .001 * double(rng() % 1200)), first, (first + 1) % width,
                                 (first + 2) % width);
        }
        CheckCircuit(circuit, width, checks);
    }
    // Local row bits on repeated targets are overwritten in order. Backwards
    // transposition of the unreduced local matrix would give a wrong answer.
    std::vector<Gate> repeated{{Fourier(2), 0, 0}, {MatrixFor(1, .37), 1}, {Fourier(3), 1, 1, 2}};
    CheckCircuit(repeated, 3, checks);

    // Logical three-qubit references remapped across packed word boundaries,
    // with an independently fixed idle bit and arbitrary nonzero start bits.
    for (size_t width : {size_t(63), size_t(64), size_t(65), size_t(129), size_t(1024)})
    {
        const std::array<size_t, 3> active{0, width / 2, width - 1};
        std::vector<Gate> logical{{Fourier(1), 0},
                                  {MatrixFor(1, .57), 1},
                                  {MatrixFor(3, 0.), 2, 1},
                                  {MatrixFor(2, -.37), 0},
                                  {Fourier(1), 2}};
        std::vector<Gate> wide;
        for (const auto &gate : logical)
            wide.emplace_back(gate.getRawOperatorMatrix(), active[gate.getQubit1()], active[gate.getQubit2()],
                              active[gate.getQubit3()]);
        for (size_t limit : {size_t(0), size_t(1), std::numeric_limits<size_t>::max()})
        {
            Simulator sim;
            sim.SetTrimValue(0.);
            sim.SetMaxDoublingsForBackwardPaths(limit);
            sim.SetCircuit(wide);
            for (size_t start : {size_t(0), size_t(5)})
                for (size_t target = 0; target < 8; ++target)
                {
                    std::vector<bool> startBits(width), targetBits(width);
                    startBits[width - 2] = targetBits[width - 2] = true;
                    for (size_t bit = 0; bit < active.size(); ++bit)
                    {
                        startBits[active[bit]] = ((start >> bit) & 1) != 0;
                        targetBits[active[bit]] = ((target >> bit) & 1) != 0;
                    }
                    Close(sim.Propagate(startBits, targetBits), DenseAmplitude(logical, 3, start, target, limit, 0.),
                          checks);
                    targetBits[width - 2] = false;
                    Close(sim.Propagate(startBits, targetBits), 0., checks);
                }
        }
    }

    Simulator lifecycle;
    lifecycle.SetTrimValue(0.);
    lifecycle.SetCircuit(structured);
    const auto start = Bits(13, 6), target = Bits(23, 6);
    const auto expected = DenseAmplitude(structured, 6, 13, 23, 0, 0.);
    auto cloned = lifecycle.Clone();
    auto copied = lifecycle;
    for (size_t limit : {size_t(0), size_t(1), std::numeric_limits<size_t>::max(), size_t(0)})
    {
        lifecycle.SetMaxDoublingsForBackwardPaths(limit);
        Close(lifecycle.Propagate(start, target), expected, checks);
        Close(cloned.Propagate(start, target), expected, checks);
        Close(copied.Propagate(start, target), expected, checks);
    }
    lifecycle.PropagateAll(structured, start);
    const auto full = lifecycle.GetAmplitudes();
    // Unsplit queries leave the stored state alone, including the early exit
    // for endpoints that differ on a qubit (5) that no gate touches.
    Close(lifecycle.Propagate(start, Bits(13 | 32, 6)), 0., checks);
    Require(lifecycle.GetAmplitudes().size() == full.size(), "Out-of-cone query discarded the forward state");
    lifecycle.SaveAmplitudes();
    lifecycle.SetMaxDoublingsForBackwardPaths(1);
    Require(lifecycle.GetAmplitudes().size() == full.size(), "Changing the cut discarded the forward state");
    for (const auto &entry : full)
        Close(lifecycle.GetAmplitudes().at(entry.first), entry.second, checks);
    lifecycle.RestoreAmplitudes();
    for (const auto &entry : full)
        Close(lifecycle.GetAmplitudes().at(entry.first), entry.second, checks);
    // A full-state request still retains amplitudes outside any previous target.
    auto dense = std::vector<Complex>(64);
    dense[13] = 1.;
    for (const auto &gate : structured)
        dense = DenseStep(dense, gate, 0., true);
    for (const auto &entry : lifecycle.GetAmplitudes())
    {
        size_t state = 0;
        for (size_t bit = 0; bit < 6; ++bit)
            state |= size_t(entry.first.get(bit)) << bit;
        Close(entry.second, dense[state], checks);
    }
    lifecycle.SetCircuit({Gate(MatrixFor(1, .91), 5)});
    Close(lifecycle.Propagate(Bits(0, 6), Bits(32, 6)), std::sin(.91 / 2), checks);
    bool rejected = false;
    try
    {
        lifecycle.SetCircuit({Gate(Fourier(1), 1024)});
    }
    catch (const std::out_of_range &)
    {
        rejected = true;
    }
    Require(rejected, "Out-of-capacity gate was not rejected");
    Close(lifecycle.Propagate(Bits(0, 6), Bits(32, 6)), std::sin(.91 / 2), checks);
    lifecycle.Reset();
    Close(lifecycle.Propagate(start, start), 1., checks);
    Close(lifecycle.Propagate(start, target), 0., checks);
    lifecycle.SetCircuit({});
    Close(lifecycle.Propagate({}, {}), 1., checks);
    rejected = false;
    try
    {
        lifecycle.Propagate(Bits(0, 3), Bits(0, 4));
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    Require(rejected, "Mismatched endpoint widths were not rejected");
    std::cout << "PASS " << checks << " path integral lightcone comparisons\n";
}
} // namespace PathIntegralLightconeRegression
} // namespace QC
