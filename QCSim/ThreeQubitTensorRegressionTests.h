#pragma once

#include "MPOSimulator.h"
#include "MPSSimulator.h"
#include "TensorNetworkRegressionTests.h"

#include <array>
#include <iostream>
#include <type_traits>

namespace QC
{
namespace ThreeQubitTensorRegression
{
using namespace TensorNetworkRegression;

// Independent dense embedding: the first argument is the least significant
// matrix bit. In particular this does not use the simulator's routing helper.
inline Matrix Embed(const Matrix &gate, const std::vector<int> &qubits, int count)
{
    const int size = 1 << count;
    Matrix result = Matrix::Zero(size, size);
    for (int col = 0; col < size; ++col)
    {
        int local = 0, outside = col;
        for (size_t bit = 0; bit < qubits.size(); ++bit)
        {
            local |= ((col >> qubits[bit]) & 1) << bit;
            outside &= ~(1 << qubits[bit]);
        }
        for (int row = 0; row < gate.rows(); ++row)
        {
            int global = outside;
            for (size_t bit = 0; bit < qubits.size(); ++bit)
                global |= ((row >> bit) & 1) << qubits[bit];
            result(global, col) = gate(row, local);
        }
    }
    return result;
}

template <bool MPO> void Evolve(Matrix &state, const Matrix &gate, const std::vector<int> &sites, int count)
{
    const Matrix op = Embed(gate, sites, count);
    if constexpr (MPO)
        state = (op * state * op.adjoint()).eval();
    else
        state = (op * state).eval();
}

template <bool MPO, class Simulator> Matrix Storage(const Simulator &simulator)
{
    if constexpr (MPO)
        return simulator.getUnnormalizedDensityMatrix();
    else
        return simulator.getRegisterStorage();
}

template <class Simulator> class LegacyMPSDispatch : public Simulator
{
  public:
    using Simulator::Simulator;
    using GateClass = typename Simulator::GateClass;
    using Index = typename Simulator::IndexType;
    int localCalls = 0, threeCalls = 0;

    void ApplyGate(const GateClass &gate, Index target, Index control = 0) override
    {
        ++localCalls;
        Simulator::ApplyGate(gate, target, control);
    }

    void ApplyGate(const GateClass &gate, Index target, Index control1, Index control2) override
    {
        ++threeCalls;
        Simulator::ApplyGate(gate, target, control1, control2);
    }
};

template <class Simulator> class LegacyMPODispatch : public Simulator
{
  public:
    using Simulator::Simulator;
    using GateClass = typename Simulator::GateClass;
    using Index = typename Simulator::IndexType;
    int localCalls = 0, threeCalls = 0;

    // Before three-site support, MPO's AppliedGate entry point dispatched to
    // this overload; preserve it for users that subclass the wrapper or impl.
    void ApplyOperator(const Gate &gate) override
    {
        ++localCalls;
        Simulator::ApplyOperator(gate);
    }

    void ApplyGate(const GateClass &gate, Index target, Index control1, Index control2) override
    {
        ++threeCalls;
        Simulator::ApplyGate(gate, target, control1, control2);
    }
};

template <bool MPO, class Simulator> void CheckLegacyDispatch()
{
    using Probe = std::conditional_t<MPO, LegacyMPODispatch<Simulator>, LegacyMPSDispatch<Simulator>>;
    Probe probe(3);
    probe.SetMultithreading(false);
    Simulator &base = probe; // Exercise virtual dispatch without overload hiding.
    base.ApplyGate(Gate(Gates::PauliXGate<>().getRawOperatorMatrix(), 0));
    base.ApplyGate(Gate(Gates::CNOTGate<>().getRawOperatorMatrix(), 1, 0));
    Require(probe.localCalls == 2 && probe.threeCalls == 0, "AppliedGate bypassed a legacy virtual override");
    base.ApplyGate(Gate(Gates::ToffoliGate<>().getRawOperatorMatrix(), 2, 1, 0));
    Require(probe.localCalls == 2 && probe.threeCalls == 1, "Applied three-site gate used the legacy overload");
    Close(probe.GetProbability(2), 0., "Virtual-dispatch compatibility probe did not apply its gates");
}

template <bool MPO, class Simulator> void CheckWideRejection()
{
    Simulator simulator(5);
    simulator.SetMultithreading(false);
    simulator.ApplyGate(Gates::HadamardGate<>(), 0);
    simulator.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    const Matrix before = Storage<MPO>(simulator);
    // A matrix on more than three qubits stays rejected. Compact gates on more qubits are
    // applied in place by the one-argument entry points (see CheckCompactGates and
    // CheckWideTables); the overloads that take a matrix and sites cannot carry them.
    const std::array<Gate, 1> unsupported{Gate(Matrix::Identity(16, 16), 0, 2, 4)};
    const std::array<Gate, 3> compact{
        Gate::MultiControlled(Gates::PauliXGate<>(), 0, {2, 3, 4}),
        Gate::Diagonal(std::vector<Complex>(16, Complex(0., 1.)), {0, 2, 3, 4}),
        Gate::MultiplexedRotation(Gates::RotationAxis::Y, std::vector<double>(8, 0.3), 0, {2, 3, 4})};
    using Interface =
        std::conditional_t<MPO, TensorNetworks::MPOSimulatorInterface, TensorNetworks::MPSSimulatorInterface>;
    Interface &api = simulator;
    for (const auto &gate : unsupported)
    {
        Throws<std::invalid_argument>([&] { api.ApplyGate(gate); }, "Tensor network accepted a wide AppliedGate");
        Throws<std::invalid_argument>([&] { api.ApplyGate(gate, 0, 2); }, "Tensor network accepted a wide local gate");
        Throws<std::invalid_argument>([&] { api.ApplyGate(gate, 0, 2, 4); },
                                      "Tensor network accepted a wide three-site gate");
        Throws<std::invalid_argument>([&] { api.ApplyGates({gate}); }, "Tensor network accepted a wide gate batch");
        if constexpr (MPO)
        {
            Throws<std::invalid_argument>([&] { api.ApplyOperator(gate); }, "MPO accepted a wide operator");
            Throws<std::invalid_argument>([&] { api.ApplyOperator(gate, 0, 2); }, "MPO accepted a wide local operator");
            Throws<std::invalid_argument>([&] { api.ApplyOperator(gate, 0, 2, 4); },
                                          "MPO accepted a wide three-site operator");
            Throws<std::invalid_argument>([&] { api.ApplyOperators({gate}); }, "MPO accepted a wide operator batch");
            Throws<std::invalid_argument>([&] { api.ApplyKrausOperators({gate}); },
                                          "MPO accepted a wide Kraus operator");
            Throws<std::invalid_argument>([&] { api.ApplyOperatorAndNormalize(gate); },
                                          "MPO accepted wide postselection");
            Throws<std::invalid_argument>([&] { api.ApplyOperatorAndNormalize(gate, 0, 2); },
                                          "MPO accepted wide local postselection");
            Throws<std::invalid_argument>([&] { api.ApplyOperatorAndNormalize(gate, 0, 2, 4); },
                                          "MPO accepted wide three-site postselection");
        }
        Close(Storage<MPO>(simulator), before, "Rejected wide gate changed tensor-network state");
        if constexpr (std::is_same_v<Simulator, TensorNetworks::MPSSimulator> ||
                      std::is_same_v<Simulator, TensorNetworks::MPOSimulator>)
            Require(simulator.getQubitsMap() == std::vector<typename Simulator::IndexType>({0, 1, 2, 3, 4}),
                    "Rejected wide gate changed tensor-network routing");
    }
    for (const Gate &controlled : compact)
    {
        Throws<std::invalid_argument>([&] { api.ApplyGate(controlled, 0, 2); },
                                      "Tensor network accepted a compact gate as a local matrix");
        Throws<std::invalid_argument>([&] { api.ApplyGate(controlled, 0, 2, 4); },
                                      "Tensor network accepted a compact gate as a three-site matrix");
        if constexpr (MPO)
        {
            Throws<std::invalid_argument>([&] { api.ApplyOperator(controlled, 0, 2); },
                                          "MPO accepted a compact gate as a local operator matrix");
            Throws<std::invalid_argument>([&] { api.ApplyKrausOperators({controlled}); },
                                          "MPO accepted a compact Kraus operator");
            Throws<std::invalid_argument>([&] { api.ApplyOperatorAndNormalize(controlled); },
                                          "MPO accepted compact postselection");
        }
    }
    Close(Storage<MPO>(simulator), before, "Rejected compact gate changed tensor-network state");
}

// The full matrix of a compact gate straight from its definition: bit i of the table or
// control index is qubits[i]. This does not use the library's own expansion.
inline Matrix CompactReference(const Gates::StructuredGate &gate, int count)
{
    const int size = 1 << count;
    Matrix result = Matrix::Zero(size, size);
    for (int col = 0; col < size; ++col)
    {
        size_t selector = 0;
        for (size_t bit = 0; bit < gate.qubits.size(); ++bit)
            selector |= size_t((col >> gate.qubits[bit]) & 1) << bit;
        if (gate.kind == Gates::StructuredGate::Kind::Diagonal)
        {
            result(col, col) = gate.diagonal[selector];
            continue;
        }
        Complex block[2][2] = {{1., 0.}, {0., 1.}};
        if (gate.kind == Gates::StructuredGate::Kind::Controlled)
        {
            const size_t wanted = gate.controlState == std::numeric_limits<size_t>::max()
                                      ? (size_t(1) << gate.qubits.size()) - 1
                                      : gate.controlState;
            if (selector == wanted)
                for (int r = 0; r < 2; ++r)
                    for (int c = 0; c < 2; ++c)
                        block[r][c] = gate.matrix[2 * r + c];
        }
        else
        {
            const double c = std::cos(gate.angles[selector] / 2), sn = std::sin(gate.angles[selector] / 2);
            if (gate.axis == Gates::RotationAxis::X)
                block[0][0] = block[1][1] = c, block[0][1] = block[1][0] = Complex(0., -sn);
            else if (gate.axis == Gates::RotationAxis::Y)
                block[0][0] = block[1][1] = c, block[0][1] = -sn, block[1][0] = sn;
            else
                block[0][0] = Complex(c, -sn), block[1][1] = Complex(c, sn);
        }
        const int bit = (col >> gate.target) & 1;
        for (int out = 0; out < 2; ++out)
            result((col & ~(1 << gate.target)) | (out << gate.target), col) = block[out][bit];
    }
    return result;
}

template <bool MPO> void EvolveCompact(Matrix &state, const Gate &gate, int count)
{
    const Matrix op = CompactReference(*gate.getStructuredGate(), count);
    if constexpr (MPO)
        state = (op * state * op.adjoint()).eval();
    else
        state = (op * state).eval();
}

// Diagonal and multiplexed tables on more than three qubits are applied in place, as the chain
// operator of the table: on scattered qubits in any table order, through every entry point, next
// to ordinary gates, with and without a bond limit. A table that couples the chain beyond the
// rank limit is rejected before anything changes, and so is, on MPS, a diagonal that is not
// unitary.
template <bool MPO, class Simulator> void CheckWideTables()
{
    const bool routed = std::is_same_v<Simulator, TensorNetworks::MPSSimulator> ||
                        std::is_same_v<Simulator, TensorNetworks::MPOSimulator>;
    const int count = MPO ? 6 : 8, size = 1 << count;
    std::mt19937 random(MPO ? 613 : 713);
    const auto pick = [&random](int below) { return int(random() % unsigned(below)); };
    std::uniform_real_distribution<double> angle(-3., 3.);
    const auto start = [&](Simulator &simulator, Matrix &reference) {
        simulator.SetMultithreading(false);
        reference = Matrix::Zero(size, MPO ? size : 1);
        reference(0, 0) = 1.;
        for (int layer = 0; layer < 3; ++layer)
            for (int site = layer % 2; site + 1 < count; site += 2)
            {
                const Matrix gate = Unitary(4, 7000 + 31 * layer + site + (MPO ? 7 : 0));
                simulator.ApplyGate(Gate(gate), site + 1, site);
                Evolve<MPO>(reference, gate, {site + 1, site}, count);
            }
        if (routed)
        {
            // scramble the logical-to-physical map: the compact gates take logical qubits
            const Matrix gate = Unitary(4, 7777);
            simulator.ApplyGate(Gate(gate), count - 1, 0);
            Evolve<MPO>(reference, gate, {count - 1, 0}, count);
            simulator.ApplyGate(Gate(gate), 1, count - 2);
            Evolve<MPO>(reference, gate, {1, count - 2}, count);
        }
    };
    // `arity` distinct qubits in a random order
    const auto someQubits = [&](int arity) {
        std::vector<size_t> qubits(count);
        for (int q = 0; q < count; ++q)
            qubits[q] = q;
        std::shuffle(qubits.begin(), qubits.end(), random);
        qubits.resize(arity);
        return qubits;
    };
    // A product of one-qubit phases: rank one.
    const auto productTable = [&](size_t arity) {
        std::vector<Complex> table(size_t(1) << arity, 1.);
        for (size_t bit = 0; bit < arity; ++bit)
        {
            const Complex phase = std::polar(1., angle(random));
            for (size_t i = 0; i < table.size(); ++i)
                if ((i >> bit) & 1)
                    table[i] *= phase;
        }
        return table;
    };
    // exp(-i sum of couplings between table bits that follow each other): rank two between them.
    const auto couplingTable = [&](size_t arity) {
        std::vector<Complex> table(size_t(1) << arity);
        std::vector<double> coupling(arity);
        for (auto &value : coupling)
            value = angle(random);
        for (size_t i = 0; i < table.size(); ++i)
        {
            double phase = 0.;
            for (size_t bit = 0; bit + 1 < arity; ++bit)
                phase += coupling[bit] * ((((i >> bit) ^ (i >> (bit + 1))) & 1) ? -1. : 1.);
            table[i] = std::polar(1., phase);
        }
        return table;
    };
    const auto randomTable = [&](size_t arity) {
        std::vector<Complex> table(size_t(1) << arity);
        for (auto &value : table)
            value = std::polar(1., angle(random));
        return table;
    };
    const auto randomAngles = [&](size_t controls) {
        std::vector<double> angles(size_t(1) << controls);
        for (auto &value : angles)
            value = angle(random);
        return angles;
    };

    for (int trial = 0; trial < 3; ++trial)
    {
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        const auto apply = [&](const Gate &gate) {
            simulator.ApplyGate(gate);
            EvolveCompact<MPO>(reference, gate, count);
        };
        const int widest = MPO ? count : count - 1;
        apply(Gate::Diagonal(productTable(size_t(widest)), someQubits(widest)));
        apply(Gate::Diagonal(couplingTable(5), someQubits(5)));
        // a table without structure on four qubits: rank up to four
        apply(Gate::Diagonal(randomTable(4), someQubits(4)));
        {
            const int site = pick(count - 1);
            const Matrix local = Unitary(4, 7100 + trial);
            simulator.ApplyGate(Gate(local), site, site + 1);
            Evolve<MPO>(reference, local, {site, site + 1}, count);
        }
        for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
        {
            const auto qubits = someQubits(4 + pick(2));
            apply(Gate::MultiplexedRotation(axis, randomAngles(qubits.size() - 1), qubits[0],
                                            std::vector<size_t>(qubits.begin() + 1, qubits.end())));
            if (trial == 1)
                // a read between the gates must leave no stale environments or cached forms
                (void)simulator.GetProbability(pick(count));
        }
        Close(Storage<MPO>(simulator), reference, "Wide diagonal and multiplexed tables on a tensor network mismatch");
    }

    // The named entry points and the adjoint.
    {
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        const auto qubits = someQubits(5);
        const auto table = couplingTable(5);
        simulator.ApplyDiagonalGate(table, qubits);
        EvolveCompact<MPO>(reference, Gate::Diagonal(table, qubits), count);
        const std::vector<size_t> controls(qubits.begin() + 1, qubits.end());
        const auto angles = randomAngles(controls.size());
        simulator.ApplyMultiplexedRX(angles, qubits[0], controls);
        EvolveCompact<MPO>(reference, Gate::MultiplexedRotation(Gates::RotationAxis::X, angles, qubits[0], controls),
                           count);
        simulator.ApplyMultiplexedRY(angles, qubits[2], {qubits[0], qubits[1], qubits[3], qubits[4]});
        EvolveCompact<MPO>(reference,
                           Gate::MultiplexedRotation(Gates::RotationAxis::Y, angles, qubits[2],
                                                     {qubits[0], qubits[1], qubits[3], qubits[4]}),
                           count);
        simulator.ApplyMultiplexedRZ(angles, qubits[4], {qubits[3], qubits[2], qubits[1], qubits[0]});
        EvolveCompact<MPO>(reference,
                           Gate::MultiplexedRotation(Gates::RotationAxis::Z, angles, qubits[4],
                                                     {qubits[3], qubits[2], qubits[1], qubits[0]}),
                           count);
        Close(Storage<MPO>(simulator), reference, "Named wide-table entry points on a tensor network mismatch");
        // The adjoint undoes the gate.
        const Matrix before = Storage<MPO>(simulator);
        const Gate gate = Gate::MultiplexedRotation(Gates::RotationAxis::Y, angles, qubits[1],
                                                    {qubits[0], qubits[2], qubits[3], qubits[4]});
        simulator.ApplyGate(gate);
        simulator.ApplyGate(gate.Adjoint());
        Close(Storage<MPO>(simulator), before, "The adjoint of a wide multiplexer did not undo it");

        // Rejections leave the state as it was.
        const Matrix kept = Storage<MPO>(simulator);
        std::vector<Complex> shrinking = table;
        shrinking[3] = 0.5;
        if constexpr (MPO)
        {
            // rho -> O rho O^dagger also for a diagonal that is not unitary
            const Gate filter = Gate::Diagonal(shrinking, qubits);
            Matrix expected = kept;
            EvolveCompact<true>(expected, filter, count);
            simulator.ApplyGate(filter);
            Close(Storage<MPO>(simulator), expected, "A wide diagonal filter on an MPO mismatch");
        }
        else
        {
            Throws<std::invalid_argument>([&] { simulator.ApplyDiagonalGate(shrinking, qubits); },
                                          "MPS applied a wide diagonal that is not unitary");
            Close(Storage<MPO>(simulator), kept, "A rejected wide diagonal changed the tensor-network state");
        }
    }

    // A table that couples the chain beyond the rank limit is rejected before anything changes.
    {
        const int wide = MPO ? 10 : 14;
        Simulator simulator(wide);
        simulator.SetMultithreading(false);
        simulator.ApplyGate(Gates::HadamardGate<>(), 0);
        simulator.ApplyGate(Gates::CNOTGate<>(), 1, 0);
        const auto bonds = simulator.getBondDimensions();
        const double probability = simulator.GetProbability(1);
        std::vector<size_t> qubits(wide);
        for (int q = 0; q < wide; ++q)
            qubits[q] = q;
        Throws<std::invalid_argument>([&] { simulator.ApplyDiagonalGate(randomTable(size_t(wide)), qubits); },
                                      "Tensor network accepted a table beyond the rank limit");
        Require(simulator.getBondDimensions() == bonds && std::abs(simulator.GetProbability(1) - probability) < 1E-12,
                "A table beyond the rank limit changed the tensor-network state");
        // The same qubits with a table of low rank are fine.
        simulator.ApplyDiagonalGate(productTable(size_t(wide)), qubits);
        Require(simulator.getBondDimensions() == bonds && std::abs(simulator.GetProbability(1) - probability) < 1E-12,
                "A product of one-qubit phases changed the bonds or the probabilities");
    }

    // With a bond limit the sweep after the gate cuts like a two-qubit gate does.
    for (int cap : {2, 5})
    {
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        simulator.setLimitBondDimension(cap);
        simulator.ApplyGate(Gate::Diagonal(couplingTable(5), someQubits(5)));
        const auto qubits = someQubits(5);
        simulator.ApplyGate(Gate::MultiplexedRotation(Gates::RotationAxis::Y, randomAngles(4), qubits[0],
                                                      std::vector<size_t>(qubits.begin() + 1, qubits.end())));
        for (const auto dimension : simulator.getBondDimensions())
            Require(dimension <= cap, "A wide table left a bond above the limit");
        // What the chain holds now evolves exactly once the limit is lifted.
        simulator.dontLimitBondDimension();
        reference = Storage<MPO>(simulator);
        Require(reference.allFinite(), "A truncated wide table left a non-finite state");
        if constexpr (!MPO)
            Require(std::abs(reference.norm() - 1.) < 1E-9, "A truncated wide table left the MPS off its norm");
        const Gate gate = Gate::Diagonal(couplingTable(4), someQubits(4));
        simulator.ApplyGate(gate);
        EvolveCompact<MPO>(reference, gate, count);
        Close(Storage<MPO>(simulator), reference, "Evolution after a truncated wide table mismatch");
    }
}

// Compact gates on tensor networks: a multi-controlled one-qubit gate is applied where its
// qubits are, with any number of controls and any control values; other compact gates on at
// most three qubits go through their matrix.
template <bool MPO, class Simulator> void CheckCompactGates()
{
    const bool routed = std::is_same_v<Simulator, TensorNetworks::MPSSimulator> ||
                        std::is_same_v<Simulator, TensorNetworks::MPOSimulator>;
    const int count = 7, size = 1 << count;
    std::mt19937 random(MPO ? 811 : 911);
    const auto pick = [&random](int below) { return int(random() % unsigned(below)); };
    const auto start = [&](Simulator &simulator, Matrix &reference) {
        simulator.SetMultithreading(false);
        reference = Matrix::Zero(size, MPO ? size : 1);
        reference(0, 0) = 1.;
        for (int layer = 0; layer < 3; ++layer)
            for (int site = layer % 2; site + 1 < count; site += 2)
            {
                const Matrix gate = Unitary(4, 1000 + 31 * layer + site + (MPO ? 7 : 0));
                simulator.ApplyGate(Gate(gate), site + 1, site);
                Evolve<MPO>(reference, gate, {site + 1, site}, count);
            }
        if (routed)
        {
            // scramble the logical-to-physical map: the compact gates take logical qubits
            const Matrix gate = Unitary(4, 1777);
            simulator.ApplyGate(Gate(gate), 6, 0);
            Evolve<MPO>(reference, gate, {6, 0}, count);
            simulator.ApplyGate(Gate(gate), 1, 5);
            Evolve<MPO>(reference, gate, {1, 5}, count);
        }
    };
    const auto controlledGate = [&](int controls, const Matrix &u) {
        std::vector<size_t> qubits(count);
        for (int q = 0; q < count; ++q)
            qubits[q] = q;
        std::shuffle(qubits.begin(), qubits.end(), random);
        const size_t target = qubits[0];
        std::vector<size_t> chosen(qubits.begin() + 1, qubits.begin() + 1 + controls);
        const size_t state = pick(3) == 0 ? std::numeric_limits<size_t>::max() : size_t(pick(1 << controls));
        return Gate::MultiControlled(u, target, chosen, state);
    };

    for (int trial = 0; trial < 4; ++trial)
    {
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        for (int step = 0; step < 10; ++step)
        {
            // 1 to 6 controls: up to every other qubit of the register
            const Gate gate = controlledGate(1 + pick(count - 1), Unitary(2, 2000 + 97 * trial + step));
            simulator.ApplyGate(gate);
            EvolveCompact<MPO>(reference, gate, count);
            const int site = pick(count - 1);
            const Matrix local = Unitary(4, 3000 + 53 * trial + step);
            simulator.ApplyGate(Gate(local), site, site + 1);
            Evolve<MPO>(reference, local, {site, site + 1}, count);
            if (step == 4)
                // a read between the gates must leave no stale environments or cached forms
                (void)simulator.GetProbability(pick(count));
        }
        Close(Storage<MPO>(simulator), reference, "Multi-controlled gates on a tensor network mismatch");
    }

    // The named entry points and the batch interface.
    {
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        const auto apply = [&](const Gate &gate) { EvolveCompact<MPO>(reference, gate, count); };
        simulator.ApplyMCX(3, {0, 1, 5, 6});
        apply(Gate::MultiControlled(Gates::PauliXGate<>(), 3, {0, 1, 5, 6}));
        simulator.ApplyMCZ(0, {2, 4, 6}, 0b101);
        apply(Gate::MultiControlled(Gates::PauliZGate<>(), 0, {2, 4, 6}, 0b101));
        simulator.ApplyMCP(0.37, 6, {0, 1, 2, 3, 4, 5});
        apply(Gate::MultiControlled(Gates::PhaseShiftGate<>(0.37), 6, {0, 1, 2, 3, 4, 5}));
        simulator.ApplyMCRY(-0.83, 2, {6, 0, 4, 1}, 0b0110);
        apply(Gate::MultiControlled(Gates::RyGate<>(-0.83), 2, {6, 0, 4, 1}, 0b0110));
        simulator.ApplyMCU(0.4, -0.3, 0.9, 5, {1, 3, 0, 6}, 0.21);
        apply(Gate::MultiControlled(Gates::UGate<>(0.4, -0.3, 0.9, 0.21), 5, {1, 3, 0, 6}));
        const Gate batch[2] = {Gate::MultiControlled(Gates::HadamardGate<>(), 1, {0, 2, 3, 4, 5, 6}, 0b010101),
                               Gate::MultiControlled(Gates::PauliYGate<>(), 4, {5, 0, 3, 1})};
        simulator.ApplyGates({batch[0], batch[1]});
        apply(batch[0]);
        apply(batch[1]);
        // the identity on the selected subspace changes nothing, not even a bond
        const auto bonds = simulator.getBondDimensions();
        simulator.ApplyGate(Gate::MultiControlled(Matrix(Matrix::Identity(2, 2)), 0, {1, 2, 3, 4, 5, 6}));
        Require(simulator.getBondDimensions() == bonds, "A controlled identity changed bond dimensions");
        Close(Storage<MPO>(simulator), reference, "Named multi-controlled entry points on a tensor network mismatch");

        // Compact gates on at most three qubits go through their matrix.
        const std::vector<size_t> local = routed ? std::vector<size_t>{5, 0, 3} : std::vector<size_t>{3, 4, 2};
        std::vector<Complex> table(8);
        for (int i = 0; i < 8; ++i)
            table[i] = std::polar(1., 0.3 * i - 0.7);
        Gate gate = Gate::Diagonal(table, local);
        simulator.ApplyGate(gate);
        apply(gate);
        simulator.ApplyDiagonalGate(std::vector<Complex>{1., Complex(0., 1.), -1., Complex(0., -1.)},
                                    {local[1], local[0]});
        apply(Gate::Diagonal({1., Complex(0., 1.), -1., Complex(0., -1.)}, {local[1], local[0]}));
        gate = Gate::MultiplexedRotation(Gates::RotationAxis::Y, {0.3, -1.1, 2.2, 0.6}, local[1], {local[2], local[0]});
        simulator.ApplyGate(gate);
        apply(gate);
        simulator.ApplyMultiplexedRZ({0.9, -0.4}, local[0], {local[1]});
        apply(Gate::MultiplexedRotation(Gates::RotationAxis::Z, {0.9, -0.4}, local[0], {local[1]}));
        gate = Gate::MultiControlled(Unitary(2, 4242), local[2], {local[0], local[1]}, 0b01);
        simulator.ApplyGate(gate);
        apply(gate);
        gate = Gate::MultiControlled(Unitary(2, 4243), local[1], {});
        simulator.ApplyGate(gate);
        apply(gate);
        Close(Storage<MPO>(simulator), reference, "Small compact gates on a tensor network mismatch");

        // Malformed gates are rejected before anything changes.
        Throws<std::invalid_argument>([&] { simulator.ApplyDiagonalGate(std::vector<Complex>(8, 1.), {0, 1, 2, 3}); },
                                      "Tensor network accepted a diagonal table of the wrong size");
        Throws<std::invalid_argument>([&] { simulator.ApplyMultiplexedRX(std::vector<double>(8, 0.2), 0, {1, 0, 3}); },
                                      "Tensor network accepted a multiplexer whose target is a control");
        Throws<std::invalid_argument>([&] { simulator.ApplyMCX(2, {0, 2, 4, 5}); },
                                      "Tensor network accepted a target among its controls");
        Throws<std::invalid_argument>([&] { simulator.ApplyMCX(2, {0, 1, 4, 9}); },
                                      "Tensor network accepted a control outside the register");
        Matrix shrinking = Matrix::Identity(2, 2);
        shrinking(1, 1) = 0.5;
        if constexpr (MPO)
        {
            // rho -> O rho O^dagger also for an operator that is not unitary
            gate = Gate::MultiControlled(shrinking, 1, {0, 2, 4, 6}, 0b0011);
            simulator.ApplyGate(gate);
            apply(gate);
            simulator.ApplyOperator(gate);
            apply(gate);
        }
        else
            Throws<std::invalid_argument>(
                [&] { simulator.ApplyGate(Gate::MultiControlled(shrinking, 1, {0, 2, 4, 6})); },
                "MPS applied a wide controlled operator that is not unitary");
        Close(Storage<MPO>(simulator), reference, "A rejected compact gate changed the tensor-network state");
    }

    if constexpr (!MPO)
    {
        // After a Kraus operator the chain keeps explicit environments; the gate ends them
        // and what follows, reads included, sees the same state.
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        Matrix k0 = Matrix::Identity(2, 2), k1 = Matrix::Zero(2, 2);
        k0(1, 1) = 0.8;
        k1(0, 1) = 0.6;
        const size_t chosen = simulator.ApplyKrausOperator({k0, k1}, 0.35, 3);
        Evolve<false>(reference, chosen ? k1 : k0, {3}, count);
        reference /= reference.norm();
        Gate gate = Gate::MultiControlled(Unitary(2, 6001), 5, {0, 2, 3, 6}, 0b1001);
        simulator.ApplyGate(gate);
        EvolveCompact<false>(reference, gate, count);
        (void)simulator.GetProbability(2);
        gate = Gate::MultiControlled(Unitary(2, 6002), 1, {4, 5, 6, 0, 2});
        simulator.ApplyGate(gate);
        EvolveCompact<false>(reference, gate, count);
        Close(Storage<false>(simulator), reference, "Multi-controlled gate after a Kraus operator mismatch");
    }

    // With a bond limit the sweep after the gate cuts like a two-qubit gate does.
    for (int cap : {2, 5})
    {
        Simulator simulator(count);
        Matrix reference;
        start(simulator, reference);
        simulator.setLimitBondDimension(cap);
        for (int step = 0; step < 4; ++step)
            simulator.ApplyGate(controlledGate(3 + pick(count - 3), Unitary(2, 5000 + 11 * cap + step)));
        for (const auto dimension : simulator.getBondDimensions())
            Require(dimension <= cap, "A multi-controlled gate left a bond above the limit");
        // What the chain holds now evolves exactly once the limit is lifted.
        simulator.dontLimitBondDimension();
        reference = Storage<MPO>(simulator);
        Require(reference.allFinite(), "A truncated multi-controlled gate left a non-finite state");
        if constexpr (!MPO)
            Require(std::abs(reference.norm() - 1.) < 1E-9,
                    "A truncated multi-controlled gate left the MPS off its norm");
        const Gate gate = controlledGate(5, Unitary(2, 5100 + cap));
        simulator.ApplyGate(gate);
        EvolveCompact<MPO>(reference, gate, count);
        const Matrix local = Unitary(4, 5200 + cap);
        simulator.ApplyGate(Gate(local), 2, 3);
        Evolve<MPO>(reference, local, {2, 3}, count);
        Close(Storage<MPO>(simulator), reference, "Evolution after a truncated multi-controlled gate mismatch");
    }
}

// Independent full-operator tensor SVD reference. This deliberately computes
// U*S, unlike the simulator's inverse-free theta*V implementation.
inline Matrix TruncateOperator(const Matrix &rho, int count, int first, int cap)
{
    const int size = 1 << count, vectorSize = size * size;
    Eigen::VectorXcd tensor(vectorSize);
    for (int bra = 0; bra < size; ++bra)
        for (int ket = 0; ket < size; ++ket)
        {
            int index = 0;
            for (int q = 0; q < count; ++q)
                index |= ((((ket >> q) & 1) + 2 * ((bra >> q) & 1)) << (2 * q));
            tensor[index] = rho(ket, bra);
        }
    const int firstRows = 1 << (2 * (first + 2));
    const Matrix rightCut = Eigen::Map<const Matrix>(tensor.data(), firstRows, vectorSize / firstRows);
    Eigen::JacobiSVD<Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV> rightSvd(rightCut);
    const int rightRank = std::min<int>(cap, rightSvd.rank());
    const Matrix leftBlock =
        rightSvd.matrixU().leftCols(rightRank) * rightSvd.singularValues().head(rightRank).asDiagonal();
    const int secondRows = firstRows / 4;
    const Matrix leftCut = Eigen::Map<const Matrix>(leftBlock.data(), secondRows, 4 * rightRank);
    Eigen::JacobiSVD<Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV> leftSvd(leftCut);
    const int leftRank = std::min<int>(cap, leftSvd.rank());
    const Matrix compressedLeft = leftSvd.matrixU().leftCols(leftRank) *
                                  leftSvd.singularValues().head(leftRank).asDiagonal() *
                                  leftSvd.matrixV().leftCols(leftRank).adjoint();
    const Matrix reconstructed = Eigen::Map<const Matrix>(compressedLeft.data(), firstRows, rightRank) *
                                 rightSvd.matrixV().leftCols(rightRank).adjoint();
    const Eigen::Map<const Eigen::VectorXcd> flattened(reconstructed.data(), vectorSize);
    Matrix result(size, size);
    for (int bra = 0; bra < size; ++bra)
        for (int ket = 0; ket < size; ++ket)
        {
            int index = 0;
            for (int q = 0; q < count; ++q)
                index |= ((((ket >> q) & 1) + 2 * ((bra >> q) & 1)) << (2 * q));
            result(ket, bra) = flattened[index];
        }
    return result;
}

inline void MPOAudit()
{
    using Simulator = TensorNetworks::MPOSimulator;
    const Matrix unitary = Unitary(8, 961);
    for (int exponent : {-1000, -400, 400, 1000})
    {
        Simulator scaled(5);
        scaled.SetMultithreading(false);
        scaled.ApplyGate(Gate(Unitary(8, 82)), 1, 2, 3);
        scaled.ApplyGate(Gates::CNOTGate<>(), 0, 1);
        scaled.ApplyGate(Gates::CNOTGate<>(), 4, 3);
        scaled.ApplyAmplitudeDamping(2, 0.23);
        const Matrix before = scaled.getUnnormalizedDensityMatrix();
        scaled.ApplyGate(Gate(Matrix(unitary * std::ldexp(1., exponent))), 1, 2, 3);
        const auto snapshot = std::dynamic_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(scaled.getState());
        for (const auto &gamma : snapshot->gammas)
            for (Eigen::Index index = 0; index < gamma.size(); ++index)
                Require(std::isfinite(gamma.data()[index].real()) && std::isfinite(gamma.data()[index].imag()),
                        "Extreme MPO gate left non-finite tensors");
        if (std::abs(exponent) < 500)
        {
            Matrix expected = before;
            Evolve<true>(expected, unitary, {1, 2, 3}, 5);
            Matrix actual = scaled.getUnnormalizedDensityMatrix();
            TensorNetworks::ThreeSite::ScalePowerOfTwo(actual, -2 * exponent);
            Close(actual, expected, "Extreme representable MPO gate lost relative precision");
        }
        scaled.ApplyGate(Gate(Matrix(unitary.adjoint() * std::ldexp(1., -exponent))), 1, 2, 3);
        Close(scaled.getUnnormalizedDensityMatrix(), before, "Reciprocal extreme MPO gates changed the state");
    }
    {
        Simulator subnormal(3);
        subnormal.SetMultithreading(false);
        for (int exponent : {-1073, 500, 573})
            subnormal.ApplyGate(Gate(Matrix(Matrix::Identity(8, 8) * std::ldexp(1., exponent))), 2, 0, 1);
        Matrix expected = Matrix::Zero(8, 8);
        expected(0, 0) = 1.;
        Close(subnormal.getUnnormalizedDensityMatrix(), expected, "Subnormal MPO gate reciprocal overflowed");
    }
    {
        Simulator zero(3);
        zero.ApplyGate(Gate(Matrix::Zero(8, 8)), 0, 1, 2);
        Close(zero.getUnnormalizedDensityMatrix(), Matrix::Zero(8, 8), "Zero MPO operator was not preserved");
        zero.ApplyGate(Gate(unitary), 0, 1, 2);
        Close(zero.getUnnormalizedDensityMatrix(), Matrix::Zero(8, 8), "Evolution resurrected zero MPO operator");
    }
    Simulator compressed(5);
    compressed.SetMultithreading(false);
    compressed.setToMixtureOfBasisStates(std::vector<std::pair<size_t, double>>{{size_t(0), 0.63}, {size_t(21), 0.37}});
    for (int q = 0; q < 4; ++q)
        compressed.ApplyGate(Gate(Unitary(4, 300 + q)), q, q + 1);
    compressed.ApplyAmplitudeDamping(0, 0.17);
    Matrix expected = compressed.getUnnormalizedDensityMatrix();
    Evolve<true>(expected, unitary, {1, 2, 3}, 5);
    expected = TruncateOperator(expected, 5, 1, 2);
    compressed.setLimitBondDimension(2);
    compressed.ApplyGate(Gate(unitary), 1, 2, 3);
    Close(compressed.getUnnormalizedDensityMatrix(), expected, "MPO two-cut truncation differs from dense SVD");
    compressed.dontLimitBondDimension();
    for (int step = 0; step < 4; ++step)
    {
        compressed.GetProbability(2); // exercises a read before invalidating its cache
        compressed.ApplyGate(Gate(unitary), 1, 2, 3);
        Evolve<true>(expected, unitary, {1, 2, 3}, 5);
        compressed.ApplyGate(Gates::RyGate<>(0.21), 0);
        Evolve<true>(expected, Gates::RyGate<>(0.21).getRawOperatorMatrix(), {0}, 5);
        compressed.ApplyGate(Gates::CNOTGate<>(), 4, 3);
        Evolve<true>(expected, Gates::CNOTGate<>().getRawOperatorMatrix(), {4, 3}, 5);
        Close(compressed.getUnnormalizedDensityMatrix(), expected, "MPO truncated center or cache became stale");
    }
}

template <int D> void CheckPhysicalKernels()
{
    std::mt19937 random(2411 + D);
    std::uniform_real_distribution<double> distribution(-0.5, 0.5);
    Matrix filter = Matrix::Identity(8, 8);
    filter(0, 1) = Complex(0.3, -0.2);
    filter.row(3).setZero();
    Matrix denseFilter = Unitary(8, 2429);
    denseFilter.row(0) = Matrix::Identity(8, 8).row(0);
    denseFilter.row(3).setZero();
    denseFilter.row(6) *= Complex(0.2, -0.4);
    const std::array<Matrix, 9> gates{Gates::ToffoliGate<>().getRawOperatorMatrix(),
                                      Gates::CCZGate<>().getRawOperatorMatrix(),
                                      Gates::CCRxGate<>(0.41).getRawOperatorMatrix(),
                                      Gates::CCHGate<>().getRawOperatorMatrix(),
                                      Unitary(8, 2423),
                                      filter,
                                      denseFilter,
                                      Matrix::Zero(8, 8),
                                      Matrix::Identity(8, 8)};
    std::vector<int> leftDimensions{1, 3, 4, 8, 17, 64};
    if constexpr (D == 4)
        leftDimensions.insert(leftDimensions.end(), {127, 128, 129, 256});
    for (int L : leftDimensions)
        for (int R : {1, 3})
        {
            Matrix initial(D * D * L, D * R);
            for (Eigen::Index i = 0; i < initial.size(); ++i)
                initial.data()[i] = Complex(distribution(random), distribution(random));
            const Eigen::VectorXd lambda = Eigen::VectorXd::LinSpaced(L, 0.1, 0.9);
            const Matrix weighted = TensorNetworks::ThreeSite::Weighted(initial, &lambda, L);
            for (Eigen::Index col = 0; col < initial.cols(); ++col)
                for (Eigen::Index row = 0; row < initial.rows(); ++row)
                    Close(weighted(row, col), initial(row, col) * lambda[row % L],
                          "Three-site block weighting mismatch");
            std::array<int, 3> sites{0, 1, 2};
            do
            {
                for (const Matrix &gate : gates)
                {
                    const Matrix ordered = Embed(gate, {sites[0], sites[1], sites[2]}, 3);
                    Matrix actual = initial, expected = initial;
                    for (int r = 0; r < R; ++r)
                        for (int l = 0; l < L; ++l)
                        {
                            Matrix block(8, D == 2 ? 1 : 8);
                            const auto position = [&](int ket, int bra) {
                                const int p1 = (ket & 1) + (D == 4 ? 2 * (bra & 1) : 0);
                                const int p2 = ((ket >> 1) & 1) + (D == 4 ? 2 * ((bra >> 1) & 1) : 0);
                                const int p3 = ((ket >> 2) & 1) + (D == 4 ? 2 * ((bra >> 2) & 1) : 0);
                                return std::pair<int, int>{l + L * (p1 + D * p2), p3 + D * r};
                            };
                            for (int bra = 0; bra < block.cols(); ++bra)
                                for (int ket = 0; ket < 8; ++ket)
                                {
                                    const auto [row, col] = position(ket, bra);
                                    block(ket, bra) = initial(row, col);
                                }
                            if constexpr (D == 2)
                                block = (ordered * block).eval();
                            else
                                block = (ordered * block * ordered.adjoint()).eval();
                            for (int bra = 0; bra < block.cols(); ++bra)
                                for (int ket = 0; ket < 8; ++ket)
                                {
                                    const auto [row, col] = position(ket, bra);
                                    expected(row, col) = block(ket, bra);
                                }
                        }
                    TensorNetworks::ThreeSite::Apply<D>(actual, gate, 0, sites[0], sites[1], sites[2], L, R);
                    Close(actual, expected, "Three-site contiguous physical kernel differs from dense oracle");
                }
            } while (std::next_permutation(sites.begin(), sites.end()));
        }
}

inline void MPOOperatorAPI()
{
    using Simulator = TensorNetworks::MPOSimulator;
    using Interface = TensorNetworks::MPOSimulatorInterface;
    Matrix filter = Unitary(8, 2011);
    for (int row = 0; row < 8; ++row)
        filter.row(row) *= 0.25 + 0.07 * row;
    std::array<int, 3> sites{0, 2, 4};
    do
    {
        Simulator raw(5), normalized(5);
        raw.SetMultithreading(false);
        normalized.SetMultithreading(false);
        raw.setToMixtureOfBasisStates(std::vector<std::pair<size_t, double>>{{0, 0.6}, {21, 0.4}});
        normalized.setToMixtureOfBasisStates(std::vector<std::pair<size_t, double>>{{0, 0.6}, {21, 0.4}});
        Matrix expected = raw.getUnnormalizedDensityMatrix();
        Evolve<true>(expected, filter, {sites[0], sites[1], sites[2]}, 5);
        Interface &rawApi = raw, &normalizedApi = normalized;
        rawApi.ApplyOperator(Gate(filter), sites[0], sites[1], sites[2]);
        normalizedApi.ApplyOperatorAndNormalize(Gate(filter, sites[0], sites[1], sites[2]));
        Close(raw.getUnnormalizedDensityMatrix(), expected, "Three-site MPO operator API mismatch");
        Close(normalized.getUnnormalizedDensityMatrix(), Matrix(expected / expected.trace()),
              "Three-site normalized MPO operator API mismatch");
        const Gate inverse(Matrix(Unitary(8, 2099)), sites[2], sites[0], sites[1]);
        rawApi.ApplyOperators({Gate(Gates::RyGate<>(0.21).getRawOperatorMatrix(), 1), inverse});
        Evolve<true>(expected, Gates::RyGate<>(0.21).getRawOperatorMatrix(), {1}, 5);
        Evolve<true>(expected, inverse.getRawOperatorMatrix(), {sites[2], sites[0], sites[1]}, 5);
        Close(raw.getUnnormalizedDensityMatrix(), expected, "Mixed-arity MPO operator batch mismatch");
    } while (std::next_permutation(sites.begin(), sites.end()));

    for (bool route : {false, true})
    {
        Simulator rejected(5);
        rejected.SetMultithreading(false);
        rejected.ApplyGate(Gates::HadamardGate<>(), 0);
        rejected.ApplyGate(Gates::CNOTGate<>(), 1, 0);
        rejected.SaveState();
        const Matrix before = rejected.getUnnormalizedDensityMatrix();
        const auto map = rejected.getQubitsMap(), dimensions = rejected.getBondDimensions();
        rejected.setLimitBondDimension(1);
        rejected.setRestoreTraceAfterTruncation(true);
        rejected.setHermitizeAfterTruncation(true);
        int callbacks = 0;
        rejected.SetBondDimensionCallback([&](const auto &) { ++callbacks; });
        Matrix projector = Matrix::Zero(8, 8);
        projector(7, 7) = 1.;
        Throws<std::runtime_error>(
            [&] { rejected.ApplyOperatorAndNormalize(Gate(projector), route ? 0 : 1, 2, route ? 4 : 3); },
            "Zero-weight three-site MPO postselection was accepted");
        Require(rejected.getQubitsMap() == map && rejected.getBondDimensions() == dimensions && callbacks == 0,
                "Failed three-site MPO postselection changed routing or callbacks");
        Close(rejected.getUnnormalizedDensityMatrix(), before, "Failed three-site MPO postselection changed state");
        rejected.dontLimitBondDimension();
        rejected.ApplyGate(Gates::PauliXGate<>(), 4);
        rejected.RestoreState();
        Close(rejected.getUnnormalizedDensityMatrix(), before, "MPO postselection overwrote the saved state");
        Throws<std::invalid_argument>([&] { rejected.ApplyOperator(Gate(Matrix::Identity(16, 16))); },
                                      "MPO accepted a four-qubit operator");
        Throws<std::invalid_argument>([&] { rejected.ApplyOperatorAndNormalize(Gate(Matrix::Identity(16, 16))); },
                                      "MPO accepted a normalized four-qubit operator");
    }
    TensorNetworks::MPOSimulatorImpl adjacent(3);
    adjacent.SetMultithreading(false);
    Interface &api = adjacent;
    api.ApplyOperator(Gate(filter, 2, 0, 1));
    Matrix expected = Matrix::Zero(8, 8);
    expected(0, 0) = 1.;
    Evolve<true>(expected, filter, {2, 0, 1}, 3);
    Close(adjacent.getUnnormalizedDensityMatrix(), expected, "Direct MPO three-site operator mismatch");
    api.ApplyOperatorAndNormalize(Gate(filter), 1, 2, 0);
    Evolve<true>(expected, filter, {1, 2, 0}, 3);
    expected /= expected.trace();
    Close(adjacent.getUnnormalizedDensityMatrix(), expected, "Direct normalized MPO three-site operator mismatch");
}

template <bool MPO> void Check()
{
    using Simulator = std::conditional_t<MPO, TensorNetworks::MPOSimulator, TensorNetworks::MPSSimulator>;
    using Interface =
        std::conditional_t<MPO, TensorNetworks::MPOSimulatorInterface, TensorNetworks::MPSSimulatorInterface>;
    using Impl = std::conditional_t<MPO, TensorNetworks::MPOSimulatorImpl, TensorNetworks::MPSSimulatorImpl>;
    CheckLegacyDispatch<MPO, Simulator>();
    CheckLegacyDispatch<MPO, Impl>();
    CheckWideRejection<MPO, Simulator>();
    CheckWideRejection<MPO, Impl>();
    constexpr int count = 5, size = 1 << count;
    const Matrix dense = Unitary(8, 831);
    const Matrix sparse = Gates::ToffoliGate<>().getRawOperatorMatrix();
    Matrix diagonal = Matrix::Identity(8, 8);
    for (int i = 0; i < 8; ++i)
        diagonal(i, i) = std::polar(1., 0.17 * i);
    const std::array<Matrix, 3> gates{dense, sparse, diagonal};
    for (bool remapped : {false, true})
    {
        std::array<int, 3> sites{0, 2, 4};
        do
        {
            Simulator simulator(count);
            simulator.SetMultithreading(false);
            Matrix reference = Matrix::Zero(size, MPO ? size : 1);
            reference(0, 0) = 1.;
            if constexpr (MPO)
            {
                simulator.setToMixtureOfBasisStates(
                    std::vector<std::pair<size_t, double>>{{size_t(0), 0.7}, {size_t(19), 0.3}});
                reference(0, 0) = 0.7;
                reference(19, 19) = 0.3;
            }
            for (int q = 0; q < count; ++q)
            {
                const Gates::RyGate<> rotation(0.24 + 0.16 * q);
                simulator.ApplyGate(rotation, q);
                Evolve<MPO>(reference, rotation.getRawOperatorMatrix(), {q}, count);
            }
            for (int q = 0; q + 1 < count; ++q)
            {
                const Gates::CNOTGate<> cx;
                simulator.ApplyGate(cx, q + 1, q);
                Evolve<MPO>(reference, cx.getRawOperatorMatrix(), {q + 1, q}, count);
            }
            if (remapped)
            {
                const Gates::SwapGate<> swap;
                simulator.ApplyGate(swap, 0, 4);
                Evolve<MPO>(reference, swap.getRawOperatorMatrix(), {0, 4}, count);
            }
            // Prime read caches before modifying three sites.
            simulator.GetProbability(2);
            const auto before = simulator.getQubitsMap();
            std::array<int, 3> physicalOrder = sites;
            std::sort(physicalOrder.begin(), physicalOrder.end(), [&](int a, int b) { return before[a] < before[b]; });
            const int middleLogical = physicalOrder[1];
            const auto middle = before[middleLogical];
            int callbacks = 0, summaries = 0, meetings = 0;
            simulator.SetBondDimensionCallback([&](const auto &dimensions) {
                Require(dimensions.size() == count - 1, "Three-site bond callback returned wrong size");
                ++callbacks;
            });
            simulator.SetBondDimensionSummaryCallback([&](auto maximum) {
                Require(maximum == simulator.getMaxBondDimension(), "Three-site summary callback was stale");
                ++summaries;
            });
            simulator.SetMeetingPositionCallback([&](const auto &) {
                ++meetings;
                return 0;
            });
            Interface &api = simulator;
            api.ApplyGate(Gate(dense), sites[0], sites[1], sites[2]);
            Evolve<MPO>(reference, dense, {sites[0], sites[1], sites[2]}, count);
            Close(Storage<MPO>(simulator), reference, "Three-site dense gate/order mismatch");
            Require(simulator.getQubitsMap()[middleLogical] == middle, "Three-site routing moved the median site");
            const auto &after = simulator.getQubitsMap();
            Require(after[physicalOrder[0]] == middle - 1 && after[physicalOrder[2]] == middle + 1,
                    "Three-site routing did not place the endpoints beside the median");
            Require(callbacks == before[physicalOrder[2]] - before[physicalOrder[0]] - 1 && callbacks == summaries &&
                        meetings == 0,
                    "Three-site routing did not preserve callback semantics");
            simulator.SetBondDimensionCallback(nullptr);
            simulator.SetBondDimensionSummaryCallback(nullptr);
            simulator.SetMeetingPositionCallback(nullptr);
            for (int index = 1; index < 3; ++index)
            {
                simulator.ApplyGate(Gate(gates[index], sites[0], sites[1], sites[2]));
                Evolve<MPO>(reference, gates[index], {sites[0], sites[1], sites[2]}, count);
                Close(Storage<MPO>(simulator), reference, "Three-site sparse/diagonal gate mismatch");
            }
            const auto snapshot = simulator.getState();
            auto clone = simulator.Clone();
            const Gate inverse(Matrix(dense.adjoint()), sites[2], sites[0], sites[1]);
            simulator.ApplyGates({inverse});
            clone->ApplyGate(inverse);
            Close(Storage<MPO>(*clone), Storage<MPO>(simulator), "Three-site clone/batch mismatch");
            simulator.setState(snapshot);
            Close(Storage<MPO>(simulator), reference, "Three-site snapshot restore mismatch");

            const auto map = simulator.getQubitsMap();
            Throws<std::invalid_argument>([&] { simulator.ApplyGate(Gate(dense), 0, 0, 4); },
                                          "Duplicate three-site target accepted");
            Throws<std::invalid_argument>([&] { simulator.ApplyGate(Gate(dense), 0, 2, count); },
                                          "Out-of-range third target accepted");
            Matrix invalid = dense;
            invalid(0, 0) = std::numeric_limits<double>::infinity();
            Throws<std::invalid_argument>([&] { simulator.ApplyGate(Gate(invalid), 0, 2, 4); },
                                          "Non-finite three-site matrix accepted");
            Require(simulator.getQubitsMap() == map, "Invalid three-site gate mutated routing");
            Close(Storage<MPO>(simulator), reference, "Invalid three-site gate mutated state");
        } while (std::next_permutation(sites.begin(), sites.end()));
    }

    // Both exposed impl and wrapper APIs accept all six local argument orders.
    std::array<int, 3> local{1, 2, 3};
    do
    {
        Impl simulator(count);
        simulator.SetMultithreading(false);
        simulator.ApplyGate(Gate(dense, local[0], local[1], local[2]));
        Matrix reference = Matrix::Zero(size, MPO ? size : 1);
        reference(0, 0) = 1.;
        Evolve<MPO>(reference, dense, {local[0], local[1], local[2]}, count);
        Close(Storage<MPO>(simulator), reference, "Direct three-site implementation order mismatch");
    } while (std::next_permutation(local.begin(), local.end()));

    {
        Simulator boundary(3);
        boundary.SetMultithreading(false);
        Matrix expected = Matrix::Zero(8, MPO ? 8 : 1);
        expected(0, 0) = 1.;
        boundary.ApplyGate(Gate(dense), 2, 0, 1);
        Evolve<MPO>(expected, dense, {2, 0, 1}, 3);
        Close(Storage<MPO>(boundary), expected, "Three-site full-chain boundary mismatch");
    }

    // Unoccupied huge columns must not set the binary scale of a tiny occupied
    // column. Test all argument orders and recover the MPO's unrepresentable
    // intermediate trace with two finite reciprocal powers of two.
    std::array<int, 3> extremeSites{0, 1, 2};
    do
    {
        Simulator supported(3);
        supported.SetMultithreading(false);
        const size_t basis = (size_t{1} << extremeSites[0]) | (size_t{1} << extremeSites[2]);
        supported.setToBasisState(basis);
        Matrix extreme = Matrix::Identity(8, 8) * std::ldexp(1., 990);
        extreme(5, 5) = std::ldexp(1., -1030);
        supported.ApplyGate(Gate(extreme), extremeSites[0], extremeSites[1], extremeSites[2]);
        if constexpr (MPO)
            for (int step = 0; step < 2; ++step)
                supported.ApplyGate(Gate(Matrix(Matrix::Identity(8, 8) * std::ldexp(1., 515))), extremeSites[0],
                                    extremeSites[1], extremeSites[2]);
        Matrix expected = Matrix::Zero(8, MPO ? 8 : 1);
        expected(basis, MPO ? basis : 0) = 1.;
        Close(Storage<MPO>(supported), expected, "Inactive huge three-qubit columns erased a tiny active coefficient");
    } while (std::next_permutation(extremeSites.begin(), extremeSites.end()));

    {
        Simulator filtered(count);
        filtered.SetMultithreading(false);
        filtered.ApplyGate(Gate(dense), 1, 2, 3);
        filtered.ApplyGate(Gates::CNOTGate<>(), 0, 1);
        filtered.ApplyGate(Gates::CNOTGate<>(), 4, 3);
        Matrix expected = Storage<MPO>(filtered);
        Matrix filter = Matrix::Identity(8, 8);
        for (int i = 0; i < 8; ++i)
            filter(i, i) = 0.1 + 0.1 * i;
        filtered.ApplyGate(Gate(filter), 1, 2, 3);
        Evolve<MPO>(expected, filter, {1, 2, 3}, count);
        if constexpr (!MPO)
            expected /= expected.norm();
        Close(Storage<MPO>(filtered), expected, "Three-site finite nonunitary update mismatch");
        filtered.ApplyGate(Gates::CNOTGate<>(), 4, 0);
        Evolve<MPO>(expected, Gates::CNOTGate<>().getRawOperatorMatrix(), {4, 0}, count);
        Close(Storage<MPO>(filtered), expected, "Three-site filter left stale outer environments");
    }

    // A nonunitary noise update changes the environments and cached spectra.
    Simulator noisy(count);
    noisy.SetMultithreading(false);
    for (int q = 0; q < count; ++q)
        noisy.ApplyGate(Gates::HadamardGate<>(), q);
    noisy.ApplyGate(Gates::CNOTGate<>(), 4, 0);
    Matrix reference = Storage<MPO>(noisy);
    Matrix k0 = Matrix::Identity(2, 2), k1 = Matrix::Zero(2, 2);
    k0(1, 1) = std::sqrt(0.7);
    k1(0, 1) = std::sqrt(0.3);
    if constexpr (MPO)
    {
        noisy.ApplyKrausOperators({k0, k1}, 2);
        Matrix first = reference, second = reference;
        Evolve<true>(first, k0, {2}, count);
        Evolve<true>(second, k1, {2}, count);
        reference = first + second;
    }
    else
    {
        const size_t selected = noisy.ApplyKrausOperator({k0, k1}, 0.37, 2);
        Evolve<false>(reference, selected == 0 ? k0 : k1, {2}, count);
        reference /= reference.norm();
    }
    for (int step = 0; step < 3; ++step)
    {
        noisy.ApplyGate(Gate(gates[step], 4, 0, 2));
        Evolve<MPO>(reference, gates[step], {4, 0, 2}, count);
        Close(Storage<MPO>(noisy), reference, "Three-site evolution after noise mismatch");
        noisy.ApplyGate(Gates::CNOTGate<>(), 1, 4);
        Evolve<MPO>(reference, Gates::CNOTGate<>().getRawOperatorMatrix(), {1, 4}, count);
        Close(Storage<MPO>(noisy), reference, "Two-site evolution after noisy three-site gate mismatch");
    }
    for (auto mode : {Simulator::TruncationMode::DiscardedWeight, Simulator::TruncationMode::RelativeToMax})
    {
        Simulator compressed(count);
        compressed.SetMultithreading(false);
        compressed.setLimitBondDimension(2);
        compressed.setLimitEntanglement(1E-6);
        compressed.setTruncationMode(mode);
        if constexpr (MPO)
        {
            compressed.setHermitizeAfterTruncation(true);
            compressed.setRestoreTraceAfterTruncation(true);
        }
        for (int step = 0; step < 4; ++step)
            compressed.ApplyGate(Gate(dense, step % count, (step + 1) % count, (step + 3) % count));
        const Matrix state = Storage<MPO>(compressed);
        Require(state.allFinite() && compressed.getMaxBondDimension() <= 2,
                "Three-site compression produced invalid tensors or exceeded bond limit");
        if constexpr (MPO)
            Close(state.trace(), 1., "Three-site MPO trace restoration failed");
        else
            Close(state.squaredNorm(), 1., "Three-site MPS normalization failed");
        // Start from the actual compressed state, then verify subsequent exact
        // evolution. This detects stale spectra/environments without assuming
        // a particular approximate SVD basis.
        Matrix subsequent = Storage<MPO>(compressed);
        compressed.dontLimitBondDimension();
        compressed.dontLimitEntanglement();
        compressed.ApplyGate(Gate(dense, 1, 3, 4));
        Evolve<MPO>(subsequent, dense, {1, 3, 4}, count);
        compressed.ApplyGate(Gates::RyGate<>(0.37), 2);
        Evolve<MPO>(subsequent, Gates::RyGate<>(0.37).getRawOperatorMatrix(), {2}, count);
        compressed.ApplyGate(Gates::CNOTGate<>(), 4, 0);
        Evolve<MPO>(subsequent, Gates::CNOTGate<>().getRawOperatorMatrix(), {4, 0}, count);
        Close(Storage<MPO>(compressed), subsequent, "Evolution after three-site compression mismatch");
        compressed.ReCanonicalize();
        Close(Storage<MPO>(compressed), subsequent, "Three-site recanonicalization changed represented state");
        const auto saved = compressed.getState();
        compressed.ApplyGate(Gates::PauliXGate<>(), 0);
        compressed.setState(saved);
        Close(Storage<MPO>(compressed), subsequent, "Restoring after three-site compression changed represented state");
        if constexpr (MPO)
            compressed.ApplyAmplitudeDamping(3, 0.2);
        else
            compressed.ApplyKrausOperator({k0, k1}, 0.2, 3);
        compressed.ApplyGate(Gate(dense, 0, 2, 4));
        Require(Storage<MPO>(compressed).allFinite(), "Noise after restored three-site compression failed");
    }
}

inline bool Run(bool mpo)
{
    try
    {
        if (mpo)
        {
            CheckPhysicalKernels<4>();
            Check<true>();
            MPOAudit();
            MPOOperatorAPI();
            CheckCompactGates<true, TensorNetworks::MPOSimulator>();
            CheckCompactGates<true, TensorNetworks::MPOSimulatorImpl>();
            CheckWideTables<true, TensorNetworks::MPOSimulator>();
            CheckWideTables<true, TensorNetworks::MPOSimulatorImpl>();
        }
        else
        {
            CheckPhysicalKernels<2>();
            Check<false>();
            CheckCompactGates<false, TensorNetworks::MPSSimulator>();
            CheckCompactGates<false, TensorNetworks::MPSSimulatorImpl>();
            CheckWideTables<false, TensorNetworks::MPSSimulator>();
            CheckWideTables<false, TensorNetworks::MPSSimulatorImpl>();
        }
        std::cout << "Three-qubit tensor routing, contractions, splits and lifecycle: passed\n";
        return true;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Three-qubit tensor regression: " << error.what() << '\n';
        return false;
    }
}
} // namespace ThreeQubitTensorRegression
} // namespace QC
