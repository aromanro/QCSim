#pragma once

#include "PathIntegral.h"
#include "QubitRegister.h"
#include <array>
#include <iostream>
#include <random>

namespace QC
{
namespace PathStructuredRegression
{
using Gate = Gates::AppliedGate<>;
using Simulator = PathIntegral::PathIntegralSimulator;
using Complex = std::complex<double>;

inline void Require(bool value, const char *message)
{
    if (!value)
        throw std::runtime_error(message);
}

inline std::vector<bool> Bits(size_t value, size_t width)
{
    std::vector<bool> bits(width);
    for (size_t q = 0; q < std::min(width, size_t{64}); ++q)
        bits[q] = ((value >> q) & 1) != 0;
    return bits;
}

inline void DenseComparisons()
{
    std::vector<Gate> gates;
    for (size_t q = 0; q < 7; ++q)
        gates.emplace_back(Gates::UGate<>(.731 + q * .13, -.29, .13 * q).getRawOperatorMatrix(), q);
    std::vector<Complex> diagonal(32);
    for (size_t i = 0; i < diagonal.size(); ++i)
        diagonal[i] = std::polar(1., -.41 + .17 * i);
    gates.push_back(Gate::Diagonal(diagonal, {6, 1, 4, 0, 3}));
    gates.push_back(Gate::MultiControlled(Gates::PauliXGate<>(), 2, {6, 1, 5, 0, 3}, 21));
    gates.push_back(Gate::MultiControlled(Gates::UGate<>(.391, -.73, .17, .31), 4, {0, 6, 2, 1}, 10));
    for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
    {
        std::vector<double> angles(32);
        for (size_t i = 0; i < angles.size(); ++i)
            angles[i] = -.117 + i * .073;
        gates.push_back(Gate::MultiplexedRotation(axis, angles, 0, {5, 3, 6, 2, 4}));
    }
    gates.push_back(Gate::Diagonal({std::polar(1., .393)}, {}));
    gates.emplace_back(Gates::HadamardGate<>().getRawOperatorMatrix(), 6);
    QubitRegister<> reference(7);
    reference.ApplyGates(gates);
    for (size_t split : {size_t{0}, size_t{1}, std::numeric_limits<size_t>::max()})
        for (bool growth : {false, true})
            for (bool cache : {false, true})
            {
                Simulator sim;
                sim.SetTrimValue(0.);
                sim.SetMaxDoublingsForBackwardPaths(split);
                sim.SetUseGrowthAwareSplit(growth);
                if (!cache)
                {
                    sim.SetCompiledGateCacheMaxBytes(0);
                    sim.SetForwardCacheMaxBytes(0);
                }
                sim.SetCircuit(gates);
                for (unsigned repeat = 0; repeat < 2; ++repeat)
                    for (size_t output = 0; output < 128; ++output)
                        Require(std::abs(sim.Propagate(Bits(output, 7)) - reference.getBasisStateAmplitude(output)) <
                                    4e-12,
                                "Structured path query disagrees with statevector");
                auto copy = sim;
                sim.Reset();
                Require(std::abs(copy.Propagate(Bits(117, 7)) - reference.getBasisStateAmplitude(117)) < 4e-12,
                        "Structured compiled path cache lost payload ownership");
                copy.SetMaxDoublingsForBackwardPaths(split == 0 ? 2 : 0);
                copy.SetUseGrowthAwareSplit(!growth);
                Require(std::abs(copy.Propagate(Bits(35, 7)) - reference.getBasisStateAmplitude(35)) < 4e-12,
                        "Structured path transpose/repartition mismatch");
            }
    Simulator all;
    all.SetTrimValue(0.);
    all.PropagateAll(gates, Bits(0, 7));
    for (size_t output = 0; output < 128; ++output)
    {
        const auto found = all.GetAmplitudes().find(PathIntegral::FastVectorBool(Bits(output, 7)));
        Require(found != all.GetAmplitudes().end() &&
                    std::abs(found->second - reference.getBasisStateAmplitude(output)) < 4e-12,
                "Structured full path frontier mismatch");
    }
    PathIntegral::AmplitudeMap stepped;
    stepped[PathIntegral::FastVectorBool(Bits(0, 7))] = 1.;
    for (const auto &gate : gates)
        all.PropagateStep(gate, stepped);
    for (const auto &entry : stepped)
    {
        size_t index = 0;
        for (size_t q = 0; q < 7; ++q)
            index |= size_t(entry.first.get(q)) << q;
        Require(std::abs(entry.second - reference.getBasisStateAmplitude(index)) < 4e-12,
                "Structured step frontier mismatch");
    }
}

inline void WideAndTiny()
{
    std::vector<size_t> controls(63);
    std::vector<bool> start(130, false);
    for (size_t q = 0; q < controls.size(); ++q)
    {
        controls[q] = q + 64;
        start[q + 64] = true;
    }
    const Gate controlled = Gate::MultiControlled(Gates::PauliXGate<>(), 129, controls);
    const Gate phase = Gate::Diagonal({1., Complex(0., 1.), -1., Complex(0., -1.)}, {129, 80});
    auto output = start;
    output[129] = true;
    for (size_t split : {size_t{0}, std::numeric_limits<size_t>::max()})
    {
        Simulator wide;
        wide.SetTrimValue(0.);
        wide.SetMaxDoublingsForBackwardPaths(split);
        wide.SetCircuit({controlled, phase});
        for (unsigned repeat = 0; repeat < 3; ++repeat)
            Require(wide.Propagate(start, output) == Complex(0., -1.),
                    "64-local-bit structured gate failed on a wide register");
        wide.PropagateAll({controlled, phase}, start);
        Require(wide.GetAmplitudes().size() == 1 && wide.QubitProbability(129) == 1.,
                "Wide structured gate expanded a permutation");
    }
    const double angle = 1e-200;
    const auto rotation = Gate::MultiControlled(Gates::RyGate<>(angle), 0, {});
    for (const auto &gate : {rotation, Gate(Gates::RyGate<>(angle).getRawOperatorMatrix(), 0)})
    {
        Simulator tiny;
        tiny.SetTrimValue(0.);
        tiny.SetCircuit({gate});
        for (unsigned repeat = 0; repeat < 3; ++repeat)
            Require(tiny.Propagate(Bits(1, 1)) == Complex(.5 * angle, 0.),
                    "Exact path mode dropped a representable tiny coefficient");
        tiny.PropagateAll({gate}, Bits(0, 1));
        const auto found = tiny.GetAmplitudes().find(PathIntegral::FastVectorBool(Bits(1, 1)));
        Require(found != tiny.GetAmplitudes().end() && found->second == Complex(.5 * angle, 0.),
                "Tiny full-frontier coefficient was lost");
    }
}

inline void SymmetricPayloads()
{
    const std::vector<Gate> symmetric{
        Gate::Diagonal({1., Complex(0., 1.), -1., Complex(0., -1.)}, {2, 0}),
        Gate::MultiplexedRotation(Gates::RotationAxis::X, {.13, -.27, .39, -.51}, 1, {2, 0}),
        Gate::MultiplexedRotation(Gates::RotationAxis::Z, {.13, -.27, .39, -.51}, 1, {2, 0}),
        Gate::MultiControlled(Gates::PauliXGate<>(), 1, {2, 0})};
    for (const auto &gate : symmetric)
    {
        const Gate shared = gate;
        Require(shared.getStructuredGate() == gate.getStructuredGate(),
                "Copying a symmetric compact gate duplicated its immutable payload");
        std::vector<Gate> circuit;
        for (size_t q = 0; q < 3; ++q)
            circuit.emplace_back(Gates::HadamardGate<>().getRawOperatorMatrix(), q);
        circuit.push_back(shared);
        QubitRegister<> reference(3);
        reference.ApplyGates(circuit);
        Simulator simulator;
        simulator.SetTrimValue(0.);
        simulator.SetCircuit(circuit);
        for (size_t split : {size_t{0}, std::numeric_limits<size_t>::max(), size_t{1}})
        {
            simulator.SetMaxDoublingsForBackwardPaths(split);
            for (size_t output = 0; output < 8; ++output)
                Require(std::abs(simulator.Propagate(Bits(output, 3)) - reference.getBasisStateAmplitude(output)) <
                            1e-14,
                        "Sharing a symmetric transposed payload changed path propagation");
        }
    }
}

template <class Exception, class Action> bool Throws(Action action)
{
    try
    {
        action();
    }
    catch (const Exception &)
    {
        return true;
    }
    catch (...)
    {
    }
    return false;
}

// A compact gate on at most three qubits must be listed exactly like the
// dense gate with the same matrix, and so run through the dense kernels.
inline void LoweredTables()
{
    using Sparse = PathIntegral::Detail::SparseGate;
    Eigen::MatrixXcd skew(2, 2);
    skew << Complex(.6, 0.), Complex(0., -.8), Complex(.8, 0.), Complex(0., .6);
    const std::vector<Gate> fewQubits{
        Gate::Diagonal({1., Complex(0., 1.)}, {4}),
        Gate::Diagonal({1., Complex(0., 1.), 0., Complex(0., -1.)}, {5, 2}),
        Gate::Diagonal({1., -1., Complex(0., 1.), 1., std::polar(1., .3), 1., 1., -1.}, {1, 6, 3}),
        Gate::MultiControlled(skew, 3, {}),
        Gate::MultiControlled(Gates::PauliXGate<>(), 2, {5}),
        Gate::MultiControlled(skew, 2, {5}, 0),
        Gate::MultiControlled(skew, 4, {0, 6}, 2),
        Gate::MultiplexedRotation(Gates::RotationAxis::Y, {.37}, 1, {}),
        Gate::MultiplexedRotation(Gates::RotationAxis::X, {0., 1.1}, 1, {4}),
        Gate::MultiplexedRotation(Gates::RotationAxis::Y, {.3, 0., std::acos(-1.), -1.7}, 6, {2, 0}),
        Gate::MultiplexedRotation(Gates::RotationAxis::Z, {.3, 0., -.9, 2.1}, 0, {3, 1})};
    for (const auto &source : fewQubits)
        for (const Gate &gate : {source, source.Adjoint(), source.Transposed()})
        {
            const auto &operation = *gate.getStructuredGate();
            const auto ordered = operation.OrderedQubits();
            const Gate dense(operation.DenseMatrix<Eigen::MatrixXcd>(), ordered[0], ordered.size() > 1 ? ordered[1] : 0,
                             ordered.size() > 2 ? ordered[2] : 0);
            for (double epsilon : {0., 1E-15})
            {
                const Sparse lowered(gate, epsilon), reference(dense, epsilon);
                Require(!lowered.structured && lowered.arity == reference.arity &&
                            lowered.counts == reference.counts && lowered.maxFanout == reference.maxFanout &&
                            lowered.diagonal == reference.diagonal && lowered.injective == reference.injective,
                        "A small compact gate was not lowered to the dense transition metadata");
                for (unsigned column = 0; column < (1u << lowered.arity); ++column)
                    for (unsigned t = 0; t < lowered.counts[column]; ++t)
                        Require(lowered.rows[column][t] == reference.rows[column][t] &&
                                    lowered.values[column][t] == reference.values[column][t],
                                "A lowered compact gate lists different transitions than its dense matrix");
            }
        }
    const Sparse wide(Gate::MultiControlled(Gates::PauliXGate<>(), 0, {1, 2, 3}), 0.),
        phase(Gate::Diagonal({Complex(0., 1.)}, {}), 0.);
    Require(wide.structured && wide.arity == 4 && phase.structured && phase.arity == 0,
            "A compact gate on four qubits, or on none, must stay compact");
    // The adjoint and the transpose of the dense matrix, from the shared tables.
    for (const auto &gate : fewQubits)
    {
        const auto matrix = gate.getStructuredGate()->DenseMatrix<Eigen::MatrixXcd>();
        Require((gate.Adjoint().getStructuredGate()->DenseMatrix<Eigen::MatrixXcd>() - matrix.adjoint()).norm() == 0. &&
                    (gate.Transposed().getStructuredGate()->DenseMatrix<Eigen::MatrixXcd>() - matrix.transpose())
                            .norm() == 0.,
                "The adjoint or transpose of a compact gate is not that of its matrix");
    }
}

// The selector of a state is gathered from whole words; it must be the
// bit-by-bit value for every qubit order, including runs that end at a word
// boundary and controls spread over several words.
inline void WordSelectors()
{
    std::mt19937_64 engine(77);
    const std::vector<std::vector<size_t>> lists{
        {},
        {5},
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
        {9, 8, 7, 3, 4, 5, 0},
        {62, 63, 64, 65, 66},
        {63, 127, 128, 0, 191, 190, 64},
        {130, 131, 132, 70, 71, 5, 133, 134, 199},
    };
    for (const auto &qubits : lists)
        for (unsigned trial = 0; trial < 50; ++trial)
        {
            std::array<uint64_t, 4> words{engine(), engine(), engine(), engine()};
            size_t expected = 0;
            for (size_t i = 0; i < qubits.size(); ++i)
                expected |= size_t((words[qubits[i] / 64] >> (qubits[i] % 64)) & 1) << i;
            const size_t wanted = qubits.empty() ? 0 : engine() % (size_t{1} << qubits.size());
            const Gates::StructuredSelector selector(qubits, wanted);
            Require(selector.Gather(words.data()) == expected && selector.Matches(words.data()) == (expected == wanted),
                    "A selector gathered from words differs from the bit-by-bit value");
            // Force a match as well: it is rare for a random state.
            for (size_t i = 0; i < qubits.size(); ++i)
            {
                const uint64_t bit = uint64_t{1} << (qubits[i] % 64);
                if ((wanted >> i) & 1)
                    words[qubits[i] / 64] |= bit;
                else
                    words[qubits[i] / 64] &= ~bit;
            }
            Require(selector.Matches(words.data()) && selector.Gather(words.data()) == wanted,
                    "A satisfied control pattern was not recognized");
        }
}

// Tables are validated and analysed once per gate, and shared by the gates
// derived from it; the derived gates must act as the operations they stand for.
inline void SharedTables()
{
    using Sparse = PathIntegral::Detail::SparseGate;
    const std::vector<size_t> controls{0, 1, 2, 3, 4, 5};
    // More angles than the prefix that usually decides the fan-out: this
    // identity table needs the values derived from all of it.
    const auto identity = Gate::MultiplexedRotation(Gates::RotationAxis::Y, std::vector<double>(64, 0.), 6, controls);
    const auto *payload = identity.getCompactPayload();
    Require(payload && !payload->Info(false) && !payload->Coefficients(false),
            "A compact gate analysed its table before any backend asked");
    const Sparse still(identity, 0.);
    const auto *info = payload->Info(false);
    Require(still.structured == payload && still.diagonal && still.injective && still.maxFanout == 1 && info &&
                !identity.isBranching() && identity.getMaxBranching() == 1 && !payload->Coefficients(false),
            "An identity multiplexer was not recognized from its derived values alone");
    auto moved = identity;
    moved.setQubits({7, 6, 5, 4, 3, 2, 1});
    for (const Gate &derived : {identity.Adjoint(), identity.Transposed(), moved, Gate(identity)})
        Require(derived.getCompactPayload()->core == payload->core && derived.getCompactPayload()->Info(false) == info &&
                    Sparse(derived, 1E-15).maxFanout == 1 && payload->Info(false) == info,
                "A derived compact gate copied or analysed the shared table again");
    std::vector<double> angles(64);
    for (size_t i = 0; i < angles.size(); ++i)
        angles[i] = .11 + .07 * double(i);
    const auto generic = Gate::MultiplexedRotation(Gates::RotationAxis::Y, angles, 6, controls);
    Require(Sparse(generic, 0.).maxFanout == 2 && generic.isBranching() &&
                !generic.getCompactPayload()->Info(false) && !generic.getCompactPayload()->Coefficients(false),
            "A multiplexer whose first angles branch should not need its whole table");
    // A frontier broader than the table keeps the cosines and sines; the
    // gates derived afterwards use the same ones.
    std::vector<Gate> circuit;
    for (size_t q = 0; q < 7; ++q)
        circuit.emplace_back(Gates::HadamardGate<>().getRawOperatorMatrix(), q);
    circuit.push_back(generic);
    Simulator broad;
    broad.SetTrimValue(0.);
    broad.PropagateAll(circuit, Bits(0, 8));
    const auto *kept = generic.getCompactPayload()->Coefficients(false);
    Require(kept && kept->size() == 64 && generic.Adjoint().getCompactPayload()->Coefficients(false) == kept,
            "A broad frontier did not keep the rotation coefficients with the gate");
    // The gate, its adjoint, a remapped copy and a transposed backward half
    // against the statevector applying the corresponding StructuredGate.
    auto remapped = generic;
    remapped.setQubits({0, 6, 5, 4, 3, 2, 1});
    std::vector<std::complex<double>> phases(128);
    for (size_t i = 0; i < phases.size(); ++i)
        phases[i] = std::polar(1., .05 * double(i * i % 97));
    const auto table = Gate::Diagonal(phases, {6, 5, 4, 3, 2, 1, 0});
    auto shifted = table;
    shifted.setQubits({1, 2, 3, 4, 5, 6, 7});
    circuit.push_back(remapped);
    circuit.push_back(table.Adjoint());
    circuit.push_back(shifted);
    circuit.push_back(Gate::MultiControlled(Gates::RyGate<>(.7), 7, {0, 3, 5, 6}, 5).Adjoint());
    circuit.push_back(generic.Adjoint());
    circuit.emplace_back(Gates::HadamardGate<>().getRawOperatorMatrix(), 6);
    Require(shifted.getCompactPayload()->core == table.getCompactPayload()->core &&
                shifted.getStructuredGate() != table.getStructuredGate() &&
                shifted.getStructuredGate()->diagonal == phases && table.Adjoint().Adjoint().getStructuredGate() ==
                                                                       table.getStructuredGate(),
            "Remapping or inverting a diagonal did not share its table");
    QubitRegister<> reference(8);
    reference.ApplyGates(circuit);
    for (size_t split : {size_t{0}, size_t{1}, std::numeric_limits<size_t>::max()})
        for (bool cache : {false, true})
        {
            Simulator sim;
            sim.SetTrimValue(0.);
            sim.SetMaxDoublingsForBackwardPaths(split);
            if (!cache)
                sim.SetCompiledGateCacheMaxBytes(0);
            sim.SetCircuit(circuit);
            for (size_t output = 0; output < 256; output += 5)
                Require(std::abs(sim.Propagate(Bits(output, 8)) - reference.getBasisStateAmplitude(output)) < 4e-12,
                        "A derived compact gate disagrees with the statevector");
        }
    // A table far wider than any frontier used here: applying it must not
    // depend on its size, so the same gate is applied many times.
    const size_t wide = 18;
    std::vector<size_t> many(wide);
    for (size_t q = 0; q < wide; ++q)
        many[q] = q;
    const auto huge = Gate::Diagonal(std::vector<std::complex<double>>(size_t{1} << wide, Complex(0., 1.)), many);
    PathIntegral::AmplitudeMap single;
    single[PathIntegral::FastVectorBool(Bits(5, wide))] = 1.;
    for (unsigned repeat = 0; repeat < 4000; ++repeat)
        broad.PropagateStep(huge, single);
    Require(single.size() == 1 && single.begin()->second == Complex(1.),
            "Repeated application of a wide diagonal lost its phase");
}

// An injective gate that moves no present state keeps the amplitude index;
// one that moves some must rebuild it. Both are checked through lookups.
inline void UnmovedKeys()
{
    std::vector<Gate> prepare;
    for (size_t q = 0; q < 6; ++q)
        prepare.emplace_back(Gates::HadamardGate<>().getRawOperatorMatrix(), q);
    const std::vector<Gate> steps{
        Gate::MultiControlled(Gates::PauliXGate<>(), 7, {6, 8, 9, 0}),               // control 6 is never set
        Gate(Gates::CNOTGate<>().getRawOperatorMatrix(), 7, 6),                       // dense, control never set
        Gate(Gates::ToffoliGate<>().getRawOperatorMatrix(), 8, 0, 6),                 // dense, one control never set
        Gate::MultiControlled(Gates::PauliXGate<>(), 6, {0, 1, 2, 3}),               // moves a sixteenth of the states
        Gate(Gates::CNOTGate<>().getRawOperatorMatrix(), 7, 6),                       // now moves those
        Gate::MultiControlled(Gates::PauliYGate<>(), 9, {7, 6, 5, 0}, 3),            // moves some, with a phase
        Gate(Gates::ToffoliGate<>().getRawOperatorMatrix(), 8, 9, 7),
    };
    Simulator sim;
    sim.SetTrimValue(0.);
    sim.PropagateAll(prepare, Bits(0, 10));
    QubitRegister<> reference(10);
    reference.ApplyGates(prepare);
    for (const auto &gate : steps)
    {
        sim.PropagateStep(gate, sim.GetAmplitudes());
        reference.ApplyGate(gate);
        size_t found = 0;
        for (size_t output = 0; output < 1024; ++output)
        {
            const auto entry = sim.GetAmplitudes().find(PathIntegral::FastVectorBool(Bits(output, 10)));
            const Complex amplitude = entry == sim.GetAmplitudes().end() ? Complex(0.) : entry->second;
            found += entry != sim.GetAmplitudes().end();
            Require(std::abs(amplitude - reference.getBasisStateAmplitude(output)) < 1e-14,
                    "A permutation step left a stale amplitude index");
        }
        Require(found == 64 && sim.GetAmplitudes().size() == 64, "A permutation step changed the support size");
    }
}

// The same misuse is reported with the same exception for a compact gate and
// for a dense one, and the single-index qubit setters behave alike.
inline void MisuseAndRemap()
{
    for (bool compactGate : {false, true})
    {
        const auto make = [&](size_t target, size_t control) {
            return compactGate ? Gate::MultiControlled(Gates::PauliXGate<>(), target, {control})
                               : Gate(Gates::CNOTGate<>().getRawOperatorMatrix(), target, control);
        };
        Simulator sim;
        PathIntegral::AmplitudeMap narrow;
        narrow[PathIntegral::FastVectorBool(Bits(0, 3))] = 1.;
        Require(Throws<std::out_of_range>([&] { sim.PropagateStep(make(5, 1), narrow); }) && narrow.size() == 1,
                "A gate qubit beyond the register must be an out_of_range error");
        Require(Throws<std::out_of_range>([&] { sim.PropagateAll({make(5, 1)}, Bits(0, 3)); }),
                "A gate qubit beyond the start state must be an out_of_range error");
        for (size_t large : {size_t{1024}, size_t{5000}})
        {
            Require(Throws<std::out_of_range>([&] { sim.SetCircuit({make(large, 1)}); }),
                    "A gate qubit beyond the capacity must be an out_of_range error when a circuit is set");
            Require(Throws<std::length_error>([&] { sim.PropagateAll({make(large, 1)}); }),
                    "A circuit wider than the capacity must fail like the construction of its start state");
        }
        // Target and control exchanged one index at a time.
        Gate swapped = make(0, 1);
        swapped.setQubit1(1);
        swapped.setQubit2(0);
        swapped.setQubit3(9);
        Require(swapped.getQubit1() == 1 && swapped.getQubit2() == 0 && swapped.getQubit3() == 9 &&
                    swapped.getQubitsNumber() == 2,
                "The single-index setters did not exchange target and control");
        PathIntegral::AmplitudeMap state;
        state[PathIntegral::FastVectorBool(Bits(1, 2))] = 1.;
        sim.PropagateStep(swapped, state);
        Require(state.size() == 1 && state.find(PathIntegral::FastVectorBool(Bits(3, 2))) != state.end(),
                "A gate remapped index by index acts on the wrong qubits");
    }
    // A compact gate left with a repeated qubit is rejected where it is used.
    Gate repeated = Gate::MultiControlled(Gates::PauliXGate<>(), 0, {1, 2, 3});
    repeated.setQubit1(2);
    Simulator sim;
    PathIntegral::AmplitudeMap state;
    state[PathIntegral::FastVectorBool(Bits(0, 6))] = 1.;
    Require(Throws<std::invalid_argument>([&] { sim.PropagateStep(repeated, state); }) &&
                Throws<std::invalid_argument>([&] { sim.SetCircuit({repeated}); }) &&
                Throws<std::invalid_argument>([&] { repeated.setQubits({1, 1, 2, 3}); }),
            "A compact gate with a repeated qubit was accepted");
    repeated.setQubit2(5);
    repeated.setQubit3(0);
    sim.PropagateStep(repeated, state);
    Require(repeated.getQubits() == std::vector<size_t>({2, 5, 0, 3}) && state.size() == 1,
            "A completed index-by-index remap was not accepted");
}

inline void Run()
{
    const double pi = std::acos(-1.);
    for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
        for (double epsilon : {0., 1E-15})
        {
            const auto monomial = Gate::MultiplexedRotation(axis, {0., pi, -pi, 2. * pi}, 0, {2, 1});
            const PathIntegral::Detail::SparseGate transitions(monomial, epsilon);
            Require(transitions.injective && transitions.maxFanout == 1,
                    "Exact multiplexed quarter turns did not retain one-path metadata");
            const auto identity = Gate::MultiplexedRotation(axis, {0., 0.}, 0, {1});
            Require(PathIntegral::Detail::SparseGate(identity, epsilon).diagonal,
                    "Identity multiplexed rotation lost diagonal metadata");
            Simulator exact;
            exact.SetTrimValue(epsilon);
            exact.PropagateAll({monomial}, Bits(4, 3));
            Require(exact.GetAmplitudes().size() == 1, "Exact multiplexed quarter turn grew the path frontier");
        }
    for (double angle : {std::nextafter(pi, 0.), std::nextafter(pi, 2. * pi), 1E-200})
    {
        const auto gate = Gate::MultiplexedRotation(Gates::RotationAxis::Y, {angle}, 0, {});
        const PathIntegral::Detail::SparseGate transitions(gate, 0.);
        Require(!transitions.injective && transitions.maxFanout == 2,
                "A nearby or tiny multiplexed angle was snapped to a permutation");
        Simulator exact;
        exact.SetTrimValue(0.);
        exact.PropagateAll({gate}, Bits(0, 1));
        Require(exact.GetAmplitudes().size() == 2, "Exact path mode dropped a nonzero rotation branch");
    }
    Eigen::MatrixXcd merging(2, 2);
    merging << 1., 1., 0., 0.;
    const auto compact = Gate::MultiControlled(merging, 0, {1});
    Require(compact.isBranching() && compact.getMaxBranching() == 2,
            "Compact controlled backward branching used columns instead of rows");
    for (size_t split : {size_t{0}, size_t{1}, std::numeric_limits<size_t>::max()})
    {
        Simulator filtered;
        filtered.SetTrimValue(0.);
        filtered.SetMaxDoublingsForBackwardPaths(split);
        filtered.SetCircuit({compact});
        Require(filtered.Propagate(Bits(2, 2), Bits(2, 2)) == Complex(1.) &&
                    filtered.Propagate(Bits(3, 2), Bits(2, 2)) == Complex(1.),
                "Nonunitary compact control disagrees between forward and backward propagation");
    }
    DenseComparisons();
    WideAndTiny();
    SymmetricPayloads();
    LoweredTables();
    WordSelectors();
    SharedTables();
    UnmovedKeys();
    MisuseAndRemap();
    std::cout << "Structured path gates, cache/split modes, wide registers and tiny coefficients passed.\n";
}
} // namespace PathStructuredRegression
} // namespace QC
