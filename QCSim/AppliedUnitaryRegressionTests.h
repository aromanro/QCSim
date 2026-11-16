#pragma once

#include "PauliGatePrecisionRegressionTests.h"
#include "QuantumGate.h"
#include <iostream>
#include <thread>
#include <utility>

namespace QC
{
namespace AppliedUnitaryRegression
{
using PauliGatePrecisionRegression::CheckTransfer;
using PauliGatePrecisionRegression::Coefficient;
using PauliGatePrecisionRegression::Require;

template <class Applied> inline void Verify(const Applied &gate)
{
    const auto descriptor = Gates::GetUnitaryGateDescriptor(gate);
    const auto &matrix = gate.getRawOperatorMatrix();
    Require(descriptor && Gates::VerifyUnitaryGateDescriptor(matrix, descriptor),
            "A typed applied gate lost its verified constructor descriptor");
    const unsigned arity = unsigned(gate.getQubitsNumber());
    const auto table = PauliDetail::MatrixTransfer(matrix, arity, descriptor);
    Require((table->unchanged & 1U) && table->offsets[1] == 1,
            "A typed applied gate branches or rescales an untouched identity");
    CheckTransfer(Eigen::MatrixXcd(matrix), arity, table);
}

template <class Gate> inline void GateAndAdjoint(const Gate &source)
{
    using Applied = Gates::AppliedGate<>;
    const Applied forward(source, 2, 0, 1);
    Verify(forward);
    const auto inverse = forward.Adjoint();
    Verify(inverse);
    Require((inverse.getRawOperatorMatrix().array() == forward.getRawOperatorMatrix().adjoint().array()).all(),
            "Applied adjoint recomputed rather than transposed its stored matrix");
    const auto twice = inverse.Adjoint();
    Require(Gates::GetUnitaryGateDescriptor(twice) == Gates::GetUnitaryGateDescriptor(forward) &&
                (twice.getRawOperatorMatrix().array() == forward.getRawOperatorMatrix().array()).all(),
            "A double adjoint changed the matrix or descriptor");
    Require(inverse.getQubit1() == 2 && inverse.getQubit2() == 0 && inverse.getQubit3() == 1,
            "An adjoint changed its physical qubit mapping");
    const Applied raw(source.getRawOperatorMatrix(), 2, 0, 1);
    Require(!Gates::GetUnitaryGateDescriptor(raw) && !Gates::GetUnitaryGateDescriptor(raw.Adjoint()),
            "Raw matrix construction inferred a unitary descriptor");

    Applied copied(inverse);
    Verify(copied);
    Applied moved(std::move(copied));
    Verify(moved);
    Applied assigned;
    assigned = moved;
    Verify(assigned);
    // Assignment to itself, by copy and by move. Through a reference, so that
    // the compiler does not warn about what the test is there to check.
    Applied &itself = assigned;
    assigned = itself;
    assigned = std::move(itself);
    Verify(assigned);
    Applied moveAssigned;
    moveAssigned = std::move(assigned);
    Verify(moveAssigned);
    moveAssigned = raw;
    Require(!Gates::GetUnitaryGateDescriptor(moveAssigned), "Raw copy assignment retained old unitary provenance");
    moveAssigned = Applied(source.getRawOperatorMatrix());
    Require(!Gates::GetUnitaryGateDescriptor(moveAssigned), "Raw move assignment retained old unitary provenance");
}

inline void MutationAndStorage()
{
    using Applied = Gates::AppliedGate<>;
    Gates::RyGate<> source(.3);
    Applied snapshot(source, 0);
    source.setOperator(Eigen::MatrixXcd(2. * Eigen::MatrixXcd::Identity(2, 2)));
    Verify(snapshot);
    Applied untrusted(source, 0);
    const auto untrustedDescriptor = Gates::GetUnitaryGateDescriptor(untrusted);
    Require(untrustedDescriptor &&
                !Gates::VerifyUnitaryGateDescriptor(untrusted.getRawOperatorMatrix(), untrustedDescriptor),
            "Typed construction trusted a previously mutated source");
    const auto table = PauliDetail::MatrixTransfer(untrusted.getRawOperatorMatrix(), 1, untrustedDescriptor);
    Require(Coefficient(*table, 0, 0) == 4., "Typed applied construction unitarized a mutated source");
    const auto badInverse = untrusted.Adjoint();
    Require(!Gates::VerifyUnitaryGateDescriptor(badInverse.getRawOperatorMatrix(),
                                                Gates::GetUnitaryGateDescriptor(badInverse)),
            "Adjoint restored trust to a mutated source");

    Gates::QuantumGateWithOp<> &base = snapshot;
    base.setOperator(Eigen::MatrixXcd(2. * Eigen::MatrixXcd::Identity(2, 2)));
    Require(!Gates::GetUnitaryGateDescriptor(snapshot), "A base-reference matrix mutation retained applied provenance");
    snapshot = Applied(Gates::RxGate<>(.3));
    base = Eigen::MatrixXcd::Identity(2, 2);
    Require(!Gates::GetUnitaryGateDescriptor(snapshot), "Base matrix assignment retained applied provenance");

    struct DerivedApplied : Applied
    {
        using Applied::AppliedGate;

        explicit DerivedApplied(const Applied &source) : Applied(source)
        {
        }

        void Mutate()
        {
            this->operatorMat(0, 0) = 2.;
        }

        void ReplaceDense(const Eigen::MatrixXcd &matrix)
        {
            this->operatorMat = matrix;
        }
    };

    DerivedApplied derived(Gates::RyGate<>(.3), 0);
    derived.Mutate();
    const Applied derivedCopy(derived);
    const Applied derivedMove(std::move(derived));
    for (const auto *gate : {&derivedCopy, &derivedMove})
        Require(
            !Gates::VerifyUnitaryGateDescriptor(gate->getRawOperatorMatrix(), Gates::GetUnitaryGateDescriptor(*gate)),
            "A derived protected mutation bypassed applied-gate verification");

    using Fixed = Eigen::Matrix<std::complex<double>, 2, 2>;
    const Gates::RxGate<Fixed> fixedSource(.3);
    const Gates::AppliedGate<Fixed> fixed(fixedSource, 3);
    Verify(fixed);
    Verify(fixed.Adjoint());
    Verify(Applied(fixedSource, 3));
    using RowMajor = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const Gates::AppliedGate<RowMajor> rowMajor(Gates::CCUGate<>(.43, -.71, .19, .23), 2, 0, 1);
    Verify(rowMajor);
    Verify(rowMajor.Adjoint());
    const auto compact = Applied::MultiplexedRotation(Gates::RotationAxis::Y, {.1, .3}, 0, {1});
    Require(compact.getStructuredGate() && compact.Adjoint().getStructuredGate() &&
                !Gates::GetUnitaryGateDescriptor(compact) && !Gates::GetUnitaryGateDescriptor(compact.Adjoint()),
            "Applied unitary storage interfered with compact gate storage");
    Applied mixed = compact;
    Require(mixed.getStructuredGate() == compact.getStructuredGate(), "Compact copies stopped sharing their payload");
    mixed = Applied(Gates::RxGate<>(.3), 0);
    Verify(mixed);
    Require(!mixed.getStructuredGate(), "A typed dense assignment retained compact storage");
    mixed = compact;
    Require(!Gates::GetUnitaryGateDescriptor(mixed) && mixed.getStructuredGate() == compact.getStructuredGate(),
            "Compact assignment retained a unitary payload");
    Applied movedCompact(std::move(mixed));
    movedCompact.setQubits({2, 3});
    Require(compact.getQubit1() == 0 && compact.getQubit2() == 1 && movedCompact.getQubit1() == 2 &&
                movedCompact.getQubit2() == 3,
            "Compact remapping changed a shared source payload");
    DerivedApplied replacedCompact(compact);
    replacedCompact.ReplaceDense(Gates::RxGate<>(.3).getRawOperatorMatrix());
    Require(!replacedCompact.getStructuredGate() && !Gates::GetUnitaryGateDescriptor(replacedCompact),
            "A derived dense replacement misinterpreted a structured payload");
    DerivedApplied emptiedTyped(Gates::RxGate<>(.3));
    emptiedTyped.ReplaceDense(Eigen::MatrixXcd());
    Require(!emptiedTyped.getStructuredGate(), "An empty dense matrix misinterpreted a unitary payload as structured");

    std::vector<Applied> retained;
    for (unsigned index = 0; index < 48; ++index)
        retained.emplace_back(Gates::UGate<>(.31 + .01 * index, -.71, .19, .23));
    for (const auto &gate : retained)
        Verify(gate);
    Applied fromThread;
    std::thread worker([&] { fromThread = Applied(Gates::CCUGate<>(.43, -.71, .19, .23), 2, 0, 1); });
    worker.join();
    Verify(fromThread);
    Verify(fromThread.Adjoint());
}

inline void CircuitIdentity()
{
    using Applied = Gates::AppliedGate<>;
    std::vector<Applied> circuit;
    for (unsigned index = 0; index < 6; ++index)
    {
        const size_t q = 3 * index;
        circuit.emplace_back(Gates::CCPhaseShiftGate<>(.3), q, q + 1, q + 2);
        circuit.emplace_back(Gates::CCUGate<>(.43, -.71, .19, .23), q, q + 1, q + 2);
        circuit.emplace_back(Gates::PhaseShiftGate<>(.3), q);
        circuit.emplace_back(Gates::RXXGate<>(.3), q, q + 1);
        circuit.emplace_back(Gates::XXPlusYYGate<>(.51, -.29), q, q + 1);
    }
    PauliPropagator simulator;
    simulator.SetNrQubits(100);
    simulator.SetCoefficientThreshold(0.);
    for (const auto &gate : circuit)
        simulator.ApplyGate(gate);
    for (auto gate = circuit.rbegin(); gate != circuit.rend(); ++gate)
        simulator.ApplyGate(gate->Adjoint());
    PauliStringXZWithCoefficient observable(100);
    observable.Z[99] = true;
    PauliStringStorage output;
    Require(simulator.ExpectationValue(observable, output) == 1. && output.size() == 1,
            "A typed AppliedGate circuit or its adjoint expanded a disjoint observable");
}

// A constructor's matrix must verify against its own descriptor wherever the
// gate is built and whatever the compiler does with the arithmetic (fused
// multiply-add contraction, constant folding): the descriptor keeps the stored
// coefficients, so nothing is evaluated a second time.
template <class Gate> inline void VerifyReproducible(const Gate &gate, const char *message)
{
    const auto descriptor = Gates::GetUnitaryGateDescriptor(gate);
    Require(descriptor && descriptor.hasBuilt &&
                Gates::VerifyUnitaryGateDescriptor(gate.getRawOperatorMatrix(), descriptor),
            message);
    auto inverse = descriptor;
    inverse.adjoint = true;
    Require(Gates::VerifyUnitaryGateDescriptor(Eigen::MatrixXcd(gate.getRawOperatorMatrix().adjoint()), inverse),
            message);
}

inline void ConstructorReproducibility()
{
    // A fixed linear congruential sequence: the same angles on every platform.
    uint64_t state = 0x9E3779B97F4A7C15ULL;
    const auto next = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return (double(state >> 11) / double(uint64_t(1) << 53)) * 6.4 - 3.2;
    };
    for (unsigned sample = 0; sample < 2000; ++sample)
    {
        const double t = next(), p = next(), l = next(), g = next();
        const Gates::UGate<> u(t, p, l, g);
        VerifyReproducible(u, "A U gate does not reproduce its own constructor matrix");
        VerifyReproducible(Gates::ControlledUGate<>(t, p, l, g), "A controlled U gate lost its constructor matrix");
        VerifyReproducible(Gates::ControlledGate<>(u), "A wrapped controlled U gate lost its constructor matrix");
        VerifyReproducible(Gates::DoublyControlledGate<>(u), "A wrapped CC-U gate lost its constructor matrix");
        VerifyReproducible(Gates::CCUGate<>(t, p, l, g), "A CCU gate lost its constructor matrix");
        VerifyReproducible(Gates::PhaseShiftGate<>(t), "A phase gate lost its constructor matrix");
        VerifyReproducible(Gates::RxGate<>(t), "An Rx gate lost its constructor matrix");
        VerifyReproducible(Gates::RyGate<>(p), "An Ry gate lost its constructor matrix");
        VerifyReproducible(Gates::RzGate<>(l), "An Rz gate lost its constructor matrix");
        VerifyReproducible(Gates::ControlledRzGate<>(g), "A controlled Rz gate lost its constructor matrix");
        VerifyReproducible(Gates::CCRyGate<>(t), "A CC-Ry gate lost its constructor matrix");
        VerifyReproducible(Gates::ControlledPhaseShiftGate<>(p), "A controlled phase gate lost its constructor matrix");
        VerifyReproducible(Gates::RZXGate<>(t), "An RZX gate lost its constructor matrix");
        VerifyReproducible(Gates::RXYGate<>(l), "An RXY gate lost its constructor matrix");
        VerifyReproducible(Gates::XXPlusYYGate<>(t, p), "An XX+YY gate lost its constructor matrix");
        VerifyReproducible(Gates::XXMinusYYGate<>(l, g), "An XX-YY gate lost its constructor matrix");
    }
    VerifyReproducible(Gates::TGate<>(), "A T gate lost its constructor matrix");
    VerifyReproducible(Gates::TDGGate<>(), "A T-dagger gate lost its constructor matrix");

    // The stored coefficients are part of the candidate: a matrix that differs
    // from them in its last bit is not the constructor's matrix.
    const Gates::UGate<> u(.43, -.71, .19, .23);
    auto descriptor = Gates::GetUnitaryGateDescriptor(u);
    Eigen::MatrixXcd changed = u.getRawOperatorMatrix();
    changed(1, 0) = {std::nextafter(changed(1, 0).real(), 2.), changed(1, 0).imag()};
    Require(!Gates::VerifyUnitaryGateDescriptor(changed, descriptor),
            "A changed matrix entry passed the verification of a stored constructor");
    // A candidate without stored coefficients still verifies by evaluation.
    descriptor = {Gates::UnitaryGateDescriptor::Kind::Rx, {.3}, 0};
    Require(!descriptor.hasBuilt && Gates::UnitaryGateMatrix(descriptor).rows() == 2 &&
                std::abs(Gates::UnitaryGateMatrix(descriptor)(0, 0) - std::cos(.15)) < 1e-15,
            "A bare descriptor no longer evaluates its matrix");

    // The practical consequence: U-family gates on other qubits leave an
    // observable alone, on every compiler.
    PauliPropagator simulator;
    simulator.SetNrQubits(40);
    simulator.SetCoefficientThreshold(0.);
    for (unsigned index = 0; index < 6; ++index)
    {
        const double t = next(), p = next(), l = next(), g = next();
        const int q = int(6 * index);
        simulator.ApplyGate(Gates::UGate<>(t, p, l, g), q);
        simulator.ApplyGate(Gates::ControlledUGate<>(t, p, l, g), q + 1, q + 2);
        simulator.ApplyGate(Gates::CCUGate<>(t, p, l, g), q + 3, q + 4, q + 5);
    }
    PauliStringXZWithCoefficient observable(40);
    observable.Z[39] = true;
    PauliStringStorage output;
    Require(simulator.ExpectationValue(observable, output) == 1. && output.size() == 1,
            "Typed U-family gates on other qubits expanded an untouched observable");
}

inline void Run()
{
    ConstructorReproducibility();
    GateAndAdjoint(Gates::TGate<>());
    GateAndAdjoint(Gates::PhaseShiftGate<>(.3));
    GateAndAdjoint(Gates::UGate<>(.43, -.71, .19, .23));
    GateAndAdjoint(Gates::ControlledUGate<>(.43, -.71, .19, .23));
    GateAndAdjoint(Gates::CCPhaseShiftGate<>(.3));
    GateAndAdjoint(Gates::CCUGate<>(.43, -.71, .19, .23));
    GateAndAdjoint(Gates::RXXGate<>(.3));
    GateAndAdjoint(Gates::XXPlusYYGate<>(.51, -.29));
    GateAndAdjoint(Gates::XXMinusYYGate<>(.51, -.29));
    GateAndAdjoint(Gates::FredkinGate<>());
    MutationAndStorage();
    CircuitIdentity();
    std::cout << "Typed AppliedGate, adjoint, mutation and storage regressions passed.\n";
}
} // namespace AppliedUnitaryRegression
} // namespace QC
