#pragma once

#include "DensityMatrix.h"

#include <iostream>

namespace QC
{
namespace StructuredGateTests
{

inline Eigen::VectorXcd State(size_t dimension, double offset)
{
    Eigen::VectorXcd result(dimension);
    for (size_t i = 0; i < dimension; ++i)
        result(i) = {std::sin(0.37 * i + offset), std::cos(0.19 * i - offset)};
    result.normalize();
    return result;
}

// Deliberately dense oracle, assembling the action on each basis vector without
// using the calculator's masks, block traversal, or structured execution code.
inline Eigen::MatrixXcd Dense(const Gates::StructuredGate &gate, size_t nrQubits)
{
    const size_t dimension = size_t{1} << nrQubits;
    Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(dimension, dimension);
    for (size_t input = 0; input < dimension; ++input)
    {
        size_t selector = 0;
        for (size_t k = 0; k < gate.qubits.size(); ++k)
            if ((input / (size_t{1} << gate.qubits[k])) % 2)
                selector += size_t{1} << k;
        if (gate.kind == Gates::StructuredGate::Kind::Diagonal)
        {
            result(input, input) = gate.diagonal[selector];
            continue;
        }
        Eigen::Matrix2cd local;
        if (gate.kind == Gates::StructuredGate::Kind::Controlled)
        {
            const size_t selected = gate.controlState == std::numeric_limits<size_t>::max()
                                        ? (size_t{1} << gate.qubits.size()) - 1
                                        : gate.controlState;
            if (selector != selected)
            {
                result(input, input) = 1.;
                continue;
            }
            local << gate.matrix[0], gate.matrix[1], gate.matrix[2], gate.matrix[3];
        }
        else if (gate.axis == Gates::RotationAxis::X)
            local = Gates::RxGate<>(gate.angles[selector]).getRawOperatorMatrix();
        else if (gate.axis == Gates::RotationAxis::Y)
            local = Gates::RyGate<>(gate.angles[selector]).getRawOperatorMatrix();
        else
            local = Gates::RzGate<>(gate.angles[selector]).getRawOperatorMatrix();
        const size_t bit = size_t{1} << gate.target;
        const size_t col = (input / bit) % 2;
        result(input & ~bit, input) = local(0, col);
        result(input | bit, input) = local(1, col);
    }
    return result;
}

inline std::vector<Gates::StructuredGate> Operations(size_t nrQubits)
{
    std::vector<Gates::StructuredGate> result;
    Gates::StructuredGate diagonal;
    diagonal.qubits = {nrQubits - 1, 0, 3, 1};
    for (size_t changed : {size_t{0}, size_t{1}, size_t{8}, size_t{9}})
    {
        diagonal.diagonal.assign(16, 1.);
        for (size_t i = 0; i < changed; ++i)
            diagonal.diagonal[i] = std::polar(1., .31 + .17 * i);
        result.push_back(diagonal);
    }
    diagonal.diagonal.clear();
    for (size_t i = 0; i < 16; ++i)
        diagonal.diagonal.push_back(std::polar(1., 0.23 * i - 0.7));
    result.push_back(diagonal);
    diagonal.qubits.clear();
    diagonal.diagonal.clear();
    for (size_t q = nrQubits; q-- > 0;)
        diagonal.qubits.push_back(q);
    for (size_t i = 0; i < (size_t{1} << nrQubits); ++i)
        diagonal.diagonal.push_back(std::polar(1., 0.13 * i + 0.2));
    result.push_back(diagonal);
    diagonal.qubits.clear();
    diagonal.diagonal = {std::polar(1., 0.42)};
    result.push_back(diagonal);
    std::vector<Eigen::MatrixXcd> matrices{
        Gates::PauliXGate<>().getRawOperatorMatrix(), Gates::PauliYGate<>().getRawOperatorMatrix(),
        Gates::PauliZGate<>().getRawOperatorMatrix(), Gates::PhaseShiftGate<>(0.37).getRawOperatorMatrix(),
        Gates::RxGate<>(0.38).getRawOperatorMatrix(), Gates::RyGate<>(-0.82).getRawOperatorMatrix(),
        Gates::RzGate<>(1.29).getRawOperatorMatrix(), Gates::UGate<>(0.3, -0.7, 0.9, 0.27).getRawOperatorMatrix()};
    for (size_t target : {size_t{0}, size_t{2}, nrQubits - 1})
    {
        for (size_t m = 0; m < matrices.size(); ++m)
        {
            Gates::StructuredGate gate;
            gate.kind = Gates::StructuredGate::Kind::Controlled;
            gate.target = target;
            for (size_t q = nrQubits; q-- > 0;)
                if (q != target && gate.qubits.size() < 5)
                    gate.qubits.push_back(q);
            gate.controlState = m % 2 ? 10 : std::numeric_limits<size_t>::max();
            const auto &matrix = matrices[m];
            gate.matrix = {matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1)};
            result.push_back(gate);
            if (m == 0)
            {
                gate.qubits.clear();
                gate.controlState = 0;
                result.push_back(gate);
            }
        }
        for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
        {
            Gates::StructuredGate gate;
            gate.kind = Gates::StructuredGate::Kind::MultiplexedRotation;
            gate.target = target;
            gate.axis = axis;
            for (size_t q = nrQubits; q-- > 0;)
                if (q != target)
                    gate.qubits.push_back(q);
            for (size_t i = 0; i < (size_t{1} << gate.qubits.size()); ++i)
                gate.angles.push_back(0.37 * i - 1.2);
            result.push_back(gate);
            gate.qubits.clear();
            gate.angles.resize(1);
            result.push_back(gate);
        }
    }
    return result;
}

template <class Simulator> void ApplyPublic(Simulator &sim, const Gates::StructuredGate &gate)
{
    if (gate.kind == Gates::StructuredGate::Kind::Diagonal)
        sim.ApplyDiagonalGate(gate.diagonal, gate.qubits);
    else if (gate.kind == Gates::StructuredGate::Kind::Controlled)
    {
        Eigen::MatrixXcd matrix(2, 2);
        matrix << gate.matrix[0], gate.matrix[1], gate.matrix[2], gate.matrix[3];
        sim.ApplyMultiControlledGate(matrix, gate.target, gate.qubits, gate.controlState);
    }
    else if (gate.axis == Gates::RotationAxis::X)
        sim.ApplyMultiplexedRX(gate.angles, gate.target, gate.qubits);
    else if (gate.axis == Gates::RotationAxis::Y)
        sim.ApplyMultiplexedRY(gate.angles, gate.target, gate.qubits);
    else
        sim.ApplyMultiplexedRZ(gate.angles, gate.target, gate.qubits);
}

template <class Simulator> bool InvalidInputs(Simulator &sim)
{
    const auto unchanged = sim.Clone();
    const auto density = unchanged->getDensityMatrix();
    const auto rejects = [&](const auto &action) {
        try
        {
            action();
        }
        catch (const std::invalid_argument &)
        {
            return (sim.getDensityMatrix() - density).norm() == 0.;
        }
        return false;
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!rejects([&] { sim.ApplyDiagonalGate({1., 1.}, {0, 1}); }) ||
        !rejects([&] { sim.ApplyDiagonalGate({1., 1., 1., 1.}, {0, 0}); }) ||
        !rejects([&] { sim.ApplyDiagonalGate({1., nan}, {0}); }) || !rejects([&] { sim.ApplyMCX(0, {0}); }) ||
        !rejects([&] { sim.ApplyMCX(0, {1, 1}); }) || !rejects([&] { sim.ApplyMCX(sim.getNrQubits(), {}); }) ||
        !rejects([&] { sim.ApplyMCX(0, {sim.getNrQubits()}); }) || !rejects([&] { sim.ApplyMCX(0, {1}, 2); }) ||
        !rejects([&] { sim.ApplyMCRY(nan, 0, {1}); }) || !rejects([&] { sim.ApplyMultiplexedRX({0.2}, 0, {1}); }) ||
        !rejects([&] { sim.ApplyMultiplexedRY({0., nan}, 0, {1}); }) ||
        !rejects([&] { sim.ApplyMultiControlledGate(Eigen::MatrixXcd::Identity(4, 4), 0, {1}); }))
        return false;
    return true;
}

struct RestoreParallelSettings
{
    size_t sv = QubitRegister<>::GetParallelMinBasisStates();
    size_t dm = DensityMatrix<>::GetParallelMinElements();
    int threads = QubitRegister<>::GetParallelThreads();

    RestoreParallelSettings()
    {
        QubitRegister<>::SetParallelMinBasisStates(1);
        DensityMatrix<>::SetParallelMinElements(1);
        QubitRegister<>::SetParallelThreads(2);
    }

    ~RestoreParallelSettings()
    {
        QubitRegister<>::SetParallelMinBasisStates(sv);
        DensityMatrix<>::SetParallelMinElements(dm);
        QubitRegister<>::SetParallelThreads(threads);
    }
};

inline bool Traversal(bool parallel)
{
    // Exercise contiguous two-pair SIMD and its bit-zero/strided fallbacks.
    QubitRegister<> wide(9, 19);
    wide.SetMultithreading(parallel);
    for (size_t target : {size_t{0}, size_t{1}, size_t{6}})
        for (bool controlZero : {false, true})
            for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
            {
                if (target == 0 && controlZero)
                    continue;
                Gates::StructuredGate operation;
                operation.kind = Gates::StructuredGate::Kind::MultiplexedRotation;
                operation.axis = axis;
                operation.target = target;
                operation.qubits = {8, 3, controlZero ? size_t{0} : size_t{5}};
                operation.angles = {.13, -.27, .39, -.51, .73, -.91, 1.17, -1.33};
                const auto input = State(512, .43);
                wide.setRegisterStorage(input);
                ApplyPublic(wide, operation);
                if ((wide.getRegisterStorage() - Dense(operation, 9) * input).norm() > 3E-12)
                    return false;
            }
    for (size_t controls : {size_t{2}, size_t{7}})
    {
        Gates::StructuredGate operation;
        operation.kind = Gates::StructuredGate::Kind::Controlled;
        operation.target = 6;
        operation.qubits = {8, 1, 5, 2, 7, 4, 3};
        operation.qubits.resize(controls);
        operation.controlState = controls == 2 ? 1 : 43;
        const auto matrix = Gates::UGate<>(.37, -.71, .29).getRawOperatorMatrix();
        operation.matrix = {matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1)};
        const auto input = State(512, .63);
        wide.setRegisterStorage(input);
        ApplyPublic(wide, operation);
        if ((wide.getRegisterStorage() - Dense(operation, 9) * input).norm() > 3E-12)
            return false;
    }
    // Cross tile boundaries with all qubits in an intentionally permuted table.
    QubitRegister<> sim(12, 41);
    sim.SetMultithreading(parallel);
    Eigen::VectorXcd expected = State(4096, 0.71);
    sim.setRegisterStorage(expected);
    const std::vector<size_t> qubits{11, 2, 9, 0, 4, 8, 3, 10, 6, 1, 7, 5};
    const auto select = [](size_t index, const std::vector<size_t> &bits) {
        size_t result = 0;
        for (size_t k = 0; k < bits.size(); ++k)
            result |= ((index >> bits[k]) & 1) << k;
        return result;
    };
    std::vector<std::complex<double>> values(4096);
    for (size_t i = 0; i < values.size(); ++i)
        values[i] = std::polar(1., 0.003 * i);
    for (size_t i = 0; i < 4096; ++i)
        expected(i) *= values[select(i, qubits)];
    sim.ApplyDiagonalGate(values, qubits);
    if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
        return false;
    for (size_t target : {size_t{0}, size_t{5}, size_t{11}})
    {
        std::vector<size_t> controls;
        for (size_t q : qubits)
            if (q != target)
                controls.push_back(q);
        std::vector<double> angles(2048);
        for (size_t i = 0; i < angles.size(); ++i)
            angles[i] = -1.4 + 0.0021 * i;
        const size_t bit = size_t{1} << target;
        for (size_t low = 0; low < 4096; ++low)
            if (!(low & bit))
            {
                const double half = angles[select(low, controls)] * 0.5;
                const auto a = expected(low), b = expected(low | bit);
                expected(low) = std::cos(half) * a - std::sin(half) * b;
                expected(low | bit) = std::sin(half) * a + std::cos(half) * b;
            }
        sim.ApplyMultiplexedRY(angles, target, controls);
        if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
            return false;
    }
    const std::vector<size_t> controls{11, 9, 6, 4, 1};
    for (size_t low = 0; low < 4096; ++low)
        if (!(low & 8) && select(low, controls) == 21)
            std::swap(expected(low), expected(low | 8));
    sim.ApplyMCX(3, controls, 21);
    if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
        return false;

    Eigen::VectorXcd storage = Eigen::VectorXcd::Constant(128, {42., -7.});
    using Strided = Eigen::Map<Eigen::VectorXcd, 0, Eigen::InnerStride<2>>;
    Strided strided(storage.data(), 64);
    Eigen::VectorXcd reference = State(64, 0.9);
    strided = reference;
    for (const auto &gate : Operations(6))
    {
        reference = (Dense(gate, 6) * reference).eval();
        gate.Validate(6);
        QubitRegisterCalculator<Strided>::ApplyStructuredGateInPlace(strided, gate, 64, parallel);
        if ((strided - reference).norm() > 3E-12)
            return false;
        if (gate.kind == Gates::StructuredGate::Kind::MultiplexedRotation)
        {
            const auto coefficients = QubitRegisterCalculator<Strided>::PrepareStructuredRotation(gate);
            reference = (Dense(gate, 6).conjugate() * reference).eval();
            QubitRegisterCalculator<Strided>::ApplyStructuredGateInPlace(strided, gate, 64, parallel, 0, true,
                                                                         &coefficients);
            if ((strided - reference).norm() > 3E-12)
                return false;
        }
    }
    for (size_t i = 1; i < 128; i += 2)
        if (storage(i) != std::complex<double>(42., -7.))
            return false;
    return true;
}

template <class Matrix> bool Density(bool parallel)
{
    DensityMatrix<Eigen::VectorXcd, Matrix> sim(6, 31);
    sim.SetMultithreading(parallel);
    const auto first = State(64, 0.1), second = State(64, 0.7);
    Eigen::MatrixXcd expected = 0.37 * first * first.adjoint() + 0.63 * second * second.adjoint();
    sim.setDensityMatrix(expected);
    for (const auto &gate : Operations(6))
    {
        const auto dense = Dense(gate, 6);
        expected = (dense * expected * dense.adjoint()).eval();
        ApplyPublic(sim, gate);
        if ((sim.getDensityMatrix() - expected).norm() > 3E-12)
            return false;
    }
    sim.SaveState();
    sim.ApplyMCU(0.3, -0.4, 0.7, 2, {5, 0, 1, 3, 4}, 0.19);
    sim.RestoreState();
    return (sim.getDensityMatrix() - expected).norm() < 3E-12 && InvalidInputs(sim);
}

// Circuit storage must use the same compact representation as direct execution.
inline bool AppliedRecording()
{
    using Applied = Gates::AppliedGate<>;

    Eigen::MatrixXcd merging(2, 2);
    merging << 1., 1., 0., 0.;
    const auto compactMerge = Applied::MultiControlled(merging, 0, {1});
    if (!compactMerge.isBranching() || compactMerge.getMaxBranching() != 2)
        return false;
    if (!Gates::CSGate<>().isDiagonal() || !Gates::CCRzGate<>(.3).isDiagonal() || !Gates::CCYGate<>().isAntidiagonal())
        return false;
    Gates::ControlledGate<> changing{Gates::SGate<>()};
    changing.SetOperation(Gates::PauliXGate<>().getRawOperatorMatrix());
    if (changing.isDiagonal() || !changing.isAntidiagonal())
        return false;
    const double pi = std::acos(-1.);

    struct Capture : StructuredGateOperations<Capture, Eigen::MatrixXcd>
    {
        Gates::StructuredGate operation;

        void ApplyStructuredGate(const Gates::StructuredGate &gate)
        {
            gate.Validate(3);
            operation = gate;
        }
    } capture;

    capture.ApplyMCRX(pi, 0, {2, 1});
    if (capture.operation.matrix[0] != std::complex<double>(0.) || Applied(capture.operation).isBranching())
        return false;
    capture.ApplyMCRX(std::nextafter(pi, 0.), 0, {2, 1});
    if (capture.operation.matrix[0] == std::complex<double>(0.) || !Applied(capture.operation).isBranching())
        return false;
    capture.ApplyMCRY(pi, 0, {2, 1});
    if (capture.operation.matrix != std::array<std::complex<double>, 4>{0., -1., 1., 0.})
        return false;
    capture.ApplyMCRZ(pi, 0, {2, 1});
    if (capture.operation.matrix[0] != std::complex<double>(0., -1.) ||
        capture.operation.matrix[3] != std::complex<double>(0., 1.))
        return false;
    capture.ApplyMCP(pi, 0, {2, 1});
    if (capture.operation.matrix[3] != std::complex<double>(-1.))
        return false;
    const Gates::RxGate<> legacyRotation(pi);
    capture.ApplyMultiControlledGate(legacyRotation, 0, {2, 1});
    if (capture.operation.matrix[0] != legacyRotation.getRawOperatorMatrix()(0, 0))
        return false;
    for (double sign : {-1., 1.})
    {
        const double angle = sign * pi;
        if (Gates::RXXGate<>(angle).getRawOperatorMatrix()(0, 0) != std::complex<double>(0.))
            return false;
        const auto swap = Gates::XXPlusYYGate<>(angle).getRawOperatorMatrix();
        if (swap(1, 1) != std::complex<double>(0.))
            return false;
        for (double direction : {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()})
            if (Gates::RXXGate<>(std::nextafter(angle, direction)).getRawOperatorMatrix()(0, 0) ==
                std::complex<double>(0.))
                return false;
        QubitRegister<> exact(1);
        exact.ApplyMultiplexedRY({angle}, 0, {});
        if (exact.getBasisStateAmplitude(0) != std::complex<double>(0.) ||
            exact.getBasisStateAmplitude(1) != std::complex<double>(sign))
            return false;
    }

    struct RecordedRegister : QubitRegister<>
    {
        using QubitRegister<>::QubitRegister;

        size_t Count() const
        {
            return computeGates.size();
        }
    };

    std::vector<Applied> circuit;
    Eigen::VectorXcd initial = State(128, 0.41), expected = initial;
    for (const auto &operation : Operations(7))
    {
        Applied gate(operation);
        if (!gate.getStructuredGate() || gate.getQubits().size() != gate.getQubitsNumber() ||
            gate.getMaxBranching() > 2)
            return false;
        Applied copy = gate;
        Applied moved = std::move(copy);
        copy = moved;
        moved = Applied();
        moved = std::move(copy);
        circuit.push_back(moved);
        expected = (Dense(operation, 7) * expected).eval();
        circuit.emplace_back(Gates::HadamardGate<>().getRawOperatorMatrix(), 3);
        expected = (Gates::HadamardGate<>().getOperatorMatrix(7, 3) * expected).eval();
    }
    RecordedRegister sim(7, 37);
    sim.setRegisterStorage(initial);
    sim.ComputeStart();
    sim.ApplyGates(circuit);
    sim.ComputeEnd();
    if (sim.Count() != circuit.size() || (sim.getRegisterStorage() - expected).norm() > 4E-12)
        return false;
    const auto clone = sim.Clone();
    clone->Uncompute();
    sim.Uncompute();
    if ((sim.getRegisterStorage() - initial).norm() > 4E-12 || (clone->getRegisterStorage() - initial).norm() > 4E-12)
        return false;
    sim.Compute();
    if ((sim.getRegisterStorage() - expected).norm() > 4E-12)
        return false;
    sim.ComputeClear();
    for (auto gate = circuit.rbegin(); gate != circuit.rend(); ++gate)
        sim.ApplyGate(gate->Adjoint());
    if ((sim.getRegisterStorage() - initial).norm() > 4E-12)
        return false;
    const auto expectation = sim.ExpectationValue(circuit);
    if (std::abs(expectation - initial.dot(expected)) > 4E-12)
        return false;
    sim.ComputeStart();
    sim.ApplyGate(Gates::PauliXGate<>(), 0);
    const Eigen::VectorXcd recordedState = sim.getRegisterStorage();
    (void)sim.ExpectationValue(circuit);
    if (sim.Count() != 1 || (sim.getRegisterStorage() - recordedState).norm() != 0.)
        return false;
    auto invalidCircuit = circuit;
    invalidCircuit.push_back(Applied::MultiControlled(Gates::PauliXGate<>(), 100, {}));
    try
    {
        (void)sim.ExpectationValue(invalidCircuit);
        return false;
    }
    catch (const std::invalid_argument &)
    {
    }
    if (sim.Count() != 1 || (sim.getRegisterStorage() - recordedState).norm() != 0.)
        return false;
    sim.ApplyGate(Gates::PauliXGate<>(), 0);
    sim.ComputeEnd();
    if (sim.Count() != 2)
        return false;
    sim.ComputeClear();

    auto control = Applied::MultiControlled(Gates::RyGate<>(0.7), 6, {4, 2, 0, 5}, 9);
    const auto original = control;
    control.setQubits({1, 2, 3, 4, 5});
    control.setQubit1(0);
    if (original.getQubit1() != 6 || control.getQubit1() != 0 || control.getQubits()[4] != 5)
        return false;
    sim.ApplyGate(control);
    auto operation = *control.getStructuredGate();
    if ((sim.getRegisterStorage() - Dense(operation, 7) * initial).norm() > 4E-12)
        return false;
    Gates::QuantumGateWithOp<> &base = control;
    try
    {
        (void)base.getRawOperatorMatrix();
        return false;
    }
    catch (const std::invalid_argument &)
    {
    }
    base.setOperator(Gates::PauliZGate<>().getRawOperatorMatrix());
    if (control.getStructuredGate() || control.getQubitsNumber() != 1 ||
        control.getStructure().kind != Gates::GateStructure::Kind::Diagonal)
        return false;

    struct MutableApplied : Applied
    {
        using Applied::Applied;

        void Replace(const Eigen::MatrixXcd &matrix)
        {
            operatorMat = matrix;
        }

        void OperatorChanged() override
        {
        }
    };

    for (bool direct : {false, true})
    {
        MutableApplied mutableGate(*original.getStructuredGate());
        if (direct)
            mutableGate.Replace(Gates::PauliXGate<>().getRawOperatorMatrix());
        else
            mutableGate.setOperator(Gates::PauliXGate<>().getRawOperatorMatrix());
        Applied copied(mutableGate);
        Applied moved(std::move(mutableGate));
        for (const auto *gate : {&copied, &moved})
            if (gate->getStructuredGate() || gate->getQubitsNumber() != 1 ||
                gate->getStructure().kind != Gates::GateStructure::Kind::Antidiagonal)
                return false;
    }

    std::vector<size_t> controls;
    for (size_t q = 0; q < 63; ++q)
        controls.push_back(q + 100);
    auto wide = Applied::MultiControlled(Gates::PauliXGate<>(), 70, controls);
    if (wide.getQubitsNumber() != 64 || wide.getMaxBranching() != 1 || wide.isBranching())
        return false;
    const auto before = sim.getRegisterStorage();
    try
    {
        sim.ApplyGate(wide);
        return false;
    }
    catch (const std::invalid_argument &)
    {
    }
    if ((sim.getRegisterStorage() - before).norm() != 0.)
        return false;
    auto diagonal = Applied::Diagonal({1., std::complex<double>(0., 1.)}, {3});
    try
    {
        Gates::AppliedGate<Eigen::Matrix2cd> fixed(*diagonal.getStructuredGate());
        return false;
    }
    catch (const std::invalid_argument &)
    {
    }
    auto mux = Applied::MultiplexedRotation(Gates::RotationAxis::X, {1E-12, 0.7}, 4, {2});
    Eigen::Index next = -1;
    if (diagonal.isBranching() || !mux.isBranchingInState(0, next) || next != 0)
        return false;
    try
    {
        auto bad = Applied::MultiControlled(Gates::PauliXGate<>(), 0, {0});
        return false;
    }
    catch (const std::invalid_argument &)
    {
    }

    // Remapping the qubits, the adjoint and the transpose share the tables of
    // the gate they come from. A backend that reads a StructuredGate gets the
    // remapped or conjugated copy, built once.
    {
        std::vector<std::complex<double>> phases(64);
        for (size_t i = 0; i < phases.size(); ++i)
            phases[i] = std::polar(1., .1 * double(i));
        const auto table = Applied::Diagonal(phases, {0, 1, 2, 3, 4, 5});
        const Gates::StructuredGate *tables = table.getStructuredGate();
        auto remapped = table;
        remapped.setQubits({6, 5, 4, 3, 2, 1});
        const auto inverse = table.Adjoint();
        if (remapped.getCompactPayload()->core != table.getCompactPayload()->core ||
            inverse.getCompactPayload()->core != table.getCompactPayload()->core ||
            &table.getCompactPayload()->core->operation != tables || inverse.Adjoint().getStructuredGate() != tables ||
            table.Transposed().getStructuredGate() != tables)
            return false;
        const Gates::StructuredGate *copy = remapped.getStructuredGate();
        if (copy == tables || copy != remapped.getStructuredGate() || copy != Applied(remapped).getStructuredGate() ||
            copy->qubits != std::vector<size_t>({6, 5, 4, 3, 2, 1}) || copy->diagonal != phases ||
            inverse.getStructuredGate()->diagonal[3] != std::conj(phases[3]) ||
            inverse.getStructuredGate()->qubits != tables->qubits)
            return false;
        const auto rotation = Applied::MultiplexedRotation(Gates::RotationAxis::Y, {.3, -.7}, 1, {0});
        if (rotation.Adjoint().getStructuredGate()->angles != std::vector<double>({-.3, .7}) ||
            rotation.Transposed().getStructuredGate()->angles != std::vector<double>({-.3, .7}) ||
            rotation.Transposed().Transposed().getStructuredGate() != rotation.getStructuredGate())
            return false;
        QubitRegister<> direct(7, 11), shared(7, 11);
        const auto start = State(128, 0.37);
        direct.setRegisterStorage(start);
        shared.setRegisterStorage(start);
        direct.ApplyStructuredGate(*copy);
        direct.ApplyStructuredGate(tables->Adjoint());
        shared.ApplyGate(remapped);
        shared.ApplyGate(inverse);
        if ((direct.getRegisterStorage() - shared.getRegisterStorage()).norm() != 0. ||
            (direct.getRegisterStorage() - Dense(tables->Adjoint(), 7) * (Dense(*copy, 7) * start)).norm() > 4E-12)
            return false;
    }

    // The single-index setters behave as for a dense gate: target and control
    // can be exchanged one index at a time, an index beyond the arity is kept
    // and ignored, and nothing is copied. A repeated qubit is rejected when the
    // gate is applied, and at once by a complete list.
    {
        auto exchanged = Applied::MultiControlled(Gates::PauliXGate<>(), 0, {1});
        const auto core = exchanged.getCompactPayload()->core;
        Applied legacy(Gates::CNOTGate<>().getRawOperatorMatrix(), 0, 1);
        for (Applied *gate : {&exchanged, &legacy})
        {
            gate->setQubit1(1);
            gate->setQubit2(0);
            gate->setQubit3(5);
            if (gate->getQubit1() != 1 || gate->getQubit2() != 0 || gate->getQubit3() != 5 ||
                gate->getQubitsNumber() != 2 || gate->getQubits() != std::vector<size_t>({1, 0}))
                return false;
        }
        if (exchanged.getCompactPayload()->core != core)
            return false;
        QubitRegister<> compactRegister(3, 5), denseRegister(3, 5);
        const auto start = State(8, 0.61);
        compactRegister.setRegisterStorage(start);
        denseRegister.setRegisterStorage(start);
        compactRegister.ApplyGate(exchanged);
        denseRegister.ApplyGate(legacy);
        if ((compactRegister.getRegisterStorage() - denseRegister.getRegisterStorage()).norm() != 0.)
            return false;
        auto repeated = Applied::MultiControlled(Gates::PauliXGate<>(), 0, {1, 2});
        repeated.setQubit1(2);
        const auto before = compactRegister.getRegisterStorage();
        try
        {
            compactRegister.ApplyGate(repeated);
            return false;
        }
        catch (const std::invalid_argument &)
        {
        }
        try
        {
            repeated.setQubits({2, 1, 2});
            return false;
        }
        catch (const std::invalid_argument &)
        {
        }
        if ((compactRegister.getRegisterStorage() - before).norm() != 0. || repeated.getQubit1() != 2)
            return false;
        repeated.setQubit3(0);
        compactRegister.ApplyGate(repeated);
        if (repeated.getQubits() != std::vector<size_t>({2, 1, 0}) ||
            (compactRegister.getRegisterStorage() - Dense(*repeated.getStructuredGate(), 3) * before).norm() > 4E-12)
            return false;
    }

    // Full-register replay on two/three qubits must retain all physical indices.
    for (size_t n : {size_t(2), size_t(3)})
    {
        QubitRegister<> shortRegister(n, 17);
        const auto first = State(size_t{1} << n, 0.2);
        shortRegister.setRegisterStorage(first);
        const auto full = Gates::CNOTGate<>().getOperatorMatrix(n, n - 1, 0);
        shortRegister.ComputeStart();
        shortRegister.ApplyOperatorMatrix(full);
        shortRegister.ApplyMCX(0, {1});
        shortRegister.ComputeEnd();
        shortRegister.Uncompute();
        if ((shortRegister.getRegisterStorage() - first).norm() > 2E-12)
            return false;
    }
    return true;
}

template <class Matrix> bool AppliedDensity()
{
    using Applied = Gates::AppliedGate<Matrix>;
    DensityMatrix<Eigen::VectorXcd, Matrix> sim(6, 31);
    const auto first = State(64, 0.2), second = State(64, 0.8);
    Eigen::MatrixXcd expected = 0.29 * first * first.adjoint() + 0.71 * second * second.adjoint();
    const Eigen::MatrixXcd initial = expected;
    sim.setDensityMatrix(expected);
    std::vector<Applied> circuit;
    for (const auto &operation : Operations(6))
    {
        circuit.emplace_back(operation);
        const auto full = Dense(operation, 6);
        expected = (full * expected * full.adjoint()).eval();
    }
    sim.ApplyGates(circuit);
    if ((sim.getDensityMatrix() - expected).norm() > 4E-12)
        return false;
    for (auto gate = circuit.rbegin(); gate != circuit.rend(); ++gate)
        sim.ApplyGate(gate->Adjoint());
    if ((sim.getDensityMatrix() - initial).norm() >= 4E-12)
        return false;
    DensityMatrix<Eigen::VectorXcd, Matrix> scaled(1, 7);
    Matrix tiny = Matrix::Zero(2, 2);
    tiny(0, 0) = 1E-300;
    tiny(1, 1) = 1.;
    scaled.setDensityMatrix(tiny);
    scaled.ApplyGate(Applied::Diagonal({1E200, 0.}, {0}));
    if (!scaled.getDensityMatrix().allFinite() || std::abs(scaled.getDensityMatrix()(0, 0) / 1E100 - 1.) > 2E-14)
        return false;
    scaled.ApplyDiagonalGate({1E-200, 0.}, {0});
    return scaled.getDensityMatrix().allFinite() && std::abs(scaled.getDensityMatrix()(0, 0) / 1E-300 - 1.) < 2E-14;
}

// The entry points that take tables from their caller apply them where they are;
// an applied gate is not validated entry by entry again. Both must give the bits
// of the StructuredGate entry point, keep every rejection, and leave a recording
// that owns its tables.
inline bool BorrowedTables()
{
    using Complex = std::complex<double>;
    const size_t n = 7, dimension = size_t{1} << n;
    for (bool parallel : {false, true})
    {
        // Diagonals in the three traversals: runs of equal factors (table qubits on
        // top), pairs (a table qubit at bit zero) with factors that are exactly one,
        // one with a negative zero, and a full permuted table.
        std::vector<Gates::StructuredGate> gates;
        const auto diagonal = [&](std::vector<size_t> qubits, double offset, bool withOnes) {
            Gates::StructuredGate gate;
            gate.qubits = std::move(qubits);
            gate.diagonal.resize(size_t{1} << gate.qubits.size());
            for (size_t i = 0; i < gate.diagonal.size(); ++i)
                gate.diagonal[i] = std::polar(1. + 0.01 * double(i % 5), offset + 0.37 * double(i));
            if (withOnes)
                for (size_t i = 0; i < gate.diagonal.size(); ++i)
                {
                    if (i % 4 < 2)
                        gate.diagonal[i] = 1.;
                    if (i % 8 == 6)
                        gate.diagonal[i] = Complex(1., -0.);
                }
            gates.push_back(std::move(gate));
        };
        diagonal({6, 5, 4}, 0.2, false);
        diagonal({0, 3}, 0.4, false);
        diagonal({3, 0}, 0.5, true);
        diagonal({0, 1, 2, 3, 4, 5, 6}, 0.6, true);
        diagonal({6, 4, 2, 0, 1, 3, 5}, 0.7, false);
        diagonal({2, 0, 6}, 0.8, true);
        for (const auto &gate : gates)
        {
            const Eigen::VectorXcd initial = State(dimension, 0.33);
            const auto apply = [&](const auto &action) {
                QubitRegister<> sim(n, 5);
                sim.SetMultithreading(parallel);
                sim.setRegisterStorage(initial);
                action(sim);
                return Eigen::VectorXcd(sim.getRegisterStorage());
            };
            const Eigen::VectorXcd structured = apply([&](auto &sim) { sim.ApplyStructuredGate(gate); });
            const Eigen::VectorXcd direct =
                apply([&](auto &sim) { sim.ApplyDiagonalGate(gate.diagonal, gate.qubits); });
            const Eigen::VectorXcd viewed =
                apply([&](auto &sim) { sim.ApplyStructuredGateView(Gates::StructuredGateView(gate)); });
            const Eigen::VectorXcd applied = apply([&](auto &sim) { sim.ApplyGate(Gates::AppliedGate<>(gate)); });
            if (std::memcmp(structured.data(), direct.data(), sizeof(Complex) * dimension) != 0 ||
                std::memcmp(structured.data(), viewed.data(), sizeof(Complex) * dimension) != 0 ||
                std::memcmp(structured.data(), applied.data(), sizeof(Complex) * dimension) != 0)
                return false;
            if ((structured - Dense(gate, n) * initial).norm() > 3E-13)
                return false;
        }

        // A recording made through a borrowed table owns a copy of it.
        {
            QubitRegister<> sim(n, 5);
            sim.SetMultithreading(parallel);
            const Eigen::VectorXcd initial = State(dimension, 0.41);
            sim.setRegisterStorage(initial);
            std::vector<Complex> values(8);
            for (size_t i = 0; i < values.size(); ++i)
                values[i] = std::polar(1., 0.3 + double(i));
            std::vector<double> angles{0.3, -0.7, 1.1, 0.2};
            std::vector<size_t> qubits{1, 4, 6}, controls{5, 2};
            sim.ComputeStart();
            sim.ApplyDiagonalGate(values, qubits);
            sim.ApplyMultiplexedRY(angles, 0, controls);
            sim.ApplyMCRX(0.4, 3, controls);
            sim.ComputeEnd();
            std::fill(values.begin(), values.end(), Complex(7., 7.));
            std::fill(angles.begin(), angles.end(), 9.);
            qubits = {0, 1, 2};
            controls = {6, 1};
            sim.Uncompute();
            if ((sim.getRegisterStorage() - initial).norm() > 3E-12)
                return false;
        }

        // An applied gate whose qubits were changed one at a time is checked where it
        // is used; its tables are not read for that.
        {
            Gates::StructuredGate gate;
            gate.qubits = {0, 1, 2};
            gate.diagonal.assign(8, Complex(0., 1.));
            const Eigen::VectorXcd initial = State(dimension, 0.52);
            const auto rejected = [&](const Gates::AppliedGate<> &bad) {
                QubitRegister<> sim(n, 5);
                sim.SetMultithreading(parallel);
                sim.setRegisterStorage(initial);
                DensityMatrix<> density(n, 5);
                density.SetMultithreading(parallel);
                const Eigen::MatrixXcd rho = initial * initial.adjoint();
                density.setDensityMatrix(rho);
                bool first = false, second = false;
                try
                {
                    sim.ApplyGate(bad);
                }
                catch (const std::invalid_argument &)
                {
                    first =
                        std::memcmp(sim.getRegisterStorage().data(), initial.data(), sizeof(Complex) * dimension) == 0;
                }
                try
                {
                    density.ApplyGate(bad);
                }
                catch (const std::invalid_argument &)
                {
                    second = (density.getDensityMatrix() - rho).norm() == 0.;
                }
                return first && second;
            };
            Gates::AppliedGate<> repeated(gate);
            repeated.setQubit1(2);
            Gates::AppliedGate<> outside(gate);
            outside.setQubit2(n);
            if (!rejected(repeated) || !rejected(outside))
                return false;
            // The same gate moved onto valid qubits is applied, and equals the dense action.
            Gates::AppliedGate<> moved(gate);
            moved.setQubit1(6);
            moved.setQubit3(4);
            QubitRegister<> sim(n, 5);
            sim.SetMultithreading(parallel);
            sim.setRegisterStorage(initial);
            sim.ApplyGate(moved);
            Gates::StructuredGate reference = gate;
            reference.qubits = {6, 1, 4};
            if ((sim.getRegisterStorage() - Dense(reference, n) * initial).norm() > 3E-13)
                return false;
        }

        // The density matrix through the same three entry points.
        {
            Gates::StructuredGate gate;
            gate.kind = Gates::StructuredGate::Kind::MultiplexedRotation;
            gate.axis = Gates::RotationAxis::X;
            gate.target = 2;
            gate.qubits = {4, 0};
            gate.angles = {0.3, -1.2, 0.8, 2.1};
            const auto first = State(32, 0.15), second = State(32, 0.65);
            const Eigen::MatrixXcd rho = 0.4 * first * first.adjoint() + 0.6 * second * second.adjoint();
            const auto apply = [&](const auto &action) {
                DensityMatrix<> sim(5, 5);
                sim.SetMultithreading(parallel);
                sim.setDensityMatrix(rho);
                action(sim);
                return Eigen::MatrixXcd(sim.getDensityMatrix());
            };
            const Eigen::MatrixXcd structured = apply([&](auto &sim) { sim.ApplyStructuredGate(gate); });
            const Eigen::MatrixXcd direct =
                apply([&](auto &sim) { sim.ApplyMultiplexedRX(gate.angles, gate.target, gate.qubits); });
            const Eigen::MatrixXcd applied = apply([&](auto &sim) { sim.ApplyGate(Gates::AppliedGate<>(gate)); });
            const auto dense = Dense(gate, 5);
            if ((structured - direct).norm() != 0. || (structured - applied).norm() != 0. ||
                (structured - dense * rho * dense.adjoint()).norm() > 3E-13)
                return false;
        }
    }
    return true;
}

// An applied compact gate is executed from its shared tables: a remapped,
// adjoint or transposed gate copies none of them, and the cosines and sines of
// a multiplexer are computed once and kept. Each such gate must act as the
// dense matrix of its own descriptor, on the statevector and on the density
// matrix, the first time and every time after.
inline bool SharedTableKernels()
{
    using Complex = std::complex<double>;
    using Applied = Gates::AppliedGate<>;
    const size_t n = 6, dimension = size_t{1} << n;
    Eigen::MatrixXcd skew(2, 2);
    skew << Complex(.6, .1), Complex(.2, -.7), Complex(-.3, .4), Complex(.5, .5);
    std::vector<Complex> diagonal(8);
    for (size_t i = 0; i < diagonal.size(); ++i)
        diagonal[i] = std::polar(1. + 0.05 * double(i), 0.4 * double(i) - 1.);
    std::vector<Applied> gates{Applied::Diagonal(diagonal, {4, 0, 2}), Applied::MultiControlled(skew, 3, {5, 1}, 2)};
    for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
        gates.push_back(Applied::MultiplexedRotation(axis, {0.3, -1.2, 0.8, 2.1}, 1, {4, 2}));
    for (bool parallel : {false, true})
        for (const Applied &source : gates)
        {
            Applied remapped = source;
            auto qubits = remapped.getQubits();
            std::rotate(qubits.begin(), qubits.begin() + 1, qubits.end());
            for (size_t &q : qubits)
                q = (q + 1) % n;
            remapped.setQubits(qubits);
            const std::vector<Applied> variants{source,
                                                source.Adjoint(),
                                                source.Transposed(),
                                                remapped,
                                                remapped.Adjoint(),
                                                remapped.Adjoint().Adjoint(),
                                                remapped.Transposed().Adjoint()};
            const Eigen::VectorXcd initial = State(dimension, 0.27);
            const Eigen::MatrixXcd rho =
                0.35 * initial * initial.adjoint() + 0.65 * State(dimension, 0.81) * State(dimension, 0.81).adjoint();
            for (const Applied &gate : variants)
            {
                const Eigen::MatrixXcd dense = Dense(*gate.getStructuredGate(), n);
                QubitRegister<> sim(n, 5);
                sim.SetMultithreading(parallel);
                sim.setRegisterStorage(initial);
                const Eigen::VectorXcd start = sim.getRegisterStorage();
                DensityMatrix<> density(n, 5);
                density.SetMultithreading(parallel);
                density.setDensityMatrix(rho);
                Eigen::VectorXcd expected = start;
                Eigen::MatrixXcd expectedRho = rho;
                // The second and third applications run on what the first one left with the tables.
                for (int repetition = 0; repetition < 3; ++repetition)
                {
                    sim.ApplyGate(gate);
                    density.ApplyGate(gate);
                    expected = (dense * expected).eval();
                    expectedRho = (dense * expectedRho * dense.adjoint()).eval();
                    if ((sim.getRegisterStorage() - expected).norm() > 1E-12 * (1. + expected.norm()) ||
                        (density.getDensityMatrix() - expectedRho).norm() > 1E-12 * (1. + expectedRho.norm()))
                        return false;
                }
                // Recorded and undone: the replay uses the adjoint gate over the same tables.
                QubitRegister<> recorded(n, 5);
                recorded.SetMultithreading(parallel);
                recorded.setRegisterStorage(initial);
                recorded.ComputeStart();
                recorded.ApplyGate(gate);
                recorded.ComputeEnd();
                if ((recorded.getRegisterStorage() - dense * start).norm() > 1E-12 * (1. + (dense * start).norm()))
                    return false;
                // Uncompute applies the adjoint, which undoes a unitary gate only.
                if ((dense * dense.adjoint() - Eigen::MatrixXcd::Identity(dimension, dimension)).norm() < 1E-9)
                {
                    recorded.Uncompute();
                    if ((recorded.getRegisterStorage() - start).norm() > 1E-12)
                        return false;
                }
            }
        }
    return true;
}

// A large table that is not in state order goes through a buffer (blocked
// traversal). Every amplitude must get the factor of its own table entry,
// whatever the permutation, for tables that leave qubits out, in one thread and
// in several; an amplitude whose factor is exactly one keeps its bits.
inline bool BlockedDiagonals()
{
    using Complex = std::complex<double>;
    const auto check = [](size_t n, const std::vector<size_t> &qubits, bool parallel, bool withOnes) {
        const size_t dimension = size_t{1} << n;
        std::vector<Complex> table(size_t{1} << qubits.size());
        for (size_t i = 0; i < table.size(); ++i)
            table[i] = std::polar(1. + 1e-3 * double(i % 7), 0.37 * double(i % 1021) + 0.1);
        if (withOnes)
            for (size_t i = 0; i < table.size(); ++i)
            {
                if (i % 5 < 2)
                    table[i] = 1.;
                if (i % 16 == 9)
                    table[i] = Complex(1., -0.);
            }
        Eigen::VectorXcd initial(dimension);
        for (size_t i = 0; i < dimension; ++i)
            initial[i] = Complex(0.3 + 1e-6 * double(i % 4093), -0.2 + 1e-6 * double(i % 2039));
        QubitRegister<> sim(n, 5);
        sim.SetMultithreading(parallel);
        sim.setRegisterStorage(initial);
        // setRegisterStorage normalizes: read back what the register holds.
        const Eigen::VectorXcd before = sim.getRegisterStorage();
        sim.ApplyDiagonalGate(table, qubits);
        const Eigen::VectorXcd &after = sim.getRegisterStorage();
        double worst = 0.;
        for (size_t i = 0; i < dimension; ++i)
        {
            size_t index = 0;
            for (size_t k = 0; k < qubits.size(); ++k)
                index |= ((i >> qubits[k]) & 1) << k;
            const Complex factor = table[index];
            if (factor == Complex(1.))
            {
                if (std::memcmp(&after[i], &before[i], sizeof(Complex)) != 0)
                    return false;
                continue;
            }
            worst = std::max(worst, std::abs(after[i] - before[i] * factor) / std::abs(before[i]));
        }
        return worst < 1E-14;
    };
    // The lowest state bits select the highest table bits in all of these.
    const auto reversed = [](size_t n) {
        std::vector<size_t> qubits(n);
        for (size_t i = 0; i < n; ++i)
            qubits[i] = n - 1 - i;
        return qubits;
    };
    // The blocked traversal starts at 2^18 table entries in one thread and at 2^20
    // when the threads share the work; the tables below that leave qubits out are
    // still that large.
    for (bool parallel : {false, true})
    {
        const size_t n = parallel ? 21 : 20;
        std::vector<size_t> rotated(n), mixed = reversed(n);
        for (size_t i = 0; i < n; ++i)
            rotated[i] = (i + n / 2) % n;
        // Reversed, with the order inside the two halves of the table changed.
        std::swap(mixed[0], mixed[3]);
        std::swap(mixed[1], mixed[6]);
        std::swap(mixed[n - 1], mixed[n - 4]);
        std::swap(mixed[n - 2], mixed[n - 7]);
        // Tables that leave out a low state bit, a middle one, or two bits.
        std::vector<size_t> withoutLow = reversed(n), withoutMiddle = reversed(n), withoutTwo = reversed(n);
        withoutLow.erase(std::find(withoutLow.begin(), withoutLow.end(), size_t{2}));
        withoutMiddle.erase(std::find(withoutMiddle.begin(), withoutMiddle.end(), n / 2));
        withoutTwo.erase(std::find(withoutTwo.begin(), withoutTwo.end(), n - 1));
        withoutTwo.erase(std::find(withoutTwo.begin(), withoutTwo.end(), size_t{5}));
        if (!check(n, reversed(n), parallel, false) || !check(n, reversed(n), parallel, true) ||
            !check(n, rotated, parallel, false) || !check(n, mixed, parallel, true) ||
            !check(n, withoutLow, parallel, false) || !check(n, withoutMiddle, parallel, true) ||
            (!parallel && !check(n, withoutTwo, parallel, false)))
            return false;
    }
    return true;
}

} // namespace StructuredGateTests
} // namespace QC

inline bool StatevectorStructuredGateTests()
{
    using namespace QC;
    StructuredGateTests::RestoreParallelSettings restore;
    if (!StructuredGateTests::AppliedRecording())
        return false;
    if (!StructuredGateTests::BorrowedTables())
        return false;
    if (!StructuredGateTests::BlockedDiagonals())
        return false;
    if (!StructuredGateTests::SharedTableKernels())
        return false;
    std::cout << "Structured statevector gates and replay..." << std::endl;
    for (bool parallel : {false, true})
    {
        if (!StructuredGateTests::Traversal(parallel))
            return false;
        QubitRegister<> sim(7, 37);
        sim.SetMultithreading(parallel);
        const Eigen::VectorXcd initial = StructuredGateTests::State(128, 0.21);
        sim.setRegisterStorage(initial);
        Eigen::VectorXcd expected = initial;
        sim.ComputeStart();
        for (const auto &gate : StructuredGateTests::Operations(7))
        {
            const auto dense = StructuredGateTests::Dense(gate, 7);
            expected = (dense * expected).eval();
            StructuredGateTests::ApplyPublic(sim, gate);
            if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
                return false;
            // Interleave the old recording representation with compact operations.
            sim.ApplyGate(Gates::HadamardGate<>(), 3);
            expected = Gates::HadamardGate<>().getOperatorMatrix(7, 3) * expected;
        }
        sim.ComputeEnd();
        const auto clone = sim.Clone();
        sim.Uncompute();
        if ((sim.getRegisterStorage() - initial).norm() > 3E-12)
            return false;
        sim.Compute();
        if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
            return false;
        clone->Uncompute();
        if ((clone->getRegisterStorage() - initial).norm() > 3E-12 || !StructuredGateTests::InvalidInputs(sim))
            return false;
        sim.ComputeClear();
        sim.Uncompute();
        if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
            return false;
        sim.ComputeStart();
        sim.ApplyMCX(0, {1, 2, 3, 4, 5, 6});
        sim.ApplyMCZ(0, {1, 2, 3, 4, 5, 6});
        sim.ApplyMCP(0.27, 0, {1, 2, 3, 4, 5, 6});
        sim.ApplyMCRX(0.72, 0, {1, 2, 3, 4, 5, 6});
        sim.ApplyMCRY(-0.35, 0, {1, 2, 3, 4, 5, 6});
        sim.ApplyMCRZ(0.47, 0, {1, 2, 3, 4, 5, 6});
        sim.ApplyMCU(0.3, 0.4, 0.5, 0, {1, 2, 3, 4, 5, 6}, -0.17);
        sim.ComputeEnd();
        sim.Uncompute();
        if ((sim.getRegisterStorage() - expected).norm() > 3E-12)
            return false;
    }
    return true;
}

inline bool DensityMatrixStructuredGateTests()
{
    QC::StructuredGateTests::RestoreParallelSettings restore;
    using RowMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    std::cout << "Structured density gates, mixed states and storage orders..." << std::endl;
    return QC::StructuredGateTests::AppliedDensity<Eigen::MatrixXcd>() &&
           QC::StructuredGateTests::AppliedDensity<RowMatrix>() &&
           QC::StructuredGateTests::Density<Eigen::MatrixXcd>(false) &&
           QC::StructuredGateTests::Density<Eigen::MatrixXcd>(true) &&
           QC::StructuredGateTests::Density<RowMatrix>(false) && QC::StructuredGateTests::Density<RowMatrix>(true);
}
