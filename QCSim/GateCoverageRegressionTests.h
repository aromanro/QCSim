#pragma once

#include "Clifford.h"
#include "ExtendedStabilizer.h"
#include "PauliPropagator.h"
#include "QuantumGate.h"
#include "QubitRegisterCalculator.h"
#include <iostream>

namespace QC
{
namespace GateCoverageRegression
{
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline Vector Apply(const Vector &state, const Matrix &matrix, const std::vector<unsigned> &qubits)
{
    Vector output = Vector::Zero(state.size());
    for (unsigned col = 0; col < unsigned(state.size()); ++col)
    {
        unsigned input = 0, base = col;
        for (unsigned bit = 0; bit < qubits.size(); ++bit)
        {
            input |= ((col >> qubits[bit]) & 1U) << bit;
            base &= ~(1U << qubits[bit]);
        }
        for (unsigned row = 0; row < unsigned(matrix.rows()); ++row)
        {
            unsigned dest = base;
            for (unsigned bit = 0; bit < qubits.size(); ++bit)
                dest |= ((row >> bit) & 1U) << qubits[bit];
            output[dest] += matrix(row, input) * state[col];
        }
    }
    return output;
}

inline void RealDenseStatevectorGates()
{
    // A nonsymmetric, dense real orthogonal matrix.
    const auto mixingGate = [](unsigned dim) {
        const Eigen::VectorXd axis = Eigen::VectorXd::LinSpaced(dim, 1., dim).normalized();
        Matrix gate = (Eigen::MatrixXd::Identity(dim, dim) - 2. * axis * axis.transpose()).cast<std::complex<double>>();
        gate.row(0).swap(gate.row(1));
        return gate;
    };
    for (unsigned n : {5U, 14U})
    {
        Vector initial(1U << n);
        for (Eigen::Index i = 0; i < initial.size(); ++i)
            initial(i) = {std::sin(.13 * i + .4), std::cos(.37 * i - .1)};
        initial.normalize();
        for (const std::vector<unsigned> &qubits : std::vector<std::vector<unsigned>>{
                 {1, 3}, {3, 1}, {1, 2, n - 1}, {n - 1, 2, 1}, {n - 1, 0, 2}})
        {
            const unsigned arity = static_cast<unsigned>(qubits.size()), dim = 1U << arity;
            std::array<size_t, 3> bits{};
            for (unsigned i = 0; i < arity; ++i)
                bits[i] = size_t{1} << qubits[i];
            Matrix gate = mixingGate(dim);
            const auto check = [&](const Vector &state, const Matrix &matrix) {
                const Vector expected = Apply(state, matrix, qubits);
                const auto structure = Gates::ClassifyGateMatrix(matrix);
                for (bool parallel : {false, true})
                {
                    Vector actual = state;
                    QubitRegisterCalculator<>::ApplyGateInPlace(actual, matrix, structure, bits, arity, state.size(),
                                                               parallel);
                    Require((actual - expected).norm() < 2e-15, "Real dense gate changed its amplitudes or target order");
                    // Compare the imaginary part separately: tiny imaginary coefficients must survive dispatch.
                    const double scale = expected.imag().cwiseAbs().maxCoeff();
                    Require((actual.imag() - expected.imag()).cwiseAbs().maxCoeff() <= 2e-14 * scale,
                            "Dense gate lost an imaginary contribution");
                }
            };
            check(initial, gate);
            if (arity == 3)
            {
                Matrix controlled = Matrix::Identity(8, 8);
                controlled.bottomRightCorner(4, 4) = mixingGate(4);
                check(initial, controlled);
            }
            gate(dim - 1, dim - 1).imag(1e-180);
            check(initial.real().cast<std::complex<double>>(), gate);
        }
    }
}

inline double Expectation(const Vector &state, const std::string &pauli)
{
    std::complex<double> result = 0.;
    for (unsigned col = 0; col < unsigned(state.size()); ++col)
    {
        unsigned row = col;
        std::complex<double> phase = 1.;
        for (unsigned q = 0; q < pauli.size(); ++q)
        {
            const bool bit = (col >> q) & 1U;
            if (pauli[q] == 'X' || pauli[q] == 'Y')
                row ^= 1U << q;
            if (pauli[q] == 'Y')
                phase *= std::complex<double>(0., bit ? -1. : 1.);
            if (pauli[q] == 'Z' && bit)
                phase = -phase;
        }
        result += std::conj(state[row]) * phase * state[col];
    }
    return result.real();
}

template <class Simulator> void CheckObservables(const Vector &reference, const Simulator &sim)
{
    const char letters[] = {'I', 'X', 'Y', 'Z'};
    for (unsigned code = 0; code < 64; ++code)
    {
        std::string pauli(3, 'I');
        for (unsigned q = 0; q < 3; ++q)
            pauli[q] = letters[(code >> (2 * q)) & 3U];
        Require(std::abs(Expectation(reference, pauli) - sim.ExpectationValue(pauli)) < 2e-11,
                "gate backend disagrees with independent dense reference");
    }
}

inline Vector Prepare(ExtendedStabilizer &extended, PauliPropagator &pauli)
{
    Vector state = Vector::Zero(8);
    state[0] = 1.;
    for (unsigned q = 0; q < 3; ++q)
    {
        const Gates::UGate<> gate(.37 + q * .21, -.29 + q * .13, .19 - q * .17);
        extended.ApplyGate(gate, q);
        pauli.ApplyGate(gate, q);
        state = Apply(state, gate.getRawOperatorMatrix(), {q});
    }
    return state;
}

template <class Gate> void CheckGeneric(const Gate &gate)
{
    const auto &matrix = gate.getRawOperatorMatrix();
    Require((matrix.adjoint() * matrix - Matrix::Identity(matrix.rows(), matrix.cols())).norm() < 2e-14,
            "new gate is not unitary");
    ExtendedStabilizer extended(3);
    PauliPropagator pauli;
    pauli.SetNrQubits(3);
    pauli.SetCoefficientThreshold(0.);
    Vector state = Prepare(extended, pauli);
    std::vector<unsigned> qubits{2};
    if (gate.getQubitsNumber() > 1)
        qubits.push_back(0);
    if (gate.getQubitsNumber() > 2)
        qubits.push_back(1);
    state = Apply(state, matrix, qubits);
    extended.ApplyGate(gate, 2, 0, 1);
    pauli.ApplyGate(gate, 2, 0, 1);
    CheckObservables(state, extended);
    CheckObservables(state, pauli);
    PauliPropagator imported;
    imported.SetNrQubits(3);
    imported.SetCoefficientThreshold(0.);
    imported.SetOperations(pauli.GetOperations());
    CheckObservables(state, imported);
}

template <class Single> void CheckControls(const Single &single)
{
    const Gates::ControlledGate<> controlled(single);
    const Gates::DoublyControlledGate<> doubly(single);
    Require((controlled.getRawOperatorMatrix().bottomRightCorner(2, 2) - single.getRawOperatorMatrix()).norm() == 0.,
            "controlled gate lost the single-qubit phase");
    Require((doubly.getRawOperatorMatrix().bottomRightCorner(2, 2) - single.getRawOperatorMatrix()).norm() == 0.,
            "doubly controlled gate lost the single-qubit phase");
    CheckGeneric(controlled);
    CheckGeneric(doubly);
    ExtendedStabilizer extended(3);
    PauliPropagator pauli;
    pauli.SetNrQubits(3);
    pauli.SetCoefficientThreshold(0.);
    Vector state = Prepare(extended, pauli);
    extended.ApplyControlledGate(single, 2, 0);
    pauli.ApplyControlledGate(single, 0, 2);
    state = Apply(state, controlled.getRawOperatorMatrix(), {2, 0});
    extended.ApplyDoublyControlledGate(single, 2, 0, 1);
    pauli.ApplyDoublyControlledGate(single, 0, 1, 2);
    state = Apply(state, doubly.getRawOperatorMatrix(), {2, 0, 1});
    CheckObservables(state, extended);
    CheckObservables(state, pauli);
}

template <class Gate, class ExtAction, class PauliAction>
void CheckNative(const Gate &gate, ExtAction extAction, PauliAction pauliAction)
{
    ExtendedStabilizer extended(3);
    PauliPropagator pauli;
    pauli.SetNrQubits(3);
    pauli.SetCoefficientThreshold(0.);
    Vector state = Prepare(extended, pauli);
    state = Apply(state, gate.getRawOperatorMatrix(), {2, 0});
    extAction(extended);
    pauliAction(pauli);
    CheckObservables(state, extended);
    CheckObservables(state, pauli);
}

inline void Native()
{
    for (double angle : {.413, -1.271, 1.57079632679489661923})
    {
#define QC_CHECK_PAIR(NAME)                                                                                            \
    CheckNative(                                                                                                       \
        Gates::NAME##Gate<>(angle), [=](auto &s) { s.Apply##NAME(2, 0, angle); },                                      \
        [=](auto &s) { s.Apply##NAME(2, 0, angle); })
        QC_CHECK_PAIR(RXX);
        QC_CHECK_PAIR(RYY);
        QC_CHECK_PAIR(RZZ);
        QC_CHECK_PAIR(RZX);
        QC_CHECK_PAIR(RXY);
#undef QC_CHECK_PAIR
        for (double beta : {0., .317, -1.123})
        {
            CheckNative(
                Gates::XXPlusYYGate<>(angle, beta), [=](auto &s) { s.ApplyXXPlusYY(2, 0, angle, beta); },
                [=](auto &s) { s.ApplyXXPlusYY(2, 0, angle, beta); });
            CheckNative(
                Gates::XXMinusYYGate<>(angle, beta), [=](auto &s) { s.ApplyXXMinusYY(2, 0, angle, beta); },
                [=](auto &s) { s.ApplyXXMinusYY(2, 0, angle, beta); });
        }
    }
    CheckNative(Gates::ECRGate<>(), [](auto &s) { s.ApplyECR(2, 0); }, [](auto &s) { s.ApplyECR(2, 0); });
    CheckNative(Gates::DCXGate<>(), [](auto &s) { s.ApplyDCX(2, 0); }, [](auto &s) { s.ApplyDCX(2, 0); });
    Require(PauliDetail::ECRTransfer()->clifford && PauliDetail::DCXTransfer()->clifford,
            "named Clifford gate unnecessarily branches");
    const auto tiny = PauliDetail::CompilePairRotation('X', 'Y', 1e-17);
    Require(tiny->maxOutputs == 2 && !tiny->clifford, "tiny rotation was rounded to a Clifford");
}

inline void Definitions()
{
    const std::complex<double> i(0., 1.);
    Matrix expected(4, 4);
    expected << 0., 1., 0., i, 1., 0., -i, 0., 0., i, 0., 1., -i, 0., 1., 0.;
    Require((Gates::ECRGate<>().getRawOperatorMatrix() - expected / std::sqrt(2.)).norm() < 1e-15,
            "ECR convention mismatch");
    expected << 1., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 1., 0.;
    Require((Gates::DCXGate<>().getRawOperatorMatrix() - expected).norm() == 0., "DCX direction mismatch");
    const double angle = .729, beta = -.331, c = std::cos(angle / 2), s = std::sin(angle / 2);
    expected << c, 0., -i * s, 0., 0., c, 0., i * s, -i * s, 0., c, 0., 0., i * s, 0., c;
    Require((Gates::RZXGate<>(angle).getRawOperatorMatrix() - expected).norm() == 0., "RZX ordering mismatch");
    expected = Matrix::Identity(4, 4);
    expected(1, 1) = expected(2, 2) = c;
    expected(1, 2) = -i * s * std::exp(-i * beta);
    expected(2, 1) = -i * s * std::exp(i * beta);
    Require((Gates::XXPlusYYGate<>(angle, beta).getRawOperatorMatrix() - expected).norm() < 1e-15,
            "XX+YY beta phase mismatch");
    expected = Matrix::Identity(4, 4);
    expected(0, 0) = expected(3, 3) = c;
    expected(0, 3) = -i * s * std::exp(-i * beta);
    expected(3, 0) = -i * s * std::exp(i * beta);
    Require((Gates::XXMinusYYGate<>(angle, beta).getRawOperatorMatrix() - expected).norm() < 1e-15,
            "XX-YY beta phase mismatch");

    ExtendedStabilizer extended(3);
    PauliPropagator pauli;
    pauli.SetNrQubits(3);
    Vector state = Prepare(extended, pauli);
    Gates::AppliedGate<> applied(Gates::CCUGate<>(angle, beta, .171, -.413).getRawOperatorMatrix(), 2, 0, 1);
    extended.ApplyGate(applied);
    pauli.ApplyGate(applied);
    state = Apply(state, applied.getRawOperatorMatrix(), {2, 0, 1});
    CheckObservables(state, extended);
    CheckObservables(state, pauli);
    Gates::QuantumGateWithOp<> &base = applied;
    base.setOperator(Gates::CCYGate<>().getRawOperatorMatrix());
    extended.ApplyGate(applied);
    pauli.ApplyGate(applied);
    state = Apply(state, applied.getRawOperatorMatrix(), {2, 0, 1});
    CheckObservables(state, extended);
    CheckObservables(state, pauli);
}

template <class Gate, class Action> void CheckClifford(const Gate &gate, Action action)
{
    for (unsigned basis = 0; basis < 8; ++basis)
    {
        Clifford::StabilizerSimulator sim(3);
        Vector state = Vector::Zero(8);
        state[basis] = 1.;
        for (unsigned q = 0; q < 3; ++q)
            if ((basis >> q) & 1U)
                sim.ApplyX(q);
        sim.ApplyH(1);
        sim.ApplyS(1);
        sim.ApplyCX(0, 1);
        state = Apply(state, Gates::HadamardGate<>().getRawOperatorMatrix(), {1});
        state = Apply(state, Gates::SGate<>().getRawOperatorMatrix(), {1});
        state = Apply(state, Gates::CNOTGate<>().getRawOperatorMatrix(), {0, 1});
        state = Apply(state, gate.getRawOperatorMatrix(), {2, 0});
        action(sim);
        CheckObservables(state, sim);
    }
}

inline bool Run()
{
    try
    {
        Definitions();
        RealDenseStatevectorGates();
        CheckGeneric(Gates::RXXGate<>(.337));
        CheckGeneric(Gates::RYYGate<>(-.733));
        CheckGeneric(Gates::RZZGate<>(.173));
        CheckGeneric(Gates::RZXGate<>(-.937));
        CheckGeneric(Gates::RXYGate<>(.727));
        CheckGeneric(Gates::XXPlusYYGate<>(.173, .733));
        CheckGeneric(Gates::XXMinusYYGate<>(.473, -.333));
        CheckGeneric(Gates::ECRGate<>());
        CheckGeneric(Gates::DCXGate<>());
        CheckGeneric(Gates::CCUGate<>(.731, -.253, .837, -.471));
        CheckControls(Gates::HadamardGate<>());
        CheckControls(Gates::HyGate<>());
        CheckControls(Gates::SGate<>());
        CheckControls(Gates::SDGGate<>());
        CheckControls(Gates::TGate<>());
        CheckControls(Gates::TDGGate<>());
        CheckControls(Gates::PhaseShiftGate<>(.333));
        CheckControls(Gates::PauliXGate<>());
        CheckControls(Gates::PauliYGate<>());
        CheckControls(Gates::PauliZGate<>());
        CheckControls(Gates::SquareRootNOTGate<>());
        CheckControls(Gates::SquareRootNOTDagGate<>());
        CheckControls(Gates::SplitterGate<>());
        CheckControls(Gates::RxGate<>(.391));
        CheckControls(Gates::RyGate<>(-.717));
        CheckControls(Gates::RzGate<>(.919));
        CheckControls(Gates::UGate<>(.331, -.119, .731));
        Native();
        CheckClifford(Gates::ECRGate<>(), [](auto &s) { s.ApplyECR(2, 0); });
        CheckClifford(Gates::DCXGate<>(), [](auto &s) { s.ApplyDCX(2, 0); });
        const double halfPi = std::acos(-1.) * .5;
        for (double angle : {-4. * halfPi, -2. * halfPi, 0., 2. * halfPi, 4. * halfPi})
            for (double beta : {-halfPi, 0., halfPi})
            {
                CheckClifford(Gates::XXPlusYYGate<>(angle, beta), [=](auto &s) { s.ApplyXXPlusYY(2, 0, angle, beta); });
                CheckClifford(Gates::XXMinusYYGate<>(angle, beta),
                              [=](auto &s) { s.ApplyXXMinusYY(2, 0, angle, beta); });
            }
#define QC_CHECK_CLIFFORD_PAIR(NAME)                                                                                   \
    CheckClifford(Gates::NAME##Gate<>(halfPi), [=](auto &s) { s.Apply##NAME(2, 0, halfPi); })
        QC_CHECK_CLIFFORD_PAIR(RXX);
        QC_CHECK_CLIFFORD_PAIR(RYY);
        QC_CHECK_CLIFFORD_PAIR(RZZ);
        QC_CHECK_CLIFFORD_PAIR(RZX);
        QC_CHECK_CLIFFORD_PAIR(RXY);
#undef QC_CHECK_CLIFFORD_PAIR
        Clifford::StabilizerSimulator invalid(3);
        bool rejected = false;
        try
        {
            invalid.ApplyRXY(0, 1, .123);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        Require(rejected && invalid.ExpectationValue("ZZZ") == 1., "non-Clifford angle changed the tableau");
        std::cout << "Expanded gate coverage regressions passed.\n";
        return true;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Expanded gate regression: " << error.what() << '\n';
        return false;
    }
}

} // namespace GateCoverageRegression
} // namespace QC
