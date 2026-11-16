#pragma once

#include "ExtendedStabilizer.h"
#include "PauliPropagator.h"
#include <iostream>

namespace QC
{
namespace PauliGatePrecisionRegression
{
inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline Eigen::MatrixXcd Pauli(unsigned code, unsigned qubits)
{
    const unsigned dimension = 1U << qubits;
    Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(dimension, dimension);
    for (unsigned col = 0; col < dimension; ++col)
    {
        unsigned row = col;
        std::complex<double> value = 1.;
        for (unsigned q = 0; q < qubits; ++q)
        {
            const unsigned local = (code >> (2 * q)) & 3U;
            if (local & 1U)
                row ^= 1U << q;
            if (local == 3)
                value *= std::complex<double>(0., (col & (1U << q)) ? -1. : 1.);
            else if (local == 2 && (col & (1U << q)))
                value = -value;
        }
        result(row, col) = value;
    }
    return result;
}

inline double Coefficient(const PauliDetail::LocalTransfer &table, unsigned input, unsigned output)
{
    for (unsigned i = table.offsets[input]; i < table.offsets[input + 1]; ++i)
        if (table.entries[i].pauli == output)
            return table.entries[i].coefficient;
    return 0.;
}

inline void CheckTransfer(const Eigen::MatrixXcd &matrix, unsigned qubits,
                          const std::shared_ptr<const PauliDetail::LocalTransfer> &table)
{
    const unsigned dimension = 1U << qubits;
    std::vector<Eigen::MatrixXcd> paulis;
    for (unsigned code = 0; code < dimension * dimension; ++code)
        paulis.push_back(Pauli(code, qubits));
    for (unsigned input = 0; input < dimension * dimension; ++input)
    {
        const Eigen::MatrixXcd expected = matrix.adjoint() * paulis[input] * matrix;
        Eigen::MatrixXcd actual = Eigen::MatrixXcd::Zero(dimension, dimension);
        for (unsigned i = table->offsets[input]; i < table->offsets[input + 1]; ++i)
            actual += table->entries[i].coefficient * paulis[table->entries[i].pauli];
        Require((expected - actual).norm() < 3e-12 * std::max(1., expected.norm()),
                "Compiled Pauli transfer differs from independent dense multiplication");
    }
}

inline void ControlledBlocks()
{
    std::vector<Eigen::MatrixXcd> targets;
    Eigen::MatrixXcd arbitrary(2, 2);
    arbitrary << std::complex<double>(.3, .7), std::complex<double>(-.4, .2), std::complex<double>(.6, -.9),
        std::complex<double>(1.1, .3);
    targets.push_back(arbitrary);
    targets.push_back(Gates::HadamardGate<>().getRawOperatorMatrix());
    targets.push_back((1. + 1e-15) * targets.back());
    targets.push_back(Gates::RyGate<>(M_PI_2 + 1e-15).getRawOperatorMatrix());
    targets.push_back(Gates::UGate<>(.731, -.235, .913, .319).getRawOperatorMatrix());
    targets.push_back(Eigen::MatrixXcd::Zero(2, 2));
    arbitrary << 1., 1., 0., 0.;
    targets.push_back(arbitrary);
    for (unsigned qubits = 1; qubits <= 3; ++qubits)
        for (unsigned target = 0; target < qubits; ++target)
            for (unsigned selected = 0; selected < (1U << qubits); ++selected)
            {
                if ((selected >> target) & 1U)
                    continue;
                for (const auto &v : targets)
                {
                    Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Identity(1U << qubits, 1U << qubits);
                    for (unsigned row = 0; row < 2; ++row)
                        for (unsigned col = 0; col < 2; ++col)
                            matrix(selected | (row << target), selected | (col << target)) = v(row, col);
                    CheckTransfer(matrix, qubits, PauliDetail::MatrixTransfer(matrix, qubits));
                }
            }
    const auto nearby = PauliDetail::MatrixTransfer(Gates::RyGate<>(M_PI_2 + 1e-15).getRawOperatorMatrix(), 1);
    Require(Coefficient(*nearby, 2, 2) < -2e-16, "A genuine near-Clifford Pauli coefficient was discarded");
    const auto scaled = PauliDetail::MatrixTransfer(targets[2], 1);
    Require(!(scaled->unchanged & 1U) && Coefficient(*scaled, 0, 0) > 1.,
            "A near-unitary scaled matrix was silently made unitary");
}

inline void Interactions()
{
    for (bool plus : {false, true})
        for (double theta : {0., 1e-12, -.731, M_PI, M_PI_2, M_PI + 1e-15})
            for (double beta : {0., .319, M_PI_2, M_PI_2 + 1e-15})
            {
                const auto table = PauliDetail::CompileXXYY(plus, theta, beta);
                const Eigen::MatrixXcd matrix = plus ? Gates::XXPlusYYGate<>(theta, beta).getRawOperatorMatrix()
                                                     : Gates::XXMinusYYGate<>(theta, beta).getRawOperatorMatrix();
                CheckTransfer(matrix, 2, table);
                Require((table->unchanged & 1U) && table->offsets[1] == 1,
                        "Unitary interaction does not preserve its identity input exactly");
            }
    const double beta = std::nextafter(M_PI_2, std::numeric_limits<double>::infinity());
    const auto nearby = PauliDetail::CompileXXYY(true, M_PI, beta);
    const double expected = std::cos(beta) * std::cos(beta);
    Require(Coefficient(*nearby, 5, 5) > 0. && std::abs(Coefficient(*nearby, 5, 5) / expected - 1.) < 1e-14,
            "XXYY lost the small complement near a quarter-turn beta");
    const auto nearZero = PauliDetail::CompileXXYY(true, M_PI, 1e-16);
    Require(Coefficient(*nearZero, 7, 7) > 0. && std::abs(Coefficient(*nearZero, 7, 7) / 1e-32 - 1.) < 1e-14,
            "XXYY lost the small complement near zero beta");
}

inline void IdentityGrowthAndCache()
{
    PauliPropagator simulator;
    simulator.SetNrQubits(100);
    simulator.SetCoefficientThreshold(0.);
    for (unsigned gate = 0; gate < 8; ++gate)
    {
        simulator.ApplyGate(Gates::CCHGate<>(), 3 * gate, 3 * gate + 1, 3 * gate + 2);
        simulator.ApplyXXPlusYY(3 * gate, 3 * gate + 1, .731, .243);
    }
    PauliStringXZWithCoefficient observable(100);
    observable.Z[99] = true;
    PauliStringStorage output;
    Require(simulator.ExpectationValue(observable, output) == 1. && output.size() == 1,
            "Disjoint unitary gates grew an untouched Pauli observable");
    auto matrix = Gates::CCHGate<>().getRawOperatorMatrix();
    const auto first = PauliDetail::MatrixTransfer(matrix, 3), second = PauliDetail::MatrixTransfer(matrix, 3);
    Require(first == second, "Repeated generic matrices did not share their immutable compilation");
    matrix(6, 6) = std::nextafter(matrix(6, 6).real(), 2.);
    const auto changed = PauliDetail::MatrixTransfer(matrix, 3);
    Require(changed != first, "Matrix mutation reused an obsolete compilation");
    Require(detail::FindCanonicalClifford(Gates::HadamardGate<>().getRawOperatorMatrix()) != nullptr,
            "Canonical H was not recognized");
    matrix = Gates::HadamardGate<>().getRawOperatorMatrix();
    matrix(0, 0) = std::nextafter(matrix(0, 0).real(), 2.);
    Require(detail::FindCanonicalClifford(matrix) == nullptr, "A modified H was snapped to a Clifford");
    matrix = Eigen::MatrixXcd::Identity(2, 2);
    matrix(1, 0) = 1e-200;
    Require(detail::FindCanonicalClifford(matrix) == nullptr,
            "A tiny matrix entry was squared away during recognition");
    const auto tiny = PauliDetail::MatrixTransfer(matrix, 1);
    Require(Coefficient(*tiny, 0, 1) != 0., "A tiny physical Pauli coefficient was snapped to identity");
}

template <class Gate> inline void CheckTypedUnitary(const Gate &gate)
{
    const auto descriptor = Gates::GetUnitaryGateDescriptor(gate);
    const auto &matrix = gate.getRawOperatorMatrix();
    const unsigned arity = unsigned(gate.getQubitsNumber());
    Require(descriptor && Gates::VerifyUnitaryGateDescriptor(matrix, descriptor),
            "A named unitary gate lost its verified constructor metadata");
    const auto table = PauliDetail::MatrixTransfer(matrix, arity, descriptor);
    Require((table->unchanged & 1U) && table->offsets[1] == 1,
            "A verified typed gate branches or rescales an untouched identity");
    CheckTransfer(matrix, arity, table);
}

template <class Gate> inline void CheckTypedFamily(const Gate &gate)
{
    CheckTypedUnitary(gate);
    CheckTypedUnitary(Gates::ControlledGate<>(gate));
    CheckTypedUnitary(Gates::DoublyControlledGate<>(gate));
}

inline void TypedUnitaryPrecision()
{
    CheckTypedFamily(Gates::HadamardGate<>());
    CheckTypedFamily(Gates::HyGate<>());
    CheckTypedFamily(Gates::SGate<>());
    CheckTypedFamily(Gates::SDGGate<>());
    CheckTypedFamily(Gates::TGate<>());
    CheckTypedFamily(Gates::TDGGate<>());
    CheckTypedFamily(Gates::PauliXGate<>());
    CheckTypedFamily(Gates::PauliYGate<>());
    CheckTypedFamily(Gates::PauliZGate<>());
    CheckTypedFamily(Gates::SquareRootNOTGate<>());
    CheckTypedFamily(Gates::SquareRootNOTDagGate<>());
    CheckTypedFamily(Gates::SplitterGate<>());
    for (double angle : {.3, .43, .731, 1e-16, M_PI_2, std::nextafter(M_PI_2, 2.)})
    {
        CheckTypedFamily(Gates::PhaseShiftGate<>(angle));
        CheckTypedFamily(Gates::RxGate<>(angle));
        CheckTypedFamily(Gates::RyGate<>(angle));
        CheckTypedFamily(Gates::RzGate<>(angle));
        CheckTypedFamily(Gates::UGate<>(angle, -.71, .19, .23));
        CheckTypedUnitary(Gates::RXXGate<>(angle));
        CheckTypedUnitary(Gates::RYYGate<>(angle));
        CheckTypedUnitary(Gates::RZZGate<>(angle));
        CheckTypedUnitary(Gates::RZXGate<>(angle));
        CheckTypedUnitary(Gates::RXYGate<>(angle));
        CheckTypedUnitary(Gates::XXPlusYYGate<>(angle, -.29));
        CheckTypedUnitary(Gates::XXMinusYYGate<>(angle, -.29));
    }
    CheckTypedUnitary(Gates::ControlledPhaseShiftGate<>(.3));
    CheckTypedUnitary(Gates::ControlledRxGate<>(.3));
    CheckTypedUnitary(Gates::ControlledRyGate<>(.3));
    CheckTypedUnitary(Gates::ControlledRzGate<>(.3));
    CheckTypedUnitary(Gates::ControlledUGate<>(.43, -.71, .19, .23));
    CheckTypedUnitary(Gates::CCUGate<>(.43, -.71, .19, .23));

    PauliPropagator simulator;
    simulator.SetNrQubits(100);
    simulator.SetCoefficientThreshold(0.);
    for (unsigned gate = 0; gate < 8; ++gate)
    {
        const int q = 3 * gate;
        simulator.ApplyGate(Gates::CCPhaseShiftGate<>(.3), q, q + 1, q + 2);
        simulator.ApplyGate(Gates::CCUGate<>(.43, -.71, .19, .23), q, q + 1, q + 2);
        simulator.ApplyDoublyControlledGate(Gates::RxGate<>(.3), q + 1, q + 2, q);
        simulator.ApplyControlledGate(Gates::UGate<>(.43, -.71, .19, .23), q + 1, q);
        simulator.ApplyGate(Gates::XXPlusYYGate<>(.51, -.29), q, q + 1);
    }
    PauliStringXZWithCoefficient observable(100);
    observable.Z[99] = true;
    PauliStringStorage output;
    Require(simulator.ExpectationValue(observable, output) == 1. && output.size() == 1,
            "Non-Clifford typed C/CC gates grow a disjoint observable");

    Gates::RxGate<> modified(.3);
    const auto descriptor = Gates::GetUnitaryGateDescriptor(modified);
    const auto original = PauliDetail::MatrixTransfer(modified.getRawOperatorMatrix(), 1, descriptor);
    modified.setOperator(Eigen::MatrixXcd(2. * Eigen::MatrixXcd::Identity(2, 2)));
    Require(!Gates::VerifyUnitaryGateDescriptor(modified.getRawOperatorMatrix(), descriptor),
            "Matrix mutation retained verified unitary semantics");
    const auto changed =
        PauliDetail::MatrixTransfer(modified.getRawOperatorMatrix(), 1, Gates::GetUnitaryGateDescriptor(modified));
    Require(changed != original && Coefficient(*changed, 0, 0) == 4., "A mutated typed filter was unitarized");
    Gates::DoublyControlledGate<> controlled(modified);
    const auto filtered =
        PauliDetail::MatrixTransfer(controlled.getRawOperatorMatrix(), 3, Gates::GetUnitaryGateDescriptor(controlled));
    CheckTransfer(controlled.getRawOperatorMatrix(), 3, filtered);
    Require(!(filtered->unchanged & 1U), "Wrapping a mutated typed matrix restored unitary semantics");
    controlled.SetOperation(Gates::HadamardGate<>());
    const auto beforeFailedSet = Gates::GetUnitaryGateDescriptor(controlled);
    bool invalidRejected = false;
    try
    {
        controlled.SetOperation(Eigen::MatrixXcd::Zero(3, 3));
    }
    catch (const std::invalid_argument &)
    {
        invalidRejected = true;
    }
    Require(invalidRejected && Gates::GetUnitaryGateDescriptor(controlled) == beforeFailedSet &&
                Gates::VerifyUnitaryGateDescriptor(controlled.getRawOperatorMatrix(), beforeFailedSet),
            "A failed controlled setter invalidated unchanged unitary metadata");
    auto altered = controlled.getRawOperatorMatrix();
    altered(0, 0) = std::nextafter(1., 2.);
    controlled.setOperator(altered);
    Require(!Gates::VerifyUnitaryGateDescriptor(controlled.getRawOperatorMatrix(),
                                                Gates::GetUnitaryGateDescriptor(controlled)),
            "Controlled provenance ignored a changed inactive sector");

    struct DerivedRotation : Gates::RxGate<>
    {
        void ChangeProtectedMatrix()
        {
            this->operatorMat(1, 0) = 1e-200;
        }
    };

    DerivedRotation derived;
    derived.ChangeProtectedMatrix();
    const auto derivedTable =
        PauliDetail::MatrixTransfer(derived.getRawOperatorMatrix(), 1, Gates::GetUnitaryGateDescriptor(derived));
    Require(!(derivedTable->unchanged & 1U) && Coefficient(*derivedTable, 0, 1) != 0.,
            "A derived mutation bypassed the unitary metadata check");

    const Gates::PhaseShiftGate<> phase(.3);
    const auto typed =
        PauliDetail::MatrixTransfer(phase.getRawOperatorMatrix(), 1, Gates::GetUnitaryGateDescriptor(phase));
    const auto raw = PauliDetail::MatrixTransfer(phase.getRawOperatorMatrix(), 1);
    // The raw entry is compiled from the matrix entries, apart from the typed
    // one. Its identity input may be rescaled by the rounding of the stored
    // matrix, but it has no other output.
    Require(typed != raw && (typed->unchanged & 1U) && raw->KeepsIdentitySupport() &&
                std::abs(Coefficient(*raw, 0, 0) - 1.) < 4e-16,
            "Raw matrix cache entries acquired typed unitary semantics or branch the identity");
    const Gates::RyGate<> nearby(std::nextafter(M_PI_2, 2.));
    const auto nearbyTable =
        PauliDetail::MatrixTransfer(nearby.getRawOperatorMatrix(), 1, Gates::GetUnitaryGateDescriptor(nearby));
    Require(Coefficient(*nearbyTable, 2, 2) != 0., "A neighboring typed rotation was snapped to a Clifford");
    const Gates::CCRxGate<> tinyRotation(1e-16);
    const auto smallTable = PauliDetail::MatrixTransfer(tinyRotation.getRawOperatorMatrix(), 3,
                                                        Gates::GetUnitaryGateDescriptor(tinyRotation));
    Require(Coefficient(*smallTable, 2, 10) > 0., "A tiny controlled diagonal-difference term was lost");
    Require(PauliDetail::InteractionTransfer('X', 'X', 0., -0.) == PauliDetail::InteractionTransfer('X', 'X', -0., 0.),
            "Signed zero split an interaction cache entry");
    Require(PauliDetail::ParameterizedTransfer(OperationType::CU, 0., -0., 0., -0.) ==
                PauliDetail::ParameterizedTransfer(OperationType::CU, -0., 0., -0., 0.),
            "Signed zero split a parameter cache entry");
    for (double angle : {M_PI, 10. * M_PI_2, -M_PI, -10. * M_PI_2})
    {
        const Gates::RxGate<> rx(angle);
        const Gates::RyGate<> ry(angle);
        const Gates::RzGate<> rz(angle);
        const double sign = angle > 0. ? 1. : -1.;
        Require(rx.getRawOperatorMatrix()(0, 0) == std::complex<double>(0.) &&
                    rx.getRawOperatorMatrix()(0, 1) == std::complex<double>(0., -sign),
                "Legacy Rx did not produce an exact encoded half turn");
        Require(ry.getRawOperatorMatrix()(0, 0) == std::complex<double>(0.) &&
                    ry.getRawOperatorMatrix()(1, 0) == std::complex<double>(sign),
                "Legacy Ry did not produce an exact encoded half turn");
        Require(rz.getRawOperatorMatrix()(0, 0) == std::complex<double>(0., -sign) &&
                    rz.getRawOperatorMatrix()(1, 1) == std::complex<double>(0., sign),
                "Legacy Rz retained a quarter-turn rounding residue");
        const Gates::ControlledRxGate<> crx(angle);
        const Gates::CCRxGate<> ccrx(angle);
        Require(Gates::AppliedGate<>(crx.getRawOperatorMatrix()).getMaxBranching() == 1 &&
                    Gates::AppliedGate<>(ccrx.getRawOperatorMatrix()).getMaxBranching() == 1,
                "Controlled quarter-turn matrices retain false branching");
        const double neighbor = std::nextafter(angle, angle > 0. ? 100. : -100.);
        Require(Gates::RxGate<>(neighbor).getRawOperatorMatrix()(0, 0) != std::complex<double>(0.) &&
                    Gates::CCRxGate<>(neighbor).getRawOperatorMatrix()(6, 6) != std::complex<double>(0.),
                "A neighboring controlled rotation was snapped to a quarter turn");
        Require(Gates::PhaseShiftGate<>(angle).getRawOperatorMatrix()(1, 1) == std::complex<double>(-1.) &&
                    Gates::ControlledPhaseShiftGate<>(angle).getRawOperatorMatrix()(3, 3) == std::complex<double>(-1.),
                "Legacy phase gate retained a half-turn rounding residue");
    }
    Require(
        (Gates::UGate<>(M_PI).getRawOperatorMatrix().array() == Gates::RyGate<>(M_PI).getRawOperatorMatrix().array())
            .all(),
        "U did not retain exact quarter-turn trigonometry");
    for (double angle : {1e16, -1e16, 1e19, -1e19})
    {
        long long turns = 0;
        Require(!detail::TryGetQuarterTurns(angle, turns), "A huge ambiguous angle was recognized as a quarter turn");
        Require(Gates::PhaseShiftGate<>(angle).getRawOperatorMatrix()(1, 1) == std::polar(1., angle),
                "A huge legacy phase was snapped instead of using ordinary trigonometry");
        const Gates::RxGate<> rotation(angle);
        Require(rotation.getRawOperatorMatrix()(0, 0) == std::complex<double>(std::cos(.5 * angle)) &&
                    rotation.getRawOperatorMatrix()(0, 1) == std::complex<double>(0., -std::sin(.5 * angle)),
                "A huge legacy rotation was snapped instead of using ordinary trigonometry");
    }
    const Gates::PauliProductRotationGate<'Q', 'X'> invalidAxis(.3);
    Require(!Gates::VerifyUnitaryGateDescriptor(invalidAxis.getRawOperatorMatrix(),
                                                Gates::GetUnitaryGateDescriptor(invalidAxis)),
            "An invalid Pauli axis acquired analytic semantics");
    CheckTransfer(invalidAxis.getRawOperatorMatrix(), 2,
                  PauliDetail::MatrixTransfer(invalidAxis.getRawOperatorMatrix(), 2,
                                              Gates::GetUnitaryGateDescriptor(invalidAxis)));
}

inline void CoherentEulerAndAdjoints()
{
    const double huge = std::numeric_limits<double>::max();
    const std::array<std::array<double, 4>, 9> angles{{{{.43, 2.3, 2.485, 2.665}},
                                                       {{.43, 2.54, 2.707, 2.863}},
                                                       {{.43, 7. * M_PI_2, 8. * M_PI_2, 0.}},
                                                       {{.43, 1., 1e-17, 0.}},
                                                       {{.43, 1e16, -1e16, 1.}},
                                                       {{0., 1e16, -1e16, 1e-300}},
                                                       {{.43, huge, huge, huge}},
                                                       {{0., 0., 0., 1e-300}},
                                                       {{std::nextafter(M_PI_2, 2.), .3, -.7, .2}}}};
    const auto checkAdjoint = [](const auto &gate) {
        const auto &matrix = gate.getRawOperatorMatrix();
        auto descriptor = Gates::GetUnitaryGateDescriptor(gate);
        const unsigned width = unsigned(gate.getQubitsNumber());
        const auto forward = PauliDetail::MatrixTransfer(matrix, width, descriptor);
        descriptor.adjoint = true;
        const Eigen::MatrixXcd inverse = matrix.adjoint();
        Require(Gates::VerifyUnitaryGateDescriptor(inverse, descriptor), "An adjoint lost its verified constructor");
        const auto table = PauliDetail::MatrixTransfer(inverse, width, descriptor);
        Require((table->unchanged & 1U) && table->offsets[1] == 1, "A verified adjoint branches an untouched identity");
        CheckTransfer(inverse, width, table);
        for (unsigned input = 0; input < (1U << (2 * width)); ++input)
            for (unsigned output = 0; output < (1U << (2 * width)); ++output)
                Require(Coefficient(*table, input, output) == Coefficient(*forward, output, input),
                        "An adjoint did not transpose its certified real Pauli transfer");
    };
    for (const auto &a : angles)
    {
        const Gates::UGate<> gate(a[0], a[1], a[2], a[3]);
        const auto &matrix = gate.getRawOperatorMatrix();
        const Eigen::MatrixXcd expected = std::polar(1., a[3]) * Gates::PhaseShiftGate<>(a[1]).getRawOperatorMatrix() *
                                          Gates::RyGate<>(a[0]).getRawOperatorMatrix() *
                                          Gates::PhaseShiftGate<>(a[2]).getRawOperatorMatrix();
        Require(matrix.allFinite() && (matrix - expected).norm() < 3e-15 &&
                    (matrix.adjoint() * matrix - Eigen::MatrixXcd::Identity(2, 2)).norm() < 3e-15,
                "Euler constructor differs from coherent phase-factor composition");
        CheckTypedFamily(gate);
        CheckTypedUnitary(Gates::ControlledUGate<>(a[0], a[1], a[2], a[3]));
        CheckTypedUnitary(Gates::CCUGate<>(a[0], a[1], a[2], a[3]));
        checkAdjoint(gate);
        checkAdjoint(Gates::CCUGate<>(a[0], a[1], a[2], a[3]));
        CheckTransfer(Gates::ControlledUGate<>(a[0], a[1], a[2], a[3]).getRawOperatorMatrix(), 2,
                      PauliDetail::CompileCU(a[0], a[1], a[2], a[3]));
    }
    const Gates::UGate<> tinyGlobal(0., 1e16, -1e16, 1e-300), tinyRelative(.43, 1., 1e-300, 0.);
    for (unsigned diagonal = 0; diagonal < 2; ++diagonal)
    {
        const auto value = tinyGlobal.getRawOperatorMatrix()(diagonal, diagonal);
        Require(value.imag() != 0. && std::abs((value.imag() / value.real()) / 1e-300 - 1.) < 3e-15,
                "Cancelling huge Euler phases erased a tiny global phase");
    }
    Require(tinyRelative.getRawOperatorMatrix()(0, 1).imag() != 0., "A tiny relative Euler phase was erased");
    // Explicit matrices retain their actual entries, including the old phase-sum artifact.
    Gates::UGate<> changed(.3, 1e16, -1e16, 1.);
    Eigen::MatrixXcd legacy(2, 2);
    const double c = std::cos(.15), s = std::sin(.15);
    legacy << std::polar(1., 1.) * c, -std::polar(1., 1. - 1e16) * s, std::polar(1., 1. + 1e16) * s,
        std::polar(1., (1. + 1e16) - 1e16) * c;
    changed.setOperator(legacy);
    Require(!Gates::VerifyUnitaryGateDescriptor(legacy, changed.getUnitaryDescriptor()),
            "An explicitly supplied Euler phase artifact was certified unitary");
    const auto raw = PauliDetail::MatrixTransfer(legacy, 1);
    const auto mutated = PauliDetail::MatrixTransfer(legacy, 1, changed.getUnitaryDescriptor());
    Require(!(mutated->unchanged & 1U), "A mutated Euler filter acquired an exact identity transfer");
    for (unsigned input = 0; input < 4; ++input)
        for (unsigned output = 0; output < 4; ++output)
            Require(Coefficient(*raw, input, output) == Coefficient(*mutated, input, output),
                    "A mutated Euler filter lost a raw transfer coefficient");
    PauliPropagator simulator;
    simulator.SetNrQubits(100);
    simulator.SetCoefficientThreshold(0.);
    for (unsigned q = 0; q < 4; ++q)
    {
        const auto &a = angles[q];
        simulator.ApplyGate(Gates::CCUGate<>(a[0], a[1], a[2], a[3]), 3 * q, 3 * q + 1, 3 * q + 2);
    }
    PauliStringXZWithCoefficient observable(100);
    observable.Z[99] = true;
    PauliStringStorage output;
    Require(simulator.ExpectationValue(observable, output) == 1. && output.size() == 1,
            "Ordinary Euler gates grow a disjoint observable");
    CheckTypedUnitary(Gates::FredkinGate<>());
    checkAdjoint(Gates::FredkinGate<>());
    ExtendedStabilizer approximate(3, ExtendedStabilizerApproximationPolicy::Approximate(0., 1));
    approximate.ApplyRy(0, .317);
    const double bound = approximate.GetApproximationErrorBound();
    approximate.SetApproximationPolicy(ExtendedStabilizerApproximationPolicy::Exact());
    approximate.ApplyGate(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), 0, 1, 2);
    approximate.ApplyGate(Gates::FredkinGate<>(), 0, 1, 2);
    Require(bound > 0. && bound < 1. && approximate.GetApproximationErrorBound() == bound,
            "Certified Euler/Fredkin gates inflated a prior approximation bound");
}

inline void MalformedMatrixRejection()
{
    struct MatrixGate
    {
        Eigen::MatrixXcd matrix;
        unsigned width;

        const Eigen::MatrixXcd &getRawOperatorMatrix() const
        {
            return matrix;
        }

        unsigned getQubitsNumber() const
        {
            return width;
        }
    };

    std::vector<MatrixGate> invalid{{Eigen::MatrixXcd::Zero(0, 0), 0},
                                    {Eigen::MatrixXcd::Zero(2, 1), 1},
                                    {Eigen::MatrixXcd::Zero(2, 3), 1},
                                    {Eigen::MatrixXcd::Zero(3, 3), 2},
                                    {Eigen::MatrixXcd::Zero(4, 8), 2},
                                    {Eigen::MatrixXcd::Zero(8, 4), 3},
                                    {Eigen::MatrixXcd::Identity(16, 16), 4},
                                    {Eigen::MatrixXcd::Identity(2, 2), 32},
                                    {Eigen::MatrixXcd::Identity(2, 2), std::numeric_limits<unsigned>::max()}};
    for (double value : {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()})
    {
        MatrixGate gate{Eigen::MatrixXcd::Identity(8, 8), 3};
        gate.matrix(7, 6) = {0., value};
        invalid.push_back(std::move(gate));
    }
    PauliPropagator pauli;
    pauli.SetNrQubits(3);
    pauli.ApplyGate(Gates::RyGate<>(.317), 0);
    const auto operations = pauli.GetOperations().size();
    ExtendedStabilizer extended(3);
    extended.ApplyRy(0, .317);
    const auto amplitudes = extended.GetFrames().front().amplitudes;
    const auto signs = extended.GetFrames().front().signs;
    const auto valid = Gates::CCRxGate<>(.37);
    const auto cached = PauliDetail::MatrixTransfer(valid.getRawOperatorMatrix(), 3, valid.getUnitaryDescriptor());
    for (unsigned controls : {3U, 32U, std::numeric_limits<unsigned>::max()})
    {
        auto descriptor = valid.getUnitaryDescriptor();
        descriptor.controls = controls;
        Require(!Gates::VerifyUnitaryGateDescriptor(valid.getRawOperatorMatrix(), descriptor),
                "An invalid control count reached descriptor dimensions or shifts");
    }
    for (const auto &gate : invalid)
    {
        const auto rejects = [](const auto &call) {
            try
            {
                call();
            }
            catch (const std::invalid_argument &)
            {
                return;
            }
            throw std::runtime_error("A malformed local matrix was accepted");
        };
        if (gate.matrix.rows() != 2 || gate.matrix.cols() != 2 || !gate.matrix.allFinite())
            rejects([&] { detail::MatrixKey(gate.matrix); });
        rejects([&] { detail::MatrixPauliSum(gate.matrix, gate.width); });
        rejects([&] { PauliDetail::MatrixTransfer(gate.matrix, gate.width); });
        rejects([&] { pauli.ApplyGate(gate, 0, 1, 2); });
        rejects([&] { extended.ApplyGate(gate, 0, 1, 2); });
        Require(pauli.GetOperations().size() == operations && extended.GetFrames().front().amplitudes == amplitudes &&
                    extended.GetFrames().front().signs == signs,
                "Malformed-matrix rejection changed a simulator");
        Require(PauliDetail::MatrixTransfer(valid.getRawOperatorMatrix(), 3, valid.getUnitaryDescriptor()) == cached,
                "Malformed-matrix rejection polluted a valid cached transfer");
    }
}

inline void ExactInverseDescriptors()
{
    using Descriptor = Gates::UnitaryGateDescriptor;
    using Kind = Descriptor::Kind;
    const auto inverse = Gates::AreUnitaryGateDescriptorsInverse;
    const auto checkPair = [&](Descriptor a, Descriptor b) {
        for (bool adjoint : {false, true})
        {
            a.adjoint = b.adjoint = adjoint;
            Require(inverse(a, b) && inverse(b, a), "An exact constructor inverse was not recognized");
            const auto first = Gates::UnitaryGateMatrix(a), second = Gates::UnitaryGateMatrix(b);
            Require(first.size() && second.size() &&
                        (second * first - Eigen::MatrixXcd::Identity(first.rows(), first.cols())).norm() < 1e-14,
                    "An inverse descriptor relation disagrees with coherent constructor matrices");
        }
        b = a;
        b.adjoint = !a.adjoint;
        Require(inverse(a, b) && inverse(b, a), "An exact descriptor-adjoint pair was not recognized");
    };
    for (unsigned controls = 0; controls < 3; ++controls)
    {
        for (Kind kind : {Kind::Identity, Kind::H, Kind::Hy, Kind::X, Kind::Y, Kind::Z})
            checkPair({kind, {}, controls}, {kind, {}, controls});
        for (auto kinds : {std::pair<Kind, Kind>{Kind::S, Kind::Sdg}, {Kind::T, Kind::Tdg}, {Kind::Sx, Kind::Sxdg}})
            checkPair({kinds.first, {}, controls}, {kinds.second, {}, controls});
        for (Kind kind : {Kind::Phase, Kind::Rx, Kind::Ry, Kind::Rz})
            for (double angle : {0., .3, M_PI_2, 1e-300, std::numeric_limits<double>::max()})
            {
                const Descriptor a{kind, {angle}, controls}, b{kind, {-angle}, controls};
                checkPair(a, b);
                auto nearby = b;
                nearby.angles[0] = std::nextafter(nearby.angles[0], std::numeric_limits<double>::infinity());
                Require(!inverse(a, nearby), "A neighboring or tiny rotation was mistaken for an exact inverse");
                nearby = a;
                nearby.controls = (controls + 1) % 3;
                nearby.adjoint = true;
                Require(!inverse(a, nearby), "Inverse metadata ignored the number of controls");
            }
        for (const auto &angles :
             {std::array<double, 4>{.43, 2.54, 2.707, 2.863}, std::array<double, 4>{.3, 1e16, -1e16, 1e-300}})
        {
            const Descriptor a{Kind::U, angles, controls};
            const Descriptor b{Kind::U, {-angles[0], -angles[2], -angles[1], -angles[3]}, controls};
            checkPair(a, b);
            for (unsigned parameter = 0; parameter < 4; ++parameter)
            {
                auto nearby = b;
                nearby.angles[parameter] =
                    std::nextafter(nearby.angles[parameter], std::numeric_limits<double>::infinity());
                Require(!inverse(a, nearby), "A neighboring Euler parameter was mistaken for an exact inverse");
            }
        }
    }
    checkPair({Kind::Fredkin}, {Kind::Fredkin});
    for (char first : {'X', 'Y', 'Z'})
        for (char second : {'X', 'Y', 'Z'})
            checkPair({Kind::PairRotation, {.37}, 0, first, second}, {Kind::PairRotation, {-.37}, 0, first, second});
    for (Kind kind : {Kind::XXPlusYY, Kind::XXMinusYY})
    {
        const Descriptor a{kind, {.37, 1e-300}}, b{kind, {-.37, 1e-300}};
        checkPair(a, b);
        auto nearby = b;
        nearby.angles[1] = std::nextafter(nearby.angles[1], 1.);
        Require(!inverse(a, nearby), "An inverse XX/YY descriptor ignored a tiny beta change");
    }
    Descriptor splitter{Kind::Splitter};
    Require(!inverse(splitter, splitter), "Splitter was incorrectly declared self-inverse");
    auto splitterAdjoint = splitter;
    splitterAdjoint.adjoint = true;
    Require(inverse(splitter, splitterAdjoint), "Splitter's exact adjoint was not recognized");
    Descriptor unused{Kind::H, {1e16, 1e-300, -.7, .23}, 0, 'Q', 'Q', true};
    Require(inverse(Descriptor{Kind::H}, unused), "Unused fixed-gate fields changed its inverse relation");
    Require(!inverse(Descriptor{}, Descriptor{}) &&
                !inverse(Descriptor{static_cast<Kind>(255)}, Descriptor{static_cast<Kind>(255)}),
            "An unsupported descriptor acquired inverse semantics");
    Require(!inverse(Descriptor{Kind::PairRotation, {.37}, 0, 'X', 'Y'},
                     Descriptor{Kind::PairRotation, {-.37}, 0, 'Y', 'X'}),
            "Pair-rotation inverse metadata exchanged its physical axes");
}

inline void CrossFamilyInverseDescriptors()
{
    using Descriptor = Gates::UnitaryGateDescriptor;
    using Kind = Descriptor::Kind;
    const auto inverse = Gates::AreUnitaryGateDescriptorsInverse;
    const auto check = [&](const Descriptor &a, const Descriptor &b, bool expected) {
        Require(inverse(a, b) == expected && inverse(b, a) == expected,
                "An exact cross-family inverse relation changed with operand order");
        if (expected)
        {
            const auto first = Gates::UnitaryGateMatrix(a), second = Gates::UnitaryGateMatrix(b);
            Require(first.size() && second.size() &&
                        (second * first - Eigen::MatrixXcd::Identity(first.rows(), first.cols())).norm() < 1e-14,
                    "A cross-family inverse disagrees with its verified gate matrices");
        }
    };
    for (unsigned controls = 0; controls < 3; ++controls)
    {
        std::vector<std::pair<Descriptor, Descriptor>> aliases;
        const std::array<double, 7> phases{{0., M_PI / 4., -M_PI / 4., M_PI_2, -M_PI_2, M_PI, -M_PI}};
        const std::array<Kind, 7> kinds{{Kind::Identity, Kind::T, Kind::Tdg, Kind::S, Kind::Sdg, Kind::Z, Kind::Z}};
        for (size_t i = 0; i < phases.size(); ++i)
            aliases.push_back({{Kind::Phase, {phases[i]}, controls}, {kinds[i], {}, controls}});
        for (double angle : {0., .731, 1e-300, std::numeric_limits<double>::max()})
            aliases.push_back({{Kind::U, {angle, 0., 0., 0.}, controls}, {Kind::Ry, {angle}, controls}});
        for (double angle : {0., .3, M_PI / 4., 1e-300, std::numeric_limits<double>::max()})
        {
            aliases.push_back({{Kind::U, {0., angle, 0., 0.}, controls}, {Kind::Phase, {angle}, controls}});
            aliases.push_back({{Kind::U, {0., 0., angle, 0.}, controls}, {Kind::Phase, {angle}, controls}});
        }
        aliases.push_back({{Kind::U, {0., M_PI / 4., 0., 0.}, controls}, {Kind::T, {}, controls}});
        aliases.push_back({{Kind::U, {0., 0., -M_PI_2, 0.}, controls}, {Kind::Sdg, {}, controls}});
        for (const auto &alias : aliases)
        {
            for (bool adjoint : {false, true})
            {
                auto a = alias.first, b = alias.second;
                a.adjoint = adjoint;
                b.adjoint = !adjoint;
                check(a, b, true);
                b.controls = (controls + 1) % 3;
                check(a, b, false);
            }
            auto a = alias.first, b = alias.second;
            switch (b.kind)
            {
            case Kind::T:
                b.kind = Kind::Tdg;
                break;
            case Kind::Tdg:
                b.kind = Kind::T;
                break;
            case Kind::S:
                b.kind = Kind::Sdg;
                break;
            case Kind::Sdg:
                b.kind = Kind::S;
                break;
            case Kind::Phase:
            case Kind::Ry:
                b.angles[0] = -b.angles[0];
                break;
            default:
                break;
            }
            check(a, b, true);
            a.adjoint = b.adjoint = true;
            check(a, b, true);
            a = alias.first;
            b = alias.second;
            b.adjoint = true;
            if (a.kind == Kind::Phase)
            {
                a.angles[0] = std::nextafter(a.angles[0], std::numeric_limits<double>::infinity());
                check(a, b, false);
            }
            else
            {
                // Every exact-zero Euler reduction prerequisite is significant,
                // even when changing it by the smallest representable amount.
                for (unsigned parameter = 0; parameter < 4; ++parameter)
                    if (alias.first.angles[parameter] == 0.)
                    {
                        a = alias.first;
                        a.angles[parameter] = std::numeric_limits<double>::denorm_min();
                        check(a, b, false);
                    }
                if (b.kind == Kind::Phase || b.kind == Kind::Ry)
                {
                    a = alias.first;
                    b.angles[0] = std::nextafter(b.angles[0], std::numeric_limits<double>::infinity());
                    check(a, b, false);
                }
            }
        }
        check({Kind::U, {0., 1e16, -1e16, 0.}, controls}, {Kind::Identity, {}, controls}, false);
        check({Kind::U, {.731, 0., 0., 2. * M_PI}, controls}, {Kind::Ry, {-.731}, controls}, false);
        check({Kind::Phase, {.3}, controls}, {Kind::Rx, {-.3}, controls}, false);
    }
}

inline void ExtendedFilters()
{
    ExtendedStabilizer canonical(12);
    for (unsigned q = 0; q < 12; ++q)
        canonical.ApplyGate(Gates::HadamardGate<>(), q);
    Require(canonical.GetFrames().front().GetFrameSize() == 1,
            "Canonical generic Cliffords expanded the stabilizer frame");
    for (unsigned q = 0; q < 6; ++q)
        canonical.ApplyControlledGate(Gates::PauliXGate<>(), 2 * q, 2 * q + 1);
    Require(canonical.GetFrames().front().GetFrameSize() == 1, "Generic controlled X expanded the stabilizer frame");

    ExtendedStabilizer nearby(1);
    nearby.ApplyRy(0, .731);
    nearby.ApplyGate(Gates::RyGate<>(-.731 + 1e-15), 0);
    Require(nearby.GetFrames().front().GetFrameSize() == 2 && std::abs(nearby.ExpectationValue("X")) > 2e-16,
            "Generic frame multiplication discarded a genuine cancellation coefficient");
    for (double scale : {2., 1e-300, 1e300, 1e-310, std::numeric_limits<double>::denorm_min()})
    {
        ExtendedStabilizer filtered(1);
        const Gates::AppliedGate<> scalar(Eigen::MatrixXcd(scale * Eigen::MatrixXcd::Identity(2, 2)));
        filtered.ApplyGate(scalar, 0);
        filtered.ApplyRy(0, .731);
        Require(std::abs(filtered.ExpectationValue("Z") - std::cos(.731)) < 2e-14,
                "Scalar filter broke normalized stabilizer expectation values");
        const auto amplitudes = filtered.GetFrames().front().amplitudes;
        bool rejected = false;
        try
        {
            filtered.ApplyGate(Gates::AppliedGate<>(Eigen::MatrixXcd::Zero(2, 2)), 0);
        }
        catch (const std::domain_error &)
        {
            rejected = true;
        }
        Require(rejected && filtered.GetFrames().front().amplitudes == amplitudes,
                "An annihilating scalar filter changed the state before rejection");
    }
    ExtendedStabilizer zero(1);
    Eigen::MatrixXcd projector = Eigen::MatrixXcd::Zero(2, 2);
    projector(1, 1) = 1.;
    bool rejected = false;
    try
    {
        zero.ApplyGate(Gates::AppliedGate<>(projector), 0);
    }
    catch (const std::domain_error &)
    {
        rejected = true;
    }
    Require(rejected && zero.ExpectationValue("Z") == 1.,
            "A state-dependent annihilating filter changed the original frame");
    ExtendedStabilizer nonunitary(1);
    nonunitary.ApplyH(0);
    Eigen::MatrixXcd filter = Eigen::MatrixXcd::Zero(2, 2);
    filter(0, 0) = 2.;
    filter(1, 1) = 1.;
    nonunitary.ApplyGate(Gates::AppliedGate<>(filter), 0);
    Require(std::abs(nonunitary.ExpectationValue("Z") - .6) < 2e-14,
            "Generic nonunitary filter did not normalize its output");
    ExtendedStabilizer illConditioned(1);
    illConditioned.setToBasisState(1);
    filter(0, 0) = 1.;
    filter(1, 1) = 1e-15;
    illConditioned.ApplyGate(Gates::AppliedGate<>(filter), 0);
    Require(illConditioned.ExpectationValue("Z") == -1., "A near-singular nonzero filter was rejected");
    ExtendedStabilizer tinyOutput(1);
    filter << 0., 0., 1e-300, 1.;
    tinyOutput.ApplyGate(Gates::AppliedGate<>(filter), 0);
    Require(tinyOutput.ExpectationValue("Z") == -1., "A tiny filtered state was mistaken for annihilation");
    const auto policy = ExtendedStabilizerApproximationPolicy::Approximate(0., 1);
    ExtendedStabilizer approximate(1, policy);
    approximate.ApplyRy(0, .731);
    const auto statistics = approximate.GetApproximationStatistics();
    rejected = false;
    try
    {
        approximate.ApplyGate(Gates::AppliedGate<>(projector), 0);
    }
    catch (const std::domain_error &)
    {
        rejected = true;
    }
    Require(rejected && approximate.ExpectationValue("Z") == 1. &&
                approximate.GetApproximationStatistics().traceDistanceErrorBound ==
                    statistics.traceDistanceErrorBound &&
                approximate.GetApproximationStatistics().pruningEvents == statistics.pruningEvents,
            "Rejected approximate filtering changed the state or pruning statistics");
    filter << 1., 0., 0., 2.;
    approximate.ApplyGate(Gates::AppliedGate<>(filter), 0);
    Require(approximate.GetApproximationPolicy().mode == policy.mode &&
                approximate.GetApproximationPolicy().maxComponents == policy.maxComponents &&
                approximate.GetApproximationStatistics().traceDistanceErrorBound == 1.,
            "Filtering changed approximation policy or understated a conditioned error bound");
}

inline void MatrixKeyPrefixCopies()
{
    detail::LocalMatrixKey empty;
    const detail::LocalMatrixKey emptyCopy(empty);
    Require(emptyCopy == empty && emptyCopy.dimension == 0 && emptyCopy.hash == empty.hash,
            "Copying a default matrix key changed its empty prefix");
    Eigen::MatrixXcd one = Eigen::MatrixXcd::Identity(2, 2), three(8, 8);
    one(1, 0) = {1e-300, -.37};
    for (Eigen::Index row = 0; row < 8; ++row)
        for (Eigen::Index col = 0; col < 8; ++col)
            three(row, col) = {double(1 + row + 8 * col), double(row - col)};
    const auto oneKey = detail::MatrixKey(one), threeKey = detail::MatrixKey(three);
    detail::LocalMatrixKey assigned(oneKey);
    const auto check = [&](const detail::LocalMatrixKey &expected) {
        Require(assigned == expected && expected == assigned && assigned.hash == expected.hash,
                "Matrix-key copy or assignment lost its active prefix or hash");
    };
    check(detail::MatrixKey(one));
    assigned = threeKey;
    check(detail::MatrixKey(three));
    assigned = oneKey;
    check(detail::MatrixKey(one));
    assigned = assigned;
    check(oneKey);
    auto fresh = detail::MatrixKey(one);
    std::fill(fresh.bits.begin() + 8, fresh.bits.end(), UINT64_C(0xAAAAAAAAAAAAAAAA));
    check(fresh); // the old 8x8 tail cannot participate in 2x2 equality
    one(0, 1) = {-0., -0.};
    check(detail::MatrixKey(one));
    assigned = empty;
    check(emptyCopy);
    assigned = threeKey;
    const detail::LocalMatrixKey copied(assigned);
    Require(copied == threeKey && copied.hash == threeKey.hash,
            "Copying a full-size matrix key omitted a matrix coefficient");

    for (Eigen::Index dimension : {2, 4, 8})
    {
        Eigen::MatrixXcd matrix(dimension, dimension);
        for (Eigen::Index row = 0; row < dimension; ++row)
            for (Eigen::Index col = 0; col < dimension; ++col)
                matrix(row, col) = {double(1 + row + dimension * col), double(2 * row - col)};
        matrix(0, 0) = {-0., -0.};
        matrix(0, 1) = {std::numeric_limits<double>::denorm_min(), -1e-300};
        matrix(1, 0) = {-std::numeric_limits<double>::max(), std::numeric_limits<double>::denorm_min()};
        matrix(1, 1) = {std::numeric_limits<double>::min(), -std::numeric_limits<double>::max()};
        const auto checkAdjoint = [&](const auto &input) {
            const auto key = detail::MatrixKey(input);
            Require(detail::AdjointMatrixKeyHash(key) == detail::MatrixKey(input.adjoint()).hash,
                    "Bitwise adjoint hashing changed a matrix coefficient or its logical position");
            auto changedTail = key;
            std::fill(changedTail.bits.begin() + 2 * dimension * dimension, changedTail.bits.end(),
                      UINT64_C(0xAAAAAAAAAAAAAAAA));
            Require(detail::AdjointMatrixKeyHash(changedTail) == detail::AdjointMatrixKeyHash(key),
                    "Adjoint hashing read beyond the matrix key's active prefix");
        };
        checkAdjoint(matrix);
        const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMajor = matrix;
        checkAdjoint(rowMajor);
        const Eigen::Index outerStride = 2 * dimension + 3;
        std::vector<std::complex<double>> storage(static_cast<size_t>(outerStride * dimension));
        Eigen::Map<Eigen::MatrixXcd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> strided(
            storage.data(), dimension, dimension, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(outerStride, 2));
        strided = matrix;
        checkAdjoint(strided);
    }
}

// The operator of a matrix on the listed qubits of a register (local bit i is qubits[i]).
inline Eigen::MatrixXcd Embed(const Eigen::MatrixXcd &matrix, const std::vector<size_t> &qubits, unsigned total)
{
    const unsigned dimension = 1U << total;
    Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(dimension, dimension);
    for (unsigned col = 0; col < dimension; ++col)
    {
        unsigned input = 0, rest = col;
        for (size_t bit = 0; bit < qubits.size(); ++bit)
        {
            input |= ((col >> qubits[bit]) & 1U) << bit;
            rest &= ~(1U << qubits[bit]);
        }
        for (unsigned row = 0; row < unsigned(matrix.rows()); ++row)
        {
            unsigned output = rest;
            for (size_t bit = 0; bit < qubits.size(); ++bit)
                output |= ((row >> bit) & 1U) << qubits[bit];
            result(output, col) += matrix(row, input);
        }
    }
    return result;
}

// Every Pauli expectation of a three-qubit simulator against a dense state.
inline void CheckExpectations(const PauliPropagator &simulator, const Eigen::VectorXcd &state, const char *message)
{
    const char letters[] = {'I', 'X', 'Z', 'Y'};
    for (unsigned code = 0; code < 64; ++code)
    {
        std::string pauli(3, 'I');
        for (unsigned q = 0; q < 3; ++q)
            pauli[q] = letters[(code >> (2 * q)) & 3U];
        const double expected = (state.adjoint() * Pauli(code, 3) * state)(0, 0).real();
        Require(std::abs(simulator.ExpectationValue(pauli) - expected) < 2e-12, message);
    }
}

inline bool SameTable(const PauliDetail::LocalTransfer &a, const PauliDetail::LocalTransfer &b)
{
    if (a.qubits != b.qubits || a.offsets != b.offsets || a.entries.size() != b.entries.size() ||
        a.unchanged != b.unchanged || a.clifford != b.clifford || a.maxOutputs != b.maxOutputs)
        return false;
    for (size_t i = 0; i < a.entries.size(); ++i)
        if (a.entries[i].pauli != b.entries[i].pauli || a.entries[i].coefficient != b.entries[i].coefficient)
            return false;
    return true;
}

// A matrix without a constructor candidate is compiled from its entries. The
// rounding of a cancellation must not become a table entry, and nothing else
// may be dropped or rescaled.
inline void RawMatrixResidue()
{
    using Applied = Gates::AppliedGate<>;
    const Gates::UGate<> u(.43, -.71, .19, .23);
    const Gates::ControlledUGate<> cu(.43, -.71, .19, .23);
    const Gates::CCUGate<> ccu(.43, -.71, .19, .23);
    const auto check = [](const auto &gate, unsigned qubits) {
        const Eigen::MatrixXcd matrix = gate.getRawOperatorMatrix();
        const auto raw = PauliDetail::MatrixTransfer(matrix, qubits);
        const auto typed = PauliDetail::MatrixTransfer(matrix, qubits, Gates::GetUnitaryGateDescriptor(gate));
        CheckTransfer(matrix, qubits, raw);
        Require(raw != typed && raw->KeepsIdentitySupport() && std::abs(Coefficient(*raw, 0, 0) - 1.) < 4e-16 &&
                    !raw->nonunitary,
                "A raw unitary matrix branches or visibly rescales the identity");
        // The same outputs as the analytic table, and nothing of rounding size.
        Require(raw->offsets == typed->offsets && raw->maxOutputs == typed->maxOutputs,
                "A raw unitary matrix has rounding residue among its table entries");
        for (size_t i = 0; i < raw->entries.size(); ++i)
            Require(raw->entries[i].pauli == typed->entries[i].pauli &&
                        std::abs(raw->entries[i].coefficient - typed->entries[i].coefficient) < 4e-16,
                    "A raw unitary matrix differs from its analytic transfer");
    };
    check(u, 1);
    check(cu, 2);
    check(ccu, 3);
    check(Gates::PhaseShiftGate<>(.3), 1);
    check(Gates::ControlledRyGate<>(-.717), 2);
    check(Gates::CCPhaseShiftGate<>(.3), 3);

    // A typed gate whose matrix was replaced keeps its candidate, fails the
    // verification and takes the same path as the raw matrix.
    Gates::RxGate<> replaced(.3);
    replaced.setOperator(Eigen::MatrixXcd(u.getRawOperatorMatrix()));
    const auto unverified =
        PauliDetail::MatrixTransfer(replaced.getRawOperatorMatrix(), 1, Gates::GetUnitaryGateDescriptor(replaced));
    Require(SameTable(*unverified, *PauliDetail::MatrixTransfer(u.getRawOperatorMatrix(), 1)),
            "An unverified typed gate did not take the raw matrix path");

    // Gates given as matrices on qubits the observable does not touch leave it alone.
    PauliPropagator simulator;
    simulator.SetNrQubits(100);
    simulator.SetCoefficientThreshold(0.);
    for (unsigned gate = 0; gate < 8; ++gate)
    {
        const size_t q = 3 * gate;
        simulator.ApplyGate(Applied(Gates::UGate<>(.3 + gate, .2, .1, .4).getRawOperatorMatrix(), q));
        simulator.ApplyGate(
            Applied(Gates::ControlledUGate<>(.5 + gate, .3, -.2, .11).getRawOperatorMatrix(), q + 1, q));
        simulator.ApplyGate(
            Applied(Gates::CCUGate<>(.7 + gate, -.33, .17, -.41).getRawOperatorMatrix(), q, q + 1, q + 2));
    }
    PauliStringXZWithCoefficient observable(100);
    observable.Z[99] = true;
    PauliStringStorage output;
    Require(std::abs(simulator.ExpectationValue(observable, output) - 1.) < 1e-14 && output.size() == 1,
            "Raw unitary matrices grew an observable they do not touch");

    // The same circuit from matrices and from typed gates reaches the same terms.
    const auto brickwork = [](bool raw, PauliStringStorage &terms) {
        PauliPropagator sim;
        sim.SetNrQubits(6);
        sim.SetCoefficientThreshold(0.);
        sim.SetStepsBetweenDeduplication(1);
        for (int layer = 0; layer < 2; ++layer)
        {
            for (int q = 0; q < 6; ++q)
            {
                const Gates::UGate<> gate(.3 + .1 * q + layer, .2 * q, -.15 * layer + .4);
                if (raw)
                    sim.ApplyGate(Applied(gate.getRawOperatorMatrix(), size_t(q)));
                else
                    sim.ApplyGate(gate, q);
            }
            for (int q = layer % 2; q + 1 < 6; q += 2)
            {
                const Gates::ControlledUGate<> gate(.5 + .07 * q, .3, -.2 * layer, .11);
                if (raw)
                    sim.ApplyGate(Applied(gate.getRawOperatorMatrix(), size_t(q), size_t(q + 1)));
                else
                    sim.ApplyGate(gate, q, q + 1);
            }
        }
        return sim.ExpectationValue(std::string("Z"), terms);
    };
    PauliStringStorage rawTerms, typedTerms;
    const double rawValue = brickwork(true, rawTerms), typedValue = brickwork(false, typedTerms);
    Require(rawTerms.size() == typedTerms.size() && std::abs(rawValue - typedValue) < 1e-13,
            "A circuit of raw matrices has more terms than the same typed circuit");

    // What is not rounding survives: a tiny product, a small cancellation
    // that exceeds its rounding, a scale and a filter.
    const Eigen::MatrixXcd identity2 = Eigen::MatrixXcd::Identity(2, 2), identity4 = Eigen::MatrixXcd::Identity(4, 4);
    const auto embedded = [&](const Eigen::MatrixXcd &single) { return Embed(single, {0}, 3); };
    const auto tiny = PauliDetail::MatrixTransfer(embedded(Gates::RxGate<>(1e-9).getRawOperatorMatrix()), 3);
    Require(std::abs(std::abs(Coefficient(*tiny, 3, 2)) / std::sin(1e-9) - 1.) < 1e-15 &&
                Coefficient(*tiny, 2, 3) == -Coefficient(*tiny, 3, 2),
            "A tiny rotation coefficient of a raw matrix was dropped");
    const auto nearby =
        PauliDetail::MatrixTransfer(embedded(Gates::RyGate<>(M_PI_2 + 1e-15).getRawOperatorMatrix()), 3);
    Require(Coefficient(*nearby, 2, 2) < -2e-16,
            "A genuine near-Clifford coefficient of an embedded matrix was dropped");
    const auto scaled = PauliDetail::MatrixTransfer(Eigen::MatrixXcd((1. + 1e-15) * ccu.getRawOperatorMatrix()), 3);
    Require(scaled->KeepsIdentitySupport() && Coefficient(*scaled, 0, 0) > 1. + 1e-15 && !(scaled->unchanged & 1U),
            "A scaled raw matrix lost its scale");
    Eigen::MatrixXcd filter = identity2;
    filter(1, 1) = .5;
    const auto filtered = PauliDetail::MatrixTransfer(filter, 1);
    Require(filtered->nonunitary && Coefficient(*filtered, 0, 0) == .625 && Coefficient(*filtered, 0, 2) == .375,
            "A raw filter lost its actual action");
    Eigen::MatrixXcd shear = identity4;
    shear(3, 0) = 1e-200;
    const auto sheared = PauliDetail::MatrixTransfer(shear, 2);
    unsigned nonzero = 0;
    for (unsigned output = 1; output < 16; ++output)
        nonzero += Coefficient(*sheared, 0, output) != 0.;
    // 1e-200 (|0><3| + |3><0|) is 1e-200 (XX - YY) / 2.
    Require(nonzero == 2 && Coefficient(*sheared, 0, 5) != 0. && Coefficient(*sheared, 0, 15) != 0.,
            "A tiny entry of a raw two-qubit matrix was dropped");

    // A gate that only rescales the identity must be applied to every term but
    // does not bring the observable onto its qubits: gates before it on those
    // qubits stay outside the cone.
    PauliPropagator cone;
    cone.SetNrQubits(8);
    cone.ApplyRY(5, .4);
    cone.ApplyGate(Applied(Eigen::MatrixXcd((1. + 1e-15) * Gates::HadamardGate<>().getRawOperatorMatrix()), 5));
    cone.ApplyRY(0, .3);
    std::vector<PauliDetail::Operation> operations;
    for (auto &op : cone.GetOperations())
        operations.push_back(PauliDetail::Operation::Import(std::move(op)));
    Require(operations.size() == 3 && IsPauliLocalGate(operations[1].type) &&
                operations[1].Transfer().KeepsIdentitySupport() && !(operations[1].Transfer().unchanged & 1U),
            "The scaled matrix was not recorded as a rescaling table");
    std::vector<PauliDetail::Term<1>> terms(1, PauliDetail::Term<1>(8));
    terms[0].Z[0] = 1;
    const auto needed = PauliDetail::BackwardCone(terms, operations, 8);
    Require(needed == std::vector<unsigned char>({0, 1, 1}), "A rescaling gate widened the backward cone");
    output.clear();
    const double rescaled = cone.ExpectationValue(std::string("Z"), output);
    Require(output.size() == 2 &&
                std::abs(rescaled - std::cos(.3) * Coefficient(operations[1].Transfer(), 0, 0)) < 1e-15,
            "A rescaling gate outside the cone was not applied");
}

// A qubit that no entry of a matrix changes (a control in any position and
// for any control value, a spectator, every qubit of a diagonal matrix) has
// the columns of the inputs with Z or Y on it copied from the columns without
// them. Every column must still be the action of the matrix, entry by entry,
// with nothing of rounding size kept and nothing else dropped.
inline void FixedQubitColumns()
{
    const auto check = [](const Eigen::MatrixXcd &matrix, unsigned qubits) {
        const auto table = PauliDetail::CompileMatrix(matrix, qubits);
        const unsigned dimension = 1U << qubits, count = dimension * dimension;
        std::vector<Eigen::MatrixXcd> paulis;
        for (unsigned code = 0; code < count; ++code)
            paulis.push_back(Pauli(code, qubits));
        for (unsigned input = 0; input < count; ++input)
        {
            const Eigen::MatrixXcd transformed = matrix.adjoint() * paulis[input] * matrix;
            for (unsigned output = 0; output < count; ++output)
            {
                const double expected = (paulis[output] * transformed).trace().real() / dimension;
                const double actual = Coefficient(*table, input, output);
                Require(std::abs(actual - expected) < 4e-15 && (actual != 0.) == (std::abs(expected) > 1e-13),
                        "A column of a matrix with unchanged qubits differs from dense multiplication");
            }
        }
    };
    // The single-qubit block on the target, selected by the given values of the other qubits.
    const auto controlled = [](const Eigen::MatrixXcd &block, unsigned qubits, unsigned target, unsigned value) {
        const unsigned dimension = 1U << qubits, others = (dimension - 1U) & ~(1U << target);
        Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Identity(dimension, dimension);
        for (unsigned row = 0; row < 2; ++row)
            for (unsigned col = 0; col < 2; ++col)
                matrix((value & others) | (row << target), (value & others) | (col << target)) = block(row, col);
        return matrix;
    };
    Eigen::MatrixXcd arbitrary(2, 2);
    arbitrary << std::complex<double>(.3, .7), std::complex<double>(-.4, .2), std::complex<double>(.6, -.9),
        std::complex<double>(1.1, .3);
    const std::vector<Eigen::MatrixXcd> blocks = {Gates::UGate<>(.43, -.71, .19, .23).getRawOperatorMatrix(),
                                                  Gates::RyGate<>(-.717).getRawOperatorMatrix(),
                                                  Gates::PhaseShiftGate<>(.3).getRawOperatorMatrix(), arbitrary};
    for (const auto &block : blocks)
        for (unsigned qubits = 2; qubits <= 3; ++qubits)
            for (unsigned target = 0; target < qubits; ++target)
                for (unsigned value = 0; value < (1U << qubits); ++value)
                    if (!((value >> target) & 1U))
                        check(controlled(block, qubits, target, value), qubits);
    // A gate on some qubits of a larger identity: the others are spectators.
    for (size_t q = 0; q < 3; ++q)
    {
        check(Embed(blocks[0], {q}, 3), 3);
        check(Embed(Gates::ControlledUGate<>(.5, .3, -.2, .11).getRawOperatorMatrix(), {q, (q + 1) % 3}, 3), 3);
    }
    // Diagonal matrices, unitary or not, change no qubit at all.
    for (unsigned qubits = 1; qubits <= 3; ++qubits)
    {
        const unsigned dimension = 1U << qubits;
        Eigen::MatrixXcd phases = Eigen::MatrixXcd::Zero(dimension, dimension), filter = phases;
        for (unsigned k = 0; k < dimension; ++k)
        {
            phases(k, k) = std::polar(1., .37 + 1.3 * k * k);
            filter(k, k) = std::complex<double>(.4 + .17 * k, .05 * k * k - .3);
        }
        check(phases, qubits);
        check(filter, qubits);
    }
    // A different block for every value of the controls.
    for (unsigned target = 0; target < 3; ++target)
    {
        Eigen::MatrixXcd multiplexed = Eigen::MatrixXcd::Zero(8, 8);
        unsigned index = 0;
        for (unsigned value = 0; value < 8; ++value)
            if (!((value >> target) & 1U))
            {
                const Eigen::MatrixXcd block =
                    Gates::UGate<>(.3 + .41 * index, .2 - .13 * index, -.5 + .29 * index, .1 * index)
                        .getRawOperatorMatrix();
                ++index;
                for (unsigned row = 0; row < 2; ++row)
                    for (unsigned col = 0; col < 2; ++col)
                        multiplexed(value | (row << target), value | (col << target)) = block(row, col);
            }
        check(multiplexed, 3);
    }
    // The copied columns of a controlled unitary are those of its analytic
    // table: the same outputs for every input, with Z or Y on a control too.
    const Gates::CCUGate<> ccu(.43, -.71, .19, .23);
    const auto raw = PauliDetail::CompileMatrix(ccu.getRawOperatorMatrix(), 3);
    const auto typed = PauliDetail::MatrixTransfer(ccu.getRawOperatorMatrix(), 3, Gates::GetUnitaryGateDescriptor(ccu));
    Require(raw->offsets == typed->offsets, "A raw controlled matrix and its analytic table have different outputs");
    for (size_t i = 0; i < raw->entries.size(); ++i)
        Require(raw->entries[i].pauli == typed->entries[i].pauli &&
                    std::abs(raw->entries[i].coefficient - typed->entries[i].coefficient) < 4e-16,
                "A raw controlled matrix differs from its analytic table");
}

// The recording paths that avoid building or hashing a controlled matrix, and
// the caches, which keep more gates than they have sets.
inline void RecordingPaths()
{
    const auto sameAsControlledGate = [](const auto &single) {
        const Gates::ControlledGate<> one(single);
        const Gates::DoublyControlledGate<> two(single);
        const auto descriptor = Gates::GetUnitaryGateDescriptor(single);
        Require(SameTable(*PauliDetail::ControlledTransfer(single.getRawOperatorMatrix(), 1, descriptor),
                          *PauliDetail::MatrixTransfer(one.getRawOperatorMatrix(), 2,
                                                       Gates::GetUnitaryGateDescriptor(one))) &&
                    SameTable(*PauliDetail::ControlledTransfer(single.getRawOperatorMatrix(), 2, descriptor),
                              *PauliDetail::MatrixTransfer(two.getRawOperatorMatrix(), 3,
                                                           Gates::GetUnitaryGateDescriptor(two))),
                "A controlled table differs from the table of the controlled matrix");
    };
    sameAsControlledGate(Gates::RzGate<>(.3));
    sameAsControlledGate(Gates::RxGate<>(M_PI_2));
    sameAsControlledGate(Gates::PhaseShiftGate<>(1e-16));
    sameAsControlledGate(Gates::TGate<>());
    sameAsControlledGate(Gates::HadamardGate<>());
    sameAsControlledGate(Gates::PauliXGate<>());
    sameAsControlledGate(Gates::SquareRootNOTGate<>());
    sameAsControlledGate(Gates::UGate<>(.43, -.71, .19, .23));
    sameAsControlledGate(Gates::AppliedGate<>(Gates::UGate<>(.43, -.71, .19, .23), 0).Adjoint());
    sameAsControlledGate(Gates::AppliedGate<>(Gates::RyGate<>(.37), 0).Adjoint());
    Eigen::MatrixXcd arbitrary(2, 2);
    arbitrary << std::complex<double>(.3, .7), std::complex<double>(-.4, .2), std::complex<double>(.6, -.9),
        std::complex<double>(1.1, .3);
    sameAsControlledGate(Gates::AppliedGate<>(arbitrary, 0));
    Gates::RxGate<> replaced(.3);
    replaced.setOperator(arbitrary);
    sameAsControlledGate(replaced);
    for (unsigned controls : {0U, 3U})
    {
        bool rejected = false;
        try
        {
            PauliDetail::ControlledTransfer(arbitrary, controls);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        Require(rejected, "An unsupported number of controls was accepted");
    }
    bool rejected = false;
    try
    {
        PauliDetail::ControlledTransfer(Eigen::MatrixXcd(Eigen::MatrixXcd::Identity(4, 4)), 1);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    Require(rejected, "A controlled table was built from a matrix that is not 2 by 2");

    // Eight gates in turn, twice: the second round finds every table again,
    // where one entry per cache set would have compiled most of them twice.
    // The keys are exact numbers, so the cache sets they fall into are the
    // same on every platform.
    constexpr unsigned gates = 8;
    std::vector<std::shared_ptr<const PauliDetail::LocalTransfer>> first;
    for (unsigned round = 0; round < 2; ++round)
        for (unsigned gate = 0; gate < gates; ++gate)
        {
            const double angle = 800. + .125 * (gate + 1);
            Eigen::MatrixXcd two = Eigen::MatrixXcd::Identity(4, 4), one = Eigen::MatrixXcd::Identity(2, 2);
            two(3, 3) = one(1, 1) = 1. + angle;
            const std::shared_ptr<const PauliDetail::LocalTransfer> tables[] = {
                PauliDetail::ParameterizedTransfer(OperationType::CRZ, angle),
                PauliDetail::InteractionTransfer('Z', 'Z', angle), PauliDetail::MatrixTransfer(two, 2),
                PauliDetail::ControlledTransfer(one, 2)};
            for (const auto &table : tables)
                if (!round)
                    first.push_back(table);
                else
                    Require(first[4 * gate + size_t(&table - tables)] == table,
                            "A table cache lost a gate of a short repeating sequence");
        }

    // A negative zero is the angle zero: the same cache entry, whichever of
    // the two is asked for first and wherever it stands among the angles.
    const auto plain = PauliDetail::ParameterizedTransfer(OperationType::U, 0., .3125, 0.);
    Require(PauliDetail::ParameterizedTransfer(OperationType::U, -0., .3125, -0.) == plain &&
                PauliDetail::ParameterizedTransfer(OperationType::U, 0., .3125, -0.) == plain &&
                PauliDetail::ParameterizedTransfer(OperationType::U, 0., .3125, 0.) == plain,
            "A negative zero angle has a cache entry of its own");
    const auto negativeFirst = PauliDetail::ParameterizedTransfer(OperationType::CU, -0., -0., .6875, -0.);
    Require(PauliDetail::ParameterizedTransfer(OperationType::CU, 0., 0., .6875, 0.) == negativeFirst &&
                PauliDetail::ParameterizedTransfer(OperationType::CU, -0., 0., .6875, -0.) == negativeFirst,
            "A negative zero angle asked for first has a cache entry of its own");
    Require(SameTable(*PauliDetail::ParameterizedTransfer(OperationType::CRZ, -0.),
                      *PauliDetail::ParameterizedTransfer(OperationType::CRZ, 0.)) &&
                (PauliDetail::ParameterizedTransfer(OperationType::CRZ, -0.)->unchanged & 0xFFFFU) == 0xFFFFU,
            "A controlled rotation by negative zero is not the identity table");

    // Two gates in turn among many others, on every cache: whatever was
    // replaced in between, each call gives the table of its own gate.
    Eigen::MatrixXcd firstMatrix = Eigen::MatrixXcd::Identity(4, 4), secondMatrix = firstMatrix;
    firstMatrix(3, 3) = std::polar(1., .4375);
    secondMatrix(3, 3) = std::polar(1., -1.1875);
    const auto firstTable = PauliDetail::CompileMatrix(firstMatrix, 2),
               secondTable = PauliDetail::CompileMatrix(secondMatrix, 2);
    for (unsigned other = 0; other < 100; ++other)
    {
        Eigen::MatrixXcd otherMatrix = firstMatrix;
        otherMatrix(3, 3) = std::polar(1., 1. + .0625 * other);
        PauliDetail::MatrixTransfer(otherMatrix, 2);
        PauliDetail::ParameterizedTransfer(OperationType::CP, 1. + .0625 * other);
        PauliDetail::InteractionTransfer('Z', 'Z', 1. + .0625 * other);
        const bool odd = other & 1U;
        const double angle = odd ? -1.1875 : .4375;
        Require(SameTable(*PauliDetail::MatrixTransfer(odd ? secondMatrix : firstMatrix, 2),
                          odd ? *secondTable : *firstTable) &&
                    SameTable(*PauliDetail::ParameterizedTransfer(OperationType::CP, angle),
                              *PauliDetail::CompileCP(angle)) &&
                    SameTable(*PauliDetail::InteractionTransfer('Z', 'Z', angle),
                              *PauliDetail::CompilePairRotation('Z', 'Z', angle)),
                "A table cache returned the table of another gate");
    }
}

// Gate objects reach the native kernels when the action is exactly theirs.
inline void NativeDispatch()
{
    using Applied = Gates::AppliedGate<>;
    // Every Pauli on the gate qubits through the recorded operation and
    // through the table of the same matrix: the same terms, bit for bit.
    const auto check = [](const auto &gate, bool expectNative, const char *name) {
        const unsigned arity = unsigned(gate.getQubitsNumber());
        const int qubits[3] = {5, 2, 7};
        PauliPropagator recorded, tabled;
        recorded.SetNrQubits(9);
        tabled.SetNrQubits(9);
        recorded.ApplyGate(gate, qubits[0], qubits[1], qubits[2]);
        const auto table =
            PauliDetail::MatrixTransfer(gate.getRawOperatorMatrix(), arity, Gates::GetUnitaryGateDescriptor(gate));
        std::vector<std::unique_ptr<Operator>> operations;
        operations.push_back(std::make_unique<OperatorLocal>(arity == 1   ? OperationType::LOCAL_1Q
                                                             : arity == 2 ? OperationType::LOCAL_2Q
                                                                          : OperationType::LOCAL_3Q,
                                                             qubits[0], qubits[1], qubits[2], table));
        tabled.SetOperations(std::move(operations));
        const auto recordedOperations = recorded.GetOperations();
        Require(recordedOperations.size() == 1, "A gate object was not recorded as one operation");
        const bool native = !IsPauliLocalGate(recordedOperations.front()->GetType());
        if (native != expectNative)
            throw std::runtime_error(
                std::string(expectNative ? "No native kernel was used for " : "A table was replaced for ") + name);
        for (unsigned label = 0; label < (1U << (2 * arity)); ++label)
        {
            PauliStringXZWithCoefficient pauli(9);
            for (unsigned q = 0; q < arity; ++q)
            {
                pauli.X[size_t(qubits[q])] = (label >> (2 * q)) & 1U;
                pauli.Z[size_t(qubits[q])] = (label >> (2 * q + 1)) & 1U;
            }
            pauli.X[0] = true; // a spectator
            pauli.Coefficient = .37;
            PauliStringStorage a, b;
            recorded.ExpectationValue(pauli, a);
            tabled.ExpectationValue(pauli, b);
            Require(a.size() == b.size(), "A native kernel and its table give different numbers of terms");
            for (const auto &term : a)
            {
                bool found = false;
                for (const auto &other : b)
                    found = found || (term.X == other.X && term.Z == other.Z && term.Coefficient == other.Coefficient);
                if (!found)
                    throw std::runtime_error(std::string("A native kernel differs from its table for ") + name);
            }
        }
    };
    check(Gates::HadamardGate<>(), true, "H");
    check(Gates::SGate<>(), true, "S");
    check(Gates::SDGGate<>(), true, "Sdg");
    check(Gates::PauliXGate<>(), true, "X");
    check(Gates::PauliYGate<>(), true, "Y");
    check(Gates::PauliZGate<>(), true, "Z");
    check(Gates::SquareRootNOTGate<>(), true, "SX");
    check(Gates::SquareRootNOTDagGate<>(), true, "SXdg");
    check(Gates::CNOTGate<>(), true, "CX");
    check(Gates::ControlledYGate<>(), true, "CY");
    check(Gates::ControlledZGate<>(), true, "CZ");
    check(Gates::SwapGate<>(), true, "Swap");
    check(Gates::iSwapGate<>(), true, "iSwap");
    check(Gates::iSwapDagGate<>(), true, "iSwapDag");
    check(Gates::ControlledGate<>(Gates::PauliXGate<>()), true, "controlled X");
    check(Applied(Gates::SGate<>(), 0).Adjoint(), true, "adjoint S");
    for (double angle : {.3, -1.271, 1e-9, M_PI_2, -M_PI_2, M_PI, 3. * M_PI_2, 9. * M_PI_2, std::nextafter(M_PI_2, 2.)})
    {
        check(Gates::RxGate<>(angle), true, "Rx");
        check(Gates::RyGate<>(angle), true, "Ry");
        check(Gates::RzGate<>(angle), true, "Rz");
        check(Gates::PhaseShiftGate<>(angle), true, "Phase");
        check(Applied(Gates::RyGate<>(angle), 0).Adjoint(), true, "adjoint Ry");
    }
    check(Gates::TGate<>(), true, "T");
    check(Gates::TDGGate<>(), true, "Tdg");
    // A matrix that is exactly a native Clifford is recorded as that kernel
    // too, also with its qubits in the other order.
    check(Applied(Gates::HadamardGate<>().getRawOperatorMatrix(), 0), true, "raw H");
    Eigen::MatrixXcd controlOnBitZero = Eigen::MatrixXcd::Zero(4, 4);
    controlOnBitZero(0, 0) = controlOnBitZero(2, 2) = controlOnBitZero(3, 1) = controlOnBitZero(1, 3) = 1.;
    check(Applied(controlOnBitZero, 0, 1), true, "CX with its control on local bit 0");
    // These have no kernel of their own; whatever is recorded must agree.
    check(Gates::UGate<>(.43, -.71, .19, .23), false, "U");
    check(Gates::RZZGate<>(.3), false, "RZZ");
    check(Gates::ECRGate<>(), false, "ECR");
    check(Gates::ControlledRzGate<>(.3), false, "CRz");
    check(Gates::CCUGate<>(.43, -.71, .19, .23), false, "CCU");
    check(Gates::ToffoliGate<>(), false, "CCX");
    check(Applied(Eigen::MatrixXcd((1. + 1e-15) * Gates::HadamardGate<>().getRawOperatorMatrix()), 0), false,
          "scaled H");
    PauliPropagator simulator;
    simulator.SetNrQubits(3);
    simulator.ApplyControlledGate(Gates::PauliZGate<>(), 2, 0);
    simulator.ApplyGate(Gates::RzGate<>(.3), 1);
    simulator.ApplyGate(Gates::HyGate<>(), 1);
    simulator.ApplyGate(Gates::SplitterGate<>(), 1);
    const auto operations = simulator.GetOperations();
    Require(operations.size() == 4 && operations[0]->GetType() == OperationType::CZ &&
                operations[1]->GetType() == OperationType::RZ,
            "A controlled Pauli or a rotation object was not recorded as its native operation");

    // How a table is recorded is remembered with its cache entry and belongs
    // to that table alone. More matrices than the cache holds pass between
    // the uses of a Clifford and of a rotation: whatever takes their entries
    // or moves them is recorded as itself, and so are they when they return.
    PauliPropagator circuit;
    circuit.SetNrQubits(3);
    circuit.SetCoefficientThreshold(0.);
    circuit.SetStepsBetweenDeduplication(1);
    Eigen::VectorXcd state = Eigen::VectorXcd::Zero(8);
    state[0] = 1.;
    const Gates::HadamardGate<> hadamard;
    const Gates::RzGate<> rotation(.3);
    size_t expectedNative = 0, expectedTables = 0;
    for (unsigned step = 0; step < 120; ++step)
    {
        const size_t q = step % 3;
        Eigen::MatrixXcd matrix;
        if (step % 5 == 0)
        {
            matrix = hadamard.getRawOperatorMatrix();
            circuit.ApplyGate(hadamard, int(q));
            ++expectedNative;
        }
        else if (step % 5 == 3)
        {
            matrix = rotation.getRawOperatorMatrix();
            circuit.ApplyGate(rotation, int(q));
            ++expectedNative;
        }
        else
        {
            matrix = Gates::UGate<>(.1 + .37 * step, .2 * step, -.3 + .11 * step).getRawOperatorMatrix();
            circuit.ApplyGate(Applied(matrix, q));
            ++expectedTables;
        }
        state = Embed(matrix, {q}, 3) * state;
    }
    size_t recordedNative = 0, recordedTables = 0;
    for (const auto &operation : circuit.GetOperations())
        ++(IsPauliLocalGate(operation->GetType()) ? recordedTables : recordedNative);
    Require(recordedNative == expectedNative && recordedTables == expectedTables,
            "A gate was recorded the way another gate of its cache entry had been");
    CheckExpectations(circuit, state, "A circuit of more matrices than the cache holds differs from its dense state");
}

// Compact gates on at most three qubits follow the matrix path.
inline void CompactGates()
{
    using Applied = Gates::AppliedGate<>;
    const std::complex<double> i(0., 1.);
    std::vector<Applied> gates;
    gates.push_back(Applied::Diagonal({std::exp(.3 * i), std::exp(-.7 * i)}, {1}));
    gates.push_back(Applied::Diagonal({1., std::exp(.3 * i), std::exp(-.2 * i), std::exp(.9 * i)}, {2, 0}));
    gates.push_back(Applied::Diagonal({1., std::exp(.1 * i), std::exp(.2 * i), std::exp(.3 * i), std::exp(.4 * i),
                                       std::exp(.5 * i), std::exp(.6 * i), -1.},
                                      {1, 2, 0}));
    const Gates::UGate<> u(.43, -.71, .19, .23);
    gates.push_back(Applied::MultiControlled(u, 2, {}));
    gates.push_back(Applied::MultiControlled(u, 2, {0}));
    gates.push_back(Applied::MultiControlled(u, 0, {2}, 0));
    gates.push_back(Applied::MultiControlled(u, 1, {2, 0}));
    gates.push_back(Applied::MultiControlled(u, 1, {0, 2}, 1));
    gates.push_back(Applied::MultiControlled(Gates::PauliXGate<>(), 0, {1}));
    gates.push_back(Applied::MultiControlled(Gates::PauliXGate<>(), 2, {0, 1}));
    for (auto axis : {Gates::RotationAxis::X, Gates::RotationAxis::Y, Gates::RotationAxis::Z})
    {
        gates.push_back(Applied::MultiplexedRotation(axis, {.37}, 1, {}));
        gates.push_back(Applied::MultiplexedRotation(axis, {.37, -.91}, 2, {0}));
        gates.push_back(Applied::MultiplexedRotation(axis, {.37, -.91, 1.3, .05}, 0, {2, 1}));
    }
    for (const auto &gate : gates)
    {
        const auto *compact = gate.getStructuredGate();
        Require(compact != nullptr, "A compact applied gate lost its compact form");
        Eigen::VectorXcd state = Eigen::VectorXcd::Zero(8);
        state[0] = 1.;
        const auto prepare = [&](PauliPropagator &simulator) {
            simulator.SetNrQubits(3);
            simulator.SetCoefficientThreshold(0.);
            for (int q = 0; q < 3; ++q)
                simulator.ApplyGate(Gates::UGate<>(.37 + q * .21, -.29 + q * .13, .19 - q * .17), q);
        };
        for (size_t q = 0; q < 3; ++q)
            state = Embed(Gates::UGate<>(.37 + q * .21, -.29 + q * .13, .19 - q * .17).getRawOperatorMatrix(), {q}, 3) *
                    state;
        const auto qubits = compact->OrderedQubits();
        state = Embed(compact->DenseMatrix<Eigen::MatrixXcd>(), qubits, 3) * state;
        PauliPropagator own, base, explicitQubits, structured;
        prepare(own);
        own.ApplyGate(gate);
        CheckExpectations(own, state, "A compact gate on its own qubits differs from its dense matrix");
        prepare(base);
        const Gates::QuantumGateWithOp<> &reference = gate;
        base.ApplyGate(reference, int(qubits[0]), qubits.size() > 1 ? int(qubits[1]) : 0,
                       qubits.size() > 2 ? int(qubits[2]) : 0);
        CheckExpectations(base, state, "A compact gate behind a base reference differs from its dense matrix");
        prepare(explicitQubits);
        explicitQubits.ApplyGate(gate, int(qubits[0]), qubits.size() > 1 ? int(qubits[1]) : 0,
                                 qubits.size() > 2 ? int(qubits[2]) : 0);
        CheckExpectations(explicitQubits, state, "A compact gate on explicit qubits differs from its dense matrix");
        prepare(structured);
        structured.ApplyStructuredGate(*compact);
        CheckExpectations(structured, state, "A structured gate differs from its dense matrix");
    }
    // A compact controlled X on two qubits is a CX: one native operation.
    PauliPropagator simulator;
    simulator.SetNrQubits(5);
    simulator.ApplyGate(Applied::MultiControlled(Gates::PauliXGate<>(), 0, {3}));
    Require(simulator.GetOperations().size() == 1 && simulator.GetOperations().front()->GetType() == OperationType::CX,
            "A compact controlled X was not recorded as CX");
    // More than three qubits has no local table: rejected, nothing recorded.
    const auto wide = Applied::MultiControlled(Gates::PauliZGate<>(), 0, {1, 2, 3});
    const Gates::QuantumGateWithOp<> &wideReference = wide;
    for (unsigned entry = 0; entry < 4; ++entry)
    {
        bool rejected = false;
        try
        {
            if (entry == 0)
                simulator.ApplyGate(wide);
            else if (entry == 1)
                simulator.ApplyGate(wide, 0, 1, 2);
            else if (entry == 2)
                simulator.ApplyGate(wideReference, 0, 1, 2);
            else
                simulator.ApplyStructuredGate(*wide.getStructuredGate());
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        Require(rejected && simulator.GetOperations().size() == 1, "A compact gate on four qubits was not rejected");
    }
}

// Outcome probabilities take <I> = 1; a matrix that is not unitary breaks
// that, so they are refused while expectation values stay available.
inline void NonunitaryProbabilities()
{
    using Applied = Gates::AppliedGate<>;
    Eigen::MatrixXcd filter = Eigen::MatrixXcd::Identity(2, 2);
    filter(1, 1) = .5;
    const auto refused = [](PauliPropagator &simulator) {
        unsigned count = 0;
        for (unsigned entry = 0; entry < 4; ++entry)
            try
            {
                if (entry == 0)
                    simulator.Probability(0);
                else if (entry == 1)
                    simulator.Sample({0, 1});
                else if (entry == 2)
                    simulator.SampleCounts({0, 1}, 32);
                else
                    simulator.Measure({0});
            }
            catch (const std::domain_error &)
            {
                ++count;
            }
        return count;
    };
    PauliPropagator simulator;
    simulator.SetNrQubits(2);
    simulator.SetSeed(7);
    simulator.ApplyH(0);
    simulator.ApplyCX(0, 1);
    Require(refused(simulator) == 0, "Probabilities of a unitary circuit were refused");
    simulator.ClearOperations();
    simulator.ApplyH(0);
    simulator.ApplyCX(0, 1);
    // A unitary given as a matrix is unitary only up to rounding: still fine.
    simulator.ApplyGate(Applied(Gates::UGate<>(.43, -.71, .19, .23).getRawOperatorMatrix(), 1));
    simulator.SaveState();
    Require(std::abs(simulator.Probability(0) + simulator.Probability(1) + simulator.Probability(2) +
                     simulator.Probability(3) - 1.) < 1e-12,
            "Probabilities after a raw unitary matrix are not normalised");
    const size_t operations = simulator.GetOperations().size();
    simulator.ApplyGate(Applied(filter, 0));
    Require(std::abs(simulator.ExpectationValue("II") - .625) < 1e-15,
            "A raw filter lost its unnormalised expectation");
    Require(refused(simulator) == 4 && simulator.GetOperations().size() == operations + 1,
            "Probabilities after a matrix that is not unitary were not refused");
    PauliPropagator imported, shared;
    imported.SetNrQubits(2);
    imported.SetOperations(simulator.GetOperations());
    shared.SetNrQubits(2);
    shared.ShareOperationsFrom(simulator);
    Require(refused(imported) == 4 && refused(shared) == 4, "A copied circuit forgot its matrix that is not unitary");
    simulator.RestoreState();
    Require(refused(simulator) == 0, "Removing the matrix that is not unitary did not restore probabilities");
    simulator.ApplyControlledGate(Eigen::MatrixXcd(2. * Eigen::MatrixXcd::Identity(2, 2)), 0, 1);
    Require(refused(simulator) == 4, "A controlled matrix that is not unitary was not noticed");
    simulator.ClearOperations();
    Require(refused(simulator) == 0, "Clearing the circuit did not restore probabilities");
}

// A raw matrix of the form a*I + b*P (a rotation about a Pauli product given
// as a plain matrix) gets its table from a and b. It must be the table of the
// general compilation up to rounding, for every Pauli that moves a bit, for
// unitary and non-unitary pairs, and only for matrices of exactly that form.
inline void RawPauliPairs()
{
    using Complex = std::complex<double>;
    const std::vector<std::pair<Complex, Complex>> pairs{
        {Complex(std::cos(.37), 0.), Complex(0., -std::sin(.37))}, // exp(-i t P)
        {Complex(std::cos(1.1) * .6, std::cos(1.1) * .8), Complex(std::sin(1.1) * .8, -std::sin(1.1) * .6)},
        {Complex(.3, -.7), Complex(.2, .45)},                       // not unitary
        {Complex(1., 0.), Complex(0., 1e-200)},                     // a tiny genuine term
        {Complex(0., 0.), Complex(.6, .8)},                         // a Pauli with a phase
        {Complex(std::sqrt(.5), 0.), Complex(0., -std::sqrt(.5))}}; // the quarter turn
    for (unsigned qubits = 1; qubits <= 3; ++qubits)
    {
        const unsigned dimension = 1U << qubits, count = dimension * dimension;
        for (unsigned code = 1; code < count; ++code)
        {
            const Eigen::MatrixXcd pauli = Pauli(code, qubits);
            const bool moves = !pauli.isDiagonal();
            for (const auto &[a, b] : pairs)
            {
                Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(dimension, dimension);
                for (unsigned row = 0; row < dimension; ++row)
                    for (unsigned col = 0; col < dimension; ++col)
                        matrix(row, col) = row == col ? a + b * pauli(row, col) : b * pauli(row, col);
                const auto pair = PauliDetail::CompilePauliPair(matrix, qubits);
                if (!moves)
                {
                    Require(!pair, "A diagonal matrix was taken for a rotation about a Pauli that moves a bit");
                    continue;
                }
                Require(pair != nullptr, "A matrix of the form a I + b P was not recognized");
                const auto general = PauliDetail::CompileMatrix(matrix, qubits);
                const double scale = std::max(1., std::norm(a) + std::norm(b));
                for (unsigned input = 0; input < count; ++input)
                    for (unsigned output = 0; output < count; ++output)
                        Require(std::abs(Coefficient(*pair, input, output) - Coefficient(*general, input, output)) <=
                                    8e-16 * scale,
                                "The table of a I + b P differs from the general compilation");
                Require(pair->nonunitary == general->nonunitary, "The two compilations disagree on unitarity");
                CheckTransfer(matrix, qubits, pair);
                // A matrix one rounding away from the form takes the general path.
                Eigen::MatrixXcd nearby = matrix;
                nearby(dimension - 1, dimension - 1) = Complex(std::nextafter(a.real(), 2.), a.imag());
                Require(!PauliDetail::CompilePauliPair(nearby, qubits),
                        "A matrix that is not exactly a I + b P took its shortcut");
            }
        }
    }
    // The public path: a raw rotation matrix leaves the identity alone and equals its analytic table.
    const auto check = [](const auto &gate, unsigned qubits) {
        const Eigen::MatrixXcd matrix = gate.getRawOperatorMatrix();
        const auto raw = PauliDetail::MatrixTransfer(matrix, qubits);
        const auto typed = PauliDetail::MatrixTransfer(matrix, qubits, Gates::GetUnitaryGateDescriptor(gate));
        Require(raw->KeepsIdentitySupport() && std::abs(Coefficient(*raw, 0, 0) - 1.) < 4e-16 && !raw->nonunitary &&
                    raw->offsets == typed->offsets,
                "A raw rotation matrix branches or rescales the identity");
        for (size_t i = 0; i < raw->entries.size(); ++i)
            Require(raw->entries[i].pauli == typed->entries[i].pauli &&
                        std::abs(raw->entries[i].coefficient - typed->entries[i].coefficient) < 4e-16,
                    "A raw rotation matrix differs from its analytic transfer");
    };
    check(Gates::RxGate<>(.41), 1);
    check(Gates::RyGate<>(-1.3), 1);
    check(Gates::RXXGate<>(.83), 2);
    check(Gates::RYYGate<>(-.29), 2);
    check(Gates::RZXGate<>(2.2), 2);
    // The tiny term survives: nothing of it is rounding.
    Eigen::MatrixXcd tiny = Eigen::MatrixXcd::Identity(2, 2);
    tiny(0, 1) = tiny(1, 0) = Complex(0., 1e-200);
    const auto table = PauliDetail::CompilePauliPair(tiny, 1);
    Require(table && Coefficient(*table, 2, 3) != 0., "A tiny genuine coefficient of a I + b P was dropped");
}

inline void Run()
{
    ControlledBlocks();
    Interactions();
    IdentityGrowthAndCache();
    TypedUnitaryPrecision();
    CoherentEulerAndAdjoints();
    MalformedMatrixRejection();
    ExactInverseDescriptors();
    CrossFamilyInverseDescriptors();
    MatrixKeyPrefixCopies();
    RawMatrixResidue();
    RawPauliPairs();
    FixedQubitColumns();
    RecordingPaths();
    NativeDispatch();
    CompactGates();
    NonunitaryProbabilities();
    ExtendedFilters();
    std::cout << "Pauli/extended gate precision, sparse transfers and filter regressions passed.\n";
}
} // namespace PauliGatePrecisionRegression
} // namespace QC
