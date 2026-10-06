#pragma once

#include "QuantumGate.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace QC
{
namespace TensorNetworkRegression
{

using Matrix = Eigen::MatrixXcd;
using Complex = std::complex<double>;
using Gate = Gates::AppliedGate<>;

inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline void Close(Complex actual, Complex expected, const char *message)
{
    Require(std::isfinite(actual.real()) && std::isfinite(actual.imag()) && std::abs(actual - expected) < 2E-10,
            message);
}

inline void Close(const Matrix &actual, const Matrix &expected, const char *message)
{
    Require(actual.rows() == expected.rows() && actual.cols() == expected.cols() && actual.allFinite() &&
                (actual - expected).cwiseAbs().maxCoeff() < 2E-10,
            message);
}

template <class Exception, class Func> void Throws(Func &&operation, const char *message)
{
    try
    {
        operation();
    }
    catch (const Exception &)
    {
        return;
    }
    throw std::runtime_error(message);
}

inline Matrix Unitary(int size, unsigned int seed)
{
    std::mt19937 random(seed);
    std::normal_distribution<double> distribution;
    Matrix matrix(size, size);
    for (int column = 0; column < size; ++column)
        for (int row = 0; row < size; ++row)
            matrix(row, column) = Complex(distribution(random), distribution(random));
    return matrix.householderQr().householderQ() * Matrix::Identity(size, size);
}

inline Matrix Pauli(char p)
{
    Matrix matrix = Matrix::Zero(2, 2);
    switch (std::toupper(static_cast<unsigned char>(p)))
    {
    case 'I':
        matrix.setIdentity();
        break;
    case 'X':
        matrix(0, 1) = matrix(1, 0) = 1.;
        break;
    case 'Y':
        matrix(0, 1) = Complex(0., -1.);
        matrix(1, 0) = Complex(0., 1.);
        break;
    case 'Z':
        matrix(0, 0) = 1.;
        matrix(1, 1) = -1.;
        break;
    default:
        throw std::invalid_argument("Invalid Pauli in test reference");
    }
    return matrix;
}

inline std::vector<Gate> PauliGates(const std::string &pauli)
{
    std::vector<Gate> gates;
    for (size_t q = 0; q < pauli.size(); ++q)
        if (std::toupper(static_cast<unsigned char>(pauli[q])) != 'I')
            gates.emplace_back(Pauli(pauli[q]), q);
    return gates;
}

inline Complex DensePauliExpectation(const Matrix &density, const std::string &pauli)
{
    Complex result = 0.;
    for (Eigen::Index input = 0; input < density.rows(); ++input)
    {
        Eigen::Index output = input;
        Complex phase = 1.;
        for (size_t q = 0; q < pauli.size(); ++q)
        {
            const Eigen::Index mask = Eigen::Index(1) << q;
            const bool one = (input & mask) != 0;
            switch (std::toupper(static_cast<unsigned char>(pauli[q])))
            {
            case 'X':
                output ^= mask;
                break;
            case 'Y':
                output ^= mask;
                phase *= Complex(0., one ? -1. : 1.);
                break;
            case 'Z':
                if (one)
                    phase = -phase;
                break;
            }
        }
        result += phase * density(input, output);
    }
    return result;
}

template <class Sim> void Prepare(Sim &simulator, int depth = 3)
{
    simulator.SetMultithreading(false);
    const Gate one(Unitary(2, 19)), two(Unitary(4, 29));
    const auto count = static_cast<Eigen::Index>(simulator.getNrQubits());
    for (int layer = 0; layer < depth; ++layer)
    {
        for (Eigen::Index q = 0; q < count; ++q)
            simulator.ApplyGate(one, q);
        for (Eigen::Index q = layer % 2; q + 1 < count; q += 2)
            simulator.ApplyGate(two, q + 1, q);
    }
}

template <class Sim, class SingleExpectation> void CheckObservableBatch(SingleExpectation singleExpectation)
{
    Sim simulator(8, 1);
    simulator.SetSeed(0x4241544348ULL);
    simulator.setLimitBondDimension(8);
    Prepare(simulator, 4);
    const std::vector<std::string> paulis = {"IIIIIIII", "YIIIIIII", "IIZIIIII", "XXYZIIII", "XXYYIIII",
                                             "XXYYIIII", "XXYYYYII", "XXYYYYIX", "XXXYYYYZ", "XXYIIIII",
                                             "IIIIIIXY", "IIIIIIYZ", "iixyziii", "iiiiiiii"};
    const auto check = [&] {
        const auto values = simulator.ExpectationValues(paulis);
        Require(values.size() == paulis.size(), "Batch changed the number of results");
        for (size_t i = 0; i < paulis.size(); ++i)
            Close(values[i], singleExpectation(simulator, paulis[i]), "Batch differs from a single observable");
    };
    check();
    simulator.MoveAtBeginningOfChain({7, 2});
    check();
    simulator.SaveState();
    simulator.ApplyGate(Gate(Unitary(2, 67)), 3);
    check();
    simulator.RestoreState();
    check();
    simulator.ReCanonicalize();
    check();
    auto clone = simulator.Clone();
    const auto cloned = clone->ExpectationValues(paulis);
    const auto original = simulator.ExpectationValues(paulis);
    for (size_t i = 0; i < original.size(); ++i)
        Close(cloned[i], original[i], "Clone changed batch observables");
    Require(simulator.ExpectationValues({}).empty(), "Empty batch returned results");
    for (const std::string invalid : {"I", "IIIIIIIII", "IIIIIII?"})
        Throws<std::invalid_argument>(
            [&] {
                simulator.ExpectationValues({"IIIIIIII", invalid});
            },
            "Invalid batch Pauli was accepted");
    check();
}

template <class Sim> Eigen::Index MaximumBond(const Sim &simulator)
{
    Eigen::Index maximum = 1;
    for (const auto dimension : simulator.getBondDimensions())
        maximum = std::max(maximum, dimension);
    return maximum;
}

template <class Sim> void CheckBondSummary()
{
    Sim simulator(8, 1);
    simulator.SetSeed(0x53554D4D415259ULL);
    simulator.SetMultithreading(false);
    Eigen::Index lastSnapshot = 1;
    int fullCalls = 0, summaryCalls = 0;
    simulator.SetBondDimensionCallback([&](const auto &dimensions) {
        lastSnapshot = 1;
        for (auto dimension : dimensions)
            lastSnapshot = std::max(lastSnapshot, dimension);
        ++fullCalls;
    });
    simulator.SetBondDimensionSummaryCallback([&](auto maximum) {
        Require(maximum == lastSnapshot && maximum == MaximumBond(simulator),
                "Summary callback differs from actual bond dimensions");
        ++summaryCalls;
    });
    const auto check = [&] {
        Require(simulator.getMaxBondDimension() == MaximumBond(simulator), "Maximum bond getter is stale");
        Require(fullCalls == summaryCalls, "Summary notification boundaries differ from full snapshots");
    };
    check();
    Prepare(simulator);
    check();
    Require(summaryCalls > 0 && simulator.getMaxBondDimension() > 1, "Summary test did not grow the bonds");
    simulator.MoveAtBeginningOfChain({7, 3});
    check();
    simulator.SaveState();
    simulator.setLimitBondDimension(2);
    simulator.Trim();
    check();
    simulator.ReCanonicalize();
    check();
    simulator.RestoreState();
    check();
    auto clone = simulator.Clone();
    clone->SetBondDimensionCallback(nullptr);
    clone->SetBondDimensionSummaryCallback([](auto) {});
    Require(clone->getMaxBondDimension() == simulator.getMaxBondDimension(), "Clone changed maximum bond");
    clone->Clear();
    Require(clone->getMaxBondDimension() == 1, "Clone clear retained stale dimensions");
    clone->RestoreState();
    Require(clone->getMaxBondDimension() == simulator.getMaxBondDimension(), "Clone restore lost summary state");
    check();
    auto state = simulator.getState();
    simulator.Clear();
    check();
    simulator.setState(state);
    check();
    simulator.Clear();
    simulator.setStateDestructive(state);
    check();
    simulator.MeasureQubit(3);
    check();
    simulator.InitOnesState();
    check();
    simulator.SetBondDimensionSummaryCallback(nullptr);
    simulator.RestoreState();
    Require(simulator.getMaxBondDimension() == MaximumBond(simulator) && simulator.getMaxBondDimension() > 1,
            "Disabled summary fallback lost entangled bond dimensions");
    Sim single(1, 1);
    single.SetBondDimensionSummaryCallback([](auto) {});
    Require(single.getMaxBondDimension() == 1, "One-site maximum bond is not one");
}

} // namespace TensorNetworkRegression
} // namespace QC
