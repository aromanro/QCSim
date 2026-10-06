#pragma once

#include "MPOSimulator.h"
#include "QuantumGate.h"

#include <iostream>
#include <numeric>
#include <type_traits>

namespace QC
{
namespace MPORegression
{

using Simulator = TensorNetworks::MPOSimulator;
using Implementation = TensorNetworks::MPOSimulatorImpl;
using Interface = TensorNetworks::MPOSimulatorInterface;
using Matrix = Interface::MatrixClass;
using Index = Interface::IndexType;

inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline void Close(const Matrix &actual, const Matrix &expected, const char *message)
{
    Require(actual.rows() == expected.rows() && actual.cols() == expected.cols() && actual.allFinite() &&
                (actual - expected).norm() < 1E-10,
            message);
}

inline void Close(std::complex<double> actual, std::complex<double> expected, const char *message)
{
    Require(std::isfinite(actual.real()) && std::isfinite(actual.imag()) && std::abs(actual - expected) < 1E-10,
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

template <class Sim> void CheckStaleEnvironments()
{
    Sim mpo(3);
    Gates::CNOTGate<> cnot;
    mpo.ApplyGate(Gates::RyGate<>(2. * std::asin(std::sqrt(0.1))), 0);
    mpo.ApplyGate(cnot, 1, 0);
    mpo.ApplyGate(cnot, 2, 1);
    const Matrix expected = mpo.getUnnormalizedDensityMatrix();
    Matrix first = Matrix::Identity(2, 2), second = first;
    first(1, 1) = std::sqrt(0.1);
    second(0, 0) = std::sqrt(0.1);
    // On sqrt(.9)|000> + sqrt(.1)|111>, each pair multiplies both
    // amplitudes by the same scalar. Conditional normalization cancels it.
    for (int i = 0; i < 40; ++i)
    {
        mpo.ApplyOperatorAndNormalize(Gates::SingleQubitGate<>(first), 1);
        mpo.ApplyOperatorAndNormalize(Gates::SingleQubitGate<>(second), 0);
    }
    Close(mpo.getUnnormalizedDensityMatrix(), expected, "Filter pairs changed the reference state");
    const auto snapshot = mpo.getState();
    mpo.ApplyGate(Gates::TwoQubitsGate<>(Matrix::Identity(4, 4)), 1, 2);
    Close(mpo.getUnnormalizedDensityMatrix(), expected, "Identity gate discarded a significant component");
    mpo.setState(snapshot);
    mpo.ReCanonicalize();
    Close(mpo.getUnnormalizedDensityMatrix(), expected, "Gauge repair discarded a significant component");
    mpo.setState(snapshot);
    mpo.ApplyKrausOperators(std::vector<Matrix>{Matrix::Identity(4, 4)}, 1, 2);
    Close(mpo.getUnnormalizedDensityMatrix(), expected, "Identity channel discarded a significant component");
    Close(mpo.Trace(), 1., "Exact evolution changed the trace");
}

inline void StaleEnvironments()
{
    CheckStaleEnvironments<Implementation>();
    CheckStaleEnvironments<Simulator>();
}

inline void WideSampling()
{
    Simulator mpo(80);
    mpo.SetSeed(0x4D504F01ULL);
    for (Index q = 0; q < 72; ++q)
        mpo.ApplyGate(Gates::HadamardGate<>(), q);
    mpo.ApplyGate(Gates::PauliXGate<>(), 78);
    int zeroCount = 0;
    for (int shot = 0; shot < 128; ++shot)
    {
        const auto sample = mpo.MeasureNoCollapse();
        Require(sample.size() == 80, "Full sampling returned an incomplete outcome");
        for (Index q = 72; q < 80; ++q)
            Require(sample.at(q) == (q == 78), "A deterministic tail bit was sampled incorrectly");
        zeroCount += !sample.at(60);
        const auto subset = mpo.MeasureNoCollapse({0, 60, 78, 79});
        Require(subset.size() == 4 && subset.at(78) && !subset.at(79), "Subset sampling biased the tail");
    }
    Require(zeroCount > 35 && zeroCount < 93, "High-index unbiased bit became deterministic");
    Close(mpo.GetProbability(60), 0.5, "Sampling collapsed the state");
    Close(mpo.Trace(), 1., "Sampling changed the trace");

    Simulator correlated(70);
    correlated.SetSeed(0x4D504F02ULL);
    for (Index q = 0; q < 62; ++q)
        correlated.ApplyGate(Gates::HadamardGate<>(), q);
    correlated.ApplyGate(Gates::CNOTGate<>(), 69, 0);
    for (int shot = 0; shot < 64; ++shot)
    {
        const auto sample = correlated.MeasureNoCollapse({0, 65, 69});
        Require(sample.at(0) == sample.at(69) && !sample.at(65), "Routed subset lost a long-range correlation");
    }

    Simulator scaled(2);
    scaled.ApplyOperator(Gates::SingleQubitGate<>(1E-100 * Matrix::Identity(2, 2)), 0);
    const auto sample = scaled.MeasureNoCollapse();
    Require(sample.size() == 2 && !sample.at(0) && !sample.at(1), "Sampling depends on the operator's overall scale");
    scaled.ApplyOperator(Gates::SingleQubitGate<>(Matrix::Zero(2, 2)), 0);
    Throws<std::runtime_error>([&] { scaled.MeasureNoCollapse(); }, "Zero operator was sampled without an error");
}

inline void LogicalOverlaps()
{
    Simulator a(3), b(3);
    a.setToBasisState(size_t{1});
    b.setToBasisState(size_t{1});
    a.MoveAtBeginningOfChain({2});
    Close(a.HilbertSchmidtOverlap(a), 1., "Reordered pure state has nonunit self-overlap");
    Close(a.HilbertSchmidtOverlap(b), 1., "Overlap ignored the caller's logical mapping");
    Close(b.HilbertSchmidtOverlap(a), 1., "Overlap ignored the other logical mapping");

    Simulator wide(24), other(24);
    Implementation physical(24);
    for (Interface *sim : std::vector<Interface *>{&wide, &other, &physical})
        sim->ApplyBitFlipNoise(0, 0.25);
    wide.MoveAtBeginningOfChain({23});
    wide.setLimitBondDimension(1);
    wide.setLimitEntanglement(0.9);
    const auto originalMap = wide.getQubitsMap();
    int callbacks = 0;
    wide.SetBondDimensionCallback([&](const auto &) { ++callbacks; });
    for (const Interface *sim : std::vector<const Interface *>{&wide, &other, &physical})
    {
        Close(wide.HilbertSchmidtOverlap(*sim), 0.625, "Large-register overlap used dense storage or truncation");
        Close(sim->HilbertSchmidtOverlap(wide), 0.625, "Reverse large-register overlap failed");
    }
    Require(callbacks == 0 && wide.getQubitsMap() == originalMap, "Overlap mutated routing or invoked callbacks");

    Simulator mixedA(4), mixedB(4);
    for (Index q = 0; q < 4; ++q)
    {
        mixedA.ApplyGate(Gates::RyGate<>(0.3 * (q + 1)), q);
        mixedB.ApplyGate(Gates::RxGate<>(0.2 * (q + 1)), q);
    }
    mixedA.ApplyGate(Gates::CNOTGate<>(), 3, 0);
    mixedB.ApplyGate(Gates::CNOTGate<>(), 0, 2);
    mixedA.ApplyAmplitudeDamping(0, 0.3);
    mixedB.ApplyDepolarizingNoise(2, 0.2);
    const Matrix rhoA = mixedA.getDensityMatrix(), rhoB = mixedB.getDensityMatrix();
    const auto expected = (rhoA.adjoint() * rhoB).trace();
    Close(mixedA.HilbertSchmidtOverlap(mixedB), expected, "Routed mixed-state overlap differs from dense reference");
    Close(mixedB.HilbertSchmidtOverlap(mixedA), std::conj(expected), "Overlap lost conjugate symmetry");
    Close(mixedA.getDensityMatrix(), rhoA, "Overlap changed its caller");
    Close(mixedB.getDensityMatrix(), rhoB, "Overlap changed its argument");

    Implementation complexA(1), complexB(1);
    auto state = std::dynamic_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(complexA.getState());
    state->gammas[0](0, 0, 1, 0) = {0., 0.25};
    complexA.setState(state);
    state->gammas[0](0, 0, 1, 0) = 0.5;
    complexB.setState(state);
    Close(complexA.HilbertSchmidtOverlap(complexB), {1., -0.125}, "Overlap conjugated the wrong operand");
}

inline void AtomicRoutingFailures()
{
    Simulator mpo(3);
    mpo.ApplyGate(Gates::HadamardGate<>(), 0);
    mpo.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    mpo.setLimitBondDimension(1);
    mpo.setKrausCompletenessCheck(Interface::KrausCompletenessCheck::Strict);
    const Matrix original = mpo.getUnnormalizedDensityMatrix();
    const auto originalMap = mpo.getQubitsMap();
    const auto originalBonds = mpo.getBondDimensions();
    auto unchanged = [&]() {
        Close(mpo.getUnnormalizedDensityMatrix(), original, "Rejected routed operation changed the operator");
        Require(mpo.getQubitsMap() == originalMap && mpo.getBondDimensions() == originalBonds,
                "Rejected routed operation changed mappings or bonds");
    };
    const Matrix incomplete = 0.5 * Matrix::Identity(4, 4);
    Throws<std::invalid_argument>([&] { mpo.ApplyKrausOperators(std::vector<Matrix>{incomplete}, 0, 2); },
                                  "Strict matrix Kraus channel was not rejected");
    unchanged();
    Throws<std::invalid_argument>(
        [&] { mpo.ApplyKrausOperators(std::vector<Gates::AppliedGate<>>{Gates::AppliedGate<>(incomplete, 0, 2)}); },
        "Strict applied Kraus channel was not rejected");
    unchanged();
    Throws<std::runtime_error>([&] { mpo.ApplyOperatorAndNormalize(Gates::TwoQubitsGate<>(Matrix::Zero(4, 4)), 0, 2); },
                               "Impossible normalized operator was not rejected");
    unchanged();
    Throws<std::runtime_error>([&] { mpo.ApplyOperatorAndNormalize(Gates::AppliedGate<>(Matrix::Zero(4, 4), 0, 2)); },
                               "Impossible normalized applied operator was not rejected");
    unchanged();
}

inline void InvalidBasisState()
{
    Simulator mpo(3);
    mpo.setToBasisState(size_t{1});
    mpo.MoveAtBeginningOfChain({2});
    const Matrix original = mpo.getDensityMatrix();
    const auto originalMap = mpo.getQubitsMap();
    Throws<std::invalid_argument>([&] { mpo.setToBasisState(size_t{8}); },
                                  "Out-of-range integer basis state was accepted");
    Close(mpo.getDensityMatrix(), original, "Invalid basis state changed the operator");
    Require(mpo.getQubitsMap() == originalMap, "Invalid basis state reset the mapping");
    Implementation physical(3);
    physical.setToBasisState(size_t{1});
    Throws<std::invalid_argument>([&] { physical.setToBasisState(size_t{8}); },
                                  "Implementation accepted invalid basis state");
    Close(physical.getBasisStateProbability(size_t{1}), 1.,
          "Implementation changed state after invalid initialization");
}

template <class Sim> void CheckDenseBounds()
{
    Sim mpo(64);
    std::vector<Index> all(64);
    std::iota(all.begin(), all.end(), Index{0});
    Throws<std::runtime_error>([&] { mpo.PartialTrace(all); }, "Dense partial trace accepted 64 retained qubits");
    Throws<std::invalid_argument>([&] { mpo.FidelityWithStatevector(Eigen::VectorXcd::Ones(1)); },
                                  "Wide register accepted an unrelated one-element statevector");
    all.resize(14);
    Throws<std::runtime_error>([&] { mpo.PartialTrace(all); }, "Partial trace exceeded the dense allocation limit");
    std::vector<bool> basis(64, false);
    basis[63] = true;
    mpo.setToBasisState(basis);
    Matrix expected = Matrix::Zero(4, 4);
    expected(1, 1) = 1.;
    Close(mpo.PartialTrace({63, 0}), expected, "Wide register's small partial trace has incorrect bit order");
}

inline void DenseBounds()
{
    CheckDenseBounds<Simulator>();
    CheckDenseBounds<Implementation>();
    Simulator shortRegister(3);
    shortRegister.setToBasisState(size_t{1});
    shortRegister.MoveAtBeginningOfChain({2});
    Eigen::VectorXcd psi = Eigen::VectorXcd::Zero(8);
    psi[1] = 1.;
    Close(shortRegister.FidelityWithStatevector(psi), 1., "Fidelity lost the logical bit mapping");
}

inline void ThresholdTrim()
{
    for (auto mode : {Interface::TruncationMode::RelativeToMax, Interface::TruncationMode::DiscardedWeight})
        for (bool cap : {false, true})
        {
            Simulator mpo(2);
            mpo.ApplyGate(Gates::RyGate<>(0.2), 0);
            mpo.ApplyGate(Gates::CNOTGate<>(), 1, 0);
            mpo.setTruncationMode(mode);
            mpo.setLimitEntanglement(mode == Interface::TruncationMode::RelativeToMax ? 0.2 : 0.03);
            if (cap)
                mpo.setLimitBondDimension(100);
            mpo.Trim();
            Require(mpo.getBondDimensions()[0] == 1, "Trim ignored a threshold without an exceeded bond cap");
            Matrix expected = Matrix::Zero(4, 4);
            expected(0, 0) = std::pow(std::cos(0.1), 2);
            Close(mpo.getUnnormalizedDensityMatrix(), expected,
                  "Threshold Trim differs from the analytic rank-one approximation");
        }
}

inline Matrix DenseApply(const Matrix &rho, const Matrix &op, size_t target, size_t control = 0)
{
    Matrix lifted = Matrix::Zero(rho.rows(), rho.cols());
    const size_t mask = (size_t{1} << target) | (op.rows() == 4 ? size_t{1} << control : 0);
    for (size_t input = 0; input < static_cast<size_t>(rho.rows()); ++input)
    {
        const size_t local = ((input >> target) & 1) | (op.rows() == 4 ? ((input >> control) & 1) << 1 : 0);
        for (size_t output = 0; output < static_cast<size_t>(op.rows()); ++output)
        {
            const size_t row =
                (input & ~mask) | ((output & 1) << target) | (op.rows() == 4 ? (output >> 1) << control : 0);
            lifted(row, input) = op(output, local);
        }
    }
    return lifted * rho * lifted.adjoint();
}

inline Matrix DensePartial(const Matrix &rho, const std::vector<Index> &keep)
{
    const size_t dim = size_t{1} << keep.size();
    Matrix result = Matrix::Zero(dim, dim);
    size_t mask = 0;
    for (Index q : keep)
        mask |= size_t{1} << q;
    for (size_t row = 0; row < static_cast<size_t>(rho.rows()); ++row)
        for (size_t col = 0; col < static_cast<size_t>(rho.cols()); ++col)
        {
            if ((row & ~mask) != (col & ~mask))
                continue;
            size_t r = 0, c = 0;
            for (size_t i = 0; i < keep.size(); ++i)
            {
                r |= ((row >> keep[i]) & 1) << i;
                c |= ((col >> keep[i]) & 1) << i;
            }
            result(r, c) += rho(row, col);
        }
    return result / rho.trace();
}

template <class Sim> void CheckSingleQubitNormalization()
{
    for (int compression = 0; compression < 3; ++compression)
        for (bool scaled : {false, true})
        {
            Sim sim(4);
            sim.SetMultithreading(false);
            sim.ApplyGate(Gates::HadamardGate<>(), 1);
            sim.ApplyGate(Gates::CNOTGate<>(), 2, 1);
            sim.ApplyOperator(Gates::SingleQubitGate<>(0.5 * Matrix::Identity(2, 2)), 0);
            sim.ReCanonicalize();
            if constexpr (std::is_same_v<Sim, Simulator>)
                sim.MoveAtBeginningOfChain({3, 1});
            if (scaled)
            {
                auto state = sim.getState();
                auto base = std::dynamic_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(state);
                base->scaleExponent += 80;
                for (Index i = 0; i < base->gammas[0].size(); ++i)
                    base->gammas[0].data()[i] *= std::ldexp(1., -80);
                sim.setState(state);
            }
            if (compression == 1)
                sim.setLimitBondDimension(16);
            else if (compression == 2)
                sim.setLimitEntanglement(0.);

            Matrix expected = sim.getUnnormalizedDensityMatrix();
            const auto bonds = sim.getBondDimensions();
            // Zero, below-epsilon, and overflowing post-operation traces must
            // leave both the represented state and tensor storage unchanged.
            const auto before = std::dynamic_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
            for (double factor : {0., 1E-9, 1E200})
            {
                const Matrix bad = factor * Matrix::Identity(2, 2);
                Throws<std::runtime_error>([&] { sim.ApplyOperatorAndNormalize(Gates::SingleQubitGate<>(bad), 2); },
                                           "Invalid single-site normalization was accepted");
                Throws<std::runtime_error>([&] { sim.ApplyOperatorAndNormalize(Gates::AppliedGate<>(bad, 2)); },
                                           "Invalid applied single-site normalization was accepted");
                const auto after = std::dynamic_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
                Require(after->scaleExponent == before->scaleExponent && sim.getBondDimensions() == bonds,
                        "Failed single-site normalization changed scale or bonds");
                for (size_t q = 0; q < before->gammas.size(); ++q)
                    Require(std::equal(before->gammas[q].data(), before->gammas[q].data() + before->gammas[q].size(),
                                       after->gammas[q].data()),
                            "Failed single-site normalization changed a tensor");
                for (size_t q = 0; q < before->lambdas.size(); ++q)
                    Require((before->lambdas[q].array() == after->lambdas[q].array()).all(),
                            "Failed single-site normalization changed Schmidt weights");
                Close(sim.getUnnormalizedDensityMatrix(), expected, "Failed single-site normalization changed rho");
            }

            Matrix filter(2, 2);
            filter << .7, std::complex<double>(.1, .2), std::complex<double>(-.2, .1), .4;
            Matrix projectOne = Matrix::Zero(2, 2);
            projectOne(1, 1) = 1.;
            const std::vector<Matrix> ops{Gates::HadamardGate<>().getRawOperatorMatrix(), filter, filter, projectOne};
            const std::vector<Index> targets{3, 2, 0, 2};
            for (size_t i = 0; i < ops.size(); ++i)
            {
                const auto previousBonds = sim.getBondDimensions();
                sim.MeasureNoCollapse(); // populate the sampling cache before changing a site
                expected = DenseApply(expected, ops[i], targets[i]);
                expected /= expected.trace();
                if (i % 2 == 0)
                    sim.ApplyOperatorAndNormalize(Gates::SingleQubitGate<>(ops[i]), targets[i]);
                else
                    sim.ApplyOperatorAndNormalize(Gates::AppliedGate<>(ops[i], targets[i]));
                Close(sim.Trace(), 1., "Single-site normalization lost the trace scale");
                Close(sim.getUnnormalizedDensityMatrix(), expected, "Single-site normalization differs from dense rho");
                Require(sim.getBondDimensions() == previousBonds, "Single-site normalization changed bond dimensions");
                if (i + 1 == ops.size())
                    Require(sim.MeasureNoCollapse().at(2), "Single-site normalization retained stale sampling data");

                // Exercise canonical weights after the unitary case, and QR
                // transport after filters, before testing another local update.
                sim.ApplyGate(Gates::CNOTGate<>(), 2, 1);
                expected = DenseApply(expected, Gates::CNOTGate<>().getRawOperatorMatrix(), 2, 1);
                Close(sim.getUnnormalizedDensityMatrix(), expected, "Single-site normalization damaged later gates");
            }
        }
}

inline void SingleQubitNormalization()
{
    CheckSingleQubitNormalization<Implementation>();
    CheckSingleQubitNormalization<Simulator>();
}

inline void LocalEnvironmentsAndQueries()
{
    Simulator sim(5);
    sim.SetMultithreading(false);
    Matrix rho = Matrix::Zero(32, 32);
    rho(0, 0) = 1.;
    const Matrix x = Gates::PauliXGate<>().getRawOperatorMatrix();
    const std::vector<Matrix> noise{std::sqrt(.93) * Matrix::Identity(2, 2), std::sqrt(.07) * x};
    for (int step = 0; step < 30; ++step)
    {
        const Index q = (step * 7 + 2) % 5, target = (q + 2) % 5;
        const Gates::RyGate<> rotation(.13 + .011 * step);
        sim.ApplyGate(rotation, q);
        rho = DenseApply(rho, rotation.getRawOperatorMatrix(), q);
        sim.ApplyKrausOperators(noise, q);
        rho = (DenseApply(rho, noise[0], q) + DenseApply(rho, noise[1], q)).eval();
        sim.ApplyGate(Gates::CNOTGate<>(), target, q);
        rho = DenseApply(rho, Gates::CNOTGate<>().getRawOperatorMatrix(), target, q);
        Close(sim.getUnnormalizedDensityMatrix(), rho, "Local QR environments changed noisy routed evolution");
        if (step % 5 == 0)
        {
            Eigen::VectorXcd psi(32);
            for (Index i = 0; i < 32; ++i)
                psi[i] = std::complex<double>(std::sin(.7 * i), std::cos(.3 * i));
            psi.normalize();
            Close(sim.FidelityWithStatevector(psi), (psi.dot(rho * psi) / rho.trace()).real(),
                  "Tensor fidelity differs from dense reference");
            for (const auto &keep :
                 {std::vector<Index>{}, std::vector<Index>{4, 0, 2}, std::vector<Index>{3, 0, 4, 1, 2}})
                Close(sim.PartialTrace(keep), DensePartial(rho, keep),
                      "Reduced density matrix ordering or traced gaps are incorrect");
            Close(sim.HermiticityResidual(), (rho - rho.adjoint()).norm(),
                  "Tensor Hermiticity residual differs from dense reference");
        }
    }
    // Compare local truncation against a fully repaired copy of exactly the same input.
    sim.setLimitBondDimension(3);
    for (int step = 0; step < 12; ++step)
    {
        sim.ApplyKrausOperators(noise, step % 5);
        auto repaired = sim.Clone();
        repaired->ReCanonicalize();
        const Index q = step % 4;
        sim.ApplyGate(Gates::CNOTGate<>(), q + 1, q);
        repaired->ApplyGate(Gates::CNOTGate<>(), q + 1, q);
        Require((sim.getUnnormalizedDensityMatrix() - repaired->getUnnormalizedDensityMatrix()).norm() < 1E-8,
                "Local truncated update differs from fully canonical reference");
    }
    Simulator wide(14);
    for (int q = 0; q < 14; ++q)
        wide.ApplyGate(Gates::HadamardGate<>(), q);
    Close(wide.FidelityWithStatevector(Eigen::VectorXcd::Ones(16384) / 128.), 1.,
          "Wide statevector fidelity is incorrect");
    Eigen::VectorXcd sparse = Eigen::VectorXcd::Zero(16384);
    sparse[8195] = 1.;
    Close(wide.FidelityWithStatevector(sparse), 1. / 16384., "Sparse statevector fidelity is incorrect");
    Simulator complexSim(6);
    Matrix complexRho = Matrix::Zero(64, 64);
    complexRho(0, 0) = 1.;
    for (Index q = 0; q < 6; ++q)
    {
        Matrix op = Gates::RyGate<>(.21 + .12 * q).getRawOperatorMatrix();
        op.row(1) *= std::exp(std::complex<double>(0., .19 + .13 * q));
        complexSim.ApplyGate(Gates::SingleQubitGate<>(op), q);
        complexRho = DenseApply(complexRho, op, q);
    }
    complexSim.ApplyGate(Gates::CNOTGate<>(), 5, 0);
    complexRho = DenseApply(complexRho, Gates::CNOTGate<>().getRawOperatorMatrix(), 5, 0);
    complexSim.ApplyKrausOperators(noise, 3);
    complexRho = (DenseApply(complexRho, noise[0], 3) + DenseApply(complexRho, noise[1], 3)).eval();
    Eigen::VectorXcd densePsi(64);
    for (Index i = 0; i < 64; ++i)
        densePsi[i] = std::complex<double>(std::sin(.37 * i), std::cos(.23 * i));
    densePsi.normalize();
    Close(complexSim.FidelityWithStatevector(densePsi), densePsi.dot(complexRho * densePsi).real(),
          "Dense complex statevector contraction lost conjugation or logical bit order");
}

inline void StableHermiticity()
{
    Implementation wide(64);
    for (double epsilon : {0., 1E-12, 3E-10, .1})
    {
        auto state = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(wide.getState());
        state->gammas[0](0, 0, 1, 0) = std::complex<double>(0., epsilon);
        wide.setState(state);
        Require(std::abs(wide.HermiticityResidual() - std::sqrt(2.) * epsilon) < 1E-13,
                "Hermiticity residual lost small raw anti-Hermitian components");
        Require(wide.IsHermitian() == (epsilon < 1E-10), "Default Hermiticity tolerance changed");
    }
    auto state = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(wide.getState());
    state->gammas[0](0, 0, 0, 0) = 0.;
    wide.setState(state);
    Close(wide.HermiticityResidual(), std::sqrt(.02), "Zero-trace raw Hermiticity query was normalized");
}

inline void SamplingCacheAndSavedStates()
{
    Simulator sim(6), fresh(6);
    sim.ApplyGate(Gates::HadamardGate<>(), 0);
    sim.ApplyGate(Gates::CNOTGate<>(), 5, 0);
    for (int variant = 0; variant < 14; ++variant)
    {
        if (variant == 1)
            sim.ApplyGate(Gates::PauliXGate<>(), 3);
        if (variant == 2)
            sim.setToBasisState(41);
        if (variant == 3)
            sim.InitOnesState();
        if (variant == 4)
            sim.Clear();
        if (variant == 5)
            sim.setToMixtureOfBasisStates(std::vector<std::pair<size_t, double>>{{3, .4}, {60, .6}});
        if (variant == 6)
        {
            sim.ApplyGate(Gates::HadamardGate<>(), 1);
            sim.MeasureQubit(1);
        }
        if (variant == 7)
        {
            auto state = fresh.getState();
            sim.setStateDestructive(state);
        }
        if (variant == 8)
            sim.ApplyKrausOperators(std::vector<Matrix>{std::sqrt(.3) * Matrix::Identity(2, 2),
                                                        std::sqrt(.7) * Gates::PauliXGate<>().getRawOperatorMatrix()},
                                    2);
        if (variant == 9)
        {
            Matrix filter = Matrix::Identity(2, 2);
            filter(1, 1) = .7;
            sim.ApplyOperatorAndNormalize(Gates::SingleQubitGate<>(filter), 4);
        }
        if (variant == 10)
            sim.MoveAtBeginningOfChain({4});
        if (variant == 11)
        {
            sim.setLimitBondDimension(1);
            sim.Trim();
        }
        if (variant == 12)
            sim.Hermitize();
        if (variant == 13)
        {
            sim.ReCanonicalize();
            sim.RestoreTrace();
        }
        sim.SetSeed(711 + variant);
        fresh.SetSeed(711 + variant);
        for (int shot = 0; shot < 32; ++shot)
        {
            fresh.setState(sim.getState()); // force a new suffix calculation in the reference
            Require(sim.MeasureNoCollapse() == fresh.MeasureNoCollapse(), "Sampling reused stale tensor environments");
        }
    }
    sim.setToBasisState(37);
    sim.SaveState();
    auto clone = sim.Clone();
    sim.setToBasisState(2);
    sim.SaveState();
    clone->setToBasisState(7);
    clone->RestoreState();
    Close(clone->getBasisStateProbability(37), 1., "Clone lost or aliased its saved state");
    clone->ApplyGate(Gates::CNOTGate<>(), 4, 0);
    clone->RestoreStateDestructive();
    Close(clone->getBasisStateProbability(37), 1., "Destructive restore lost saved canonical metadata");
    // Private snapshots remain valid even when copied through the C++ copy constructor.
    Simulator copied = sim;
    copied.setToBasisState(13);
    copied.RestoreStateDestructive();
    sim.RestoreState();
    Close(sim.getBasisStateProbability(2), 1., "Destructive restore mutated another simulator's saved snapshot");
    sim.dontLimitBondDimension();
    sim.ApplyGate(Gates::HadamardGate<>(), 0);
    sim.ApplyKrausOperators(std::vector<Matrix>{std::sqrt(.9) * Matrix::Identity(2, 2),
                                                std::sqrt(.1) * Gates::PauliXGate<>().getRawOperatorMatrix()},
                            1);
    sim.ApplyGate(Gates::CNOTGate<>(), 4, 1);
    sim.SaveState();
    const Matrix mixedSaved = sim.getDensityMatrix();
    clone = sim.Clone();
    clone->Clear();
    clone->RestoreState();
    clone->ApplyGate(Gates::CNOTGate<>(), 5, 2);
    Close(clone->getDensityMatrix(), DenseApply(mixedSaved, Gates::CNOTGate<>().getRawOperatorMatrix(), 5, 2),
          "Saved mixed canonical metadata damaged subsequent gates");
}

inline void TransactionalNotifications()
{
    Simulator sim(5);
    sim.ApplyGate(Gates::HadamardGate<>(), 0);
    sim.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    const Matrix before = sim.getDensityMatrix();
    int notifications = 0;
    sim.SetBondDimensionCallback([&](const auto &) { ++notifications; });
    Throws<std::runtime_error>([&] { sim.ApplyOperatorAndNormalize(Gates::TwoQubitsGate<>(Matrix::Zero(4, 4)), 4, 0); },
                               "Zero routed operator was accepted");
    Require(notifications == 0, "Failed transaction emitted bond notifications");
    Close(sim.getDensityMatrix(), before, "Failed routed operation changed the state");
    sim.ApplyOperatorAndNormalize(Gates::CNOTGate<>(), 4, 0);
    Require(notifications == 1, "Successful routed transaction did not emit one committed notification");
    const auto expected = DenseApply(before, Gates::CNOTGate<>().getRawOperatorMatrix(), 4, 0);
    Close(sim.getDensityMatrix(), expected, "Routing rollback damaged a later normalized operation");
}

inline void PatchExceptionRecovery()
{
    class ThrowingRepair : public Implementation
    {
      public:
        ThrowingRepair() : Implementation(2)
        {
        }

        bool fail = true;

        void ReCanonicalize() override
        {
            if (fail)
            {
                fail = false;
                throw std::runtime_error("injected repair failure");
            }
            Implementation::ReCanonicalize();
        }
    };

    ThrowingRepair sim;
    Throws<std::runtime_error>([&] { sim.Hermitize(); }, "Exception injection failed");
    sim.Clear();
    sim.setRestoreTraceAfterTruncation(true);
    sim.setLimitBondDimension(1);
    sim.ApplyGate(Gates::RyGate<>(.2), 0);
    sim.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    Close(sim.Trace(), 1., "Exception permanently disabled post-truncation patches");
    ThrowingRepair nested;
    nested.setLimitBondDimension(1);
    nested.setHermitizeAfterTruncation(true);
    nested.ApplyGate(Gates::RyGate<>(.2), 0);
    Throws<std::runtime_error>([&] { nested.ApplyGate(Gates::CNOTGate<>(), 1, 0); },
                               "Nested patch exception injection failed");
    nested.Clear();
    nested.setHermitizeAfterTruncation(false);
    nested.setRestoreTraceAfterTruncation(true);
    nested.ApplyGate(Gates::RyGate<>(.2), 0);
    nested.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    Close(nested.Trace(), 1., "Exception left the outer post-truncation guard enabled");
}

inline void NormalizationAndNoOpTrim()
{
    Implementation sim(3);
    sim.setToBasisState(5);
    const Matrix before = sim.getDensityMatrix();
    Throws<std::invalid_argument>([&] { sim.setToQubitState(3); },
                                  "Out-of-range single-qubit initialization was accepted");
    Close(sim.getDensityMatrix(), before, "Rejected initialization cleared the state");
    auto state = sim.getState();
    sim.setState(state); // deliberately unknown canonical form
    sim.setLimitBondDimension(10);
    sim.Trim();
    const auto after = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
    const auto original = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(state);
    for (size_t q = 0; q < after->gammas.size(); ++q)
        Require(std::equal(after->gammas[q].data(), after->gammas[q].data() + after->gammas[q].size(),
                           original->gammas[q].data()),
                "No-op cap-only Trim changed tensor storage");
    sim.ApplyOperator(Gates::SingleQubitGate<>(Matrix::Zero(2, 2)), 0);
    Throws<std::runtime_error>([&] { sim.PartialTrace({0}); }, "Partial trace accepted zero trace");
    Throws<std::runtime_error>([&] { sim.HilbertSchmidtOverlap(sim); }, "Tensor overlap accepted zero trace");
    Throws<std::runtime_error>([&] { sim.FidelityWithStatevector(Eigen::VectorXcd::Ones(8)); },
                               "Fidelity accepted zero trace");
    Throws<std::runtime_error>([&] { sim.GetProbability(0); }, "Probability accepted zero trace");
    Throws<std::runtime_error>([&] { sim.ExpectationValue("III"); }, "Expectation accepted zero trace");
    Throws<std::runtime_error>([&] { sim.getBasisStateProbability(0); }, "Basis probability accepted zero trace");
    Close(sim.HermiticityResidual(), 0., "Raw zero operator residual failed");
}

template <class Sim> void CheckBasisProbabilityNormalization()
{
    Sim sim(1);
    auto state = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
    for (const auto trace : {std::complex<double>{}, std::complex<double>{-1., 0.}, std::complex<double>{0., 1.},
                             std::complex<double>{std::numeric_limits<double>::epsilon(), 0.}})
    {
        state->gammas[0].setZero();
        state->gammas[0](0, 0, 0, 0) = trace;
        sim.setState(state);
        Throws<std::runtime_error>([&] { sim.getBasisStateProbability(size_t{0}); },
                                   "Integer basis probability accepted an unsafe trace");
        Throws<std::runtime_error>([&] { sim.getBasisStateProbability(std::vector<bool>{false}); },
                                   "Vector basis probability accepted an unsafe trace");
    }
    // Finite tensor entries can still overflow when contracted into a trace.
    state->gammas[0](0, 0, 0, 0) = std::numeric_limits<double>::max();
    state->gammas[0](0, 1, 1, 0) = std::numeric_limits<double>::max();
    sim.setState(state);
    Require(!std::isfinite(sim.Trace().real()), "Non-finite trace test did not overflow");
    Throws<std::runtime_error>([&] { sim.getBasisStateProbability(size_t{0}); },
                               "Integer basis probability accepted a non-finite trace");
    Throws<std::runtime_error>([&] { sim.getBasisStateProbability(std::vector<bool>{false}); },
                               "Vector basis probability accepted a non-finite trace");
    state->gammas[0](0, 0, 0, 0) = 3.;
    state->gammas[0](0, 1, 1, 0) = 1.;
    sim.setState(state);
    for (size_t basis = 0; basis < 2; ++basis)
    {
        const double expected = basis == 0 ? .75 : .25;
        Close(sim.getBasisStateProbability(basis), expected,
              "Integer basis probability did not normalize by the trace");
        Close(sim.getBasisStateProbability(std::vector<bool>{basis != 0}), expected,
              "Vector basis probability did not normalize by the trace");
    }
}

inline void BasisProbabilityNormalization()
{
    CheckBasisProbabilityNormalization<Implementation>();
    CheckBasisProbabilityNormalization<Simulator>();
}

inline void FidelityMappings()
{
    Simulator sim(6);
    for (Index q = 0; q < 6; ++q)
    {
        sim.ApplyGate(Gates::RyGate<>(.31 + .17 * q), q);
        sim.ApplyGate(Gates::RxGate<>(.27 + .11 * q), q);
    }
    sim.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    sim.ApplyGate(Gates::CNOTGate<>(), 4, 3);
    sim.ApplyAmplitudeDamping(2, .23);
    sim.ApplyOperator(Gates::SingleQubitGate<>(.7 * Matrix::Identity(2, 2)), 0);
    const Matrix expectedRho = sim.getDensityMatrix();
    for (int mapping = 0; mapping < 3; ++mapping)
    {
        if (mapping == 1)
            sim.MoveAtBeginningOfChain({5}); // a non-self-inverse permutation
        if (mapping == 2)
            for (Index q = 0; q < 6; ++q)
                sim.MoveAtBeginningOfChain({q}); // reversal equals its inverse
        const auto originalMap = sim.getQubitsMap();
        for (Index count : {0, 1, 2, 32, 33, 64})
        {
            Eigen::VectorXcd psi = Eigen::VectorXcd::Zero(64);
            for (Index i = 0; i < count; ++i)
                psi[(13 * i + 7) % 64] = std::complex<double>(std::sin(.37 * (i + 1)), std::cos(.23 * (i + 1)));
            if (count)
                psi.normalize();
            Close(sim.FidelityWithStatevector(psi), psi.dot(expectedRho * psi).real(),
                  "Fidelity lost sparse amplitudes, conjugation, normalization or logical mapping");
        }
        Throws<std::invalid_argument>([&] { sim.FidelityWithStatevector(Eigen::VectorXcd::Ones(63)); },
                                      "Fidelity mapping bypassed statevector dimension validation");
        Require(sim.getQubitsMap() == originalMap, "Fidelity changed the logical mapping");
        Close(sim.getDensityMatrix(), expectedRho, "Fidelity changed the density matrix");
    }
    sim.ApplyOperator(Gates::SingleQubitGate<>(Matrix::Zero(2, 2)), 0);
    Throws<std::runtime_error>([&] { sim.FidelityWithStatevector(Eigen::VectorXcd::Zero(64)); },
                               "Mapped sparse fidelity accepted a zero-trace operator");
}

inline void ComplexSliceContractions()
{
    auto makeState = [](const std::vector<Index> &bonds, double phase) {
        auto state = std::make_shared<TensorNetworks::MPOSimulatorBaseState>();
        for (size_t q = 0; q + 1 < bonds.size(); ++q)
        {
            Interface::TensorType gamma(bonds[q], 2, 2, bonds[q + 1]);
            for (Index i = 0; i < gamma.size(); ++i)
                gamma.data()[i] =
                    .03 * std::complex<double>(std::sin((i + 1.) * (q + 1.) + phase), std::cos(.7 * (i + 1.) + phase));
            gamma(0, 0, 0, 0) += 1.;
            gamma(0, 1, 1, 0) += 1.;
            state->gammas.push_back(std::move(gamma));
            if (q + 2 < bonds.size())
                state->lambdas.push_back(Interface::LambdaType::Ones(bonds[q + 1]));
        }
        return state;
    };
    Implementation a(4), b(4);
    a.setState(makeState({1, 2, 3, 2, 1}, .2));
    b.setState(makeState({1, 3, 2, 4, 1}, .8));
    const Matrix rawA = a.getUnnormalizedDensityMatrix(), rawB = b.getUnnormalizedDensityMatrix();
    const auto expectedOverlap = (rawA.adjoint() * rawB).trace() / (std::conj(rawA.trace()) * rawB.trace());
    for (bool multithreading : {false, true})
    {
        a.SetMultithreading(multithreading);
        b.SetMultithreading(multithreading);
        Close(a.TraceOfSquare(), (rawA * rawA).trace(), "Strided trace-of-square contraction changed ket/bra ordering");
        Close(b.TraceOfSquare(), (rawB * rawB).trace(), "Strided trace-of-square contraction lost complex values");
        Close(a.HilbertSchmidtOverlap(b), expectedOverlap, "Strided overlap lost conjugation or mismatched bond sizes");
        Close(b.HilbertSchmidtOverlap(a), std::conj(expectedOverlap), "Strided overlap lost conjugate symmetry");
    }
}

inline bool Run()
{
    std::cout << "\nMPO regression tests (sampling seeds 0x4D504F01, 0x4D504F02)" << std::endl;
    const std::pair<const char *, void (*)()> tests[] = {
        {"stale canonical environments", StaleEnvironments},
        {"wide and scaled sampling", WideSampling},
        {"logical and large-register overlaps", LogicalOverlaps},
        {"atomic routed failures", AtomicRoutingFailures},
        {"invalid basis initialization", InvalidBasisState},
        {"dense-query bounds", DenseBounds},
        {"threshold-only trimming", ThresholdTrim},
        {"single-qubit normalization transactions", SingleQubitNormalization},
        {"local environments and tensor queries", LocalEnvironmentsAndQueries},
        {"stable raw Hermiticity", StableHermiticity},
        {"sampling cache and saved-state cloning", SamplingCacheAndSavedStates},
        {"transactional routing notifications", TransactionalNotifications},
        {"patch exception recovery", PatchExceptionRecovery},
        {"normalization and no-op trimming", NormalizationAndNoOpTrim},
        {"basis probability normalization", BasisProbabilityNormalization},
        {"sparse and dense fidelity mappings", FidelityMappings},
        {"complex strided slice contractions", ComplexSliceContractions}};
    bool passed = true;
    for (const auto &test : tests)
    {
        try
        {
            test.second();
            std::cout << "PASS: " << test.first << std::endl;
        }
        catch (const std::exception &error)
        {
            std::cout << "FAIL: " << test.first << ": " << error.what() << std::endl;
            passed = false;
        }
    }
    return passed;
}

} // namespace MPORegression
} // namespace QC
