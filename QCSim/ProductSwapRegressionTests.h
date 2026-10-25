#pragma once

#include <MPOSimulator.h>
#include <MPSSimulator.h>
#include <QuantumGate.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace QC
{
namespace ProductSwapRegression
{

using Matrix = Eigen::MatrixXcd;
using MPS = TensorNetworks::MPSSimulatorImpl;
using MPO = TensorNetworks::MPOSimulatorImpl;

inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline void Close(const Matrix &actual, const Matrix &expected, const char *message)
{
    Require(actual.rows() == expected.rows() && actual.cols() == expected.cols() && actual.allFinite(), message);
    const double scale = expected.cwiseAbs().maxCoeff();
    Require((actual - expected).cwiseAbs().maxCoeff() / (scale > 0. ? scale : 1.) < 5E-10, message);
}

inline Eigen::Index SwappedIndex(Eigen::Index i, int first, int second)
{
    return ((i >> first) & 1) == ((i >> second) & 1) ? i
                                                     : i ^ ((Eigen::Index{1} << first) | (Eigen::Index{1} << second));
}

inline Matrix DenseSwap(const Matrix &state, int first, int second)
{
    Matrix result(state.rows(), state.cols());
    for (Eigen::Index row = 0; row < state.rows(); ++row)
        for (Eigen::Index column = 0; column < state.cols(); ++column)
            result(row, column) =
                state(SwappedIndex(row, first, second), state.cols() == 1 ? 0 : SwappedIndex(column, first, second));
    return result;
}

template <class Sim> void Product(Sim &sim)
{
    sim.SetMultithreading(false);
    for (Eigen::Index q = 0; q < static_cast<Eigen::Index>(sim.getNrQubits()); ++q)
    {
        sim.ApplyGate(Gates::RyGate<>(.31 + .23 * q), q);
        sim.ApplyGate(Gates::RxGate<>(.17 + .19 * q), q);
    }
}

template <class Sim> void Compression(Sim &sim, int mode)
{
    if (mode == 0)
        return;
    sim.setLimitBondDimension(mode == 1 ? 1 : 2);
    sim.setTruncationMode(mode == 2 ? Sim::TruncationMode::RelativeToMax : Sim::TruncationMode::DiscardedWeight);
    sim.setLimitEntanglement(.03);
}

inline void CheckMPS()
{
    const Gates::SwapGate<> swap;
    // This deliberately untagged gate exercises the original contraction/SVD path.
    const Gates::TwoQubitsGate<> generalSwap(swap.getRawOperatorMatrix());
    Require(!generalSwap.isSwapGate(), "MPS reference must use the general gate path");
    for (int first : {0, 2, 4})
        for (bool reverse : {false, true})
            for (int mode : {0, 1, 2, 3})
                for (double weight : {.125, 1., 3.})
                {
                    MPS sim(6, 1);
                    Product(sim);
                    auto state = std::static_pointer_cast<TensorNetworks::MPSSimulatorBaseState>(sim.getState());
                    // Nonunit middle and left Schmidt values, independently scaled B
                    // tensors, and an unnormalized input. Lambdas are already in the Bs.
                    state->lambdas[first][0] = weight;
                    if (first != 0)
                        state->lambdas[first - 1][0] = weight;
                    for (Eigen::Index p = 0; p < 2; ++p)
                    {
                        state->gammas[first].data()[p] *= weight;
                        state->gammas[first + 1].data()[p] *= 1.7;
                    }
                    sim.setState(state);
                    sim.EnableBondDimensionSummary(true);
                    Compression(sim, mode);
                    MPS reference = sim;
                    const int target = first + (reverse ? 0 : 1), control = first + (reverse ? 1 : 0);
                    Matrix dense = DenseSwap(sim.getRegisterStorage(), first, first + 1);
                    const double expectedNorm = weight * 1.7 * (first == 0 ? 1. : weight);
                    dense /= expectedNorm;
                    sim.ApplyGate(swap, target, control);
                    reference.ApplyGate(generalSwap, target, control);
                    Close(sim.getRegisterStorage(), dense, "MPS product SWAP lost normalization or a scalar weight");
                    Close(sim.getRegisterStorage(), reference.getRegisterStorage(),
                          "MPS product SWAP differs from SVD");
                    state = std::static_pointer_cast<TensorNetworks::MPSSimulatorBaseState>(sim.getState());
                    Require(state->lambdas[first].size() == 1 && std::abs(state->lambdas[first][0] - 1.) < 1E-12 &&
                                sim.getMaxBondDimension() == 1,
                            "MPS product SWAP left invalid Schmidt weights or dimension summary");
                    for (int q : {4, 0, 3, 1, 2})
                    {
                        sim.ApplyGate(Gates::CNOTGate<>(), q + 1, q);
                        reference.ApplyGate(Gates::CNOTGate<>(), q + 1, q);
                        Close(sim.getRegisterStorage(), reference.getRegisterStorage(),
                              "MPS product SWAP damaged subsequent compressed gates");
                        Require(sim.getBondDimensions() == reference.getBondDimensions(), "MPS ranks changed");
                    }
                }

    // A rank-one middle bond does not make sites with entangled outer bonds independent.
    MPS entangled(4, 1);
    entangled.ApplyGate(Gates::HadamardGate<>(), 0);
    entangled.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    entangled.ApplyGate(Gates::HadamardGate<>(), 2);
    entangled.ApplyGate(Gates::CNOTGate<>(), 3, 2);
    Require(entangled.getBondDimensions() == std::vector<Eigen::Index>({2, 1, 2}), "Invalid MPS fallback fixture");
    const Matrix expected = DenseSwap(entangled.getRegisterStorage(), 1, 2);
    entangled.ApplyGate(swap, 1, 2);
    Close(entangled.getRegisterStorage(), expected, "MPS skipped an entangled outer bond");
    Require(entangled.getBondDimension(1) == 4, "MPS SWAP failed to grow an entangled cut");

    MPS repeated(6, 1);
    Product(repeated);
    const Matrix original = repeated.getRegisterStorage();
    for (int iteration = 0; iteration < 100; ++iteration)
        for (int first = 0; first < 5; ++first)
        {
            repeated.ApplyGate(swap, first, first + 1);
            repeated.ApplyGate(swap, first + 1, first);
        }
    Close(repeated.getRegisterStorage(), original, "Repeated MPS product SWAPs drifted");
    repeated.setToBasisState(size_t{1});
    repeated.ApplyGate(swap, 0, 1);
    const auto sample = repeated.MeasureNoCollapse();
    Require(!sample.at(0) && sample.at(1), "MPS product SWAP damaged sampling");
}

class MPOProbe : public MPO
{
  public:
    using MPO::MPO;

    CanonicalMetadata Metadata() const
    {
        return GetCanonicalMetadata();
    }
};

inline void CheckMPO()
{
    const Gates::SwapGate<> swap;
    const std::vector<Matrix> swapChannel{swap.getRawOperatorMatrix()};
    for (int first : {0, 2, 4})
        for (bool reverse : {false, true})
            for (bool canonical : {false, true})
                for (int mode : {0, 1, 2, 3})
                    for (double weight : {.125, 1., 3.})
                    {
                        MPOProbe sim(6, 1);
                        Product(sim);
                        for (int q = 0; q < 6; ++q)
                            sim.ApplyDepolarizingNoise(q, .11 + .03 * q);
                        sim.ApplyOperator(Gates::SingleQubitGate<>(std::sqrt(weight) * Matrix::Identity(2, 2)), 0);
                        if (canonical)
                            sim.ReCanonicalize();
                        auto state = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
                        Require(state->lambdas[first].size() == 1, "Mixed MPO fixture is not a product");
                        if (canonical)
                            Require(std::abs(state->lambdas[first][0] - 1.) > .01,
                                    "MPO fixture must exercise a nonunit singular value");
                        sim.EnableBondDimensionSummary(true);
                        Compression(sim, mode);
                        sim.setHermitizeAfterTruncation(true);
                        sim.setRestoreTraceAfterTruncation(true);
                        MPO reference = sim;
                        const Matrix before = sim.getUnnormalizedDensityMatrix();
                        const int target = first + (reverse ? 0 : 1), control = first + (reverse ? 1 : 0);
                        // Fill the sampling cache before changing the representation.
                        sim.MeasureNoCollapse();
                        sim.ApplyGate(swap, target, control);
                        reference.ApplyKrausOperators(swapChannel, target, control);
                        Close(sim.getUnnormalizedDensityMatrix(), DenseSwap(before, first, first + 1),
                              "MPO product SWAP lost trace, scale or ket/bra ordering");
                        Close(sim.getUnnormalizedDensityMatrix(), reference.getUnnormalizedDensityMatrix(),
                              "MPO product SWAP differs from the general channel");
                        Require(sim.getMaxBondDimension() == 1, "MPO product SWAP damaged the dimension summary");
                        const auto metadata = sim.Metadata();
                        Require(metadata.valid == canonical &&
                                    (canonical || (metadata.first == first && metadata.last == first)),
                                "MPO product SWAP lost canonical-center metadata");
                        state = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
                        Require(std::abs(Eigen::Map<Eigen::Vector4cd>(state->gammas[first + 1].data()).norm() - 1.) <
                                    1E-12,
                                "MPO product SWAP left a noncanonical right tensor");
                        sim.ApplyGate(swap, control, target);
                        reference.ApplyKrausOperators(swapChannel, control, target);
                        Close(sim.getUnnormalizedDensityMatrix(), before, "Two MPO product SWAPs changed the operator");
                        for (int q : {4, 0, 3, 1, 2})
                        {
                            sim.ApplyGate(Gates::CNOTGate<>(), q + 1, q);
                            reference.ApplyGate(Gates::CNOTGate<>(), q + 1, q);
                            Close(sim.getUnnormalizedDensityMatrix(), reference.getUnnormalizedDensityMatrix(),
                                  "MPO product SWAP damaged subsequent compressed gates");
                            Require(sim.getBondDimensions() == reference.getBondDimensions(), "MPO ranks changed");
                        }
                    }

    // Arbitrary imported gauges, complex operators, and explicit scale exponents.
    // Relative comparisons also catch loss of very small (nonzero) operators.
    for (int exponent : {-800, -400, 0, 400, 800})
        for (int first : {0, 2})
        {
            MPO sim(4, 1);
            Product(sim);
            auto state = std::static_pointer_cast<TensorNetworks::MPOSimulatorBaseState>(sim.getState());
            for (int p = 0; p < 4; ++p)
            {
                state->gammas[first].data()[p] *= std::complex<double>(.3, .4) * std::ldexp(1., exponent);
                state->gammas[first + 1].data()[p] *= 2.;
            }
            for (auto &lambda : state->lambdas)
                lambda[0] = .37; // deliberately stale: setState must prepare its environments
            state->scaleExponent = 23;
            sim.setState(state);
            MPO reference = sim;
            const Matrix expected = DenseSwap(sim.getUnnormalizedDensityMatrix(), first, first + 1);
            sim.ApplyGate(swap, first, first + 1);
            reference.ApplyKrausOperators(swapChannel, first, first + 1);
            Close(sim.getUnnormalizedDensityMatrix(), expected, "MPO product SWAP lost an imported scale");
            Close(sim.getUnnormalizedDensityMatrix(), reference.getUnnormalizedDensityMatrix(),
                  "Scaled MPO product SWAP differs from SVD");
            sim.ApplyGate(Gates::CNOTGate<>(), 2, 1);
            reference.ApplyGate(Gates::CNOTGate<>(), 2, 1);
            Close(sim.getUnnormalizedDensityMatrix(), reference.getUnnormalizedDensityMatrix(),
                  "Scaled MPO product SWAP damaged a later gate");
        }

    MPO entangled(4, 1);
    entangled.ApplyGate(Gates::HadamardGate<>(), 0);
    entangled.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    entangled.ApplyGate(Gates::HadamardGate<>(), 2);
    entangled.ApplyGate(Gates::CNOTGate<>(), 3, 2);
    Require(entangled.getBondDimensions() == std::vector<Eigen::Index>({4, 1, 4}), "Invalid MPO fallback fixture");
    const Matrix swappedEntangled = DenseSwap(entangled.getUnnormalizedDensityMatrix(), 1, 2);
    entangled.ApplyGate(swap, 1, 2);
    Close(entangled.getUnnormalizedDensityMatrix(), swappedEntangled, "MPO skipped an entangled outer bond");
    Require(entangled.getBondDimension(1) == 16, "MPO SWAP failed to grow an entangled cut");

    MPO repeated(6, 1);
    Product(repeated);
    for (int q = 0; q < 6; ++q)
        repeated.ApplyDepolarizingNoise(q, .17);
    repeated.ReCanonicalize();
    const Matrix original = repeated.getUnnormalizedDensityMatrix();
    for (int iteration = 0; iteration < 100; ++iteration)
        for (int first = 0; first < 5; ++first)
        {
            repeated.ApplyGate(swap, first, first + 1);
            repeated.ApplyGate(swap, first + 1, first);
        }
    Close(repeated.getUnnormalizedDensityMatrix(), original, "Repeated mixed MPO product SWAPs drifted");

    MPO cached(3, 1);
    cached.setToBasisState(size_t{1});
    cached.MeasureNoCollapse();
    cached.ApplyGate(swap, 0, 1);
    const auto sample = cached.MeasureNoCollapse();
    Require(!sample.at(0) && sample.at(1), "MPO product SWAP reused a stale sampling cache");
    cached.ApplyOperator(Gates::SingleQubitGate<>(Matrix::Zero(2, 2)), 0);
    cached.ApplyGate(swap, 0, 1);
    Require(cached.getUnnormalizedDensityMatrix().norm() == 0., "MPO zero-operator fallback failed");

    // Physical reordering of copies for mapped overlap shares the shortcut.
    TensorNetworks::MPOSimulator a(6, 1), b(6, 1);
    Product(a);
    Product(b);
    for (int q = 0; q < 6; ++q)
    {
        a.ApplyDepolarizingNoise(q, .07 * q);
        b.ApplyAmplitudeDamping(q, .03 * q);
    }
    a.ApplyGate(swap, 0, 5);
    b.ApplyGate(swap, 1, 4);
    const Matrix rhoA = a.getUnnormalizedDensityMatrix(), rhoB = b.getUnnormalizedDensityMatrix();
    const auto overlap = (rhoA.adjoint() * rhoB).trace();
    Require(std::abs(a.HilbertSchmidtOverlap(b) - overlap) < 1E-11, "MPO product overlap reordering failed");
    Close(a.getUnnormalizedDensityMatrix(), rhoA, "Overlap changed the caller");
    Close(b.getUnnormalizedDensityMatrix(), rhoB, "Overlap changed the argument");
}

inline bool Run(bool mpo)
{
    try
    {
        if (mpo)
            CheckMPO();
        else
            CheckMPS();
        std::cout << (mpo ? "MPO" : "MPS") << " product SWAP regressions passed" << std::endl;
        return true;
    }
    catch (const std::exception &error)
    {
        std::cout << (mpo ? "MPO" : "MPS") << " product SWAP regression: " << error.what() << std::endl;
        return false;
    }
}

} // namespace ProductSwapRegression
} // namespace QC
