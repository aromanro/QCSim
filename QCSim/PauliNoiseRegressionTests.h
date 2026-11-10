#pragma once
#include "PauliNativeRegressionTests.h"
#include "PauliPropagator.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

// Noise channels of the Pauli propagator, checked against a dense density
// matrix evolved in the Schrodinger picture. The propagator applies the
// adjoint channel to the observable, so one recorded circuit must reproduce
// the ensemble averages of the mixed state exactly.

namespace QC
{
namespace PauliNoiseRegression
{
using PauliNativeRegression::Embed;
using PauliNativeRegression::Matrix;
using PauliNativeRegression::Pauli;
using PauliNativeRegression::Require;
constexpr double pi = PauliNativeRegression::pi;

inline void Close(double a, double b, const char *message, double eps = 1e-10)
{
    Require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) < eps, message);
}

// Public channel arguments index a qubit as I, X, Y, Z. Local labels (and the
// Pauli() helper) use I=0, X=1, Z=2, Y=3. Kept separate from the code under test.
inline unsigned LabelOfIndex(unsigned index, unsigned qubits)
{
    const unsigned code[4] = {0, 1, 3, 2};
    unsigned label = 0;
    for (unsigned q = 0; q < qubits; ++q)
        label |= code[(index >> (2 * q)) & 3] << (2 * q);
    return label;
}

inline std::string PauliString(unsigned label, int width)
{
    std::string s(static_cast<size_t>(width), 'I');
    for (int q = 0; q < width; ++q)
        s[q] = "IXZY"[(label >> (2 * q)) & 3];
    return s;
}

struct Dense
{
    int width;
    Matrix rho;

    explicit Dense(int n) : width(n), rho(Matrix::Zero(Eigen::Index(1) << n, Eigen::Index(1) << n))
    {
        rho(0, 0) = 1.;
    }

    void Unitary(const Matrix &u, const std::array<int, 3> &q, int arity)
    {
        const Matrix full = Embed(u, q, arity, width);
        rho = (full * rho * full.adjoint()).eval();
    }

    void Channel(const std::vector<Matrix> &kraus, const std::array<int, 3> &q, int arity)
    {
        Matrix next = Matrix::Zero(rho.rows(), rho.cols());
        for (const auto &k : kraus)
        {
            const Matrix full = Embed(k, q, arity, width);
            next += full * rho * full.adjoint();
        }
        rho = next;
    }

    double Expectation(unsigned label) const
    {
        return (rho * Pauli(label, static_cast<unsigned>(width))).trace().real();
    }

    double Probability(size_t outcome) const
    {
        return rho(Eigen::Index(outcome), Eigen::Index(outcome)).real();
    }
};

// ---- reference Kraus operators ----

inline std::vector<Matrix> PauliKraus(const std::vector<double> &probabilities, unsigned qubits)
{
    std::vector<Matrix> kraus;
    for (unsigned index = 0; index < probabilities.size(); ++index)
        if (probabilities[index] > 0.)
            kraus.push_back(std::sqrt(probabilities[index]) * Pauli(LabelOfIndex(index, qubits), qubits));
    return kraus;
}

inline std::vector<Matrix> PhaseDampingKraus(double gamma)
{
    Matrix k0 = Matrix::Zero(2, 2), k1 = Matrix::Zero(2, 2);
    k0(0, 0) = 1.;
    k0(1, 1) = std::sqrt(1. - gamma);
    k1(1, 1) = std::sqrt(gamma);
    return {k0, k1};
}

inline std::vector<Matrix> DampingKraus(double gamma, double excited)
{
    const double ground = 1. - excited;
    Matrix k0 = Matrix::Zero(2, 2), k1 = Matrix::Zero(2, 2), k2 = Matrix::Zero(2, 2), k3 = Matrix::Zero(2, 2);
    k0(0, 0) = std::sqrt(ground);
    k0(1, 1) = std::sqrt(ground * (1. - gamma));
    k1(0, 1) = std::sqrt(ground * gamma);
    k2(0, 0) = std::sqrt(excited * (1. - gamma));
    k2(1, 1) = std::sqrt(excited);
    k3(1, 0) = std::sqrt(excited * gamma);
    return {k0, k1, k2, k3};
}

inline std::vector<Matrix> ResetKraus(bool toOne)
{
    Matrix k0 = Matrix::Zero(2, 2), k1 = Matrix::Zero(2, 2);
    k0(toOne ? 1 : 0, 0) = 1.;
    k1(toOne ? 1 : 0, 1) = 1.;
    return {k0, k1};
}

// A random trace-preserving channel: E_k = A_k (sum_j A_j^dagger A_j)^(-1/2).
inline std::vector<Matrix> RandomKraus(std::mt19937_64 &rng, unsigned qubits, int count)
{
    std::normal_distribution<double> normal;
    const Eigen::Index dim = Eigen::Index(1) << qubits;
    std::vector<Matrix> kraus(static_cast<size_t>(count), Matrix::Zero(dim, dim));
    Matrix sum = Matrix::Zero(dim, dim);
    for (auto &k : kraus)
    {
        for (Eigen::Index r = 0; r < dim; ++r)
            for (Eigen::Index c = 0; c < dim; ++c)
                k(r, c) = std::complex<double>(normal(rng), normal(rng));
        sum += k.adjoint() * k;
    }
    const Matrix inverseRoot = Eigen::SelfAdjointEigenSolver<Matrix>(sum).operatorInverseSqrt();
    for (auto &k : kraus)
        k = (k * inverseRoot).eval();
    return kraus;
}

// The observable sum_k E_k^dagger P E_k on the channel qubits.
inline Matrix Adjoint(const std::vector<Matrix> &kraus, const Matrix &p)
{
    Matrix result = Matrix::Zero(p.rows(), p.cols());
    for (const auto &k : kraus)
        result += k.adjoint() * p * k;
    return result;
}

inline size_t comparisons = 0;

// Every Pauli on the channel qubits, propagated through the one recorded
// operation of sim, against the dense adjoint. q lists the channel qubits.
// This covers outputs with zero expectation on |0>, and spectator qubits.
inline void CheckAction(PauliPropagator &sim, const std::vector<Matrix> &kraus, const std::array<int, 3> &q,
                        unsigned qubits, const char *message)
{
    const int width = sim.GetNrQubits();
    for (unsigned input = 0; input < (1U << (2 * qubits)); ++input)
    {
        PauliStringXZWithCoefficient p(static_cast<size_t>(width));
        for (unsigned j = 0; j < qubits; ++j)
        {
            p.X[q[j]] = (input >> (2 * j)) & 1;
            p.Z[q[j]] = (input >> (2 * j + 1)) & 1;
        }
        // A spectator factor must pass through untouched. The callers leave
        // the last qubit of the register free.
        const int spectator = width - 1;
        const bool useSpectator = spectator != q[0] && (qubits == 1 || spectator != q[1]);
        if (useSpectator)
            p.X[spectator] = p.Z[spectator] = true;
        PauliStringStorage output;
        sim.ExpectationValue(std::move(p), output);
        Matrix got = Matrix::Zero(Eigen::Index(1) << qubits, Eigen::Index(1) << qubits);
        for (const auto &term : output)
        {
            if (term.Coefficient == 0.)
                continue;
            unsigned label = 0;
            for (unsigned j = 0; j < qubits; ++j)
                label |= unsigned(term.X[q[j]]) << (2 * j) | unsigned(term.Z[q[j]]) << (2 * j + 1);
            for (int other = 0; other < width; ++other)
            {
                if (other == q[0] || (qubits == 2 && other == q[1]))
                    continue;
                const bool expected = useSpectator && other == spectator;
                Require(bool(term.X[other]) == expected && bool(term.Z[other]) == expected,
                        "channel touched a spectator");
            }
            got += term.Coefficient * Pauli(label, qubits);
        }
        const Matrix want = Adjoint(kraus, Pauli(input, qubits));
        Require((got - want).norm() < 1e-11, message);
        ++comparisons;
    }
}

inline void Actions()
{
    std::mt19937_64 rng(20261110);
    std::uniform_real_distribution<double> uniform(0., 1.);
    // Qubit placements: one word, across a word boundary, reversed order.
    const std::array<std::array<int, 3>, 4> placements{{{{0, 1, 0}}, {{2, 0, 0}}, {{0, 65, 0}}, {{129, 3, 0}}}};
    for (const auto &q : placements)
    {
        const int width = std::max(q[0], q[1]) + 2;
        const auto fresh = [&] {
            PauliPropagator sim;
            sim.SetNrQubits(width);
            return sim;
        };
        for (int repeat = 0; repeat < 4; ++repeat)
        {
            // One-qubit Pauli channel, asymmetric in X, Y and Z.
            {
                const double px = 0.3 * uniform(rng), py = 0.3 * uniform(rng), pz = 0.3 * uniform(rng);
                auto sim = fresh();
                sim.ApplyPauliChannel(q[0], px, py, pz);
                CheckAction(sim, PauliKraus({1. - px - py - pz, px, py, pz}, 1), q, 1, "one-qubit Pauli channel");
            }
            // Two-qubit Pauli channel, not symmetric under exchange of its qubits.
            {
                std::array<double, 16> probabilities;
                double total = 0.;
                for (auto &p : probabilities)
                    total += p = uniform(rng) * uniform(rng);
                for (auto &p : probabilities)
                    p /= total;
                auto sim = fresh();
                sim.ApplyTwoQubitPauliChannel(q[0], q[1], probabilities);
                CheckAction(sim, PauliKraus(std::vector<double>(probabilities.begin(), probabilities.end()), 2), q, 2,
                            "two-qubit Pauli channel");
            }
            {
                const double gamma = uniform(rng);
                auto sim = fresh();
                sim.ApplyPhaseDamping(q[0], gamma);
                CheckAction(sim, PhaseDampingKraus(gamma), q, 1, "phase damping");
            }
            {
                const double gamma = uniform(rng), excited = repeat ? uniform(rng) : 0.;
                auto sim = fresh();
                sim.ApplyAmplitudeDamping(q[0], gamma, excited);
                CheckAction(sim, DampingKraus(gamma, excited), q, 1, "amplitude damping");
            }
            // Thermal relaxation, including t1 < t2 <= 2 t1. The reference
            // composes damping with the extra dephasing that reaches t2.
            {
                const double t1 = 50. + 100. * uniform(rng), t2 = t1 * (0.2 + 1.8 * uniform(rng)),
                             duration = 80. * uniform(rng);
                const double excited = repeat % 2 ? 0.3 * uniform(rng) : 0.;
                auto sim = fresh();
                sim.ApplyThermalRelaxation(q[0], duration, t1, t2, excited);
                const double gamma = -std::expm1(-duration / t1);
                const double extra = std::exp(-duration / t2) / std::sqrt(1. - gamma);
                auto kraus = DampingKraus(gamma, excited);
                std::vector<Matrix> composed;
                for (const auto &flip : PauliKraus({0.5 * (1. + extra), 0., 0., 0.5 * (1. - extra)}, 1))
                    for (const auto &k : kraus)
                        composed.push_back(flip * k);
                CheckAction(sim, composed, q, 1, "thermal relaxation");
            }
            for (unsigned qubits : {1U, 2U})
            {
                const auto kraus = RandomKraus(rng, qubits, 1 + repeat);
                auto sim = fresh();
                sim.ApplyKrausChannel(kraus, q[0], q[1]);
                CheckAction(sim, kraus, q, qubits, "Kraus channel");
            }
            for (bool toOne : {false, true})
            {
                auto sim = fresh();
                sim.ApplyReset(q[0], toOne);
                CheckAction(sim, ResetKraus(toOne), q, 1, "reset");
            }
        }
        // Boundary strengths: identity channels and complete depolarization.
        for (double p : {0., 0.75, 1.})
        {
            auto sim = fresh();
            sim.ApplyDepolarizingNoise(q[0], p);
            CheckAction(sim, PauliKraus({1. - p, p / 3., p / 3., p / 3.}, 1), q, 1, "depolarizing boundary");
            auto two = fresh();
            two.ApplyTwoQubitDepolarizingNoise(q[0], q[1], p);
            std::vector<double> probabilities(16, p / 15.);
            probabilities[0] = 1. - p;
            CheckAction(two, PauliKraus(probabilities, 2), q, 2, "two-qubit depolarizing boundary");
        }
        for (double p : {0., 1.})
        {
            auto sim = fresh();
            sim.ApplyBitFlipNoise(q[0], p);
            CheckAction(sim, PauliKraus({1. - p, p, 0., 0.}, 1), q, 1, "bit flip boundary");
            auto phase = fresh();
            phase.ApplyPhaseFlipNoise(q[0], p);
            CheckAction(phase, PauliKraus({1. - p, 0., 0., p}, 1), q, 1, "phase flip boundary");
            auto damping = fresh();
            damping.ApplyAmplitudeDamping(q[0], p);
            CheckAction(damping, DampingKraus(p, 0.), q, 1, "amplitude damping boundary");
        }
    }
    std::cout << "PASS " << comparisons << " dense adjoint-channel comparisons (1Q/2Q, word boundaries, spectators)\n";
}

// How a channel is recorded: a diagonal transfer matrix must not branch,
// however the channel was specified.
inline void Recording()
{
    using T = OperationType;
    const auto type = [](const PauliPropagator &sim) { return sim.GetOperations().back()->GetType(); };
    PauliPropagator sim;
    sim.SetNrQubits(3);
    sim.ApplyDepolarizingNoise(0, 0.1);
    Require(type(sim) == T::PAULI_CHANNEL_1Q, "depolarizing is a diagonal channel");
    sim.ApplyPhaseDamping(1, 0.3);
    Require(type(sim) == T::PAULI_CHANNEL_1Q, "phase damping is a diagonal channel");
    sim.ApplyKrausChannel(PhaseDampingKraus(0.3), 1);
    Require(type(sim) == T::PAULI_CHANNEL_1Q, "Kraus phase damping was not recognized as diagonal");
    sim.ApplyKrausChannel(PauliKraus({0.7, 0.1, 0.05, 0.15}, 1), 2);
    Require(type(sim) == T::PAULI_CHANNEL_1Q, "Kraus Pauli channel was not recognized as diagonal");
    std::vector<double> two(16, 0.02);
    two[0] = 0.7;
    sim.ApplyKrausChannel(PauliKraus(two, 2), 2, 0);
    Require(type(sim) == T::PAULI_CHANNEL_2Q, "two-qubit Kraus Pauli channel was not recognized as diagonal");
    sim.ApplyTwoQubitDepolarizingNoise(0, 2, 0.2);
    Require(type(sim) == T::PAULI_CHANNEL_2Q, "two-qubit depolarizing is a diagonal channel");
    sim.ApplyAmplitudeDamping(0, 0.2);
    Require(type(sim) == T::CHANNEL_1Q, "amplitude damping is a local table");
    // At the symmetric population the shift to the identity vanishes.
    sim.ApplyAmplitudeDamping(0, 0.2, 0.5);
    Require(type(sim) == T::PAULI_CHANNEL_1Q, "symmetric damping is diagonal");
    sim.ApplyThermalRelaxation(1, 10., 100., 150.);
    Require(type(sim) == T::CHANNEL_1Q, "thermal relaxation is a local table");
    sim.ApplyKrausChannel(DampingKraus(0.2, 0.), 1);
    Require(type(sim) == T::CHANNEL_1Q, "Kraus damping is a local table");
    sim.ApplyReset(2);
    Require(type(sim) == T::RESET, "reset type");
    Require(sim.GetOperations().size() == 11, "one recorded operation per channel");

    // Every table maps the identity to exactly itself: a channel away from the
    // observable must leave both its value and its cone alone.
    for (const auto &op : sim.GetOperations())
        if (IsPauliLocalGate(op->GetType()))
            Require(static_cast<const OperatorLocal &>(*op).GetTransfer()->unchanged & 1, "identity column not pinned");

    // A unitary given as one Kraus operator keeps its Clifford structure.
    PauliPropagator clifford;
    clifford.SetNrQubits(1);
    Matrix h(2, 2);
    h << 1., 1., 1., -1.;
    clifford.ApplyKrausChannel(std::vector<Matrix>{h / std::sqrt(2.)}, 0);
    Require(static_cast<const OperatorLocal &>(*clifford.GetOperations().back()).GetTransfer()->clifford,
            "Hadamard Kraus operator lost its Clifford table");
    Close(clifford.ExpectationValue("X"), 1., "Hadamard channel");

    // Validation.
    const auto rejects = [](auto &&call) {
        try
        {
            call();
        }
        catch (const std::invalid_argument &)
        {
            return true;
        }
        catch (const std::out_of_range &)
        {
            return true;
        }
        return false;
    };
    PauliPropagator bad;
    bad.SetNrQubits(2);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Require(rejects([&] { bad.ApplyPauliChannel(0, 0.5, 0.4, 0.2); }), "probabilities above one accepted");
    Require(rejects([&] { bad.ApplyPauliChannel(0, -0.1, 0., 0.); }), "negative probability accepted");
    Require(rejects([&] { bad.ApplyPauliChannel(0, nan, 0., 0.); }), "NaN probability accepted");
    Require(rejects([&] { bad.ApplyDepolarizingNoise(2, 0.1); }), "channel outside the register accepted");
    Require(rejects([&] { bad.ApplyTwoQubitDepolarizingNoise(1, 1, 0.1); }), "repeated channel qubit accepted");
    Require(rejects([&] { bad.ApplyTwoQubitPauliChannel(0, 1, std::array<double, 16>{}); }),
            "unnormalized channel accepted");
    Require(rejects([&] { bad.ApplyPhaseDamping(0, 1.5); }), "phase damping above one accepted");
    Require(rejects([&] { bad.ApplyAmplitudeDamping(0, 0.5, -0.1); }), "negative population accepted");
    Require(rejects([&] { bad.ApplyThermalRelaxation(0, 1., 10., 21.); }), "T2 > 2 T1 accepted");
    Require(rejects([&] { bad.ApplyThermalRelaxation(0, -1., 10., 10.); }), "negative duration accepted");
    Require(rejects([&] { bad.ApplyThermalRelaxation(0, 1., 0., 10.); }), "zero T1 accepted");
    Require(rejects([&] { bad.ApplyReset(-1); }), "reset outside the register accepted");
    Require(rejects([&] { bad.ApplyKrausChannel(std::vector<Matrix>{}, 0); }), "empty channel accepted");
    Require(rejects([&] { bad.ApplyKrausChannel(std::vector<Matrix>{Matrix::Identity(3, 3)}, 0); }),
            "3x3 Kraus accepted");
    Require(rejects([&] { bad.ApplyKrausChannel(std::vector<Matrix>{0.9 * Matrix::Identity(2, 2)}, 0); }),
            "trace-decreasing channel accepted");
    Require(
        rejects([&] { bad.ApplyKrausChannel(std::vector<Matrix>{Matrix::Identity(2, 2), Matrix::Identity(4, 4)}, 0); }),
        "mixed Kraus sizes accepted");
    Matrix infinite = Matrix::Identity(2, 2);
    infinite(0, 1) = std::numeric_limits<double>::infinity();
    Require(rejects([&] { bad.ApplyKrausChannel(std::vector<Matrix>{infinite}, 0); }), "nonfinite Kraus accepted");
    Require(bad.GetOperations().empty(), "a rejected channel was recorded");
    bad.ApplyThermalRelaxation(0, 1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    Close(bad.ExpectationValue("ZI"), 1., "infinite relaxation times are the identity");
    std::cout << "PASS channel recording, diagonal detection, pinned identity and validation\n";
}

struct Step
{
    int kind;
    std::array<int, 3> q;
    double a, b, c;
    std::vector<Matrix> kraus;
    std::array<double, 16> probabilities;
};

inline std::vector<Step> RandomCircuit(std::mt19937_64 &rng, int width, int length)
{
    std::uniform_real_distribution<double> uniform(0., 1.);
    std::uniform_int_distribution<int> qubit(0, width - 1), kind(0, 14);
    std::vector<Step> steps;
    for (int i = 0; i < length; ++i)
    {
        Step s{};
        s.kind = kind(rng);
        s.q[0] = qubit(rng);
        do
            s.q[1] = qubit(rng);
        while (s.q[1] == s.q[0]);
        s.a = uniform(rng);
        s.b = uniform(rng);
        s.c = uniform(rng);
        if (s.kind == 11)
            s.kraus = RandomKraus(rng, 1, 3);
        else if (s.kind == 12)
            s.kraus = RandomKraus(rng, 2, 2);
        else if (s.kind == 8)
        {
            double total = 0.;
            for (auto &p : s.probabilities)
                total += p = uniform(rng) * uniform(rng);
            // Mostly the identity, as a noise channel is.
            for (auto &p : s.probabilities)
                p *= 0.3 / total;
            s.probabilities[0] += 0.7;
        }
        steps.push_back(std::move(s));
    }
    return steps;
}

inline void Record(PauliPropagator &sim, Dense *dense, const Step &s)
{
    const Matrix x = Pauli(1, 1), z = Pauli(2, 1);
    const auto one = [&](const Matrix &u) {
        if (dense)
            dense->Unitary(u, s.q, 1);
    };
    const auto channel = [&](const std::vector<Matrix> &kraus, int arity) {
        if (dense)
            dense->Channel(kraus, s.q, arity);
    };
    const double angle = 4. * pi * (s.a - 0.5);
    switch (s.kind)
    {
    case 0:
        sim.ApplyH(s.q[0]);
        one((x + z) / std::sqrt(2.));
        break;
    case 1:
        sim.ApplyRX(s.q[0], angle);
        one(PauliNativeRegression::Rotation(1, angle));
        break;
    case 2:
        sim.ApplyRY(s.q[0], angle);
        one(PauliNativeRegression::Rotation(3, angle));
        break;
    case 3:
        sim.ApplyRZ(s.q[0], angle);
        one(PauliNativeRegression::Rotation(2, angle));
        break;
    case 4:
    case 5:
        // Embed puts local bit 0 on q[0]: the target, with q[1] the control.
        sim.ApplyCX(s.q[1], s.q[0]);
        if (dense)
            dense->Unitary(PauliNativeRegression::Controlled(x), s.q, 2);
        break;
    case 6:
        sim.ApplyPauliChannel(s.q[0], 0.2 * s.a, 0.2 * s.b, 0.2 * s.c);
        channel(PauliKraus({1. - 0.2 * (s.a + s.b + s.c), 0.2 * s.a, 0.2 * s.b, 0.2 * s.c}, 1), 1);
        break;
    case 7:
        sim.ApplyDepolarizingNoise(s.q[0], s.a < 0.1 ? 0.75 : 0.3 * s.a);
        {
            const double p = s.a < 0.1 ? 0.75 : 0.3 * s.a;
            channel(PauliKraus({1. - p, p / 3., p / 3., p / 3.}, 1), 1);
        }
        break;
    case 8:
        sim.ApplyTwoQubitPauliChannel(s.q[0], s.q[1], s.probabilities);
        channel(PauliKraus(std::vector<double>(s.probabilities.begin(), s.probabilities.end()), 2), 2);
        break;
    case 9:
        sim.ApplyPhaseDamping(s.q[0], s.a);
        channel(PhaseDampingKraus(s.a), 1);
        break;
    case 10:
        sim.ApplyAmplitudeDamping(s.q[0], s.a, s.b < 0.5 ? 0. : 0.4 * s.c);
        channel(DampingKraus(s.a, s.b < 0.5 ? 0. : 0.4 * s.c), 1);
        break;
    case 11:
        sim.ApplyKrausChannel(s.kraus, s.q[0]);
        channel(s.kraus, 1);
        break;
    case 12:
        sim.ApplyKrausChannel(s.kraus, s.q[0], s.q[1]);
        channel(s.kraus, 2);
        break;
    case 13:
        sim.ApplyReset(s.q[0], s.a < 0.3);
        channel(ResetKraus(s.a < 0.3), 1);
        break;
    default:
        sim.ApplyTwoQubitDepolarizingNoise(s.q[0], s.q[1], 0.4 * s.a);
        {
            std::vector<double> probabilities(16, 0.4 * s.a / 15.);
            probabilities[0] = 1. - 0.4 * s.a;
            channel(PauliKraus(probabilities, 2), 2);
        }
        break;
    }
}

// Random noisy circuits: every Pauli expectation and every basis probability
// against the dense mixed state, for each execution configuration.
inline void Circuits()
{
    std::mt19937_64 rng(977);
    size_t checked = 0;
    for (int trial = 0; trial < 24; ++trial)
    {
        const int width = 2 + trial % 3;
        const auto steps = RandomCircuit(rng, width, 14 + trial);
        for (int config = 0; config < 4; ++config)
        {
            PauliPropagator sim;
            sim.SetNrQubits(width);
            if (config == 1)
                sim.SetStepsBetweenDeduplication(1);
            if (config == 2)
                sim.SetStepsBetweenDeduplication(5);
            if (config == 3)
            {
                sim.EnableParallel(3);
                sim.SetParallelThreshold(1);
                sim.SetBatchSize(1);
                sim.SetStepsBetweenDeduplication(7);
            }
            Dense dense(width);
            for (const auto &s : steps)
                Record(sim, &dense, s);
            Close(dense.rho.trace().real(), 1., "reference lost trace");
            for (unsigned label = 0; label < (1U << (2 * width)); ++label)
            {
                Close(sim.ExpectationValue(PauliString(label, width)), dense.Expectation(label), "noisy expectation");
                ++checked;
            }
            for (size_t outcome = 0; outcome < (size_t(1) << width); ++outcome)
            {
                Close(sim.Probability(outcome), dense.Probability(outcome), "noisy probability", 1e-9);
                ++checked;
            }
            // Legacy operators round trip the recorded program.
            PauliPropagator copy;
            copy.SetNrQubits(width);
            copy.SetOperations(sim.GetOperations());
            const auto probe = PauliString(unsigned(trial * 7 + 5) & ((1U << (2 * width)) - 1), width);
            Close(copy.ExpectationValue(probe), sim.ExpectationValue(probe), "legacy round trip changed a channel",
                  1e-12);
            PauliPropagator shared;
            shared.SetNrQubits(width);
            shared.ShareOperationsFrom(sim);
            Close(shared.ExpectationValue(probe), sim.ExpectationValue(probe), "shared program changed a channel",
                  1e-12);
        }
    }
    std::cout << "PASS " << checked << " noisy expectations and probabilities against a dense mixed state\n";
}

// A subclass keeps its virtual Apply (an imported exact type would use the
// packed kernel), so it follows the unbatched custom-operator path.
struct SlowPauliChannel : OperatorPauliChannel
{
    using OperatorPauliChannel::OperatorPauliChannel;

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<SlowPauliChannel>(*this);
    }
};

struct SlowReset : OperatorReset
{
    using OperatorReset::OperatorReset;

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<SlowReset>(*this);
    }
};

// The packed kernels, the non-branching runs they share with Clifford gates,
// and the fused one-qubit blocks, against the per-operation virtual path.
inline void Engine()
{
    std::mt19937_64 rng(4242);
    for (int width : {7, 66, 130})
    {
        // Wide registers put the active qubits on both sides of a word boundary.
        const int base = width - 7;
        PauliPropagator fast;
        fast.SetNrQubits(width);
        std::uniform_real_distribution<double> uniform(0., 1.);
        for (int layer = 0; layer < 4; ++layer)
            for (int j = 0; j < 7; ++j)
            {
                const int q = base + j, next = base + (j + 1) % 7;
                fast.ApplyRY(q, 0.3 + 0.1 * uniform(rng));
                fast.ApplyDepolarizingNoise(q, 0.02 + 0.05 * uniform(rng));
                fast.ApplyCX(q, next);
                fast.ApplyTwoQubitDepolarizingNoise(q, next, 0.03);
                fast.ApplyH(q);
                fast.ApplyPhaseDamping(q, 0.1 * uniform(rng));
                if ((layer * 7 + j) % 9 == 4)
                    fast.ApplyReset(q, j % 2 == 0);
                if ((layer * 7 + j) % 11 == 3)
                    fast.ApplyDepolarizingNoise(next, 0.75);
            }
        auto operations = fast.GetOperations();
        for (auto &op : operations)
        {
            if (IsPauliDiagonalChannel(op->GetType()))
            {
                const auto &channel = static_cast<const OperatorPauliChannel &>(*op);
                if (op->GetNrQubits() == 1)
                    op = std::make_unique<SlowPauliChannel>(op->GetQubit(0), channel.GetEigenvalues());
                else
                    op = std::make_unique<SlowPauliChannel>(op->GetQubit(0), op->GetQubit(1), channel.GetEigenvalues());
            }
            else if (op->GetType() == OperationType::RESET)
                op = std::make_unique<SlowReset>(op->GetQubit(0), static_cast<const OperatorReset &>(*op).IsToOne());
        }
        PauliPropagator slow;
        slow.SetNrQubits(width);
        slow.SetOperations(std::move(operations));
        for (int config = 0; config < 5; ++config)
        {
            for (auto *sim : {&fast, &slow})
            {
                sim->SetStepsBetweenDeduplication(config == 1 ? 1 : config >= 2 ? 6 : std::numeric_limits<int>::max());
                // Truncated runs must take the same decisions on both paths.
                sim->SetCoefficientThreshold(config == 3 ? 1e-4 : 0.);
                sim->SetPauliWeightThreshold(config == 4 ? 4 : std::numeric_limits<size_t>::max());
                if (config == 2)
                {
                    sim->EnableParallel(3);
                    sim->SetParallelThreshold(1);
                    sim->SetBatchSize(1);
                }
                else
                    sim->DisableParallel();
            }
            for (int probe = 0; probe < 12; ++probe)
            {
                std::string pauli(static_cast<size_t>(width), 'I');
                for (int j = 0; j < 7; ++j)
                    pauli[base + j] = "IXYZ"[(probe * 5 + j * (probe + 1)) % 4];
                Close(fast.ExpectationValue(pauli), slow.ExpectationValue(pauli), "packed channel kernel", 1e-12);
            }
        }
    }

    // Large exact frontier: the rotations of one qubit, with the channels (and
    // a reset) between them, fuse into one table. The final rotation layer is
    // visited first and grows the frontier past the fusion size; the leading
    // Hadamards keep the other qubits unfinished meanwhile. Deduplicating at
    // every operation makes each index a checkpoint, which prevents fusion.
    PauliPropagator fused;
    const int width = 16;
    fused.SetNrQubits(width);
    for (int q = 0; q < width - 2; ++q)
        fused.ApplyH(q);
    for (int q = width - 2; q < width; ++q)
    {
        if (q == width - 2)
            fused.ApplyReset(q, true);
        fused.ApplyRX(q, 0.3);
        fused.ApplyDepolarizingNoise(q, 0.05);
        fused.ApplyRZ(q, -0.7);
        fused.ApplyPhaseDamping(q, 0.2);
        fused.ApplyRY(q, 0.9);
        fused.ApplyPauliChannel(q, 0.01, 0.02, 0.03);
    }
    for (int q = 0; q < width; ++q)
        fused.ApplyCX(q, (q + 1) % width);
    for (int q = 0; q < width; ++q)
        fused.ApplyRY(q, 0.4 + 0.01 * q);
    PauliPropagator unfused;
    unfused.SetNrQubits(width);
    unfused.ShareOperationsFrom(fused);
    unfused.SetStepsBetweenDeduplication(1);
    for (int workers : {0, 3})
    {
        if (workers)
            fused.EnableParallel(workers);
        for (const char *pauli : {"ZZZZZZZZZZZZZZZZ", "XZXZXZXZXZXZXZXZ", "XXXXXXXXXXXXXXXX"})
            Close(fused.ExpectationValue(std::string(pauli)), unfused.ExpectationValue(std::string(pauli)),
                  "fused channel block", 1e-11);
    }

    // A channel on an idle wire is outside the cone of the observable.
    PauliPropagator idle;
    idle.SetNrQubits(70);
    idle.ApplyRY(0, 0.4);
    idle.ApplyCX(0, 1);
    const double before = idle.ExpectationValue("ZI");
    PauliStringStorage reference;
    idle.ExpectationValue(std::string("ZI"), reference);
    idle.ApplyAmplitudeDamping(69, 0.5);
    idle.ApplyKrausChannel(RandomKraus(rng, 2, 3), 40, 68);
    idle.ApplyReset(33, true);
    idle.ApplyDepolarizingNoise(2, 0.75);
    Require(idle.ExpectationValue("ZI") == before, "idle channel changed an unrelated observable");
    PauliStringStorage expansion;
    idle.ExpectationValue(std::string("ZI"), expansion);
    Require(reference.size() == 2 && expansion.size() == reference.size(), "idle channel added terms");
    std::cout << "PASS packed kernels, non-branching runs, fused blocks and idle-wire cones\n";
}

// Frontier sizes. Diagonal channels never add a term; a zero eigenvalue and a
// reset remove terms; damping adds one term per damped Z, merged by
// deduplication when the labels coincide.
inline void Frontier()
{
    const auto terms = [](PauliPropagator &sim, const std::string &pauli) {
        PauliStringStorage expansion;
        sim.ExpectationValue(pauli, expansion);
        size_t count = 0;
        for (const auto &term : expansion)
            count += term.Coefficient != 0.;
        return count;
    };
    const int width = 6;
    const auto rotated = [&] {
        PauliPropagator sim;
        sim.SetNrQubits(width);
        for (int q = 0; q < width; ++q)
            sim.ApplyRY(q, 0.5);
        return sim;
    };
    const std::string z(static_cast<size_t>(width), 'Z');
    auto noiseless = rotated();
    const size_t ideal = terms(noiseless, z);
    Require(ideal == 64, "unexpected noiseless frontier");
    // Diagonal channels on both sides of the rotations: visited on the single
    // input term and on the full frontier.
    PauliPropagator diagonal;
    diagonal.SetNrQubits(width);
    const auto channels = [&] {
        for (int q = 0; q < width; ++q)
        {
            diagonal.ApplyDepolarizingNoise(q, 0.1);
            diagonal.ApplyPhaseDamping(q, 0.2);
            diagonal.ApplyTwoQubitDepolarizingNoise(q, (q + 1) % width, 0.05);
        }
    };
    channels();
    for (int q = 0; q < width; ++q)
        diagonal.ApplyRY(q, 0.5);
    channels();
    Require(terms(diagonal, z) == ideal, "a diagonal channel changed the frontier size");
    // Complete depolarization before the rotations acts on the whole frontier:
    // every term carries X or Z on that qubit and is removed.
    PauliPropagator erased;
    erased.SetNrQubits(width);
    erased.ApplyDepolarizingNoise(0, 0.75);
    for (int q = 0; q < width; ++q)
        erased.ApplyRY(q, 0.5);
    Require(terms(erased, z) == 0, "erased observable kept terms");
    Require(terms(erased, "IZZZZZ") == 32, "zero eigenvalue removed the wrong terms");
    // Reset after the rotations: Z becomes the identity, nothing is added.
    auto reset = rotated();
    reset.ApplyReset(0);
    Require(terms(reset, z) == 32, "reset changed the frontier size");
    Close(reset.ExpectationValue(z), std::pow(std::cos(0.5), 5), "reset expectation");

    // One layer of damping on a weight-w Z string: 2^w terms from one.
    for (int layers : {1, 3})
        for (bool dedup : {false, true})
        {
            PauliPropagator damped;
            damped.SetNrQubits(width);
            if (dedup)
                damped.SetStepsBetweenDeduplication(1);
            for (int layer = 0; layer < layers; ++layer)
                for (int q = 0; q < width; ++q)
                    damped.ApplyAmplitudeDamping(q, 0.1);
            // Without deduplication each layer leaves one more identity branch
            // per qubit: (layers + 1)^w. Merged, the labels are the 2^w subsets.
            const size_t expected = dedup ? 64 : static_cast<size_t>(std::pow(layers + 1, width) + 0.5);
            Require(terms(damped, z) == expected, "damping frontier");
            Close(damped.ExpectationValue(z), 1., "damped ground state");
        }
    std::cout << "PASS frontier sizes: diagonal channels, zero eigenvalues, reset and damping\n";
}

// Terminal sampling, probabilities and recorded measurements of a mixed state.
inline void Sampling()
{
    std::mt19937_64 rng(31337);
    const int width = 4;
    const auto steps = RandomCircuit(rng, width, 40);
    for (int config = 0; config < 3; ++config)
    {
        PauliPropagator sim;
        sim.SetNrQubits(width);
        sim.SetSeed(1234 + config);
        // 0: Z-moment table. 1: prefix tree only. 2: no persistent cache.
        if (config == 1)
            sim.SetSamplingCacheMaxNodes(8);
        if (config == 2)
            sim.SetSamplingCacheMaxNodes(0);
        Dense dense(width);
        for (const auto &s : steps)
            Record(sim, &dense, s);
        std::vector<int> qubits{2, 0, 3, 1};
        const size_t shots = config == 2 ? 3000 : 40000;
        const auto counts = sim.SampleCounts(qubits, shots);
        std::vector<double> frequency(size_t(1) << width);
        size_t total = 0;
        for (const auto &[bits, count] : counts)
        {
            size_t outcome = 0;
            for (size_t i = 0; i < bits.size(); ++i)
                if (bits[i])
                    outcome |= size_t(1) << qubits[i];
            frequency[outcome] += double(count) / shots;
            total += count;
        }
        Require(total == shots, "sampled shot count");
        for (size_t outcome = 0; outcome < frequency.size(); ++outcome)
        {
            const double p = dense.Probability(outcome);
            Require(std::abs(frequency[outcome] - p) < 5. * std::sqrt(p * (1. - p) / shots) + 2e-3,
                    "sampled distribution of a mixed state");
            Close(sim.Probability(outcome), p, "probability after sampling", 1e-9);
        }
        // Sampling neither consumed nor disturbed the recorded program.
        Close(sim.ExpectationValue("ZIIZ"), dense.Expectation(2 | (2 << 6)), "expectation after sampling");
    }

    // A recorded measurement of a mixed state: the projector with the measured
    // probability is the exact conditional state.
    for (int trial = 0; trial < 6; ++trial)
    {
        PauliPropagator sim;
        sim.SetNrQubits(3);
        sim.SetSeed(99 + trial);
        Dense dense(3);
        for (const auto &s : RandomCircuit(rng, 3, 20))
            Record(sim, &dense, s);
        const bool one = sim.Measure({1})[0];
        Matrix projector = Matrix::Zero(2, 2);
        projector(one ? 1 : 0, one ? 1 : 0) = 1.;
        dense.Channel({projector}, {{1, 0, 0}}, 1);
        const double probability = dense.rho.trace().real();
        Require(probability > 1e-12, "measured an impossible outcome");
        dense.rho /= probability;
        for (unsigned label = 0; label < 64; ++label)
            Close(sim.ExpectationValue(PauliString(label, 3)), dense.Expectation(label), "post-measurement mixed state",
                  1e-9);
        // Noise after the measurement still acts on the conditional state.
        sim.ApplyAmplitudeDamping(1, 0.3);
        dense.Channel(DampingKraus(0.3, 0.), {{1, 0, 0}}, 1);
        Close(sim.ExpectationValue("IZI"), dense.Expectation(2 << 2), "channel after a measurement", 1e-9);
    }

    // Reset is deterministic: no outcome is drawn and the seed is irrelevant.
    PauliPropagator a, b;
    for (auto *sim : {&a, &b})
    {
        sim->SetNrQubits(2);
        sim->ApplyH(0);
        sim->ApplyCX(0, 1);
        sim->ApplyReset(0);
        sim->ApplyH(0);
    }
    a.SetSeed(1);
    b.SetSeed(2);
    Require(a.GetOperations().size() == 4, "reset recorded more than one operation");
    for (const char *pauli : {"XI", "IZ", "XZ", "ZZ"})
        Require(a.ExpectationValue(pauli) == b.ExpectationValue(pauli), "reset depended on the seed");
    Close(a.ExpectationValue("XI"), 1., "reset then Hadamard");
    Close(a.ExpectationValue("IZ"), 0., "partner of a reset Bell pair is maximally mixed");

    // Checkpoints and clones keep channels.
    PauliPropagator checkpoint;
    checkpoint.SetNrQubits(2);
    checkpoint.ApplyRY(0, 0.8);
    checkpoint.ApplyAmplitudeDamping(0, 0.25);
    checkpoint.SaveState();
    const double saved = checkpoint.ExpectationValue("ZI");
    checkpoint.ApplyDepolarizingNoise(0, 0.5);
    checkpoint.ApplyReset(0, true);
    Close(checkpoint.ExpectationValue("ZI"), -1., "reset to one");
    checkpoint.RestoreState();
    Require(checkpoint.ExpectationValue("ZI") == saved, "restore lost a channel");
    std::cout << "PASS mixed-state sampling, probabilities, measurements, reset and checkpoints\n";
}

inline void Run()
{
    Actions();
    Recording();
    Circuits();
    Engine();
    Frontier();
    Sampling();
}
} // namespace PauliNoiseRegression
} // namespace QC
