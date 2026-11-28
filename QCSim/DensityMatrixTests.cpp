#include <array>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "StructuredGateRegressionTests.h"
#include "Tests.h"

#include "DensityMatrix.h"
#include "QubitRegister.h"

#include <unsupported/Eigen/KroneckerProduct>

#define _USE_MATH_DEFINES
#include <math.h>

#define DM_NR_QUBITS_LIMIT 6

static void DM_FillOneQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates)
{
    gates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::HyGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::SGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::SDGGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::TGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::TDGGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::PhaseShiftGate<>>(0.38 * M_PI));
    gates.emplace_back(std::make_shared<QC::Gates::PauliXGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::PauliYGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::PauliZGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTDagGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::RxGate<>>(M_PI / 3));
    gates.emplace_back(std::make_shared<QC::Gates::RyGate<>>(M_PI / 7));
    gates.emplace_back(std::make_shared<QC::Gates::RzGate<>>(M_PI / 5));
    gates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 3, M_PI / 5));
}

static void DM_FillTwoQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates)
{
    gates.emplace_back(std::make_shared<QC::Gates::SwapGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::iSwapDagGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::DecrementGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledYGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledZGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledHadamardGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledPhaseShiftGate<>>(M_PI / 3));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledUGate<>>(M_PI / 3, M_PI / 7));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledRxGate<>>(M_PI / 5));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledRyGate<>>(M_PI / 3));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledRzGate<>>(M_PI / 7));
}

static Eigen::MatrixXcd DM_FourierMatrix(Eigen::Index dimension)
{
    Eigen::MatrixXcd fourier(dimension, dimension);
    const double scale = 1. / std::sqrt(static_cast<double>(dimension));
    for (Eigen::Index row = 0; row < dimension; ++row)
        for (Eigen::Index col = 0; col < dimension; ++col)
            fourier(row, col) =
                scale * std::polar(1., 2. * M_PI * static_cast<double>(row * col) / static_cast<double>(dimension));
    return fourier;
}

static void DM_FillThreeQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates)
{
    gates.emplace_back(std::make_shared<QC::Gates::ToffoliGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::FredkinGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::CCZGate<>>());

    // A dense complex gate exercises the generic three-qubit kernel rather than a specialized path.
    gates.emplace_back(std::make_shared<QC::Gates::ThreeQubitsGate<>>(DM_FourierMatrix(8)));
}

template <class Callable>
static bool DM_ExpectInvalidArgument(Callable &&callable, const char *description, const char *messagePart = nullptr)
{
    try
    {
        callable();
    }
    catch (const std::invalid_argument &ex)
    {
        if (messagePart && std::string(ex.what()).find(messagePart) == std::string::npos)
        {
            std::cout << description << " reported an unrelated error: " << ex.what() << std::endl;
            return false;
        }
        return true;
    }
    catch (const std::exception &ex)
    {
        std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
        return false;
    }

    std::cout << description << " did not throw std::invalid_argument" << std::endl;
    return false;
}

template <class Callable> static bool DM_ExpectDomainError(Callable &&callable, const char *description)
{
    try
    {
        callable();
    }
    catch (const std::domain_error &)
    {
        return true;
    }
    catch (const std::exception &ex)
    {
        std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
        return false;
    }

    std::cout << description << " did not throw std::domain_error" << std::endl;
    return false;
}

// checks that a density matrix equals the outer product |psi><psi| of the statevector simulator
static bool CompareWithStatevector(const QC::DensityMatrix<> &dm, const QC::QubitRegister<> &reg, int nrQubits,
                                   double eps = 1E-9)
{
    const Eigen::MatrixXcd expected = reg.getDensityMatrix();
    const Eigen::MatrixXcd &got = dm.getDensityMatrix();

    if (expected.rows() != got.rows() || expected.cols() != got.cols())
    {
        std::cout << "Density matrix dimensions mismatch for " << nrQubits << " qubits" << std::endl;
        return false;
    }

    for (int i = 0; i < expected.rows(); ++i)
        for (int j = 0; j < expected.cols(); ++j)
            if (!approxEqual(expected(i, j), got(i, j), eps))
            {
                std::cout << "Density matrix element (" << i << ", " << j << ") mismatch for " << nrQubits
                          << " qubits: " << expected(i, j) << " (statevector) vs " << got(i, j) << " (density matrix)"
                          << std::endl;
                return false;
            }

    return true;
}

// applies random unitary circuits to both a statevector simulator and the density matrix simulator
// and checks that rho == |psi><psi| after every gate
static bool DensityMatrixUnitaryTest()
{
    std::cout << "\nDensity matrix simulator - unitary circuits compared against the statevector simulator"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(25, 50);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        for (int t = 0; t < 10; ++t)
        {
            QC::DensityMatrix<> dm(nrQubits);
            QC::QubitRegister<> reg(nrQubits);

            const int lim = nrGatesDistr(gen);

            for (int i = 0; i < lim; ++i)
            {
                const int g = gateDistr(gen);
                const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
                if (twoQubitsGate && nrQubits < 2)
                    continue;

                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = qubit1 + 1;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

                dm.ApplyGate(*gates[g], qubit1, qubit2);
                reg.ApplyGate(*gates[g], qubit1, qubit2);

                if (!CompareWithStatevector(dm, reg, nrQubits))
                    return false;
            }

            // sanity: a pure state must have unit trace and purity 1
            if (!approxEqual(dm.Trace().real(), 1., 1E-9) || !approxEqual(dm.Trace().imag(), 0., 1E-9))
            {
                std::cout << "Trace is not 1 for a pure state: " << dm.Trace() << std::endl;
                return false;
            }

            if (!approxEqual(dm.Purity(), 1., 1E-9))
            {
                std::cout << "Purity is not 1 for a pure state: " << dm.Purity() << std::endl;
                return false;
            }

            if (!dm.IsHermitian(1E-9))
            {
                std::cout << "Density matrix is not Hermitian" << std::endl;
                return false;
            }
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

// Independently embeds a small operator in the full Hilbert space. Local bit zero corresponds to
// qubits[0], local bit one to qubits[1], and so on, matching the simulator's gate convention.
static Eigen::MatrixXcd DM_EmbedOperator(const Eigen::MatrixXcd &op, const std::vector<size_t> &qubits, size_t nrQubits)
{
    const size_t nrBasisStates = 1ULL << nrQubits;
    size_t operatedMask = 0;
    for (const size_t qubit : qubits)
        operatedMask |= 1ULL << qubit;
    const size_t unaffectedMask = (nrBasisStates - 1) ^ operatedMask;

    Eigen::MatrixXcd embedded =
        Eigen::MatrixXcd::Zero(static_cast<Eigen::Index>(nrBasisStates), static_cast<Eigen::Index>(nrBasisStates));
    for (size_t row = 0; row < nrBasisStates; ++row)
        for (size_t col = 0; col < nrBasisStates; ++col)
        {
            if ((row & unaffectedMask) != (col & unaffectedMask))
                continue;

            size_t localRow = 0;
            size_t localCol = 0;
            for (size_t localQubit = 0; localQubit < qubits.size(); ++localQubit)
            {
                localRow |= ((row >> qubits[localQubit]) & 1ULL) << localQubit;
                localCol |= ((col >> qubits[localQubit]) & 1ULL) << localQubit;
            }
            embedded(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) =
                op(static_cast<Eigen::Index>(localRow), static_cast<Eigen::Index>(localCol));
        }

    return embedded;
}

static Eigen::MatrixXcd DM_DeterministicMixedState(size_t nrQubits)
{
    const Eigen::Index dimension = static_cast<Eigen::Index>(1ULL << nrQubits);
    Eigen::VectorXcd psi1(dimension);
    Eigen::VectorXcd psi2(dimension);
    for (Eigen::Index i = 0; i < dimension; ++i)
    {
        psi1(i) = std::complex<double>(static_cast<double>(i + 1), static_cast<double>(i % 3) - 1.);
        psi2(i) =
            std::complex<double>(static_cast<double>((2 * i + 1) % 7) - 3., static_cast<double>((3 * i + 2) % 5) - 2.);
    }
    psi1.normalize();
    psi2.normalize();
    return 0.37 * (psi1 * psi1.adjoint()) + 0.63 * (psi2 * psi2.adjoint());
}

// Cross-checks the column/row gate kernels against an independently embedded dense operator on a
// genuinely mixed state. Both specialized gate objects and recorded AppliedGate objects are used.
static bool DensityMatrixDenseGateReferenceTest()
{
    std::cout << "\nDensity matrix simulator - gates compared against independently embedded dense operators"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> oneQubitGates;
    oneQubitGates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
    oneQubitGates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 5., M_PI / 7.));

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> twoQubitGates;
    twoQubitGates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
    twoQubitGates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
    twoQubitGates.emplace_back(std::make_shared<QC::Gates::iSwapDagGate<>>());
    twoQubitGates.emplace_back(std::make_shared<QC::Gates::TwoQubitsGate<>>(DM_FourierMatrix(4)));

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> threeQubitGates;
    DM_FillThreeQubitGates(threeQubitGates);

    constexpr size_t nrQubits = 4;
    const Eigen::MatrixXcd initial = DM_DeterministicMixedState(nrQubits);
    auto checkPlacement = [&](const QC::Gates::QuantumGateWithOp<> &gate, const std::vector<size_t> &qubits) -> bool {
        const Eigen::MatrixXcd embedded = DM_EmbedOperator(gate.getRawOperatorMatrix(), qubits, nrQubits);
        const Eigen::MatrixXcd expected = embedded * initial * embedded.adjoint();

        QC::DensityMatrix<> direct(nrQubits);
        direct.setDensityMatrix(initial);
        direct.ApplyGate(gate, qubits[0], qubits.size() > 1 ? qubits[1] : 0, qubits.size() > 2 ? qubits[2] : 0);

        QC::DensityMatrix<> recorded(nrQubits);
        recorded.setDensityMatrix(initial);
        recorded.ApplyGate(QC::Gates::AppliedGate<>(gate.getRawOperatorMatrix(), qubits[0],
                                                    qubits.size() > 1 ? qubits[1] : 0,
                                                    qubits.size() > 2 ? qubits[2] : 0));

        if ((direct.getDensityMatrix() - expected).norm() > 1E-9 ||
            (recorded.getDensityMatrix() - expected).norm() > 1E-9)
        {
            std::cout << "Dense gate reference mismatch for a " << qubits.size() << "-qubit gate on";
            for (const size_t qubit : qubits)
                std::cout << " " << qubit;
            std::cout << std::endl;
            return false;
        }
        return true;
    };

    for (const auto &gate : oneQubitGates)
        for (size_t q1 = 0; q1 < nrQubits; ++q1)
            if (!checkPlacement(*gate, {q1}))
                return false;

    for (const auto &gate : twoQubitGates)
        for (size_t q1 = 0; q1 < nrQubits; ++q1)
            for (size_t q2 = 0; q2 < nrQubits; ++q2)
                if (q1 != q2 && !checkPlacement(*gate, {q1, q2}))
                    return false;

    for (const auto &gate : threeQubitGates)
        for (size_t q1 = 0; q1 < nrQubits; ++q1)
            for (size_t q2 = 0; q2 < nrQubits; ++q2)
                for (size_t q3 = 0; q3 < nrQubits; ++q3)
                    if (q1 != q2 && q1 != q3 && q2 != q3 && !checkPlacement(*gate, {q1, q2, q3}))
                        return false;

    std::cout << "Success" << std::endl;
    return true;
}

// Exercises every specialized three-qubit path plus a dense complex 8x8 gate. Every ordered
// placement is checked, including non-adjacent placements in registers larger than three qubits.
static bool DensityMatrixThreeQubitGateTest()
{
    std::cout << "\nDensity matrix simulator - three qubit gates compared against the statevector simulator"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillThreeQubitGates(gates);

    QC::Gates::HadamardGate<> h;
    QC::Gates::RyGate<> ry(M_PI / 7.);
    QC::Gates::CNOTGate<> cnot;

    for (int nrQubits = 3; nrQubits <= 5; ++nrQubits)
    {
        for (const auto &gate : gates)
        {
            for (int q1 = 0; q1 < nrQubits; ++q1)
                for (int q2 = 0; q2 < nrQubits; ++q2)
                    for (int q3 = 0; q3 < nrQubits; ++q3)
                    {
                        if (q1 == q2 || q1 == q3 || q2 == q3)
                            continue;

                        QC::DensityMatrix<> dm(nrQubits);
                        QC::QubitRegister<> reg(nrQubits);
                        for (int q = 0; q < nrQubits; ++q)
                        {
                            if (q % 2 == 0)
                            {
                                dm.ApplyGate(h, q);
                                reg.ApplyGate(h, q);
                            }
                            else
                            {
                                dm.ApplyGate(ry, q);
                                reg.ApplyGate(ry, q);
                            }
                        }
                        for (int q = 1; q < nrQubits; ++q)
                        {
                            dm.ApplyGate(cnot, q, q - 1);
                            reg.ApplyGate(cnot, q, q - 1);
                        }

                        dm.ApplyGate(*gate, q1, q2, q3);
                        reg.ApplyGate(*gate, q1, q2, q3);
                        if (!CompareWithStatevector(dm, reg, nrQubits))
                            return false;
                    }
        }
    }

    std::cout << "Success" << std::endl;
    return true;
}

// checks the noise channels: trace preservation, that noise turns a pure state into a mixed one,
// and the analytic amplitude damping result
static bool DensityMatrixChannelsTest()
{
    std::cout << "\nDensity matrix simulator - quantum channels" << std::endl;

    // amplitude damping on |1>: rho' = [[gamma, 0], [0, 1 - gamma]]
    {
        const double gamma = 0.3;
        QC::DensityMatrix<> dm(1);
        QC::Gates::PauliXGate<> x;
        dm.ApplyGate(x, 0); // |0> -> |1>
        dm.ApplyAmplitudeDamping(0, gamma);

        const auto &rho = dm.getDensityMatrix();
        if (!approxEqual(rho(0, 0), std::complex<double>(gamma, 0), 1E-9) ||
            !approxEqual(rho(1, 1), std::complex<double>(1. - gamma, 0), 1E-9) ||
            !approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
            !approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
        {
            std::cout << "Amplitude damping result is wrong:\n" << rho << std::endl;
            return false;
        }
    }

    // dephasing of |+> suppresses the off diagonal coherences while keeping the populations
    {
        QC::DensityMatrix<> dm(1);
        QC::Gates::HadamardGate<> h;
        dm.ApplyGate(h, 0); // |+>

        dm.ApplyPhaseDamping(0, 1.0); // full dephasing

        const auto &rho = dm.getDensityMatrix();
        if (!approxEqual(rho(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho(1, 1), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
            !approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
        {
            std::cout << "Full dephasing result is wrong:\n" << rho << std::endl;
            return false;
        }

        if (approxEqual(dm.Purity(), 1., 1E-9))
        {
            std::cout << "Dephased state should be mixed but purity is 1" << std::endl;
            return false;
        }
    }

    // depolarizing, bit flip and phase flip must all preserve the trace
    {
        std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> oneQubitGates;
        DM_FillOneQubitGates(oneQubitGates);

        std::uniform_int_distribution gateDistr(0, static_cast<int>(oneQubitGates.size()) - 1);

        for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
        {
            std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

            QC::DensityMatrix<> dm(nrQubits);

            for (int i = 0; i < 20; ++i)
            {
                const int q = qubitDistr(gen);
                dm.ApplyGate(*oneQubitGates[gateDistr(gen)], q);
                dm.ApplyDepolarizingNoise(q, 0.1);
                dm.ApplyBitFlipNoise(q, 0.05);
                dm.ApplyPhaseFlipNoise(q, 0.05);

                if (!approxEqual(dm.Trace().real(), 1., 1E-9) || !approxEqual(dm.Trace().imag(), 0., 1E-9))
                {
                    std::cout << "A channel did not preserve the trace: " << dm.Trace() << std::endl;
                    return false;
                }

                if (!dm.IsHermitian(1E-9))
                {
                    std::cout << "A channel broke Hermiticity" << std::endl;
                    return false;
                }
            }
        }
    }

    // reset must always bring a qubit to |0>
    {
        QC::DensityMatrix<> dm(1);
        QC::Gates::HadamardGate<> h;
        dm.ApplyGate(h, 0);
        dm.ApplyReset(0);

        if (!approxEqual(dm.GetQubitProbability(0), 0., 1E-9))
        {
            std::cout << "Reset did not bring the qubit to |0>: p1 = " << dm.GetQubitProbability(0) << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

static bool DensityMatrixGenericChannelTest()
{
    std::cout << "\nDensity matrix simulator - generic Kraus channels against explicit dense evolution" << std::endl;

    // One-qubit amplitude damping embedded in the middle of a three-qubit mixed state.
    {
        QC::DensityMatrix<> dm(3);
        const Eigen::MatrixXcd before = DM_DeterministicMixedState(3);
        dm.setDensityMatrix(before);

        constexpr double gamma = 0.29;
        Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(2, 2);
        k0(0, 0) = 1.;
        k0(1, 1) = std::sqrt(1. - gamma);
        Eigen::MatrixXcd k1 = Eigen::MatrixXcd::Zero(2, 2);
        k1(0, 1) = std::sqrt(gamma);

        const Eigen::MatrixXcd fullK0 = DM_EmbedOperator(k0, {1}, 3);
        const Eigen::MatrixXcd fullK1 = DM_EmbedOperator(k1, {1}, 3);
        const Eigen::MatrixXcd expected = fullK0 * before * fullK0.adjoint() + fullK1 * before * fullK1.adjoint();
        dm.ApplyChannel({k0, k1}, 1);
        if ((dm.getDensityMatrix() - expected).norm() > 1E-10)
        {
            std::cout << "Generic one-qubit Kraus evolution does not match the explicit matrix result" << std::endl;
            return false;
        }
    }

    // A correlated channel on non-adjacent qubits zero and two.
    {
        QC::DensityMatrix<> dm(3);
        const Eigen::MatrixXcd before = DM_DeterministicMixedState(3);
        dm.setDensityMatrix(before);
        QC::Gates::CNOTGate<> cnot;
        constexpr double p = 0.37;
        const Eigen::MatrixXcd k0 = std::sqrt(1. - p) * Eigen::MatrixXcd::Identity(4, 4);
        const Eigen::MatrixXcd k1 = std::sqrt(p) * cnot.getRawOperatorMatrix();
        const Eigen::MatrixXcd fullK0 = DM_EmbedOperator(k0, {0, 2}, 3);
        const Eigen::MatrixXcd fullK1 = DM_EmbedOperator(k1, {0, 2}, 3);
        const Eigen::MatrixXcd expected = fullK0 * before * fullK0.adjoint() + fullK1 * before * fullK1.adjoint();

        dm.ApplyChannel({k0, k1}, 0, 2);
        if ((dm.getDensityMatrix() - expected).norm() > 1E-10 ||
            !approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-10) || !dm.IsHermitian(1E-10))
        {
            std::cout << "Generic non-adjacent two-qubit Kraus evolution does not match the explicit matrix result"
                      << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;
    return true;
}

// The d x d blocks of a random real or complex isometry form a trace-preserving Kraus family.
static std::vector<Eigen::MatrixXcd> DM_RandomKraus(Eigen::Index dimension, Eigen::Index count, std::mt19937_64 &random,
                                                    bool real = false)
{
    std::normal_distribution<double> normal;
    Eigen::MatrixXcd wide(count * dimension, dimension);
    for (Eigen::Index r = 0; r < wide.rows(); ++r)
        for (Eigen::Index c = 0; c < wide.cols(); ++c)
            wide(r, c) = std::complex<double>(normal(random), real ? 0. : normal(random));
    const Eigen::MatrixXcd isometry = Eigen::HouseholderQR<Eigen::MatrixXcd>(wide).householderQ() *
                                      Eigen::MatrixXcd::Identity(wide.rows(), wide.cols());
    std::vector<Eigen::MatrixXcd> kraus;
    for (Eigen::Index k = 0; k < count; ++k)
        kraus.emplace_back(isometry.block(k * dimension, 0, dimension, dimension));
    return kraus;
}

// sum_k conj(E_k) (x) E_k on a block indexed row + d * column.
static Eigen::MatrixXcd DM_Superoperator(const std::vector<Eigen::MatrixXcd> &kraus)
{
    const Eigen::Index d = kraus.front().rows();
    Eigen::MatrixXcd superoperator = Eigen::MatrixXcd::Zero(d * d, d * d);
    for (const auto &e : kraus)
        for (Eigen::Index output = 0; output < d * d; ++output)
            for (Eigen::Index input = 0; input < d * d; ++input)
                superoperator(output, input) += e(output % d, input % d) * std::conj(e(output / d, input / d));
    return superoperator;
}

// Channel kernels against sum_k E rho E^dagger with independently embedded operators:
// complex, real and Pauli maps, every target position (index bit zero takes a different vector path),
// reversed and non-adjacent pairs, below and above the parallel threshold.
static bool DensityMatrixSuperoperatorTest()
{
    std::cout << "\nDensity matrix simulator - local channels against explicit dense evolution" << std::endl;

    std::mt19937_64 random(0x5EED5EEDULL);
    const auto pauli = [](int id) {
        Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(2, 2);
        if (id == 0)
            m.setIdentity();
        else if (id == 1)
            m << 0., 1., 1., 0.;
        else if (id == 2)
            m << 0., std::complex<double>(0., -1.), std::complex<double>(0., 1.), 0.;
        else
            m << 1., 0., 0., -1.;
        return m;
    };

    for (const size_t nrQubits : {size_t{3}, size_t{8}})
    {
        const Eigen::MatrixXcd before = DM_DeterministicMixedState(nrQubits);
        const auto check = [&](const std::vector<Eigen::MatrixXcd> &kraus, const std::vector<size_t> &qubits,
                               const char *description) {
            Eigen::MatrixXcd expected = Eigen::MatrixXcd::Zero(before.rows(), before.cols());
            for (const auto &e : kraus)
            {
                const Eigen::MatrixXcd full = DM_EmbedOperator(e, qubits, nrQubits);
                expected += full * before * full.adjoint();
            }
            if (qubits.size() == 3)
            {
                using RowMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
                QC::DensityMatrix<Eigen::VectorXcd, RowMatrix> rowMajor(nrQubits);
                rowMajor.setDensityMatrix(before);
                rowMajor.ApplyChannel(std::vector<RowMatrix>(kraus.begin(), kraus.end()), qubits[0], qubits[1],
                                      qubits[2]);
                if (!rowMajor.getDensityMatrix().isApprox(expected, 1E-10))
                {
                    std::cout << description << " differs from the reference with row-major storage" << std::endl;
                    return false;
                }
            }
            for (int entry = 0; entry < (qubits.size() == 3 ? 1 : 2); ++entry)
            {
                QC::DensityMatrix<> dm(nrQubits);
                dm.setDensityMatrix(before);
                if (qubits.size() == 3)
                    dm.ApplyChannel(kraus, qubits[0], qubits[1], qubits[2]);
                else if (entry == 0)
                    dm.ApplyChannel(kraus, qubits[0], qubits.size() > 1 ? qubits[1] : 0);
                else
                    dm.ApplySuperoperator(DM_Superoperator(kraus), qubits[0], qubits.size() > 1 ? qubits[1] : 0);
                const double error = (dm.getDensityMatrix() - expected).norm();
                if (error > 1E-10 || !approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-10) ||
                    !dm.IsHermitian(1E-10))
                {
                    std::cout << description << (entry == 0 ? " (Kraus)" : " (superoperator)") << " on " << nrQubits
                              << " qubits, first target " << qubits[0] << ": error " << error << std::endl;
                    return false;
                }
            }
            return true;
        };

        for (size_t q = 0; q < nrQubits; ++q)
        {
            constexpr double gamma = 0.23;
            Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(2, 2), k1 = Eigen::MatrixXcd::Zero(2, 2);
            k0(0, 0) = 1.;
            k0(1, 1) = std::sqrt(1. - gamma);
            k1(0, 1) = std::sqrt(gamma);
            if (!check(DM_RandomKraus(2, 3, random), {q}, "Random one-qubit channel") ||
                !check({k0, k1}, {q}, "Amplitude damping") ||
                !check({std::sqrt(0.7) * pauli(0), std::sqrt(0.1) * pauli(1), std::sqrt(0.15) * pauli(2),
                        std::sqrt(0.05) * pauli(3)},
                       {q}, "One-qubit Pauli channel"))
                return false;
        }

        std::vector<Eigen::MatrixXcd> pauliPair;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                pauliPair.emplace_back(std::sqrt(a == 0 && b == 0 ? 0.85 : 0.01) *
                                       Eigen::kroneckerProduct(pauli(b), pauli(a)).eval());
        const size_t last = nrQubits - 1;
        for (const auto &pair :
             std::vector<std::vector<size_t>>{{0, 1}, {1, 0}, {0, last}, {last, 0}, {1, 2}, {last, 1}})
            if (!check(DM_RandomKraus(4, 5, random), pair, "Random two-qubit channel") ||
                !check(pauliPair, pair, "Two-qubit Pauli channel"))
                return false;
        if (!check(DM_RandomKraus(4, 2, random, true), {last, 0}, "Dense real two-qubit channel") ||
            !check(DM_RandomKraus(4, 2, random, true), {last, 1}, "Dense real two-qubit channel") ||
            !check(DM_RandomKraus(8, 3, random), {1, last, 0}, "Dense complex three-qubit channel") ||
            !check(DM_RandomKraus(8, 2, random, true), {last, 0, 1}, "Dense real three-qubit channel") ||
            !check(DM_RandomKraus(8, 1, random), {0, last, 1}, "Single three-qubit Kraus operator"))
            return false;
        if (nrQubits >= 5 && !check(DM_RandomKraus(8, 3, random, true), {last, 2, 3}, "Spread three-qubit channel"))
            return false;
    }

    // This trace-preserving map has a positive, diagonally dominant Choi matrix. On |+><+|,
    // the first three terms in this channel's off-diagonal row cancel exactly.
    QC::DensityMatrix<> ordered(1);
    const Eigen::MatrixXcd plus = Eigen::MatrixXcd::Constant(2, 2, 0.5);
    const double residual = std::ldexp(1., -60);
    Eigen::MatrixXcd map(4, 4);
    map.row(0) << .5, 1. / 64, 1. / 64, .5;
    map.row(1) << 1. / 16, 1. / 16, -1. / 8, residual;
    map.row(2) << 1. / 16, -1. / 8, 1. / 16, residual;
    map.row(3) << .5, -1. / 64, -1. / 64, .5;
    ordered.setDensityMatrix(plus);
    ordered.ApplySuperoperator(map, 0);
    if (ordered.getDensityMatrix()(1, 0) != 0.5 * residual || ordered.getDensityMatrix()(0, 1) != 0.5 * residual)
    {
        std::cout << "A dense channel lost a representable coherence after cancellation" << std::endl;
        return false;
    }
    // ApplySuperoperator also accepts finite linear maps. A reordered sum of
    // three positive terms would overflow before the negative term is added.
    const double large = 1.7E308;
    map.row(1) << large, large, -large, large;
    ordered.setDensityMatrix(plus);
    ordered.ApplySuperoperator(map, 0);
    if (ordered.getDensityMatrix()(1, 0) != large)
    {
        std::cout << "A finite dense local map overflowed during summation" << std::endl;
        return false;
    }

    QC::DensityMatrix<> dm(2);
    const Eigen::MatrixXcd identity4 = Eigen::MatrixXcd::Identity(4, 4);
    Eigen::MatrixXcd nonFinite = identity4;
    nonFinite(2, 1) = std::numeric_limits<double>::quiet_NaN();
    if (!DM_ExpectInvalidArgument([&] { dm.ApplySuperoperator(Eigen::MatrixXcd::Identity(8, 8), 0); },
                                  "Three-index superoperator") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplySuperoperator(Eigen::MatrixXcd::Identity(4, 16), 0); },
                                  "Non-square superoperator") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplySuperoperator(nonFinite, 0); }, "Non-finite superoperator") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplySuperoperator(identity4, 2); }, "Out-of-range superoperator qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplySuperoperator(Eigen::MatrixXcd::Identity(16, 16), 1, 1); },
                                  "Duplicate superoperator qubits"))
        return false;

    std::cout << "Success" << std::endl;
    return true;
}

// Reduced density matrix of `qubits` (the first is the low local bit) by explicit sums over the other bits.
static Eigen::MatrixXcd DM_ReducedDensity(const Eigen::VectorXcd &psi, const std::vector<size_t> &qubits)
{
    const size_t dimension = size_t{1} << qubits.size();
    size_t mask = 0;
    for (const size_t q : qubits)
        mask |= size_t{1} << q;
    const auto local = [&](size_t state) {
        size_t index = 0;
        for (size_t k = 0; k < qubits.size(); ++k)
            index |= ((state >> qubits[k]) & 1) << k;
        return index;
    };
    Eigen::MatrixXcd rho = Eigen::MatrixXcd::Zero(dimension, dimension);
    for (size_t row = 0; row < static_cast<size_t>(psi.size()); ++row)
        for (size_t l = 0; l < dimension; ++l)
        {
            size_t col = row & ~mask;
            for (size_t k = 0; k < qubits.size(); ++k)
                if (l & (size_t{1} << k))
                    col |= size_t{1} << qubits[k];
            rho(local(row), l) += psi(row) * std::conj(psi(col));
        }
    return rho;
}

// A statevector takes one Kraus operator per trajectory step, with probabilities read from the reduced
// density matrix of the targets. Enumerating the branches of a channel, each weighted by its probability,
// must rebuild the channel exactly: sum_k p_k |psi_k><psi_k| = sum_k E_k rho E_k^dagger.
static bool StatevectorKrausTrajectoryTest()
{
    std::cout << "\nStatevector simulator - Kraus trajectory steps against the density matrix" << std::endl;

    // The low-level calculator also accepts an explicit number of complete local groups.
    // A partial final panel must leave the NaN padding out of the reduction.
    for (const size_t groups : {63u, 64u, 65u, 127u, 128u})
    {
        Eigen::VectorXcd local(8);
        for (int i = 0; i < 8; ++i)
            local[i] = std::complex<double>(i + 1., 2. * i - 3.);
        Eigen::VectorXcd storage =
            Eigen::VectorXcd::Constant(8 * (groups + 63), std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 0; i < groups; ++i)
            storage.segment<8>(8 * i) = local;
        const auto rho = QC::QubitRegisterCalculator<>::LocalDensityMatrix<3>(storage, {1, 2, 4}, 8 * groups, false);
        for (int r = 0; r < 8; ++r)
            for (int c = 0; c < 8; ++c)
                if (rho[8 * r + c] != double(groups) * local[r] * std::conj(local[c]))
                {
                    std::cout << "Incorrect partial reduced-density panel for " << groups << " groups" << std::endl;
                    return false;
                }
    }

    // All sixteen output amplitudes have squared magnitude below one subnormal unit.
    // Together they have a representable norm, for both placements of the paired SIMD lanes.
    const double tinyScale = std::ldexp(1., -537);
    for (size_t first : {size_t{0}, size_t{3}})
        for (const auto coefficient : {std::complex<double>(tinyScale, 0.), std::complex<double>(0., tinyScale)})
        {
            QC::QubitRegister<> reg(4);
            reg.setToEqualSuperposition();
            const Eigen::MatrixXcd map = coefficient * Eigen::MatrixXcd::Identity(4, 4);
            const double weight = reg.ApplyLocalMapNorm(map, first, 1);
            if (weight != std::numeric_limits<double>::denorm_min() ||
                reg.getBasisStateAmplitude(0) != coefficient * 0.25)
            {
                std::cout << "A local map lost a representable subnormal squared norm" << std::endl;
                return false;
            }
        }

    // Some reduction blocks contain only zeros, and each square in the others rounds to zero.
    // The map must distinguish their representable total from exact cancellation.
    {
        constexpr size_t count = size_t{1} << 16;
        Eigen::VectorXcd input = Eigen::VectorXcd::Zero(count);
        input.tail(count / 2).setConstant(std::complex<double>(1. / 256., -1. / 256.));
        const Eigen::MatrixXcd tiny = std::complex<double>(0., tinyScale) * Eigen::MatrixXcd::Identity(4, 4);
        Eigen::MatrixXcd cancel = Eigen::MatrixXcd::Zero(4, 4);
        cancel(0, 0) = 1.;
        cancel(0, 1) = -1.;
        for (bool parallel : {false, true})
        {
            QC::QubitRegister<> reg(16);
            reg.SetMultithreading(parallel);
            reg.setRegisterStorage(input);
            if (reg.ApplyLocalMapNorm(tiny, 2, 5) != std::numeric_limits<double>::denorm_min())
            {
                std::cout << "A local map lost a subnormal norm distributed across blocks" << std::endl;
                return false;
            }
            reg.setRegisterStorage(input);
            if (reg.ApplyLocalMapNorm(cancel, 0, 5) != 0. ||
                (reg.getRegisterStorage().array() != std::complex<double>(0.)).any())
            {
                std::cout << "An exactly cancelling local map produced a nonzero result" << std::endl;
                return false;
            }
        }
    }

    std::mt19937_64 random(0xC0FFEEULL);
    std::normal_distribution<double> normal;
    const auto randomState = [&](size_t nrQubits) {
        Eigen::VectorXcd psi(Eigen::Index{1} << nrQubits);
        for (Eigen::Index i = 0; i < psi.size(); ++i)
            psi(i) = std::complex<double>(normal(random), normal(random));
        psi.normalize();
        return psi;
    };

    // The reduced density matrix, below and above the parallel threshold, for every target order.
    for (const size_t nrQubits : {size_t{4}, size_t{15}})
    {
        const Eigen::VectorXcd psi = randomState(nrQubits);
        QC::QubitRegister<> reg(nrQubits);
        reg.setRegisterStorage(psi);
        const size_t last = nrQubits - 1;
        for (const auto &qubits : std::vector<std::vector<size_t>>{
                 {0}, {1}, {last}, {0, 1}, {1, 0}, {last, 0}, {2, last}, {0, 1, 2}, {last, 0, 2}, {2, 3, 1}})
        {
            const Eigen::MatrixXcd expected = DM_ReducedDensity(psi, qubits);
            const auto evaluate = [&] {
                return qubits.size() == 1   ? reg.LocalDensityMatrix(qubits[0])
                       : qubits.size() == 2 ? reg.LocalDensityMatrix(qubits[0], qubits[1])
                                            : reg.LocalDensityMatrix(qubits[0], qubits[1], qubits[2]);
            };
            const auto checkPopulations = [&](const Eigen::MatrixXcd &rho) {
                const auto same = [&](const auto &populations) {
                    // The compiler can contract the scalar population and density reductions differently.
                    for (size_t i = 0; i < populations.size(); ++i)
                        if (!approxEqual(populations[i], rho(i, i).real(), 1E-14))
                            return false;
                    return true;
                };
                if (qubits.size() == 1)
                    return same(reg.LocalPopulations(qubits[0]));
                if (qubits.size() == 2)
                    return same(reg.LocalPopulations(qubits[0], qubits[1]));
                return same(reg.LocalPopulations(qubits[0], qubits[1], qubits[2]));
            };
            const Eigen::MatrixXcd parallel = evaluate();
            const bool parallelPopulations = checkPopulations(parallel);
            reg.SetMultithreading(false);
            const Eigen::MatrixXcd serial = evaluate();
            const bool serialPopulations = checkPopulations(serial);
            reg.SetMultithreading(true);
            if ((parallel - expected).norm() > 1E-12 || parallel != serial || !parallelPopulations ||
                !serialPopulations)
            {
                std::cout << "Reduced density matrix of " << qubits.size() << " qubits is wrong on " << nrQubits
                          << " qubits (error " << (parallel - expected).norm()
                          << "), or depends on the number of threads" << std::endl;
                return false;
            }
        }
    }

    constexpr size_t nrQubits = 4;
    const Eigen::VectorXcd psi = randomState(nrQubits);
    const Eigen::MatrixXcd rho = psi * psi.adjoint();
    constexpr double gamma = 0.31;
    Eigen::MatrixXcd damp0 = Eigen::MatrixXcd::Zero(2, 2), damp1 = Eigen::MatrixXcd::Zero(2, 2);
    damp0(0, 0) = 1.;
    damp0(1, 1) = std::sqrt(1. - gamma);
    damp1(0, 1) = std::sqrt(gamma);

    const auto check = [&](const std::vector<Eigen::MatrixXcd> &kraus, const std::vector<size_t> &qubits,
                           const char *description) {
        const size_t second = qubits.size() > 1 ? qubits[1] : 0;
        QC::DensityMatrix<> dm(nrQubits);
        dm.setDensityMatrix(rho);
        const bool triple = qubits.size() == 3;
        if (triple)
            dm.ApplyChannel(kraus, qubits[0], second, qubits[2]);
        else
            dm.ApplyChannel(kraus, qubits[0], second);
        const auto draw = [&](QC::QubitRegister<> &reg, double sample) {
            return triple ? reg.ApplyKrausOperator(kraus, sample, qubits[0], second, qubits[2])
                          : reg.ApplyKrausOperator(kraus, sample, qubits[0], second);
        };

        std::vector<double> weights;
        std::vector<Eigen::VectorXcd> branches;
        double total = 0.;
        for (const auto &e : kraus)
        {
            branches.emplace_back(DM_EmbedOperator(e, qubits, nrQubits) * psi);
            weights.push_back(branches.back().squaredNorm());
            total += weights.back();
        }

        Eigen::MatrixXcd mixture = Eigen::MatrixXcd::Zero(rho.rows(), rho.cols());
        double before = 0.;
        for (size_t k = 0; k < kraus.size(); ++k)
        {
            // A draw in the middle of this operator's share must select it.
            QC::QubitRegister<> reg(nrQubits);
            reg.setRegisterStorage(psi);
            reg.ComputeStart();
            const size_t chosen = draw(reg, (before + 0.5 * weights[k]) / total);
            reg.ComputeEnd();
            before += weights[k];
            const Eigen::VectorXcd expected = branches[k] / std::sqrt(weights[k]);
            if (chosen != k || (reg.getRegisterStorage() - expected).norm() > 1E-12 ||
                std::abs(reg.getRegisterStorage().norm() - 1.) > 1E-12)
            {
                std::cout << description << ": operator " << k << " was not selected or applied correctly (selected "
                          << chosen << ")" << std::endl;
                return false;
            }
            mixture += weights[k] * (reg.getRegisterStorage() * reg.getRegisterStorage().adjoint());
            reg.setRegisterStorage(psi);
            reg.Compute();
            if ((reg.getRegisterStorage() - expected).norm() > 1E-12)
            {
                std::cout << description << ": recorded branch did not replay" << std::endl;
                return false;
            }
            if (qubits.size() >= 2)
            {
                reg.setRegisterStorage(psi);
                const double norm = triple ? reg.ApplyLocalMapNorm(kraus[k], qubits[0], second, qubits[2])
                                           : reg.ApplyLocalMapNorm(kraus[k], qubits[0], second);
                if (std::abs(norm - weights[k]) > 1E-12 || (reg.getRegisterStorage() - branches[k]).norm() > 1E-12)
                {
                    std::cout << description << ": fused application and norm differ from the reference" << std::endl;
                    return false;
                }
            }
        }
        if ((mixture - dm.getDensityMatrix()).norm() > 1E-12)
        {
            std::cout << description << ": the weighted branches do not rebuild the channel (error "
                      << (mixture - dm.getDensityMatrix()).norm() << ")" << std::endl;
            return false;
        }

        // The ends of the unit interval select the first and the last operator with any weight.
        QC::QubitRegister<> first(nrQubits), final(nrQubits);
        first.setRegisterStorage(psi);
        final.setRegisterStorage(psi);
        if (draw(first, 0.) != 0 || draw(final, std::nextafter(1., 0.)) != kraus.size() - 1)
        {
            std::cout << description << ": the ends of the unit interval select the wrong operators" << std::endl;
            return false;
        }
        return true;
    };

    for (size_t q = 0; q < nrQubits; ++q)
        if (!check(DM_RandomKraus(2, 3, random), {q}, "Random one-qubit Kraus set") ||
            !check({damp0, damp1}, {q}, "Amplitude damping"))
            return false;
    for (const auto &pair : std::vector<std::vector<size_t>>{{0, 1}, {1, 0}, {0, 3}, {3, 1}, {2, 3}})
        if (!check(DM_RandomKraus(4, 5, random), pair, "Random two-qubit Kraus set"))
            return false;
    if (!check(DM_RandomKraus(4, 2, random, true), {3, 1}, "Dense real two-qubit Kraus set") ||
        !check(DM_RandomKraus(8, 3, random, true), {3, 2, 1}, "Dense real three-qubit Kraus set") ||
        !check(DM_RandomKraus(8, 2, random), {2, 3, 0}, "Dense complex three-qubit Kraus set"))
        return false;
    // Opposite relative phases cause exact destructive interference in different spectator sectors.
    // The sectors lie in separate reduction blocks. The two dense maps have relative weights 4:1.
    {
        using Complex = std::complex<double>;
        const std::array<Complex, 4> phase{1., Complex(0., 1.), -1., Complex(0., -1.)};
        Eigen::MatrixXcd op(8, 8);
        for (int r = 0; r < 8; ++r)
            for (int c = 0; c < 8; ++c)
                op(r, c) = phase[(r + c) % 4] / 8.;
        const std::vector<Eigen::MatrixXcd> operators{op, op.conjugate()};
        for (const std::array<size_t, 3> qubits : {std::array<size_t, 3>{0, 12, 2}, {2, 0, 12}, {12, 2, 0}, {1, 12, 3}})
        {
            const auto offset = [&](unsigned local) {
                size_t index = 0;
                for (unsigned b = 0; b < 3; ++b)
                    if (local & (1u << b))
                        index |= size_t{1} << qubits[b];
                return index;
            };
            const size_t mask = offset(7);
            const size_t other = (size_t{1} << 15) | ((mask + 1) & ~mask);
            Eigen::VectorXcd input = Eigen::VectorXcd::Zero(1 << 16);
            input[0] = 0.25;
            input[offset(5)] = Complex(0., 0.25);
            input[other] = 0.5;
            input[other | offset(5)] = Complex(0., -0.5);
            for (size_t k = 0; k < 2; ++k)
            {
                QC::QubitRegister<> reg(16);
                reg.setRegisterStorage(input);
                reg.ComputeStart();
                const size_t selected =
                    reg.ApplyKrausOperator(operators, k == 0 ? 0.79 : 0.81, qubits[0], qubits[1], qubits[2]);
                reg.ComputeEnd();
                Eigen::VectorXcd expected = Eigen::VectorXcd::Zero(input.size());
                for (unsigned r = 0; r < 8; ++r)
                    expected[(k == 0 ? other : 0) | offset(r)] =
                        (k == 0 ? phase[r % 4] : std::conj(phase[r % 4])) / std::sqrt(8.);
                if (selected != k || !reg.getRegisterStorage().isApprox(expected, 1E-12))
                {
                    std::cout << "Dense complex Kraus interference selected or normalized the wrong branch"
                              << std::endl;
                    return false;
                }
                reg.setRegisterStorage(input);
                reg.Compute();
                if (!reg.getRegisterStorage().isApprox(expected, 1E-12))
                {
                    std::cout << "Dense complex Kraus interference did not replay" << std::endl;
                    return false;
                }
            }
        }
    }
    // Even a tiny imaginary Kraus coefficient can cancel a real one through the state's coherence.
    {
        const double scale = 1E-140;
        Eigen::MatrixXcd filter = Eigen::MatrixXcd::Zero(4, 4);
        filter(0, 0) = scale;
        filter(0, 1) = std::complex<double>(0., scale);
        Eigen::VectorXcd input = Eigen::VectorXcd::Zero(4);
        input[0] = 1.;
        input[2] = std::complex<double>(0., 1.);
        QC::QubitRegister<> reg(2);
        reg.setRegisterStorage(input);
        if (reg.ApplyKrausOperator({filter, scale * Eigen::MatrixXcd::Identity(4, 4)}, 0., 1, 0) != 1)
            return false;
    }
    {
        Eigen::MatrixXcd stay = Eigen::MatrixXcd::Identity(8, 8), jump = Eigen::MatrixXcd::Zero(8, 8);
        stay(5, 5) = 0.5;
        jump(2, 5) = std::complex<double>(0., std::sqrt(0.75));
        if (!check({stay, jump}, {1, 3, 0}, "Three-qubit population-weighted channel"))
            return false;
    }
    // The two contributing columns almost cancel. Recompute the small weight from the amplitudes
    // and replay both the 8x8 operator and its normalization on the original targets.
    {
        Eigen::MatrixXcd filter = Eigen::MatrixXcd::Zero(8, 8);
        filter(3, 0) = 0.5;
        filter(3, 7) = -0.5;
        Eigen::VectorXcd input = Eigen::VectorXcd::Zero(16);
        input[0] = 1.;
        input[13] = 1. - 1E-9;
        input.normalize();
        QC::QubitRegister<> reg(4);
        reg.setRegisterStorage(input);
        reg.ComputeStart();
        reg.SampleKrausOperator({filter}, 3, 0, 2);
        reg.ComputeEnd();
        const Eigen::VectorXcd expected = Eigen::VectorXcd::Unit(16, 9);
        if (!reg.getRegisterStorage().isApprox(expected, 1E-12))
            return false;
        reg.setRegisterStorage(input);
        reg.Compute();
        if (!reg.getRegisterStorage().isApprox(expected, 1E-12))
            return false;
    }
    // More operators than the stack buffer holds.
    std::vector<Eigen::MatrixXcd> many;
    for (const auto &e : DM_RandomKraus(2, 3, random))
        for (int copy = 0; copy < 7; ++copy)
            many.emplace_back(e / std::sqrt(7.));
    if (!check(many, {2}, "Kraus set with 21 operators"))
        return false;

    // A sparse two-qubit channel, including complex phases and a jump between different columns. Each branch
    // must preserve the rest of an entangled register and renormalize even if the state came from a filter.
    {
        Eigen::MatrixXcd stay = Eigen::MatrixXcd::Identity(4, 4), jump = Eigen::MatrixXcd::Zero(4, 4);
        stay(1, 1) = std::sqrt(0.25);
        stay(2, 2) = std::sqrt(0.75);
        jump(3, 1) = std::complex<double>(0., std::sqrt(0.75));
        jump(1, 2) = -0.5;
        if (!check({stay, jump}, {3, 0}, "Two-qubit population-weighted channel"))
            return false;
        QC::QubitRegister<> reg(nrQubits);
        Eigen::VectorXcd scaled = 3. * psi;
        reg.setRegisterStorageFastNoNormalize(scaled);
        const Eigen::VectorXcd jumped = DM_EmbedOperator(jump, {3, 0}, nrQubits) * psi;
        if (reg.ApplyKrausOperator({stay, jump}, std::nextafter(1., 0.), 3, 0) != 1 ||
            (reg.getRegisterStorage() - jumped.normalized()).norm() > 1E-12)
        {
            std::cout << "A population-weighted branch on an unnormalized register is wrong" << std::endl;
            return false;
        }
    }

    // A nonzero overlap in a row, however small, can decide the branch at a boundary.
    {
        QC::QubitRegister<> reg(1);
        reg.ApplyGate(QC::Gates::HadamardGate<>(), 0);
        Eigen::MatrixXcd projector = Eigen::MatrixXcd::Zero(2, 2);
        projector(0, 0) = 1.;
        Eigen::MatrixXcd tilted = projector;
        tilted(0, 1) = 1E-12;
        if (reg.ApplyKrausOperator({tilted, projector}, 0.5, 0) != 0)
        {
            std::cout << "A nearly diagonal Kraus effect lost its coherence" << std::endl;
            return false;
        }
    }

    // Squaring the coefficient first would overflow although the image and branch weight are finite.
    {
        QC::QubitRegister<> reg(1);
        Eigen::VectorXcd input(2);
        input << 1., 1E-150;
        reg.setRegisterStorageFastNoNormalize(input);
        Eigen::MatrixXcd gain = Eigen::MatrixXcd::Zero(2, 2);
        gain(0, 1) = 1E200;
        reg.ApplyKrausOperator({gain}, 0.5, 0);
        if (std::abs(reg.getBasisStateProbability(0) - 1.) > 1E-12)
        {
            std::cout << "A sparse Kraus operator on a small amplitude was normalized incorrectly" << std::endl;
            return false;
        }
    }

    // In a projective measurement a positive branch can be smaller than the squared-weight range.
    // Splitting that projector into four equal Kraus operators makes even its amplitude underflow.
    for (const double amplitude :
         {1.2 * std::ldexp(1., -537), std::ldexp(1., -540), std::numeric_limits<double>::denorm_min()})
    {
        Eigen::MatrixXcd projector = Eigen::MatrixXcd::Zero(8, 8);
        projector(7, 7) = 1.;
        const Eigen::MatrixXcd other = Eigen::MatrixXcd::Identity(8, 8) - projector;
        const size_t count = amplitude == std::numeric_limits<double>::denorm_min() ? 4 : 1;
        std::vector<Eigen::MatrixXcd> family(count, projector / std::sqrt(double(count)));
        family.emplace_back(other);
        Eigen::VectorXcd input = Eigen::VectorXcd::Zero(8);
        input(0) = 1.;
        input(7) = amplitude;
        QC::QubitRegister<> reg(3);
        reg.setRegisterStorage(input);
        reg.ComputeStart();
        // The first branch weighs about 1.44 denorm_min in the first case.
        const double draw = amplitude > std::ldexp(1., -537) ? std::numeric_limits<double>::denorm_min() : 0.;
        if (reg.ApplyKrausOperator(family, draw, 2, 0, 1) != 0 ||
            std::abs(reg.getBasisStateAmplitude(7) - 1.) > 1E-13 ||
            std::abs(reg.getRegisterStorage().squaredNorm() - 1.) > 1E-13)
            return false;
        reg.ComputeEnd();
        reg.setRegisterStorage(input);
        reg.Compute();
        if (std::abs(reg.getBasisStateAmplitude(7) - 1.) > 1E-13)
            return false;
    }
    // The selected matrix times its normalization factor would overflow on an unoccupied column.
    // Recording must still act on that column if the same gate sequence is replayed on another input.
    {
        Eigen::MatrixXcd filter = Eigen::MatrixXcd::Zero(8, 8);
        filter(0, 0) = 1E-150;
        filter(7, 7) = 1E200;
        QC::QubitRegister<> reg(3);
        reg.ComputeStart();
        reg.ApplyKrausOperator({filter}, 0.5, 0, 1, 2);
        reg.ComputeEnd();
        if (!reg.getRegisterStorage().allFinite() || std::abs(reg.getBasisStateAmplitude(0) - 1.) > 1E-13)
            return false;
        Eigen::VectorXcd probe = Eigen::VectorXcd::Zero(8);
        probe(7) = 1E-300;
        reg.setRegisterStorageFastNoNormalize(probe);
        reg.Compute();
        if (std::abs(reg.getBasisStateAmplitude(7) / 1E50 - 1.) > 1E-13)
            return false;
    }
    // A zero population can be an empty column or a nonzero amplitude whose square underflows.
    // Put the only such amplitude in the last reduction block, including the odd AVX group.
    for (size_t q : {size_t{0}, size_t{4}})
        for (bool dense : {false, true})
            for (bool occupied : {false, true})
            {
                constexpr size_t n = 16, q2 = 7, q3 = 13;
                const size_t spectator =
                    ((size_t{1} << n) - 1) ^ ((size_t{1} << q) | (size_t{1} << q2) | (size_t{1} << q3));
                Eigen::VectorXcd input = Eigen::VectorXcd::Zero(size_t{1} << n);
                input[0] = 1.;
                input[spectator | (size_t{1} << q) | (size_t{1} << q3)] = {-0., occupied ? -1E-200 : -0.};
                Eigen::MatrixXcd stay = Eigen::MatrixXcd::Zero(8, 8), gain = stay;
                stay(0, 0) = 1.;
                gain(6, 5) = 1E200;
                if (dense)
                    gain(6, 3) = {0.25, 0.125};
                QC::QubitRegister<> reg(n);
                reg.setRegisterStorageFastNoNormalize(input);
                const size_t chosen = reg.ApplyKrausOperator({stay, gain}, 0.75, q, q2, q3);
                const size_t expectedIndex = occupied ? spectator | (size_t{1} << q2) | (size_t{1} << q3) : 0;
                const std::complex<double> expected = occupied ? std::complex<double>(0., -1.) : 1.;
                if (chosen != size_t(occupied) ||
                    std::abs(reg.getBasisStateAmplitude(expectedIndex) - expected) > 1E-13)
                    return false;
            }
    // The draw lies exactly at the end of the first branch and inside the next, however small it is.
    {
        const Eigen::MatrixXcd halfI = 0.5 * Eigen::MatrixXcd::Identity(8, 8);
        Eigen::MatrixXcd halfP = Eigen::MatrixXcd::Zero(8, 8);
        halfP(7, 7) = 0.5;
        Eigen::VectorXcd input = Eigen::VectorXcd::Zero(8);
        input(0) = 1.;
        input(7) = std::ldexp(1., -540);
        QC::QubitRegister<> reg(3);
        reg.setRegisterStorage(input);
        if (reg.ApplyKrausOperator({halfI, halfP, halfI - halfP, halfI, halfI}, 0.25, 0, 2, 1) != 1 ||
            std::abs(reg.getBasisStateAmplitude(7) - 1.) > 1E-13)
            return false;
    }
    // An earlier tiny weight moves the midpoint strictly inside the first of two equal large branches.
    {
        std::vector<Eigen::MatrixXcd> projectors(4, Eigen::MatrixXcd::Zero(4, 4));
        for (int i = 0; i < 4; ++i)
            projectors[i](i, i) = 1.;
        Eigen::VectorXcd input(4);
        input << std::ldexp(1., -530), 1., 1., 0.;
        QC::QubitRegister<> reg(2);
        reg.setRegisterStorage(input);
        if (reg.ApplyKrausOperator(projectors, 0.5, 0, 1) != 1 || std::abs(reg.getBasisStateAmplitude(1) - 1.) > 1E-13)
            return false;
    }
    // A small common scale forces the exceptional path. Cancel first and normalize the actual image.
    {
        const double delta = std::ldexp(1., -35), scale = std::ldexp(1., -600);
        Eigen::VectorXcd input = Eigen::VectorXcd::Zero(8), expected = Eigen::VectorXcd::Zero(8);
        input(0) = 0.5;
        input(7) = 0.5 * (1. - delta);
        Eigen::MatrixXcd filter = Eigen::MatrixXcd::Zero(8, 8);
        filter(0, 0) = 1.5 * scale;
        filter(0, 7) = -1.5 * scale;
        filter(2, 0) = delta * scale;
        expected(0) = 1.5;
        expected(2) = 1.;
        expected.normalize();
        QC::QubitRegister<> reg(3);
        Eigen::VectorXcd original = input;
        reg.setRegisterStorageFastNoNormalize(original);
        reg.ComputeStart();
        reg.ApplyKrausOperator({filter}, 0.5, 0, 1, 2);
        reg.ComputeEnd();
        if (!reg.getRegisterStorage().isApprox(expected, 1E-12))
            return false;
        reg.setRegisterStorageFastNoNormalize(input);
        reg.Compute();
        if (!reg.getRegisterStorage().isApprox(expected, 1E-12))
            return false;
    }
    // An operator without weight on the state is never selected: |0> under a Z-basis measurement channel.
    {
        Eigen::MatrixXcd p0 = Eigen::MatrixXcd::Zero(2, 2), p1 = Eigen::MatrixXcd::Zero(2, 2);
        p0(0, 0) = 1.;
        p1(1, 1) = 1.;
        QC::QubitRegister<> reg(2);
        if (reg.ApplyKrausOperator({p1, p0}, 0., 0) != 1 ||
            reg.ApplyKrausOperator({p0, p1}, std::nextafter(1., 0.), 1) != 0 ||
            std::abs(reg.getBasisStateProbability(0) - 1.) > 1E-14)
        {
            std::cout << "An operator with no probability was selected" << std::endl;
            return false;
        }
        // The register's own stream gives the same kind of step.
        reg.SetSeed(11);
        reg.ApplyGate(QC::Gates::HadamardGate<>(), 0);
        size_t ones = 0;
        for (int trial = 0; trial < 4000; ++trial)
        {
            QC::QubitRegister<> copy(2);
            copy.setRegisterStorage(reg.getRegisterStorage());
            copy.SetSeed(1000 + trial);
            ones += copy.SampleKrausOperator({p0, p1}, 0);
        }
        if (ones < 1850 || ones > 2150)
        {
            std::cout << "Sampled measurement channel gave " << ones << " ones in 4000 trials" << std::endl;
            return false;
        }
    }

    // A branch whose weight is what a cancellation in the reduced density matrix left: the entries of the matrix
    // are known to their own rounding, far above such a weight, and the state must still leave with unit norm and
    // in the state the operator projects on.
    {
        Eigen::MatrixXcd minus(2, 2), bell = Eigen::MatrixXcd::Zero(4, 4), cnot = Eigen::MatrixXcd::Zero(4, 4);
        minus << 0.5, -0.5, -0.5, 0.5; // |-><-|
        bell(0, 0) = bell(3, 3) = 0.5; // the projector on (|00> - |11>) / sqrt(2)
        bell(0, 3) = bell(3, 0) = -0.5;
        cnot(0, 0) = cnot(1, 1) = cnot(2, 3) = cnot(3, 2) = 1.; // the control is the high bit
        for (const size_t qubits : {size_t{2}, size_t{3}, size_t{10}})
            for (const double delta : {1E-3, 1E-5, 1E-6})
                for (const bool pair : {false, true})
                {
                    // Qubit 0 is |+> turned by delta, so <-|psi> = -sin(delta); with qubit 1 flipped along, the
                    // pair is as far from the projector's state.
                    QC::QubitRegister<> reg(qubits);
                    const double half = std::acos(-1.) / 4. + delta;
                    Eigen::MatrixXcd turn(2, 2);
                    turn << std::cos(half), -std::sin(half), std::sin(half), std::cos(half);
                    reg.ApplyGate(QC::Gates::AppliedGate<>(turn, 0));
                    if (pair)
                        reg.ApplyGate(QC::Gates::AppliedGate<>(cnot, 1, 0));
                    for (size_t q = pair ? 2 : 1; q < qubits; ++q)
                        reg.ApplyGate(QC::Gates::HadamardGate<>(), q);
                    if (pair)
                        reg.ApplyKrausOperator({bell}, 0.5, 0, 1);
                    else
                        reg.ApplyKrausOperator({minus}, 0.5, 0);
                    const Eigen::VectorXcd state = reg.getRegisterStorage();
                    // Every amplitude with qubit 0 set is minus its partner's (qubit 1 set as well, for the pair).
                    double mismatch = 0.;
                    const size_t flip = pair ? 3 : 1;
                    for (size_t i = 0; i < static_cast<size_t>(state.size()); ++i)
                        if ((i & flip) == 0)
                            mismatch = std::max(mismatch, std::abs(state[i] + state[i | flip]));
                    if (std::abs(state.squaredNorm() - 1.) > 1E-10 || mismatch > 1E-10)
                    {
                        std::cout << "A Kraus operator of weight " << delta * delta << " on " << qubits
                                  << " qubits left squared norm " << state.squaredNorm() << " and a mismatch of "
                                  << mismatch << std::endl;
                        return false;
                    }
                }
    }

    // Weights that the reduced density matrix cannot hold at all (its rounding is above them, and may leave them
    // negative) are read from the amplitudes before anything is drawn: an operator that has a weight is applied, and
    // one that has none is refused with the register as it was. A recorded step replays to the same state.
    {
        Eigen::MatrixXcd minus(2, 2);
        minus << 0.5, -0.5, -0.5, 0.5; // |-><-|
        const auto turned = [](double delta) {
            const double half = std::acos(-1.) / 4. + delta; // |+> turned by delta: <-|psi> = -sin(delta)
            Eigen::MatrixXcd turn(2, 2);
            turn << std::cos(half), -std::sin(half), std::sin(half), std::cos(half);
            return turn;
        };
        for (const size_t qubits : {size_t{1}, size_t{4}, size_t{10}})
            for (const double delta : {1E-8, 3E-9, 1E-9, 1E-10, 1E-12})
            {
                QC::QubitRegister<> reg(qubits);
                reg.ApplyGate(QC::Gates::AppliedGate<>(turned(delta), 0));
                for (size_t q = 1; q < qubits; ++q)
                    reg.ApplyGate(QC::Gates::HadamardGate<>(), q);
                const Eigen::VectorXcd before = reg.getRegisterStorage();
                reg.ComputeStart();
                try
                {
                    reg.ApplyKrausOperator({minus}, 0.5, 0);
                }
                catch (const std::exception &error)
                {
                    std::cout << "A Kraus operator of weight " << delta * delta << " on " << qubits
                              << " qubits was refused: " << error.what() << std::endl;
                    return false;
                }
                reg.ComputeEnd();
                const Eigen::VectorXcd direct = reg.getRegisterStorage();
                double mismatch = 0.;
                for (size_t i = 0; i < static_cast<size_t>(direct.size()); i += 2)
                    mismatch = std::max(mismatch, std::abs(direct[i] + direct[i | 1]));
                Eigen::VectorXcd restored = before;
                reg.setRegisterStorageFastNoNormalize(restored);
                reg.Compute();
                const double replayed = (reg.getRegisterStorage() - direct).norm();
                if (std::abs(direct.squaredNorm() - 1.) > 1E-12 || mismatch > 1E-12 || replayed != 0.)
                {
                    std::cout << "A Kraus operator of weight " << delta * delta << " on " << qubits
                              << " qubits left squared norm " << direct.squaredNorm() << ", a mismatch of " << mismatch
                              << " and a replay off by " << replayed << std::endl;
                    return false;
                }
            }
        // Qubit 0 exactly in |+>: the operator has no weight.
        std::mt19937_64 angles(3);
        std::uniform_real_distribution<double> angle(0., 2. * std::acos(-1.));
        for (int trial = 0; trial < 200; ++trial)
        {
            const size_t qubits = 2 + trial % 9;
            QC::QubitRegister<> reg(qubits);
            for (size_t q = 1; q < qubits; ++q)
                reg.ApplyGate(QC::Gates::AppliedGate<>(turned(angle(angles)), q));
            reg.ApplyGate(QC::Gates::HadamardGate<>(), 0);
            const Eigen::VectorXcd before = reg.getRegisterStorage();
            if (!DM_ExpectDomainError([&] { reg.ApplyKrausOperator({minus}, 0.5, 0); },
                                      "Kraus operator with no weight on the state"))
                return false;
            if ((reg.getRegisterStorage() - before).norm() != 0.)
            {
                std::cout << "A refused Kraus operator changed the register" << std::endl;
                return false;
            }
            // In a set, it is never the one drawn, whatever the rounding of the density matrix made of its weight.
            Eigen::MatrixXcd plus(2, 2);
            plus << 0.5, 0.5, 0.5, 0.5;
            if (reg.ApplyKrausOperator({minus, plus}, 0., 0) != 1 ||
                reg.ApplyKrausOperator({plus, minus}, std::nextafter(1., 0.), 0) != 0)
            {
                std::cout << "A Kraus operator with no weight on the state was drawn" << std::endl;
                return false;
            }
        }
    }

    // A weight read from the amplitudes that falls below its cutoff counts as none, although it is positive: the
    // message then says what was found and what it was compared with. The state (1/2, 1/2, 1/2, 1/2) and the row
    // (1/2, -1/2 + 2^-49) leave 2^-50 on each pair, exactly, which is at eight units of rounding.
    {
        QC::QubitRegister<> exact(2);
        const Eigen::VectorXcd half = Eigen::VectorXcd::Constant(4, 0.5);
        exact.setRegisterStorage(half);
        Eigen::MatrixXcd faint(2, 2);
        faint << 0.5, -0.5 + std::ldexp(1., -49), 0., 0.;
        std::string message;
        try
        {
            exact.ApplyKrausOperator({faint}, 0.5, 0);
        }
        catch (const std::domain_error &error)
        {
            message = error.what();
        }
        if (message.find("largest weight found") == std::string::npos ||
            (exact.getRegisterStorage() - half).norm() != 0.)
        {
            std::cout << "A Kraus weight below its cutoff was not reported as found, or the register changed: \""
                      << message << "\"" << std::endl;
            return false;
        }
    }

    QC::QubitRegister<> reg(3);
    const Eigen::MatrixXcd i2 = Eigen::MatrixXcd::Identity(2, 2), i4 = Eigen::MatrixXcd::Identity(4, 4);
    if (!DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({}, 0.5, 0); }, "Empty Kraus set") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({Eigen::MatrixXcd::Identity(8, 8)}, 0.5, 0); },
                                  "Missing targets for an 8x8 Kraus operator") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({Eigen::MatrixXcd::Identity(8, 8)}, 0.5, 1, 0, 1); },
                                  "Duplicate three-qubit Kraus targets") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({Eigen::MatrixXcd::Identity(8, 8)}, 0.5, 1, 0, 3); },
                                  "Out-of-range third Kraus target") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({i4}, 0.5, 0, 1, 2); }, "4x4 Kraus on three targets") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({i2, i4}, 0.5, 0); }, "Mixed Kraus dimensions") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({i2}, 0.5, 3); }, "Out-of-range Kraus qubit") ||
        !DM_ExpectInvalidArgument([&] { reg.ApplyKrausOperator({i4}, 0.5, 1, 1); }, "Duplicate Kraus qubits") ||
        !DM_ExpectInvalidArgument([&] { reg.LocalDensityMatrix(0, 1, 3); }, "Out-of-range density-matrix qubit") ||
        !DM_ExpectDomainError([&] { reg.ApplyKrausOperator({Eigen::MatrixXcd::Zero(2, 2)}, 0.5, 0); },
                              "Kraus set without weight"))
        return false;

    std::cout << "Success" << std::endl;
    return true;
}

static bool DensityMatrixInitializationAndValidationTest()
{
    std::cout << "\nDensity matrix simulator - initialization and invalid input handling" << std::endl;

    if (!DM_ExpectInvalidArgument([] { QC::DensityMatrix<> invalid(0); }, "Zero-qubit construction") ||
        !DM_ExpectInvalidArgument([] { QC::DensityMatrix<> invalid(std::numeric_limits<size_t>::digits); },
                                  "Oversized-qubit construction"))
        return false;
    try
    {
        QC::DensityMatrix<> overflow(std::numeric_limits<size_t>::digits / 2);
        std::cout << "Overflowing density-matrix storage construction did not throw" << std::endl;
        return false;
    }
    catch (const std::length_error &)
    {
    }

    QC::DensityMatrix<> dm(2);
    if (!dm.GetMultithreading())
    {
        std::cout << "Default multithreading state should be true" << std::endl;
        return false;
    }
    dm.SetMultithreading(false);
    if (dm.GetMultithreading())
    {
        std::cout << "SetMultithreading(false) failed" << std::endl;
        return false;
    }
    dm.SetMultithreading(true);
    if (!dm.GetMultithreading())
    {
        std::cout << "SetMultithreading(true) failed" << std::endl;
        return false;
    }

    QC::Gates::PauliXGate<> x;
    QC::Gates::CNOTGate<> cnot;
    QC::Gates::ToffoliGate<> toffoli;
    const Eigen::MatrixXcd original = dm.getDensityMatrix();

    Eigen::MatrixXcd malformed = Eigen::MatrixXcd::Identity(3, 3);
    QC::Gates::AppliedGate<> malformedGate(malformed, 0);
    Eigen::MatrixXcd tooLarge = Eigen::MatrixXcd::Identity(16, 16);
    QC::Gates::AppliedGate<> tooLargeGate(tooLarge, 0, 1);
    Eigen::MatrixXcd nonFiniteGateMatrix = Eigen::MatrixXcd::Identity(2, 2);
    nonFiniteGateMatrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
    QC::Gates::AppliedGate<> nonFiniteGate(nonFiniteGateMatrix, 0);

    if (!DM_ExpectInvalidArgument([&] { dm.ApplyGate(x, 2); }, "Out-of-range gate qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyGate(cnot, 0, 0); }, "Duplicate two-qubit gate arguments") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyGate(toffoli, 0, 1, 1); }, "Duplicate three-qubit gate arguments") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyGate(malformedGate); }, "Non-power-of-two gate matrix") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyGate(tooLargeGate); }, "Unsupported four-qubit gate") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyGate(nonFiniteGate); }, "Non-finite gate matrix") ||
        !DM_ExpectInvalidArgument([&] { dm.GetQubitProbability(2); }, "Out-of-range probability qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.MeasureQubit(2); }, "Out-of-range measured qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.DephaseMeasure(2); }, "Out-of-range dephased qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.MeasureNoCollapse(std::set<size_t>{0, 2}); }, "Out-of-range subset qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.RepeatedMeasure(1, 0, 1); }, "Reversed measurement range") ||
        !DM_ExpectInvalidArgument([&] { dm.RepeatedMeasureUnordered(0, 2, 1); }, "Out-of-range measurement range"))
        return false;

    if (!DM_ExpectInvalidArgument([&] { dm.ApplyBitFlipNoise(0, -0.1); }, "Negative noise probability") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyAmplitudeDamping(0, 1.1); }, "Noise probability above one") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyPhaseDamping(0, std::numeric_limits<double>::quiet_NaN()); },
                                  "NaN noise probability"))
        return false;

    const Eigen::MatrixXcd i2 = Eigen::MatrixXcd::Identity(2, 2);
    const Eigen::MatrixXcd i4 = Eigen::MatrixXcd::Identity(4, 4);
    const Eigen::MatrixXcd i8 = Eigen::MatrixXcd::Identity(8, 8);
    Eigen::MatrixXcd nonFiniteKraus = i2;
    nonFiniteKraus(1, 1) = std::numeric_limits<double>::infinity();
    if (!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({}, 0); }, "Empty Kraus channel") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({0.5 * i2}, 0); }, "Non trace-preserving Kraus channel") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({i2, i4}, 0); }, "Mismatched Kraus dimensions") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({i4}, 0, 0); }, "Duplicate two-qubit channel arguments") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({i2}, 2); }, "Out-of-range channel qubit") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({i8}, 0); }, "Missing targets for an 8x8 channel",
                                  "three target qubits") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({Eigen::MatrixXcd::Identity(3, 3)}, 0); },
                                  "3x3 channel operator", "dimensions") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({Eigen::MatrixXcd::Zero(4, 3)}, 0); }, "4x3 channel operator",
                                  "dimensions") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({i8}, 0, 1, 1); }, "Duplicate three-qubit channel targets") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({i4}, 0, 1, 0); }, "4x4 channel on three targets") ||
        !DM_ExpectInvalidArgument([&] { dm.ApplyChannel({nonFiniteKraus}, 0); }, "Non-finite Kraus operator"))
        return false;

    Eigen::VectorXcd wrongVector = Eigen::VectorXcd::Zero(3);
    Eigen::VectorXcd zeroVector = Eigen::VectorXcd::Zero(4);
    Eigen::VectorXcd nonFiniteVector = Eigen::VectorXcd::Ones(4);
    nonFiniteVector(2) = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXcd wrongMatrix = Eigen::MatrixXcd::Identity(2, 2);
    if (!DM_ExpectInvalidArgument([&] { dm.setFromStatevector(wrongVector); }, "Wrong statevector dimension") ||
        !DM_ExpectInvalidArgument([&] { dm.setFromStatevector(zeroVector); }, "Zero statevector") ||
        !DM_ExpectInvalidArgument([&] { dm.setFromStatevector(nonFiniteVector); }, "Non-finite statevector") ||
        !DM_ExpectInvalidArgument([&] { dm.setDensityMatrix(wrongMatrix); }, "Wrong density-matrix dimension") ||
        !DM_ExpectInvalidArgument([&] { dm.ExpectationValue(wrongMatrix); }, "Wrong observable dimension") ||
        !DM_ExpectInvalidArgument([&] { dm.ExpectationValue(std::string("IA")); }, "Invalid Pauli operator") ||
        !DM_ExpectInvalidArgument([&] { dm.setToBasisState(4); }, "Out-of-range basis state") ||
        !DM_ExpectInvalidArgument([&] { dm.setToMixtureOfBasisStates({}); }, "Empty basis-state mixture") ||
        !DM_ExpectInvalidArgument([&] { dm.setToMixtureOfBasisStates({{0, 0.}, {7, 1.}}); },
                                  "Mixture without a valid positive weight") ||
        !DM_ExpectInvalidArgument([&] { dm.setToMixtureOfBasisStates({{0, std::numeric_limits<double>::infinity()}}); },
                                  "Non-finite mixture weight"))
        return false;

    if ((dm.getDensityMatrix() - original).norm() > 1E-12)
    {
        std::cout << "An invalid operation changed the density matrix" << std::endl;
        return false;
    }
    dm.ApplyChannel({i4}, 1); // The omitted second target defaults to qubit zero.
    if (dm.getDensityMatrix() != original)
        return false;

    // Non-normalized finite vectors are accepted and normalized into a unit-trace pure state.
    Eigen::VectorXcd psi(4);
    psi << 1., std::complex<double>(0., 2.), -1., 0.5;
    dm.setFromStatevector(psi);
    if (!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-12) || !approxEqual(dm.Purity(), 1., 1E-12))
    {
        std::cout << "setFromStatevector did not normalize a finite non-zero vector" << std::endl;
        return false;
    }

    // Duplicate valid states accumulate; non-positive and out-of-range entries are ignored.
    dm.setToMixtureOfBasisStates({{1, 1.}, {1, 2.}, {3, 1.}, {0, -4.}, {7, 9.}});
    if (!approxEqual(dm.getBasisStateProbability(1), 0.75, 1E-12) ||
        !approxEqual(dm.getBasisStateProbability(3), 0.25, 1E-12) ||
        !approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-12) || !approxEqual(dm.Purity(), 0.625, 1E-12))
    {
        std::cout << "Basis-state mixture normalization or duplicate accumulation is wrong" << std::endl;
        return false;
    }

    // All closed-interval endpoints are valid for every predefined noise channel.
    for (const double probability : {0., 1.})
    {
        QC::DensityMatrix<> boundary(1);
        boundary.ApplyGate(QC::Gates::HadamardGate<>(), 0);
        boundary.ApplyBitFlipNoise(0, probability);
        boundary.ApplyPhaseFlipNoise(0, probability);
        boundary.ApplyDepolarizingNoise(0, probability);
        boundary.ApplyAmplitudeDamping(0, probability);
        boundary.ApplyPhaseDamping(0, probability);
        if (!boundary.getDensityMatrix().allFinite() || !boundary.IsHermitian(1E-10) ||
            !approxEqual(boundary.Trace(), std::complex<double>(1., 0.), 1E-10))
        {
            std::cout << "A noise channel failed at probability endpoint " << probability << std::endl;
            return false;
        }
    }

    QC::DensityMatrix<> cleared(1);
    cleared.Clear();
    if (!DM_ExpectDomainError([&] { cleared.MeasureNoCollapse(); }, "Sampling an empty state"))
        return false;

    QC::DensityMatrix<> negativePopulation(1);
    Eigen::MatrixXcd invalidRho = Eigen::MatrixXcd::Zero(2, 2);
    invalidRho(0, 0) = 1.1;
    invalidRho(1, 1) = -0.1;
    negativePopulation.setDensityMatrix(invalidRho);
    if (!DM_ExpectDomainError([&] { negativePopulation.RepeatedMeasure(2); }, "Sampling a negative population"))
        return false;

    QC::DensityMatrix<> complexPopulation(1);
    invalidRho.setZero();
    invalidRho(0, 0) = std::complex<double>(1., 0.01);
    complexPopulation.setDensityMatrix(invalidRho);
    if (!DM_ExpectDomainError([&] { complexPopulation.MeasureQubit(0); }, "Measuring a non-real population"))
        return false;

    std::cout << "Success" << std::endl;
    return true;
}

// checks that single qubit measurement probabilities match the statevector simulator
static bool DensityMatrixMeasurementTest()
{
    std::cout << "\nDensity matrix simulator - measurement probabilities compared against the statevector simulator"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(15, 30);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        for (int t = 0; t < 10; ++t)
        {
            QC::DensityMatrix<> dm(nrQubits);
            QC::QubitRegister<> reg(nrQubits);

            const int lim = nrGatesDistr(gen);
            for (int i = 0; i < lim; ++i)
            {
                const int g = gateDistr(gen);
                const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
                if (twoQubitsGate && nrQubits < 2)
                    continue;

                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = qubit1 + 1;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

                dm.ApplyGate(*gates[g], qubit1, qubit2);
                reg.ApplyGate(*gates[g], qubit1, qubit2);
            }

            for (int q = 0; q < nrQubits; ++q)
                if (!approxEqual(dm.GetQubitProbability(q), reg.GetQubitProbability(q), 1E-9))
                {
                    std::cout << "Qubit " << q << " probability mismatch for " << nrQubits
                              << " qubits: " << reg.GetQubitProbability(q) << " (statevector) vs "
                              << dm.GetQubitProbability(q) << " (density matrix)" << std::endl;
                    return false;
                }
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

// reference Pauli string expectation value <psi|P|psi> computed directly from the statevector,
// P = (x) P_i where character i acts on qubit i (qubit 0 is the least significant bit)
static std::complex<double> StatevectorPauliExpectation(const QC::QubitRegister<> &reg, const std::string &pauliString)
{
    const size_t nrQubits = reg.getNrQubits();
    const size_t nrBasisStates = reg.getNrBasisStates();

    size_t flipMask = 0;
    size_t signMask = 0;
    size_t yCount = 0;
    for (size_t i = 0; i < nrQubits; ++i)
    {
        const size_t bit = 1ULL << i;
        switch (toupper(static_cast<unsigned char>(pauliString[i])))
        {
        case 'I':
            break;
        case 'X':
            flipMask |= bit;
            break;
        case 'Y':
            flipMask |= bit;
            signMask |= bit;
            ++yCount;
            break;
        case 'Z':
            signMask |= bit;
            break;
        }
    }

    std::complex<double> res = 0.;
    for (size_t k = 0; k < nrBasisStates; ++k)
    {
        size_t s = signMask & k;
        bool negative = false;
        while (s)
        {
            negative = !negative;
            s &= s - 1;
        }

        // <psi|P|psi> = sum_k conj(psi_{k^flip}) * phase(k) * psi_k
        const std::complex<double> term =
            std::conj(reg.getBasisStateAmplitude(k ^ flipMask)) * reg.getBasisStateAmplitude(k);
        res += negative ? -term : term;
    }

    static const std::complex<double> iPow[4] = {{1., 0.}, {0., 1.}, {-1., 0.}, {0., -1.}};
    res *= iPow[yCount & 3];

    return res;
}

// builds the full 2^N x 2^N Pauli string matrix (Kronecker product, qubit 0 is the least significant bit)
static Eigen::MatrixXcd DM_BuildPauliMatrix(const std::string &pauliString)
{
    Eigen::Matrix2cd I;
    I << 1., 0., 0., 1.;
    Eigen::Matrix2cd X;
    X << 0., 1., 1., 0.;
    Eigen::Matrix2cd Y;
    Y << 0., std::complex<double>(0., -1.), std::complex<double>(0., 1.), 0.;
    Eigen::Matrix2cd Z;
    Z << 1., 0., 0., -1.;

    auto single = [&](char c) -> Eigen::Matrix2cd {
        switch (toupper(static_cast<unsigned char>(c)))
        {
        case 'X':
            return X;
        case 'Y':
            return Y;
        case 'Z':
            return Z;
        default:
            return I;
        }
    };

    // qubit 0 is the least significant bit, so it must be the rightmost factor of the Kronecker product
    Eigen::MatrixXcd result = single(pauliString[0]);
    for (size_t i = 1; i < pauliString.size(); ++i)
        result = Eigen::kroneckerProduct(single(pauliString[i]), result).eval();

    return result;
}

static bool DensityMatrixPauliExpectationTest()
{
    std::cout
        << "\nDensity matrix simulator - Pauli string expectation values compared against the statevector simulator"
        << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(15, 30);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    const char paulis[4] = {'I', 'X', 'Y', 'Z'};
    std::uniform_int_distribution pauliDistr(0, 3);

    for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        for (int t = 0; t < 10; ++t)
        {
            QC::DensityMatrix<> dm(nrQubits);
            QC::QubitRegister<> reg(nrQubits);

            const int lim = nrGatesDistr(gen);
            for (int i = 0; i < lim; ++i)
            {
                const int g = gateDistr(gen);
                const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
                if (twoQubitsGate && nrQubits < 2)
                    continue;

                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = qubit1 + 1;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

                dm.ApplyGate(*gates[g], qubit1, qubit2);
                reg.ApplyGate(*gates[g], qubit1, qubit2);
            }

            // try a few random Pauli strings
            for (int s = 0; s < 5; ++s)
            {
                std::string pauliString(nrQubits, 'I');
                for (int q = 0; q < nrQubits; ++q)
                    pauliString[q] = paulis[pauliDistr(gen)];

                const std::complex<double> expected = StatevectorPauliExpectation(reg, pauliString);
                const std::complex<double> got = dm.ExpectationValue(pauliString);

                // also cross check against the full matrix operator overload
                const std::complex<double> viaMatrix = dm.ExpectationValue(DM_BuildPauliMatrix(pauliString));

                if (!approxEqual(got, expected, 1E-9) || !approxEqual(got, viaMatrix, 1E-9))
                {
                    std::cout << "Pauli string " << pauliString << " expectation mismatch for " << nrQubits
                              << " qubits: " << expected << " (statevector) vs " << got << " (density matrix) vs "
                              << viaMatrix << " (full matrix)" << std::endl;
                    return false;
                }
            }
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

// checks the sampling statistics of MeasureQubit match the analytic GetQubitProbability,
// and that DephaseMeasure only kills the coherences between the two measurement sectors
static bool DensityMatrixMeasurementStatisticsTest()
{
    std::cout << "\nDensity matrix simulator - measurement sampling statistics and non selective measurement"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(15, 25);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    const int nrMeasurements = 20000;

    for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        // build a fixed random circuit so both the reference probabilities and the sampled ones use the same state
        std::vector<QC::Gates::AppliedGate<>> circuit;
        const int lim = nrGatesDistr(gen);
        for (int i = 0; i < lim; ++i)
        {
            const int g = gateDistr(gen);
            const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
            if (twoQubitsGate && nrQubits < 2)
                continue;

            int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
            int qubit2 = qubit1 + 1;
            if (twoQubitsGate && dist_bool(gen))
                std::swap(qubit1, qubit2);

            circuit.emplace_back(gates[g]->getRawOperatorMatrix(), qubit1, twoQubitsGate ? qubit2 : 0);
        }

        // analytic single qubit probabilities of the prepared state
        QC::DensityMatrix<> ref(nrQubits);
        for (const auto &gate : circuit)
            ref.ApplyGate(gate);

        std::vector<double> expectedP1(nrQubits);
        for (int q = 0; q < nrQubits; ++q)
            expectedP1[q] = ref.GetQubitProbability(q);

        // sample each qubit independently by re-preparing and measuring just that qubit
        for (int q = 0; q < nrQubits; ++q)
        {
            int ones = 0;
            for (int m = 0; m < nrMeasurements; ++m)
            {
                QC::DensityMatrix<> dm(nrQubits);
                dm.SetRandomSeed(0xD3A51000ULL +
                                 static_cast<uint64_t>(nrQubits * nrMeasurements + q * nrMeasurements + m));
                for (const auto &gate : circuit)
                    dm.ApplyGate(gate);

                if (dm.MeasureQubit(q) == 1)
                    ++ones;
            }

            const double sampled = static_cast<double>(ones) / nrMeasurements;
            if (std::abs(sampled - expectedP1[q]) > 0.05)
            {
                std::cout << "Measurement statistics mismatch for qubit " << q << " with " << nrQubits
                          << " qubits: " << expectedP1[q] << " (analytic) vs " << sampled << " (sampled)" << std::endl;
                return false;
            }
        }
    }

    // DephaseMeasure keeps the populations but removes the coherences between the measured sectors
    {
        QC::DensityMatrix<> dm(1);
        QC::Gates::HadamardGate<> h;
        dm.ApplyGate(h, 0); // |+>: rho = 1/2 [[1,1],[1,1]]

        const double p1Before = dm.GetQubitProbability(0);
        dm.DephaseMeasure(0);
        const auto &rho = dm.getDensityMatrix();

        if (!approxEqual(dm.GetQubitProbability(0), p1Before, 1E-9))
        {
            std::cout << "DephaseMeasure changed the population: " << p1Before << " -> " << dm.GetQubitProbability(0)
                      << std::endl;
            return false;
        }

        if (!approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
            !approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
        {
            std::cout << "DephaseMeasure did not remove the off diagonal coherences:\n" << rho << std::endl;
            return false;
        }

        if (!approxEqual(rho(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho(1, 1), std::complex<double>(0.5, 0), 1E-9))
        {
            std::cout << "DephaseMeasure changed the diagonal populations:\n" << rho << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

// Verifies the state update, not just measurement frequencies: rho -> P_m rho P_m / p_m.
static bool DensityMatrixCollapseTest()
{
    std::cout << "\nDensity matrix simulator - projective measurement collapse" << std::endl;

    auto checkProjection = [](QC::DensityMatrix<> &dm, size_t measuredQubit, uint64_t seed, const char *description,
                              size_t &result) -> bool {
        const Eigen::MatrixXcd before = dm.getDensityMatrix();
        dm.SetRandomSeed(seed);
        result = dm.MeasureQubit(measuredQubit);
        const size_t mask = 1ULL << measuredQubit;

        double probability = 0.;
        for (size_t state = 0; state < dm.getNrBasisStates(); ++state)
            if (((state & mask) != 0) == (result != 0))
                probability += before(state, state).real();

        Eigen::MatrixXcd expected = before;
        for (size_t row = 0; row < dm.getNrBasisStates(); ++row)
            for (size_t col = 0; col < dm.getNrBasisStates(); ++col)
                if (((row & mask) != 0) != (result != 0) || ((col & mask) != 0) != (result != 0))
                    expected(row, col) = 0.;
                else
                    expected(row, col) /= probability;

        if ((dm.getDensityMatrix() - expected).norm() > 1E-10 ||
            !approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-10) || !dm.IsHermitian(1E-10))
        {
            std::cout << description << " collapse does not equal P_m rho P_m / p_m" << std::endl;
            return false;
        }
        return true;
    };

    QC::Gates::HadamardGate<> h;
    QC::Gates::RyGate<> ry(M_PI / 5.);
    QC::Gates::CNOTGate<> cnot;
    size_t ignoredResult = 0;

    // Pure separable state.
    QC::DensityMatrix<> pure(2);
    pure.ApplyGate(h, 0);
    pure.ApplyGate(ry, 1);
    if (!checkProjection(pure, 1, 0x5011DULL, "Pure-state", ignoredResult))
        return false;

    // Mixed state containing both classical populations and coherences.
    QC::DensityMatrix<> mixed(3);
    mixed.ApplyGate(h, 0);
    mixed.ApplyGate(cnot, 2, 0);
    mixed.ApplyGate(h, 1);
    mixed.ApplyAmplitudeDamping(2, 0.23);
    mixed.ApplyPhaseDamping(1, 0.31);
    if (approxEqual(mixed.Purity(), 1., 1E-10) ||
        !checkProjection(mixed, 1, 0xC011A95EULL, "Mixed-state", ignoredResult))
        return false;

    // Entangled Bell state: verify both the projected matrix and the subsequent correlation.
    QC::DensityMatrix<> bell(2);
    bell.ApplyGate(h, 0);
    bell.ApplyGate(cnot, 1, 0);
    size_t q0 = 0;
    if (!checkProjection(bell, 0, 0xBE11ULL, "Entangled-state", q0))
        return false;
    const size_t q1 = bell.MeasureQubit(1);
    if (q0 != q1)
    {
        std::cout << "Sequential measurements broke an entangled-qubit correlation" << std::endl;
        return false;
    }

    std::cout << "Success" << std::endl;
    return true;
}

// checks the general Hermitian-observable ExpectationValue overload against a direct Tr(rho O)
// computation, verifies a few analytic Pauli expectation edge cases, that setFromStatevector and
// setDensityMatrix reproduce the reference matrix element, and that ExpectationValue throws on a
// malformed Pauli string
static bool DensityMatrixExpectationExtrasTest()
{
    std::cout << "\nDensity matrix simulator - general observable expectation, edge cases and error handling"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(15, 30);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        for (int t = 0; t < 10; ++t)
        {
            QC::DensityMatrix<> dm(nrQubits);
            QC::QubitRegister<> reg(nrQubits);

            const int lim = nrGatesDistr(gen);
            for (int i = 0; i < lim; ++i)
            {
                const int g = gateDistr(gen);
                const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
                if (twoQubitsGate && nrQubits < 2)
                    continue;

                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = qubit1 + 1;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

                dm.ApplyGate(*gates[g], qubit1, qubit2);
                reg.ApplyGate(*gates[g], qubit1, qubit2);
            }

            // build a random Hermitian observable O = M + M^dagger and compare the fast overload with Tr(rho O)
            const Eigen::Index dim = static_cast<Eigen::Index>(dm.getNrBasisStates());
            Eigen::MatrixXcd M = Eigen::MatrixXcd::Random(dim, dim);
            const Eigen::MatrixXcd O = M + M.adjoint();

            const std::complex<double> got = dm.ExpectationValue(O);
            std::complex<double> expected = 0.;
            for (Eigen::Index row = 0; row < dim; ++row)
                for (Eigen::Index col = 0; col < dim; ++col)
                    expected += dm.getDensityMatrix()(row, col) * O(col, row);

            if (!approxEqual(got, expected, 1E-9))
            {
                std::cout << "General observable expectation mismatch for " << nrQubits << " qubits: " << expected
                          << " (Tr(rho O)) vs " << got << " (ExpectationValue)" << std::endl;
                return false;
            }

            // setFromStatevector must reproduce |psi><psi| of the statevector simulator
            QC::DensityMatrix<> dmFromPsi(nrQubits);
            dmFromPsi.setFromStatevector(reg.getRegisterStorage());
            if ((dmFromPsi.getDensityMatrix() - dm.getDensityMatrix()).norm() > 1E-9)
            {
                std::cout << "setFromStatevector does not reproduce the evolved density matrix for " << nrQubits
                          << " qubits" << std::endl;
                return false;
            }

            // setDensityMatrix round trips the state
            QC::DensityMatrix<> dmCopy(nrQubits);
            dmCopy.setDensityMatrix(dm.getDensityMatrix());
            if ((dmCopy.getDensityMatrix() - dm.getDensityMatrix()).norm() > 1E-12)
            {
                std::cout << "setDensityMatrix does not round trip the state for " << nrQubits << " qubits"
                          << std::endl;
                return false;
            }

            // the all identity Pauli string expectation must equal the trace (1)
            const std::string allI(nrQubits, 'I');
            if (!approxEqual(dm.ExpectationValue(allI), std::complex<double>(1., 0.), 1E-9))
            {
                std::cout << "Identity Pauli string expectation is not the trace for " << nrQubits
                          << " qubits: " << dm.ExpectationValue(allI) << std::endl;
                return false;
            }
        }
    }

    // analytic single qubit edge cases
    {
        // <+|X|+> = 1, <+|Z|+> = 0
        QC::DensityMatrix<> dm(1);
        QC::Gates::HadamardGate<> h;
        dm.ApplyGate(h, 0);
        if (!approxEqual(dm.ExpectationValue(std::string("X")), std::complex<double>(1., 0.), 1E-9) ||
            !approxEqual(dm.ExpectationValue(std::string("Z")), std::complex<double>(0., 0.), 1E-9))
        {
            std::cout << "Single qubit |+> Pauli expectations are wrong: <X>=" << dm.ExpectationValue(std::string("X"))
                      << " <Z>=" << dm.ExpectationValue(std::string("Z")) << std::endl;
            return false;
        }
    }

    {
        // <1|Z|1> = -1
        QC::DensityMatrix<> dm(1);
        QC::Gates::PauliXGate<> x;
        dm.ApplyGate(x, 0);
        if (!approxEqual(dm.ExpectationValue(std::string("Z")), std::complex<double>(-1., 0.), 1E-9))
        {
            std::cout << "Single qubit |1> <Z> is not -1: " << dm.ExpectationValue(std::string("Z")) << std::endl;
            return false;
        }
    }

    // error handling: a Pauli string of the wrong length must throw
    {
        QC::DensityMatrix<> dm(2);
        bool threw = false;
        try
        {
            dm.ExpectationValue(std::string("X")); // only 1 char for 2 qubits
        }
        catch (const std::invalid_argument &)
        {
            threw = true;
        }

        if (!threw)
        {
            std::cout << "ExpectationValue did not throw on a Pauli string of the wrong length" << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

// samples full-register and contiguous-subregister distributions using all four RepeatedMeasure
// overloads and compares the empirical histograms against the exact diagonal populations
static bool DensityMatrixSamplingTest()
{
    std::cout << "\nDensity matrix simulator - repeated sampling against the diagonal populations" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(15, 25);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    const int nrSamples = 100000;

    for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        // prepare a fixed random state (unitary circuit plus some noise so it is mixed)
        QC::DensityMatrix<> dm(nrQubits);
        const int lim = nrGatesDistr(gen);
        for (int i = 0; i < lim; ++i)
        {
            const int g = gateDistr(gen);
            const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
            if (twoQubitsGate && nrQubits < 2)
                continue;

            int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
            int qubit2 = qubit1 + 1;
            if (twoQubitsGate && dist_bool(gen))
                std::swap(qubit1, qubit2);

            dm.ApplyGate(*gates[g], qubit1, qubit2);
        }
        dm.ApplyDepolarizingNoise(qubitDistr(gen), 0.3);
        dm.SetRandomSeed(0x5A6D1100ULL + static_cast<uint64_t>(nrQubits));

        // exact populations and a snapshot used to verify sampling does not change the state
        const size_t nrBasisStates = dm.getNrBasisStates();
        std::vector<double> expected(nrBasisStates);
        for (size_t s = 0; s < nrBasisStates; ++s)
            expected[s] = dm.getBasisStateProbability(s);
        const Eigen::MatrixXcd stateBeforeSampling = dm.getDensityMatrix();

        auto checkHistogram = [nrSamples](const auto &counts, const std::vector<double> &probabilities,
                                          const std::string &description) -> bool {
            size_t total = 0;
            for (const auto &measurement : counts)
                total += measurement.second;
            if (total != static_cast<size_t>(nrSamples))
            {
                std::cout << description << " returned " << total << " samples instead of " << nrSamples << std::endl;
                return false;
            }

            for (size_t s = 0; s < probabilities.size(); ++s)
            {
                const auto it = counts.find(s);
                const size_t count = it == counts.end() ? 0 : it->second;
                const double sampled = static_cast<double>(count) / nrSamples;
                if (std::abs(sampled - probabilities[s]) > 0.02)
                {
                    std::cout << description << " mismatch for outcome " << s << ": " << probabilities[s]
                              << " (exact) vs " << sampled << " (sampled)" << std::endl;
                    return false;
                }
            }

            return true;
        };

        const auto orderedCounts = dm.RepeatedMeasure(nrSamples);
        if (!checkHistogram(orderedCounts, expected, "Ordered full-register sampling"))
            return false;

        const auto unorderedCounts = dm.RepeatedMeasureUnordered(nrSamples);
        if (!checkHistogram(unorderedCounts, expected, "Unordered full-register sampling"))
            return false;

        // sample a contiguous subregister and compare against its exact marginal. Starting at qubit one
        // when possible also checks that the returned outcomes are shifted and packed correctly.
        {
            const size_t firstQubit = nrQubits > 1 ? 1 : 0;
            const size_t secondQubit = static_cast<size_t>(nrQubits - 1);
            const size_t firstPartMask = (1ULL << firstQubit) - 1;
            const size_t measuredPartMask = (1ULL << (secondQubit + 1)) - 1 - firstPartMask;
            const size_t subregisterStates = 1ULL << (secondQubit - firstQubit + 1);
            std::vector<double> expectedSubregister(subregisterStates, 0.);
            for (size_t s = 0; s < nrBasisStates; ++s)
                expectedSubregister[(s & measuredPartMask) >> firstQubit] += expected[s];

            const auto orderedSubregisterCounts = dm.RepeatedMeasure(firstQubit, secondQubit, nrSamples);
            if (!checkHistogram(orderedSubregisterCounts, expectedSubregister, "Ordered subregister sampling"))
                return false;

            const auto unorderedSubregisterCounts = dm.RepeatedMeasureUnordered(firstQubit, secondQubit, nrSamples);
            if (!checkHistogram(unorderedSubregisterCounts, expectedSubregister, "Unordered subregister sampling"))
                return false;
        }

        if (!dm.RepeatedMeasure(0).empty() || !dm.RepeatedMeasureUnordered(0).empty() ||
            !dm.RepeatedMeasure(0, 0, 0).empty() || !dm.RepeatedMeasureUnordered(0, 0, 0).empty())
        {
            std::cout << "Zero-shot sampling returned a non-empty histogram" << std::endl;
            return false;
        }

        if ((dm.getDensityMatrix() - stateBeforeSampling).norm() > 1E-12)
        {
            std::cout << "Repeated sampling changed the state for " << nrQubits << " qubits" << std::endl;
            return false;
        }

        std::cout << "." << std::flush;
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

static bool DensityMatrixStateSaveRestoreCloneTest()
{
    std::cout << "\nDensity matrix simulator - state save, restore and clone" << std::endl;

    QC::DensityMatrix<> dm(3);
    QC::Gates::HadamardGate<> h;
    QC::Gates::CNOTGate<> cnot;
    dm.ApplyGate(h, 0);
    dm.ApplyGate(cnot, 2, 0);
    dm.ApplyDepolarizingNoise(1, 0.3);

    const Eigen::MatrixXcd savedState = dm.getDensityMatrix();
    dm.SaveState();

    dm.setToBasisState(7);
    const Eigen::MatrixXcd stateAtClone = dm.getDensityMatrix();
    auto clone = dm.Clone();
    if ((clone->getDensityMatrix() - stateAtClone).norm() > 1E-12)
    {
        std::cout << "Density matrix clone did not copy the current state" << std::endl;
        return false;
    }

    // A regular restore keeps the snapshot available for later restores.
    dm.RestoreState();
    if ((dm.getDensityMatrix() - savedState).norm() > 1E-12)
    {
        std::cout << "Density matrix RestoreState did not restore the snapshot" << std::endl;
        return false;
    }
    dm.setToBasisState(7);
    dm.RestoreState();
    if ((dm.getDensityMatrix() - savedState).norm() > 1E-12)
    {
        std::cout << "Density matrix RestoreState consumed the snapshot" << std::endl;
        return false;
    }

    // A destructive restore consumes the original's snapshot.
    dm.setToBasisState(7);
    dm.RestoreStateDestructive();
    if ((dm.getDensityMatrix() - savedState).norm() > 1E-12)
    {
        std::cout << "Density matrix RestoreStateDestructive did not restore the snapshot" << std::endl;
        return false;
    }
    dm.setToBasisState(0);
    const Eigen::MatrixXcd stateWithoutSnapshot = dm.getDensityMatrix();
    dm.RestoreStateDestructive();
    if ((dm.getDensityMatrix() - stateWithoutSnapshot).norm() > 1E-12)
    {
        std::cout << "Density matrix destructive restore did not consume the snapshot" << std::endl;
        return false;
    }

    // The clone owns an independent copy of both the current state and the saved snapshot.
    clone->RestoreState();
    if ((clone->getDensityMatrix() - savedState).norm() > 1E-12)
    {
        std::cout << "Density matrix clone did not copy the saved snapshot" << std::endl;
        return false;
    }
    clone->setToBasisState(7);
    if ((dm.getDensityMatrix() - stateWithoutSnapshot).norm() > 1E-12)
    {
        std::cout << "Changing a density matrix clone changed the original" << std::endl;
        return false;
    }
    clone->RestoreStateDestructive();
    if ((clone->getDensityMatrix() - savedState).norm() > 1E-12)
    {
        std::cout << "Density matrix clone's destructive restore failed" << std::endl;
        return false;
    }

    // A clone is a complete simulator snapshot, including the random-engine position.
    QC::DensityMatrix<> sampler(1);
    sampler.ApplyGate(h, 0);
    sampler.SetRandomSeed(0xC10EULL);
    (void)sampler.MeasureNoCollapse(); // advance the engine before cloning
    auto samplerClone = sampler.Clone();
    if (sampler.RepeatedMeasure(2000) != samplerClone->RepeatedMeasure(2000))
    {
        std::cout << "Density matrix clone did not preserve the random-engine state" << std::endl;
        return false;
    }

    std::cout << "Success" << std::endl;
    return true;
}

// smallest eigenvalue of a Hermitian matrix (used to check positive semidefiniteness)
static double DM_SmallestEigenvalue(const Eigen::MatrixXcd &rho)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(rho);
    return solver.eigenvalues().minCoeff();
}

// unitary + noise evolution: unit trace, hermiticity, positive semidefiniteness, purity bounds,
// plus a few analytic limiting cases (full depolarizing, full amplitude damping, reset idempotence)
static bool DensityMatrixInvariantsTest()
{
    std::cout << "\nDensity matrix simulator - physics invariants (trace, hermiticity, positivity, purity)"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    DM_FillOneQubitGates(gates);
    DM_FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(20, 40);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);
    std::uniform_int_distribution noiseDistr(0, 5);
    std::uniform_real_distribution<double> probDistr(0.05, 0.5);

    for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

        for (int t = 0; t < 10; ++t)
        {
            QC::DensityMatrix<> dm(nrQubits);

            const int lim = nrGatesDistr(gen);
            for (int i = 0; i < lim; ++i)
            {
                const int g = gateDistr(gen);
                const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
                if (twoQubitsGate && nrQubits < 2)
                    continue;

                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = qubit1 + 1;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

                dm.ApplyGate(*gates[g], qubit1, qubit2);

                // sprinkle in a random noise channel every few gates
                if (i % 4 == 3)
                {
                    const int q = qubitDistr(gen);
                    const double p = probDistr(gen);
                    switch (noiseDistr(gen))
                    {
                    case 0:
                        dm.ApplyBitFlipNoise(q, p);
                        break;
                    case 1:
                        dm.ApplyPhaseFlipNoise(q, p);
                        break;
                    case 2:
                        dm.ApplyDepolarizingNoise(q, p);
                        break;
                    case 3:
                        dm.ApplyAmplitudeDamping(q, p);
                        break;
                    case 4:
                        dm.ApplyPhaseDamping(q, p);
                        break;
                    default:
                        dm.ApplyReset(q);
                        break;
                    }
                }
            }

            const auto &rho = dm.getDensityMatrix();

            // trace must remain 1 (real)
            if (!approxEqual(dm.Trace().real(), 1., 1E-9) || !approxEqual(dm.Trace().imag(), 0., 1E-9))
            {
                std::cout << "Trace is not 1 after noisy evolution: " << dm.Trace() << std::endl;
                return false;
            }

            // must stay Hermitian
            if (!dm.IsHermitian(1E-9))
            {
                std::cout << "Density matrix is not Hermitian after noisy evolution" << std::endl;
                return false;
            }

            // must stay positive semidefinite (allow a tiny negative numerical slack)
            const double minEig = DM_SmallestEigenvalue(rho);
            if (minEig < -1E-9)
            {
                std::cout << "Density matrix is not positive semidefinite, smallest eigenvalue: " << minEig
                          << std::endl;
                return false;
            }

            // purity must lie within [1/2^N, 1]
            const double purity = dm.Purity();
            const double minPurity = 1. / static_cast<double>(dm.getNrBasisStates());
            if (purity < minPurity - 1E-9 || purity > 1. + 1E-9)
            {
                std::cout << "Purity " << purity << " out of bounds [" << minPurity << ", 1] for " << nrQubits
                          << " qubits" << std::endl;
                return false;
            }
        }
    }

    // analytic limiting cases on a single qubit
    {
        // depolarizing at p = 3/4 sends any state to the maximally mixed I/2
        // (for the {I, X, Y, Z} Kraus parametrization the fixed point is p = 3/4, not p = 1)
        QC::DensityMatrix<> dm(1);
        QC::Gates::HadamardGate<> h;
        dm.ApplyGate(h, 0);
        dm.ApplyDepolarizingNoise(0, 0.75);

        const auto &rho = dm.getDensityMatrix();
        if (!approxEqual(rho(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho(1, 1), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
            !approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
        {
            std::cout << "Depolarizing at p = 3/4 did not produce the maximally mixed state:\n" << rho << std::endl;
            return false;
        }

        if (!approxEqual(dm.Purity(), 0.5, 1E-9))
        {
            std::cout << "Maximally mixed single qubit purity is not 1/2: " << dm.Purity() << std::endl;
            return false;
        }
    }

    {
        // full amplitude damping (gamma = 1) sends |1> to |0>
        QC::DensityMatrix<> dm(1);
        QC::Gates::PauliXGate<> x;
        dm.ApplyGate(x, 0); // |1>
        dm.ApplyAmplitudeDamping(0, 1.);

        const auto &rho = dm.getDensityMatrix();
        if (!approxEqual(rho(0, 0), std::complex<double>(1., 0), 1E-9) ||
            !approxEqual(rho(1, 1), std::complex<double>(0., 0), 1E-9))
        {
            std::cout << "Full amplitude damping did not relax |1> to |0>:\n" << rho << std::endl;
            return false;
        }
    }

    {
        // reset is idempotent: applying it twice equals applying it once
        QC::DensityMatrix<> dm1(2);
        QC::DensityMatrix<> dm2(2);
        QC::Gates::HadamardGate<> h;
        QC::Gates::CNOTGate<> cnot;
        for (auto *dm : {&dm1, &dm2})
        {
            dm->ApplyGate(h, 0);
            dm->ApplyGate(cnot, 1, 0);
        }

        dm1.ApplyReset(0);
        dm2.ApplyReset(0);
        dm2.ApplyReset(0); // second reset must be a no-op

        const auto &r1 = dm1.getDensityMatrix();
        const auto &r2 = dm2.getDensityMatrix();
        if ((r1 - r2).norm() > 1E-9)
        {
            std::cout << "Reset is not idempotent" << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

static bool DensityMatrixNewFeaturesTest()
{
    std::cout << "\nDensity matrix simulator - bitvector initialization, partial trace, overlap and fidelity"
              << std::endl;

    // 1. Bitvector initialization overloads
    {
        QC::DensityMatrix<> dm(2);
        dm.setToBasisState(std::vector<bool>{true, false}); // |01> state (qubit 0 = 1, qubit 1 = 0) => state 1
        if (!approxEqual(dm.getBasisStateProbability(1), 1.0, 1E-9))
        {
            std::cout << "setToBasisState(vector<bool>) failed" << std::endl;
            return false;
        }

        dm.setToMixtureOfBasisStates(
            std::vector<std::pair<std::vector<bool>, double>>{{{false, false}, 0.4}, {{true, true}, 0.6}});
        if (!approxEqual(dm.getBasisStateProbability(0), 0.4, 1E-9) ||
            !approxEqual(dm.getBasisStateProbability(3), 0.6, 1E-9))
        {
            std::cout << "setToMixtureOfBasisStates(vector<bool>) failed" << std::endl;
            return false;
        }
    }

    // 2. Partial trace: Bell state (|00> + |11>)/sqrt(2) -> trace out qubit 1 -> reduced state on qubit 0 is I/2
    {
        QC::DensityMatrix<> dm(2);
        QC::Gates::HadamardGate<> h;
        QC::Gates::CNOTGate<> cnot;
        dm.ApplyGate(h, 0);
        dm.ApplyGate(cnot, 1, 0);

        const Eigen::MatrixXcd rho0 = dm.PartialTrace({0}); // keep qubit 0
        if (!approxEqual(rho0(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho0(1, 1), std::complex<double>(0.5, 0), 1E-9) ||
            !approxEqual(rho0(0, 1), std::complex<double>(0, 0), 1E-9))
        {
            std::cout << "PartialTrace of Bell state failed:\n" << rho0 << std::endl;
            return false;
        }
    }

    // 3. Overlap and Fidelity
    {
        QC::DensityMatrix<> dm1(1);
        QC::DensityMatrix<> dm2(1);
        QC::Gates::HadamardGate<> h;
        dm1.ApplyGate(h, 0); // |+>
        dm2.ApplyGate(h, 0); // |+>

        if (!approxEqual(dm1.HilbertSchmidtOverlap(dm2).real(), 1.0, 1E-9))
        {
            std::cout << "HilbertSchmidtOverlap of identical pure states is not 1" << std::endl;
            return false;
        }

        QC::QubitRegister<> reg(1);
        reg.ApplyGate(h, 0);
        if (!approxEqual(dm1.FidelityWithStatevector(reg.getRegisterStorage()), 1.0, 1E-9))
        {
            std::cout << "FidelityWithStatevector with identical state is not 1" << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;
    return true;
}

template <class Matrix> class DM_ForcedMeasurement : public QC::DensityMatrix<Eigen::VectorXcd, Matrix>
{
  public:
    using QC::DensityMatrix<Eigen::VectorXcd, Matrix>::DensityMatrix;

    // Exercise a rare branch deterministically without searching billions of seeds.
    void SelectOutcome(size_t outcome, double probability)
    {
        this->uniformZeroOne = outcome == 0 ? std::uniform_real_distribution<double>(0., probability * 0.5)
                                            : std::uniform_real_distribution<double>(1. - probability * 0.5, 1.);
    }
};

template <class Matrix> static bool DM_RareCollapseChecks(bool parallel)
{
    using Complex = std::complex<double>;
    for (size_t qubit = 0; qubit < 3; ++qubit)
        for (size_t outcome = 0; outcome < 2; ++outcome)
        {
            const double probability = 1E-8;
            const size_t first = outcome << qubit;
            const size_t otherQubit = (qubit + 1) % 3;
            const size_t second = first | (size_t{1} << otherQubit);
            Eigen::VectorXcd branch = Eigen::VectorXcd::Zero(8);
            branch[first] = std::sqrt(0.3);
            branch[second] = Complex(0., std::sqrt(0.7));
            const Matrix expected = branch * branch.adjoint();
            Matrix rho = probability * expected;
            rho(first ^ (size_t{1} << qubit), first ^ (size_t{1} << qubit)) = 1. - probability;

            // Model the anti-Hermitian residual from separate ket/bra gate
            // passes. It is below validation tolerance before collapse.
            rho(first, first) += Complex(0., 1E-17);
            rho(second, second) -= Complex(0., 1E-17);
            const Complex error(2E-17, 3E-17);
            rho(first, second) += error;
            rho(second, first) -= std::conj(error);

            DM_ForcedMeasurement<Matrix> dm(3, 1);
            dm.SetMultithreading(parallel);
            dm.SetRandomSeed(99);
            dm.setDensityMatrix(rho);
            dm.SelectOutcome(outcome, probability);
            if (dm.MeasureQubit(qubit) != outcome || !dm.getDensityMatrix().allFinite() ||
                (dm.getDensityMatrix() - expected).norm() > 1E-12 || !dm.IsHermitian(1E-14))
            {
                std::cout << "Rare collapse amplified non-Hermitian roundoff" << std::endl;
                return false;
            }
            // An off-diagonal residual also becomes a population after a gate.
            dm.ApplyGate(QC::Gates::HadamardGate<Matrix>(), otherQubit);
            dm.MeasureNoCollapse();
            dm.MeasureQubit(otherQubit);
        }

    // Positive probabilities smaller than the old 1E-20 cutoff must still
    // normalize, including when computing a reciprocal would overflow.
    for (const double probability : {1E-30, 1E-310})
    {
        Matrix rho = Matrix::Zero(2, 2);
        rho(0, 0) = probability;
        rho(1, 1) = 1.;
        DM_ForcedMeasurement<Matrix> dm(1, 1);
        dm.SetMultithreading(parallel);
        dm.SetRandomSeed(99);
        dm.setDensityMatrix(rho);
        dm.SelectOutcome(0, probability);
        if (dm.MeasureQubit(0) != 0 || !dm.getDensityMatrix().allFinite() || dm.getBasisStateProbability(0) != 1.)
        {
            std::cout << "Tiny-probability collapse lost normalization" << std::endl;
            return false;
        }
    }
    return true;
}

template <class Matrix> static bool DM_TiledCollapseChecks(bool parallel)
{
    using Complex = std::complex<double>;
    // Cover partial tiles, a full tile, multiple tiles, and every placement of
    // the measured bit relative to the contiguous storage dimension.
    for (size_t qubits : {1, 4, 5, 6, 8})
    {
        const size_t dimension = size_t{1} << qubits;
        Eigen::VectorXcd psi(dimension);
        for (size_t i = 0; i < dimension; ++i)
            psi(i) = Complex(std::sin(double(i) + 0.3), std::cos(0.7 * double(i)));
        psi.normalize();
        Matrix initial = 0.8 * (psi * psi.adjoint());
        initial.diagonal().array() += 0.2 / double(dimension);
        // Exercise projection of both diagonal and off-diagonal roundoff.
        for (size_t row = 0; row < dimension; ++row)
        {
            initial(row, row) += Complex(0., 1E-17);
            for (size_t column = row + 1; column < dimension; ++column)
            {
                const Complex error(2E-17, -3E-17);
                initial(row, column) += error;
                initial(column, row) -= std::conj(error);
            }
        }
        DM_ForcedMeasurement<Matrix> dm(qubits, 1);
        dm.SetMultithreading(parallel);
        for (size_t qubit = 0; qubit < qubits; ++qubit)
            for (size_t outcome = 0; outcome < 2; ++outcome)
            {
                const size_t mask = size_t{1} << qubit;
                double probability = 0.;
                for (size_t i = 0; i < dimension; ++i)
                    if (((i & mask) != 0) == (outcome != 0))
                        probability += initial(i, i).real();
                Matrix expected = Matrix::Zero(dimension, dimension);
                for (size_t row = 0; row < dimension; ++row)
                    for (size_t column = 0; column < dimension; ++column)
                        if (((row & mask) != 0) == (outcome != 0) && ((column & mask) != 0) == (outcome != 0))
                            expected(row, column) =
                                (0.5 * initial(row, column) + 0.5 * std::conj(initial(column, row))) / probability;
                dm.setDensityMatrix(initial);
                dm.SetRandomSeed(99);
                dm.SelectOutcome(outcome, probability);
                if (dm.MeasureQubit(qubit) != outcome || !dm.getDensityMatrix().allFinite() ||
                    (dm.getDensityMatrix() - expected).norm() > 1E-14 || !dm.IsHermitian(1E-14))
                {
                    std::cout << "Tiled collapse differs from the projected reference: " << qubits << " qubits, q"
                              << qubit << ", outcome " << outcome << std::endl;
                    return false;
                }
            }
    }
    return true;
}

static bool DensityMatrixRareCollapseTest()
{
    std::cout << "\nDensity matrix simulator - rare measurement collapse" << std::endl;

    struct RestoreThreshold
    {
        size_t value = QC::DensityMatrix<>::GetParallelMinElements();
        int team = QC::DensityMatrix<>::GetParallelThreads();
#ifdef _OPENMP
        int threads = omp_get_max_threads();
#endif

        ~RestoreThreshold()
        {
            QC::DensityMatrix<>::SetParallelMinElements(value);
            QC::DensityMatrix<>::SetParallelThreads(team);
#ifdef _OPENMP
            omp_set_num_threads(threads);
#endif
        }
    } restoreThreshold;

    QC::DensityMatrix<>::SetParallelMinElements(1);
    using RowMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    for (bool parallel : {false, true})
        if (!DM_RareCollapseChecks<Eigen::MatrixXcd>(parallel) || !DM_RareCollapseChecks<RowMatrix>(parallel) ||
            !DM_TiledCollapseChecks<Eigen::MatrixXcd>(parallel) || !DM_TiledCollapseChecks<RowMatrix>(parallel))
            return false;

#ifdef _OPENMP
    // Tile sizing must also work with odd team sizes and with fewer physical
    // tile rows than requested workers; no particular machine is assumed.
    // A fixed team, since the automatic one has two threads for these small matrices.
    for (int threads : {1, 2, 3, 7})
    {
        omp_set_num_threads(threads);
        QC::DensityMatrix<>::SetParallelThreads(threads);
        if (!DM_TiledCollapseChecks<Eigen::MatrixXcd>(true) || !DM_TiledCollapseChecks<RowMatrix>(true))
            return false;
    }
    omp_set_num_threads(restoreThreshold.threads);
    QC::DensityMatrix<>::SetParallelThreads(restoreThreshold.team);
    bool nestedResults[2]{};
#pragma omp parallel for num_threads(2)
    for (int worker = 0; worker < 2; ++worker)
        nestedResults[worker] =
            DM_TiledCollapseChecks<Eigen::MatrixXcd>(true) && DM_TiledCollapseChecks<RowMatrix>(true);
    if (!nestedResults[0] || !nestedResults[1])
        return false;
#endif

    // Reproduce actual gate roundoff as well as the constructed residual above.
    DM_ForcedMeasurement<Eigen::MatrixXcd> dm(2, 1);
    dm.SetMultithreading(false);
    dm.SetRandomSeed(99);
    dm.setToBasisState(1);
    dm.ApplyGate(QC::Gates::HadamardGate<>(), 1);
    std::mt19937_64 random(99);
    std::uniform_real_distribution<double> angle(-3., 3.);
    std::vector<Eigen::MatrixXcd> gates;
    for (size_t i = 0; i < 100; ++i)
    {
        const double theta = angle(random), phi = angle(random), lambda = angle(random), gamma = angle(random);
        QC::Gates::UGate<> gate(theta, phi, lambda, gamma);
        gates.push_back(gate.getRawOperatorMatrix());
        dm.ApplyGate(gate, 0);
    }
    for (auto it = gates.rbegin(); it != gates.rend(); ++it)
        dm.ApplyGate(QC::Gates::AppliedGate<>(it->adjoint().eval(), 0));
    dm.ApplyGate(QC::Gates::RxGate<>(2E-4), 0);
    dm.SelectOutcome(0, 1E-8);
    if (dm.MeasureQubit(0) != 0 || !dm.IsHermitian(1E-14))
        return false;
    dm.MeasureQubit(1);
    std::cout << "Success" << std::endl;
    return true;
}

// From DensityMatrix<>::OneQubitSuperoperatorQubits (dense, 11 by default) and
// AntidiagonalSuperoperatorQubits (5) on, non-diagonal one-qubit gates are applied as their
// superoperator in one pass instead of two. Checks that path at the default thresholds against the
// statevector, for column- and row-major storage.
static bool DensityMatrixOneQubitSuperoperatorTest()
{
    std::cout << "\nDensity matrix simulator - one-qubit gates as superoperators on large registers" << std::endl;

    using RowMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> oneQubitGates, gates;
    DM_FillOneQubitGates(oneQubitGates);
    gates = oneQubitGates;
    DM_FillTwoQubitGates(gates);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    // The default thresholds, as fixed sizes: an override that raises them keeps this test bounded
    // (three density matrices of 4^n entries) and checks whichever path then applies.
    for (const size_t size : {size_t{5}, size_t{11}})
    {
        const int nrQubits = static_cast<int>(std::max<size_t>(2, size));
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        QC::DensityMatrix<> dm(nrQubits);
        QC::DensityMatrix<Eigen::VectorXcd, RowMatrix> rowDm(nrQubits);
        QC::QubitRegister<> reg(nrQubits);
        const auto apply = [&](const QC::Gates::QuantumGateWithOp<> &gate, int qubit1, int qubit2) {
            dm.ApplyGate(gate, qubit1, qubit2);
            // A row-major register takes row-major gates: the same matrix, whose structure is classified again.
            rowDm.ApplyGate(QC::Gates::AppliedGate<RowMatrix>(RowMatrix(gate.getRawOperatorMatrix()), qubit1, qubit2));
            reg.ApplyGate(gate, qubit1, qubit2);
        };

        // Every one-qubit gate at least once (dense, antidiagonal and diagonal), then a random circuit.
        for (const auto &gate : oneQubitGates)
            apply(*gate, qubitDistr(gen), 0);
        for (int i = 0; i < 40; ++i)
        {
            const auto &gate = *gates[gateDistr(gen)];
            if (gate.getQubitsNumber() == 2)
            {
                int qubit1 = qubitDistr2(gen), qubit2 = qubit1 + 1;
                if (dist_bool(gen))
                    std::swap(qubit1, qubit2);
                apply(gate, qubit1, qubit2);
            }
            else
                apply(gate, qubitDistr(gen), 0);
        }

        if (!CompareWithStatevector(dm, reg, nrQubits))
            return false;
        const Eigen::MatrixXcd &columns = dm.getDensityMatrix();
        const RowMatrix &rows = rowDm.getDensityMatrix();
        for (Eigen::Index r = 0; r < columns.rows(); ++r)
            for (Eigen::Index c = 0; c < columns.cols(); ++c)
                if (!approxEqual(columns(r, c), rows(r, c), 1E-9))
                {
                    std::cout << "Row-major density matrix differs at (" << r << ", " << c << ") for " << nrQubits
                              << " qubits" << std::endl;
                    return false;
                }
    }

    std::cout << "Success" << std::endl;
    return true;
}

// Gates need not be unitary. A superoperator coefficient is a product of two operator entries and can
// overflow, or round to zero, where the two passes stay in range: the gate then takes the two passes.
// An antidiagonal and a dense operator with an entry of 1e200, each at the size where it would take the
// superoperator, then a sequence whose middle operator has entries of 1e-200.
static bool DensityMatrixNonUnitaryOverflowTest()
{
    std::cout << "\nDensity matrix simulator - non-unitary one-qubit operators with large entries" << std::endl;

    struct Case
    {
        int nrQubits;
        Eigen::MatrixXcd op;
        bool startInOne;  // prepare |1> on qubit 0 instead of |0>
        Eigen::Index out; // the basis state the operator maps it to
    };

    Eigen::MatrixXcd antidiagonal(2, 2), dense(2, 2);
    antidiagonal << 0., 1e200, 1., 0.; // |0> -> |1>
    dense << 1e200, 1., 1., 0.;        // |1> -> |0>
    const std::array<Case, 2> cases{{{5, antidiagonal, false, 1}, {11, dense, true, 0}}};

    for (const auto &c : cases)
    {
        QC::DensityMatrix<> dm(c.nrQubits);
        if (c.startInOne)
            dm.ApplyGate(QC::Gates::PauliXGate<>(), 0);
        dm.ApplyGate(QC::Gates::AppliedGate<>(c.op, 0), 0);
        const Eigen::MatrixXcd &rho = dm.getDensityMatrix();
        if (!rho.allFinite())
        {
            std::cout << "A finite non-unitary operator gave a non-finite density matrix for " << c.nrQubits
                      << " qubits" << std::endl;
            return false;
        }
        Eigen::MatrixXcd expected = Eigen::MatrixXcd::Zero(rho.rows(), rho.cols());
        expected(c.out, c.out) = 1.;
        if ((rho - expected).cwiseAbs().maxCoeff() > 1E-12)
        {
            std::cout << "A non-unitary operator with large entries gave the wrong density matrix for " << c.nrQubits
                      << " qubits" << std::endl;
            return false;
        }
    }

    // Underflow: 1e100 X, 1e-200 X, 1e100 X multiply to X. The coefficient 1e-200 * 1e-200 of the
    // middle one rounds to zero, while the two passes scale the 1e200 population down in two steps.
    Eigen::MatrixXcd x(2, 2);
    x << 0., 1., 1., 0.;
    for (const int nrQubits : {5, 11})
    {
        QC::DensityMatrix<> dm(nrQubits);
        for (const double scale : {1e100, 1e-200, 1e100})
            dm.ApplyGate(QC::Gates::AppliedGate<>(Eigen::MatrixXcd(scale * x), 0), 0);
        const Eigen::MatrixXcd &rho = dm.getDensityMatrix();
        Eigen::MatrixXcd expected = Eigen::MatrixXcd::Zero(rho.rows(), rho.cols());
        expected(1, 1) = 1.;
        if (!rho.allFinite() || (rho - expected).cwiseAbs().maxCoeff() > 1E-12)
        {
            std::cout << "Operators with very small and very large entries gave the wrong density matrix for "
                      << nrQubits << " qubits: population of |1> " << rho(1, 1).real() << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;
    return true;
}

bool DensityMatrixTests()
{
    std::cout << "\nDensity matrix simulator tests" << std::endl;

    const std::mt19937 savedGenerator = gen;
    gen.seed(0xD35A17U);
    const bool result =
        DensityMatrixStructuredGateTests() && DensityMatrixUnitaryTest() && DensityMatrixDenseGateReferenceTest() &&
        DensityMatrixThreeQubitGateTest() && DensityMatrixChannelsTest() && DensityMatrixGenericChannelTest() &&
        DensityMatrixSuperoperatorTest() && StatevectorKrausTrajectoryTest() &&
        DensityMatrixInitializationAndValidationTest() && DensityMatrixMeasurementTest() &&
        DensityMatrixPauliExpectationTest() && DensityMatrixInvariantsTest() &&
        DensityMatrixMeasurementStatisticsTest() && DensityMatrixCollapseTest() && DensityMatrixRareCollapseTest() &&
        DensityMatrixExpectationExtrasTest() && DensityMatrixSamplingTest() &&
        DensityMatrixStateSaveRestoreCloneTest() && DensityMatrixNewFeaturesTest() &&
        DensityMatrixOneQubitSuperoperatorTest() && DensityMatrixNonUnitaryOverflowTest();
    gen = savedGenerator;
    return result;
}
