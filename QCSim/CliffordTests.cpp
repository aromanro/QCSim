#include "Clifford.h"
#include "ExtendedStabilizer.h"
#include "Tests.h"

#include "PauliStringXZCoeff.h"
#include "QubitRegister.h"

#include <cstdlib>
#include <future>

std::shared_ptr<QC::Gates::QuantumGateWithOp<>> GetTwoQubitsGate(int code)
{
    switch (code)
    {
    case 9:
        return std::make_shared<QC::Gates::CNOTGate<>>();
    case 10:
        return std::make_shared<QC::Gates::ControlledYGate<>>();
    case 11:
        return std::make_shared<QC::Gates::ControlledZGate<>>();
    case 12:
        return std::make_shared<QC::Gates::SwapGate<>>();
    case 13:
        return std::make_shared<QC::Gates::iSwapGate<>>();
    case 14:
        return std::make_shared<QC::Gates::iSwapDagGate<>>();
    }

    return nullptr;
}

std::shared_ptr<QC::Gates::QuantumGateWithOp<>> GetGate(int code, double param)
{
    switch (code)
    {
    case 0:
        return std::make_shared<QC::Gates::HadamardGate<>>();
    case 1:
        return std::make_shared<QC::Gates::SGate<>>();
    case 2:
        return std::make_shared<QC::Gates::SDGGate<>>();
    case 3:
        return std::make_shared<QC::Gates::PauliXGate<>>();
    case 4:
        return std::make_shared<QC::Gates::PauliYGate<>>();
    case 5:
        return std::make_shared<QC::Gates::PauliZGate<>>();
    case 6:
        return std::make_shared<QC::Gates::SquareRootNOTGate<>>();
    case 7:
        return std::make_shared<QC::Gates::SquareRootNOTDagGate<>>();
    case 8:
        return std::make_shared<QC::Gates::HyGate<>>();
    // for non-clifford
    case 15:
        return std::make_shared<QC::Gates::RxGate<>>(param);
    case 16:
        return std::make_shared<QC::Gates::RyGate<>>(param);
    case 17:
        return std::make_shared<QC::Gates::RzGate<>>(param);
    default:
        return GetTwoQubitsGate(code);
    }

    return nullptr;
}

void ApplyTwoQubitsGate(QC::Clifford::StabilizerSimulator &simulator, int code, int qubit1, int qubit2)
{
    switch (code)
    {
    case 9:
        simulator.ApplyCX(qubit1, qubit2);
        break;
    case 10:
        simulator.ApplyCY(qubit1, qubit2);
        break;
    case 11:
        simulator.ApplyCZ(qubit1, qubit2);
        break;
    case 12:
        simulator.ApplySwap(qubit1, qubit2);
        break;
    case 13:
        simulator.ApplyISwap(qubit1, qubit2);
        break;
    case 14:
        simulator.ApplyISwapDag(qubit1, qubit2);
        break;
    }
}

void ApplyGate(QC::Clifford::StabilizerSimulator &simulator, int code, int qubit1, int qubit2)
{
    switch (code)
    {
    case 0:
        simulator.ApplyH(qubit1);
        break;
    case 1:
        simulator.ApplyS(qubit1);
        break;
    case 2:
        simulator.ApplySdg(qubit1);
        break;
    case 3:
        simulator.ApplyX(qubit1);
        break;
    case 4:
        simulator.ApplyY(qubit1);
        break;
    case 5:
        simulator.ApplyZ(qubit1);
        break;
    case 6:
        simulator.ApplySx(qubit1);
        break;
    case 7:
        simulator.ApplySxDag(qubit1);
        break;
    case 8:
        simulator.ApplyK(qubit1);
        break;
    default:
        ApplyTwoQubitsGate(simulator, code, qubit1, qubit2);
        break;
    }
}

void PrintTwoQubitsGate(int code, int qubit1, int qubit2)
{
    switch (code)
    {
    case 9:
        std::cout << "CX " << qubit1 << " " << qubit2 << std::endl;
        break;
    case 10:
        std::cout << "CY " << qubit1 << " " << qubit2 << std::endl;
        break;
    case 11:
        std::cout << "CZ " << qubit1 << " " << qubit2 << std::endl;
        break;
    case 12:
        std::cout << "SWAP " << qubit1 << " " << qubit2 << std::endl;
        break;
    case 13:
        std::cout << "iSWAP " << qubit1 << " " << qubit2 << std::endl;
        break;
    case 14:
        std::cout << "iSWAPDag " << qubit1 << " " << qubit2 << std::endl;
        break;
    }
}

void PrintGate(int code, int qubit1, int qubit2)
{
    switch (code)
    {
    case 0:
        std::cout << "H " << qubit1 << std::endl;
        break;
    case 1:
        std::cout << "S " << qubit1 << std::endl;
        break;
    case 2:
        std::cout << "SDG " << qubit1 << std::endl;
        break;
    case 3:
        std::cout << "X " << qubit1 << std::endl;
        break;
    case 4:
        std::cout << "Y " << qubit1 << std::endl;
        break;
    case 5:
        std::cout << "Z " << qubit1 << std::endl;
        break;
    case 6:
        std::cout << "SX " << qubit1 << std::endl;
        break;
    case 7:
        std::cout << "SXDG " << qubit1 << std::endl;
        break;
    case 8:
        std::cout << "K " << qubit1 << std::endl;
        break;
    default:
        PrintTwoQubitsGate(code, qubit1, qubit2);
    }
}

void ConstructCircuit(size_t nrQubits, std::vector<int> &gates, std::vector<size_t> &qubits1,
                      std::vector<size_t> &qubits2, std::uniform_int_distribution<int> &gateDistr,
                      std::uniform_int_distribution<int> &qubitDistr)
{
    for (int i = 0; i < static_cast<int>(gates.size()); ++i)
    {
        gates[i] = gateDistr(gen);

        qubits1[i] = qubitDistr(gen);
        qubits2[i] = qubitDistr(gen);

        if (qubits2[i] == qubits1[i])
            qubits2[i] = (qubits1[i] + 1) % static_cast<int>(nrQubits);

        if (dist_bool(gen))
            std::swap(qubits1[i], qubits2[i]);
    }
}

void ExecuteCircuit(size_t nrShots, size_t nrQubits, const std::vector<int> &gates, const std::vector<size_t> &qubits1,
                    const std::vector<size_t> &qubits2, std::unordered_map<size_t, int> &results1,
                    std::unordered_map<size_t, int> &results2)
{
    size_t remainingCounts = nrShots;
    size_t nrThreads = QC::QubitRegisterCalculator<>::GetNumberOfThreads();
    nrThreads = std::min(nrThreads, std::max<size_t>(remainingCounts, 1ULL));

    std::vector<std::future<void>> tasks(nrThreads);

    const size_t cntPerThread = static_cast<size_t>(ceil(static_cast<double>(remainingCounts) / nrThreads));

    std::mutex resultsMutex;

    for (size_t th = 0; th < nrThreads; ++th)
    {
        const size_t curCnt = std::min(cntPerThread, remainingCounts);
        remainingCounts -= curCnt;

        tasks[th] = std::async(std::launch::async,
                               [&gates, &qubits1, &qubits2, &results1, &results2, curCnt, nrQubits, &resultsMutex]() {
                                   for (int i = 0; i < static_cast<int>(curCnt); ++i)
                                   {
                                       QC::QubitRegister qubitRegister(nrQubits);
                                       QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

                                       for (int j = 0; j < static_cast<int>(gates.size()); ++j)
                                       {
                                           ApplyGate(cliffordSim, gates[j], qubits1[j], qubits2[j]);
                                           const auto gateptr = GetGate(gates[j]);
                                           qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
                                       }

                                       // now do the measurements
                                       size_t val1 = 0;
                                       size_t val2 = 0;
                                       for (int q = 0; q < static_cast<int>(nrQubits); ++q)
                                       {
                                           val1 <<= 1;
                                           val2 <<= 1;

                                           if (cliffordSim.MeasureQubit(q))
                                               val1 |= 1;
                                           if (qubitRegister.MeasureQubit(q))
                                               val2 |= 1;
                                       }

                                       {
                                           const std::lock_guard lock(resultsMutex);
                                           ++results1[val1];
                                           ++results2[val2];
                                       }
                                   }
                               });
    }

    for (size_t i = 0; i < nrThreads; ++i)
        tasks[i].get();
}

bool CheckProbability(QC::QubitRegister<> &qubitRegister, QC::Clifford::StabilizerSimulator &cliffordSim)
{
    const double probThreshold = 1E-10;
    size_t nrQubits = qubitRegister.getNrQubits();

    for (size_t q = 0; q < nrQubits; ++q)
    {
        const double prob1 = cliffordSim.GetQubitProbability(q);
        const double prob2 = qubitRegister.GetQubitProbability(q);

        if (std::abs(prob1 - prob2) > probThreshold)
        {
            std::cout << "\nFailed qubits probabilities" << std::endl;
            std::cout << "Probability statevector: " << prob2 << ", Probability stabilizer: " << prob1 << std::endl;
            return false;
        }
    }

    return true;
}

bool CheckAllStatesProbability(QC::QubitRegister<> &qubitRegister, QC::Clifford::StabilizerSimulator &cliffordSim)
{
    const double probThreshold = 1E-10;
    size_t nrQubits = qubitRegister.getNrQubits();
    const size_t nrStates = 1ULL << nrQubits;

    for (size_t state = 0; state < nrStates; ++state)
    {
        const double prob1 = cliffordSim.getBasisStateProbability(state);
        const double prob2 = qubitRegister.getBasisStateProbability(state);

        if (std::abs(prob1 - prob2) > probThreshold)
        {
            std::cout << "\nFailed states probabilities" << std::endl;
            std::cout << "Probability statevector: " << prob2 << ", Probability stabilizer: " << prob1
                      << ", State: " << state << std::endl;
            return false;
        }
    }

    return true;
}

bool CheckMeasurements(QC::Clifford::StabilizerSimulator &cliffordSim)
{
    size_t nrQubits = cliffordSim.getNrQubits();
    for (size_t q = 0; q < nrQubits; ++q)
    {
        const bool res1 = cliffordSim.MeasureQubit(q);
        const bool res2 = cliffordSim.MeasureQubit(q);
        if (res1 != res2)
        {
            std::cout << "\nFailed qubits measurements" << std::endl;
            std::cout << "Measurement 1: " << res2 << ", Measurement 2: " << res1 << ", for qubit: " << q << std::endl;
            return false;
        }
    }
    return true;
}

namespace
{
class CliffordRegressionSimulator : public QC::Clifford::StabilizerSimulator
{
  public:
    using StabilizerSimulator::StabilizerSimulator;

    std::mt19937_64 RandomEngine() const
    {
        return gen;
    }

    bool DistributionIsPrepared() const
    {
        return (validDistributions & FullDistribution) != 0;
    }

    bool SameTableau(const CliffordRegressionSimulator &other) const
    {
        if (getNrQubits() != other.getNrQubits())
            return false;
        for (size_t q = 0; q < getNrQubits(); ++q)
        {
            if (!(inverseZ[q] == other.inverseZ[q]) || inverseZ[q].PhaseSign != other.inverseZ[q].PhaseSign ||
                !(inverseX[q] == other.inverseX[q]) || inverseX[q].PhaseSign != other.inverseX[q].PhaseSign)
                return false;
        }
        return true;
    }
};

bool CliffordCheck(bool condition, const char *message)
{
    if (!condition)
        std::cout << "Clifford regression failed: " << message << std::endl;
    return condition;
}

template <class Exception, class F> bool CliffordRejects(F &&f)
{
    try
    {
        f();
    }
    catch (const Exception &)
    {
        return true;
    }
    catch (...)
    {
        return false;
    }
    return false;
}

bool CliffordPackedRegression()
{
    using Row = QC::Clifford::detail::TableauRow<false>;
    using ConstRow = QC::Clifford::detail::TableauRow<true>;
    static_assert(std::is_copy_constructible_v<Row> && std::is_copy_constructible_v<ConstRow>);
    static_assert(!std::is_copy_assignable_v<Row> && !std::is_move_assignable_v<Row>);
    static_assert(!std::is_copy_assignable_v<ConstRow> && !std::is_move_assignable_v<ConstRow>);
    QC::Clifford::detail::PackedTableau signedRows(2, 65);
    signedRows[0].X[64] = true;
    signedRows[1].CopyFrom(signedRows[0]);
    signedRows[1].PhaseSign = true;
    if (!CliffordCheck(!(signedRows[0] == signedRows[1]), "Row equality includes sign"))
        return false;
    signedRows[0].CopyFrom(signedRows[1]);
    if (!CliffordCheck(signedRows[0] == signedRows[1], "Explicit row copy includes sign"))
        return false;
    std::mt19937 rng(9173);
    for (size_t n : {0, 1, 2, 63, 64, 65, 127, 128, 129, 1025})
        for (unsigned repeat = 0; repeat < 128; ++repeat)
        {
            QC::Clifford::detail::PackedTableau rows(2, n);
            auto left = rows[0], right = rows[1];
            int phase = 0;
            std::vector<bool> expectedX(n), expectedZ(n);
            for (size_t q = 0; q < n; ++q)
            {
                const int x = rng() & 1, z = rng() & 1, rx = rng() & 1, rz = rng() & 1;
                left.X[q] = x;
                left.Z[q] = z;
                right.X[q] = rx;
                right.Z[q] = rz;
                expectedX[q] = x != rx;
                expectedZ[q] = z != rz;
                phase += QC::PauliStringXZWithSign::g(x, z, rx, rz);
            }
            left.PhaseSign = bool(rng() & 1);
            right.PhaseSign = bool(rng() & 1);
            phase += 2 * int(bool(left.PhaseSign) != bool(right.PhaseSign));
            const unsigned extra = (phase & 1) ? 1 : 0;
            phase += extra;
            left.Multiply(right, extra);
            if (!CliffordCheck(bool(left.PhaseSign) == ((phase & 2) != 0), "Packed product phase"))
                return false;
            for (size_t q = 0; q < n; ++q)
                if (!CliffordCheck(left.X[q] == expectedX[q] && left.Z[q] == expectedZ[q], "Packed product bits"))
                    return false;
            left.Multiply(left);
            if (!CliffordCheck(!left.HasX() && !left.PhaseSign, "Aliased packed square"))
                return false;
            for (size_t w = 0; w < left.Words(); ++w)
                if (!CliffordCheck(left.Z.words[w] == 0, "Packed square is identity"))
                    return false;
        }
    return true;
}

bool CliffordDistributionRegression()
{
    CliffordRegressionSimulator cached(4);
    QC::QubitRegister<> cachedRef(4);
    cached.ApplyH(0);
    cached.ApplyCX(1, 0);
    cached.ApplyH(2);
    cachedRef.ApplyGate(QC::Gates::HadamardGate<>(), 0);
    cachedRef.ApplyGate(QC::Gates::CNOTGate<>(), 1, 0);
    cachedRef.ApplyGate(QC::Gates::HadamardGate<>(), 2);
    for (int gate : {1, 2, 5, 11, 3, 4})
        for (size_t q = 0; q < 4; ++q)
        {
            cached.AllProbabilities();
            ApplyGate(cached, gate, int(q), int((q + 1) % 4));
            cachedRef.ApplyGate(*GetGate(gate), q, (q + 1) % 4);
            if (!CliffordCheck(cached.DistributionIsPrepared() && CheckAllStatesProbability(cachedRef, cached),
                               "Diagonal and Pauli gates preserve or patch prepared probabilities"))
                return false;
            CliffordRegressionSimulator cold(cached);
            cached.SetSeed(781);
            cold.SetSeed(781);
            if (!CliffordCheck(cached.SampleBasisStates(16) == cold.SampleBasisStates(16),
                               "Patched cache matches a fresh distribution with the same seed"))
                return false;
        }
    std::mt19937 rng(78381);
    for (size_t n = 2; n <= 7; ++n)
    {
        CliffordRegressionSimulator sim(n);
        QC::QubitRegister<> ref(n);
        sim.SetMultithreading(false);
        ref.SetMultithreading(false);
        for (size_t step = 0; step < 90; ++step)
        {
            const size_t a = rng() % n, b = (a + 1 + rng() % (n - 1)) % n;
            const int gate = rng() % 15;
            ApplyGate(sim, gate, int(a), int(b));
            ref.ApplyGate(*GetGate(gate), a, b);
            const auto probabilities = sim.AllProbabilities();
            for (size_t state = 0; state < probabilities.size(); ++state)
                if (!CliffordCheck(approxEqual(probabilities[state], ref.getBasisStateProbability(state)),
                                   "Affine distribution vs state vector"))
                    return false;
            const CliffordRegressionSimulator before(sim);
            sim.SetSeed(77);
            const auto shots = sim.SampleBasisStates(32);
            sim.SetSeed(77);
            if (!CliffordCheck(shots == sim.SampleBasisStates(32) && sim.SameTableau(before),
                               "Terminal sampling preserves state and seeded repeatability"))
                return false;
            for (const auto &shot : shots)
            {
                size_t state = 0;
                for (size_t q = 0; q < n; ++q)
                    if (shot[q])
                        state |= size_t(1) << q;
                if (!CliffordCheck(probabilities[state] > 0.0, "Samples stay in affine support"))
                    return false;
            }
        }
    }
    CliffordRegressionSimulator bell(65);
    bell.ApplyH(0);
    for (size_t q = 1; q < 65; ++q)
        bell.ApplyCX(q, 0);
    bell.SetSeed(71);
    size_t ones = 0;
    for (const auto &shot : bell.SampleBasisStates(2048))
    {
        ones += shot[0];
        for (size_t q = 1; q < 65; ++q)
            if (!CliffordCheck(shot[q] == shot[0], "Sampling preserves wide GHZ correlations"))
                return false;
    }
    return CliffordCheck(ones > 850 && ones < 1200, "Terminal sampler frequencies");
}

bool CliffordGenericGateRegression()
{
    const auto check = [](const QC::Gates::QuantumGateWithOp<> &gate) {
        CliffordRegressionSimulator actual(3);
        QC::QubitRegister<> reference(3);
        actual.ApplyH(0);
        actual.ApplyS(0);
        actual.ApplyCX(2, 0);
        actual.ApplyH(1);
        reference.ApplyGate(QC::Gates::HadamardGate<>(), 0);
        reference.ApplyGate(QC::Gates::SGate<>(), 0);
        reference.ApplyGate(QC::Gates::CNOTGate<>(), 2, 0);
        reference.ApplyGate(QC::Gates::HadamardGate<>(), 1);
        actual.GetQubitProbability(2);
        const QC::Gates::AppliedGate<> applied(gate.getRawOperatorMatrix(), 2, 0, 1);
        actual.ApplyGate(applied);
        reference.ApplyGate(gate, 2, 0, 1);
        const auto &state = reference.getRegisterStorage();
        for (unsigned code = 0; code < 64; ++code)
        {
            std::string pauli(3, 'I');
            const char axes[] = {'I', 'X', 'Y', 'Z'};
            for (unsigned q = 0; q < 3; ++q)
                pauli[q] = axes[(code >> (2 * q)) & 3U];
            std::complex<double> expected = 0.;
            for (unsigned input = 0; input < 8; ++input)
            {
                unsigned output = input;
                std::complex<double> coefficient = 1.;
                for (unsigned q = 0; q < 3; ++q)
                {
                    const bool bit = (input >> q) & 1U;
                    if (pauli[q] == 'X' || pauli[q] == 'Y')
                        output ^= 1U << q;
                    if (pauli[q] == 'Y')
                        coefficient *= std::complex<double>(0., bit ? -1. : 1.);
                    else if (pauli[q] == 'Z' && bit)
                        coefficient = -coefficient;
                }
                expected += std::conj(state[output]) * coefficient * state[input];
            }
            if (!CliffordCheck(std::abs(expected.real() - actual.ExpectationValue(pauli)) < 3e-14,
                               "Generic Clifford gate matches statevector Pauli expectations"))
                return false;
        }
        return CliffordCheck(std::abs(actual.GetQubitProbability(2) - reference.GetQubitProbability(2)) < 3e-14,
                             "Generic Clifford gate refreshes cached probabilities");
    };
    for (int code = 0; code <= 14; ++code)
        if (!check(*GetGate(code)))
            return false;
    if (!check(QC::Gates::ECRGate<>()) || !check(QC::Gates::DCXGate<>()) || !check(QC::Gates::SplitterGate<>()))
        return false;
    for (const auto &matrix : {QC::Gates::DCXGate<>().getRawOperatorMatrix().adjoint().eval(),
                               QC::Gates::SplitterGate<>().getRawOperatorMatrix().adjoint().eval()})
        if (!check(QC::Gates::AppliedGate<>(matrix, 0, 1, 2)))
            return false;

    // All 24 projective one-qubit Cliffords, encoded without a sequence of
    // rounded Hadamard products: eight monomial and sixteen full matrices.
    const auto phase = [](unsigned power) { return QC::detail::LocalPauliSum::IPower(power); };
    for (unsigned a = 0; a < 4; ++a)
    {
        for (unsigned flip = 0; flip < 2; ++flip)
        {
            Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(2, 2);
            matrix(flip, 0) = 1.;
            matrix(1 ^ flip, 1) = phase(a);
            if (!check(QC::Gates::AppliedGate<>(matrix, 0)))
                return false;
        }
        for (unsigned b = 0; b < 4; ++b)
        {
            Eigen::MatrixXcd matrix = QC::Gates::HadamardGate<>().getRawOperatorMatrix();
            for (unsigned row = 0; row < 2; ++row)
                for (unsigned col = 0; col < 2; ++col)
                    matrix(row, col) = QC::detail::CliffordPhase(matrix(row, col), a * row + b * col);
            if (!check(QC::Gates::AppliedGate<>(matrix, 0)))
                return false;
        }
    }

    // Exercise every generator image independently, rather than only one
    // initial stabilizer state. ConjugatePauli returns U^dagger P U.
    const auto checkImages = [](const Eigen::MatrixXcd &matrix,
                                std::shared_ptr<const QC::detail::LocalGate> compiled = {}) {
        if (!compiled)
            compiled = QC::detail::FindCanonicalCliffordOwned(matrix);
        if (!CliffordCheck(compiled != nullptr, "Canonical local Clifford is recognized"))
            return false;
        const unsigned dimension = static_cast<unsigned>(matrix.rows());
        const auto pauliMatrix = [dimension](unsigned x, unsigned z) {
            Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(dimension, dimension);
            const auto coefficient = QC::detail::LocalPauliSum::IPower(QC::Clifford::detail::Popcount(x & z));
            for (unsigned col = 0; col < dimension; ++col)
                result(col ^ x, col) = (QC::Clifford::detail::Popcount(col & z) & 1U) ? -coefficient : coefficient;
            return result;
        };
        for (unsigned x = 0; x < dimension; ++x)
            for (unsigned z = 0; z < dimension; ++z)
            {
                bool negative;
                const auto image = compiled->ConjugatePauli({uint8_t(x), uint8_t(z)}, negative);
                const Eigen::MatrixXcd expected = matrix.adjoint() * pauliMatrix(x, z) * matrix;
                const Eigen::MatrixXcd actual = (negative ? -1. : 1.) * pauliMatrix(image.x, image.z);
                if (!CliffordCheck((expected - actual).cwiseAbs().maxCoeff() < 2e-15,
                                   "Synthesized Clifford matches every dense Pauli image"))
                    return false;
            }
        return true;
    };
    // The fixed gates retain their original 1/sqrt(2) bits. Quarter-angle
    // constructors use sqrt(.5); recognition accepts either complete encoding.
    constexpr double quarter = 1.57079632679489661923;
    constexpr double legacyAmplitude = 0x1.6a09e667f3bccp-1;
    constexpr double constructorAmplitude = 0x1.6a09e667f3bcdp-1;
    std::array<Eigen::MatrixXcd, 3> legacyMatrices{Eigen::MatrixXcd(2, 2), Eigen::MatrixXcd(2, 2),
                                                   Eigen::MatrixXcd(2, 2)};
    legacyMatrices[0] << legacyAmplitude, legacyAmplitude, legacyAmplitude, -legacyAmplitude;
    legacyMatrices[1] << legacyAmplitude, std::complex<double>(0., -legacyAmplitude),
        std::complex<double>(0., legacyAmplitude), -legacyAmplitude;
    legacyMatrices[2] << legacyAmplitude, std::complex<double>(0., legacyAmplitude),
        std::complex<double>(0., legacyAmplitude), legacyAmplitude;
    if (!CliffordCheck(legacyMatrices[0] == QC::Gates::HadamardGate<>().getRawOperatorMatrix() &&
                           legacyMatrices[1] == QC::Gates::HyGate<>().getRawOperatorMatrix() &&
                           legacyMatrices[2] == QC::Gates::SplitterGate<>().getRawOperatorMatrix(),
                       "Fixed gate constructor coefficients retain their legacy binary encodings"))
        return false;
    for (unsigned qubits = 1; qubits <= 3; ++qubits)
        for (unsigned q = 0; q < qubits; ++q)
        {
            for (const auto &legacy : legacyMatrices)
            {
                const QC::Gates::AppliedGate<> raw(legacy, 0);
                const auto embedded =
                    static_cast<const QC::Gates::QuantumGateWithOp<> &>(raw).getOperatorMatrix(qubits, q).eval();
                if (!checkImages(embedded) || !check(QC::Gates::AppliedGate<>(embedded, 0, 1, 2)))
                    return false;
            }
            for (const auto &embedded : {QC::Gates::RxGate<>(quarter).getOperatorMatrix(qubits, q).eval(),
                                         QC::Gates::RyGate<>(quarter).getOperatorMatrix(qubits, q).eval()})
                if (!checkImages(embedded) || !check(QC::Gates::AppliedGate<>(embedded, 0, 1, 2)))
                    return false;
        }
    for (int theta : {-1, 1})
        for (unsigned phi = 0; phi < 4; ++phi)
            for (unsigned lambda = 0; lambda < 4; ++lambda)
                for (unsigned gamma = 0; gamma < 4; ++gamma)
                {
                    const QC::Gates::UGate<> gate(theta * quarter, phi * quarter, lambda * quarter, gamma * quarter);
                    if (!checkImages(gate.getRawOperatorMatrix()) || !check(gate))
                        return false;
                }
    for (unsigned first = 0; first < 3; ++first)
        for (unsigned second = 0; second < 3; ++second)
            if (first != second)
            {
                const auto embedded = QC::Gates::ECRGate<>().getOperatorMatrix(3, first, second).eval();
                if (!checkImages(embedded) || !check(QC::Gates::AppliedGate<>(embedded, 0, 1, 2)))
                    return false;
            }
    const Eigen::MatrixXcd swap = QC::Gates::SwapGate<>().getRawOperatorMatrix();
    const Eigen::MatrixXcd ecr = QC::Gates::ECRGate<>().getRawOperatorMatrix();
    const auto tensor = [](const Eigen::MatrixXcd &low, const Eigen::MatrixXcd &high) {
        Eigen::MatrixXcd result(low.rows() * high.rows(), low.cols() * high.cols());
        for (Eigen::Index row = 0; row < result.rows(); ++row)
            for (Eigen::Index col = 0; col < result.cols(); ++col)
                result(row, col) = low(row % low.rows(), col % low.cols()) * high(row / low.rows(), col / low.cols());
        return result;
    };
    for (const auto &matrix : {(swap * ecr).eval(), (swap * ecr * swap).eval(),
                               tensor(QC::Gates::RxGate<>(quarter).getRawOperatorMatrix(),
                                      QC::Gates::SquareRootNOTGate<>().getRawOperatorMatrix()),
                               tensor(ecr, QC::Gates::SquareRootNOTGate<>().getRawOperatorMatrix())})
        if (!checkImages(matrix) || !check(QC::Gates::AppliedGate<>(matrix, 0, 1, 2)))
            return false;
    QC::ExtendedStabilizer fromU(8), nativeH(8), typedH(8);
    for (unsigned q = 0; q < 8; ++q)
    {
        fromU.ApplyGate(QC::Gates::UGate<>(quarter, 0., 2 * quarter), q);
        nativeH.ApplyH(q);
        typedH.ApplyGate(QC::Gates::HadamardGate<>(), q);
        if (!CliffordCheck(fromU.GetFrames().front().GetFrameSize() == 1 &&
                               nativeH.GetFrames().front().GetFrameSize() == 1 &&
                               typedH.GetFrames().front().GetFrameSize() == 1,
                           "Quarter-angle U gates preserve a single Extended stabilizer component"))
            return false;
    }
    for (unsigned basis = 0; basis < 256; ++basis)
        if (!CliffordCheck(
                std::abs(fromU.getBasisStateProbability(basis) - nativeH.getBasisStateProbability(basis)) < 1e-15 &&
                    std::abs(fromU.getBasisStateProbability(basis) - typedH.getBasisStateProbability(basis)) < 1e-15,
                "Eight disjoint U-encoded Hadamards match native and typed Hadamards"))
            return false;
    for (unsigned q = 0; q < 3; ++q)
    {
        const Eigen::MatrixXcd embedded = QC::Gates::HadamardGate<>().getOperatorMatrix(3, q);
        if (!check(QC::Gates::AppliedGate<>(embedded, 0, 1, 2)) || !checkImages(embedded))
            return false;
    }
    // SX/SXdg use the exact dyadic eighth-phase family. Embedding them or
    // composing with signed permutations must not force a generic expansion.
    for (unsigned qubits = 1; qubits <= 3; ++qubits)
        for (unsigned q = 0; q < qubits; ++q)
            for (bool dagger : {false, true})
            {
                const Eigen::MatrixXcd embedded = dagger
                                                      ? QC::Gates::SquareRootNOTDagGate<>().getOperatorMatrix(qubits, q)
                                                      : QC::Gates::SquareRootNOTGate<>().getOperatorMatrix(qubits, q);
                std::vector<Eigen::MatrixXcd> matrices{
                    embedded, QC::Gates::PauliZGate<>().getOperatorMatrix(qubits, q) * embedded,
                    QC::Gates::SGate<>().getOperatorMatrix(qubits, q) * embedded};
                if (qubits > 1)
                    matrices.push_back(QC::Gates::CNOTGate<>().getOperatorMatrix(qubits, (q + 1) % qubits, q) *
                                       embedded);
                for (const auto &matrix : matrices)
                    if (!checkImages(matrix) || !check(QC::Gates::AppliedGate<>(matrix, 0, 1, 2)))
                        return false;
            }
    std::array<unsigned, 3> permutation{0, 1, 2};
    do
    {
        Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(8, 8);
        for (unsigned input = 0; input < 8; ++input)
        {
            unsigned output = 0;
            for (unsigned q = 0; q < 3; ++q)
                output |= ((input >> q) & 1U) << permutation[q];
            matrix(output, input) = 1.;
        }
        if (!check(QC::Gates::AppliedGate<>(matrix, 0, 1, 2)) || !checkImages(matrix))
            return false;
    } while (std::next_permutation(permutation.begin(), permutation.end()));

    // Canonical H layers sandwiched by independently generated signed affine
    // permutations cover all support sizes and nonlocal generator images.
    const Eigen::MatrixXcd heldMatrix = QC::Gates::HadamardGate<>().getOperatorMatrix(3, 1);
    const auto heldCompilation = QC::detail::FindCanonicalCliffordOwned(heldMatrix);
    std::mt19937 random(7751);
    for (unsigned trial = 0; trial < 96; ++trial)
    {
        const unsigned hMask = trial & 7U;
        const double amplitude = 1. / std::sqrt(double(1U << QC::Clifford::detail::Popcount(hMask)));
        Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(8, 8);
        for (unsigned row = 0; row < 8; ++row)
            for (unsigned col = 0; col < 8; ++col)
                if (!((row ^ col) & ~hMask))
                    matrix(row, col) =
                        (QC::Clifford::detail::Popcount(row & col & hMask) & 1U) ? -amplitude : amplitude;
        for (unsigned side = 0; side < 2; ++side)
            for (unsigned step = 0; step < 12; ++step)
            {
                const unsigned target = random() % 3, control = (target + 1 + random() % 2) % 3;
                if (random() & 1U)
                {
                    for (unsigned index = 0; index < 8; ++index)
                        if (((index >> control) & 1U) && !((index >> target) & 1U))
                        {
                            if (side)
                                matrix.col(index).swap(matrix.col(index ^ (1U << target)));
                            else
                                matrix.row(index).swap(matrix.row(index ^ (1U << target)));
                        }
                }
                else
                    for (unsigned row = 0; row < 8; ++row)
                        for (unsigned col = 0; col < 8; ++col)
                            if (((side ? col : row) >> target) & 1U)
                                matrix(row, col) = QC::detail::CliffordPhase(matrix(row, col), 1);
            }
        if (!checkImages(matrix) || !check(QC::Gates::AppliedGate<>(matrix, 0, 1, 2)))
            return false;
        if (trial < 8)
        {
            // Cover support sizes 1, 2, 4 and 8 with all four balanced global
            // phases, constructing their canonical entries directly.
            const double balanced = 1. / std::sqrt(double(2U << QC::Clifford::detail::Popcount(hMask)));
            Eigen::MatrixXcd eighth = Eigen::MatrixXcd::Zero(8, 8);
            for (unsigned row = 0; row < 8; ++row)
                for (unsigned col = 0; col < 8; ++col)
                {
                    const auto value = matrix(row, col);
                    if (value.real() != 0.)
                        eighth(row, col) = std::copysign(1., value.real()) * std::complex<double>(balanced, balanced);
                    else if (value.imag() != 0.)
                        eighth(row, col) = std::copysign(1., value.imag()) * std::complex<double>(-balanced, balanced);
                }
            for (unsigned turn = 0; turn < 4; ++turn)
            {
                if (!checkImages(eighth) || !check(QC::Gates::AppliedGate<>(eighth, 0, 1, 2)))
                    return false;
                for (Eigen::Index index = 0; index < eighth.size(); ++index)
                    eighth.data()[index] = QC::detail::CliffordPhase(eighth.data()[index], 1);
            }
        }
    }
    if (!checkImages(heldMatrix, heldCompilation))
        return false;
    // Retained compilations outlive both bounded-cache eviction and the
    // producing thread's thread_local caches.
    const auto fromOtherThread = std::async(std::launch::async, [heldMatrix] {
                                     return QC::detail::FindCanonicalCliffordOwned(heldMatrix);
                                 }).get();
    if (!checkImages(heldMatrix, fromOtherThread))
        return false;
    for (int turn = -4; turn <= 4; ++turn)
    {
        const double angle = turn * quarter;
        if (!check(QC::Gates::RxGate<>(angle)) || !check(QC::Gates::RyGate<>(angle)) ||
            !check(QC::Gates::RzGate<>(angle)) || !check(QC::Gates::PhaseShiftGate<>(angle)) ||
            !check(QC::Gates::RXXGate<>(angle)) || !check(QC::Gates::RYYGate<>(angle)) ||
            !check(QC::Gates::RZZGate<>(angle)) || !check(QC::Gates::RZXGate<>(angle)) ||
            !check(QC::Gates::RXYGate<>(angle)))
            return false;
    }
    // Constructor matrices reduce quarter turns modulo eight, including the
    // largest odd quarter turn supported by the exact-angle helper.
    for (int turn : {-1048575, -11, -9, 9, 11, 1048575})
    {
        const double angle = turn * quarter;
        if (!check(QC::Gates::RxGate<>(angle)) || !check(QC::Gates::RyGate<>(angle)) ||
            !check(QC::Gates::RzGate<>(angle)) || !check(QC::Gates::RXXGate<>(angle)) ||
            !check(QC::Gates::RYYGate<>(angle)) || !check(QC::Gates::RZZGate<>(angle)) ||
            !check(QC::Gates::RZXGate<>(angle)) || !check(QC::Gates::RXYGate<>(angle)))
            return false;
    }
    CliffordRegressionSimulator actual(3);
    actual.ApplyH(0);
    const CliffordRegressionSimulator before(actual);
    // Either complete encoding is valid, but no entry may select a different
    // family. Neighbors outside the two accepted values and tiny new entries
    // must retain the generic path and be rejected by the Clifford-only API.
    const Eigen::MatrixXcd highH = QC::Gates::UGate<>(quarter, 0., 2 * quarter).getRawOperatorMatrix();
    const Eigen::MatrixXcd lowH = legacyMatrices[0];
    for (bool balanced : {false, true})
    {
        Eigen::MatrixXcd canonical =
            balanced ? tensor(highH, QC::Gates::SquareRootNOTGate<>().getRawOperatorMatrix()) : highH;
        const double low = legacyAmplitude * (balanced ? .5 : 1.);
        const double high = constructorAmplitude * (balanced ? .5 : 1.);
        for (unsigned variant = 0; variant < 7; ++variant)
        {
            Eigen::MatrixXcd perturbed = canonical;
            const double magnitude = variant == 0   ? std::nextafter(low, 0.)
                                     : variant == 1 ? std::nextafter(high, 1.)
                                     : variant == 2 ? low
                                                    : high;
            for (Eigen::Index index = 0; index < perturbed.size(); ++index)
            {
                const auto value = perturbed.data()[index];
                perturbed.data()[index] = {value.real() == 0. ? 0. : std::copysign(magnitude, value.real()),
                                           value.imag() == 0. ? 0. : std::copysign(magnitude, value.imag())};
            }
            if (variant == 2 || variant == 3)
            {
                const double opposite = variant == 2 ? high : low;
                const auto first = perturbed(0, 0);
                perturbed(0, 0) = {std::copysign(opposite, first.real()),
                                   balanced ? std::copysign(opposite, first.imag()) : 0.};
            }
            else if (variant == 4)
                perturbed(0, 0) = {high, balanced ? low : 1e-310};
            else if (variant == 5)
                perturbed *= 1. + 1e-12;
            else if (variant == 6)
                perturbed(0, 0) = -perturbed(0, 0);
            if (!CliffordCheck(QC::detail::FindCanonicalClifford(perturbed) == nullptr &&
                                   CliffordRejects<std::invalid_argument>(
                                       [&] { actual.ApplyGate(QC::Gates::AppliedGate<>(perturbed, 0, 1, 2)); }) &&
                                   actual.SameTableau(before) && checkImages(canonical),
                               "Exact recognition rejects mixed encodings, outer neighbors and invalid phases"))
                return false;
        }
    }
    // A scaled raw matrix that is bit-identical to a constructor encoding is
    // indistinguishable from that encoding. This does not accept other rounded
    // products: even products of canonical H matrices can be noncanonical.
    if (!CliffordCheck(lowH * std::nextafter(1., 2.) == highH && checkImages(highH),
                       "Bit-identical scaled aliases share the accepted constructor encoding"))
        return false;
    for (const auto &matrix : {(lowH * lowH).eval(), (highH * highH).eval(), tensor(lowH, lowH), tensor(highH, highH),
                               tensor(tensor(lowH, lowH), lowH), tensor(tensor(highH, highH), highH)})
        if (!CliffordCheck(QC::detail::FindCanonicalClifford(matrix) == nullptr,
                           "Arbitrary rounded products are not snapped to canonical Clifford matrices"))
            return false;
    Eigen::MatrixXcd altered = QC::Gates::HadamardGate<>().getRawOperatorMatrix();
    altered(0, 0) = std::nextafter(altered(0, 0).real(), 1.);
    if (!CliffordCheck(
            CliffordRejects<std::invalid_argument>([&] { actual.ApplyGate(QC::Gates::TGate<>(), 0); }) &&
                CliffordRejects<std::invalid_argument>(
                    [&] { actual.ApplyGate(QC::Gates::RxGate<>(quarter + 1e-15), 0); }) &&
                CliffordRejects<std::invalid_argument>(
                    [&] { actual.ApplyGate(QC::Gates::AppliedGate<>(altered, 0)); }) &&
                CliffordRejects<std::invalid_argument>([&] { actual.ApplyGate(QC::Gates::CNOTGate<>(), 0, 0); }) &&
                CliffordRejects<std::invalid_argument>(
                    [&] { actual.ApplyGate(QC::Gates::ToffoliGate<>(), 0, 1, 2); }) &&
                CliffordRejects<std::invalid_argument>(
                    [&] { actual.ApplyGate(QC::Gates::FredkinGate<>(), 0, 1, 2); }) &&
                actual.SameTableau(before),
            "Generic Clifford gate rejects non-Clifford and invalid inputs before mutation"))
        return false;
    Eigen::MatrixXcd embedded = QC::Gates::HadamardGate<>().getOperatorMatrix(3, 1);
    const auto original = QC::detail::FindCanonicalCliffordOwned(embedded);
    for (unsigned variant = 0; variant < 3; ++variant)
    {
        Eigen::MatrixXcd perturbed = embedded;
        if (variant == 0)
            perturbed(0, 0) = std::nextafter(perturbed(0, 0).real(), 1.);
        else if (variant == 1)
            perturbed(0, 1) = 1e-310;
        else
            perturbed *= 1. + 1e-12;
        if (!CliffordCheck(QC::detail::FindCanonicalClifford(perturbed) == nullptr &&
                               CliffordRejects<std::invalid_argument>(
                                   [&] { actual.ApplyGate(QC::Gates::AppliedGate<>(perturbed, 0, 1, 2)); }) &&
                               actual.SameTableau(before) && QC::detail::FindCanonicalClifford(embedded) != nullptr &&
                               checkImages(embedded, original),
                           "Exact synthesis rejects perturbations and scaled filters without cache poisoning"))
            return false;
    }
    const Eigen::MatrixXcd dyadic = QC::Gates::SquareRootNOTGate<>().getOperatorMatrix(3, 0);
    for (unsigned variant = 0; variant < 10; ++variant)
    {
        Eigen::MatrixXcd perturbed = dyadic;
        switch (variant)
        {
        case 0:
            perturbed *= std::nextafter(1., 2.);
            break;
        case 1:
            perturbed(0, 0) = {std::nextafter(.5, 1.), .5};
            break;
        case 2:
            perturbed(0, 0) = {.5, std::nextafter(.5, 1.)};
            break;
        case 3:
            perturbed(0, 2) = 1e-310;
            break;
        case 4:
            perturbed(0, 0) = 1. / std::sqrt(2.);
            break;
        case 5:
            perturbed(0, 1) = 1. / std::sqrt(2.);
            break;
        case 6:
            perturbed *= 1. + 1e-12;
            break;
        case 7:
            perturbed *= std::sqrt(2.);
            break;
        case 8:
            perturbed = std::complex<double>(.5, .5) * Eigen::MatrixXcd::Identity(8, 8);
            break;
        default:
            perturbed(0, 1) = 0.;
            break;
        }
        if (!CliffordCheck(QC::detail::FindCanonicalClifford(perturbed) == nullptr &&
                               CliffordRejects<std::invalid_argument>(
                                   [&] { actual.ApplyGate(QC::Gates::AppliedGate<>(perturbed, 0, 1, 2)); }) &&
                               actual.SameTableau(before) && checkImages(dyadic),
                           "Eighth-phase recognition preserves global normalization and exact coefficients"))
            return false;
    }
    return true;
}

bool CliffordSingleRotationRegression()
{
    constexpr double quarter = 1.57079632679489661923;
    using Simulator = QC::Clifford::StabilizerSimulator;
    using ApplyRotation = void (Simulator::*)(size_t, double);
    const ApplyRotation rotations[] = {&Simulator::ApplyRx, &Simulator::ApplyRy, &Simulator::ApplyRz,
                                       &Simulator::ApplyPhase};
    for (unsigned axis = 0; axis < 4; ++axis)
        for (int turns = -8; turns <= 8; ++turns)
            for (unsigned input = 0; input < 6; ++input)
            {
                CliffordRegressionSimulator actual(1);
                Eigen::Vector2cd state = Eigen::Vector2cd::Zero();
                state[0] = 1.;
                if (input & 1U)
                {
                    actual.ApplyX(0);
                    state = QC::Gates::PauliXGate<>().getRawOperatorMatrix() * state;
                }
                if (input >= 2)
                {
                    actual.ApplyH(0);
                    state = QC::Gates::HadamardGate<>().getRawOperatorMatrix() * state;
                }
                if (input >= 4)
                {
                    actual.ApplyS(0);
                    state = QC::Gates::SGate<>().getRawOperatorMatrix() * state;
                }
                // Populate the probability cache before both diagonal and
                // non-diagonal updates, including the Pauli bit-flip cases.
                actual.GetQubitProbability(0);
                const double angle = turns * quarter;
                Eigen::Matrix2cd matrix;
                if (axis == 0)
                    matrix = QC::Gates::RxGate<>(angle).getRawOperatorMatrix();
                else if (axis == 1)
                    matrix = QC::Gates::RyGate<>(angle).getRawOperatorMatrix();
                else if (axis == 2)
                    matrix = QC::Gates::RzGate<>(angle).getRawOperatorMatrix();
                else
                    matrix = QC::Gates::PhaseShiftGate<>(angle).getRawOperatorMatrix();
                state = matrix * state;
                (actual.*rotations[axis])(0, angle);
                const auto coherence = std::conj(state[0]) * state[1];
                if (!CliffordCheck(std::abs(actual.ExpectationValue("X") - 2. * coherence.real()) < 3e-14 &&
                                       std::abs(actual.ExpectationValue("Y") - 2. * coherence.imag()) < 3e-14 &&
                                       std::abs(actual.ExpectationValue("Z") -
                                                (std::norm(state[0]) - std::norm(state[1]))) < 3e-14 &&
                                       std::abs(actual.GetQubitProbability(0) - std::norm(state[1])) < 3e-14,
                                   "Single-qubit Clifford rotations and cached probabilities"))
                    return false;
            }
    for (auto apply : rotations)
    {
        CliffordRegressionSimulator actual(1);
        actual.ApplyH(0);
        actual.GetQubitProbability(0);
        const CliffordRegressionSimulator before(actual);
        for (double invalid : {std::nextafter(quarter, 0.), std::nextafter(quarter, 2.),
                               std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
            if (!CliffordCheck(CliffordRejects<std::invalid_argument>([&] { (actual.*apply)(0, invalid); }) &&
                                   actual.SameTableau(before) && actual.GetQubitProbability(0) == .5,
                               "Reject non-Clifford rotations before mutation"))
                return false;
        if (!CliffordCheck(CliffordRejects<std::out_of_range>([&] { (actual.*apply)(1, 0.); }) &&
                               actual.SameTableau(before),
                           "Validate even identity rotations"))
            return false;
    }
    return true;
}

// Compact gates on at most three qubits are recognized through their dense
// matrices; the named gates keep their short circuits; a matrix that is not a
// Clifford is rejected before the tableau or its caches change.
bool CliffordCompactGateRegression()
{
    using Applied = QC::Gates::AppliedGate<>;
    constexpr double quarter = 1.57079632679489661923;
    const std::complex<double> i(0., 1.);
    const Eigen::MatrixXcd x = QC::Gates::PauliXGate<>().getRawOperatorMatrix();
    const auto prepare = [](CliffordRegressionSimulator &sim) {
        sim.ApplyH(0);
        sim.ApplyS(0);
        sim.ApplyCX(2, 0);
        sim.ApplyH(1);
        sim.ApplyCX(3, 1);
        sim.ApplySx(3);
    };
    const auto same = [&](const Applied &compact, const auto &native) {
        CliffordRegressionSimulator actual(4), expected(4);
        prepare(actual);
        prepare(expected);
        // Fill the probability cache first: the gate has to refresh it.
        for (size_t q = 0; q < 4; ++q)
            actual.GetQubitProbability(q);
        actual.ApplyGate(compact);
        native(expected);
        bool ok = actual.SameTableau(expected);
        for (size_t q = 0; q < 4; ++q)
            ok = ok && actual.GetQubitProbability(q) == expected.GetQubitProbability(q);
        return ok;
    };
    if (!CliffordCheck(same(Applied::Diagonal({1., 1., 1., -1.}, {0, 2}), [](auto &s) { s.ApplyCZ(0, 2); }),
                       "Compact CZ diagonal") ||
        !CliffordCheck(same(Applied::Diagonal({1., i}, {1}), [](auto &s) { s.ApplyS(1); }), "Compact S diagonal") ||
        !CliffordCheck(same(Applied::MultiControlled(x, 3, {1}), [](auto &s) { s.ApplyCX(3, 1); }),
                       "Compact X with one control") ||
        !CliffordCheck(same(Applied::MultiControlled(x, 3, {1}, 0),
                            [](auto &s) {
                                s.ApplyX(1);
                                s.ApplyCX(3, 1);
                                s.ApplyX(1);
                            }),
                       "Compact X with a control on zero") ||
        !CliffordCheck(same(Applied::MultiControlled(QC::Gates::HadamardGate<>().getRawOperatorMatrix(), 2, {}),
                            [](auto &s) { s.ApplyH(2); }),
                       "Compact gate without controls") ||
        !CliffordCheck(same(Applied::MultiplexedRotation(QC::Gates::RotationAxis::Z, {quarter, -quarter}, 0, {2}),
                            [=](auto &s) { s.ApplyRZZ(0, 2, quarter); }),
                       "Compact multiplexed Z rotation at quarter turns") ||
        !CliffordCheck(same(Applied::MultiplexedRotation(QC::Gates::RotationAxis::X, {2 * quarter, 0.}, 1, {0}),
                            [](auto &s) {
                                s.ApplyX(0);
                                s.ApplyCX(1, 0);
                                s.ApplyX(0);
                                s.ApplyS(0);
                            }),
                       "Compact multiplexed X rotation by a half turn"))
        return false;

    CliffordRegressionSimulator actual(4);
    prepare(actual);
    for (size_t q = 0; q < 4; ++q)
        actual.GetQubitProbability(q);
    const CliffordRegressionSimulator before(actual);
    const auto rejects = [&](const Applied &gate) {
        return CliffordRejects<std::invalid_argument>([&] { actual.ApplyGate(gate); }) && actual.SameTableau(before) &&
               actual.GetQubitProbability(3) == before.GetQubitProbability(3);
    };
    if (!CliffordCheck(rejects(Applied::MultiControlled(x, 3, {0, 1})), "Compact Toffoli is not a Clifford") ||
        !CliffordCheck(rejects(Applied::Diagonal({1., std::polar(1., .5 * quarter)}, {0})),
                       "Compact T diagonal is not a Clifford") ||
        !CliffordCheck(rejects(Applied::MultiplexedRotation(QC::Gates::RotationAxis::Y, {.3, quarter}, 1, {0})),
                       "Compact multiplexed rotation by a generic angle is not a Clifford") ||
        !CliffordCheck(rejects(Applied::Diagonal(std::vector<std::complex<double>>(16, 1.), {0, 1, 2, 3})),
                       "Compact gates on four qubits are rejected") ||
        !CliffordCheck(rejects(Applied::MultiControlled(x, 3, {0, 1, 2})), "Compact X with three controls is rejected"))
        return false;

    // The named gates keep the catalog's circuits, whatever recognizes them first.
    const auto operations = [](const Eigen::MatrixXcd &matrix) {
        const auto compiled = QC::detail::FindCanonicalCliffordOwned(matrix);
        return compiled ? compiled->Cliffords().size() : size_t(1000);
    };
    if (!CliffordCheck(operations(QC::Gates::CNOTGate<>().getRawOperatorMatrix()) == 1 &&
                           operations(QC::Gates::SquareRootNOTGate<>().getRawOperatorMatrix()) == 1 &&
                           operations(QC::Gates::PauliYGate<>().getRawOperatorMatrix()) == 1 &&
                           operations(QC::Gates::ECRGate<>().getRawOperatorMatrix()) == 4 &&
                           operations(QC::Gates::RxGate<>(quarter).getRawOperatorMatrix()) == 1,
                       "Named Clifford gates keep their short circuits") ||
        !CliffordCheck(QC::detail::FindCanonicalClifford(QC::Gates::TGate<>().getRawOperatorMatrix()) == nullptr &&
                           QC::detail::FindCanonicalClifford(QC::Gates::ControlledRzGate<>(.3).getRawOperatorMatrix()) ==
                               nullptr &&
                           QC::detail::FindCanonicalClifford(QC::Gates::ToffoliGate<>().getRawOperatorMatrix()) ==
                               nullptr &&
                           QC::detail::FindCanonicalClifford(QC::Gates::HadamardGate<>().getOperatorMatrix(3, 1)) !=
                               nullptr,
                       "Recognition order changed which matrices are Cliffords"))
        return false;
    // More matrices than one cache set holds: each lookup still answers for its own matrix.
    std::vector<Eigen::MatrixXcd> matrices;
    for (int turns = -8; turns <= 8; ++turns)
    {
        matrices.push_back(QC::Gates::RxGate<>(turns * quarter).getRawOperatorMatrix());
        matrices.push_back(QC::Gates::RZXGate<>(turns * quarter).getRawOperatorMatrix());
        matrices.push_back(QC::Gates::RzGate<>(turns * quarter + .01).getRawOperatorMatrix());
        matrices.push_back(QC::Gates::RXYGate<>(turns * quarter).getOperatorMatrix(3, 2, 0));
    }
    for (unsigned round = 0; round < 3; ++round)
        for (size_t index = 0; index < matrices.size(); ++index)
        {
            const auto *compiled = QC::detail::FindCanonicalClifford(matrices[index]);
            const bool clifford = index % 4 != 2;
            if (!CliffordCheck((compiled != nullptr) == clifford, "Cached recognition answered for another matrix"))
                return false;
            if (!compiled)
                continue;
            const Eigen::MatrixXcd &matrix = matrices[index];
            const unsigned dimension = static_cast<unsigned>(matrix.rows());
            for (unsigned z = 1; z < dimension; z <<= 1)
            {
                bool negative;
                const auto image = compiled->ConjugatePauli({0, uint8_t(z)}, negative);
                Eigen::MatrixXcd pauli = Eigen::MatrixXcd::Zero(dimension, dimension), expected = pauli;
                for (unsigned col = 0; col < dimension; ++col)
                {
                    pauli(col, col) = (col & z) ? -1. : 1.;
                    const auto coefficient =
                        QC::detail::LocalPauliSum::IPower(QC::Clifford::detail::Popcount(image.x & image.z));
                    expected(col ^ image.x, col) =
                        (negative ? -1. : 1.) *
                        ((QC::Clifford::detail::Popcount(col & image.z) & 1U) ? -coefficient : coefficient);
                }
                if (!CliffordCheck((matrix.adjoint() * pauli * matrix - expected).cwiseAbs().maxCoeff() < 2e-15,
                                   "Cached recognition returned another matrix's circuit"))
                    return false;
            }
        }
    return true;
}

bool CliffordBoundaryRegression()
{
    CliffordRegressionSimulator tiny(2);
    tiny.ApplyH(0);
    if (!CliffordCheck(tiny.ContainsBasisState(1) && tiny.Log2BasisStateProbability(1) == -1.0 &&
                           !tiny.ContainsBasisState(2) && tiny.Log2BasisStateProbability(2) == -INFINITY &&
                           !tiny.ContainsBasisState(4) && tiny.Log2BasisStateProbability(4) == -INFINITY,
                       "Support and logarithmic probability bounds"))
        return false;
    CliffordRegressionSimulator underflow(1076);
    for (size_t q = 0; q < 1074; ++q)
        underflow.ApplyH(q);
    if (!CliffordCheck(underflow.getBasisStateProbability(0) == std::ldexp(1.0, -1074) &&
                           underflow.Log2BasisStateProbability(0) == -1074.0,
                       "Smallest double basis probability"))
        return false;
    underflow.ApplyH(1074);
    std::vector<bool> wideBits(1076);
    if (!CliffordCheck(underflow.getBasisStateProbability(0) == 0.0 && underflow.ContainsBasisState(0) &&
                           underflow.ContainsBasisState(wideBits) &&
                           underflow.Log2BasisStateProbability(wideBits) == -1075.0,
                       "Supported outcomes remain distinguishable after underflow"))
        return false;
    wideBits[1075] = true;
    if (!CliffordCheck(!underflow.ContainsBasisState(wideBits) &&
                           underflow.Log2BasisStateProbability(wideBits) == -INFINITY &&
                           CliffordRejects<std::invalid_argument>([&] { tiny.ContainsBasisState(wideBits); }) &&
                           CliffordRejects<std::invalid_argument>([&] { tiny.Log2BasisStateProbability(wideBits); }),
                       "Unsupported outcomes and malformed bit vectors"))
        return false;
    CliffordRegressionSimulator sim(2), zero(2);
    if (!CliffordCheck(
            CliffordRejects<std::out_of_range>([&] { sim.ApplyH(2); }) &&
                CliffordRejects<std::out_of_range>([&] { sim.ApplySwap(2, 2); }) &&
                CliffordRejects<std::out_of_range>([&] { sim.MeasureQubit(2); }) &&
                CliffordRejects<std::out_of_range>([&] { sim.GetQubitProbability(2); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.ApplyCX(0, 0); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.ApplyCY(0, 0); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.ApplyCZ(0, 0); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.ApplyISwap(0, 0); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.ApplyISwapDag(0, 0); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.ExpectationValue("IIZ"); }) &&
                CliffordRejects<std::invalid_argument>(
                    [&] { sim.getBasisStateProbability(std::vector<bool>{false}); }) &&
                CliffordRejects<std::invalid_argument>([&] { sim.getBasisStateProbability(std::vector<bool>(3)); }) &&
                CliffordRejects<std::runtime_error>([&] { sim.ExpectationValue("Q"); }) && sim.SameTableau(zero),
            "Invalid arguments must be rejected before mutation"))
        return false;
    if (!CliffordCheck(sim.getBasisStateProbability(4) == 0.0 && sim.getBasisStateProbability(0) == 1.0 &&
                           sim.ExpectationValue("Z") == 1.0,
                       "Basis index bounds and short Pauli strings"))
        return false;
    CliffordRegressionSimulator empty(0);
    if (!CliffordCheck(empty.AllProbabilities() == std::vector<double>{1.0} &&
                           empty.getBasisStateProbability(1) == 0.0 && empty.ExpectationValue("") == 1.0,
                       "Zero-qubit probabilities"))
        return false;
    if (sizeof(size_t) == 4)
    {
        CliffordRegressionSimulator wide(32);
        if (!CliffordCheck(CliffordRejects<std::length_error>([&] { wide.AllProbabilities(); }),
                           "Reject overflowing probability count"))
            return false;
    }
    if (!CliffordCheck(CliffordRejects<std::length_error>(
                           [] { QC::Clifford::StabilizerSimulator tooLarge(std::numeric_limits<size_t>::max()); }),
                       "Tableau allocation overflow guard"))
        return false;
    CliffordRegressionSimulator seeded(65);
    seeded.SetSeed(984731);
    seeded.ApplyH(0);
    seeded.SaveState();
    const auto engine = seeded.RandomEngine();
    const CliffordRegressionSimulator snapshot(seeded);
    CliffordRegressionSimulator moved(std::move(seeded));
    if (!CliffordCheck(moved.RandomEngine() == engine && moved.SameTableau(snapshot),
                       "Move constructor preserves RNG and state"))
        return false;
    CliffordRegressionSimulator assigned(2);
    assigned.getBasisStateProbability(0); // Populate a cache of a different size before the move.
    assigned = std::move(moved);
    if (!CliffordCheck(assigned.RandomEngine() == engine && assigned.SameTableau(snapshot),
                       "Move assignment preserves RNG and state"))
        return false;
    assigned.Reset();
    assigned.RestoreState();
    if (!CliffordCheck(assigned.SameTableau(snapshot), "Move preserves saved state"))
        return false;
    QC::PauliStringXZWithSign row(65);
    row.X[64] = true;
    row.PhaseSign = true;
    QC::PauliStringXZWithSign movedRow(std::move(row));
    return CliffordCheck(movedRow.X[64] && movedRow.PhaseSign && row.X.empty() && row.Z.empty(),
                         "Pauli move transfers buffers");
}

bool CliffordCountsRegression()
{
    std::mt19937 actions(77521);
    for (size_t n : {2, 7, 65, 129})
    {
        CliffordRegressionSimulator sim(n);
        sim.SetMultithreading(false);
        sim.SaveState(); // Sampling must not replace this zero-state snapshot.
        for (size_t step = 0; step < 90; ++step)
        {
            const size_t a = actions() % n, b = (a + 1 + actions() % (n - 1)) % n;
            ApplyGate(sim, int(actions() % 15), int(a), int(b));
            if (step % 9 != 0)
                continue;
            for (size_t width : {size_t(1), size_t(3), n + 2})
            {
                std::vector<size_t> qubits(width);
                for (size_t q = 0; q < width; ++q)
                    qubits[q] = (n - 1 - q % n);
                const CliffordRegressionSimulator before(sim);
                CliffordRegressionSimulator ref(sim);
                ref.SaveState();
                ref.SetSeed(12391);
                sim.SetSeed(12391);
                std::unordered_map<std::vector<bool>, size_t> expected;
                std::vector<bool> firstExpected;
                for (size_t shot = 0; shot < 31; ++shot)
                {
                    std::vector<bool> bits(width);
                    for (size_t q = 0; q < width; ++q)
                        bits[q] = ref.MeasureQubit(qubits[q]);
                    if (shot == 0)
                        firstExpected = bits;
                    ++expected[bits];
                    ref.RestoreState();
                }
                const auto actual = sim.SampleCountsMany(qubits, 31);
                if (!CliffordCheck(actual == expected && sim.SameTableau(before) &&
                                       sim.RandomEngine() == ref.RandomEngine(),
                                   "Marginal samples and RNG match sequential measurement"))
                    return false;
                sim.SetSeed(12391);
                if (!CliffordCheck(sim.SampleCountsMany(qubits, 31) == expected,
                                   "Warm marginal cache repeats seeded counts"))
                    return false;
                CliffordRegressionSimulator single(sim);
                single.SetSeed(12391);
                const auto one = single.SampleCountsMany(qubits, 1);
                if (!CliffordCheck(one.size() == 1 && one.begin()->first == firstExpected && one.begin()->second == 1,
                                   "Cold one-shot paths match sequential measurement"))
                    return false;
                if (width <= std::numeric_limits<size_t>::digits)
                {
                    std::unordered_map<size_t, size_t> packed;
                    for (const auto &item : expected)
                    {
                        size_t key = 0;
                        for (size_t q = 0; q < width; ++q)
                            if (item.first[q])
                                key |= size_t(1) << q;
                        packed[key] = item.second;
                    }
                    sim.SetSeed(12391);
                    if (!CliffordCheck(sim.SampleCounts(qubits, 31) == packed,
                                       "Packed counts preserve requested bit order"))
                        return false;
                }
            }
        }
        sim.RestoreState();
        CliffordRegressionSimulator zero(n);
        if (!CliffordCheck(sim.SameTableau(zero), "Counts preserve the caller's saved state"))
            return false;
    }

    CliffordRegressionSimulator sim(129);
    sim.ApplyH(0);
    sim.ApplyCX(128, 0);
    sim.ApplyX(128);
    std::vector<size_t> all(129);
    for (size_t q = 0; q < all.size(); ++q)
        all[q] = 128 - q;
    for (int gate : {1, 2, 5, 11, 3, 4, 0, 9})
    {
        sim.SampleCountsMany(all, 8);
        ApplyGate(sim, gate, 0, 128);
        CliffordRegressionSimulator cold(sim);
        sim.SetSeed(711);
        cold.SetSeed(711);
        if (!CliffordCheck(sim.SampleCountsMany(all, 1) == cold.SampleCountsMany(all, 1) &&
                               sim.RandomEngine() == cold.RandomEngine(),
                           "Cold single-shot fallback matches prepared marginal sampling"))
            return false;
        sim.SetSeed(912);
        cold.SetSeed(912);
        if (!CliffordCheck(sim.SampleCountsMany(all, 17) == cold.SampleCountsMany(all, 17),
                           "Marginal cache updates after diagonal, Pauli and non-diagonal gates"))
            return false;
    }
    // Non-contiguous selection, a high physical index, and repeated outputs.
    const std::vector<size_t> repeated{128, 0, 128, 64};
    sim.SampleCountsMany(repeated, 8);
    sim.ApplyX(128);
    CliffordRegressionSimulator cold(sim);
    sim.SetSeed(918);
    cold.SetSeed(918);
    if (!CliffordCheck(sim.SampleCountsMany(repeated, 64) == cold.SampleCountsMany(repeated, 64),
                       "Pauli cache patch handles repeated selected qubits"))
        return false;
    const auto engine = sim.RandomEngine();
    if (!CliffordCheck(sim.SampleCounts({}, 3).empty() && sim.SampleCountsMany({0}, 0).empty() &&
                           CliffordRejects<std::out_of_range>([&] { sim.SampleCounts({129}, 2); }) &&
                           CliffordRejects<std::out_of_range>([&] { sim.SampleCountsMany({129}, 2); }) &&
                           CliffordRejects<std::invalid_argument>([&] { sim.SampleCounts(all, 2); }) &&
                           sim.RandomEngine() == engine,
                       "Counts validate input before consuming randomness"))
        return false;
    CliffordRegressionSimulator bits(129);
    std::vector<size_t> word(std::numeric_limits<size_t>::digits);
    for (size_t q = 0; q < word.size(); ++q)
    {
        word[q] = 128 - q;
        bits.ApplyX(word[q]);
    }
    const auto packedWord = bits.SampleCounts(word, 5);
    if (!CliffordCheck(packedWord.size() == 1 && packedWord.at(std::numeric_limits<size_t>::max()) == 5,
                       "Packed counts include the highest size_t bit"))
        return false;
    // A nontrivial zero-preserving Clifford changes the logical encoding of
    // a GHZ state. This exercises Y factors and signs in the rank-one shortcut.
    for (size_t n : {7, 65, 129})
        for (int rep = 0; rep < 12; ++rep)
        {
            CliffordRegressionSimulator encoded(n);
            encoded.SetMultithreading(false);
            for (size_t k = 0; k < 16 * n; ++k)
            {
                const size_t a = actions() % n, b = (a + 1 + actions() % (n - 1)) % n;
                switch (actions() % 4)
                {
                case 0:
                    encoded.ApplyS(a);
                    break;
                case 1:
                    encoded.ApplyCZ(a, b);
                    break;
                case 2:
                    encoded.ApplyCX(a, b);
                    break;
                default:
                    encoded.ApplyZ(a);
                    break;
                }
            }
            encoded.ApplyH(0);
            for (size_t q = 1; q < n; ++q)
                encoded.ApplyCX(q, 0);
            for (size_t q = 0; q < n; ++q)
                if (actions() & 1)
                    encoded.ApplyX(q);
            CliffordRegressionSimulator reference(encoded);
            encoded.SetSeed(730 + rep);
            reference.SetSeed(730 + rep);
            std::vector<size_t> selected(n + 2);
            std::vector<bool> expected(n + 2);
            for (size_t q = 0; q < selected.size(); ++q)
            {
                selected[q] = n - 1 - q % n;
                expected[q] = reference.MeasureQubit(selected[q]);
            }
            const auto actual = encoded.SampleCountsMany(selected, 1);
            if (!CliffordCheck(actual.size() == 1 && actual.begin()->first == expected &&
                                   encoded.RandomEngine() == reference.RandomEngine(),
                               "Rank-one sampling phase and RNG"))
                return false;
        }
    return true;
}

void PrepareCliffordRegressionState(CliffordRegressionSimulator &sim, bool cnotOnly = false)
{
    std::mt19937 rng(73);
    const size_t n = sim.getNrQubits();
    sim.SetMultithreading(false);
    for (size_t k = 0; k < 32 * n; ++k)
    {
        const size_t a = rng() % n;
        size_t b = rng() % n;
        if (a == b)
            b = (a + 1) % n;
        const unsigned gate = cnotOnly ? 2 : rng() % 3;
        if (gate == 0)
            sim.ApplyH(a);
        else if (gate == 1)
            sim.ApplyS(a);
        else if (n > 1)
            sim.ApplyCX(a, b);
    }
}

bool CliffordResetRegression()
{
    for (size_t n : {0, 1, 7, 65, 1025})
    {
        CliffordRegressionSimulator sim(n), zero(n);
        PrepareCliffordRegressionState(sim);
        if (n != 0)
        {
            sim.ApplyH(0);
            sim.ApplyS(0);
            sim.ApplyX(0);
        }
        const CliffordRegressionSimulator before(sim);
        sim.SaveState();
        for (int repeat = 0; repeat < 2; ++repeat)
        {
            sim.Reset();
            if (!CliffordCheck(sim.SameTableau(zero), "Reset must clear every X/Z bit and sign"))
                return false;
            for (size_t q = 0; q < n; ++q)
                if (!CliffordCheck(sim.GetQubitProbability(q) == 0.0, "Reset must produce |0>"))
                    return false;
        }
        sim.RestoreState();
        if (!CliffordCheck(sim.SameTableau(before), "Reset must preserve the saved state"))
            return false;
    }
    return true;
}

bool CliffordSwapRegression()
{
    // All two-qubit Paulis, including Y and both signs, must just permute.
    const char labels[] = {'I', 'Z', 'X', 'Y'};
    for (unsigned pauli = 0; pauli < 16; ++pauli)
        for (bool sign : {false, true})
        {
            QC::PauliStringXZWithSign p(2);
            for (size_t q = 0; q < 2; ++q)
            {
                p.X[q] = (pauli >> (2 * q + 1)) & 1;
                p.Z[q] = (pauli >> (2 * q)) & 1;
            }
            p.PhaseSign = sign;
            QC::PauliStringXZWithCoefficient weighted(2);
            weighted.X = p.X;
            weighted.Z = p.Z;
            weighted.Coefficient = -0.375;
            p.ApplySwap(0, 1);
            weighted.ApplySwap(0, 1);
            std::string expected;
            expected += labels[(pauli >> 2) & 3];
            expected += labels[pauli & 3];
            if (!CliffordCheck(p.PauliStringXZ::ToString() == expected && p.PhaseSign == sign &&
                                   weighted.PauliStringXZ::ToString() == expected && weighted.Coefficient == -0.375,
                               "SWAP must preserve signs and coefficients"))
                return false;
        }

    for (size_t n : {2, 65, 1025})
    {
        CliffordRegressionSimulator initial(n);
        PrepareCliffordRegressionState(initial);
        for (bool parallel : {false, true})
        {
            CliffordRegressionSimulator actual(initial), expected(initial);
            actual.SetMultithreading(parallel);
            for (size_t a : {size_t(0), n / 2, n - 1})
            {
                const size_t b = (a + 1) % n;
                actual.ApplySwap(a, b);
                expected.ApplyCX(a, b);
                expected.ApplyCX(b, a);
                expected.ApplyCX(a, b);
                if (!CliffordCheck(actual.SameTableau(expected), "SWAP must equal three CNOTs"))
                    return false;
                actual.ApplySwap(a, a);
                if (!CliffordCheck(actual.SameTableau(expected), "SWAP on one index must be the identity"))
                    return false;
            }
        }
    }
    return true;
}

bool CliffordMeasurementRegression()
{
    CliffordRegressionSimulator single(1);
    single.ApplyS(0);
    single.ApplyH(0);
    const bool outcome = single.MeasureQubit(0); // Previously asserted on the paired destabilizer.
    if (!CliffordCheck(single.MeasureQubit(0) == outcome && single.GetQubitProbability(0) == double(outcome),
                       "Measurement after S,H must collapse consistently"))
        return false;

    for (size_t n : {2, 65, 1025})
    {
        // Two GHZ groups have four equally likely basis outcomes. S,H ensures
        // the probability path also encounters a pivot with an X destabilizer.
        CliffordRegressionSimulator groups(n);
        groups.SetMultithreading(false);
        groups.ApplyS(0);
        groups.ApplyH(0);
        groups.ApplyS(1);
        groups.ApplyH(1);
        for (size_t q = 2; q < n; ++q)
            groups.ApplyCX(q, q % 2);
        const CliffordRegressionSimulator before(groups);
        for (bool parallel : {false, true})
        {
            groups.SetMultithreading(parallel);
            for (unsigned pattern = 0; pattern < 4; ++pattern)
            {
                std::vector<bool> bits(n);
                for (size_t q = 0; q < n; ++q)
                    bits[q] = (pattern >> (q % 2)) & 1;
                if (!CliffordCheck(groups.getBasisStateProbability(bits) == 0.25 && groups.SameTableau(before),
                                   "Basis probabilities must preserve the tableau and support"))
                    return false;
                if (n > 2)
                {
                    bits[n - 1] = !bits[n - 1];
                    if (!CliffordCheck(groups.getBasisStateProbability(bits) == 0.0 && groups.SameTableau(before),
                                       "Impossible basis states must have zero probability"))
                        return false;
                }
            }
        }
    }

    std::mt19937 rng(75);
    for (size_t n : {1024, 1025, 1031, 2049})
        for (int repeat = 0; repeat < 100; ++repeat)
        {
            QC::PauliStringXZWithSign p(n), identity(n);
            for (size_t q = 0; q < n; ++q)
            {
                p.X[q] = rng() & 1;
                p.Z[q] = rng() & 1;
            }
            p.PhaseSign = (rng() & 1) != 0;
            const auto source = p;
            p.Multiply(source, true);
            if (!CliffordCheck(p == identity && !p.PhaseSign, "P * P must be identity at word boundaries"))
                return false;
            p = source;
            p.Multiply(p, true);
            if (!CliffordCheck(p == identity && !p.PhaseSign, "Aliased P * P must be identity"))
                return false;
        }

    for (size_t n : {63, 64, 65, 511, 512, 513, 1023, 1024, 1025, 1031, 2048, 2049, 4097})
    {
        CliffordRegressionSimulator initial(n);
        PrepareCliffordRegressionState(initial);
        for (unsigned seed : {74, 119})
        {
            CliffordRegressionSimulator serial(initial), parallel(initial);
            serial.SetSeed(seed);
            parallel.SetSeed(seed);
            parallel.SetMultithreading(true);
            for (size_t q : {size_t(0), n / 2, n - 1, size_t(1)})
            {
                const bool result = serial.MeasureQubit(q);
                if (!CliffordCheck(parallel.MeasureQubit(q) == result && serial.SameTableau(parallel),
                                   "Parallel measurement must match the complete serial tableau"))
                    return false;
                if (!CliffordCheck(parallel.MeasureQubit(q) == result && serial.SameTableau(parallel),
                                   "Repeated measurement must be deterministic"))
                    return false;
                serial.ApplyH((q + 1) % n);
                parallel.ApplyH((q + 1) % n);
            }
        }
    }

    CliffordRegressionSimulator deterministic(1025);
    PrepareCliffordRegressionState(deterministic, true);
    deterministic.SetMultithreading(true);
    for (size_t q : {0, 512, 1024})
        if (!CliffordCheck(deterministic.GetQubitProbability(q) == 0.0 && !deterministic.MeasureQubit(q),
                           "Dense deterministic measurements must preserve |0>"))
            return false;
    return true;
}
} // namespace

namespace
{
bool CliffordMixedRegression()
{
    size_t circuits = 0, measurements = 0, checks = 0, expectations = 0;
    std::mt19937 actions(998877);
    for (size_t n = 2; n <= 8; ++n)
        for (size_t rep = 0; rep < 100; ++rep)
        {
            CliffordRegressionSimulator s(n);
            s.SetSeed(987654321 + n * 100 + rep);
            s.SetMultithreading(false);
            QC::QubitRegister<> ref(n);
            ref.SetMultithreading(false);
            s.SaveState();
            ref.SaveState();
            bool saved = true;
            for (size_t step = 0; step < 160; ++step)
            {
                const auto op = actions() % 100;
                const size_t q = actions() % n, other = (q + 1 + actions() % (n - 1)) % n;
                if (op < 72)
                {
                    const int code = actions() % 15;
                    ApplyGate(s, code, int(q), int(other));
                    ref.ApplyGate(*GetGate(code), q, other);
                }
                else if (op < 87)
                {
                    const bool value = s.MeasureQubit(q);
                    auto v = ref.getRegisterStorage();
                    for (size_t b = 0; b < size_t(v.size()); ++b)
                        if (bool((b >> q) & 1) != value)
                            v[b] = 0.;
                    if (v.squaredNorm() < 1e-10)
                    {
                        std::cout << "impossible measurement\n";
                        return false;
                    }
                    ref.setRegisterStorage(v);
                    ++measurements;
                    if (s.MeasureQubit(q) != value)
                        return false;
                }
                else if (op < 90)
                {
                    s.Reset();
                    ref.Reset();
                }
                else if (op < 93)
                {
                    s.SaveState();
                    ref.SaveState();
                    saved = true;
                }
                else if (op < 96 && saved)
                {
                    s.RestoreState();
                    ref.RestoreState();
                }
                else if (saved)
                {
                    s.RestoreSavedStateDestructive();
                    ref.RestoreStateDestructive();
                    saved = false;
                }

                for (int sample = 0; sample < 2; ++sample)
                {
                    std::vector<QC::Gates::AppliedGate<>> p;
                    std::string pauli;
                    ConstructPauliString(n, pauli, p);
                    if (!approxEqual(s.ExpectationValue(pauli), ref.ExpectationValue(p), 1e-8))
                    {
                        std::cout << "expectation mismatch n=" << n << " rep=" << rep << " step=" << step
                                  << " op=" << op << '\n';
                        return false;
                    }
                    ++expectations;
                }
                if (step % 5 == 0)
                {
                    const CliffordRegressionSimulator original(s);
                    if (!CheckProbability(ref, s) || !CheckAllStatesProbability(ref, s))
                    {
                        std::cout << "probability mismatch n=" << n << " rep=" << rep << " step=" << step << '\n';
                        return false;
                    }
                    if (!s.SameTableau(original))
                    {
                        std::cout << "query changed tableau\n";
                        return false;
                    }
                    auto cloned = s.Clone();
                    for (size_t b = 0; b < (size_t(1) << n); ++b)
                        if (!approxEqual(cloned->getBasisStateProbability(b), ref.getBasisStateProbability(b), 1e-8))
                            return false;
                    ++checks;
                }
            }
            ++circuits;
        }
    std::cout << "PASS " << circuits << " mixed circuits, 160 operations each, " << measurements
              << " forced-reference measurements, " << expectations << " Pauli expectations, " << checks
              << " full probability/marginal/clone/tableau checks, with reset/save/restore\n";
    return true;
}
} // namespace

bool CliffordRegressionTests()
{
    std::cout << "\nClifford gates, inverse tableau, probability, sampling and boundary regressions" << std::endl;
    if (!CliffordGenericGateRegression() || !CliffordSingleRotationRegression() || !CliffordMixedRegression() ||
        !CliffordPackedRegression() || !CliffordDistributionRegression() || !CliffordBoundaryRegression() ||
        !CliffordCountsRegression() || !CliffordResetRegression() || !CliffordSwapRegression() ||
        !CliffordMeasurementRegression() || !CliffordCompactGateRegression())
        return false;
    std::cout << "Success" << std::endl;
    return true;
}

bool CliffordSimulatorTests()
{
    if (!CliffordRegressionTests())
        return false;

    const size_t nrTests = 10;
    const size_t nrShots = 100000;
    const double errorThreshold = 0.01;

    std::uniform_int_distribution gateDistr(0, 14);
    std::uniform_int_distribution nrGatesDistr(5, 20);

    std::cout << "\nClifford gates simulator" << std::endl;

    for (size_t nrQubits = 4; nrQubits < 12; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);

        for (size_t t = 0; t < nrTests; ++t)
        {
            std::unordered_map<size_t, int> results1;
            std::unordered_map<size_t, int> results2;

            // generate random gates, creating a circuit, then apply the random circuits on both simulators
            const size_t nrGates = nrGatesDistr(gen);
            std::vector<int> gates(nrGates);
            std::vector<size_t> qubits1(nrGates);
            std::vector<size_t> qubits2(nrGates);

            ConstructCircuit(nrQubits, gates, qubits1, qubits2, gateDistr, qubitDistr);

            QC::QubitRegister qubitRegister(nrQubits);
            QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

            for (int j = 0; j < static_cast<int>(gates.size()); ++j)
            {
                ApplyGate(cliffordSim, gates[j], qubits1[j], qubits2[j]);
                const auto gateptr = GetGate(gates[j]);
                qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
            }

            if (!CheckProbability(qubitRegister, cliffordSim))
                return false;

            // another way of testing is now available
            if (!CheckAllStatesProbability(qubitRegister, cliffordSim))
                return false;

            // applying the measurement again on the same qubit should give the same result
            if (!CheckMeasurements(cliffordSim))
                return false;

            ExecuteCircuit(nrShots, nrQubits, gates, qubits1, qubits2, results1, results2);

            // check to see if the results are close enough
            for (const auto val : results1)
            {
                if (results2.find(val.first) == results2.end())
                    continue;

                if (std::abs(static_cast<double>(val.second) - results2[val.first]) / nrShots > errorThreshold)
                {
                    std::cout << "\nFailed" << std::endl;
                    std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
                    std::cout << "Result 1: " << static_cast<double>(val.second) / nrShots
                              << ", Result 2: " << static_cast<double>(results2[val.first]) / nrShots << std::endl;
                    return false;
                }
            }

            for (const auto val : results2)
            {
                if (results1.find(val.first) == results1.end())
                    continue;

                if (std::abs(static_cast<double>(val.second) - results1[val.first]) / nrShots > errorThreshold)
                {
                    std::cout << "\nFailed" << std::endl;
                    std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
                    std::cout << "Result 1: " << static_cast<double>(results1[val.first]) / nrShots
                              << ", Result 2: " << static_cast<double>(val.second) / nrShots << std::endl;
                    return false;
                }
            }

            std::cout << '.';
        }
    }

    std::cout << std::endl;

    std::cout << "Success" << std::endl;

    return true;
}

void ConstructPauliString(size_t nrQubits, std::string &pauliStr, std::vector<QC::Gates::AppliedGate<>> &expGates)
{
    static const QC::Gates::PauliXGate xgate;
    static const QC::Gates::PauliYGate ygate;
    static const QC::Gates::PauliZGate zgate;

    std::uniform_int_distribution pauliDistr(0, 3);

    for (int j = 0; j < static_cast<int>(nrQubits); ++j)
    {
        const int p = pauliDistr(gen);
        switch (p)
        {
        case 0:
            pauliStr += 'I';
            break;
        case 1:
            pauliStr += 'X';
            expGates.emplace_back(xgate.getRawOperatorMatrix(), j);
            break;
        case 2:
            pauliStr += 'Y';
            expGates.emplace_back(ygate.getRawOperatorMatrix(), j);
            break;
        case 3:
            pauliStr += 'Z';
            expGates.emplace_back(zgate.getRawOperatorMatrix(), j);
            break;
        }
    }
}

bool CliffordExpectationValuesTests()
{
    const size_t nrTests = 100;

    std::uniform_int_distribution gateDistr(0, 14);
    std::uniform_int_distribution nrGatesDistr(5, 20);

    std::cout << "\nClifford expectation values" << std::endl;

    for (size_t nrQubits = 2; nrQubits < 20; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);

        for (size_t t = 0; t < nrTests; ++t)
        {
            std::unordered_map<size_t, int> results1;
            std::unordered_map<size_t, int> results2;

            // generate random gates, creating a circuit, then apply the random circuits on both simulators
            const size_t nrGates = nrGatesDistr(gen);
            std::vector<int> gates(nrGates);
            std::vector<size_t> qubits1(nrGates);
            std::vector<size_t> qubits2(nrGates);

            ConstructCircuit(nrQubits, gates, qubits1, qubits2, gateDistr, qubitDistr);

            QC::QubitRegister qubitRegister(nrQubits);
            QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

            for (int j = 0; j < static_cast<int>(gates.size()); ++j)
            {
                ApplyGate(cliffordSim, gates[j], qubits1[j], qubits2[j]);
                const auto gateptr = GetGate(gates[j]);
                qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
            }

            std::vector<QC::Gates::AppliedGate<>> expGates;
            expGates.reserve(nrQubits);
            std::string pauliStr;

            ConstructPauliString(nrQubits, pauliStr, expGates);

            const auto exp1 = qubitRegister.ExpectationValue(expGates);
            const auto exp2 = cliffordSim.ExpectationValue(pauliStr);
            if (!approxEqual(exp1, exp2, 1E-7))
            {
                std::cout << std::endl
                          << "Expectation values are not equal for stabilizer and statevector simulator for "
                          << nrQubits << " qubits, values: " << exp2 << ", " << exp1 << std::endl;

                std::cout << "Pauli string: " << pauliStr << std::endl;

                std::cout << "Circuit:" << std::endl;
                for (int j = 0; j < static_cast<int>(gates.size()); ++j)
                    PrintGate(gates[j], qubits1[j], qubits2[j]);

                return false;
            }
        }
        std::cout << '.';
    }

    std::cout << std::endl;

    std::cout << "Success" << std::endl;

    return true;
}
