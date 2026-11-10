#include "Tests.h"

#include "DensityMatrix.h"
#include "NonCliffordTestGates.h"
#include "PauliPropagator.h"
#include "PauliPropagatorRegressionTests.h"
#include "QubitRegister.h"

void ApplyTwoQubitsGate(QC::PauliPropagator &simulator, int code, int qubit1, int qubit2)
{
    switch (code)
    {
    case 9:
        simulator.ApplyCX(qubit2, qubit1);
        break;
    case 10:
        simulator.ApplyCY(qubit2, qubit1);
        break;
    case 11:
        simulator.ApplyCZ(qubit2, qubit1);
        break;
    case 12:
        simulator.ApplySWAP(qubit1, qubit2);
        break;
    case 13:
        simulator.ApplyISWAP(qubit1, qubit2);
        break;
    case 14:
        simulator.ApplyISWAPDG(qubit1, qubit2);
        break;
    }
}

// The codes of NonCliffordTestGates.h, where qubit1 is the target and qubit2
// the control; the simulator takes controls first.
void ApplyGate(QC::PauliPropagator &simulator, int code, int qubit1, int qubit2, double angle = 0.0, int qubit3 = 0)
{
    const auto angles = DeriveGateAngles(angle);
    switch (code)
    {
    case 0:
        simulator.ApplyH(qubit1);
        break;
    case 1:
        simulator.ApplyS(qubit1);
        break;
    case 2:
        simulator.ApplySDG(qubit1);
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
        simulator.ApplySX(qubit1);
        break;
    case 7:
        simulator.ApplySXDG(qubit1);
        break;
    case 8:
        simulator.ApplyK(qubit1);
        break;
    case 15:
        simulator.ApplyRX(qubit1, angle);
        break;
    case 16:
        simulator.ApplyRY(qubit1, angle);
        break;
    case 17:
        simulator.ApplyRZ(qubit1, angle);
        break;
    case CodeU:
        simulator.ApplyU(qubit1, angles.theta, angles.phi, angles.lambda);
        break;
    case CodeCU:
        simulator.ApplyCU(qubit2, qubit1, angles.theta, angles.phi, angles.lambda, angles.gamma);
        break;
    case CodeCRx:
        simulator.ApplyCRX(qubit2, qubit1, angle);
        break;
    case CodeCRy:
        simulator.ApplyCRY(qubit2, qubit1, angle);
        break;
    case CodeCRz:
        simulator.ApplyCRZ(qubit2, qubit1, angle);
        break;
    case CodeCP:
        simulator.ApplyCP(qubit2, qubit1, angle);
        break;
    case CodeCS:
        simulator.ApplyCS(qubit2, qubit1);
        break;
    case CodeCSdg:
        simulator.ApplyCSDAG(qubit2, qubit1);
        break;
    case CodeCSx:
        simulator.ApplyCSX(qubit2, qubit1);
        break;
    case CodeCSxDag:
        simulator.ApplyCSXDAG(qubit2, qubit1);
        break;
    case CodeCH:
        simulator.ApplyCH(qubit2, qubit1);
        break;
    case CodeCCX:
        simulator.ApplyCCX(qubit2, qubit3, qubit1);
        break;
    case CodeCSwap:
        simulator.ApplyCSwap(qubit3, qubit1, qubit2);
        break;
    default:
        ApplyTwoQubitsGate(simulator, code, qubit1, qubit2);
        break;
    }
}

bool CheckResults(int shots, int nrQubits, std::unordered_map<size_t, size_t> &sampledResultsPauli,
                  const std::unordered_map<size_t, size_t> &sampledResultsStatevector, bool measured = false)
{
    for (const auto &kv : sampledResultsStatevector)
    {
        const size_t state = kv.first;
        const size_t countStatevector = kv.second;
        const size_t countPauli = sampledResultsPauli[state];
        const double freqStatevector = static_cast<double>(countStatevector) / static_cast<double>(shots);
        const double freqPauli = static_cast<double>(countPauli) / static_cast<double>(shots);

        if (std::abs(freqStatevector - freqPauli) > 0.1)
        {
            std::cout << std::endl
                      << "Sampled frequencies - " << (measured ? "by measurement" : "by sampling")
                      << " - are not equal for pauli propagator and statevector simulator for " << nrQubits
                      << " qubits, measurement " << state << ", values: " << freqPauli << ", " << freqStatevector
                      << std::endl;
            return false;
        }
    }

    return true;
}

bool CheckProbability(size_t nrQubits, QC::PauliPropagator &pauliSimulator, QC::QubitRegister<> &qubitRegister)
{
    for (size_t q = 0; q < nrQubits; ++q)
    {
        const auto p0 = qubitRegister.GetQubitProbability(q);
        const auto p1 = 1. - pauliSimulator.Probability0(static_cast<int>(q));
        if (std::abs(p0 - p1) > 1E-7)
        {
            std::cout << std::endl
                      << "Probability values are not equal for pauli propagator and statevector simulator for "
                      << nrQubits << " qubits, qubit " << q << ", values: " << p1 << ", " << p0 << std::endl;
            return false;
        }
    }

    return true;
}

// Without clifford, a fifth of the gates become rotations or other supported
// non-Clifford gates; a third of their angles are multiples of pi/2.
void ExecuteCircuit(QC::QubitRegister<> &qubitRegister, QC::PauliPropagator &pauliSimulator, std::vector<int> &gates,
                    std::vector<size_t> &qubits1, std::vector<size_t> &qubits2,
                    std::uniform_real_distribution<double> &angleDistr, std::bernoulli_distribution &boolDistr,
                    bool clifford)
{
    const size_t nrQubits = qubitRegister.getNrQubits();
    std::uniform_int_distribution<int> quarterTurns(-4, 4);
    for (int j = 0; j < static_cast<int>(gates.size()); ++j)
    {
        double angle = 0.0;
        size_t qubit3 = 0;
        if (!clifford && boolDistr(gen))
        {
            gates[j] = RandomNonCliffordCode(gen, nrQubits);
            angle = quarterTurns(gen) % 3 == 0 ? 0.5 * M_PI * quarterTurns(gen) : angleDistr(gen);
            if (gates[j] >= CodeCCX)
                qubit3 = ThirdQubit(nrQubits, qubits1[j], qubits2[j]);
        }

        ApplyStatevectorGate(qubitRegister, gates[j], qubits1[j], qubits2[j], angle, qubit3);

        ApplyGate(pauliSimulator, gates[j], (int)qubits1[j], (int)qubits2[j], angle, (int)qubit3);
    }
}

std::unordered_map<size_t, size_t> CollectSampledResults(QC::PauliPropagator &pauliSimulator,
                                                         std::vector<int> &measQubits, int shots, int nrQubits,
                                                         std::mt19937 &g)
{
    std::unordered_map<size_t, size_t> sampledResultsPauli;
    for (int j = 0; j < shots; ++j)
    {
        std::shuffle(measQubits.begin(), measQubits.end(), g);
        auto res = pauliSimulator.Sample(measQubits);

        size_t result = 0;
        for (size_t q = 0; q < static_cast<size_t>(nrQubits); ++q)
        {
            if (res[q])
            {
                size_t qubit = measQubits[q];
                result |= (1ULL << qubit);
            }
        }

        ++sampledResultsPauli[result];
    }
    return sampledResultsPauli;
}

std::unordered_map<size_t, size_t> CollectMeasuredResults(QC::PauliPropagator &pauliSimulator,
                                                          std::vector<int> &measQubits, int shots, int nrQubits,
                                                          std::mt19937 &g)
{
    std::unordered_map<size_t, size_t> sampledResultsPauli;
    pauliSimulator.SaveState();
    for (int j = 0; j < shots; ++j)
    {
        std::shuffle(measQubits.begin(), measQubits.end(), g);
        auto res = pauliSimulator.Measure(measQubits);

        size_t result = 0;
        for (size_t q = 0; q < static_cast<size_t>(nrQubits); ++q)
        {
            if (res[q])
            {
                size_t qubit = measQubits[q];
                result |= (1ULL << qubit);
            }
        }

        ++sampledResultsPauli[result];

        pauliSimulator.RestoreState();
    }
    return sampledResultsPauli;
}

bool TestPauliPropagatorCorNC(bool clifford = true)
{
    std::cout << "\nPauli Propagator tests with " << (clifford ? "Clifford" : "non-Clifford") << " gates" << std::endl;
    const size_t nrTests = 100;
    const size_t maxQubits = 20;

    std::uniform_int_distribution gateDistr(0, 14);
    std::uniform_int_distribution nrGatesDistr(5, 20);
    std::uniform_real_distribution angleDistr(-2. * M_PI, 2. * M_PI);
    std::bernoulli_distribution boolDistr(0.2);

    QC::PauliPropagator pauliSimulator;
    pauliSimulator.EnableParallel();

    for (size_t nrQubits = 2; nrQubits < maxQubits; ++nrQubits)
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

            pauliSimulator.SetNrQubits(static_cast<int>(nrQubits));

            QC::QubitRegister qubitRegister(nrQubits);

            ExecuteCircuit(qubitRegister, pauliSimulator, gates, qubits1, qubits2, angleDistr, boolDistr, clifford);

            std::vector<QC::Gates::AppliedGate<>> expGates;
            expGates.reserve(nrQubits);
            std::string pauliStr;

            ConstructPauliString(nrQubits, pauliStr, expGates);

            const auto exp1 = qubitRegister.ExpectationValue(expGates);
            const auto exp2 = pauliSimulator.ExpectationValue(pauliStr);
            if (!approxEqual(exp1, exp2, 1E-7))
            {
                std::cout << std::endl
                          << "Expectation values are not equal for pauli propagator and statevector simulator for "
                          << nrQubits << " qubits, values: " << exp2 << ", " << exp1 << std::endl;

                std::cout << "Pauli string: " << pauliStr << std::endl;

                return false;
            }

            if (!CheckProbability(nrQubits, pauliSimulator, qubitRegister))
                return false;

            // for big number of qubits sampling on the Pauli propagator becomes too slow
            if (nrQubits < 9 && t < 10) // also limit the number of tests, this can become very slow
            {
                std::vector<int> measQubits(nrQubits);
                std::iota(measQubits.begin(), measQubits.end(), 0);

                int shots = 1000;
                auto sampledResultsStatevector = qubitRegister.RepeatedMeasureUnordered(shots);

                auto sampledResultsPauli = CollectSampledResults(pauliSimulator, measQubits, shots, (int)nrQubits, gen);
                if (!CheckResults(shots, (int)nrQubits, sampledResultsPauli, sampledResultsStatevector, true))
                    return false;

                // now do the same but with measurements!
                sampledResultsPauli = CollectMeasuredResults(pauliSimulator, measQubits, shots, (int)nrQubits, gen);
                if (!CheckResults(shots, (int)nrQubits, sampledResultsPauli, sampledResultsStatevector, true))
                    return false;

                // this is costly so check only a few values
                for (size_t state = 0; state < 10; ++state)
                {
                    const auto pStatevector = qubitRegister.getBasisStateProbability(state);
                    const auto pPauli = pauliSimulator.Probability(state);
                    if (std::abs(pStatevector - pPauli) > 1E-5)
                    {
                        std::cout
                            << std::endl
                            << "Probability values are not equal for pauli propagator and statevector simulator for "
                            << nrQubits << " qubits, state " << state << ", values: " << pPauli << ", " << pStatevector
                            << std::endl;
                        return false;
                    }
                }
            }

            pauliSimulator.ClearOperations();
        }
        std::cout << '.';
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

// Noise channels. The propagator applies the adjoint channel to the observable;
// the density matrix applies the channel to the state. Both simulators must
// agree on every expectation value and basis probability, including for a
// two-qubit Kraus channel that is not symmetric in its qubits (which fixes the
// qubit order of ApplyKrausChannel against DensityMatrix::ApplyChannel).
bool TestPauliPropagatorNoise()
{
    std::cout << "\nPauli Propagator noise channels against the density matrix" << std::endl;
    using Matrix = Eigen::MatrixXcd;
    std::uniform_int_distribution gateDistr(0, 17);
    std::uniform_int_distribution channelDistr(0, 8);
    std::uniform_int_distribution pauliDistr(0, 3);
    std::uniform_real_distribution angleDistr(-2. * M_PI, 2. * M_PI);
    std::uniform_real_distribution probDistr(0., 1.);
    std::bernoulli_distribution noiseDistr(0.4);

    const auto damping = [](double gamma) {
        Matrix k0 = Matrix::Zero(2, 2), k1 = Matrix::Zero(2, 2);
        k0(0, 0) = 1.;
        k0(1, 1) = std::sqrt(1. - gamma);
        k1(0, 1) = std::sqrt(gamma);
        return std::vector<Matrix>{k0, k1};
    };
    const auto bitFlip = [](double p) {
        Matrix x = Matrix::Zero(2, 2);
        x(0, 1) = x(1, 0) = 1.;
        return std::vector<Matrix>{std::sqrt(1. - p) * Matrix::Identity(2, 2), std::sqrt(p) * x};
    };
    // second (x) first, with the first qubit in the low bit of the index.
    const auto tensor = [](const Matrix &second, const Matrix &first) {
        Matrix result = Matrix::Zero(4, 4);
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                result(r, c) = second(r >> 1, c >> 1) * first(r & 1, c & 1);
        return result;
    };

    for (size_t nrQubits = 2; nrQubits <= 6; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);
        for (size_t t = 0; t < 40; ++t)
        {
            QC::PauliPropagator pauliSimulator;
            pauliSimulator.SetNrQubits(static_cast<int>(nrQubits));
            if (t % 2)
                pauliSimulator.SetStepsBetweenDeduplication(3);
            QC::DensityMatrix<> densityMatrix(nrQubits);
            for (int step = 0; step < 25; ++step)
            {
                const int q1 = qubitDistr(gen);
                int q2 = qubitDistr(gen);
                while (q2 == q1)
                    q2 = qubitDistr(gen);
                if (!noiseDistr(gen))
                {
                    const int code = gateDistr(gen);
                    const double angle = angleDistr(gen);
                    const auto gate = GetGate(code, angle);
                    densityMatrix.ApplyGate(*gate, q1, q2);
                    ApplyGate(pauliSimulator, code, q1, q2, angle);
                    continue;
                }
                const double p = probDistr(gen);
                switch (channelDistr(gen))
                {
                case 0:
                    densityMatrix.ApplyBitFlipNoise(q1, p);
                    pauliSimulator.ApplyBitFlipNoise(q1, p);
                    break;
                case 1:
                    densityMatrix.ApplyPhaseFlipNoise(q1, p);
                    pauliSimulator.ApplyPhaseFlipNoise(q1, p);
                    break;
                case 2:
                    densityMatrix.ApplyDepolarizingNoise(q1, p);
                    pauliSimulator.ApplyDepolarizingNoise(q1, p);
                    break;
                case 3:
                    densityMatrix.ApplyAmplitudeDamping(q1, p);
                    pauliSimulator.ApplyAmplitudeDamping(q1, p);
                    break;
                case 4:
                    densityMatrix.ApplyPhaseDamping(q1, p);
                    pauliSimulator.ApplyPhaseDamping(q1, p);
                    break;
                case 5:
                    densityMatrix.ApplyReset(q1);
                    pauliSimulator.ApplyReset(q1);
                    break;
                case 6: {
                    const auto kraus = damping(p);
                    densityMatrix.ApplyChannel(kraus, q1);
                    pauliSimulator.ApplyKrausChannel(kraus, q1);
                    break;
                }
                case 7: {
                    // Damping on the first qubit and a bit flip on the second:
                    // exchanging the qubits gives a different channel.
                    std::vector<Matrix> kraus;
                    for (const auto &second : bitFlip(0.5 * p))
                        for (const auto &first : damping(p))
                            kraus.push_back(tensor(second, first));
                    densityMatrix.ApplyChannel(kraus, q1, q2);
                    pauliSimulator.ApplyKrausChannel(kraus, q1, q2);
                    break;
                }
                default: {
                    // Thermal relaxation with T1 < T2 <= 2 T1: damping followed
                    // by the dephasing that brings the coherence down to exp(-t/T2).
                    const double t1 = 100., t2 = 100. * (1. + p), duration = 60. * probDistr(gen);
                    const double gamma = -std::expm1(-duration / t1);
                    const double extra = std::exp(-duration / t2) / std::sqrt(1. - gamma);
                    densityMatrix.ApplyAmplitudeDamping(q1, gamma);
                    densityMatrix.ApplyPhaseFlipNoise(q1, 0.5 * (1. - extra));
                    pauliSimulator.ApplyThermalRelaxation(q1, duration, t1, t2);
                    break;
                }
                }
            }

            for (int k = 0; k < 40; ++k)
            {
                std::string pauliStr(nrQubits, 'I');
                for (auto &c : pauliStr)
                    c = "IXYZ"[pauliDistr(gen)];
                const double expected = densityMatrix.ExpectationValue(pauliStr).real();
                const double value = pauliSimulator.ExpectationValue(pauliStr);
                if (!approxEqual(expected, value, 1E-9))
                {
                    std::cout << std::endl
                              << "Noisy expectation values are not equal for pauli propagator and density matrix for "
                              << nrQubits << " qubits, values: " << value << ", " << expected << std::endl;
                    std::cout << "Pauli string: " << pauliStr << std::endl;
                    return false;
                }
            }
            for (size_t state = 0; state < (size_t(1) << nrQubits); ++state)
            {
                const double expected = densityMatrix.getBasisStateProbability(state);
                const double value = pauliSimulator.Probability(state);
                if (std::abs(expected - value) > 1E-8)
                {
                    std::cout << std::endl
                              << "Noisy probabilities are not equal for pauli propagator and density matrix for "
                              << nrQubits << " qubits, state " << state << ", values: " << value << ", " << expected
                              << std::endl;
                    return false;
                }
            }
        }
        std::cout << '.';
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool TestPauliPropagator()
{
    try
    {
        QC::PauliRegression::Run();
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        return false;
    }

    return TestPauliPropagatorNoise() && TestPauliPropagatorCorNC(true) && TestPauliPropagatorCorNC(false);
}
