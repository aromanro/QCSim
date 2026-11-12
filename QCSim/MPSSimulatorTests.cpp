#include <iostream>
#include <iterator>
#include <limits>
#include <memory>

#include "Tests.h"

#include "MPSSimulator.h"
#include "ProductSwapRegressionTests.h"
#include "QubitRegister.h"
#include "TensorNetworkRegressionTests.h"

#include <vector>

#define _USE_MATH_DEFINES
#include <future>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NR_QUBITS_LIMIT 9

void FillOneQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates)
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
    gates.emplace_back(std::make_shared<QC::Gates::SplitterGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::RxGate<>>(M_PI / 3));
    gates.emplace_back(std::make_shared<QC::Gates::RyGate<>>(M_PI / 7));
    gates.emplace_back(std::make_shared<QC::Gates::RzGate<>>(M_PI / 5));
    gates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 3, M_PI / 5));
}

void FillTwoQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates)
{
    gates.emplace_back(std::make_shared<QC::Gates::SwapGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::DecrementGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledYGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledZGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledHadamardGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledSquareRootNOTGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledSquareRootNOTDagGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledPhaseGate<>>());
    gates.emplace_back(std::make_shared<QC::Gates::ControlledPhaseShiftGate<>>(M_PI / 3));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledUGate<>>(M_PI / 3, M_PI / 7));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledRxGate<>>(M_PI / 5));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledRyGate<>>(M_PI / 3));
    gates.emplace_back(std::make_shared<QC::Gates::ControlledRzGate<>>(M_PI / 7));
}

bool OneAndTwoQubitGatesTest()
{
    std::cout << "\nMPS simulator state test for both one and two qubit gates" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(25, 50);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        for (int t = 0; t < 10; ++t)
        {
#ifdef _DEBUG
            std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl
                      << std::endl
                      << std::endl;
#endif

            QC::TensorNetworks::MPSSimulatorImpl mps(nrQubits);
            QC::QubitRegister reg(nrQubits);

            const int lim = nrGatesDistr(gen);

            for (int i = 0; i < lim; ++i)
            {
                const int gate = gateDistr(gen);
                const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = qubit1 + 1;

                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

#ifdef _DEBUG
                if (twoQubitsGate)
                    std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2
                              << std::endl;
#endif

                mps.ApplyGate(*gates[gate], qubit1, qubit2);
                reg.ApplyGate(*gates[gate], qubit1, qubit2);

                // now check the results, they should be the same
                const auto &regState = reg.getRegisterStorage();
                auto mpsState = mps.getRegisterStorage(); // this one is computed, returns value, not reference, not
                                                          // stored elsewhere

                // QC::QubitRegister regNorm(nrQubits);
                // regNorm.setRegisterStorage(mpsState);
                // mpsState = regNorm.getRegisterStorage();

                for (int s = 0; s < regState.size(); ++s)
                {
                    if (!approxEqual(regState[s], mpsState[s], 1E-3))
                    {
                        std::cout << "State " << s << " simulation test failed for the MPS simulator for " << nrQubits
                                  << " qubits" << std::endl;

                        std::cout << "Probability for the different states: " << std::norm(regState[s]) << " vs "
                                  << std::norm(mpsState[s]) << std::endl;

                        std::cout << "Reg state:\n" << regState << std::endl;
                        std::cout << "Reg state normalization: " << regState.norm() << std::endl;

                        std::cout << "MPS state:\n" << mpsState << std::endl;
                        std::cout << "MPS state normalization: " << mpsState.norm() << std::endl;

                        std::cout << std::endl;
                        for (int q = 0; q < nrQubits; ++q)
                            std::cout << "Qubit " << q << " reg probability: " << reg.GetQubitProbability(q)
                                      << " vs mps: " << mps.GetProbability(q, false) << std::endl;

                        std::cout << std::endl;
                        for (int state = 0; state < regState.size(); ++state)
                            std::cout << "State " << state << " reg probability: " << std::norm(regState[state])
                                      << " vs mps: " << std::norm(mpsState[state]) << std::endl;

                        return false;
                    }
                }
            }

#ifdef _DEBUG
            std::cout << "Test passed: " << t << std::endl;
#endif
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

std::vector<std::shared_ptr<QC::Gates::AppliedGate<>>> GenerateRandomCircuitWithGates(
    const std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates, int minGates, int maxGates, int nrQubits)
{
    std::vector<std::shared_ptr<QC::Gates::AppliedGate<>>> circuit;

    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);
    std::uniform_int_distribution nrGatesDistr(minGates, maxGates);

    std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
    std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

    const size_t lim = nrGatesDistr(gen);

    for (size_t i = 0; i < lim; ++i)
    {
        const int gate = gateDistr(gen);
        const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;

        auto appliedGate = std::make_shared<QC::Gates::AppliedGate<>>(gates[gate]->getRawOperatorMatrix());

        int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
        int qubit2 = qubit1 + 1;
        if (twoQubitsGate && dist_bool(gen))
            std::swap(qubit1, qubit2);

        appliedGate->setQubit1(qubit1);
        appliedGate->setQubit2(qubit2);

        circuit.push_back(std::move(appliedGate));
    }

    return circuit;
}

std::vector<std::shared_ptr<QC::Gates::AppliedGate<>>> GenerateRandomCircuitWithGatesNoAdjacent(
    const std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> &gates, int minGates, int maxGates, int nrQubits)
{
    std::vector<std::shared_ptr<QC::Gates::AppliedGate<>>> circuit;

    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);
    std::uniform_int_distribution nrGatesDistr(minGates, maxGates);

    std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
    std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

    const size_t lim = nrGatesDistr(gen);

    for (size_t i = 0; i < lim; ++i)
    {
        const int gate = gateDistr(gen);
        const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;

        auto appliedGate = std::make_shared<QC::Gates::AppliedGate<>>(gates[gate]->getRawOperatorMatrix());

        int qubit1 = qubitDistr(gen);
        int qubit2 = (qubit1 + 1 + qubitDistr2(gen)) % nrQubits;
        if (twoQubitsGate && dist_bool(gen))
            std::swap(qubit1, qubit2);

        appliedGate->setQubit1(qubit1);
        appliedGate->setQubit2(qubit2);

        circuit.push_back(std::move(appliedGate));
    }

    return circuit;
}

bool TestMeasurementsWithOneAndTwoQubitGatesCircuits()
{
    std::cout << "\nMPS simulator measurements test with circuits with both one and two qubit gates" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    const int nrMeasurements = 10000;

    for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        for (int t = 0; t < 5; ++t)
        {
            const auto circuit = GenerateRandomCircuitWithGates(gates, 25, 50, nrQubits);
            for (int c = 0; c < 3; ++c)
            {
                std::unordered_map<std::vector<bool>, int> measurementsRegMap;
                std::unordered_map<std::vector<bool>, int> measurementsMPSMap;

                for (int t2 = 0; t2 < nrMeasurements; ++t2)
                {
                    QC::TensorNetworks::MPSSimulatorImpl mps(nrQubits);
                    QC::QubitRegister reg(nrQubits);

                    for (const auto &gate : circuit)
                    {
                        mps.ApplyGate(*gate);
                        reg.ApplyGate(*gate);
                    }

                    std::vector<bool> measurementsReg(nrQubits);
                    std::vector<bool> measurementsMPS(nrQubits);

                    for (int q = 0; q < nrQubits; ++q)
                    {
                        measurementsReg[q] = reg.MeasureQubit(q);
                        measurementsMPS[q] = mps.MeasureQubit(q);
                    }
                    ++measurementsRegMap[measurementsReg];
                    ++measurementsMPSMap[measurementsMPS];
                }

                std::cout << "." << std::flush;

                for (const auto &[key, value] : measurementsRegMap)
                {
                    const double dif = abs((static_cast<double>(measurementsMPSMap[key] - value)) / nrMeasurements);
                    if (dif > 0.05)
                    {
                        std::cout << "Measurements test failed for the MPS simulator for " << nrQubits << " qubits"
                                  << std::endl;
                        std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
                        std::cout << "Difference: " << dif << std::endl;

                        std::cout << "Reg measurements:\n";
                        for (const auto &[key2, value2] : measurementsRegMap)
                        {
                            for (const auto &b : key2)
                                std::cout << b << " ";
                            std::cout << " : " << static_cast<double>(value2) / nrMeasurements << std::endl;
                        }

                        std::cout << "MPS measurements:\n";
                        for (const auto &[key2, value2] : measurementsMPSMap)
                        {
                            for (const auto &b : key2)
                                std::cout << b << " ";
                            std::cout << " : " << static_cast<double>(value2) / nrMeasurements << std::endl;
                        }

                        return false;
                    }
                }
            }
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool OneAndTwoQubitGatesTestMapped()
{
    std::cout << "\nMPS swapping/mapped simulator state test for both one and two qubit gates" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(25, 50);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        for (int t = 0; t < 10; ++t)
        {
#ifdef _DEBUG
            std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl
                      << std::endl
                      << std::endl;
#endif

            QC::TensorNetworks::MPSSimulator mps(nrQubits);
            QC::QubitRegister reg(nrQubits);

            const int lim = nrGatesDistr(gen);

            for (int i = 0; i < lim; ++i)
            {
                const int gate = gateDistr(gen);
                const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
                int qubit1 = qubitDistr(gen);
                int qubit2 = (qubit1 + 1 + qubitDistr2(gen)) % nrQubits;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

#ifdef _DEBUG
                if (twoQubitsGate)
                    std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2
                              << std::endl;
#endif

                mps.ApplyGate(*gates[gate], qubit1, qubit2);
                reg.ApplyGate(*gates[gate], qubit1, qubit2);

                // now check the results, they should be the same
                const auto &regState = reg.getRegisterStorage();
                auto mpsState = mps.getRegisterStorage(); // this one is computed, returns value, not reference, not
                                                          // stored elsewhere

                // QC::QubitRegister regNorm(nrQubits);
                // regNorm.setRegisterStorage(mpsState);
                // mpsState = regNorm.getRegisterStorage();

                for (int s = 0; s < regState.size(); ++s)
                {
                    if (!approxEqual(regState[s], mpsState[s], 1E-3))
                    {
                        std::cout << "State " << s << " simulation test failed for the MPS simulator for " << nrQubits
                                  << " qubits" << std::endl;

                        std::cout << "Probability for the different states: " << std::norm(regState[s]) << " vs "
                                  << std::norm(mpsState[s]) << std::endl;

                        std::cout << "Reg state:\n" << regState << std::endl;
                        std::cout << "Reg state normalization: " << regState.norm() << std::endl;

                        std::cout << "MPS state:\n" << mpsState << std::endl;
                        std::cout << "MPS state normalization: " << mpsState.norm() << std::endl;

                        std::cout << std::endl;
                        for (int q = 0; q < nrQubits; ++q)
                            std::cout << "Qubit " << q << " reg probability: " << reg.GetQubitProbability(q)
                                      << " vs mps: " << mps.GetProbability(q, false) << std::endl;

                        std::cout << std::endl;
                        for (int state = 0; state < regState.size(); ++state)
                            std::cout << "State " << state << " reg probability: " << std::norm(regState[state])
                                      << " vs mps: " << std::norm(mpsState[state]) << std::endl;

                        return false;
                    }
                }
            }

#ifdef _DEBUG
            std::cout << "Test passed: " << t << std::endl;
#endif
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool NumericalRankStabilityTestMPS()
{
    std::cout << "\nMPS simulator numerical rank stability test" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    // This fixed seed contains a mapped circuit that creates roundoff-only
    // Schmidt values and subsequently routes gates across them. With the Vidal
    // form (dividing by the lambdas) and without a numerical rank floor, the
    // pseudoinverse amplified that SVD noise into percent-level state-vector errors.
    // The implementation does not divide by the lambdas anymore (Hastings' method),
    // but this remains a useful precision regression test.
    std::mt19937 stabilityGenerator(544);
    std::bernoulli_distribution stabilityBool;
    std::uniform_int_distribution nrGatesDistribution(25, 50);
    std::uniform_int_distribution gateDistribution(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistribution(0, nrQubits - 1);
        std::uniform_int_distribution secondQubitDistribution(0, nrQubits - 2);

        for (int trial = 0; trial < 10; ++trial)
        {
            QC::TensorNetworks::MPSSimulator mps(nrQubits);
            QC::QubitRegister<> reg(nrQubits);
            const int nrGates = nrGatesDistribution(stabilityGenerator);

            for (int i = 0; i < nrGates; ++i)
            {
                const int gate = gateDistribution(stabilityGenerator);
                const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
                int qubit1 = qubitDistribution(stabilityGenerator);
                int qubit2 = (qubit1 + 1 + secondQubitDistribution(stabilityGenerator)) % nrQubits;
                if (twoQubitsGate && stabilityBool(stabilityGenerator))
                    std::swap(qubit1, qubit2);

                mps.ApplyGate(*gates[gate], qubit1, qubit2);
                reg.ApplyGate(*gates[gate], qubit1, qubit2);

                const auto &reference = reg.getRegisterStorage();
                const auto actual = mps.getRegisterStorage();
                for (Eigen::Index state = 0; state < reference.size(); ++state)
                    if (!approxEqual(reference[state], actual[state], 1E-10))
                    {
                        std::cout << "MPS numerical-rank filtering did not prevent roundoff amplification for "
                                  << nrQubits << " qubits" << std::endl;
                        return false;
                    }
            }
        }
    }

    std::cout << "Success" << std::endl;
    return true;
}

void PrintMeasurements(int nrMeasurements, const std::unordered_map<std::vector<bool>, int> &measurements)
{
    for (const auto &[k, v] : measurements)
    {
        for (const auto &b : k)
            std::cout << b << " ";
        std::cout << " : " << static_cast<double>(v) / nrMeasurements << std::endl;
    }
}

bool CheckMeasurements(int nrQubits, int nrMeasurements, std::unordered_map<std::vector<bool>, int> &measurementsRegMap,
                       std::unordered_map<std::vector<bool>, int> &measurementsMPSMap,
                       std::unordered_map<std::vector<bool>, int> &measurementsMPSMapOpt,
                       std::unordered_map<std::vector<bool>, int> &measurementsMPSMapAll)
{
    static const double threshold = 0.05;

    for (const auto &[key, value] : measurementsRegMap)
    {
        const double dif1 = abs((static_cast<double>(measurementsMPSMap[key] - value)) / nrMeasurements);
        const double dif2 = abs((static_cast<double>(measurementsMPSMapOpt[key] - value)) / nrMeasurements);
        const double dif3 = abs((static_cast<double>(measurementsMPSMapAll[key] - value)) / nrMeasurements);
        if (dif1 > threshold || dif2 > threshold || dif3 > threshold)
        {
            std::cout << "Measurements test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;
            std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
            std::cout << "Difference 1: " << dif1 << " Difference 2: " << dif2 << " Difference 3: " << dif3
                      << std::endl;

            std::cout << "Reg measurements:\n";
            PrintMeasurements(nrMeasurements, measurementsRegMap);

            if (dif1 > threshold)
            {
                std::cout << "MPS measurements:\n";
                PrintMeasurements(nrMeasurements, measurementsMPSMap);
            }

            if (dif2 > threshold)
            {
                std::cout << "MPS optimized measurements:\n";
                PrintMeasurements(nrMeasurements, measurementsMPSMapOpt);
            }

            if (dif3 > threshold)
            {
                std::cout << "MPS 'no collapse' measurements:\n";
                PrintMeasurements(nrMeasurements, measurementsMPSMapAll);
            }

            return false;
        }
    }

    return true;
}

bool TestMappedMeasurementsWithOneAndTwoQubitGatesCircuits()
{
    std::cout << "\nMPS swapping/mapped simulator measurements test with circuits with both one and two qubit gates"
              << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    const int nrMeasurements = 10000;

    for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        for (int t = 0; t < 5; ++t)
        {
            const auto circuit = GenerateRandomCircuitWithGatesNoAdjacent(gates, 25, 50, nrQubits);
            for (int c = 0; c < 3; ++c)
            {
                std::unordered_map<std::vector<bool>, int> measurementsRegMap;
                std::unordered_map<std::vector<bool>, int> measurementsMPSMap;
                std::unordered_map<std::vector<bool>, int> measurementsMPSMapOpt;
                std::unordered_map<std::vector<bool>, int> measurementsMPSMapAll;

                size_t remainingCounts = nrMeasurements;
                size_t nrThreads = QC::QubitRegisterCalculator<>::GetNumberOfThreads();
                nrThreads = std::min(nrThreads, std::max<size_t>(remainingCounts, 1ULL));

                std::vector<std::future<void>> tasks(nrThreads);

                const size_t cntPerThread = static_cast<size_t>(ceil(static_cast<double>(remainingCounts) / nrThreads));

                std::mutex resultsMutex;

                for (size_t i = 0; i < nrThreads; ++i)
                {
                    const size_t curCnt = std::min(cntPerThread, remainingCounts);
                    remainingCounts -= curCnt;

                    tasks[i] = std::async(std::launch::async, [&circuit, &measurementsRegMap, &measurementsMPSMap,
                                                               &measurementsMPSMapOpt, &measurementsMPSMapAll, curCnt,
                                                               nrQubits, &resultsMutex]() {
                        QC::TensorNetworks::MPSSimulator mps(nrQubits);
                        QC::TensorNetworks::MPSSimulator mpsOpt(nrQubits);
                        QC::TensorNetworks::MPSSimulator mpsAll(nrQubits);
                        QC::QubitRegister reg(nrQubits);

                        std::vector<bool> measurementsReg(nrQubits);
                        std::vector<bool> measurementsMPS(nrQubits);
                        std::vector<bool> measurementsMPSOpt(nrQubits);

                        std::set<Eigen::Index> qubits;
                        for (int q = 0; q < nrQubits; ++q)
                            qubits.insert(q);

                        for (size_t j = 0; j < curCnt; ++j)
                        {
                            for (const auto &gate : circuit)
                            {
                                mps.ApplyGate(*gate);
                                mpsOpt.ApplyGate(*gate);
                                if (j == 0)
                                    mpsAll.ApplyGate(
                                        *gate); // apply it only the first iteration, measure no-collapse is used
                                reg.ApplyGate(*gate);
                            }

                            auto measRes = mpsOpt.MeasureQubits(qubits);
                            for (int q = 0; q < nrQubits; ++q)
                            {
                                measurementsReg[q] = reg.MeasureQubit(q);
                                measurementsMPS[q] = mps.MeasureQubit(q);
                                measurementsMPSOpt[q] = measRes[q];
                            }
                            auto measurementsMPSAll = mpsAll.MeasureNoCollapse();
                            std::vector<bool> measurementsMPSAllVec(nrQubits);
                            for (const auto &[qubit, val] : measurementsMPSAll)
                                measurementsMPSAllVec[qubit] = val;

                            {
                                const std::lock_guard lock(resultsMutex);
                                ++measurementsRegMap[measurementsReg];
                                ++measurementsMPSMap[measurementsMPS];
                                ++measurementsMPSMapOpt[measurementsMPSOpt];
                                ++measurementsMPSMapAll[measurementsMPSAllVec];
                            }

                            reg.setToBasisState(0);
                            mps.setToBasisState(0);
                            mpsOpt.setToBasisState(0);
                            // mpsAll.setToBasisState(0); // no need to reset if only no-collapse measurements are used
                            // and the gates are applied only the first iteration
                        }
                    });
                }

                for (size_t i = 0; i < nrThreads; ++i)
                    tasks[i].get();

                std::cout << "." << std::flush;

                const bool res = CheckMeasurements(nrQubits, nrMeasurements, measurementsRegMap, measurementsMPSMap,
                                                   measurementsMPSMapOpt, measurementsMPSMapAll);

                if (!res)
                    return false;
            }
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

size_t GenerateRandomState(size_t nrQubits)
{
    assert(nrQubits > 0);

    size_t state = dist_bool(gen) ? 1 : 0;
    for (size_t i = 1; i < nrQubits; ++i)
    {
        state <<= 1;
        state |= dist_bool(gen) ? 1 : 0;
    }

    return state;
}

bool OneAndTwoQubitGatesTestMappedRandomAmplitudes()
{
    std::cout << "\nMPS swapping/mapped simulator random state test for both one and two qubit gates" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(25, 50);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        for (int t = 0; t < 10; ++t)
        {
#ifdef _DEBUG
            std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl
                      << std::endl
                      << std::endl;
#endif

            QC::TensorNetworks::MPSSimulator mps(nrQubits);
            QC::QubitRegister reg(nrQubits);

            const int lim = nrGatesDistr(gen);

            for (int i = 0; i < lim; ++i)
            {
                const int gate = gateDistr(gen);
                const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
                int qubit1 = qubitDistr(gen);
                int qubit2 = (qubit1 + 1 + qubitDistr2(gen)) % nrQubits;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

#ifdef _DEBUG
                if (twoQubitsGate)
                    std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2
                              << std::endl;
#endif

                mps.ApplyGate(*gates[gate], qubit1, qubit2);
                reg.ApplyGate(*gates[gate], qubit1, qubit2);

                for (int s = 0; s < 5; ++s)
                {
                    const size_t state = GenerateRandomState(nrQubits);

                    const auto regAmpl = reg.getBasisStateAmplitude(state);
                    const auto mpsAmpl = mps.getBasisStateAmplitude(state);

                    if (!approxEqual(regAmpl, mpsAmpl, 1E-3))
                    {
                        std::cout << "State " << state << " simulation test failed for the MPS simulator for "
                                  << nrQubits << " qubits" << std::endl;

                        std::cout << "Amplitude for the different states: " << regAmpl << " vs " << mpsAmpl
                                  << std::endl;

                        return false;
                    }

                    const auto regProb = reg.getBasisStateProbability(state);
                    const auto mpsProb = mps.getBasisStateProbability(state);

                    if (!approxEqual(regProb, mpsProb, 1E-7))
                    {
                        std::cout << "State " << state << " simulation test failed for the MPS simulator for "
                                  << nrQubits << " qubits" << std::endl;

                        std::cout << "Probability for the different states: " << regProb << " vs " << mpsProb
                                  << std::endl;

                        return false;
                    }
                }
            }

#ifdef _DEBUG
            std::cout << "Test passed: " << t << std::endl;
#endif
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool CheckQubitsProbability()
{
    std::cout << "\nMPS swapping/mapped simulator qubits probability test for both one and two qubit gates"
              << std::endl;
    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        for (int t = 0; t < 10; ++t)
        {
#ifdef _DEBUG
            std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl
                      << std::endl
                      << std::endl;
#endif

            QC::TensorNetworks::MPSSimulator mps(nrQubits);
            QC::QubitRegister reg(nrQubits);

            const auto circ = GenerateRandomCircuitWithGatesNoAdjacent(gates, 50, 150, nrQubits);
            for (const auto &gate : circ)
            {
                mps.ApplyGate(*gate);
                reg.ApplyGate(*gate);
            }

            // check qubits probabilities
            for (int q = 0; q < nrQubits; ++q)
            {
                const auto regProb = reg.GetQubitProbability(q);
                const auto mpsProb = mps.GetProbability(q, false);
                if (!approxEqual(regProb, mpsProb, 1E-3))
                {
                    std::cout << "Qubit " << q << " probability test failed for the MPS simulator for " << nrQubits
                              << " qubits" << std::endl;
                    std::cout << "Probability for the qubit in statevector: " << regProb << " vs mps: " << mpsProb
                              << std::endl;
                    return false;
                }
            }
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool StateSimulationTest()
{
    // for longer tests:
    /*
    for (int i = 0; i < 30; ++i)
    {
        if (!(OneAndTwoQubitGatesTestMappedRandomAmplitudes() && OneAndTwoQubitGatesTestMapped() &&
    TestMappedMeasurementsWithOneAndTwoQubitGatesCircuits())) return false;
    }
    */

    return /*OneQubitGatesTest() && OneAndTwoQubitGatesTest() &&
        TestMeasurementsWithOneQubitGatesCircuits() && TestMeasurementsWithOneAndTwoQubitGatesCircuits() &&*/
        OneAndTwoQubitGatesTestMappedRandomAmplitudes() && OneAndTwoQubitGatesTestMapped() &&
        TestMappedMeasurementsWithOneAndTwoQubitGatesCircuits() && CheckQubitsProbability();
}

bool matchRandomExpectationValues(QC::QubitRegister<> &reg, QC::TensorNetworks::MPSSimulator &mps)
{
    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gatesOpExp;
    FillOneQubitGates(gatesOpExp);

    const auto nrQubits = reg.getNrQubits();
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gatesOpExp.size()) - 1);
    std::uniform_int_distribution nrOpsDistr(1, static_cast<int>(nrQubits));

    std::uniform_int_distribution qbitDistr(0, static_cast<int>(nrQubits) - 1);

    for (int s = 0; s < 5; ++s)
    {
        // generate a random Pauli string and check its expectation value
        std::vector<QC::Gates::AppliedGate<>> expGates;
        expGates.reserve(nrQubits);
        std::string pauliString;

        const auto nrOps = nrOpsDistr(gen);
        for (int q = 0; q < nrOps; ++q)
        {
            auto gate = gateDistr(gen);
            const auto qubit = qbitDistr(gen);

            expGates.emplace_back(gatesOpExp[gate]->getRawOperatorMatrix(), qubit);
        }

        const auto exp1 = reg.ExpectationValue(expGates);
        const auto exp2 = mps.ExpectationValue(expGates);

        if (!approxEqual(exp1, exp2, 1E-7))
        {
            std::cout << std::endl
                      << "Expectation values are not equal for MPS and statevector simulator for " << nrQubits
                      << " qubits, values: " << exp1 << ", " << exp2 << std::endl;

            return false;
        }
    }

    return true;
}

bool checkExpectationValuesMPS()
{
    std::cout << "\nTesting MPS expectation values against statevector..." << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(25, 35);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        for (int t = 0; t < 10; ++t)
        {
#ifdef _DEBUG
            std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl
                      << std::endl
                      << std::endl;
#endif

            QC::TensorNetworks::MPSSimulator mps(nrQubits);
            QC::QubitRegister reg(nrQubits);

            const int lim = nrGatesDistr(gen);

            for (int i = 0; i < lim; ++i)
            {
                const int gate = gateDistr(gen);
                const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
                int qubit1 = qubitDistr(gen);
                int qubit2 = (qubit1 + 1 + qubitDistr2(gen)) % nrQubits;
                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

#ifdef _DEBUG
                if (twoQubitsGate)
                    std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2
                              << std::endl;
#endif

                mps.ApplyGate(*gates[gate], qubit1, qubit2);
                reg.ApplyGate(*gates[gate], qubit1, qubit2);
            }

            if (!matchRandomExpectationValues(reg, mps))
                return false;
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool TrimTestMPS()
{
    std::cout << "\nMPS simulator Trim test" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    std::uniform_int_distribution nrGatesDistr(20, 40);
    std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

    // a no-op two qubit gate: applying it with truncation enabled on a bond does exactly what Trim
    // does on that bond (contract the two neighbour sites, apply no gate, SVD with truncation)
    const QC::Gates::TwoQubitsGate<> identityTwoQubitGate(Eigen::MatrixXcd::Identity(4, 4));

    for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
        std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

        for (int t = 0; t < 5; ++t)
        {
            // build a random circuit so it can be replayed identically on several simulators
            std::vector<QC::Gates::AppliedGate<>> circuit;
            const int lim = nrGatesDistr(gen);
            for (int i = 0; i < lim; ++i)
            {
                const int gate = gateDistr(gen);
                const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
                int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
                int qubit2 = twoQubitsGate ? qubit1 + 1 : qubit1;

                if (twoQubitsGate && dist_bool(gen))
                    std::swap(qubit1, qubit2);

                circuit.emplace_back(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
            }

            const int chi = 1 + (t % 4); // trim down to a bond dimension between 1 and 4

            // reference: apply the circuit without any limit (exact MPS), then lower the limit and
            // truncate with the no-op two qubit gate exactly the bonds Trim would touch, that is
            // only the ones that exceed the (newly lowered) limit. Trimming a bond that is already
            // within the limit is a state preserving re-gauging that would still change the result
            // of subsequent truncations, so the reference must skip those bonds just like Trim does.
            QC::TensorNetworks::MPSSimulatorImpl mpsRef(nrQubits);
            for (const auto &g : circuit)
                mpsRef.ApplyGate(g);
            mpsRef.setLimitBondDimension(chi);
            const auto refBondDims = mpsRef.getBondDimensions();
            for (int q = 0; q < nrQubits - 1; ++q)
                if (refBondDims[q] > chi)
                    mpsRef.ApplyGate(identityTwoQubitGate, q, q + 1);

            // the one under test: apply the same circuit without any limit, then lower the limit and Trim
            QC::TensorNetworks::MPSSimulatorImpl mpsTrim(nrQubits);
            for (const auto &g : circuit)
                mpsTrim.ApplyGate(g);
            mpsTrim.setLimitBondDimension(chi);
            mpsTrim.Trim();

            // after trimming, every bond dimension must be within the (newly lowered) limit
            for (const auto bd : mpsTrim.getBondDimensions())
                if (bd > chi)
                {
                    std::cout << "Trim did not reduce a bond dimension to the limit (" << bd << " > " << chi << ") for "
                              << nrQubits << " qubits" << std::endl;
                    return false;
                }

            // Trim must produce the same state as the equivalent no-op two qubit gate truncations
            const auto refState = mpsRef.getRegisterStorage();
            const auto trimState = mpsTrim.getRegisterStorage();

            for (int s = 0; s < refState.size(); ++s)
                if (!approxEqual(refState[s], trimState[s], 1E-3))
                {
                    std::cout << "Trim state " << s << " differs from the reference truncation for " << nrQubits
                              << " qubits" << std::endl;
                    std::cout << "Reference: " << refState[s] << " vs Trim: " << trimState[s] << std::endl;
                    return false;
                }

            // trimming again should be a no-op since all bonds are already within the limit
            mpsTrim.Trim();
            const auto trimState2 = mpsTrim.getRegisterStorage();
            for (int s = 0; s < trimState.size(); ++s)
                if (!approxEqual(trimState[s], trimState2[s], 1E-3))
                {
                    std::cout << "Trim is not idempotent for " << nrQubits << " qubits" << std::endl;
                    return false;
                }
        }
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}

bool TruncationModeTestMPS()
{
    using TruncationMode = QC::TensorNetworks::MPSSimulatorInterface::TruncationMode;
    using MPSState = QC::TensorNetworks::MPSSimulatorBaseState;

    std::cout << "\nMPS simulator truncation mode test" << std::endl;

    // The default must be DiscardedWeight (Qiskit Aer's / ITensor's convention) - a deliberate
    // default-behavior change from what every earlier version of this simulator did.
    {
        QC::TensorNetworks::MPSSimulatorImpl defaultMps(2);
        if (defaultMps.getTruncationMode() != TruncationMode::DiscardedWeight)
        {
            std::cout << "Default truncation mode is not DiscardedWeight" << std::endl;
            return false;
        }
    }

    // getter/setter round trip for both modes
    {
        QC::TensorNetworks::MPSSimulatorImpl mps(2);
        if (!mps.setTruncationMode(TruncationMode::RelativeToMax) ||
            mps.getTruncationMode() != TruncationMode::RelativeToMax)
        {
            std::cout << "Truncation mode round trip failed for RelativeToMax" << std::endl;
            return false;
        }

        if (!mps.setTruncationMode(TruncationMode::DiscardedWeight) ||
            mps.getTruncationMode() != TruncationMode::DiscardedWeight)
        {
            std::cout << "Truncation mode round trip failed for DiscardedWeight" << std::endl;
            return false;
        }
    }

    // Build two Bell pairs (qubits 0-1 and 2-3), each with two well-separated singular values
    // (1/sqrt2, 1/sqrt2) far above any threshold used below, then join qubits 1 and 2 with a
    // generic two qubit unitary. Because the Bell-pair bonds are never close to the threshold,
    // they are truncated identically (untouched) in every run below, so the joining gate's theta
    // matrix - and therefore its raw SVD spectrum - is identical between the exact reference and
    // every truncation-mode run. This isolates the truncation SELECTION step (the only thing this
    // test wants to exercise) from the SVD itself.
    const QC::Gates::HadamardGate<> hGate;
    const QC::Gates::CNOTGate<> cnotGate;
    // A fixed, hand-picked (not block-diagonal / not control-structured) complex matrix, turned
    // into a genuinely generic unitary via QR - deliberately avoiding gates built from
    // "controlled-X" structure, which turned out to leave this particular Bell-pair setup at
    // Schmidt rank 2 regardless of rotation angle.
    Eigen::MatrixXcd seed(4, 4);
    seed << std::complex<double>(0.3, 0.7), std::complex<double>(-0.5, 0.2), std::complex<double>(0.4, -0.6),
        std::complex<double>(0.1, 0.9), std::complex<double>(-0.8, 0.1), std::complex<double>(0.6, 0.5),
        std::complex<double>(0.2, 0.3), std::complex<double>(-0.4, 0.7), std::complex<double>(0.5, -0.3),
        std::complex<double>(0.9, -0.1), std::complex<double>(-0.6, 0.4), std::complex<double>(0.2, 0.2),
        std::complex<double>(-0.2, 0.6), std::complex<double>(0.3, -0.8), std::complex<double>(0.7, 0.1),
        std::complex<double>(-0.5, -0.3);
    const Eigen::MatrixXcd joinGate = Eigen::HouseholderQR<Eigen::MatrixXcd>(seed).householderQ();

    auto buildCircuit = [&](QC::TensorNetworks::MPSSimulatorImpl &mps) {
        mps.ApplyGate(hGate, 0);
        mps.ApplyGate(cnotGate, 1, 0);
        mps.ApplyGate(hGate, 2);
        mps.ApplyGate(cnotGate, 3, 2);
        mps.ApplyGate(QC::Gates::AppliedGate<>(joinGate, 2, 1));
    };

    constexpr Eigen::Index bondIndex = 1; // the bond between qubit 1 and qubit 2, joined last above

    QC::TensorNetworks::MPSSimulatorImpl mpsRef(4);
    buildCircuit(mpsRef);

    const auto refState = std::static_pointer_cast<MPSState>(mpsRef.getState());
    const Eigen::VectorXd spectrum = refState->lambdas[bondIndex]; // descending-sorted, per Eigen's SVD

    if (spectrum.size() < 3)
    {
        std::cout << "Test setup did not produce a rich enough spectrum (" << spectrum.size() << " singular values)"
                  << std::endl;
        return false;
    }

    const double threshold = 0.3;

    // Independently compute, from the exact spectrum, how many singular values each mode should
    // keep - deliberately not calling MPSSimulatorImpl::ComputeCompressedRank, so this checks the
    // SPECIFICATION (as documented on MPSSimulatorInterface::TruncationMode), not just that the
    // implementation agrees with itself.
    Eigen::Index expectedRelative;
    {
        const double cutoffValue = threshold * spectrum[0];
        Eigen::Index i = spectrum.size() - 1;
        while (i >= 0 && spectrum[i] < cutoffValue)
            --i;
        expectedRelative = std::max<Eigen::Index>(i + 1, 1);
    }

    Eigen::Index expectedDiscardedWeight;
    {
        const double total = spectrum.squaredNorm();
        double discarded = 0.;
        Eigen::Index keep = spectrum.size();
        for (Eigen::Index i = spectrum.size() - 1; i > 0; --i)
        {
            const double sq = spectrum[i] * spectrum[i];
            if ((discarded + sq) / total >= threshold)
                break;

            discarded += sq;
            keep = i;
        }
        expectedDiscardedWeight = keep;
    }

    if (expectedRelative == expectedDiscardedWeight)
    {
        std::cout << "Test setup does not discriminate between the two truncation modes at threshold " << threshold
                  << " (both would keep " << expectedRelative << ")" << std::endl;
        return false;
    }

    // RelativeToMax explicitly requested must match the independently-computed expectation - this
    // is the regression test for "RelativeToMax still reproduces the original behavior".
    {
        QC::TensorNetworks::MPSSimulatorImpl mpsRel(4);
        mpsRel.setTruncationMode(TruncationMode::RelativeToMax);
        mpsRel.setLimitEntanglement(threshold);
        buildCircuit(mpsRel);

        const auto bondDims = mpsRel.getBondDimensions();
        if (bondDims[bondIndex] != expectedRelative)
        {
            std::cout << "RelativeToMax kept " << bondDims[bondIndex] << " singular values, expected "
                      << expectedRelative << std::endl;
            return false;
        }
    }

    // DiscardedWeight explicitly requested must match the independently-computed expectation.
    {
        QC::TensorNetworks::MPSSimulatorImpl mpsWeight(4);
        mpsWeight.setTruncationMode(TruncationMode::DiscardedWeight);
        mpsWeight.setLimitEntanglement(threshold);
        buildCircuit(mpsWeight);

        const auto bondDims = mpsWeight.getBondDimensions();
        if (bondDims[bondIndex] != expectedDiscardedWeight)
        {
            std::cout << "DiscardedWeight kept " << bondDims[bondIndex] << " singular values, expected "
                      << expectedDiscardedWeight << std::endl;
            return false;
        }
    }

    // No explicit mode set: must match DiscardedWeight (the default), not RelativeToMax - this is
    // the regression test guarding the default-flip decision.
    {
        QC::TensorNetworks::MPSSimulatorImpl mpsDefault(4);
        mpsDefault.setLimitEntanglement(threshold);
        buildCircuit(mpsDefault);

        const auto bondDims = mpsDefault.getBondDimensions();
        if (bondDims[bondIndex] != expectedDiscardedWeight)
        {
            std::cout << "Default truncation mode did not match DiscardedWeight (kept " << bondDims[bondIndex]
                      << ", expected " << expectedDiscardedWeight << ")" << std::endl;
            return false;
        }
    }

    std::cout << "Success" << std::endl;

    return true;
}

// Regression test for a real bug found during review: MPSSimulator::Clone() copies
// limitSize/limitEntanglement/chi/singularValueThreshold but, before this test existed,
// silently forgot to copy truncationMode - a clone would always come back at the default
// (DiscardedWeight) even if the original had explicitly opted into RelativeToMax.
bool CloneTestMPS()
{
    using TruncationMode = QC::TensorNetworks::MPSSimulatorInterface::TruncationMode;

    std::cout << "\nMPS simulator clone preserves truncation mode" << std::endl;

    QC::TensorNetworks::MPSSimulator mps(2);
    if (!mps.setTruncationMode(TruncationMode::RelativeToMax))
    {
        std::cout << "Failed to set RelativeToMax on the original simulator" << std::endl;
        return false;
    }

    const auto cloned = mps.Clone();
    if (cloned->getTruncationMode() != TruncationMode::RelativeToMax)
    {
        std::cout << "Clone did not preserve the RelativeToMax truncation mode (got DiscardedWeight instead)"
                  << std::endl;
        return false;
    }

    std::cout << "Success" << std::endl;

    return true;
}

static bool WideBasisInitializationTestMPS()
{
    std::cout << "\nMPS simulator - wide vector basis-state initialization" << std::endl;

    const size_t digits = std::numeric_limits<size_t>::digits;
    const size_t nrQubits = digits + 5;
    std::vector<bool> stateBits(nrQubits, false);
    stateBits[0] = true;
    stateBits[digits - 1] = true;
    stateBits[digits] = true;
    stateBits[nrQubits - 1] = true;
    const std::vector<bool> state = std::move(stateBits);

    auto checkSimulator = [&state, digits, nrQubits](QC::TensorNetworks::MPSSimulatorInterface &simulator,
                                                     const char *simulatorName) {
        simulator.setToBasisState(state);

        for (size_t q = 0; q < nrQubits; ++q)
        {
            const double expectedOneProbability = state[q] ? 1. : 0.;
            if (!approxEqual(simulator.GetProbability(static_cast<Eigen::Index>(q), false), expectedOneProbability,
                             1E-12))
            {
                std::cout << simulatorName << " initialized qubit " << q << " incorrectly" << std::endl;
                return false;
            }
        }

        std::vector<bool> matchingState = state;
        if (!approxEqual(simulator.getBasisStateAmplitude(matchingState), std::complex<double>(1., 0.), 1E-12))
        {
            std::cout << simulatorName << " did not initialize the requested wide basis state" << std::endl;
            return false;
        }

        std::vector<bool> differentState = state;
        differentState[digits] = !differentState[digits];
        if (!approxEqual(simulator.getBasisStateAmplitude(differentState), std::complex<double>(0., 0.), 1E-12))
        {
            std::cout << simulatorName << " initialized more than one basis amplitude" << std::endl;
            return false;
        }

        const std::vector<bool> shortState{true, false, true};
        simulator.setToBasisState(shortState);
        if (!approxEqual(simulator.GetProbability(0, false), 1., 1E-12) ||
            !approxEqual(simulator.GetProbability(1, false), 0., 1E-12) ||
            !approxEqual(simulator.GetProbability(2, false), 1., 1E-12) ||
            !approxEqual(simulator.GetProbability(static_cast<Eigen::Index>(nrQubits - 1), false), 0., 1E-12))
        {
            std::cout << simulatorName << " did not zero-extend a short basis vector" << std::endl;
            return false;
        }

        return true;
    };

    QC::TensorNetworks::MPSSimulatorImpl directSimulator(nrQubits);
    QC::TensorNetworks::MPSSimulator mappedSimulator(nrQubits);
    if (!checkSimulator(directSimulator, "Direct MPS simulator") ||
        !checkSimulator(mappedSimulator, "Mapped MPS simulator"))
        return false;

    const size_t allBits = std::numeric_limits<size_t>::max();
    for (const size_t scalarQubitCount : {digits, digits + 1})
    {
        auto checkScalarSimulator = [allBits, digits,
                                     scalarQubitCount](QC::TensorNetworks::MPSSimulatorInterface &simulator,
                                                       const char *simulatorName) {
            simulator.setToBasisState(allBits);

            for (size_t q = 0; q < scalarQubitCount; ++q)
            {
                const double expectedOneProbability = q < digits ? 1. : 0.;
                if (!approxEqual(simulator.GetProbability(static_cast<Eigen::Index>(q), false), expectedOneProbability,
                                 1E-12))
                {
                    std::cout << simulatorName << " initialized size_t bit " << q << " incorrectly for "
                              << scalarQubitCount << " qubits" << std::endl;
                    return false;
                }
            }

            return true;
        };

        QC::TensorNetworks::MPSSimulatorImpl directScalarSimulator(scalarQubitCount);
        QC::TensorNetworks::MPSSimulator mappedScalarSimulator(scalarQubitCount);
        if (!checkScalarSimulator(directScalarSimulator, "Direct MPS simulator") ||
            !checkScalarSimulator(mappedScalarSimulator, "Mapped MPS simulator"))
            return false;
    }

    QC::TensorNetworks::MPSSimulator mappedValidationSimulator(3);
    mappedValidationSimulator.setToBasisState(std::vector<bool>{true, false, true});
    mappedValidationSimulator.MoveAtBeginningOfChain({2});
    const auto stateBeforeInvalidCall =
        std::dynamic_pointer_cast<QC::TensorNetworks::MPSSimulatorState>(mappedValidationSimulator.getState());
    const std::vector<Eigen::Index> identityMap{0, 1, 2};
    if (!stateBeforeInvalidCall || stateBeforeInvalidCall->qubitsMap == identityMap)
    {
        std::cout << "MPS test setup did not create a non-identity logical qubit mapping" << std::endl;
        return false;
    }

    bool oversizedStateRejected = false;
    try
    {
        mappedValidationSimulator.setToBasisState(std::vector<bool>(4, false));
    }
    catch (const std::invalid_argument &)
    {
        oversizedStateRejected = true;
    }
    catch (...)
    {
    }

    const auto stateAfterInvalidCall =
        std::dynamic_pointer_cast<QC::TensorNetworks::MPSSimulatorState>(mappedValidationSimulator.getState());
    if (!oversizedStateRejected || !stateAfterInvalidCall ||
        stateAfterInvalidCall->qubitsMap != stateBeforeInvalidCall->qubitsMap ||
        stateAfterInvalidCall->qubitsMapInv != stateBeforeInvalidCall->qubitsMapInv ||
        !approxEqual(mappedValidationSimulator.GetProbability(0, false), 1., 1E-12) ||
        !approxEqual(mappedValidationSimulator.GetProbability(1, false), 0., 1E-12) ||
        !approxEqual(mappedValidationSimulator.GetProbability(2, false), 1., 1E-12))
    {
        std::cout << "Rejected MPS basis initialization changed the state or logical qubit mapping" << std::endl;
        return false;
    }

    std::cout << "Success" << std::endl;
    return true;
}

// Checks the canonical form the implementation relies on (the sites hold B = Gamma * lambda, see MPSSimulatorBase):
// for each site sum_s B_s B_s^dagger = I (right canonical) and sum_s B_s^dagger lambda_left^2 B_s = lambda_right^2
// (together they mean that the lambdas are the Schmidt values)
static bool CheckCanonicalFormMPS(const QC::TensorNetworks::MPSSimulatorImpl &mps, double tolerance,
                                  double &maxDeviation)
{
    const auto state = std::static_pointer_cast<QC::TensorNetworks::MPSSimulatorBaseState>(mps.getState());
    const auto &gammas = state->gammas;
    const auto &lambdas = state->lambdas;

    maxDeviation = 0;

    for (size_t q = 0; q < gammas.size(); ++q)
    {
        const auto &B = gammas[q];
        const Eigen::Index dimLeft = B.dimension(0);
        const Eigen::Index dimRight = B.dimension(2);

        const Eigen::VectorXd leftLambda2 =
            q == 0 ? Eigen::VectorXd::Ones(1) : Eigen::VectorXd(lambdas[q - 1].cwiseAbs2());
        const Eigen::VectorXd rightLambda2 =
            q + 1 == gammas.size() ? Eigen::VectorXd::Ones(1) : Eigen::VectorXd(lambdas[q].cwiseAbs2());

        if (leftLambda2.size() != dimLeft || rightLambda2.size() != dimRight)
        {
            std::cout << "Bond dimensions of site " << q << " don't match the lambdas" << std::endl;
            maxDeviation = std::numeric_limits<double>::infinity();
            return false;
        }

        Eigen::MatrixXcd rightCondition = Eigen::MatrixXcd::Zero(dimLeft, dimLeft);
        Eigen::MatrixXcd leftCondition = Eigen::MatrixXcd::Zero(dimRight, dimRight);

        for (int s = 0; s < 2; ++s)
        {
            Eigen::MatrixXcd Bs(dimLeft, dimRight);
            for (Eigen::Index j = 0; j < dimRight; ++j)
                for (Eigen::Index i = 0; i < dimLeft; ++i)
                    Bs(i, j) = B(i, s, j);

            rightCondition += Bs * Bs.adjoint();
            leftCondition += Bs.adjoint() * leftLambda2.asDiagonal() * Bs;
        }

        const Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(dimLeft, dimLeft);
        const Eigen::MatrixXcd expectedLeft = rightLambda2.cast<std::complex<double>>().asDiagonal();
        maxDeviation = std::max({maxDeviation, (rightCondition - identity).cwiseAbs().maxCoeff(),
                                 (leftCondition - expectedLeft).cwiseAbs().maxCoeff()});
    }

    return maxDeviation <= tolerance;
}

// These fixtures start from the B product, independently of the stored Schmidt metadata.
static std::shared_ptr<QC::TensorNetworks::MPSSimulatorBaseState> CorrelatedStateMPS(
    int half, const std::vector<double> &probabilities)
{
    using Impl = QC::TensorNetworks::MPSSimulatorImpl;
    auto state = std::make_shared<QC::TensorNetworks::MPSSimulatorBaseState>();
    for (int q = 0; q < 2 * half; ++q)
    {
        const Eigen::Index L = q < half ? Eigen::Index(1) << q : Eigen::Index(1) << (2 * half - q);
        const Eigen::Index R = q < half ? 2 * L : L / 2;
        Impl::GammaType b(L, 2, R);
        b.setZero();
        for (Eigen::Index l = 0; l < L; ++l)
            if (q < half)
                for (int p = 0; p < 2; ++p)
                {
                    const Eigen::Index r = l + p * L;
                    b(l, p, r) = q == half - 1 ? std::polar(std::sqrt(probabilities[r]), 0.17 * r) : 1.;
                }
            else
                b(l, l % 2, l / 2) = 1.;
        state->gammas.push_back(std::move(b));
        if (q + 1 < 2 * half)
            state->lambdas.push_back(Impl::LambdaType::Constant(R, 1. / std::sqrt(double(R))));
    }
    return state;
}

bool ReCanonicalizeRegressionTestMPS()
{
    using Impl = QC::TensorNetworks::MPSSimulatorImpl;
    using State = QC::TensorNetworks::MPSSimulatorBaseState;
    using Mode = Impl::TruncationMode;
    std::cout << "\nMPS simulator recanonicalization: stale weights, gauge repair, and compression" << std::endl;
    const auto require = [](bool condition, const char *message) {
        if (!condition)
            throw std::runtime_error(message);
    };
    try
    {
        for (int fixture = 0; fixture < 4; ++fixture)
        {
            auto state = std::make_shared<State>();
            state->gammas = {Impl::GammaType(1, 2, 2), Impl::GammaType(2, 2, 2), Impl::GammaType(2, 2, 1)};
            for (auto &b : state->gammas)
                b.setZero();
            for (int p = 0; p < 2; ++p)
            {
                state->gammas[0](0, p, p) = 1. / std::sqrt(2.);
                state->gammas[1](p, p, p) = 1.;
                state->gammas[2](p, p, 0) = 1.;
            }
            state->lambdas = {Impl::LambdaType::Constant(2, 1. / std::sqrt(2.)),
                              Impl::LambdaType::Constant(2, 1. / std::sqrt(2.))};
            if (fixture == 1 || fixture == 2)
                state->lambdas[0] << 1., fixture == 1 ? 1E-20 : 0.;
            if (fixture == 3)
            {
                state->gammas[0](0, 0, 0) *= 1E-20;
                state->gammas[1](0, 0, 0) *= 1E20;
            }
            Impl sim(3);
            sim.SetMultithreading(false);
            sim.setState(state);
            const Eigen::VectorXcd expected = sim.getRegisterStorage().normalized();
            for (int repeat = 0; repeat < 2; ++repeat)
            {
                sim.ReCanonicalize();
                double deviation = 0.;
                require((sim.getRegisterStorage() - expected).norm() < 1E-12, "Gauge repair changed the GHZ state");
                require(CheckCanonicalFormMPS(sim, 1E-12, deviation), "GHZ canonical form was not restored");
            }
        }

        // Rank 16 exercises BDCSVD; complex amplitudes check conjugation and tensor layouts.
        for (int half : {1, 4})
        {
            const int rank = 1 << half;
            std::vector<double> weights(rank);
            for (int i = 0; i < rank; ++i)
                weights[i] = 2. * (i + 1) / (rank * (rank + 1));
            Impl sim(2 * half);
            sim.setState(CorrelatedStateMPS(half, weights));
            const Eigen::VectorXcd expected = sim.getRegisterStorage().normalized();
            sim.ReCanonicalize();
            double deviation = 0.;
            require((sim.getRegisterStorage() - expected).norm() < 1E-11, "Repair changed a complex correlated state");
            require(CheckCanonicalFormMPS(sim, 1E-11, deviation), "Correlated state canonical form was not restored");
            require(sim.getBondDimensions()[half - 1] == rank, "Repair lost a nonzero Schmidt direction");
        }

        struct Case
        {
            Mode mode;
            int cap;
            double threshold;
            int expectedRank;
        };

        const Case cases[] = {{Mode::RelativeToMax, 0, -1., 4},    {Mode::RelativeToMax, 2, -1., 2},
                              {Mode::RelativeToMax, 0, .5, 2},     {Mode::RelativeToMax, 3, .5, 2},
                              {Mode::RelativeToMax, 2, .2, 2},     {Mode::DiscardedWeight, 0, .12, 2},
                              {Mode::DiscardedWeight, 3, .12, 2},  {Mode::DiscardedWeight, 2, .02, 2},
                              {Mode::DiscardedWeight, 1, .001, 1}, {Mode::DiscardedWeight, 0, 0., 4}};
        for (const auto &c : cases)
        {
            Impl sim(4);
            sim.setState(CorrelatedStateMPS(2, {.64, .25, .10, .01}));
            Eigen::VectorXcd expected = sim.getRegisterStorage();
            for (int i = c.expectedRank; i < 4; ++i)
                expected[i + 4 * i] = 0.;
            expected.normalize();
            sim.setTruncationMode(c.mode);
            if (c.cap)
                sim.setLimitBondDimension(c.cap);
            if (c.threshold >= 0.)
                sim.setLimitEntanglement(c.threshold);
            sim.ReCanonicalize();
            double deviation = 0.;
            require(sim.getBondDimensions()[1] == c.expectedRank, "The stronger compression setting did not win");
            for (const auto dim : sim.getBondDimensions())
                require(!c.cap || dim <= c.cap, "Canonical refresh regrew a bond beyond the cap");
            require((sim.getRegisterStorage() - expected).norm() < 1E-11,
                    "Compression retained the wrong Schmidt components");
            require(CheckCanonicalFormMPS(sim, 1E-11, deviation), "Compression left stale canonical metadata");
            if (c.mode == Mode::RelativeToMax && c.threshold >= 0.)
            {
                const auto compressed = std::static_pointer_cast<State>(sim.getState());
                for (const auto &lambda : compressed->lambdas)
                    require(lambda.size() == 1 || lambda.tail(1)[0] >= c.threshold * lambda[0] - 1E-12,
                            "Final Schmidt spectrum violates the relative cutoff");
            }
        }

        // One application discards .04, keeping rank 3. Reapplying .06 during refresh
        // would discard a further .04/.96 and incorrectly reduce the rank to 2.
        // The large components occupy 00 and 11, keeping the outer bonds balanced.
        Impl budget(4);
        budget.setState(CorrelatedStateMPS(2, {.46, .04, .04, .46}));
        budget.setLimitEntanglement(.06);
        budget.ReCanonicalize();
        require(budget.getBondDimensions()[1] == 3, "Refresh repeatedly spent the discarded-weight budget");
        double deviation = 0.;
        require(CheckCanonicalFormMPS(budget, 1E-11, deviation), "Budget test left stale canonical metadata");
        require(std::abs(budget.getRegisterStorage().squaredNorm() - 1.) < 1E-12,
                "Compression did not normalize the state");

        QC::TensorNetworks::MPSSimulator wrapped(3);
        wrapped.ApplyGate(QC::Gates::HadamardGate<>(), 0);
        wrapped.ApplyGate(QC::Gates::CNOTGate<>(), 2, 0);
        wrapped.setLimitBondDimension(1);
        wrapped.ReCanonicalize();
        for (const auto dim : wrapped.getBondDimensions())
            require(dim == 1, "Wrapper did not enforce the cap");

        Impl single(1);
        single.ApplyGate(QC::Gates::HadamardGate<>(), 0);
        const Eigen::VectorXcd before = single.getRegisterStorage();
        single.setLimitBondDimension(1);
        single.setLimitEntanglement(.9);
        single.ReCanonicalize();
        require((before - single.getRegisterStorage()).norm() < 1E-12, "Single-site repair changed the state");
    }
    catch (const std::exception &error)
    {
        std::cout << error.what() << std::endl;
        return false;
    }
    std::cout << "Success" << std::endl;
    return true;
}

bool CanonicalFormTestMPS()
{
    std::cout << "\nMPS simulator canonical form test (after gates, measurements and ReCanonicalize)" << std::endl;

    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    constexpr double tolerance = 1E-10;

    for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
    {
        std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

        for (int t = 0; t < 5; ++t)
        {
            double deviation = 0;

            // exact simulation: the canonical form must hold after the gates
            QC::TensorNetworks::MPSSimulatorImpl mps(nrQubits);
            for (const auto &gate : GenerateRandomCircuitWithGates(gates, 25, 50, nrQubits))
                mps.ApplyGate(*gate);

            if (!CheckCanonicalFormMPS(mps, tolerance, deviation))
            {
                std::cout << "Canonical form broken after applying gates for " << nrQubits
                          << " qubits, deviation: " << deviation << std::endl;
                return false;
            }

            // ... and after measurements (the collapse is propagated along the chain)
            mps.MeasureQubit(qubitDistr(gen));
            if (!CheckCanonicalFormMPS(mps, tolerance, deviation))
            {
                std::cout << "Canonical form broken after a measurement for " << nrQubits
                          << " qubits, deviation: " << deviation << std::endl;
                return false;
            }

            std::set<Eigen::Index> measuredQubits;
            for (int q = 0; q < nrQubits; q += 2)
                measuredQubits.insert(q);
            for (const auto &gate : GenerateRandomCircuitWithGates(gates, 10, 20, nrQubits))
                mps.ApplyGate(*gate);
            mps.MeasureQubits(measuredQubits);
            if (!CheckCanonicalFormMPS(mps, tolerance, deviation))
            {
                std::cout << "Canonical form broken after measuring several qubits for " << nrQubits
                          << " qubits, deviation: " << deviation << std::endl;
                return false;
            }

            // truncated simulation: the canonical form holds only approximately, ReCanonicalize must restore it exactly
            // without changing the (normalized) state
            QC::TensorNetworks::MPSSimulatorImpl mpsTrunc(nrQubits);
            mpsTrunc.setLimitBondDimension(2);
            for (const auto &gate : GenerateRandomCircuitWithGates(gates, 25, 50, nrQubits))
                mpsTrunc.ApplyGate(*gate);

            const Eigen::VectorXcd stateBefore = mpsTrunc.getRegisterStorage().normalized();
            mpsTrunc.ReCanonicalize();
            if (!CheckCanonicalFormMPS(mpsTrunc, tolerance, deviation))
            {
                std::cout << "ReCanonicalize did not restore the canonical form for " << nrQubits
                          << " qubits, deviation: " << deviation << std::endl;
                return false;
            }

            const Eigen::VectorXcd stateAfter = mpsTrunc.getRegisterStorage();
            for (Eigen::Index s = 0; s < stateBefore.size(); ++s)
                if (!approxEqual(stateBefore[s], stateAfter[s], tolerance))
                {
                    std::cout << "ReCanonicalize changed the truncated state for " << nrQubits << " qubits"
                              << std::endl;
                    return false;
                }
        }
    }

    std::cout << "Success" << std::endl;
    return true;
}

bool MultithreadingSettingTestMPS()
{
    std::cout << "\nMPS simulator multithreading setting test" << std::endl;

    QC::TensorNetworks::MPSSimulatorImpl impl(2);
    QC::TensorNetworks::MPSSimulator mps(2);
    if (!impl.GetMultithreading() || !mps.GetMultithreading())
    {
        std::cout << "Multithreading is not enabled by default" << std::endl;
        return false;
    }

    impl.SetMultithreading(false);
    mps.SetMultithreading(false);
    if (impl.GetMultithreading() || mps.GetMultithreading() || mps.Clone()->GetMultithreading())
    {
        std::cout << "Disabling multithreading did not stick (or was not cloned)" << std::endl;
        return false;
    }

    // the same circuit with and without multithreading, big enough for Eigen to parallelize the SVDs (bond dimensions
    // up to 64)
    std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
    FillOneQubitGates(gates);
    FillTwoQubitGates(gates);

    constexpr int nrQubits = 12;
    const auto circuit = GenerateRandomCircuitWithGatesNoAdjacent(gates, 300, 300, nrQubits);

#ifdef _OPENMP
    const int threadsBefore = omp_get_max_threads();
#endif

    QC::TensorNetworks::MPSSimulator mpsMultithreaded(nrQubits);
    QC::TensorNetworks::MPSSimulator mpsSingleThreaded(nrQubits);
    mpsSingleThreaded.SetMultithreading(false);

    for (const auto &gate : circuit)
    {
        mpsMultithreaded.ApplyGate(*gate);
        mpsSingleThreaded.ApplyGate(*gate);
    }

#ifdef _OPENMP
    if (omp_get_max_threads() != threadsBefore)
    {
        std::cout << "The OpenMP number of threads was not restored after single threaded operations" << std::endl;
        return false;
    }
#endif

    const auto state1 = mpsMultithreaded.getRegisterStorage();
    const auto state2 = mpsSingleThreaded.getRegisterStorage();
    for (Eigen::Index s = 0; s < state1.size(); ++s)
        if (!approxEqual(state1[s], state2[s], 1E-12))
        {
            std::cout << "Single threaded and multithreaded simulations differ" << std::endl;
            return false;
        }

    std::cout << "Success" << std::endl;
    return true;
}

bool MeetingPositionFallbackTestMPS()
{
    for (const bool optimal : {false, true})
        for (const int position : {-1, 99})
        {
            QC::TensorNetworks::MPSSimulator mps(4);
            QC::QubitRegister<> reference(4);
            mps.SetUseOptimalMeetingPosition(optimal);
            int callbacks = 0;
            mps.SetMeetingPositionCallback([&](const auto &) {
                ++callbacks;
                return position;
            });
            QC::Gates::HadamardGate<> h;
            QC::Gates::CNOTGate<> cx;
            mps.ApplyGate(h, 0);
            reference.ApplyGate(h, 0);
            mps.ApplyGate(cx, 3, 0);
            reference.ApplyGate(cx, 3, 0);
            if (callbacks != 1 || (mps.getRegisterStorage() - reference.getRegisterStorage()).norm() > 1E-10)
                return false;
        }
    return true;
}

// A SWAP gate is applied by exchanging the two logical qubits in the map: no routing, no tensor
// update and no bond notification, while the simulated state still matches the statevector.
static bool SwapRelabelTestMPS()
{
    std::cout << "\nMPS simulator SWAP relabelling test" << std::endl;

    constexpr int nrQubits = 6;
    QC::TensorNetworks::MPSSimulator mps(nrQubits);
    QC::QubitRegister<> reference(nrQubits);

    QC::Gates::HadamardGate<> h;
    QC::Gates::CNOTGate<> cx;
    QC::Gates::RyGate<> ry(M_PI / 5);
    QC::Gates::SwapGate<> swap;

    const auto applyBoth = [&](const auto &gate, int qubit, int controllingQubit = 0) {
        mps.ApplyGate(gate, qubit, controllingQubit);
        reference.ApplyGate(gate, qubit, controllingQubit);
    };

    applyBoth(h, 0);
    applyBoth(cx, 1, 0);
    applyBoth(ry, 3);
    applyBoth(cx, 4, 3);
    applyBoth(h, 5);

    int meetingCallbacks = 0, bondCallbacks = 0;
    mps.SetMeetingPositionCallback([&](const auto &) {
        ++meetingCallbacks;
        return -1;
    });
    mps.SetBondDimensionCallback([&](const auto &) { ++bondCallbacks; });

    const auto bondsBefore = mps.getBondDimensions();
    auto expectedMap = mps.getQubitsMap();
    std::swap(expectedMap[0], expectedMap[5]);
    std::swap(expectedMap[1], expectedMap[3]);

    applyBoth(swap, 0, 5); // chain ends
    applyBoth(swap, 3, 1); // non-adjacent, both entangled

    if (meetingCallbacks != 0 || bondCallbacks != 0 || mps.getBondDimensions() != bondsBefore ||
        mps.getQubitsMap() != expectedMap)
    {
        std::cout << "A SWAP gate routed qubits or changed the tensors instead of relabelling them" << std::endl;
        return false;
    }

    if ((mps.getRegisterStorage() - reference.getRegisterStorage()).norm() > 1E-10)
    {
        std::cout << "State mismatch after relabelling SWAP gates" << std::endl;
        return false;
    }

    // later gates, including routed ones, act on the relabelled qubits
    applyBoth(cx, 5, 1);
    applyBoth(ry, 0);
    applyBoth(cx, 2, 4);
    applyBoth(swap, 2, 3);
    applyBoth(cx, 0, 3);

    if ((mps.getRegisterStorage() - reference.getRegisterStorage()).norm() > 1E-10)
    {
        std::cout << "State mismatch after gates following relabelling SWAP gates" << std::endl;
        return false;
    }

    std::cout << "Success" << std::endl;

    return true;
}

// Long chains and invalid input. Sampling keeps its prefix row normalized: the joint probability
// of the prefix used to underflow past ~1074 unbiased qubits, after which every outcome was 1.
// Out-of-range qubits throw instead of indexing the qubit map, and non-finite gates are rejected
// (the SVD of a non-finite block fails and must not be read).
static bool LongChainSamplingAndValidationTestMPS()
{
    std::cout << "\nMPS simulator long-chain sampling and input validation" << std::endl;

    using IndexType = QC::TensorNetworks::MPSSimulatorInterface::IndexType;
    constexpr int nrQubits = 1300;
    QC::TensorNetworks::MPSSimulator mps(nrQubits, 17);
    QC::Gates::HadamardGate<> h;
    for (int q = 0; q < nrQubits; ++q)
        mps.ApplyGate(h, q);
    int ones = 0, sampled = 0;
    for (int shot = 0; shot < 20; ++shot)
    {
        const auto outcome = mps.MeasureNoCollapse();
        for (int q = 1100; q < nrQubits; ++q)
        {
            ones += outcome.at(q) ? 1 : 0;
            ++sampled;
        }
    }
    const double fraction = static_cast<double>(ones) / sampled;
    if (std::abs(fraction - 0.5) > 0.05)
    {
        std::cout << "Sampling qubits 1100+ of a 1300-qubit |+> state gave a one-fraction of " << fraction << std::endl;
        return false;
    }

    QC::TensorNetworks::MPSSimulator smalls(3, 5);
    const auto expectInvalid = [](auto &&callable, const char *description) {
        try
        {
            callable();
        }
        catch (const std::invalid_argument &)
        {
            return true;
        }
        catch (const std::exception &ex)
        {
            std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
            return false;
        }
        std::cout << description << " did not throw std::invalid_argument" << std::endl;
        return false;
    };
    Eigen::MatrixXcd z = Eigen::MatrixXcd::Identity(2, 2);
    z(1, 1) = -1.;
    Eigen::MatrixXcd nonFinite = Eigen::MatrixXcd::Identity(4, 4);
    nonFinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
    const QC::Gates::TwoQubitsGate<> nonFiniteGate(nonFinite);
    if (!expectInvalid(
            [&] {
                smalls.MeasureQubits({0, 99});
            },
            "MeasureQubits with an invalid qubit") ||
        !expectInvalid([&] { smalls.MeasureNoCollapse(std::set<IndexType>{99}); },
                       "MeasureNoCollapse with an invalid qubit") ||
        !expectInvalid([&] { smalls.MoveAtBeginningOfChain({99}); }, "MoveAtBeginningOfChain with an invalid qubit") ||
        !expectInvalid([&] { smalls.ExpectationValue({QC::Gates::AppliedGate<>(z, 99)}); },
                       "ExpectationValue with an invalid qubit") ||
        !expectInvalid([&] { smalls.ApplyGate(nonFiniteGate, 1, 0); }, "A non-finite two-qubit gate"))
        return false;
    if (std::abs(smalls.getBasisStateAmplitude(0) - std::complex<double>(1., 0.)) > 1E-12)
    {
        std::cout << "Rejected MPS operations changed the state" << std::endl;
        return false;
    }

    std::cout << "Success" << std::endl;
    return true;
}

static void ObservableCacheBudgetTestMPS()
{
    using namespace QC::TensorNetworkRegression;
    using Simulator = QC::TensorNetworks::MPSSimulatorImpl;
    // Five nested Bell pairs carry a bond of 32 through a run of |0> sites.
    // Construct their exact right-canonical tensors directly, avoiding large
    // SVDs and exponentially large dense references in this cache regression.
    constexpr Eigen::Index pairs = 5, bond = 1 << pairs;
    const size_t environmentBytes = sizeof(Matrix) + bond * bond * sizeof(Complex);
    const auto capacity = QC::TensorNetworks::ObservableCacheBytes / environmentBytes;
    for (const bool exceedBudget : {false, true})
    {
        const Eigen::Index middle = static_cast<Eigen::Index>(capacity) + (exceedBudget ? 8 : -8);
        const Eigen::Index count = 2 * pairs + middle;
        auto state = std::make_shared<QC::TensorNetworks::MPSSimulatorBaseState>();
        size_t prefixBytes = 0;
        Eigen::Index left = 1;
        for (Eigen::Index q = 0; q < count; ++q)
        {
            const Eigen::Index right = q < pairs ? 2 * left : q < pairs + middle ? left : left / 2;
            Simulator::GammaType gamma(left, 2, right);
            gamma.setZero();
            if (q < pairs)
                for (Eigen::Index l = 0; l < left; ++l)
                    for (Eigen::Index p = 0; p < 2; ++p)
                        gamma(l, p, 2 * l + p) = 1. / std::sqrt(2.);
            else if (q < pairs + middle)
                for (Eigen::Index l = 0; l < left; ++l)
                    gamma(l, 0, l) = 1.;
            else
                for (Eigen::Index r = 0; r < right; ++r)
                    for (Eigen::Index p = 0; p < 2; ++p)
                        gamma(2 * r + p, p, r) = 1.;
            state->gammas.push_back(std::move(gamma));
            if (q + 1 < count)
                state->lambdas.push_back(Simulator::LambdaType::Constant(right, 1. / std::sqrt(double(right))));
            if (q + 2 < count)
                prefixBytes += sizeof(Matrix) + right * right * sizeof(Complex);
            left = right;
        }
        Require((prefixBytes > QC::TensorNetworks::ObservableCacheBytes) == exceedBudget,
                "MPS fixture did not straddle the observable cache budget");
        Simulator simulator(count, 1);
        simulator.SetMultithreading(false);
        simulator.setState(state);
        const std::string parity(count, 'Z');
        std::string crossed = parity, unmatched = parity;
        crossed.back() = 'X';
        unmatched[count - 2] = 'I';
        const auto values = simulator.ExpectationValues({parity, crossed, unmatched, parity});
        Require(values.size() == 4, "Cache-boundary batch changed the number of results");
        // Each Bell pair has <ZZ>=1, <ZX>=0, <ZI>=0; the middle sites have <Z>=1.
        Close(values[0], 1., "MPS cache boundary changed Bell-pair parity");
        Close(values[1], 0., "MPS cache boundary changed crossed Bell-pair correlation");
        Close(values[2], 0., "MPS cache boundary changed an unmatched Bell-pair observable");
        Close(values[3], 1., "MPS cache boundary changed a duplicate result");
    }
}

static bool ObservableAndSummaryRegressionTestMPS()
{
    using namespace QC::TensorNetworkRegression;
    using Simulator = QC::TensorNetworks::MPSSimulator;
    std::cout << "\nMPS batch observables, cache budget and bond summaries" << std::endl;
    try
    {
        CheckObservableBatch<Simulator>([](Simulator &simulator, const std::string &pauli) {
            return simulator.ExpectationValue(PauliGates(pauli));
        });
        CheckBondSummary<Simulator>();
        ObservableCacheBudgetTestMPS();
        Simulator simulator(6, 1);
        Prepare(simulator, 4);
        QC::QubitRegister<> dense(6);
        Eigen::VectorXcd storage = simulator.getRegisterStorage();
        dense.setRegisterStorageFastNoNormalize(storage);
        const std::vector<Gate> operators = {Gate(Pauli('X'), 2), Gate(Pauli('Y'), 2), Gate(Unitary(2, 711), 4),
                                             Gate(Pauli('Z'), 0)};
        Close(simulator.ExpectationValue(operators), dense.ExpectationValue(operators),
              "General MPS observable changed repeated operator order");
        std::cout << "Success" << std::endl;
        return true;
    }
    catch (const std::exception &error)
    {
        std::cout << "MPS observable/summary regression: " << error.what() << std::endl;
        return false;
    }
}

// The Kraus steps and the gauge mode they start: probabilities on a chain whose gates have truncated, reads that
// must leave the state as it is, small Schmidt values, branches of near-zero weight, and initializations that are
// refused. Every value is compared with the dense state, the product of the tensors, which uses nothing cached.
static bool KrausGaugeModeRegressionTestMPS()
{
    using namespace QC::TensorNetworkRegression;
    using Impl = QC::TensorNetworks::MPSSimulatorImpl;
    std::cout << "\nMPS Kraus steps: truncated chains, reads, small values, refused initializations" << std::endl;
    try
    {
        const auto probabilityOfOne = [](const Eigen::VectorXcd &psi, int qubit) {
            double p = 0.;
            for (Eigen::Index i = 0; i < psi.size(); ++i)
                if ((i >> qubit) & 1)
                    p += std::norm(psi[i]);
            return p / psi.squaredNorm();
        };
        Matrix hadamard(2, 2), cnot = Matrix::Zero(4, 4), zero = Matrix::Zero(2, 2), one = Matrix::Zero(2, 2);
        hadamard << 1., 1., 1., -1.;
        hadamard /= std::sqrt(2.);
        cnot(0, 0) = cnot(1, 1) = cnot(2, 3) = cnot(3, 2) = 1.; // Gate(cnot, target, control)
        zero(0, 0) = 1.;
        one(1, 1) = 1.;
        const Matrix identity2 = Matrix::Identity(2, 2), identity4 = Matrix::Identity(4, 4);
        const auto turn = [](double angle) {
            Matrix rotation(2, 2);
            rotation << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
            return rotation;
        };

        // After gates that truncate, the canonical form no longer gives the environments of a site.
        {
            const int qubits = 10;
            Impl truncated(qubits);
            truncated.setLimitBondDimension(8);
            unsigned int seed = 100;
            for (int layer = 0; layer < 12; ++layer)
            {
                for (int q = 0; q < qubits; ++q)
                    truncated.ApplyGate(Gate(Unitary(2, ++seed), q));
                for (int q = layer % 2; q + 1 < qubits; q += 2)
                    truncated.ApplyGate(Gate(Unitary(4, ++seed), q, q + 1));
            }
            Require(truncated.getMaxBondDimension() == 8, "The MPS test chain did not reach its bond limit");
            const Eigen::VectorXcd psi = truncated.getRegisterStorage();
            for (int q = 0; q < qubits; ++q)
            {
                Impl read = truncated;
                Close(read.LocalDensityMatrix(q)(1, 1), probabilityOfOne(psi, q),
                      "MPS local density matrix is stale after truncating gates");
                Impl projected = truncated;
                projected.ApplyKrausOperator({one}, 0.5, q);
                Close(projected.getRegisterStorage().squaredNorm(), 1.,
                      "MPS Kraus step on a truncated chain left the norm off");
            }
            // The draw follows the state: just below and just above the probability of the first operator.
            const double p0 = 1. - probabilityOfOne(psi, 3);
            Impl below = truncated, above = truncated;
            Require(below.ApplyKrausOperator({zero, one}, p0 - 1E-6, 3) == 0 &&
                        above.ApplyKrausOperator({zero, one}, p0 + 1E-6, 3) == 1,
                    "MPS Kraus draw on a truncated chain does not follow the state");
            // Channels between truncating gates keep it exact.
            Impl noisy(qubits);
            noisy.setLimitBondDimension(8);
            Matrix noJump = Matrix::Zero(2, 2), jump = Matrix::Zero(2, 2);
            noJump(0, 0) = 1.;
            noJump(1, 1) = std::sqrt(0.95);
            jump(0, 1) = std::sqrt(0.05);
            std::mt19937_64 random(5);
            std::uniform_real_distribution<double> uniform;
            for (int layer = 0; layer < 12; ++layer)
            {
                for (int q = 0; q < qubits; ++q)
                    noisy.ApplyGate(Gate(Unitary(2, ++seed), q));
                for (int q = layer % 2; q + 1 < qubits; q += 2)
                    noisy.ApplyGate(Gate(Unitary(4, ++seed), q, q + 1));
                if (layer >= 6) // the first six layers truncate on the canonical path
                    for (int q = 0; q < qubits; ++q)
                        noisy.ApplyKrausOperator({noJump, jump}, uniform(random), q);
            }
            const Eigen::VectorXcd mixed = noisy.getRegisterStorage();
            Close(mixed.squaredNorm(), 1., "MPS trajectory through truncating gates left the norm off");
            for (int q = 0; q < qubits; ++q)
                Close(noisy.LocalDensityMatrix(q)(1, 1), probabilityOfOne(mixed, q),
                      "MPS local density matrix in gauge mode");
        }

        // Restoring the canonical form for a read applies none of the cuts the user asked for.
        {
            const double gamma = 0.999, expected = (1. - gamma) / (2. - gamma);
            Matrix noJump = Matrix::Zero(2, 2), jump = Matrix::Zero(2, 2);
            noJump(0, 0) = 1.;
            noJump(1, 1) = std::sqrt(1. - gamma);
            jump(0, 1) = std::sqrt(gamma);
            const auto prepare = [&](Impl &mps) {
                mps.setLimitEntanglement(0.1);
                mps.setTruncationMode(QC::TensorNetworks::MPSSimulatorInterface::TruncationMode::RelativeToMax);
                mps.ApplyGate(Gate(hadamard, 0));
                mps.ApplyGate(Gate(cnot, 1, 0));
                Require(mps.ApplyKrausOperator({noJump, jump}, 0., 1) == 0 && mps.InGaugeMode(),
                        "MPS damping step did not take its first operator");
            };
            const auto p11 = [](const Impl &mps) {
                const Eigen::VectorXcd psi = mps.getRegisterStorage();
                return std::norm(psi[3]) / psi.squaredNorm();
            };
            Impl saved(2), probed(2), observed(2), sampled(2), compressed(2);
            prepare(saved);
            const auto state = saved.getState();
            Close(p11(saved), expected, "Saving an MPS state compressed it");
            Require(saved.getMaxBondDimension() == 2 && !saved.InGaugeMode(),
                    "Saving an MPS state did not restore the canonical form as it was");
            prepare(probed);
            Close(probed.GetProbability(1, false), expected, "MPS probability after a Kraus step");
            Close(p11(probed), expected, "Reading an MPS probability compressed the state");
            prepare(observed);
            Close(observed.ExpectationValues({"ZI"})[0], 1. - 2. * expected,
                  "MPS expectation value after a Kraus step");
            Close(p11(observed), expected, "Reading an MPS expectation value compressed the state");
            prepare(sampled);
            sampled.MeasureNoCollapse();
            Close(p11(sampled), expected, "Sampling an MPS without collapse compressed the state");
            // Asked for by name, the compression is applied.
            prepare(compressed);
            compressed.ReCanonicalize();
            Close(p11(compressed), 0., "ReCanonicalize did not apply the relative cutoff");
        }

        // A small Schmidt value survives a gate in gauge mode: at the ends of a chain, and between environments that
        // were accumulated but are well conditioned.
        {
            const double angle = 1E-8;
            Impl pair(2);
            pair.ApplyGate(Gate(turn(angle), 0));
            pair.ApplyGate(Gate(cnot, 1, 0));
            pair.ApplyKrausOperator({identity2}, 0.5, 0);
            pair.ApplyGate(Gate(identity4, 1, 0));
            Require(pair.getMaxBondDimension() == 2 &&
                        std::abs(std::abs(pair.getRegisterStorage()[3]) / std::sin(angle) - 1.) < 1E-6,
                    "A gate in MPS gauge mode dropped a small amplitude");

            Impl chain(4);
            chain.ApplyGate(Gate(hadamard, 0));
            chain.ApplyGate(Gate(cnot, 1, 0));
            chain.ApplyGate(Gate(hadamard, 3));
            chain.ApplyGate(Gate(cnot, 2, 3));
            Matrix weak = Matrix::Identity(4, 4); // qubit 2 turned by a small angle when qubit 1 is set
            weak.block(2, 2, 2, 2) = turn(angle);
            chain.ApplyGate(Gate(weak, 2, 1));
            const Eigen::VectorXcd before = chain.getRegisterStorage();
            chain.ApplyKrausOperator({identity2}, 0.5, 0);
            chain.ApplyKrausOperator({identity2}, 0.5, 3);
            chain.ApplyGate(Gate(identity4, 2, 1));
            Require(chain.InGaugeMode() && (chain.getRegisterStorage() - before).norm() < 1E-13,
                    "A gate between accumulated MPS environments dropped a small Schmidt value");
        }

        // A branch whose weight is what a cancellation in the density matrix left still leaves a state of unit norm.
        {
            Matrix minus(2, 2);
            minus << 0.5, -0.5, -0.5, 0.5; // |-><-|
            for (const double delta : {1E-3, 1E-5, 1E-6})
                for (const int site : {0, 1, 2})
                {
                    Impl mps(3);
                    for (int q = 0; q < 3; ++q)
                        mps.ApplyGate(Gate(Unitary(2, 40 + q), q));
                    mps.ApplyGate(Gate(Unitary(4, 50), 0, 1));
                    mps.ApplyGate(Gate(Unitary(4, 51), 1, 2));
                    // Put the site in |+> turned by delta, whatever the others hold: <-|psi> = -sin(delta).
                    mps.ApplyKrausOperator({zero}, 0.5, site);
                    mps.ApplyGate(Gate(turn(std::acos(-1.) / 4. + delta), site));
                    mps.ApplyKrausOperator({minus}, 0.5, site);
                    Close(mps.getRegisterStorage().squaredNorm(), 1.,
                          "MPS Kraus step on a branch of near-zero weight left the norm off");
                }
        }

        // An initialization that is refused leaves the chain, and what describes it, as they were.
        {
            // A damping step that keeps the chain in gauge mode: the state is (|00> + sqrt(1/2) |11>) / sqrt(3/2),
            // which the lambdas of the Bell pair no longer describe.
            Matrix halved = Matrix::Zero(2, 2);
            halved(0, 0) = 1.;
            halved(1, 1) = std::sqrt(0.5);
            const auto prepare = [&](Impl &mps) {
                mps.ApplyGate(Gate(hadamard, 0));
                mps.ApplyGate(Gate(cnot, 1, 0));
                mps.ApplyKrausOperator({halved}, 0.5, 0);
                Require(mps.InGaugeMode(), "The MPS test chain is not in gauge mode");
            };
            Impl refused(2), ignored(2), accepted(2);
            prepare(refused);
            Throws<std::invalid_argument>([&] { refused.setToBasisState(std::vector<bool>(3, false)); },
                                          "An MPS basis state with too many bits was accepted");
            Close(refused.GetProbability(1, true), 2. / 3., "A refused MPS initialization corrupted the probabilities");
            prepare(ignored);
            ignored.setToBasisState(size_t{4});
            Close(ignored.GetProbability(1, true), 2. / 3.,
                  "An ignored MPS initialization corrupted the probabilities");
            prepare(accepted);
            Require(accepted.InGaugeMode(), "The MPS test chain is not in gauge mode");
            accepted.setToBasisState(size_t{2});
            Require(!accepted.InGaugeMode(), "An accepted MPS initialization kept the gauge mode");
            Close(accepted.GetProbability(1, false), 1., "MPS basis state after a Kraus step");
            Close(accepted.GetProbability(0, true), 1., "MPS basis state after a Kraus step");
        }

        // Without truncation nothing is cut: a weak GHZ state keeps its 1e-8 amplitude through an identity Kraus
        // operator and an identity gate in every order, and the densities agree with the tensors afterwards.
        {
            const double angle = 1E-8;
            for (const int krausQubit : {0, 1, 2})
                for (const int gateLeft : {0, 1})
                {
                    Impl ghz(3);
                    ghz.ApplyGate(Gate(turn(angle), 0));
                    ghz.ApplyGate(Gate(cnot, 1, 0));
                    ghz.ApplyGate(Gate(cnot, 2, 1));
                    ghz.ApplyKrausOperator({identity2}, 0.5, krausQubit);
                    ghz.ApplyGate(Gate(identity4, gateLeft + 1, gateLeft));
                    const Eigen::VectorXcd psi = ghz.getRegisterStorage();
                    Require(std::abs(std::abs(psi[7]) / std::sin(angle) - 1.) < 1E-6 && ghz.getBondDimension(0) == 2 &&
                                ghz.getBondDimension(1) == 2,
                            "A gate in MPS gauge mode dropped the small amplitude of a GHZ state");
                    for (int q = 0; q < 3; ++q)
                    {
                        const double p = probabilityOfOne(psi, q);
                        Require(std::abs(ghz.LocalDensityMatrix(q)(1, 1).real() - p) < 1E-6 * p,
                                "An MPS environment is stale after a gate in gauge mode");
                    }
                }
        }

        // A branch of near-zero weight where the environment of the site is not diagonal: qubit 0 of a Bell pair is
        // projected on a generic state, which leaves qubit 1 in it, and qubit 1 meets the projector on the state at
        // delta from the orthogonal one. Its weight, sin(delta)^2, is found down to the rounding of the tensors, and
        // the state leaves with unit norm; the projector on the orthogonal state itself has no weight, is refused,
        // and leaves the state as it was.
        {
            const auto projector = [](double c, double s) {
                Matrix p(2, 2);
                p << c * c, c * s, c * s, s * s;
                return p;
            };
            for (const double alpha : {0.3, 1.1})
            {
                const auto prepare = [&](Impl &mps) {
                    mps.ApplyGate(Gate(hadamard, 0));
                    mps.ApplyGate(Gate(cnot, 1, 0));
                    mps.ApplyKrausOperator({projector(std::cos(alpha), std::sin(alpha))}, 0.5, 0);
                };
                for (const double delta : {1E-3, 1E-6, 1E-8, 1E-9, 1E-11})
                {
                    Impl mps(2);
                    prepare(mps);
                    mps.ApplyKrausOperator({projector(-std::sin(alpha + delta), std::cos(alpha + delta))}, 0.5, 1);
                    const Eigen::VectorXcd psi = mps.getRegisterStorage();
                    Close(psi.squaredNorm(), 1.,
                          "MPS Kraus step on a near-zero branch between projected sites left the norm off");
                    // qubit 0 stays in the state it was projected on
                    Require(std::abs(probabilityOfOne(psi, 0) - std::sin(alpha) * std::sin(alpha)) < 1E-9,
                            "MPS Kraus step on a near-zero branch changed the other qubit");
                }
                Impl orthogonal(2);
                prepare(orthogonal);
                const Eigen::VectorXcd before = orthogonal.getRegisterStorage();
                Throws<std::domain_error>(
                    [&] { orthogonal.ApplyKrausOperator({projector(-std::sin(alpha), std::cos(alpha))}, 0.5, 1); },
                    "An MPS Kraus operator with no weight was applied");
                Require((orthogonal.getRegisterStorage() - before).norm() < 1E-13,
                        "A refused MPS Kraus operator changed the state");
            }
        }

        // Reads on a chain whose gates have truncated, with no Kraus operator anywhere: a probability, an
        // expectation value and a sample are taken from the environments of the sites, computed from the tensors,
        // and the chain goes on afterwards as if it had not been read.
        {
            const int qubits = 8;
            Impl chain(qubits);
            chain.setLimitBondDimension(4);
            unsigned int seed = 900;
            for (int layer = 0; layer < 10; ++layer)
            {
                for (int q = 0; q < qubits; ++q)
                    chain.ApplyGate(Gate(Unitary(2, ++seed), q));
                for (int q = layer % 2; q + 1 < qubits; q += 2)
                    chain.ApplyGate(Gate(Unitary(4, ++seed), q, q + 1));
            }
            Require(chain.getMaxBondDimension() == 4, "The MPS test chain did not reach its bond limit");
            const Eigen::VectorXcd psi = chain.getRegisterStorage();
            for (int q = 0; q < qubits; ++q)
            {
                const double p = probabilityOfOne(psi, q);
                Impl probed = chain;
                Close(probed.GetProbability(q, false), p, "MPS probability read from a truncated chain");
                Impl observed = chain;
                std::string pauli(qubits, 'I');
                pauli[q] = 'Z';
                Close(observed.ExpectationValues({pauli})[0], 1. - 2. * p,
                      "MPS expectation value read from a truncated chain");
            }
            Impl sampled = chain;
            sampled.SetSeed(17);
            const int shots = 20000;
            std::vector<int> ones(qubits, 0);
            for (int shot = 0; shot < shots; ++shot)
                for (const auto &outcome : sampled.MeasureNoCollapse())
                    ones[outcome.first] += outcome.second ? 1 : 0;
            for (int q = 0; q < qubits; ++q)
            {
                const double p = probabilityOfOne(psi, q);
                Require(std::abs(double(ones[q]) / shots - p) < 5. * std::sqrt(p * (1. - p) / shots) + 1E-3,
                        "Samples of a truncated MPS do not follow its state");
            }
            // Measurements with collapse: the outcomes follow the probability of the state.
            {
                const int trials = 2000;
                int measuredOnes = 0;
                for (int trial = 0; trial < trials; ++trial)
                {
                    Impl measured = chain;
                    measured.SetSeed(3000 + trial);
                    measuredOnes += measured.MeasureQubit(5) ? 1 : 0;
                }
                const double p = probabilityOfOne(psi, 5);
                Require(std::abs(double(measuredOnes) / trials - p) < 5. * std::sqrt(p * (1. - p) / trials) + 1E-3,
                        "Measurements of a truncated MPS do not follow its state");
            }
            // Reads leave the state as it was (up to its norm, which the first read puts right), and the gates
            // that follow give what they give on a chain that was never read.
            Impl read = chain, unread = chain, densityFirst = chain;
            const Eigen::VectorXcd normalized = psi / psi.norm();
            read.GetProbability(3, false);
            read.ExpectationValues({std::string(qubits, 'Z')});
            read.MeasureNoCollapse();
            read.LocalDensityMatrix(4);
            // On its own copy, and nothing else read before the gates: the first read of a chain is the one that
            // computes the environments, and each kind of read has to leave them as kept for reads only.
            densityFirst.LocalDensityMatrix(4);
            const auto apart = [](const Eigen::VectorXcd &a, const Eigen::VectorXcd &b) {
                return (a / a.norm() - b / b.norm()).cwiseAbs().maxCoeff();
            };
            Require(apart(read.getRegisterStorage(), normalized) < 1E-12 &&
                        std::abs(read.getRegisterStorage().squaredNorm() - 1.) < 1E-12,
                    "Reading a truncated MPS changed its state");
            for (int layer = 0; layer < 2; ++layer)
                for (int q = layer % 2; q + 1 < qubits; q += 2)
                {
                    const Matrix u = Unitary(4, ++seed);
                    read.ApplyGate(Gate(u, q, q + 1));
                    unread.ApplyGate(Gate(u, q, q + 1));
                    densityFirst.ApplyGate(Gate(u, q, q + 1));
                }
            Require(apart(read.getRegisterStorage(), unread.getRegisterStorage()) < 1E-12,
                    "Gates after a read of a truncated MPS differ from those on a chain that was not read");
            Require(apart(densityFirst.getRegisterStorage(), unread.getRegisterStorage()) < 1E-12,
                    "Gates after the local density matrix of a truncated MPS differ from those on a chain that was not "
                    "read");
        }

        // A measurement on a chain that holds Kraus operators is the channel of the two projectors, and splits the
        // site again where the outcome left it: nothing sweeps the chain, the bonds of the measured qubit lose the
        // rank the collapse took, and the state is the projected one.
        {
            Impl ghz(4);
            ghz.ApplyGate(Gate(hadamard, 0));
            for (int q = 0; q < 3; ++q)
                ghz.ApplyGate(Gate(cnot, q + 1, q));
            Matrix noJump = Matrix::Zero(2, 2), jump = Matrix::Zero(2, 2);
            noJump(0, 0) = 1.;
            noJump(1, 1) = std::sqrt(0.9);
            jump(0, 1) = std::sqrt(0.1);
            Require(ghz.ApplyKrausOperator({noJump, jump}, 0., 2) == 0 && ghz.InGaugeMode(),
                    "MPS damping step did not take its first operator");
            Eigen::VectorXcd expected = ghz.getRegisterStorage();
            ghz.SetSeed(5);
            const bool outcome = ghz.MeasureQubit(1);
            for (Eigen::Index i = 0; i < expected.size(); ++i)
                if (((i >> 1) & 1) != (outcome ? 1 : 0))
                    expected[i] = 0.;
            expected /= expected.norm();
            const Eigen::VectorXcd collapsed = ghz.getRegisterStorage();
            Require(ghz.InGaugeMode() && ghz.getBondDimension(0) == 1 && ghz.getBondDimension(1) == 1,
                    "A measurement in MPS gauge mode did not split the measured site again");
            Close(collapsed.squaredNorm(), 1., "A measurement in MPS gauge mode left the norm off");
            Require(std::abs(std::abs(expected.dot(collapsed)) - 1.) < 1E-12,
                    "A measurement in MPS gauge mode did not project the state");
            // The other qubits follow: a GHZ state collapses as a whole.
            for (int q = 0; q < 4; ++q)
                Close(ghz.GetProbability(q, false), outcome ? 1. : 0.,
                      "MPS probabilities after a measurement in gauge mode");
            // A jump is a projection too.
            Impl jumped(3);
            for (int q = 0; q < 3; ++q)
                jumped.ApplyGate(Gate(Unitary(2, 70 + q), q));
            jumped.ApplyGate(Gate(Unitary(4, 80), 0, 1));
            jumped.ApplyGate(Gate(Unitary(4, 81), 1, 2));
            Eigen::VectorXcd reference = jumped.getRegisterStorage();
            Require(jumped.ApplyKrausOperator({noJump, jump}, std::nextafter(1., 0.), 1) == 1,
                    "MPS damping step did not take its jump");
            for (Eigen::Index i = 0; i < reference.size(); ++i)
                reference[i] = (i >> 1) & 1 ? Complex(0.) : reference[i | 2]; // |0><1| on qubit 1
            reference /= reference.norm();
            const Eigen::VectorXcd after = jumped.getRegisterStorage();
            // The qubit that jumped is in |0>; its neighbours may still be entangled through it.
            Require(jumped.InGaugeMode(), "A jump in MPS gauge mode restored the canonical form");
            Require(std::abs(after.squaredNorm() - 1.) < 1E-12 && std::abs(std::abs(reference.dot(after)) - 1.) < 1E-12,
                    "A jump in MPS gauge mode did not leave the projected state");
            Close(jumped.GetProbability(1, true), 1., "MPS probability of the qubit that jumped");
            for (const int q : {0, 2})
                Close(jumped.GetProbability(q, false), probabilityOfOne(reference, q),
                      "MPS probabilities after a jump in gauge mode");
        }

        // The same steps against a dense reference, operation by operation: random circuits of weak entanglers
        // (amplitudes down to 1e-13), unitaries, identity and damping channels, measurements in a generic basis and
        // projectors on states nearly orthogonal to the one just projected on. The reference applies the operator
        // the chain reports having drawn and normalizes; the two states must agree after every step, to the
        // rounding of the state, or to what a near-null branch leaves of the digits of the amplitudes it subtracts.
        {
            const auto apply1 = [](Eigen::VectorXcd &psi, const Matrix &e, int q) {
                const Eigen::Index bit = Eigen::Index{1} << q;
                for (Eigen::Index i = 0; i < psi.size(); ++i)
                    if (!(i & bit))
                    {
                        const Complex a0 = psi[i], a1 = psi[i | bit];
                        psi[i] = e(0, 0) * a0 + e(0, 1) * a1;
                        psi[i | bit] = e(1, 0) * a0 + e(1, 1) * a1;
                    }
            };
            // `low` is the low bit of the matrix index
            const auto apply2 = [](Eigen::VectorXcd &psi, const Matrix &e, int low, int high) {
                const Eigen::Index bl = Eigen::Index{1} << low, bh = Eigen::Index{1} << high;
                for (Eigen::Index i = 0; i < psi.size(); ++i)
                    if (!(i & bl) && !(i & bh))
                    {
                        const Eigen::Index index[4] = {i, i | bl, i | bh, i | bl | bh};
                        Complex a[4], out[4];
                        for (int c = 0; c < 4; ++c)
                            a[c] = psi[index[c]];
                        for (int r = 0; r < 4; ++r)
                        {
                            out[r] = 0.;
                            for (int c = 0; c < 4; ++c)
                                out[r] += e(r, c) * a[c];
                        }
                        for (int c = 0; c < 4; ++c)
                            psi[index[c]] = out[c];
                    }
            };
            const auto distance = [](const Eigen::VectorXcd &a, const Eigen::VectorXcd &b) {
                const Complex overlap = b.dot(a);
                const Complex phase = std::abs(overlap) > 0. ? overlap / std::abs(overlap) : Complex(1.);
                return (a - phase * b).cwiseAbs().maxCoeff();
            };
            unsigned int unitarySeed = 5000;
            for (unsigned int trial = 0; trial < 40; ++trial)
            {
                std::mt19937_64 random(100 + trial);
                std::uniform_real_distribution<double> uniform;
                const int qubits = 3 + static_cast<int>(trial % 4);
                Impl mps(qubits);
                Eigen::VectorXcd reference = Eigen::VectorXcd::Zero(Eigen::Index{1} << qubits);
                reference[0] = 1.;
                double tolerance = 1E-13;
                const auto kraus = [&](const std::vector<Matrix> &operators, int q) {
                    const size_t drawn = mps.ApplyKrausOperator(operators, uniform(random), q);
                    apply1(reference, operators[drawn], q);
                    reference /= reference.norm();
                };
                const double angles[] = {1E-3, 1E-6, 1E-8, 1E-11, 1E-13};
                for (int step = 0; step < 40; ++step)
                {
                    const int q = static_cast<int>(random() % qubits), pair = static_cast<int>(random() % (qubits - 1));
                    switch (random() % 9)
                    {
                    case 7: // a measurement: the reference follows the outcome
                    {
                        const bool outcome = mps.MeasureQubit(q);
                        apply1(reference, outcome ? one : zero, q);
                        reference /= reference.norm();
                        break;
                    }
                    case 0: // a weak entangler
                    {
                        const Matrix weak = turn(angles[random() % 5]);
                        mps.ApplyGate(Gate(weak, pair));
                        apply1(reference, weak, pair);
                        mps.ApplyGate(Gate(cnot, pair + 1, pair));
                        apply2(reference, cnot, pair + 1, pair);
                        break;
                    }
                    case 1: {
                        const Matrix u = Unitary(2, ++unitarySeed);
                        mps.ApplyGate(Gate(u, q));
                        apply1(reference, u, q);
                        break;
                    }
                    case 2: {
                        const Matrix u = random() % 2 ? Unitary(4, ++unitarySeed) : identity4;
                        mps.ApplyGate(Gate(u, pair, pair + 1));
                        apply2(reference, u, pair, pair + 1);
                        break;
                    }
                    case 3:
                        kraus({identity2}, q);
                        break;
                    case 4:
                    case 5: {
                        const double gammas[] = {1E-6, 1E-3, 0.05, 0.3};
                        const double gamma = gammas[random() % 4];
                        Matrix noJump = Matrix::Zero(2, 2), jump = Matrix::Zero(2, 2);
                        noJump(0, 0) = 1.;
                        noJump(1, 1) = std::sqrt(1. - gamma);
                        jump(0, 1) = std::sqrt(gamma);
                        kraus({noJump, jump}, q);
                        break;
                    }
                    case 6: // a measurement in a generic basis
                    {
                        const Eigen::VectorXcd v = Unitary(2, ++unitarySeed).col(0);
                        const Matrix p = v * v.adjoint();
                        kraus({p, Matrix(identity2 - p)}, q);
                        break;
                    }
                    default: // the site is projected, then meets a projector at delta from the orthogonal state
                    {
                        const double deltas[] = {1E-2, 1E-5, 1E-8};
                        const double delta = deltas[random() % 3];
                        const Matrix basis = Unitary(2, ++unitarySeed);
                        const Eigen::VectorXcd v = basis.col(0),
                                               w = std::cos(delta) * basis.col(1) + std::sin(delta) * basis.col(0);
                        kraus({Matrix(v * v.adjoint())}, q);
                        kraus({Matrix(w * w.adjoint())}, q);
                        tolerance = std::max(tolerance, 1E-13 / delta);
                        break;
                    }
                    }
                    const Eigen::VectorXcd state = mps.getRegisterStorage();
                    Close(state.squaredNorm(), 1., "An MPS step left the norm off, against a dense reference");
                    Require(distance(state, reference) < tolerance, "An MPS step differs from the dense reference");
                    // A read between the steps, taken from the chain itself.
                    Require(std::abs(mps.GetProbability(q, false) - probabilityOfOne(reference, q)) <
                                std::max(1E-12, 10. * tolerance),
                            "An MPS probability differs from the dense reference");
                }
            }
        }

        // An operator that has rank one only to rounding, diag(1, 1e-15), leaves a site in one state only if the site
        // holds enough of what it keeps. On |1> it leaves |1>, and on 1e-10 |0> + |1> it leaves both amplitudes, now
        // 1e-5 apart: the site must not be split again around |0> alone. Alone in a chain and between two others.
        for (const int qubits : {1, 3})
            for (const double slight : {0., 1E-10})
            {
                const int site = qubits == 1 ? 0 : 1;
                const double rest = std::sqrt(1. - slight * slight);
                Matrix prepare(2, 2), nearlyRankOne = Matrix::Zero(2, 2);
                prepare << slight, -rest, rest, slight;
                nearlyRankOne(0, 0) = 1.;
                nearlyRankOne(1, 1) = 1E-15;
                Impl mps(qubits);
                // the state of each qubit, the site last changed by the operator
                std::vector<Eigen::Vector2cd> states(qubits, Eigen::Vector2cd(1., 0.));
                for (int q = 0; q < qubits; ++q)
                    if (q != site)
                    {
                        const double angle = 0.7 + 0.2 * q;
                        mps.ApplyGate(Gate(turn(angle), q));
                        states[q] = Eigen::Vector2cd(std::cos(angle), std::sin(angle));
                    }
                mps.ApplyGate(Gate(prepare, site));
                states[site] = Eigen::Vector2cd(slight, 1E-15 * rest).normalized();
                mps.ApplyKrausOperator({nearlyRankOne}, 0.5, site);
                const Eigen::VectorXcd state = mps.getRegisterStorage();
                double worst = std::abs(state.squaredNorm() - 1.);
                for (Eigen::Index i = 0; i < state.size(); ++i)
                {
                    Complex expected = 1.;
                    for (int q = 0; q < qubits; ++q)
                        expected *= states[q][(i >> q) & 1];
                    worst = std::max(worst, std::abs(state[i] - expected));
                }
                Require(worst < 1E-12, "An MPS operator of rank one only to rounding erased what it left of the state");
            }

        // A chain of one qubit has no bond for the sweeps to normalize through: after an operator that leaves almost
        // nothing of the state, the tensor is normalized by itself. The operator |0><w| leaves |0> whatever the angle
        // of w from the state orthogonal to the qubit.
        for (const double delta : {1E-6, 1E-9, 1E-12})
        {
            Impl lone(1);
            lone.ApplyGate(Gate(turn(0.6), 0));
            Matrix nearlyNull = Matrix::Zero(2, 2);
            nearlyNull(0, 0) = std::sin(0.6 + delta);
            nearlyNull(0, 1) = -std::cos(0.6 + delta);
            lone.ApplyKrausOperator({nearlyNull}, 0.5, 0);
            const Eigen::VectorXcd state = lone.getRegisterStorage();
            Require(std::abs(state.squaredNorm() - 1.) < 1E-14 && std::abs(std::abs(state[0]) - 1.) < 1E-14 &&
                        std::abs(state[1]) < 1E-14,
                    "A near-null operator on an MPS of one qubit left the norm off");
        }

        // A weight that is read from the tensors and falls below its cutoff counts as none, and the message says what
        // was found. Two qubits in |+>, the row (1/2, -1/2 + 2^-49): what it leaves is at eight units of rounding.
        {
            Impl plus(2);
            plus.ApplyGate(Gate(hadamard, 0));
            plus.ApplyGate(Gate(hadamard, 1));
            Matrix faint(2, 2);
            faint << 0.5, -0.5 + std::ldexp(1., -49), 0., 0.;
            std::string message;
            try
            {
                plus.ApplyKrausOperator({faint}, 0.5, 0);
            }
            catch (const std::domain_error &error)
            {
                message = error.what();
            }
            Require(message.find("largest weight found") != std::string::npos,
                    "An MPS Kraus weight below its cutoff was not reported as found");
            Close(plus.getRegisterStorage().squaredNorm(), 1., "A refused MPS Kraus operator left the norm off");
        }

        std::cout << "Success" << std::endl;
        return true;
    }
    catch (const std::exception &error)
    {
        std::cout << "MPS Kraus/gauge regression: " << error.what() << std::endl;
        return false;
    }
}

bool MPSSimulatorTests()
{
    std::cout << "\nMPS Simulator Tests" << std::endl;

    /*
    QC::TensorNetworks::MPSSimulatorImpl mps(2);

    QC::Gates::PauliXGate xGate;
    QC::Gates::HadamardGate hGate;
    QC::Gates::SwapGate swapGate;

    mps.ApplyGate(xGate, 0);
    mps.ApplyGate(hGate, 1);

    mps.print();

    mps.ApplyGate(swapGate, 0, 1);

    mps.print();
    */

    /*
    for (int i = 0; i < 100; ++i)
    {
        if (!checkExpectationValuesMPS())
            return false;
    }
    */

    return QC::ProductSwapRegression::Run(false) && ObservableAndSummaryRegressionTestMPS() &&
           KrausGaugeModeRegressionTestMPS() && MeetingPositionFallbackTestMPS() && SwapRelabelTestMPS() &&
           WideBasisInitializationTestMPS() && StateSimulationTest() && NumericalRankStabilityTestMPS() &&
           checkExpectationValuesMPS() && TrimTestMPS() && TruncationModeTestMPS() && CloneTestMPS() &&
           ReCanonicalizeRegressionTestMPS() && CanonicalFormTestMPS() && MultithreadingSettingTestMPS() &&
           LongChainSamplingAndValidationTestMPS();
}
