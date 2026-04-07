#include "Tests.h"

#include "QubitRegister.h"
#include "PathIntegral.h"

bool TestPathIntegral()
{
    std::cout << "\nPath Integral tests" << std::endl;

    const size_t nrQubits = 10;
    const size_t nrBasisStates = 1ULL << nrQubits;
    const size_t nrTests = 100;

    std::uniform_int_distribution gateDistr(0, 17);
    std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);
    std::uniform_int_distribution nrGatesDistr(1, 30);

    QC::PathIntegral::PathIntegralSimulator piSim;

    for (size_t t = 0; t < nrTests; ++t)
    {
        const size_t nrGates = nrGatesDistr(gen);
        std::vector<int> gates(nrGates);
        std::vector<size_t> qubits1(nrGates);
        std::vector<size_t> qubits2(nrGates);

        ConstructCircuit(nrQubits, gates, qubits1, qubits2, gateDistr, qubitDistr);

        QC::QubitRegister qubitRegister(nrQubits);

        std::vector<QC::Gates::AppliedGate<>> circuit;
        circuit.reserve(nrGates);

        for (size_t j = 0; j < nrGates; ++j)
        {
            auto gateptr = GetGate(gates[j]);
            qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
            circuit.emplace_back(gateptr->getRawOperatorMatrix(), qubits1[j], qubits2[j]);
        }

        std::vector<bool> endState(nrQubits);
        for (size_t state = 0; state < nrBasisStates; ++state)
        {
            for (size_t i = 0; i < nrQubits; ++i)
                endState[i] = ((state >> i) & 1) == 1;

            const auto regAmpl = qubitRegister.getBasisStateAmplitude(state);
            const auto piAmpl = piSim.Propagate(circuit, endState);

            if (!approxEqual(regAmpl, piAmpl, 1E-7))
            {
                std::cout << std::endl << "Amplitude mismatch for path integral and statevector simulator for state " << state << " with " << nrQubits << " qubits" << std::endl;
                std::cout << "Statevector: " << regAmpl << " vs Path integral: " << piAmpl << std::endl;
                return false;
            }
        }

        std::cout << '.';
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}
