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
    std::uniform_int_distribution nrGatesDistr(20, 40);

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

        piSim.SetCircuit(circuit);

        std::vector<bool> endState(nrQubits);
        for (size_t state = 0; state < nrBasisStates; ++state)
        {
            for (size_t i = 0; i < nrQubits; ++i)
                endState[i] = ((state >> i) & 1) == 1;

            const auto regAmpl = qubitRegister.getBasisStateAmplitude(state);

            const auto piAmpl = piSim.Propagate(endState);

            if (!approxEqual(regAmpl, piAmpl, 1E-7))
            {
                std::cout << std::endl << "Amplitude mismatch for path integral and statevector simulator for state " << state << " with " << nrQubits << " qubits" << std::endl;
                std::cout << "Statevector: " << regAmpl << " vs Path integral: " << piAmpl << std::endl;
                return false;
            }
        }

        // the above went with the paths from end state towards the 'middle' and from start state to meet them in the middle,
        // for qubits probabilities must go from start to end with all of them
        piSim.PropagateAll(circuit);

        for (size_t q = 0; q < nrQubits; ++q)
        {
			const auto qubitProb = qubitRegister.GetQubitProbability(q);
			const auto piQubitProb = piSim.QubitProbability(q);

            if (!approxEqual(qubitProb, piQubitProb, 1E-7))
            {
                std::cout << std::endl << "Qubit probability mismatch for path integral and statevector simulator for qubit " << q << " with " << nrQubits << " qubits" << std::endl;
                std::cout << "Statevector: " << qubitProb << " vs Path integral: " << piQubitProb << std::endl;
				return false;
            }
        }

        std::cout << '.';
    }

    std::cout << "\nSuccess" << std::endl;

    return true;
}
