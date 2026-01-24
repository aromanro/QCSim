#include "Tests.h"

#include "QubitRegister.h"
#include "PauliPropagator.h"


void ApplyTwoQubitsGate(QC::PauliPropagator& simulator, int code, int qubit1, int qubit2)
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


void ApplyGate(QC::PauliPropagator& simulator, int code, int qubit1, int qubit2)
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
	default:
		ApplyTwoQubitsGate(simulator, code, qubit1, qubit2);
		break;
	}
}

bool TestPauliPropagator()
{
	std::cout << "\nPauli Propagator tests" << std::endl;
	const size_t nrTests = 100;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(5, 20);
	
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

			QC::PauliPropagator pauliSimulator;
			pauliSimulator.SetNrQubits(static_cast<int>(nrQubits));

			QC::QubitRegister qubitRegister(nrQubits);
			for (int j = 0; j < static_cast<int>(gates.size()); ++j)
			{
				const auto gateptr = GetGate(gates[j]);
				qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);

				ApplyGate(pauliSimulator, gates[j], qubits1[j], qubits2[j]);
			}

			std::vector<QC::Gates::AppliedGate<>> expGates;
			expGates.reserve(nrQubits);
			std::string pauliStr;

			ConstructPauliString(nrQubits, pauliStr, expGates);

			const auto exp1 = qubitRegister.ExpectationValue(expGates);
			const auto exp2 = pauliSimulator.ExpectationValue(pauliStr);
			if (!approxEqual(exp1, exp2, 1E-7))
			{
				std::cout << std::endl << "Expectation values are not equal for pauli propagator and statevector simulator for " << nrQubits << " qubits, values: " << exp2 << ", " << exp1 << std::endl;

				std::cout << "Pauli string: " << pauliStr << std::endl;

				return false;
			}

			for (size_t q = 0; q < nrQubits; ++q)
			{
				const auto p0 = qubitRegister.GetQubitProbability(q);
				const auto p1 = 1. - pauliSimulator.Probability0(static_cast<int>(q));
				if (std::abs(p0 - p1) > 1E-7)
				{
					std::cout << std::endl << "Probability values are not equal for pauli propagator and statevector simulator for " << nrQubits << " qubits, qubit " << q << ", values: " << p1 << ", " << p0 << std::endl;
					return false;
				}
			}

			// for big number of qubits sampling on the Pauli propagator becomes too slow
			if (nrQubits < 9 && t < 10) // also limit the number of tests, this can become very slow
			{
				std::vector<int> measQubits(nrQubits);
				std::iota(measQubits.begin(), measQubits.end(), 0);
				
				int shots = 1000;
				auto sampledResultsStatevector = qubitRegister.RepeatedMeasureUnordered(shots);
				
				std::unordered_map<size_t, size_t> sampledResultsPauli;	
				for (int j = 0; j < shots; ++j)
				{
					std::shuffle(measQubits.begin(), measQubits.end(), gen);
					auto res = pauliSimulator.Sample(measQubits);

					size_t result = 0;
					for (size_t q = 0; q < nrQubits; ++q)
					{
						if (res[q])
						{
							size_t qubit = measQubits[q];
							result |= (1ULL << qubit);
						}
					}

					++sampledResultsPauli[result];
				}

				// now check the results
				for (const auto& kv : sampledResultsStatevector)
				{
					const size_t state = kv.first;
					const size_t countStatevector = kv.second;
					const size_t countPauli = sampledResultsPauli[state];
					const double freqStatevector = static_cast<double>(countStatevector) / static_cast<double>(shots);
					const double freqPauli = static_cast<double>(countPauli) / static_cast<double>(shots);

					if (std::abs(freqStatevector - freqPauli) > 0.1)
					{
						std::cout << std::endl << "Sampled frequencies are not equal for pauli propagator and statevector simulator for " << nrQubits << " qubits, measurement " << state << ", values: " << freqPauli << ", " << freqStatevector << std::endl;
						return false;
					}
				}

				// now do the same but with measurements!
				sampledResultsPauli.clear();

				pauliSimulator.SaveState();
				for (int j = 0; j < shots; ++j)
				{
					std::shuffle(measQubits.begin(), measQubits.end(), gen);
				
					auto res = pauliSimulator.Measure(measQubits);

					size_t result = 0;
					for (size_t q = 0; q < nrQubits; ++q)
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

				// now check the results
				for (const auto& kv : sampledResultsStatevector)
				{
					const size_t state = kv.first;
					const size_t countStatevector = kv.second;
					const size_t countPauli = sampledResultsPauli[state];
					const double freqStatevector = static_cast<double>(countStatevector) / static_cast<double>(shots);
					const double freqPauli = static_cast<double>(countPauli) / static_cast<double>(shots);

					if (std::abs(freqStatevector - freqPauli) > 0.1)
					{
						std::cout << std::endl << "Sampled - by measurement - frequencies are not equal for pauli propagator and statevector simulator for " << nrQubits << " qubits, measurement " << state << ", values: " << freqPauli << ", " << freqStatevector << std::endl;
						return false;
					}
				}
			}
		}
		std::cout << '.';
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}