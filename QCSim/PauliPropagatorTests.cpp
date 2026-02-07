#include "Tests.h"

#include "QubitRegister.h"
#include "PauliPropagator.h"


void ApplyTwoQubitsGate(QC::PauliPropagator& simulator, int code, int qubit1, int qubit2)
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


void ApplyGate(QC::PauliPropagator& simulator, int code, int qubit1, int qubit2, double angle = 0.0)
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
	case 15:
		simulator.ApplyRX(qubit1, angle);
		break;
	case 16:
		simulator.ApplyRY(qubit1, angle);
		break;
	case 17:
		simulator.ApplyRZ(qubit1, angle);
		break;
	default:
		ApplyTwoQubitsGate(simulator, code, qubit1, qubit2);
		break;
	}
}

bool CheckResults(int shots, int nrQubits, std::unordered_map<size_t, size_t>& sampledResultsPauli, const std::unordered_map<size_t, size_t>& sampledResultsStatevector, bool measured = false)
{
	for (const auto& kv : sampledResultsStatevector)
	{
		const size_t state = kv.first;
		const size_t countStatevector = kv.second;
		const size_t countPauli = sampledResultsPauli[state];
		const double freqStatevector = static_cast<double>(countStatevector) / static_cast<double>(shots);
		const double freqPauli = static_cast<double>(countPauli) / static_cast<double>(shots);

		if (std::abs(freqStatevector - freqPauli) > 0.1)
		{
			std::cout << std::endl << "Sampled frequencies - " << (measured ? "by measurement" : "by sampling") << " - are not equal for pauli propagator and statevector simulator for " << nrQubits << " qubits, measurement " << state << ", values: " << freqPauli << ", " << freqStatevector << std::endl;
			return false;
		}
	}

	return true;
}

bool CheckProbability(int nrQubits, QC::PauliPropagator& pauliSimulator, QC::QubitRegister<>& qubitRegister)
{
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

	return true;
}

void ExecuteCircuit(QC::QubitRegister<>& qubitRegister, QC::PauliPropagator& pauliSimulator, std::vector<int>& gates, std::vector<size_t>& qubits1, std::vector<size_t>& qubits2, std::uniform_real_distribution<double>& angleDistr, std::bernoulli_distribution& boolDistr, std::uniform_int_distribution<int>& rotationGateDistr, bool clifford)
{
	for (int j = 0; j < static_cast<int>(gates.size()); ++j)
	{
		auto gateptr = GetGate(gates[j]);
		double angle = 0.0;
		if (!clifford && boolDistr(gen))
		{
			gates[j] = rotationGateDistr(gen);
			angle = angleDistr(gen);
			gateptr = GetGate(gates[j], angle);
		}

		qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);

		ApplyGate(pauliSimulator, gates[j], (int)qubits1[j], (int)qubits2[j], angle);
	}
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
	std::uniform_int_distribution rotationGateDistr(15, 17); // RX, RY, RZ

	QC::PauliPropagator pauliSimulator;
	
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

			ExecuteCircuit(qubitRegister, pauliSimulator, gates, qubits1, qubits2, angleDistr, boolDistr, rotationGateDistr, clifford);

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

			if (!CheckProbability(nrQubits, pauliSimulator, qubitRegister))
				return false;

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

				if (!CheckResults(shots, (int)nrQubits, sampledResultsPauli, sampledResultsStatevector, true))
					return false;

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
				if (!CheckResults(shots, (int)nrQubits, sampledResultsPauli, sampledResultsStatevector, true))
					return false;
			}

			pauliSimulator.ClearOperations();
		}
		std::cout << '.';
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}


bool TestPauliPropagator()
{
	return TestPauliPropagatorCorNC(true) && TestPauliPropagatorCorNC(false);
}
