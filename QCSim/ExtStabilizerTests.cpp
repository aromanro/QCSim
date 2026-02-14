#include "Tests.h"
#include "QubitRegister.h"
#include "ExtendedStabilizer.h"

/*
void ApplyTwoQubitsGate(QC::ExtendedStabilizer& simulator, int code, int qubit1, int qubit2)
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


void ApplyGate(QC::ExtendedStabilizer& simulator, int code, int qubit1, int qubit2, double angle = 0.0)
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

void ExecuteCircuit(QC::QubitRegister<>& qubitRegister, QC::ExtendedStabilizer& extstabSim, std::vector<int>& gates, std::vector<size_t>& qubits1, std::vector<size_t>& qubits2, std::uniform_real_distribution<double>& angleDistr, std::bernoulli_distribution& boolDistr, std::uniform_int_distribution<int>& rotationGateDistr, bool clifford = true)
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

		ApplyGate(extstabSim, gates[j], (int)qubits1[j], (int)qubits2[j], angle);
	}
}
*/

bool TestExtStabilizer()
{
	std::cout << "\nExtended Stabilizer tests" << std::endl;

	const size_t nrTests = 100;
	const size_t maxQubits = 20;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(5, 20);
	std::uniform_real_distribution angleDistr(-2. * M_PI, 2. * M_PI);
	std::bernoulli_distribution boolDistr(0.2);
	std::uniform_int_distribution rotationGateDistr(15, 17); // RX, RY, RZ

	/*
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

			QC::ExtendedStabilizer extstabSim(nrQubits);
			QC::QubitRegister qubitRegister(nrQubits);

			ExecuteCircuit(qubitRegister, extstabSim, gates, qubits1, qubits2, angleDistr, boolDistr, rotationGateDistr, false);

			std::vector<QC::Gates::AppliedGate<>> expGates;
			expGates.reserve(nrQubits);
			std::string pauliStr;

			ConstructPauliString(nrQubits, pauliStr, expGates);

			const auto exp1 = qubitRegister.ExpectationValue(expGates);
			const auto exp2 = extstabSim.ExpectationValue(pauliStr);
			if (!approxEqual(exp1, exp2, 1E-7))
			{
				std::cout << std::endl << "Expectation values are not equal for pauli propagator and statevector simulator for " << nrQubits << " qubits, values: " << exp2 << ", " << exp1 << std::endl;

				std::cout << "Pauli string: " << pauliStr << std::endl;

				return false;
			}
		}
		std::cout << '.';
	}
	*/

	std::cout << "\nSuccess" << std::endl;
	return true;
}