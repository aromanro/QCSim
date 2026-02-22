#include "Tests.h"
#include "QubitRegister.h"
#include "ExtendedStabilizer.h"


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
	/*
	case 15:
		simulator.ApplyRx(qubit1, angle);
		break;
	case 16:
		simulator.ApplyRy(qubit1, angle);
		break;
	case 17:
		simulator.ApplyRz(qubit1, angle);
		break;
	*/	
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


static bool TestExtStabilizerMeasurements()
{
	std::cout << "\nExtended Stabilizer measurement tests" << std::endl;

	const size_t nrShots = 100000;
	const double errorThreshold = 0.01;
	const size_t nrTests = 10;
	const size_t maxQubits = 8;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(50, 100);
	std::uniform_real_distribution angleDistr(-2. * M_PI, 2. * M_PI);
	std::bernoulli_distribution boolDistr(0.2);
	std::uniform_int_distribution rotationGateDistr(15, 17);

	for (size_t nrQubits = 2; nrQubits < maxQubits; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);

		for (size_t t = 0; t < nrTests; ++t)
		{
			const size_t nrGates = nrGatesDistr(gen);
			std::vector<int> gates(nrGates);
			std::vector<size_t> qubits1(nrGates);
			std::vector<size_t> qubits2(nrGates);

			ConstructCircuit(nrQubits, gates, qubits1, qubits2, gateDistr, qubitDistr);

			// test 1: repeated measurement on the same qubit must be consistent
			{
				QC::ExtendedStabilizer sim(nrQubits);
				QC::QubitRegister<> reg(nrQubits);
				ExecuteCircuit(reg, sim, gates, qubits1, qubits2, angleDistr, boolDistr, rotationGateDistr);

				for (size_t q = 0; q < nrQubits; ++q)
				{
					const bool res1 = sim.Measure(q);
					const bool res2 = sim.Measure(q);
					if (res1 != res2)
					{
						std::cout << std::endl << "Repeated measurement inconsistency for qubit " << q << " with " << nrQubits << " qubits" << std::endl;
						return false;
					}
				}
			}

			// test 2: sampling comparison against statevector
			{
				std::unordered_map<size_t, size_t> stabResults;
				std::unordered_map<size_t, size_t> svResults;

				QC::QubitRegister<> reg(nrQubits);
				QC::ExtendedStabilizer sim(nrQubits);

				ExecuteCircuit(reg, sim, gates, qubits1, qubits2, angleDistr, boolDistr, rotationGateDistr);

				svResults = reg.RepeatedMeasureUnordered(nrShots);

				sim.SaveState();
				for (size_t shot = 0; shot < nrShots; ++shot)
				{
					size_t stabVal = 0;
					for (size_t q = 0; q < nrQubits; ++q)
						if (sim.Measure(q)) stabVal |= 1 << q;

					++stabResults[stabVal];

					sim.RestoreState();
				}

				for (const auto& val : svResults)
				{
					const double svFreq = static_cast<double>(val.second) / nrShots;
					const double stabFreq = stabResults.count(val.first) ? static_cast<double>(stabResults[val.first]) / nrShots : 0.0;

					if (std::abs(svFreq - stabFreq) > errorThreshold)
					{
						std::cout << std::endl << "Measurement distribution mismatch for " << nrQubits << " qubits, state " << val.first
							<< ": statevector " << svFreq << ", stabilizer " << stabFreq << std::endl;
						std::cout << "Might fail due to randomness of measurements" << std::endl;
						return false;
					}
				}

				for (const auto& val : stabResults)
				{
					if (svResults.count(val.first)) continue;

					const double stabFreq = static_cast<double>(val.second) / nrShots;
					if (stabFreq > errorThreshold)
					{
						std::cout << std::endl << "Stabilizer produced state " << val.first << " with frequency " << stabFreq
							<< " but statevector never did, for " << nrQubits << " qubits" << std::endl;
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

bool TestExtStabilizer()
{
	std::cout << "\nExtended Stabilizer tests" << std::endl;

	// basic sanity check: repeated rotations must compose correctly
	/*
	{
		constexpr size_t nrQubits = 1;
		QC::ExtendedStabilizer extstabSim(nrQubits);
		QC::QubitRegister<> qubitRegister(nrQubits);

		const double theta = M_PI / 2.0;
		auto rx = GetGate(15, theta);
		qubitRegister.ApplyGate(*rx, 0, 0);
		extstabSim.ApplyRx(0, theta);
		qubitRegister.ApplyGate(*rx, 0, 0);
		extstabSim.ApplyRx(0, theta);

		const double p1 = qubitRegister.GetQubitProbability(0);
		const double p2 = extstabSim.GetQubitProbability(0);
		if (!approxEqual(p1, p2, 1E-10))
		{
			std::cout << "\nRotation composition failed for Rx on 1 qubit, values: " << p1 << ", " << p2 << std::endl;
			return false;
		}
	}
	*/

	const size_t nrTests = 100;
	const size_t maxQubits = 20;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(50, 100);
	std::uniform_real_distribution angleDistr(-2. * M_PI, 2. * M_PI);
	std::bernoulli_distribution boolDistr(0.2);
	std::uniform_int_distribution rotationGateDistr(15, 17); // RX, RY, RZ

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

			ExecuteCircuit(qubitRegister, extstabSim, gates, qubits1, qubits2, angleDistr, boolDistr, rotationGateDistr);

			for (size_t q = 0; q < nrQubits; ++q)
			{
				double p1 = qubitRegister.GetQubitProbability(q);
				double p2 = extstabSim.GetQubitProbability(q);
				if (!approxEqual(p1, p2, 1E-5))
				{
					std::cout << std::endl << "Probabilities are not equal for statevector and stabilizer simulator for " << nrQubits << " qubits, values: " << p1 << ", " << p2 << std::endl;
					return false;
				}
			}

			
			std::vector<QC::Gates::AppliedGate<>> expGates;
			expGates.reserve(nrQubits);
			std::string pauliStr;

			ConstructPauliString(nrQubits, pauliStr, expGates);

			const auto exp1 = qubitRegister.ExpectationValue(expGates).real();
			const auto exp2 = extstabSim.ExpectationValue(pauliStr);
			if (!approxEqual(exp1, exp2, 1E-7))
			{
				std::cout << std::endl << "Expectation values are not equal for statevector and stabilizer simulator for " << nrQubits << " qubits, values: " << exp1 << ", " << exp2 << std::endl;

				std::cout << "Pauli string: " << pauliStr << std::endl;

				return false;
			}

			//extstabSim.frames.begin()->Print();
			//std::cout << "****************************************************************" << std::endl;
		}
		std::cout << '.';
	}

	std::cout << "\nSuccess" << std::endl;

	if (!TestExtStabilizerMeasurements())
		return false;

	return true;
}