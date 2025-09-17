#include "Clifford.h"
#include "Tests.h"

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
	}

	return nullptr;
}


std::shared_ptr<QC::Gates::QuantumGateWithOp<>> GetGate(int code)
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
	default:
		return GetTwoQubitsGate(code);
	}

	return nullptr;
}

void ApplyTwoQubitsGate(QC::Clifford::StabilizerSimulator& simulator, int code, int qubit1, int qubit2)
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
	}
}


void ApplyGate(QC::Clifford::StabilizerSimulator& simulator, int code, int qubit1, int qubit2)
{
	switch(code)
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


void ConstructCircuit(size_t nrQubits, std::vector<int>& gates, std::vector<size_t>& qubits1, std::vector<size_t>& qubits2, std::uniform_int_distribution<int>& gateDistr, std::uniform_int_distribution<int>& qubitDistr)
{
	for (int i = 0; i < static_cast<int>(gates.size()); ++i)
	{
		gates[i] = gateDistr(gen);

		qubits1[i] = qubitDistr(gen);
		qubits2[i] = qubitDistr(gen);

		if (qubits2[i] == qubits1[i])
			qubits2[i] = (qubits1[i] + 1) % static_cast<int>(nrQubits);

		if (dist_bool(gen)) std::swap(qubits1[i], qubits2[i]);
	}
}

void ExecuteCircuit(size_t nrShots, size_t nrQubits, const std::vector<int>& gates, const std::vector<size_t>& qubits1, const std::vector<size_t>& qubits2, std::unordered_map<size_t, int>& results1, std::unordered_map<size_t, int>& results2)
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

		tasks[th] = std::async(std::launch::async, [&gates, &qubits1, &qubits2, &results1, &results2, curCnt, nrQubits, &resultsMutex]()
			{
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

						if (cliffordSim.MeasureQubit(q)) val1 |= 1;
						if (qubitRegister.MeasureQubit(q)) val2 |= 1;
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

bool CheckProbability(QC::QubitRegister<>& qubitRegister, QC::Clifford::StabilizerSimulator& cliffordSim)
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

bool CheckAllStatesProbability(QC::QubitRegister<>& qubitRegister, QC::Clifford::StabilizerSimulator& cliffordSim)
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
			std::cout << "Probability statevector: " << prob2 << ", Probability stabilizer: " << prob1 << ", State: " << state << std::endl;
			return false;
		}
	}

	return true;
}

bool CheckMeasurements(QC::Clifford::StabilizerSimulator& cliffordSim)
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

bool CliffordSimulatorTests()
{
	const size_t nrTests = 10;
	const size_t nrShots = 100000;
	const double errorThreshold = 0.01;

	std::uniform_int_distribution gateDistr(0, 13);
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
				if (results2.find(val.first) == results2.end()) continue;

				if (std::abs(static_cast<double>(val.second) - results2[val.first]) / nrShots > errorThreshold)
				{
					std::cout << "\nFailed" << std::endl;
					std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
					std::cout << "Result 1: " << static_cast<double>(val.second) / nrShots << ", Result 2: " << static_cast<double>(results2[val.first]) / nrShots << std::endl;
					return false;
				}
			}

			for (const auto val : results2)
			{
				if (results1.find(val.first) == results1.end()) continue;

				if (std::abs(static_cast<double>(val.second) - results1[val.first]) / nrShots > errorThreshold)
				{
					std::cout << "\nFailed" << std::endl;
					std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
					std::cout << "Result 1: " << static_cast<double>(results1[val.first]) / nrShots << ", Result 2: " << static_cast<double>(val.second) / nrShots << std::endl;
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

void ConstructPauliString(size_t nrQubits, std::string& pauliStr, std::vector<QC::Gates::AppliedGate<>>& expGates)
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

	std::uniform_int_distribution gateDistr(0, 13);
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
				std::cout << std::endl << "Expectation values are not equal for stabilizer and statevector simulator for " << nrQubits << " qubits, values: " << exp1 << ", " << exp2 << std::endl;

				return false;
			}
		}
		std::cout << '.';
	}

	std::cout << std::endl;

	std::cout << "Success" << std::endl;

	return true;
}
