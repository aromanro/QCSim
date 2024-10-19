#include "Clifford.h"
#include "Tests.h"

#include "QubitRegister.h"

#include <cstdlib>
#include <future>


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
		return std::make_shared<QC::Gates::CNOTGate<>>();
	case 9:
		return std::make_shared<QC::Gates::ControlledYGate<>>();
	case 10:
		return std::make_shared<QC::Gates::ControlledZGate<>>();
	case 11:
		return std::make_shared<QC::Gates::SwapGate<>>();
	}

	return nullptr;
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
		simulator.ApplyCX(qubit1, qubit2);
		break;
	case 9:
		simulator.ApplyCY(qubit1, qubit2);
		break;
	case 10:
		simulator.ApplyCZ(qubit1, qubit2);
		break;
	case 11:
		simulator.ApplySwap(qubit1, qubit2);
		break;
	}
}

bool CliffordSimulatorTests()
{
	const size_t nrTests = 100;
	const size_t nrShots = 50000;
	const size_t nrGates = 200;
	const size_t nrQubits = 6;
	const double errorThreshold = 0.02;

	std::uniform_int_distribution gateDistr(0, 11);
	std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);

	std::cout << "\nClifford gates simulator" << std::endl;

	for (size_t t = 0; t < nrTests; ++t)
	{
		std::unordered_map<size_t, int> results1;
		std::unordered_map<size_t, int> results2;

		// generate random gates, creating a circuit, then apply the random circuits on both simulators
		std::vector<int> gates(nrGates);
		std::vector<size_t> qubits1(nrGates);
		std::vector<size_t> qubits2(nrGates);

		for (int i = 0; i < nrGates; ++i)
		{
			gates[i] = gateDistr(gen);

			qubits1[i] = qubitDistr(gen);
			qubits2[i] = qubitDistr(gen);

			if (qubits2[i] == qubits1[i])
				qubits2[i] = (qubits1[i] + 1) % nrQubits;

			if (dist_bool(gen)) std::swap(qubits1[i], qubits2[i]);
		}


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
					for (int i = 0; i < curCnt; ++i)
					{
						QC::QubitRegister qubitRegister(nrQubits);
						QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

						for (int j = 0; j < gates.size(); ++j)
						{
							auto gate = gates[j];
							const auto qubit1 = qubits1[j];
							const auto qubit2 = qubits2[j];

							ApplyGate(cliffordSim, gate, qubit1, qubit2);
							const auto gateptr = GetGate(gate);
							qubitRegister.ApplyGate(*gateptr, qubit1, qubit2);
						}

						// now do the measurements
						size_t val1 = 0;
						size_t val2 = 0;
						for (int q = 0; q < nrQubits; ++q)
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
	std::cout << std::endl;

	std::cout << "Success" << std::endl;

	return true;
}