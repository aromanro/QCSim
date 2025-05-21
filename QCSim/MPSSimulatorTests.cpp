#include <iostream>
#include <iterator>
#include <memory>


#include "Tests.h"

#include "QubitRegister.h"
#include "MPSSimulator.h"

#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>
#include <future>


void FillOneQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::HyGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SDGGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::TGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::TDGGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PhaseShiftGate<>>(0.5 * M_PI));
	gates.emplace_back(std::make_shared<QC::Gates::PauliXGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PauliYGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PauliZGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTDagGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SplitterGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::RxGate<>>(M_PI / 4));
	gates.emplace_back(std::make_shared<QC::Gates::RyGate<>>(M_PI / 4));
	gates.emplace_back(std::make_shared<QC::Gates::RzGate<>>(M_PI / 4));
	gates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 4, M_PI / 8));
}

void FillTwoQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
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
	gates.emplace_back(std::make_shared<QC::Gates::ControlledUGate<>>(M_PI / 4, M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRxGate<>>(M_PI / 5));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRyGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRzGate<>>(M_PI / 4));
}

bool OneQubitGatesTest()
{
	std::cout << "\nMPS simulator state test for one qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(50, 100);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 1; nrQubits < 9; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPSSimulatorImpl mps(nrQubits);
			QC::QubitRegister reg(nrQubits);

			const int lim = nrGatesDistr(gen);

			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const int qubit = qubitDistr(gen);

				mps.ApplyGate(*gates[gate], qubit);
				reg.ApplyGate(*gates[gate], qubit);
			}

			// now check the results, they should be the same
			const auto& regState = reg.getRegisterStorage();
			const auto mpsState = mps.getRegisterStorage(); // this one is computed, returns value, not reference, not stored elsewhere

			for (int i = 0; i < regState.size(); ++i)
			{
				if (!approxEqual(regState[i], mpsState[i]))
				{
					std::cout << "State simulation test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;
					std::cout << "Reg state:\n" << regState << std::endl;
					std::cout << "MPS state:\n" << mpsState << std::endl;
					return false;
				}
			}
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

bool OneAndTwoQubitGatesTest()
{
	std::cout << "\nMPS simulator state test for both one and two qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);
	FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);


	for (int nrQubits = 2; nrQubits < 9; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
#ifdef _DEBUG
			std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl << std::endl << std::endl;
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

				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

#ifdef _DEBUG
				if (twoQubitsGate) std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2 << std::endl;
#endif

				mps.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);

				// now check the results, they should be the same
				const auto& regState = reg.getRegisterStorage();
				auto mpsState = mps.getRegisterStorage(); // this one is computed, returns value, not reference, not stored elsewhere

				//QC::QubitRegister regNorm(nrQubits);
				//regNorm.setRegisterStorage(mpsState);
				//mpsState = regNorm.getRegisterStorage();

				for (int s = 0; s < regState.size(); ++s)
				{
					if (!approxEqual(regState[s], mpsState[s], 1E-7))
					{
						std::cout << "State " << s << " simulation test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;

						std::cout << "Probability for the different states: " << std::norm(regState[s]) << " vs " << std::norm(mpsState[s]) << std::endl;

						std::cout << "Reg state:\n" << regState << std::endl;
						std::cout << "Reg state normalization: " << regState.norm() << std::endl;

						std::cout << "MPS state:\n" << mpsState << std::endl;
						std::cout << "MPS state normalization: " << mpsState.norm() << std::endl;

						std::cout << std::endl;
						for (int q = 0; q < nrQubits; ++q)
							std::cout << "Qubit " << q << " reg probability: " << reg.GetQubitProbability(q) << " vs mps: " << mps.GetProbability(q, false) << std::endl;

						std::cout << std::endl;
						for (int state = 0; state < regState.size(); ++state)
							std::cout << "State " << state << " reg probability: " << std::norm(regState[state]) << " vs mps: " << std::norm(mpsState[state]) << std::endl;

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

std::vector<std::shared_ptr<QC::Gates::AppliedGate<>>> GenerateRandomCircuitWithGates(const std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates, int minGates, int maxGates, int nrQubits)
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
		if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

		appliedGate->setQubit1(qubit1);
		appliedGate->setQubit2(qubit2);
		
		circuit.emplace_back(std::move(appliedGate));
	}

	return circuit;
}

std::vector<std::shared_ptr<QC::Gates::AppliedGate<>>> GenerateRandomCircuitWithGatesNoAdjacent(const std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates, int minGates, int maxGates, int nrQubits)
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
		if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

		appliedGate->setQubit1(qubit1);
		appliedGate->setQubit2(qubit2);

		circuit.emplace_back(std::move(appliedGate));
	}

	return circuit;
}

bool TestMeasurementsWithOneQubitGatesCircuits()
{
	std::cout << "\nMPS simulator measurements test with circuits with one qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);

	const int nrMeasurements = 10000;

	for (int nrQubits = 1; nrQubits < 9; ++nrQubits)
	{
		for (int t = 0; t < 5; ++t)
		{
			const auto circuit = GenerateRandomCircuitWithGates(gates, 5, 10, nrQubits);
			for (int c = 0; c < 3; ++c)
			{
				std::unordered_map<std::vector<bool>, int> measurementsRegMap;
				std::unordered_map<std::vector<bool>, int> measurementsMPSMap;

				for (int i = 0; i < nrMeasurements; ++i)
				{
					QC::TensorNetworks::MPSSimulatorImpl mps(nrQubits);
					QC::QubitRegister reg(nrQubits);

					for (const auto& gate : circuit)
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

				std::cout << ".";

				for (const auto& [key, value] : measurementsRegMap)
				{
					const double dif = abs((static_cast<double>(measurementsMPSMap[key] - value)) / nrMeasurements);
					if (dif > 0.05)
					{
						std::cout << "Measurements test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;
						std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
						std::cout << "Difference: " << dif << std::endl;

						std::cout << "Reg measurements:\n";
						for (const auto& [key, value] : measurementsRegMap)
						{
							for (const auto& b : key)
								std::cout << b << " ";
							std::cout << " : " << static_cast<double>(value) / nrMeasurements << std::endl;
						}

						std::cout << "MPS measurements:\n";
						for (const auto& [key, value] : measurementsMPSMap)
						{
							for (const auto& b : key)
								std::cout << b << " ";
							std::cout << " : " << static_cast<double>(value) / nrMeasurements << std::endl;
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

bool TestMeasurementsWithOneAndTwoQubitGatesCircuits()
{
	std::cout << "\nMPS simulator measurements test with circuits with both one and two qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);
	FillTwoQubitGates(gates);

	const int nrMeasurements = 10000;

	for (int nrQubits = 2; nrQubits < 9; ++nrQubits)
	{
		for (int t = 0; t < 5; ++t)
		{
			const auto circuit = GenerateRandomCircuitWithGates(gates, 25, 50, nrQubits);
			for (int c = 0; c < 3; ++c)
			{
				std::unordered_map<std::vector<bool>, int> measurementsRegMap;
				std::unordered_map<std::vector<bool>, int> measurementsMPSMap;

				for (int t = 0; t < nrMeasurements; ++t)
				{
					QC::TensorNetworks::MPSSimulatorImpl mps(nrQubits);
					QC::QubitRegister reg(nrQubits);

					for (const auto& gate : circuit)
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

				std::cout << ".";

				for (const auto& [key, value] : measurementsRegMap)
				{
					const double dif = abs((static_cast<double>(measurementsMPSMap[key] - value)) / nrMeasurements);
					if (dif > 0.05)
					{
						std::cout << "Measurements test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;
						std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
						std::cout << "Difference: " << dif << std::endl;

						std::cout << "Reg measurements:\n";
						for (const auto& [key, value] : measurementsRegMap)
						{
							for (const auto& b : key)
								std::cout << b << " ";
							std::cout << " : " << static_cast<double>(value) / nrMeasurements << std::endl;
						}

						std::cout << "MPS measurements:\n";
						for (const auto& [key, value] : measurementsMPSMap)
						{
							for (const auto& b : key)
								std::cout << b << " ";
							std::cout << " : " << static_cast<double>(value) / nrMeasurements << std::endl;
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

	for (int nrQubits = 3; nrQubits < 9; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
#ifdef _DEBUG
			std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl << std::endl << std::endl;
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
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

#ifdef _DEBUG
				if (twoQubitsGate) std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2 << std::endl;
#endif

				mps.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);

				// now check the results, they should be the same
				const auto& regState = reg.getRegisterStorage();
				auto mpsState = mps.getRegisterStorage(); // this one is computed, returns value, not reference, not stored elsewhere

				//QC::QubitRegister regNorm(nrQubits);
				//regNorm.setRegisterStorage(mpsState);
				//mpsState = regNorm.getRegisterStorage();

				for (int s = 0; s < regState.size(); ++s)
				{
					if (!approxEqual(regState[s], mpsState[s], 1E-7))
					{
						std::cout << "State " << s << " simulation test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;

						std::cout << "Probability for the different states: " << std::norm(regState[s]) << " vs " << std::norm(mpsState[s]) << std::endl;

						std::cout << "Reg state:\n" << regState << std::endl;
						std::cout << "Reg state normalization: " << regState.norm() << std::endl;

						std::cout << "MPS state:\n" << mpsState << std::endl;
						std::cout << "MPS state normalization: " << mpsState.norm() << std::endl;

						std::cout << std::endl;
						for (int q = 0; q < nrQubits; ++q)
							std::cout << "Qubit " << q << " reg probability: " << reg.GetQubitProbability(q) << " vs mps: " << mps.GetProbability(q, false) << std::endl;

						std::cout << std::endl;
						for (int state = 0; state < regState.size(); ++state)
							std::cout << "State " << state << " reg probability: " << std::norm(regState[state]) << " vs mps: " << std::norm(mpsState[state]) << std::endl;

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

void PrintMeasurements(int nrMeasurements, const std::unordered_map<std::vector<bool>, int>& measurements)
{
	for (const auto& [k, v] : measurements)
	{
		for (const auto& b : k)
			std::cout << b << " ";
		std::cout << " : " << static_cast<double>(v) / nrMeasurements << std::endl;
	}
}

bool CheckMeasurements(int nrQubits, int nrMeasurements, std::unordered_map<std::vector<bool>, int>& measurementsRegMap, std::unordered_map<std::vector<bool>, int>& measurementsMPSMap, std::unordered_map<std::vector<bool>, int>& measurementsMPSMapOpt, std::unordered_map<std::vector<bool>, int>& measurementsMPSMapAll)
{
	static const double threshold = 0.05;

	for (const auto& [key, value] : measurementsRegMap)
	{
		const double dif1 = abs((static_cast<double>(measurementsMPSMap[key] - value)) / nrMeasurements);
		const double dif2 = abs((static_cast<double>(measurementsMPSMapOpt[key] - value)) / nrMeasurements);
		const double dif3 = abs((static_cast<double>(measurementsMPSMapAll[key] - value)) / nrMeasurements);
		if (dif1 > threshold || dif2 > threshold || dif3 > threshold)
		{
			std::cout << "Measurements test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;
			std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
			std::cout << "Difference 1: " << dif1 << " Difference 2: " << dif2 << " Difference 3: " << dif3 << std::endl;

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
	std::cout << "\nMPS swapping/mapped simulator measurements test with circuits with both one and two qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);
	FillTwoQubitGates(gates);

	const int nrMeasurements = 10000;

	for (int nrQubits = 3; nrQubits < 9; ++nrQubits)
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

					tasks[i] = std::async(std::launch::async, [&circuit, &measurementsRegMap, &measurementsMPSMap, &measurementsMPSMapOpt, &measurementsMPSMapAll, curCnt, nrQubits, &resultsMutex]()
						{
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
								for (const auto& gate : circuit)
								{
									mps.ApplyGate(*gate);
									mpsOpt.ApplyGate(*gate);
									if (j == 0) mpsAll.ApplyGate(*gate); // apply it only the first iteration, measure no-collapse is used
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
								
								{
									const std::lock_guard lock(resultsMutex);
									++measurementsRegMap[measurementsReg];
									++measurementsMPSMap[measurementsMPS];
									++measurementsMPSMapOpt[measurementsMPSOpt];
									++measurementsMPSMapAll[measurementsMPSAll];
								}

								reg.setToBasisState(0);
								mps.setToBasisState(0);
								mpsOpt.setToBasisState(0);
								//mpsAll.setToBasisState(0); // no need to reset if only no-collapse measurements are used and the gates are applied only the first iteration
							}
						});
				}

				for (size_t i = 0; i < nrThreads; ++i)
					tasks[i].get();

				std::cout << ".";

				const bool res = CheckMeasurements(nrQubits, nrMeasurements, measurementsRegMap, measurementsMPSMap, measurementsMPSMapOpt, measurementsMPSMapAll);

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

	for (int nrQubits = 3; nrQubits < 9; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
#ifdef _DEBUG
			std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl << std::endl << std::endl;
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
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

#ifdef _DEBUG
				if (twoQubitsGate) std::cout << "Applying two qubit gate " << gate << " on qubits " << qubit1 << " and " << qubit2 << std::endl;
#endif

				mps.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);

				for (int s = 0; s < 5; ++s)
				{
					const size_t state = GenerateRandomState(nrQubits);

					const auto regAmpl = reg.getBasisStateAmplitude(state);
					const auto mpsAmpl = mps.getBasisStateAmplitude(state);

					if (!approxEqual(regAmpl, mpsAmpl, 1E-7))
					{
						std::cout << "State " << state << " simulation test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;

						std::cout << "Amplitude for the different states: " << regAmpl << " vs " << mpsAmpl << std::endl;

						return false;
					}

					const auto regProb = reg.getBasisStateProbability(state);
					const auto mpsProb = mps.getBasisStateProbability(state);

					if (!approxEqual(regProb, mpsProb, 1E-7))
					{
						std::cout << "State " << state << " simulation test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;

						std::cout << "Probability for the different states: " << regProb << " vs " << mpsProb << std::endl;

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

bool StateSimulationTest()
{
	return /*OneQubitGatesTest() && OneAndTwoQubitGatesTest() &&
		TestMeasurementsWithOneQubitGatesCircuits() && TestMeasurementsWithOneAndTwoQubitGatesCircuits() &&*/
		OneAndTwoQubitGatesTestMappedRandomAmplitudes() && OneAndTwoQubitGatesTestMapped() && TestMappedMeasurementsWithOneAndTwoQubitGatesCircuits();
}

bool MPSSimulatorTests()
{
	std::cout << "\nMPS Simulator Tests" << std::endl;

	/*
	QC::Gates::PauliXGate xGate;
	QC::Gates::PauliYGate yGate;
	QC::Gates::HadamardGate hGate;
	QC::Gates::CNOTGate cnotGate;
	QC::Gates::SwapGate swapGate;
	QC::Gates::ControlledPhaseShiftGate phaseShiftGate(0.5);

	int nrQubits = 2;
	int nrStates = 1 << nrQubits;

	QC::TensorNetworks::MPSSimulator mps(nrQubits);
	QC::QubitRegister reg(nrQubits);

	mps.ApplyGate(xGate, 0);
	reg.ApplyGate(xGate, 0);

	mps.ApplyGate(hGate, 1);
	reg.ApplyGate(hGate, 1);

	mps.ApplyGate(phaseShiftGate, 0, 1);
	reg.ApplyGate(phaseShiftGate, 0, 1);

	mps.ApplyGate(hGate, 1);
	reg.ApplyGate(hGate, 1);

	mps.ApplyGate(cnotGate, 0, 1);
	reg.ApplyGate(cnotGate, 0, 1);

	for (int i = 0; i < nrStates; ++i)
	{
		const auto regAmpl = reg.getBasisStateAmplitude(i);
		const auto mpsAmpl = mps.getBasisStateAmplitude(i);
		std::cout << "State " << i << " reg: " << regAmpl << " mps: " << mpsAmpl << std::endl;
	}

	mps.print();
	*/

	return StateSimulationTest();
}