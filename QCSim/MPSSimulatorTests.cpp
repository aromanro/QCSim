#include <iostream>
#include <iterator>
#include <limits>
#include <memory>


#include "Tests.h"

#include "QubitRegister.h"
#include "MPSSimulator.h"

#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>
#include <future>


#define NR_QUBITS_LIMIT 9

void FillOneQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
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
					if (!approxEqual(regState[s], mpsState[s], 1E-3))
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
		
		circuit.push_back(std::move(appliedGate));
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
						for (const auto& [key2, value2] : measurementsRegMap)
						{
							for (const auto& b : key2)
								std::cout << b << " ";
							std::cout << " : " << static_cast<double>(value2) / nrMeasurements << std::endl;
						}

						std::cout << "MPS measurements:\n";
						for (const auto& [key2, value2] : measurementsMPSMap)
						{
							for (const auto& b : key2)
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
					if (!approxEqual(regState[s], mpsState[s], 1E-3))
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

bool NumericalRankStabilityTestMPS()
{
	std::cout << "\nMPS simulator numerical rank stability test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);
	FillTwoQubitGates(gates);

	// This fixed seed contains a mapped circuit that creates roundoff-only
	// Schmidt values and subsequently routes gates across them. Without a
	// numerical rank floor, the Vidal pseudoinverse amplifies that SVD noise
	// into percent-level state-vector errors.
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
				if (twoQubitsGate && stabilityBool(stabilityGenerator)) std::swap(qubit1, qubit2);

				mps.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);

				const auto& reference = reg.getRegisterStorage();
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
								std::vector<bool> measurementsMPSAllVec(nrQubits);
								for (const auto& [qubit, val] : measurementsMPSAll)
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

	for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT; ++nrQubits)
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

					if (!approxEqual(regAmpl, mpsAmpl, 1E-3))
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


bool CheckQubitsProbability()
{
	std::cout << "\nMPS swapping/mapped simulator qubits probability test for both one and two qubit gates" << std::endl;
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
			std::cout << "\n\n\nTest no: " << t << " for " << nrQubits << " qubits" << std::endl << std::endl << std::endl;
#endif

			QC::TensorNetworks::MPSSimulator mps(nrQubits);
			QC::QubitRegister reg(nrQubits);

			const auto circ = GenerateRandomCircuitWithGatesNoAdjacent(gates, 50, 150, nrQubits);
			for (const auto& gate : circ)
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
					std::cout << "Qubit " << q << " probability test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;
					std::cout << "Probability for the qubit in statevector: " << regProb << " vs mps: " << mpsProb << std::endl;
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
		if (!(OneAndTwoQubitGatesTestMappedRandomAmplitudes() && OneAndTwoQubitGatesTestMapped() && TestMappedMeasurementsWithOneAndTwoQubitGatesCircuits()))
			return false;
	}
	*/

	return /*OneQubitGatesTest() && OneAndTwoQubitGatesTest() &&
		TestMeasurementsWithOneQubitGatesCircuits() && TestMeasurementsWithOneAndTwoQubitGatesCircuits() &&*/
		OneAndTwoQubitGatesTestMappedRandomAmplitudes() && OneAndTwoQubitGatesTestMapped() && TestMappedMeasurementsWithOneAndTwoQubitGatesCircuits() && CheckQubitsProbability();
}

bool matchRandomExpectationValues(QC::QubitRegister<>& reg, QC::TensorNetworks::MPSSimulator& mps)
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
			std::cout << std::endl << "Expectation values are not equal for MPS and statevector simulator for " << nrQubits << " qubits, values: " << exp1 << ", " << exp2 << std::endl;
	
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

				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				circuit.emplace_back(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
			}

			const int chi = 1 + (t % 4); // trim down to a bond dimension between 1 and 4

			// reference: apply the circuit without any limit (exact MPS), then lower the limit and
			// truncate with the no-op two qubit gate exactly the bonds Trim would touch, that is
			// only the ones that exceed the (newly lowered) limit. Trimming a bond that is already
			// within the limit is a state preserving re-gauging that would still change the result
			// of subsequent truncations, so the reference must skip those bonds just like Trim does.
			QC::TensorNetworks::MPSSimulatorImpl mpsRef(nrQubits);
			for (const auto& g : circuit) mpsRef.ApplyGate(g);
			mpsRef.setLimitBondDimension(chi);
			const auto refBondDims = mpsRef.getBondDimensions();
			for (int q = 0; q < nrQubits - 1; ++q)
				if (refBondDims[q] > chi)
					mpsRef.ApplyGate(identityTwoQubitGate, q, q + 1);

			// the one under test: apply the same circuit without any limit, then lower the limit and Trim
			QC::TensorNetworks::MPSSimulatorImpl mpsTrim(nrQubits);
			for (const auto& g : circuit) mpsTrim.ApplyGate(g);
			mpsTrim.setLimitBondDimension(chi);
			mpsTrim.Trim();

			// after trimming, every bond dimension must be within the (newly lowered) limit
			for (const auto bd : mpsTrim.getBondDimensions())
				if (bd > chi)
				{
					std::cout << "Trim did not reduce a bond dimension to the limit (" << bd << " > " << chi << ") for " << nrQubits << " qubits" << std::endl;
					return false;
				}

			// Trim must produce the same state as the equivalent no-op two qubit gate truncations
			const auto refState = mpsRef.getRegisterStorage();
			const auto trimState = mpsTrim.getRegisterStorage();

			for (int s = 0; s < refState.size(); ++s)
				if (!approxEqual(refState[s], trimState[s], 1E-3))
				{
					std::cout << "Trim state " << s << " differs from the reference truncation for " << nrQubits << " qubits" << std::endl;
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
		if (!mps.setTruncationMode(TruncationMode::RelativeToMax) || mps.getTruncationMode() != TruncationMode::RelativeToMax)
		{
			std::cout << "Truncation mode round trip failed for RelativeToMax" << std::endl;
			return false;
		}

		if (!mps.setTruncationMode(TruncationMode::DiscardedWeight) || mps.getTruncationMode() != TruncationMode::DiscardedWeight)
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
	seed << std::complex<double>(0.3, 0.7), std::complex<double>(-0.5, 0.2), std::complex<double>(0.4, -0.6), std::complex<double>(0.1, 0.9),
		std::complex<double>(-0.8, 0.1), std::complex<double>(0.6, 0.5), std::complex<double>(0.2, 0.3), std::complex<double>(-0.4, 0.7),
		std::complex<double>(0.5, -0.3), std::complex<double>(0.9, -0.1), std::complex<double>(-0.6, 0.4), std::complex<double>(0.2, 0.2),
		std::complex<double>(-0.2, 0.6), std::complex<double>(0.3, -0.8), std::complex<double>(0.7, 0.1), std::complex<double>(-0.5, -0.3);
	const Eigen::MatrixXcd joinGate = Eigen::HouseholderQR<Eigen::MatrixXcd>(seed).householderQ();

	auto buildCircuit = [&](QC::TensorNetworks::MPSSimulatorImpl& mps)
	{
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
		std::cout << "Test setup did not produce a rich enough spectrum (" << spectrum.size() << " singular values)" << std::endl;
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
		while (i >= 0 && spectrum[i] < cutoffValue) --i;
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
			if ((discarded + sq) / total >= threshold) break;

			discarded += sq;
			keep = i;
		}
		expectedDiscardedWeight = keep;
	}

	if (expectedRelative == expectedDiscardedWeight)
	{
		std::cout << "Test setup does not discriminate between the two truncation modes at threshold "
			<< threshold << " (both would keep " << expectedRelative << ")" << std::endl;
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
			std::cout << "RelativeToMax kept " << bondDims[bondIndex] << " singular values, expected " << expectedRelative << std::endl;
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
			std::cout << "DiscardedWeight kept " << bondDims[bondIndex] << " singular values, expected " << expectedDiscardedWeight << std::endl;
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

	auto checkSimulator = [&state, digits, nrQubits](QC::TensorNetworks::MPSSimulatorInterface& simulator, const char* simulatorName)
		{
			simulator.setToBasisState(state);

			for (size_t q = 0; q < nrQubits; ++q)
			{
				const double expectedOneProbability = state[q] ? 1. : 0.;
				if (!approxEqual(simulator.GetProbability(static_cast<Eigen::Index>(q), false), expectedOneProbability, 1E-12))
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

			const std::vector<bool> shortState{ true, false, true };
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
	for (const size_t scalarQubitCount : { digits, digits + 1 })
	{
		auto checkScalarSimulator = [allBits, digits, scalarQubitCount](QC::TensorNetworks::MPSSimulatorInterface& simulator, const char* simulatorName)
			{
				simulator.setToBasisState(allBits);

				for (size_t q = 0; q < scalarQubitCount; ++q)
				{
					const double expectedOneProbability = q < digits ? 1. : 0.;
					if (!approxEqual(simulator.GetProbability(static_cast<Eigen::Index>(q), false), expectedOneProbability, 1E-12))
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
	mappedValidationSimulator.setToBasisState(std::vector<bool>{ true, false, true });
	mappedValidationSimulator.MoveAtBeginningOfChain({ 2 });
	const auto stateBeforeInvalidCall = std::dynamic_pointer_cast<QC::TensorNetworks::MPSSimulatorState>(mappedValidationSimulator.getState());
	const std::vector<Eigen::Index> identityMap{ 0, 1, 2 };
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
	catch (const std::invalid_argument&)
	{
		oversizedStateRejected = true;
	}
	catch (...)
	{
	}

	const auto stateAfterInvalidCall = std::dynamic_pointer_cast<QC::TensorNetworks::MPSSimulatorState>(mappedValidationSimulator.getState());
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

	return WideBasisInitializationTestMPS() && StateSimulationTest() && NumericalRankStabilityTestMPS() && checkExpectationValuesMPS() && TrimTestMPS() && TruncationModeTestMPS();
}


