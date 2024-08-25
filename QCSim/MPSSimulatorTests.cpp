#include <iostream>
#include <iterator>
#include <memory>


#include "Tests.h"

#include "QubitRegister.h"
#include "MPSSimulator.h"

#define _USE_MATH_DEFINES
#include <math.h>


void FillOneQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::HyGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PhaseGate<>>());
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

	std::uniform_int_distribution nrGatesDistr(25, 50);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);


	for (int nrQubits = 1; nrQubits < 7; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPSSimulator mps(nrQubits);
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

	std::cout << "Success" << std::endl;

	return true;
}

bool OneAndTwoQubitGatesTest()
{
	std::cout << "\nMPS simulator state test for both one and two qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGates(gates);
	FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(50, 100);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);


	for (int nrQubits = 2; nrQubits < 3/*7*/; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPSSimulator mps(nrQubits);
			QC::QubitRegister reg(nrQubits);

			const int lim = nrGatesDistr(gen);

			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;

				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				//if (twoQubitsGate) std::cout << "Applying two qubit gate on qubits " << qubit1 << " and " << qubit2 << std::endl;

				mps.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);
			}

			// now check the results, they should be the same
			const auto& regState = reg.getRegisterStorage();
			auto mpsState = mps.getRegisterStorage(); // this one is computed, returns value, not reference, not stored elsewhere

			//QC::QubitRegister regNorm(nrQubits);
			//regNorm.setRegisterStorage(mpsState);
			//mpsState = regNorm.getRegisterStorage();

			for (int i = 0; i < regState.size(); ++i)
			{
				if (!approxEqual(regState[i], mpsState[i]))
				{
					std::cout << "State simulation test failed for the MPS simulator for " << nrQubits << " qubits" << std::endl;

					std::cout << "Probability for different state: " << std::norm(regState[i]) << " vs " << std::norm(mpsState[i]) << std::endl;

					std::cout << "Reg state:\n" << regState << std::endl;
					std::cout << "Reg state normalization: " << (regState.adjoint() * regState)(0).real() << std::endl;

					std::cout << "MPS state:\n" << mpsState << std::endl;

					std::cout << "MPS state normalization: " << (mpsState.adjoint() * mpsState)(0).real() << std::endl;
					return false;
				}
			}
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

bool StateSimulationTest()
{
	return OneQubitGatesTest() && OneAndTwoQubitGatesTest();
}

bool MPSSimulatorTests()
{
	std::cout << "\nMPS Simulator Tests" << std::endl;

	return StateSimulationTest();
}