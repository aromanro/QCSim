#include <iostream>
#include <iterator>


#include "Tests.h"

#include "MPSSimulator.h"

bool StateSimulationTest()
{
	QC::Gates::PauliXGate xgate;
	QC::Gates::PauliYGate ygate;
	QC::Gates::HadamardGate hgate;
	QC::Gates::CNOTGate cnotgate;
	QC::Gates::ControlledYGate cygate;

	QC::TensorNetworks::MPSSimulator mps(2);

	mps.print();

	mps.ApplyGate(xgate, 0);
	mps.ApplyGate(hgate, 0);
	mps.ApplyGate(cnotgate, 1, 0);
	mps.ApplyGate(ygate, 1);
	mps.ApplyGate(cnotgate, 1, 0);
	mps.ApplyGate(hgate, 1);
	mps.ApplyGate(cygate, 0, 1);

	std::cout << "After applying gates:" << std::endl;

	mps.print();

	Eigen::VectorXcd v = mps.getRegisterStorage();

	std::cout << "State vector: " << v << std::endl;

	std::cout << "Probability of getting 0 for qubit 0: " << mps.GetProbability0(0) << std::endl;

	std::cout << "Measured 0: " << !mps.Measure(0) << std::endl;

	// should give 50% probability of measuring 1
	for (int shift = 0; shift < 4; ++shift)
	{
		int cnt = 0;
		for (int i = 0; i < 100; ++i)
		{
			QC::TensorNetworks::MPSSimulator mpsl(5);
			
			mpsl.ApplyGate(xgate, shift);
			mpsl.ApplyGate(hgate, shift);
			mpsl.ApplyGate(cnotgate, shift + 1, shift);
			mpsl.ApplyGate(ygate, shift + 1);
			mpsl.ApplyGate(cnotgate, shift + 1, shift);
			mpsl.ApplyGate(hgate, shift + 1);
			mpsl.ApplyGate(cygate, shift, shift + 1);
			//std::cout << "Probability of getting 0 for qubit 0: " << mpsl.GetProbability0(shift) << std::endl;
			if (mpsl.Measure(shift))
				++cnt;
		}

		std::cout << "Measured 1 for 100 runs: " << cnt << std::endl;
	}

	return true;
}

bool MPSSimulatorTests()
{
	std::cout << "\nMPS Simulator Tests" << std::endl;

	StateSimulationTest();

	return true;
}