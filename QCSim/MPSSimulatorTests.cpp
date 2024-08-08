#include <iostream>
#include <iterator>


#include "Tests.h"

#include "MPSSimulator.h"

bool StateSimulationTest()
{
	QC::TensorNetworks::MPSSimulator mps(2);

	mps.print();

	QC::Gates::PauliXGate xgate;
	mps.ApplyGate(xgate, 0);

	std::cout << "After applying X gate:" << std::endl;
	mps.print();

	QC::Gates::HadamardGate hgate;
	mps.ApplyGate(hgate, 0);

	std::cout << "After applying H gate:" << std::endl;
	mps.print();

	QC::Gates::PauliYGate ygate;
	mps.ApplyGate(ygate, 0);

	std::cout << "After applying Y gate:" << std::endl;
	mps.print();

	QC::Gates::CNOTGate cnotgate;
	mps.ApplyGate(cnotgate, 1, 0);

	std::cout << "After applying cnot gate:" << std::endl;
	mps.print();

	QC::Gates::ControlledYGate cygate;
	mps.ApplyGate(hgate, 0);
	mps.ApplyGate(cygate, 0, 1);
	
	std::cout << "After applying h and cy:" << std::endl;
	mps.print();

	std::cout << "Probability of getting 0 for qubit 0: " << mps.GetProbability0(0) << std::endl;

	std::cout << "Measured 0: " << !mps.Measure(0) << std::endl;

	

	// it's the same circuit as above, but shifted, on a larger simulator
	// should give 50% probability of measuring 1
	for (int shift = 0; shift < 3; ++shift)
	{
		int cnt = 0;
		for (int i = 0; i < 100; ++i)
		{
			QC::TensorNetworks::MPSSimulator mpsl(4);
			mpsl.ApplyGate(xgate, shift);
			mpsl.ApplyGate(hgate, shift);
			mpsl.ApplyGate(ygate, shift);
			mpsl.ApplyGate(cnotgate, shift + 1, shift);
			mpsl.ApplyGate(hgate, shift);
			mpsl.ApplyGate(cygate, shift, shift + 1);

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