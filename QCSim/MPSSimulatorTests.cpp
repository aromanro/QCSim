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
	mps.ApplyGate(cygate, 1, 0);
	
	std::cout << "After applying h and cy:" << std::endl;
	mps.print();

	return true;
}

bool MPSSimulatorTests()
{
	std::cout << "\nMPS Simulator Tests" << std::endl;

	StateSimulationTest();

	return true;
}