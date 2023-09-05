#include "Tests.h"

#include "IsingModel.h"


bool IsingModelTest()
{
	std::cout << "\nTesting QAOA on the Ising Model" << std::endl;

	QC::QubitRegister reg(3);
	Models::QAOAIsingSubalgorithm ising;


	ising.Set({ {0, 1, 1.5}, { 1, 2, 2 }, { 0, 2, -3} }, { 1.25, 3.25, -2.5 });

	for (unsigned int state = 0; state < reg.getNrBasisStates(); ++state)
		std::cout << "State: " << state << " Energy: " << ising.Energy(state) << std::endl;
	std::cout << std::endl;

	// set a fixed start instead of picking a random one
	ising.SetGamma1Start(4);
	ising.SetBeta1Start(3);
	ising.SetGamma2Start(0.3);
	ising.SetBeta2Start(1);

	ising.SetP(2);
	ising.SetMixing(true);

	unsigned int state = ising.Execute(reg);

	unsigned int realMinState = ising.GetMinEnergyState();
	if (state != realMinState)
	{
		std::cout << "Error: Min energy state found: " << state << " is not the same as the real one: " << realMinState << std::endl;
		std::cout << "Due of the stochastic nature, this is expected to fail sometimes!" << std::endl;
		return false;
	}

	std::cout << "Min energy state: " << state << " Energy: " << ising.Energy(state) << std::endl;

	return true;
}

bool MaxCutTest()
{
	std::cout << "\nTesting QAOA on the max cut problem (for min cut just change sign)" << std::endl;

	QC::QubitRegister reg(6);
	Models::QAOAIsingSubalgorithm maxcut;

	maxcut.Set({ {0, 1, -1}, { 0, 2, -1 }, { 0, 5, -1}, { 1, 2, -1 }, { 1, 3, -1 }, { 2, 3, -1 }, { 2, 4, -1 }, { 3, 5, -1 } }, { 0, 0, 0, 0, 0, 0 });

	/*
	for (unsigned int state = 0; state < reg.getNrBasisStates(); ++state)
		std::cout << "State: " << state << " Energy: " << maxcut.Energy(state) << std::endl;
	std::cout << std::endl;
	*/

	auto minStates = maxcut.GetMinEnergyStates();

	/*
	std::cout << "Min energy states: " << std::endl;
	for (auto state : states)
		std::cout << state << " Energy: " << maxcut.Energy(state) << std::endl;
	*/

	// set a fixed start instead of picking a random one
	maxcut.SetGamma1Start(6);
	maxcut.SetBeta1Start(0.7);
	maxcut.SetGamma2Start(4.2);
	maxcut.SetBeta2Start(0.5);

	maxcut.SetDeltaE(0.0005);

	maxcut.SetP(2);
	//maxcut.SetMixing(true);

	unsigned int state = maxcut.Execute(reg);


	if (minStates.find(state) == minStates.end())
	{
		std::cout << "Error: Min energy state found: " << state << " is not among the real ones" << std::endl;
		std::cout << "Due of the stochastic nature, this is expected to fail sometimes!" << std::endl;
		return false;
	}

	std::cout << "Min energy state: " << state << " Energy: " << maxcut.Energy(state) << std::endl;

	return true;
}


bool IsingTests()
{
	return IsingModelTest() && MaxCutTest();
}