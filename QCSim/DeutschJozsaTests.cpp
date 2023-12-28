#include <iostream>
#include <iterator>
#include <map>

#include "QubitRegister.h"
#include "DeutschJozsa.h"

bool DeutschWithGatesTests()
{
	std::cout << "Deutsch's algorithm using gates for oracle..." << std::endl;
	DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<> deutschWithGates;
	for (int i = 0; i < 30; ++i)
	{
		if (i < 10)
			deutschWithGates.setFunction(DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<>::FunctionType::constantZero);
		else if (i < 20)
			deutschWithGates.setFunction(DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<>::FunctionType::constantOne);
		else
			deutschWithGates.setFunction(DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<>::FunctionType::balanced);

		unsigned int state = deutschWithGates.Execute();

		if (i < 20)
		{
			// constant
			if (!deutschWithGates.WasConstantResult(state))
			{
				std::cout << "Expected constant result, got balanced" << std::endl;
				return false;
			}
			else std::cout << "Constant " << (i < 10 ? "zero" : "one") << " function set, got constant, ok" << std::endl;
		}
		else
		{
			// balanced
			if (deutschWithGates.WasConstantResult(state))
			{
				std::cout << "Expected balanced result, got constant" << std::endl;
				return false;
			}
			else std::cout << "Balanced function set, got balanced, ok" << std::endl;
		}
	}

	return true;
}

bool DeutschTests()
{
	// first the simplest case, Deutsch's algorithm (N=3)
	std::cout << "First, Deutsch's algorithm..." << std::endl;
	DeutschJozsa::DeutschJozsaAlgorithm<> deutsch;
	for (int i = 0; i < 15; ++i)
	{
		if (i < 5)
			deutsch.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantZero);
		else if (i < 10)
			deutsch.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantOne);
		else
			deutsch.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::balanced);

		unsigned int state = deutsch.Execute();

		if (i < 10)
		{
			// constant
			if (!deutsch.WasConstantResult(state))
			{
				std::cout << "Expected constant result, got balanced" << std::endl;
				return false;
			}
			else std::cout << "Constant " << (i < 5 ? "zero" : "one") << " function set, got constant, ok" << std::endl;
		}
		else
		{
			// balanced
			if (deutsch.WasConstantResult(state))
			{
				std::cout << "Expected balanced result, got constant" << std::endl;
				return false;
			}
			else std::cout << "Balanced function set, got balanced, ok" << std::endl;
		}
	}

	return DeutschWithGatesTests();
}

bool DeutschJozsaWithGatesTests()
{
	std::cout << "Now, Deutsch-Jozsa with the oracle made out of gates..." << std::endl;
	DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<> deutschJozsaWithGates(7);
	for (int i = 0; i < 30; ++i)
	{
		if (i < 10)
			deutschJozsaWithGates.setFunction(DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<>::FunctionType::constantZero);
		else if (i < 20)
			deutschJozsaWithGates.setFunction(DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<>::FunctionType::constantOne);
		else
			deutschJozsaWithGates.setFunction(DeutschJozsa::DeutschJozsaAlgorithmWithGatesOracle<>::FunctionType::balanced);

		unsigned int state = deutschJozsaWithGates.Execute();

		if (i < 20)
		{
			// constant
			if (!deutschJozsaWithGates.WasConstantResult(state))
			{
				std::cout << "Expected constant result, got balanced" << std::endl;
				return false;
			}
			else std::cout << "Constant " << (i < 10 ? "zero" : "one") << " function set, got constant, ok" << std::endl;
		}
		else
		{
			// balanced
			if (deutschJozsaWithGates.WasConstantResult(state))
			{
				std::cout << "Expected balanced result, got constant" << std::endl;
				return false;
			}
			else std::cout << "Balanced function set, got balanced, ok" << std::endl;
		}
	}

	return true;
}

bool DeutschJozsaTests()
{
	std::cout << "\nTesting Deutsch-Jozsa..." << std::endl;

	if (!DeutschTests()) return false;

	std::cout << "Now, general Deutsch-Jozsa..." << std::endl;
	DeutschJozsa::DeutschJozsaAlgorithm<> deutschJozsa(9);
	for (int i = 0; i < 30; ++i)
	{
		if (i < 10)
			deutschJozsa.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantZero);
		else if (i < 20)
			deutschJozsa.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantOne);
		else
			deutschJozsa.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::balanced);

		unsigned int state = deutschJozsa.Execute();

		if (i < 20)
		{
			// constant
			if (!deutschJozsa.WasConstantResult(state))
			{
				std::cout << "Expected constant result, got balanced" << std::endl;
				return false;
			}
			else std::cout << "Constant " << (i < 10 ? "zero" : "one") << " function set, got constant, ok" << std::endl;
		}
		else
		{
			// balanced
			if (deutschJozsa.WasConstantResult(state))
			{
				std::cout << "Expected balanced result, got constant" << std::endl;
				return false;
			}
			else std::cout << "Balanced function set, got balanced, ok" << std::endl;
		}
	}

	return DeutschJozsaWithGatesTests();
}
