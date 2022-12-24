#include "Tests.h"

#include "QuantumEraser.h"
#include "GeneralElitzurVaidmanBomb.h"

#include <iostream>
#include <map>


bool QuantumEraserTests()
{
	std::cout << "\nTesting the quantum eraser..." << std::endl;
	std::cout << "\nWithout erasing:" << std::endl;

	Paradoxes::QuantumEraser quantumEraser;

	std::map<int, int> measurements;
	const int nrMeasurements = 100000;

	for (int i = 0; i < nrMeasurements; ++i)
	{
		const unsigned int state = quantumEraser.Execute();
		++measurements[state];
	}

	for (auto m : measurements)
		std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

	int expected = nrMeasurements / 4;
	const double allowedError = expected * 0.05;
	if (!approxEqual(measurements[0], expected, allowedError) || !approxEqual(measurements[1], expected, allowedError)
		|| !approxEqual(measurements[2], expected, allowedError) || !approxEqual(measurements[3], expected, allowedError))
		return false;

	measurements.clear();
	
	std::cout << "\nWith erasing:" << std::endl;

	quantumEraser.setEraser();

	for (int i = 0; i < nrMeasurements; ++i)
	{
		const unsigned int state = quantumEraser.Execute();
		++measurements[state];
	}

	for (auto m : measurements)
		std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

	expected = nrMeasurements / 2;
	if (!approxEqual(measurements[0], expected, allowedError) || !approxEqual(measurements[1], 0)
		|| !approxEqual(measurements[2], 0) || !approxEqual(measurements[3], expected, allowedError))
		return false;

	return true;
}

bool ElitzurVaidmanBombTests()
{
	std::cout << "\nTesting the Elitzur-Vaidman bomb tester/interaction free measurement..." << std::endl;

	const unsigned int maxStages = 4;
	Paradoxes::GeneralElitzurVaidmanBomb ElitzurVaidmanBomb(maxStages); // 4 stages max, 5 qubits

	std::map<int, int> measurements;
	const int nrMeasurements = 10000;

	const double incr = 0.05 * M_PI;
	for (unsigned int stages = 1; stages <= maxStages; ++stages)
	{
		std::cout << "\nWith " << stages << " stage" << (stages == 1 ? ":" : "s:") << std::endl;

		ElitzurVaidmanBomb.setStages(stages);

		for (double theta = incr; theta < 0.9 * M_PI; theta += incr)
		{
			ElitzurVaidmanBomb.setTheta(theta);

			for (int i = 0; i < nrMeasurements; ++i)
			{
				const unsigned int state = ElitzurVaidmanBomb.Execute();

				++measurements[state];
			}

			const double result = static_cast<double>(measurements[0]) / nrMeasurements / (1. - static_cast<double>(measurements[1]) / nrMeasurements);
			const double expected = ElitzurVaidmanBomb.TheoreticalEfficiency();
			std::cout << "Theta: " << theta << " 'Experimental' efficiency: " << result << " Theoretical efficiency: " << expected << std::endl;

			if (!approxEqual(result, expected, 0.05))
				return false;

			measurements.clear();
		}
	}

	return true;
}


bool ParadoxesTests()
{
	std::cout << "\nTesting simulations of quantum paradoxes..." << std::endl;

	bool res = QuantumEraserTests();
	if (res) res = ElitzurVaidmanBombTests();

	return res;
}