#include "Tests.h"

#include "QuantumEraser.h"

#include <iostream>
#include <map>


bool QuantumEraserTests()
{
	std::cout << "\nTesting the quantum eraser..." << std::endl;
	std::cout << "\nWithout erasing:" << std::endl;

	

	Paradoxes::QuantumEraser quantumEraser;

	std::map<int, int> measurements;
	const int nrMeasurements = 100;

	for (int i = 0; i < nrMeasurements; ++i)
	{
		const unsigned int state = quantumEraser.Execute();
		++measurements[state];
	}

	for (auto m : measurements)
		std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

	int expected = nrMeasurements / 4;
	double allowedError = expected / 4.;
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
	allowedError = expected / 8.;
	if (!approxEqual(measurements[0], expected, allowedError) || !approxEqual(measurements[1], 0)
		|| !approxEqual(measurements[2], 0) || !approxEqual(measurements[3], expected, allowedError))
		return false;

	return true;
}


bool ParadoxesTests()
{
	std::cout << "\nTesting simulations of quantum paradoxes..." << std::endl;

	bool res = QuantumEraserTests();

	return res;
}