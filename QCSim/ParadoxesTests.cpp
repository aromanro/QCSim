#include <iostream>
#include <iterator>
#include <map>

#include "Tests.h"

#include "QuantumEraser.h"
#include "GeneralElitzurVaidmanBomb.h"
#include "HardyParadox.h"



bool QuantumEraserTests()
{
	std::cout << "\nTesting the quantum eraser..." << std::endl;
	std::cout << "\nWithout erasing:" << std::endl;

	Paradoxes::QuantumEraser quantumEraser;

	const int nrMeasurements = 100000;

	auto measurements = quantumEraser.ExecuteWithMultipleMeasurements(nrMeasurements);
	for (auto m : measurements)
		std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

	double expected = nrMeasurements / 4.;
	const double allowedError = expected * 0.05;

	if (!approxEqual(measurements[0], expected, allowedError) || !approxEqual(measurements[1], expected, allowedError)
		|| !approxEqual(measurements[2], expected, allowedError) || !approxEqual(measurements[3], expected, allowedError))
		return false;

	std::cout << "\nWith erasing:" << std::endl;

	quantumEraser.setEraser();

	measurements = quantumEraser.ExecuteWithMultipleMeasurements(nrMeasurements);
	for (const auto& m : measurements)
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

	const size_t maxStages = 4;
	Paradoxes::GeneralElitzurVaidmanBomb ElitzurVaidmanBomb(maxStages); // 4 stages max, 5 qubits

	const int nrMeasurements = 10000;

	const double incr = 0.05 * M_PI;
	for (size_t stages = 1; stages <= maxStages; ++stages)
	{
		std::cout << "\nWith " << stages << " stage" << (stages == 1 ? ":" : "s:") << std::endl;

		ElitzurVaidmanBomb.setStages(stages);

		for (double theta = incr; theta < 0.9 * M_PI; theta += incr)
		{
			ElitzurVaidmanBomb.setTheta(theta);

			auto measurements = ElitzurVaidmanBomb.ExecuteWithMultipleMeasurements(nrMeasurements);
			
			const double result = static_cast<double>(measurements[0]) / (nrMeasurements - static_cast<double>(measurements[1]));
			const double expected = ElitzurVaidmanBomb.TheoreticalEfficiency();
			std::cout << "Theta: " << theta << " 'Experimental' efficiency: " << result << " Theoretical efficiency: " << expected << std::endl;

			if (!approxEqual(result, expected, 0.08))
				return false;
		}
	}

	return true;
}


bool HardyParadoxTests()
{
	std::cout << "\nTesting Hardy's Paradox..." << std::endl;

	const int nrMeasurements = 100000;

	double gammamax = 0;
	double theta0max = 0;
	double theta1max = 0;

	Paradoxes::HardyParadox HardysParadox;

	const double incr = 0.025 * M_PI;
	size_t cnt1 = 0;
	size_t cnt2 = 0;
	for (double theta0 = incr; theta0 < M_PI - incr; theta0 += incr)
	{
		HardysParadox.setTheta0(theta0);
		
		for (double theta1 = incr; theta1 < M_PI - incr; theta1 += incr)
		{
			HardysParadox.setTheta1(theta1);
			auto measurements = HardysParadox.ExecuteWithMultipleMeasurements(nrMeasurements);

			const double p0 = static_cast<double>(measurements[0]) / nrMeasurements;
			const double p1 = static_cast<double>(measurements[1]) / nrMeasurements;
			const double p2 = static_cast<double>(measurements[2]) / nrMeasurements;
			const double p3 = static_cast<double>(measurements[3]) / nrMeasurements;
			const double result = p0 / (p0 + p1 + p2 + p3);
			const double expected = HardysParadox.TheoreticalGamma();

			if (result > gammamax)
			{
				gammamax = result;
				theta0max = theta0;
				theta1max = theta1;
			}

			if (cnt1 % 4 == 0 && cnt2 % 4 == 0) // don't display so many results
				std::cout << "Theta0: " << theta0 / M_PI << " pi, theta1: " << theta1 / M_PI << " pi, 'Experimental' non-locality probability: " << result << " Theoretical non-locality probability: " << expected << std::endl;
			++cnt2;

			if (!approxEqual(result, expected, 0.08))
				return false;
		}
		++cnt1;
	}

	std::cout << "\nTheoretical max is about 0.09 at theta0=theta1=0.575 pi, the 'measured' one found is " << gammamax << " at theta0: " << theta0max / M_PI << " pi, and theta1: " << theta1max / M_PI << " pi" << std::endl;

	return true;
}


bool ParadoxesTests()
{
	std::cout << "\nTesting simulations of quantum paradoxes..." << std::endl;

	return QuantumEraserTests() && ElitzurVaidmanBombTests() && HardyParadoxTests();
}