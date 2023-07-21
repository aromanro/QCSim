#include "Tests.h"

#include "PhaseEstimation.h"

#include <iostream>
#include <vector>

bool PhaseEstimationTests()
{
	std::cout << "\nTesting phase estimation algorithm..." << std::endl;

	const int nrMeasurements = 1000;
	QC::Gates::PauliXGate x;

	std::vector<double> phases{ 1. / 4., 1. / 3., 1. / 2., 2. / 3., 3. / 4. };

	for (int p = 0; p < phases.size(); ++p)
	{
		std::vector<int> measurements(8, 0);
		const double realPhase = phases[p];
		QC::Gates::PhaseShiftGate phaseShift(realPhase * 2 * M_PI);

		for (int i = 0; i < nrMeasurements; ++i)
		{
			QC::QubitRegister reg(4);
			reg.setToBasisState(0);
			reg.ApplyGate(x, 3); // flip the last qubit

			QC::SubAlgo::PhaseEstimation<> phaseEstimation(phaseShift.getRawOperatorMatrix(), 4, 3);

			const unsigned int res = phaseEstimation.Execute(reg);
			++measurements[res];
		}

		int mx = -1;
		int cnt = 0;
		for (int i = 0; i < 8; ++i)
		{
			if (measurements[i] > cnt)
			{
				mx = i;
				cnt = measurements[i];
			}
		}

		const double estimatedPhase = static_cast<double>(mx) / (1ULL << 3);

		std::cout << "Real phase: " << realPhase << ", Estimated: " << estimatedPhase << std::endl;

		if (!approxEqual(realPhase, estimatedPhase, 0.05)) return false;
	}

	return true;
}
