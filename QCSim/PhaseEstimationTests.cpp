#include "Tests.h"

#include "PhaseEstimation.h"

#include <iostream>
#include <vector>

bool PhaseEstimationTests()
{
	std::cout << "\nTesting phase estimation algorithm..." << std::endl;

	const int nrMeasurements = 10000;
	

	std::vector<double> phases{ 1. / 4., 1. / 3., 1. / 2., 2. / 3., 3. / 4. };


	for (int p = 0; p < phases.size(); ++p)
	{
		const double realPhase = phases[p];
		QC::Gates::PhaseShiftGate phaseShift(realPhase * 2 * M_PI);

		QC::SubAlgo::PhaseEstimation<> phaseEstimation(phaseShift.getRawOperatorMatrix(), 4, 3);

		QC::QubitRegister reg(4);

		//QC::Gates::PauliXGate x;
		//reg.setToBasisState(0);
		//reg.ApplyGate(x, 3); // flip the last qubit

		// either the commented above, or this:
		reg.setToQubitState(3);

		const auto measurements = phaseEstimation.ExecuteWithMultipleMeasurements(reg, nrMeasurements);

		int mx = -1;
		int prevMx = -1;
		int cnt = 0;
		int prevCnt = 0;
		for (const auto& m : measurements)
		{
			if (m.second > cnt)
			{
				prevMx = mx;
				prevCnt = cnt;

				mx = m.first;
				cnt = m.second;
			}
			else if (m.second > prevCnt)
			{
				prevMx = m.first;
				prevCnt = m.second;
			}
		}

		const double estimatedPhase = static_cast<double>(mx) / (1ULL << 3);

		std::cout << "Real phase: " << realPhase << ", Estimated: " << estimatedPhase << ", Better estimation: " << phaseEstimation.getPhase(mx, prevMx, static_cast<double>(cnt) / nrMeasurements)  << std::endl;

		if (!approxEqual(realPhase, estimatedPhase, 0.05)) return false;
	}

	return true;
}
