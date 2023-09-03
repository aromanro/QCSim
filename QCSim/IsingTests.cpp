#include "Tests.h"

#include "IsingModel.h"


bool IsingTests()
{
	std::cout << "\nTesting QAOA on the Ising Model" << std::endl;

	QC::QubitRegister reg(3);
	Models::IsingSubalgorithm ising;


	ising.Set({ {0, 1, 1.5}, { 1, 2, 2 }, { 0, 2, -3} }, { 1.25, 3.25, -2.5 });

	for (unsigned int state = 0; state < reg.getNrBasisStates(); ++state)
		std::cout << "State: " << state << " Energy: " << ising.Energy(state) << std::endl;
	std::cout << std::endl;

	double Eold = std::numeric_limits<double>::max();
	double Emin = Eold;
	double E = 0; // whatever
	double deltaE = 0.001;

	double gamma = 1;
	double beta = 1;

	double optGamma = gamma;
	double optBeta = beta;

	bool first = true;
	int i = 0;

	int p = 1;
	do {
		Eold = E;

		if (first) first = false;
		else std::tie(gamma, beta) = ising.GradientDescentStep(reg, gamma, beta, p);

		ising.Exec(reg, gamma, beta, p);
		E = ising.EnergyExpectationValue(reg);
		
		if (E < Emin)
		{
			Emin = E;
			optGamma = gamma;
			optBeta = beta;
		}

		if (i % 10 == 0)
			std::cout << "E: " << E << " gamma: " << gamma << " beta: " << beta << std::endl;

		++i;
	} while (abs(E-Eold) > deltaE);

	ising.Exec(reg, gamma, beta, p);
	auto res = reg.RepeatedMeasure(100000);

	Emin = std::numeric_limits<double>::max();
	unsigned int state = 0;
	for (const auto& r : res)
	{
		E = ising.Energy(r.first) * static_cast<double>(r.second) / 100000.;
		if (E < Emin)
		{
			Emin = E;
			state = r.first;
		}
	}

	unsigned int realMinState = ising.GetMinEnergyState();
	if (state != realMinState)
	{
		std::cout << "Error: Min energy state found: " << state << " is not the same as the real one: " << realMinState << std::endl;
		return false;
	}

	std::cout << "Min energy state: " << state << " Energy: " << ising.Energy(state) << std::endl;

	return true;
}