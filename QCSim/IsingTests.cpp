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

	int p = 1;
	for (int i = 0; abs(E - Eold) > deltaE && i < 100000; ++i) 
	{
		Eold = E;

		// this is more like stochastic gradient descent, so it would benefit from
		// the methods used in the machine learning project (momentum/nesterov/adagrad/rmsprop/adam/whatever) 
		// I won't bother, for the curious the methods are implemented there (also in the python repo, the dft notebook)

		std::tie(gamma, beta) = ising.GradientDescentStep(reg, gamma, beta, p);

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
	}

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
		std::cout << "Due of the stochastic nature, this is expected to fail sometimes!" << std::endl;
		return false;
	}

	std::cout << "Min energy state: " << state << " Energy: " << ising.Energy(state) << std::endl;

	return true;
}