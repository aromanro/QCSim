#include "Tests.h"

#include "QubitRegister.h"
#include "PauliStringSimulation.h"
#include "Schrodinger.h"

#include <iostream>
#include <map>

void PrintState(const Eigen::VectorXcd& regVals, unsigned int nrBasisStates)
{
	for (unsigned int i = 0; i < nrBasisStates; ++i)
		std::cout << regVals[i] << std::endl;
}

void PrintStates(const Eigen::VectorXcd& regVals, const Eigen::VectorXcd& regValsEx, unsigned int nrBasisStates)
{
	std::cout << "Exact:" << std::endl;
	PrintState(regValsEx, nrBasisStates);

	std::cout << "*************************************" << std::endl;

	std::cout << "Simulation:" << std::endl;
	PrintState(regVals, nrBasisStates);
}

bool PauliSimultationTests()
{
	std::cout << "\nTesting Pauli strings hamiltonian simulations..." << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_num_terms(1, 5);
	std::uniform_int_distribution<> dist_op(0, 2);
	std::uniform_real_distribution<double> dist_coeff(0.1, 10.);

	const double simTime = 1. / (2. * M_PI);
	const unsigned int nrSteps = 100;

	for (int i = 0; i < 10; ++i)
	{
		std::cout << "Generating and trying Hamiltonian... ";

		QuantumSimulation::PauliDecomposedHamiltonianSimulation sim(3, simTime, nrSteps);
		sim.setToBasisState(0);

		QuantumSimulation::PauliStringSimulation term;

		const int numTerms = dist_num_terms(gen);

		std::cout << numTerms << " term" << (numTerms > 1 ? "s" : ".") << "... ";

		for (int t = 0; t < numTerms; ++t)
		{
			for (int q = 0; q < 3; ++q)
			{
				const QuantumSimulation::PauliStringSimulation<>::PauliOp trm = (QuantumSimulation::PauliStringSimulation<>::PauliOp)dist_op(gen);
				term.setOperatorForQubit(q, trm);
			}

			sim.AddTerm(-dist_coeff(gen), term);
		}

		Eigen::MatrixXcd evOp = sim.getEvolutionOperator(simTime);
		Eigen::VectorXcd regValsEx = sim.getRegisterStorage();
		regValsEx = evOp * regValsEx;

		sim.AdjustPhaseAndNormalize(regValsEx);

		sim.setToBasisState(0);

		sim.Execute();

		sim.AdjustPhaseAndNormalize();
		Eigen::VectorXcd regVals = sim.getRegisterStorage();

		const double fidelity = sim.stateFidelity(regValsEx);
		std::cout << " Simulation result with fidelity: " << fidelity;

		if (fidelity < 0.999)
		{
			std::cout << "\nFidelity too low, failure!" << std::endl;

			PrintStates(regVals, regValsEx, sim.getNrBasisStates());

			return false;
		}

		for (unsigned int i = 0; i < sim.getNrBasisStates(); ++i)
		{
			if (!approxEqual(regVals(i), regValsEx(i), 0.5)) // in some circumstances some values can differ quite a bit but the fidelity is still high
			{
				std::cout << "\nValue from simulation does not match the value from the 'exact' computation!" << std::endl;

				PrintStates(regVals, regValsEx, sim.getNrBasisStates());

				return false;
			}
		}

		std::cout << " OK!" << std::endl;
	}

	return true;
}

bool SchrodingerSimulationTests()
{
	std::cout << "\nTesting Schrodinger simulations, this is going to take a while..." << std::endl;

	const double time = 0.2;
	const double dx = 0.1;
	const int nrSteps = 50;

	QuantumSimulation::SchrodingerSimulation schrSim(9, time, dx, nrSteps);
	QuantumSimulation::SchrodingerSimulation finDifSim(9, time, dx, nrSteps);

	const unsigned int nrStates = schrSim.getNrBasisStates();
	const double len = dx * (nrStates - 1);

	const unsigned int halfPotWidth = 25;
	
	double k = nrStates * dx / time;
	k *= 4;

	const double potential = 60;

	const unsigned int startPos = nrStates / 4;
	const double sigma = 25;

	std::cout << "First, with a potential barrier..." << std::endl;

	schrSim.setConstantPotentialInTheMiddle(potential, halfPotWidth);
	finDifSim.setConstantPotentialInTheMiddle(potential, halfPotWidth);

	schrSim.setGaussian(startPos, sigma, k);
	finDifSim.setGaussian(startPos, sigma, k);

	schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_start.dat");
	finDifSim.getRegister().writeToFile("c:\\temp\\findif_start.dat");

	for (unsigned int i = 0; i < 15; ++i)
	{
		schrSim.Execute();
		finDifSim.solveWithFiniteDifferences();

		schrSim.AdjustPhaseAndNormalize();
		finDifSim.AdjustPhaseAndNormalize();

		const double fidelity = schrSim.stateFidelity(finDifSim.getRegisterStorage());

		std::cout << "Part " << i << " Fidelity: " << fidelity << std::endl;

		schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_end.dat");
		finDifSim.getRegister().writeToFile("c:\\temp\\findif_end.dat");

		if (fidelity < 0.5) 
		{
			std::cout << "\nFidelity too low, failure!" << std::endl;
			return false;
		}
	}

	std::cout << "Second, with a potential well..." << std::endl;

	schrSim.setConstantPotentialInTheMiddle(-0.35 * potential, halfPotWidth);
	finDifSim.setConstantPotentialInTheMiddle(-0.35 * potential, halfPotWidth);

	schrSim.setGaussian(startPos, sigma, k);
	finDifSim.setGaussian(startPos, sigma, k);

	for (unsigned int i = 0; i < 15; ++i)
	{
		schrSim.Execute();
		finDifSim.solveWithFiniteDifferences();

		schrSim.AdjustPhaseAndNormalize();
		finDifSim.AdjustPhaseAndNormalize();

		const double fidelity = schrSim.stateFidelity(finDifSim.getRegisterStorage());

		std::cout << "Part " << i << " Fidelity: " << fidelity << std::endl;

		schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_well_end.dat");
		finDifSim.getRegister().writeToFile("c:\\temp\\findif_well_end.dat");

		if (fidelity < 0.5)
		{
			std::cout << "\nFidelity too low, failure!" << std::endl;
			return false;
		}
	}
	
	return true;
}

bool SimulationTests()
{
	std::cout << "\nTesting quantum simulations..." << std::endl;

	bool res = PauliSimultationTests();
	if (res) res = SchrodingerSimulationTests();

	return res;
}
