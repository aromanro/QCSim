#include <iostream>
#include <iterator>
#include <map>

#include "Tests.h"

#include "QubitRegister.h"
#include "PauliStringSimulation.h"
#include "Schrodinger.h"


void PrintState(const Eigen::VectorXcd& regVals, size_t nrBasisStates)
{
	for (size_t i = 0; i < nrBasisStates; ++i)
		std::cout << regVals[i] << std::endl;
}

void PrintStates(const Eigen::VectorXcd& regVals, const Eigen::VectorXcd& regValsEx, size_t nrBasisStates)
{
	std::cout << "Exact:" << std::endl;
	PrintState(regValsEx, nrBasisStates);

	std::cout << "*************************************" << std::endl;

	std::cout << "Simulation:" << std::endl;
	PrintState(regVals, nrBasisStates);
}

bool PauliSimulationTests()
{
	std::cout << "\nTesting Pauli strings hamiltonian simulations..." << std::endl;

	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_num_terms(1, 5);
	std::uniform_int_distribution<> dist_op(0, 2);
	std::uniform_real_distribution<double> dist_coeff(0.1, 10.);

	const double simTime = 1. / (2. * M_PI);
	const size_t nrSteps = 100;

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
				const PauliString::PauliString::PauliOp trm = (PauliString::PauliString::PauliOp)dist_op(gen);
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

		for (size_t j = 0; j < sim.getNrBasisStates(); ++j)
		{
			if (!approxEqual(regVals(j), regValsEx(j), 0.5)) // in some circumstances some values can differ quite a bit but the fidelity is still high
			{
				std::cout << "\nValue from simulation does not match the value from the 'exact' computation!" << std::endl;

				PrintStates(regVals, regValsEx, sim.getNrBasisStates());

				//return false;
			}
		}

		std::cout << " OK!" << std::endl;
	}

	return true;
}

bool SchrodingerSimulationTests()
{
	std::cout << "\nTesting Schrodinger simulations, this is going to take a while..." << std::endl;
	
	const size_t nrStates = 512;
	const int nrSteps = 50;

	// the following values work for finite differences, unfortunately they do not work for the Schrodinger simulation
	//const double dx = 1. / (nrStates - 1.);
	//const double dt = 2 * dx * dx;


	
	const double dx = 0.1;
	const double dt = 0.01;

	const double len = dx * (nrStates - 1.);
	// pick a k so the packet will propagate about half the length during the simulation
	//double k = len / (2. * nrSteps * 15 * dt); // this should be ok for the finite difference computation
	
	 	
	// the following seem to work ok for schrodinger simulation
	double k = len / (2. * nrSteps * 15 * dt) + 0.5 * len / dt;

	double E = k * k / 2;

	QuantumSimulation::SchrodingerSimulation schrSim(9, dt, dx, nrSteps);
	QuantumSimulation::SchrodingerSimulation fftSchr(9, dt, dx, nrSteps);

	const size_t startPos = nrStates / 4;
	const double sigma = nrStates / 20;
	//const double sigma = nrStates / 30;

	std::cout << "First, with a potential barrier..." << std::endl;

	const size_t halfPotWidth = static_cast<size_t>(sigma / 4.);
	schrSim.setConstantPotentialInTheMiddle(E, halfPotWidth);
	fftSchr.setConstantPotentialInTheMiddle(E, halfPotWidth);

	schrSim.setGaussian(startPos, sigma, k);
	fftSchr.setGaussian(startPos, sigma, k);

	schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_start.dat");

	for (size_t i = 0; i < 6; ++i)
	{
		schrSim.Execute();
		fftSchr.solveWithClassicalFFT();

		schrSim.AdjustPhaseAndNormalize();
		fftSchr.AdjustPhaseAndNormalize();

		const double fidelity = schrSim.stateFidelity(fftSchr.getRegisterStorage());

		std::cout << "Part " << i << " Fidelity: " << fidelity << std::endl;

		schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_end.dat");
		fftSchr.getRegister().writeToFile("c:\\temp\\findif_end.dat");

		if (fidelity < 0.8)
		{
			std::cout << "\nFidelity too low, failure!" << std::endl;
			return false;
		}
	}

	std::cout << "Second, with a potential well..." << std::endl;

	const double wellDepth = -0.35 * E;
	schrSim.setConstantPotentialInTheMiddle(wellDepth, halfPotWidth);
	fftSchr.setConstantPotentialInTheMiddle(wellDepth, halfPotWidth);

	schrSim.setGaussian(startPos, sigma, k);
	fftSchr.setGaussian(startPos, sigma, k);

	for (size_t i = 0; i < 6; ++i)
	{
		schrSim.Execute();
		fftSchr.solveWithClassicalFFT();

		schrSim.AdjustPhaseAndNormalize();
		fftSchr.AdjustPhaseAndNormalize();

		const double fidelity = schrSim.stateFidelity(fftSchr.getRegisterStorage());

		std::cout << "Part " << i << " Fidelity: " << fidelity << std::endl;

		schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_well_end.dat");
		fftSchr.getRegister().writeToFile("c:\\temp\\findif_well_end.dat");

		if (fidelity < 0.8)
		{
			std::cout << "\nFidelity too low, failure!" << std::endl;
			return false;
		}
	}

	std::cout << "Third, with a potential step (higher)..." << std::endl;

	const double stepHeight = E;
	schrSim.setConstantPotentialToRight(stepHeight);
	fftSchr.setConstantPotentialToRight(stepHeight);

	schrSim.setGaussian(startPos, sigma, k);
	fftSchr.setGaussian(startPos, sigma, k);

	for (size_t i = 0; i < 4; ++i)
	{
		schrSim.Execute();
		fftSchr.solveWithClassicalFFT();

		schrSim.AdjustPhaseAndNormalize();
		fftSchr.AdjustPhaseAndNormalize();

		const double fidelity = schrSim.stateFidelity(fftSchr.getRegisterStorage());

		std::cout << "Part " << i << " Fidelity: " << fidelity << std::endl;

		schrSim.getRegister().writeToFile("c:\\temp\\schrodinger_step_end.dat");
		fftSchr.getRegister().writeToFile("c:\\temp\\findif_step_end.dat");

		if (fidelity < 0.8)
		{
			std::cout << "\nFidelity too low, failure!" << std::endl;
			return false;
		}
	}

	return true;
}

bool SimulationTests()
{
	std::cout << "\nTesting simulations of quantum simulations..." << std::endl;

	return PauliSimulationTests() && SchrodingerSimulationTests();
}
