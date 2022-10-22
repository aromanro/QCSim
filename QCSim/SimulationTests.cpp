#include "Tests.h"

#include "QubitRegister.h"
#include "PauliStringSimulation.h"

#include <iostream>
#include <map>


void Normalize(Eigen::VectorXcd& regVals, unsigned int nrBasisStates)
{
	std::complex<double> v = abs(regVals[0]) > 1E-5 ? regVals[0] : std::complex<double>(1, 0);
	std::complex<double> accum(0, 0);
	for (unsigned int i = 0; i < nrBasisStates; ++i)
	{
		regVals[i] /= v;

		accum += regVals[i] * std::conj(regVals[i]);
	}
	const double norm = 1. / sqrt(accum.real());
	for (unsigned int i = 0; i < nrBasisStates; ++i)
		regVals[i] *= norm;
}

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

bool SimulationTests()
{
	std::cout << "\nTesting quantum simulations..." << std::endl;

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

		std::complex<double> v = abs(regValsEx[0]) > 1E-5 ? regValsEx[0] : std::complex<double>(1, 0);
		std::complex<double> accum(0);
		for (unsigned int i = 0; i < sim.getNrBasisStates(); ++i)
		{
			regValsEx[i] /= v;

			accum += regValsEx[i] * std::conj(regValsEx[i]);
		}
		double norm = 1. / sqrt(accum.real());
		for (unsigned int i = 0; i < sim.getNrBasisStates(); ++i)
			regValsEx[i] *= norm;

		sim.setToBasisState(0);

		sim.Execute();

		Eigen::VectorXcd regVals = sim.getRegisterStorage();

		Normalize(regVals, sim.getNrBasisStates());

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
