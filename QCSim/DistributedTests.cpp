#include <iostream>
#include <iterator>
#include <map>

#include "Tests.h"

#include "DistributedCNOT.h"
#include "TeleportedCNOT.h"


bool distributedCNOTTest()
{
	std::cout << "\nTesting distributed CNOT..." << std::endl;


	Distributed::DistributedCNOT distCnot;

	// for non-distributed CNOT
	QC::Gates::CNOTGate cnot;
	QC::QubitRegister reg(2);

	for (unsigned int ctrlQubit = 0; ctrlQubit <= 1; ++ctrlQubit)
		for (unsigned int targetQubit = 0; targetQubit <= 1; ++targetQubit)
		{
			// distributed CNOT
			distCnot.setToBasisState((targetQubit << 3) | ctrlQubit);
			const unsigned int measurements = distCnot.Execute() << 1;
			const Eigen::VectorXcd distStorage = distCnot.getRegisterStorage();

			// now do the same thing but not distributed
			reg.setToBasisState((targetQubit << 1) | ctrlQubit);
			reg.ApplyGate(cnot, 1, 0);
			const Eigen::VectorXcd res = reg.getRegisterStorage();
			
			// now check them to have the same results
			for (int i = 0; i < 4; ++i)
			{
				if (!approxEqual(distStorage((i & 2) << 2 | measurements | (i & 1)), res(i)))
				{
					std::cout << "Error: distributed CNOT failed for ctrl = " << ctrlQubit << ", target = " << targetQubit << ", results were different than for the local CNOT" << std::endl;

					return false;
				}

				// this should not happen, but just in case... to catch bugs
				for (int q1 = 0; q1 < 2; ++q1)
					for (int q2 = 0; q2 < 2; ++q2)
					{
						const unsigned int state = ((q2 << 1) | q1) << 1;
						if (state == measurements)
							continue;

						if (!approxEqual(distStorage(((i & 2) << 2) | state | (i & 1)), 0.0))
						{
							std::cout << "Error: distributed CNOT failed for ctrl = " << ctrlQubit << ", target = " << targetQubit << ", amplitude that should be zero is not" << std::endl;
							return false;
						}
					}
			}
		}

	for (int i = 0; i < 16; ++i)
	{
		std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> gamma(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> delta(dist_ampl(gen), dist_ampl(gen));

		distCnot.Clear();
		distCnot.setRawAmplitude(0, alpha);
		distCnot.setRawAmplitude(1, beta);
		distCnot.setRawAmplitude(0x8, gamma);
		distCnot.setRawAmplitude(0x9, delta);
		distCnot.Normalize();
		const unsigned int measurements = distCnot.Execute() << 1;

		reg.Clear();
		reg.setRawAmplitude(0, alpha);
		reg.setRawAmplitude(1, beta);
		reg.setRawAmplitude(2, gamma);
		reg.setRawAmplitude(3, delta);
		reg.Normalize();
		reg.ApplyGate(cnot, 1, 0);

		const std::complex<double> v1 = distCnot.getBasisStateAmplitude(measurements);
		const std::complex<double> v2 = distCnot.getBasisStateAmplitude(measurements | 1);
		const std::complex<double> v3 = distCnot.getBasisStateAmplitude(0x8 | measurements);
		const std::complex<double> v4 = distCnot.getBasisStateAmplitude(0x8 | measurements | 1);

		if (!approxEqual(v1, reg.getBasisStateAmplitude(0)) || !approxEqual(v2, reg.getBasisStateAmplitude(1)) ||
			!approxEqual(v3, reg.getBasisStateAmplitude(2)) || !approxEqual(v4, reg.getBasisStateAmplitude(3)))
		{
			std::cout << "Error: distributed CNOT failed for random state" << std::endl;

			std::cout << "Measured: " << measurements << std::endl << std::endl;
			reg.displayRegister();
			std::cout << std::endl << std::endl;
			distCnot.displayRegister();

			return false;
		}
	}

	std::cout << "ok" << std::endl;	

	return true;
}

bool teleportedCNOTTest()
{
	std::cout << "\nTesting teleported CNOT..." << std::endl;


	Distributed::TeleportedCNOT teleportedCnot;

	// for non-distributed CNOT
	QC::Gates::CNOTGate cnot;
	QC::QubitRegister reg(2);

	for (unsigned int ctrlQubit = 0; ctrlQubit <= 1; ++ctrlQubit)
		for (unsigned int targetQubit = 0; targetQubit <= 1; ++targetQubit)
		{
			// teleported CNOT
			teleportedCnot.setToBasisState(targetQubit | (ctrlQubit << 5));
			const unsigned int measurements = teleportedCnot.Execute();
			const Eigen::VectorXcd distStorage = teleportedCnot.getRegisterStorage();

			// now do the same thing but not teleported
			reg.setToBasisState(targetQubit | (ctrlQubit << 1));
			reg.ApplyGate(cnot, 0, 1);
			const Eigen::VectorXcd res = reg.getRegisterStorage();

			// now check them to have the same results
			for (int i = 0; i < 4; ++i)
			{
				if (!approxEqual(distStorage((i << 2) | measurements), res(i)))
				{
					std::cout << "Error: teleported CNOT failed for ctrl = " << ctrlQubit << ", target = " << targetQubit << ", results were different than for the local CNOT" << std::endl;

					
					std::cout << "Measured: " << measurements << std::endl << std::endl;
					reg.displayRegister();
					std::cout << std::endl << std::endl;
					teleportedCnot.displayRegister();
					

					return false;
				}
			}
		}

	for (int i = 0; i < 16; ++i)
	{
		std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> gamma(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> delta(dist_ampl(gen), dist_ampl(gen));

		teleportedCnot.Clear();
		teleportedCnot.setRawAmplitude(0, alpha);
		teleportedCnot.setRawAmplitude(1, beta);
		teleportedCnot.setRawAmplitude(0x20, gamma);
		teleportedCnot.setRawAmplitude(0x21, delta);
		teleportedCnot.Normalize();
		const unsigned int measurements = teleportedCnot.Execute();
		
		reg.Clear();
		reg.setRawAmplitude(0, alpha);
		reg.setRawAmplitude(1, beta);
		reg.setRawAmplitude(2, gamma);
		reg.setRawAmplitude(3, delta);
		reg.Normalize();
		reg.ApplyGate(cnot, 0, 1);

		const std::complex<double> v1 = teleportedCnot.getBasisStateAmplitude(measurements);
		const std::complex<double> v2 = teleportedCnot.getBasisStateAmplitude(measurements | 0x4);
		const std::complex<double> v3 = teleportedCnot.getBasisStateAmplitude(measurements | 0x8);
		const std::complex<double> v4 = teleportedCnot.getBasisStateAmplitude(measurements | 0xC);
		
		if (!approxEqual(v1, reg.getBasisStateAmplitude(0)) || !approxEqual(v2, reg.getBasisStateAmplitude(1)) ||
			!approxEqual(v3, reg.getBasisStateAmplitude(2)) || !approxEqual(v4, reg.getBasisStateAmplitude(3))) 
		{
			std::cout << "Error: teleported CNOT failed for random state" << std::endl;

			std::cout << "Measured: " << measurements << std::endl << std::endl;
			reg.displayRegister();
			std::cout << std::endl << std::endl;
			teleportedCnot.displayRegister();

			return false;
		}
	}

	std::cout << "ok" << std::endl;

	return true;
}

bool distributedTests()
{
	std::cout << "\nTesting distributed quantum computing..." << std::endl;

	return distributedCNOTTest() && teleportedCNOTTest();
}