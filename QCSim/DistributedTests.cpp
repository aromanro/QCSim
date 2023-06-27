#include "Tests.h"

#include "DistributedCNOT.h"
#include "TeleportedCNOT.h"

#include <iostream>
#include <map>


bool distributedCNOTTest()
{
	std::cout << "\nTesting distributed CNOT..." << std::endl;


	Distributed::DistributedCNOT distCnot;

	// for non-distributed CNOT
	QC::Gates::CNOTGate cnot;
	QC::QubitRegister reg(2);

	// TODO: maybe instead of using basis states, use random states... but for now this should suffice
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
				if (!approxEqual(distStorage((i & 2) << 2 | measurements | i & 1), res(i)))
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

						if (!approxEqual(distStorage(((i & 2) << 2) | state | i & 1), 0.0))
						{
							std::cout << "Error: distributed CNOT failed for ctrl = " << ctrlQubit << ", target = " << targetQubit << ", amplitude that should be zero is not" << std::endl;
							return false;
						}
					}
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

	// TODO: maybe instead of using basis states, use random states... but for now this should suffice
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

	std::cout << "ok" << std::endl;

	return true;
}

bool distributedTests()
{
	std::cout << "\nTesting distributed quantum computing..." << std::endl;

	return distributedCNOTTest() && teleportedCNOTTest();
}