#include <iostream>
#include <iterator>

#include "Tests.h"
#include "Teleportation.h"
#include "EntanglementSwapping.h"


bool TeleportTests()
{
	std::cout << "Teleporting with teleport sub-algorithm:" << std::endl;
	QC::SubAlgo::Teleport<> teleport;
	QC::QubitRegister<> reg;

	for (int i = 0; i < 30; ++i)
	{
		if (dist_bool(gen))
		{
			reg.setToBasisState(1); // set the qubit to teleport to 'up'
			teleport.Execute(reg);
			const unsigned int state = reg.MeasureQubit(2);

			std::cout << "Teleported 1, measured: " << state << std::endl;

			if (state != 1) return false;
		}
		else
		{
			reg.setToBasisState(0); // set the qubit to teleport to 'down'
			teleport.Execute(reg);
			const unsigned int state = reg.MeasureQubit(2);

			std::cout << "Teleported 0, measured: " << state << std::endl;

			if (state != 0) return false;
		}
	}

	std::cout << "Now testing teleporting a generic state with teleport sub-algorithm:" << std::endl;

	for (int i = 0; i < 16; ++i)
	{
		// generate and normalize
		std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));

		const double norm = sqrt(std::norm(alpha) + std::norm(beta));
		alpha /= norm;
		beta /= norm;

		reg.Clear();
		reg.setRawAmplitude(0, alpha);
		reg.setRawAmplitude(1, beta);
		reg.Normalize();

		std::cout << "Teleporting " << alpha << "|0> + " << beta << "|1>";

		const unsigned int classicalBits = teleport.Execute(reg);
		std::cout << " Measured values for the two qubits: " << classicalBits;

		// how is the whole thing looking before Bob's measurement?
		const std::complex<double> receivedAlpha = reg.getBasisStateAmplitude(classicalBits);
		const std::complex<double> receivedBeta = reg.getBasisStateAmplitude(0x4 | classicalBits);
		std::cout << "... Teleported state: " << receivedAlpha << "|0> + " << receivedBeta << "|1>" << std::endl;

		if (!approxEqual(alpha, receivedAlpha) || !approxEqual(beta, receivedBeta)) return false;
	}

	return true;
}

bool EntanglementSwappingTests()
{
	std::cout << "Entanglement swapping:" << std::endl;
	Teleportation::EntanglementSwapping entSwap;
	QC::QubitRegister reg(2);
	QC::BellState bellState;
	std::complex<double> zero(0, 0);

	for (int i = 0; i < 16; ++i)
	{
		const bool s1 = dist_bool(gen);
		const bool s2 = dist_bool(gen);

		std::cout << "Entanglement swapping the Bell" << (s1 ? "1" : "0") << (s2 ? "1" : "0") << " state..." << std::endl;

		// this is the entangled pair to be swapped, should end up in the qubits 5 and 6
		// in reg is the state to be compared with
		bellState.setBellState(reg, 0, 1, s1, s2);

		entSwap.SetBellStateToBeSent(s1, s2);
		const unsigned int meas = entSwap.Teleport(i < 8 ? false : true);

		const Eigen::VectorXcd& regStorage = entSwap.getRegisterStorage();

		for (unsigned int state = 0; state < regStorage.size(); ++state)
		{
			if ((state & 0xf) == meas)
			{
				const unsigned int redState = state >> 4;
				if (!approxEqual(regStorage[state], reg.getBasisStateAmplitude(redState)))
				{
					std::cout << "Entanglement swapping failed: Expecting " << reg.getBasisStateAmplitude(redState) << " but got " << regStorage[state] << " for state: " << redState << std::endl;
					return false;
				}

				std::cout << "Sent " << reg.getBasisStateAmplitude(redState) << " received " << regStorage[state] << " for state: " << redState << std::endl;
			}
			else if (!approxEqual(regStorage[state], zero))
			{
				std::cout << "Surprise: Expecting zero!" << std::endl;
				return false;
			}
		}
		std::cout << std::endl;
	}

	return true;
}


bool TeleportationTests()
{
	std::cout << "\nTesting teleportation..." << std::endl;

	Teleportation::QuantumTeleportationRealization<> qt;

	for (int i = 0; i < 30; ++i)
	{
		if (dist_bool(gen))
		{
			qt.SetState(0, 1); // set the qubit to teleport to 'up'
			const unsigned int state = qt.Execute();

			std::cout << "Teleported 1, measured: " << state << std::endl;

			if (state != 1) return false;
		}
		else
		{
			qt.SetState(1, 0); // set the qubit to teleport to 'down'
			const unsigned int state = qt.Execute();

			std::cout << "Teleported 0, measured: " << state << std::endl;

			if (state != 0) return false;
		}
	}

	std::cout << "Now testing teleporting a generic state:" << std::endl;

	for (int i = 0; i < 16; ++i)
	{
		// generate and normalize
		std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));

		const double norm = sqrt(std::norm(alpha) + std::norm(beta));
		alpha /= norm;
		beta /= norm;

		qt.SetState(alpha, beta);
		std::cout << "Teleporting " << alpha << "|0> + " << beta << "|1>";

		const unsigned int classicalBits = qt.Teleport(i < 8 ? false : true); // also test sending explicitely the two classical bits for half the tests, although it should not make a difference
		std::cout << " Measured values for the two qubits: " << classicalBits;

		// how is the whole thing looking before Bob's measurement?
		const std::complex<double> receivedAlpha = qt.getBasisStateAmplitude(classicalBits);
		const std::complex<double> receivedBeta = qt.getBasisStateAmplitude(0x4 | classicalBits);
		std::cout << "... Teleported state: " << receivedAlpha << "|0> + " << receivedBeta << "|1>" << std::endl;

		if (!approxEqual(alpha, receivedAlpha) || !approxEqual(beta, receivedBeta)) return false;
	}

	return TeleportTests() && EntanglementSwappingTests();
}
