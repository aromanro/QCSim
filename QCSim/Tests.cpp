#include "Tests.h"

#include "QubitRegister.h"

#include "GroverAlgorithm.h"
#include "ShorAlgorithm.h"
#include "Teleportation.h"
#include "SuperdenseCoding.h"
#include "CheckCHSHInequality.h"
#include "BernsteinVazirani.h"
#include "QuantumCryptograpy.h"
#include "DeutschJozsa.h"


#include <iostream>
#include <map>


bool ShorTests()
{
	std::cout << "\nTesting Shor..." << std::endl;

	for (int i = 0; i < 10;)
	{
		Shor::ShorAlgorithm shorAlgo;//(21, 14, 9)
		unsigned int p1;
		unsigned int p2;
		bool res = shorAlgo.factorize(p1, p2);
		if (p2 > p1) std::swap(p1, p2);

		std::cout << (res ? "Quantum algo: " : "Classical algo: ") << p1 << " " << p2 << std::endl;

		if (p1 != 5 && p2 != 3) return false;

		if (res) ++i;
	}

	return true;
}

bool GroverTests()
{
	std::cout << "\nTesting Grover..." << std::endl;

	bool res = true;

	std::map<int, int> measurements;
	const int nrMeasurements = 100;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist(0, 0xf);

	for (int i = 0; i < 5; ++i)
	{
		const int ans = dist(gen);

		std::cout << "Testing for answer: " << ans << std::endl;

		measurements.clear();
		Grover::GroverAlgorithm galgo(8);
		galgo.setCorrectQuestionState(ans);
		for (int i = 0; i < nrMeasurements; ++i)
		{
			unsigned int state = galgo.Execute();
			++measurements[state];
		}

		bool found = false;
		for (auto m : measurements)
		{
			std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

			if (m.first == ans)
			{
				found = true;
				if (static_cast<double>(m.second) / nrMeasurements < 0.9)
					res = false;
			}
		}

		if (!found) res = false;
	}

	return res;
}


bool TeleportationTests()
{
	std::cout << "\nTesting teleportation..." << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_bool(0, 1);

	Teleportation::QuantumTeleportationRealization qt;

	for (int i = 0; i < 30; ++i)
	{
		if (dist_bool(gen))
		{
			qt.SetState(0, 1); // set the qubit to teleport to 'up'
			unsigned int state = qt.Execute();

			std::cout << "Teleported 1, measured: " << state << std::endl;

			if (state != 1) return false;
		}
		else
		{
			qt.SetState(1, 0); // set the qubit to teleport to 'down'
			unsigned int state = qt.Execute();

			std::cout << "Teleported 0, measured: " << state << std::endl;

			if (state != 0) return false;
		}
	}

	std::uniform_real_distribution<> dist_ampl(-1., 1.);

	std::cout << "Now testing teleporting a generic state:" << std::endl;

	for (int i = 0; i < 16; ++i)
	{
		// generate and normalize
		std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
		std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));

		const double norm = sqrt((alpha * std::conj(alpha) + beta * std::conj(beta)).real());
		alpha /= norm;
		beta /= norm;

		qt.SetState(alpha, beta);
		std::cout << "Teleporting " << alpha << "|0> + " << beta << "|1>";

		unsigned int classicalBits = qt.Teleport(i < 8 ? false : true); // also test sending explicitely the two classical bits for half the tests, although it should not make a difference
		std::cout << " Measured values for the two qubits: " << classicalBits;

		// how is the whole thing looking before Bob's measurement?
		std::complex<double> receivedAlpha = qt.getBasisStateAmplitude(classicalBits);
		std::complex<double> receivedBeta = qt.getBasisStateAmplitude(0x4 | classicalBits);
		std::cout << "... Teleported state: " << receivedAlpha << "|0> + " << receivedBeta << "|1>" << std::endl;

		if (!approxEqual(alpha, receivedAlpha) || !approxEqual(beta, receivedBeta)) return false;
	}

	return true;
}

bool SuperdenseCodingTests()
{
	std::cout << "\nTesting superdense coding..." << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_bool(0, 1);

	Coding::SuperdenseCoding coding;

	for (int i = 0; i < 30; ++i)
	{
		// bit1 is on position 0, bit2 is on position 1
		const bool bit1 = dist_bool(gen) == 1;
		const bool bit2 = dist_bool(gen) == 1;
		coding.SetBits(bit1, bit2);
		unsigned int classicalBits = coding.Execute();

		const bool recvbit1 = (classicalBits & 1) != 0;
		const bool recvbit2 = (classicalBits & 2) != 0;
		std::cout << "Sent " << (bit1 ? "1" : "0") << (bit2 ? "1" : "0") << " using superdense coding, received: " << (recvbit1 ? "1" : "0") << (recvbit2 ? "1" : "0") << std::endl;

		if (bit1 != recvbit1 || bit2 != recvbit2) return false;
	}

	return true;
}

bool BellInequalitiesTests()
{
	std::cout << "\nTesting CHSH inequality..." << std::endl;

	static const double expected = 2. * sqrt(2.);

	BellInequalities::CheckCHSHInequality test;
	bool separateMeasurements = false;
	for (int i = 0; i < 20; ++i)
	{
		if (i >= 10) separateMeasurements = true;
		test.ResetStatistics();
		for (int i = 0; i < 100000; ++i)
			test.Check(separateMeasurements);

		const double val = test.getValue();
		std::cout << (separateMeasurements ? "Measurements separated" : "Measurements together") << " Value: " << val << (val <= 2. ? " Inequality obeyed (you won the lottery!)" : " Inequality violated") << std::endl;

		if (val <= 2.) return false;
		else if (!approxEqual(val, expected, 0.1))
		{
			std::cout << "Expected value to be about 2*sqrt(2), got a value that's too different" << std::endl;
			return false;
		}
	}

	return true;
}

bool BernsteinVaziraniTests()
{
	std::cout << "\nTesting Bernstein-Vazirani..." << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_bool(0, 1);

	std::cout << "Three qubits:" << std::endl;
	BernsteinVazirani::BernsteinVaziraniAlgorithm bv;

	unsigned int b0 = 1;
	unsigned int b1 = 2;
	unsigned int b2 = 4;

	for (int i = 0; i < 10; ++i)
	{
		const bool bit0 = dist_bool(gen) == 1;
		const bool bit1 = dist_bool(gen) == 1;
		const bool bit2 = dist_bool(gen) == 1;

		unsigned int str = 0;
		if (bit0) str |= b0;
		if (bit1) str |= b1;
		if (bit2) str |= b2;

		std::cout << "Setting string to: " << str << std::endl;
		bv.setString(str);
		unsigned int state = bv.Execute();

		if (state != str)
		{
			std::cout << "Failed, obtained a different string: " << state << std::endl;
			return false;
		}
	}

	std::cout << "Six qubits:" << std::endl;
	unsigned int b3 = 8;
	unsigned int b4 = 16;
	unsigned int b5 = 32;

	BernsteinVazirani::BernsteinVaziraniAlgorithm bvBig(6);
	for (int i = 0; i < 30; ++i)
	{
		const bool bit0 = dist_bool(gen) == 1;
		const bool bit1 = dist_bool(gen) == 1;
		const bool bit2 = dist_bool(gen) == 1;
		const bool bit3 = dist_bool(gen) == 1;
		const bool bit4 = dist_bool(gen) == 1;
		const bool bit5 = dist_bool(gen) == 1;

		unsigned int str = 0;
		if (bit0) str |= b0;
		if (bit1) str |= b1;
		if (bit2) str |= b2;
		if (bit3) str |= b3;
		if (bit4) str |= b4;
		if (bit5) str |= b5;

		std::cout << "Setting string to: " << str << std::endl;
		bvBig.setString(str);
		unsigned int state = bvBig.Execute();

		if (state != str)
		{
			std::cout << "Failed, obtained a different string: " << state << std::endl;
			return false;
		}
	}

	return true;
}

bool QuantumCryptograpyTests()
{
	std::cout << "\nTesting BB84 protocol..." << std::endl;

	QuantumCryptograpy::BB84Protocol bb84;

	for (int i = 0; i < 30; ++i)
	{
		const bool eavesdrop = i > 9;
		const bool randomEavesdrop = i > 19;

		bb84.setEavesdropping(eavesdrop);
		bb84.setRandomEavesdropping(randomEavesdrop);

		std::cout << "Transmission with" << (eavesdrop ? " eavesdropping" : "out eavesdropping") << (eavesdrop ? (randomEavesdrop ? " randomly" : " each time") : "") << std::endl;

		const unsigned int match = bb84.Execute();
		if (match)
		{
			if (eavesdrop)
			{
				std::cout << "Match, but eavesdropped" << std::endl;
				return false;
			}

			// additional check, this is only done for tests, the key remains a secret in 'real life'
			if (!bb84.compareKeys())
			{
				std::cout << "Despite a match for checked bits, the keys do not match" << std::endl;
				return false;
			}
		}
		else
		{
			// mismatch
			if (!eavesdrop)
			{
				std::cout << "Mismatch, but not eavesdropped" << std::endl;
				return false;
			}
		}
		//std::cout << "Key length: " << bb84.getReceivedKey().size() << std::endl;
	}

	return true;
}

bool DeutschTests()
{
	// first the simplest case, Deutsch's algorithm (N=3)
	std::cout << "First, Deutsch's algorithm..." << std::endl;
	DeutschJozsa::DeutschJozsaAlgorithm deutsch;
	for (int i = 0; i < 15; ++i)
	{
		if (i < 5)
			deutsch.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantZero);
		else if (i < 10)
			deutsch.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantOne);
		else
			deutsch.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::balanced);

		unsigned int state = deutsch.Execute();

		if (i < 10)
		{
			// constant
			if (!deutsch.WasConstantResult(state))
			{
				std::cout << "Expected constant result, got balanced" << std::endl;
				return false;
			}
			else std::cout << "Constant " << (i < 5 ? "zero" : "one") << " function set, got constant, ok" << std::endl;
		}
		else
		{
			// balanced
			if (deutsch.WasConstantResult(state))
			{
				std::cout << "Expected balanced result, got constant" << std::endl;
				return false;
			}
			else std::cout << "Balanced function set, got balanced, ok" << std::endl;
		}
	}

	return true;
}

bool DeutschJozsaTests()
{
	std::cout << "\nTesting Deutsch-Jozsa..." << std::endl;

	if (!DeutschTests()) return false;

	std::cout << "Now, general Deutsch-Jozsa..." << std::endl;
	DeutschJozsa::DeutschJozsaAlgorithm deutschJozsa(9);
	for (int i = 0; i < 30; ++i)
	{
		if (i < 10)
			deutschJozsa.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantZero);
		else if (i < 20)
			deutschJozsa.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::constantOne);
		else
			deutschJozsa.setFunction(DeutschJozsa::DeutschJozsaAlgorithm<>::FunctionType::balanced);

		unsigned int state = deutschJozsa.Execute();

		if (i < 20)
		{
			// constant
			if (!deutschJozsa.WasConstantResult(state))
			{
				std::cout << "Expected constant result, got balanced" << std::endl;
				return false;
			}
			else std::cout << "Constant " << (i < 10 ? "zero" : "one") << " function set, got constant, ok" << std::endl;
		}
		else
		{
			// balanced
			if (deutschJozsa.WasConstantResult(state))
			{
				std::cout << "Expected balanced result, got constant" << std::endl;
				return false;
			}
			else std::cout << "Balanced function set, got balanced, ok" << std::endl;
		}
	}

	return true;
}

bool tests()
{
	std::cout << "\nTests\n";

	bool res = registerMeasurementsTests();
	if (res) res = DeutschJozsaTests();
	if (res) res = BernsteinVaziraniTests();
	if (res) res = GroverTests();
	if (res) res = ShorTests();
	if (res) res = TeleportationTests();
	if (res) res = SuperdenseCodingTests();
	if (res) res = QuantumCryptograpyTests();
	if (res) res = BellInequalitiesTests();

	if (res)
	{
		std::cout << "\nTesting simulations of quantum simulations..." << std::endl;

		res = SimulationTests();
	}

	std::cout << "\nTests " << (res ? "succeeded" : "failed") << std::endl;

	return res;
}