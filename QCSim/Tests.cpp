#include <iostream>
#include <map>

#include "QubitRegister.h"

#include "GroverAlgorithm.h"
#include "ShorAlgorithm.h"
#include "SuperdenseCoding.h"
#include "CheckCHSHInequality.h"
#include "BernsteinVazirani.h"
#include "QuantumCryptograpy.h"
#include "DeutschJozsa.h"
#include "SimonAlgorithm.h"

#define TESTS_CPP_ 1
#include "Tests.h"

std::random_device rd;
std::mt19937 gen(rd());

std::bernoulli_distribution dist_bool;
std::uniform_real_distribution<> dist_ampl(-1., 1.);

bool ShorTests()
{
	std::cout << "\nTesting Shor..." << std::endl;

	/*
	{
		// check f
		std::map<int, int> measurements;
		const int nrMeasurements = 100;

		Shor::ShorAlgorithm shorAlgo;

		// with 2 (default) f should be 1, 2, 4, 8 for performing period finding
		//shorAlgo.setA(7); // with 7 should be 1, 7, 4, 13
		shorAlgo.setA(4); // should be 1, 4
		//shorAlgo.setA(8); // 1, 8, 4, 2
		//shorAlgo.setA(11); // 1, 11
		//shorAlgo.setA(13); // 1, 13, 4, 7
		//shorAlgo.setA(14); // 1, 14

		std::map<int, int> fmeasurements;

		for (int i = 0; i < nrMeasurements; ++i)
		{
			unsigned int state = shorAlgo.Execute();
			++measurements[state & 0x7]; //only the l bits matter

			// this should work only if Execute calls full register measurement, otherwise use the commented part below
			++fmeasurements[state >> 3]; //save the f register measurements for debugging

			//state = shorAlgo.reg.Measure(3, 6);
			//++fmeasurements[state];
		}

		std::cout << "f register (with 4 should be 1, 4)" << std::endl;
		for (auto m : fmeasurements)
			std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

		std::cout << "x register" << std::endl;
		for (auto m : measurements)
			std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
	}
	*/
	
	for (int i = 0; i < 10;)
	{
		Shor::ShorAlgorithm<> shorAlgo;//(21, 14, 9)
		unsigned int p1;
		unsigned int p2;
		bool res = shorAlgo.factorize(p1, p2);
		if (p2 > p1) std::swap(p1, p2);

		std::cout << (res ? "Quantum algo: " : "Classical algo: ") << p1 << " " << p2 << std::endl;

		if (p1 != 5 && p2 != 3) return false;

		if (res) ++i;
	}

	// this is slow, so it stays commented out, it's only for occasional tests
	/*
	for (;;) {
		Shor::ShorAlgorithm shorAlgo(21, 14, 9);
		unsigned int p1;
		unsigned int p2;
		bool res = shorAlgo.factorize(p1, p2);
		if (p2 > p1) std::swap(p1, p2);

		std::cout << (res ? "Quantum algo: " : "Classical algo: ") << p1 << " " << p2 << std::endl;

		if (p1 != 7 && p2 != 3) return false;

		if (res) break; // managed to factorize using the quantum algo
	}
	*/

	return true;
}

bool GroverTests()
{
	std::cout << "\nTesting Grover..." << std::endl;

	bool res = true;

	std::map<int, int> measurements;
	const int nrMeasurements = 100;

	std::uniform_int_distribution<> dist(0, 0xf);

	for (int i = 0; i < 5; ++i)
	{
		const int ans = dist(gen);

		std::cout << "Testing for answer: " << ans << std::endl;

		measurements.clear();
		Grover::GroverAlgorithm<> galgo(8);
		galgo.setCorrectQuestionState(ans);
		for (int j = 0; j < nrMeasurements; ++j)
		{
			const unsigned int state = galgo.Execute();
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



bool SuperdenseCodingTests()
{
	std::cout << "\nTesting superdense coding..." << std::endl;

	Coding::SuperdenseCoding<> coding;

	for (int i = 0; i < 30; ++i)
	{
		// bit1 is on position 0, bit2 is on position 1
		const bool bit1 = dist_bool(gen);
		const bool bit2 = dist_bool(gen);
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

	BellInequalities::CheckCHSHInequality<> test;
	bool separateMeasurements = false;
	for (int i = 0; i < 20; ++i)
	{
		if (i >= 10) separateMeasurements = true;
		test.ResetStatistics();
		for (int j = 0; j < 100000; ++j)
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

	std::cout << "Three qubits:" << std::endl;
	BernsteinVazirani::BernsteinVaziraniAlgorithm<> bv;

	unsigned int b0 = 1;
	unsigned int b1 = 2;
	unsigned int b2 = 4;

	for (int i = 0; i < 10; ++i)
	{
		const bool bit0 = dist_bool(gen);
		const bool bit1 = dist_bool(gen);
		const bool bit2 = dist_bool(gen);

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

	BernsteinVazirani::BernsteinVaziraniAlgorithm<> bvBig(6);
	for (int i = 0; i < 30; ++i)
	{
		const bool bit0 = dist_bool(gen);
		const bool bit1 = dist_bool(gen);
		const bool bit2 = dist_bool(gen);
		const bool bit3 = dist_bool(gen);
		const bool bit4 = dist_bool(gen);
		const bool bit5 = dist_bool(gen);

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

	QuantumCryptograpy::BB84Protocol<> bb84;

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
	DeutschJozsa::DeutschJozsaAlgorithm<> deutsch;
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
	DeutschJozsa::DeutschJozsaAlgorithm<> deutschJozsa(9);
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


bool SimonTests()
{
	std::cout << "\nTesting Simon..." << std::endl;

	for (unsigned int nrQubits = 2; nrQubits <= 4; ++nrQubits)
	{
		std::cout << "Nr of qubits: " << nrQubits << std::endl;

		// also test function generation
		Simon::Oracle<> oracle;

		Simon::SimonAlgorithm<> simonAlgorithm(nrQubits);

		const unsigned int lim = (1 << nrQubits) - 1;
		for (unsigned int functionString = 0; functionString <= lim; ++functionString)
		{
			std::cout << "Trying with string: " << functionString << "...";

			// function testing

			oracle.setString(functionString, nrQubits);
			if (!oracle.checkFunction())
			{
				std::cout << "\n Something is wrong with function contruction, check the oracle code" << std::endl;

				return false;
			}

			// ***************

			simonAlgorithm.setString(functionString);

			const unsigned int res = simonAlgorithm.Execute();

			if (res != functionString)
			{
				std::cout << "\n Result different that the set string: " << res << std::endl;

				return false;
			}

			std::cout << " ok" << std::endl;
		}
	}

	return true;
}

bool basicTests()
{
	bool res = registerMeasurementsTests();
	if (res) res = checkGates();
	if (res) res = ErrorCorrectionTests();
	if (res) res = BellInequalitiesTests();

	return res;
}

bool tests()
{
	/*
	QC::Gates::HadamardGate g;
	std::cout << g.getRawOperatorMatrix() << std::endl << std::endl;
	std::cout << g.getOperatorMatrix(2) << std::endl << std::endl;
	std::cout << g.getOperatorMatrix(3) << std::endl << std::endl;
	std::cout << g.getOperatorMatrix(2, 1) << std::endl << std::endl;
	std::cout << g.getOperatorMatrix(3, 1) << std::endl << std::endl;
	std::cout << g.getOperatorMatrix(3, 2) << std::endl << std::endl;
	exit(0);
	*/

	std::cout << "\nTests\n";

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	bool res = basicTests() && quantumAdderTests() && DeutschJozsaTests();
	if (res) res = SimonTests() && BernsteinVaziraniTests() && GroverTests();
	if (res) res = PhaseEstimationTests() && ShorTests() && TeleportationTests();
	if (res) res = SuperdenseCodingTests() && QuantumCryptograpyTests() && SimulationTests();
	if (res) res = ParadoxesTests() && GamesTests() && distributedTests();

	auto dif = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();

	std::cout << "\nTesting took: " << dif / 1000. << " seconds!" << std::endl << "\nTests " << (res ? "succeeded" : "failed") << std::endl;

	return res;
}