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
#include "QuantumCountingAlgorithm.h"

#define TESTS_CPP_ 1
#include "Tests.h"

std::random_device rd;
std::mt19937 gen(rd());

std::bernoulli_distribution dist_bool;
std::uniform_real_distribution<> dist_ampl(-1., 1.);


bool ShorTestsWithoutTensorProduct()
{
	std::cout << "\nTesting Shor without a tensor product..." << std::endl;

	/*
	{
		// check f
		std::map<int, int> measurements;
		const int nrMeasurements = 100;

		Shor::ShorAlgorithmWithoutTensorProduct shorAlgo;

		// with 2 (default) f should be 1, 2, 4, 8 for performing period finding
		shorAlgo.setA(7); // with 7 should be 1, 7, 4, 13
		//shorAlgo.setA(4); // should be 1, 4
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

		std::cout << "f register (with 4 should be 1, 4, with 7 should be 1, 7, 4, 13)" << std::endl;
		for (auto m : fmeasurements)
			std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

		std::cout << "x register (with 7 should be 4, 0, 2, 6)" << std::endl;
		for (auto m : measurements)
			std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
	}
	*/

	for (int i = 0; i < 10;)
	{
		Shor::ShorAlgorithmWithoutTensorProduct<> shorAlgo;
		unsigned int p1;
		unsigned int p2;
		bool res = shorAlgo.factorize(p1, p2);
		if (p2 > p1) std::swap(p1, p2);

		std::cout << (res ? "Quantum algo: " : "Classical algo: ") << p1 << " " << p2 << std::endl;

		if (p1 != 5 && p2 != 3) return false;

		if (res) ++i;
	}


	// this is slow, so it stays commented out, it's only for occasional tests
	// WARNING: needs 22 qubits, I did not test it yet
	/*
	for (;;) {
		Shor::ShorAlgorithmWithoutTensorProduct shorAlgo(21, 9, 5);
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

bool ShorTests()
{
	std::cout << "\nTesting Shor (with tensor product)..." << std::endl;

	/*
	{
		// check f
		std::map<int, int> measurements;
		const int nrMeasurements = 100;

		Shor::ShorAlgorithm shorAlgo;

		// with 2 (default) f should be 1, 2, 4, 8 for performing period finding
		shorAlgo.setA(7); // with 7 should be 1, 7, 4, 13
		//shorAlgo.setA(4); // should be 1, 4
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

		std::cout << "f register (with 4 should be 1, 4, with 7 should be 1, 7, 4, 13)" << std::endl;
		for (auto m : fmeasurements)
			std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

		std::cout << "x register (with 7 should be 4, 0, 2, 6)" << std::endl;
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

	return ShorTestsWithoutTensorProduct();
}


bool GroverWithGatesTests()
{
	std::cout << "\nTesting Grover with the oracle made out of gates..." << std::endl;

	const int nrMeasurements = 500;

	for (int nrQubits = 4; nrQubits <= 6; ++nrQubits)
	{
		std::cout << nrQubits << " qubits" << std::endl;

		std::uniform_int_distribution<> dist(0, (1 << nrQubits) - 1);
		Grover::GroverAlgorithmWithGatesOracle<> galgo(nrQubits);

		for (int i = 0; i < 5; ++i)
		{
			const unsigned int ans = dist(gen);

			std::cout << "Testing for answer: " << ans << std::endl;
			
			galgo.setCorrectQuestionState(ans);
			const auto measurements = galgo.ExecuteWithMultipleMeasurements(nrMeasurements);

			bool found = false;
			for (const auto& m : measurements)
			{
				std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

				if (m.first == ans)
				{
					found = true;
					if (static_cast<double>(m.second) / nrMeasurements < 0.9)
						return false;
				}
			}

			if (!found) return false;
		}
	}

	return true;
}


bool GroverTests()
{
	std::cout << "\nTesting Grover..." << std::endl;

	const int nrMeasurements = 500;

	for (int nrQubits = 4; nrQubits <= 6; ++nrQubits)
	{
		std::cout << nrQubits << " qubits" << std::endl;

		std::uniform_int_distribution<> dist(0, (1 << nrQubits) - 1);
		Grover::GroverAlgorithm<> galgo(nrQubits);

		for (int i = 0; i < 5; ++i)
		{
			const unsigned int ans = dist(gen);

			std::cout << "Testing for answer: " << ans << std::endl;

			galgo.setCorrectQuestionState(ans);
			const auto measurements = galgo.ExecuteWithMultipleMeasurements(nrMeasurements);

			bool found = false;
			for (const auto& m : measurements)
			{
				std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

				if (m.first == ans)
				{
					found = true;
					if (static_cast<double>(m.second) / nrMeasurements < 0.9)
						return false;
				}
			}

			if (!found) return false;
		}
	}

	return GroverWithGatesTests();
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


bool BernsteinVaziraniWithGatesTests()
{
	std::cout << "\nTesting Bernstein-Vazirani with the oracle made out of gates..." << std::endl;

	std::cout << "Three qubits:" << std::endl;
	BernsteinVazirani::BernsteinVaziraniAlgorithmWithGatesOracle<> bv;

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

	BernsteinVazirani::BernsteinVaziraniAlgorithmWithGatesOracle<> bvBig(6);
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

	return BernsteinVaziraniWithGatesTests();
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

bool SimonWithGatesTests()
{
	std::cout << "\nTesting Simon using gates for oracle..." << std::endl;

	for (unsigned int nrQubits = 2; nrQubits <= 4; ++nrQubits)
	{
		std::cout << "Nr of qubits: " << nrQubits << std::endl;

		// also test function generation
		Simon::SimonFunction func;

		Simon::SimonAlgorithmWithGatesOracle<> simonAlgorithm(nrQubits);

		const unsigned int lim = (1 << nrQubits) - 1;
		for (unsigned int functionString = 0; functionString <= lim; ++functionString)
		{
			std::cout << "Trying with string: " << functionString << "...";

			// function testing

			func.setString(functionString, nrQubits);
			if (!func.checkFunction())
			{
				std::cout << "\n Something is wrong with function construction, check the function code" << std::endl;

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


bool SimonTests()
{
	std::cout << "\nTesting Simon..." << std::endl;

	for (unsigned int nrQubits = 2; nrQubits <= 4; ++nrQubits)
	{
		std::cout << "Nr of qubits: " << nrQubits << std::endl;

		// also test function generation
		Simon::SimonFunction func;

		Simon::SimonAlgorithm<> simonAlgorithm(nrQubits);

		const unsigned int lim = (1 << nrQubits) - 1;
		for (unsigned int functionString = 0; functionString <= lim; ++functionString)
		{
			std::cout << "Trying with string: " << functionString << "...";

			// function testing

			func.setString(functionString, nrQubits);
			if (!func.checkFunction())
			{
				std::cout << "\n Something is wrong with function construction, check the function code" << std::endl;

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

	return SimonWithGatesTests();
}

bool CountingTests()
{
	std::cout << "\nTesting Quantum Counting..." << std::endl;

	unsigned int nrGroverQubits = 4;
	unsigned int nrPrecisionQubits = 6;
	unsigned int nrMeasurements = 10000;

	unsigned int nrGroverStates = 1 << nrGroverQubits;
	for (unsigned int nrMarked = 0; nrMarked <= nrGroverStates; ++nrMarked)
	{
		// pick 'nrMarked' states at random:		
		std::vector<unsigned int> states(nrGroverStates);
		std::iota(states.begin(), states.end(), 0);
		std::shuffle(states.begin(), states.end(), std::default_random_engine(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count())));
		states.resize(nrMarked);

		QuantumCounting::QuantumCountingAlgorithm<> quantumCountingAlgorithm(nrPrecisionQubits, nrGroverQubits);
		quantumCountingAlgorithm.SetMarkedStates(states);

		const auto res = quantumCountingAlgorithm.ExecuteWithMultipleMeasurements(nrMeasurements);

		// get the result with the most measurements:
		unsigned int nrMeasured = 0;
		unsigned int state = 0;
		for (const auto& v : res)
		{
			if (v.second > nrMeasured)
			{
				nrMeasured = v.second;
				state = v.first;
			}
		}

		//std::cout << "Measured " << res.size() << " states, most probable state: " << state << " probability: " << static_cast<double>(nrMeasured)/nrMeasurements << std::endl;

		unsigned int approxCnt = quantumCountingAlgorithm.GetCountForState(state);

		std::cout << "Nr of marked states: " << nrMarked << ", approx count: " << approxCnt << ", real theta: " << quantumCountingAlgorithm.GetCorrectThetaForMarkedStates() << ", approx theta: " << quantumCountingAlgorithm.GetThetaForState(state) << std::endl;

		if (approxCnt != nrMarked) return false;
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
	std::cout << "\nTests\n";

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	bool res = basicTests() && quantumAdderTests() && DeutschJozsaTests();
	if (res) res = SimonTests() && BernsteinVaziraniTests() && GroverTests();
	if (res) res = PhaseEstimationTests() && ShorTests() && TeleportationTests();
	if (res) res = SuperdenseCodingTests() && QuantumCryptograpyTests() && SimulationTests();
	if (res) res = ParadoxesTests() && GamesTests() && distributedTests();
	if (res) res = CountingTests();

	auto dif = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();

	std::cout << "\nTesting took: " << dif / 1000. << " seconds!" << std::endl << "\nTests " << (res ? "succeeded" : "failed") << std::endl;

	return res;
}