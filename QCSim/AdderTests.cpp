#include "Tests.h"
#include "QuantumAdder.h"
#include "DraperAdder.h"

#include <iostream>
#include <map>

extern std::uniform_int_distribution<int> dist_bool;

bool quantumHalfAdderTests()
{
	std::cout << "\nTwo qubits half-adder..." << std::endl;
	QC::QubitRegister regThreeQubits;
	Adders::TwoQubitsHalfAdder halfAdder(0, 1, 2);

	// 00
	std::cout << "Adding 0 + 0...";
	regThreeQubits.setToBasisState(0);
	halfAdder.Execute(regThreeQubits);
	unsigned int res = regThreeQubits.Measure();
	if (res & 1)
	{
		std::cout << " Half-adder measured 1 on the first qubit when adding 0 + 0" << std::endl;
		return false;
	}
	res >>= 1;
	if (res != 0)
	{
		std::cout << " Half-adder result: " << res << " for adding 0 + 0" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	// 01
	std::cout << "Adding 0 + 1...";
	regThreeQubits.setToBasisState(1);
	halfAdder.Execute(regThreeQubits);
	res = regThreeQubits.Measure();
	if ((res & 1) == 0)
	{
		std::cout << " Half-adder measured 0 on the first qubit when adding 1 + 0" << std::endl;
		return false;
	}
	res >>= 1;
	if (res != 1)
	{
		std::cout << " Half-adder result: " << res << " for adding 1 + 0" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	// 10
	std::cout << "Adding 1 + 0...";
	regThreeQubits.setToBasisState(2);
	halfAdder.Execute(regThreeQubits);
	res = regThreeQubits.Measure();
	if (res & 1)
	{
		std::cout << " Half-adder measured 1 on the first qubit when adding 0 + 1" << std::endl;
		return false;
	}
	res >>= 1;
	if (res != 1)
	{
		std::cout << " Half-adder result: " << res << " for adding 0 + 1" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	// 11
	std::cout << "Adding 1 + 1...";
	regThreeQubits.setToBasisState(3);
	halfAdder.Execute(regThreeQubits);
	res = regThreeQubits.Measure();
	if ((res & 1) == 0)
	{
		std::cout << " Half-adder measured 0 on the first qubit when adding 1 + 1" << std::endl;
		return false;
	}
	res >>= 1;
	if (res != 2)
	{
		std::cout << " Half-adder result: " << res << " for adding 1 + 1" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	return true;
}

bool quantumFullAdderTests()
{
	QC::QubitRegister regFourQubits(4);
	Adders::TwoQubitsFullAdder fullAdder(0, 1, 2, 3);
	std::cout << "Testing full adder..." << std::endl;

	// 00
	std::cout << "Adding 0 + 0...";
	regFourQubits.setToBasisState(0);
	fullAdder.Execute(regFourQubits);
	unsigned int res = regFourQubits.Measure();
	if (res & 3)
	{
		std::cout << " Full-adder altered the qubits when adding 0 + 0: " << (res & 3) << std::endl;
		return false;
	}
	res >>= 2;
	if (res != 0)
	{
		std::cout << " Full-adder result: " << res << " for adding 0 + 0" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	// 01
	std::cout << "Adding 0 + 1...";
	regFourQubits.setToBasisState(1);
	fullAdder.Execute(regFourQubits);
	res = regFourQubits.Measure();
	if ((res & 3) != 1)
	{
		std::cout << " Full-adder altered the qubits when adding 1 + 0: " << (res & 3) << std::endl;
		return false;
	}
	res >>= 2;
	if (res != 1)
	{
		std::cout << " Full-adder result: " << res << " for adding 1 + 0" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	// 10
	std::cout << "Adding 1 + 0...";
	regFourQubits.setToBasisState(2);
	fullAdder.Execute(regFourQubits);
	res = regFourQubits.Measure();
	if ((res & 3) != 2)
	{
		std::cout << " Full-adder altered the qubits when adding 0 + 1: " << (res & 3) << std::endl;
		return false;
	}
	res >>= 2;
	if (res != 1)
	{
		std::cout << " Full-adder result: " << res << " for adding 1 + 0" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	// 11
	std::cout << "Adding 1 + 1...";
	regFourQubits.setToBasisState(3);
	fullAdder.Execute(regFourQubits);
	res = regFourQubits.Measure();
	if ((res & 3) != 3)
	{
		std::cout << " Full-adder altered the qubits when adding 1 + 1: " << (res & 3) << std::endl;
		return false;
	}
	res >>= 2;
	if (res != 2)
	{
		std::cout << " Full-adder result: " << res << " for adding 1 + 1" << std::endl;
		return false;
	}
	std::cout << " ok" << std::endl;

	return true;
}

bool NQubitsAdderTests()
{
	std::cout << "NQubitsAdder, adding 3-qubit values..." << std::endl;

	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_nr(0, 7);


	Adders::NQubitsAdderAlgorithm threeQubitsAdder;

	for (int i = 0; i < 10; ++i)
	{
		const unsigned int n1 = dist_nr(gen);
		unsigned int n2 = dist_nr(gen);
		std::cout << "Computing " << n1 << "+" << n2 << "...";

		const unsigned int expected = n1 + n2;
		n2 <<= 3;
		n2 |= n1;
		threeQubitsAdder.setToBasisState(n2);
		unsigned int res = threeQubitsAdder.Execute();
		if ((res & 0x3f) != n2)
		{
			std::cout << " Adder altered the qubits, the input qubits are now: " << (res & 0x3f) << std::endl;
			return false;
		}
		res >>= 6;
		if (res != expected)
		{
			std::cout << " Adder result wrong: " << res << std::endl;
			return false;
		}
		std::cout << " ok" << std::endl;
	}

	return true;
}


bool SimpleDrapperAdderTests()
{
	const unsigned int nQubits = 3;

	std::cout << "Draper adder, adding " << nQubits << "-qubit values..." << std::endl;

	const unsigned int mask = (1 << nQubits) - 1;

	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_nr1(0, 3);
	std::uniform_int_distribution<> dist_nr2(0, 4);

	Adders::DraperAdder adder(nQubits);

	for (int i = 0; i < 20; ++i)
	{
		unsigned int n1 = dist_nr1(gen);
		unsigned int n2 = dist_nr2(gen);
		if (dist_bool(gen)) std::swap(n1, n2); // this allows having the bigger values (if the ones from distributions are not equal) have equal probability in both registers

		std::cout << "Computing " << n1 << "+" << n2 << "...";

		const unsigned int expected = n1 + n2;

		n2 <<= nQubits;
		n2 |= n1;

		std::map<unsigned int, unsigned int> measurements;
		int failures = 0;

		for (int t = 0; t < 100; ++t)
		{
			adder.setToBasisState(n2);
			unsigned int res = adder.Execute();

			if ((res & mask) != n1)
			{
				std::cout << " Adder altered the first qubits, the result is: " << res << std::endl;
				return false;
			}
			res >>= nQubits;

			if (res != expected)
				++failures;

			++measurements[res];
		}

		unsigned int mostFreqRes = measurements.begin()->first;
		unsigned int freqMax = measurements.begin()->second;
		for (auto& v : measurements)
		{
			unsigned int res = v.first;
			unsigned int freq = v.second;
			if (freq > freqMax)
			{
				freqMax = freq;
				mostFreqRes = res;
			}
		}

		if (mostFreqRes != expected)
		{
			std::cout << " Adder result wrong, number of failures/100 tries: " << failures << " Most frequent result: " << mostFreqRes << std::endl;
			std::cout << "All results: " << std::endl;

			for (auto& v : measurements)
				std::cout << v.first << ", " << v.second << " times" << std::endl;

			return false;
		}

		std::cout << " ok, failures/100 tries: " << failures << std::endl;
	}

	return true;
}


bool DrapperAdderWithCarryTests()
{
	const unsigned int nQubits = 3;

	std::cout << "Draper adder with carry, adding " << nQubits << "-qubit values..." << std::endl;

	const unsigned int mask = (1 << nQubits) - 1;

	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist_nr(0, 7);

	Adders::DraperAdderWithCarry adder(nQubits);

	for (int i = 0; i < 30; ++i)
	{
		unsigned int n1 = dist_nr(gen);
		unsigned int n2 = dist_nr(gen);

		std::cout << "Computing " << n1 << "+" << n2 << "...";

		const unsigned int expected = n1 + n2;

		n2 <<= nQubits;
		n2 |= n1;

		std::map<unsigned int, unsigned int> measurements;
		int failures = 0;

		for (int t = 0; t < 100; ++t)
		{
			adder.setToBasisState(n2);
			unsigned int res = adder.Execute();

			if ((res & mask) != n1)
			{
				std::cout << " Adder altered the first qubits, the result is: " << res << std::endl;
				return false;
			}
			res >>= nQubits;

			if (res != expected)
				++failures;

			++measurements[res];
		}

		unsigned int mostFreqRes = measurements.begin()->first;
		unsigned int freqMax = measurements.begin()->second;
		for (auto& v : measurements)
		{
			unsigned int res = v.first;
			unsigned int freq = v.second;
			if (freq > freqMax)
			{
				freqMax = freq;
				mostFreqRes = res;
			}
		}

		if (mostFreqRes != expected)
		{
			std::cout << " Adder result wrong, number of failures/100 tries: " << failures << " Most frequent result: " << mostFreqRes << std::endl;
			std::cout << "All results: " << std::endl;

			for (auto& v : measurements)
				std::cout << v.first << ", " << v.second << " times" << std::endl;

			return false;
		}

		std::cout << " ok, failures/100 tries: " << failures << std::endl;
	}

	return true;
}

bool DraperAdderTests()
{
	return SimpleDrapperAdderTests() && DrapperAdderWithCarryTests();
}

bool quantumAdderTests()
{
	std::cout << "\nTesting quantum adders..." << std::endl;

	return quantumHalfAdderTests() && quantumFullAdderTests() && NQubitsAdderTests() && DraperAdderTests();
}
