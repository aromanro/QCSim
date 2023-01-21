#include "Tests.h"
#include "QuantumAdder.h"

#include <iostream>
#include <map>

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

bool quantumAdderTests()
{
	std::cout << "\nTesting quantum adders..." << std::endl;

	if (!quantumHalfAdderTests() || !quantumFullAdderTests()) return false;

	std::cout << "Adding 3-qubit values..." << std::endl;

	std::random_device rd;
	std::mt19937 gen(rd());
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
