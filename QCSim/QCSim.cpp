// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iterator>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Tests.h"

int ExecuteTests(const std::string& arg1)
{
	if (arg1 == "0") return 0;
	else if (arg1 == "1")
		return tests(0) ? 0 : -1;
	else if (arg1 == "2")
		return tests(1) ? 0 : -1;
	else if (arg1 == "3")
		return tests(2) ? 0 : -1;
	else if (arg1 == "4")
		return tests(3) ? 0 : -1;
	else if (arg1 == "5")
		return tests(4) ? 0 : -1;
	else
	{
		std::cout << "Unknown command line argument, options are 0, 1, 2, 3, 4" << std::endl;
		return -1;
	}

	return 0;
}

int main(int argc, char* argv[])
{
	if (argc > 1)
	{
		const std::string arg1 = argv[1];

		return ExecuteTests(arg1);
	}

	while (true)
	{
		std::cout << "\n\nPick an option:" << std::endl;
		std::cout << "0 to exit" << std::endl;
		std::cout << "1 to run statevector tests" << std::endl;
		std::cout << "2 to run MPS simulator tests" << std::endl;
		std::cout << "3 to run Clifford simulator tests" << std::endl;
		std::cout << "4 to run Pauli propagator tests" << std::endl;
		std::cout << "5 to run all tests" << std::endl;
		std::cout << "Command: ";
		std::string dummy;
		getline(std::cin, dummy);
		if (dummy.empty()) continue;
		const char c = dummy.size() == 1 ? dummy[0] : 'x';

		if (c == '0') break;
		else if (c >= '1' && c <= '5') tests(c - '1');
		else
			std::cout << "Sorry, unknown command" << std::endl;
	}
	
	return 0;
}





