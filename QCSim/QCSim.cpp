// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iterator>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Tests.h"


int main(int argc, char* argv[])
{
	if (argc > 1)
	{
		const std::string arg1 = argv[1];

		if (arg1 == "0" && !tests(0)) return -1;
		else if (arg1 == "1" && !tests(1)) return -1;
		else if (arg1 == "2" && !tests(2)) return -1;
		else if (arg1 == "3" && !tests(3)) return -1;
		else if (arg1 == "4" && !tests(4)) return -1;
		else
		{
			std::cout << "Unknown command line argument, options are 0, 1, 2, 3, 4" << std::endl;
			return -1;
		}

		return 0;
	}

	while (true)
	{
		std::cout << "\n\nPick an option:" << std::endl;
		std::cout << "0 to exit" << std::endl;
		std::cout << "1 to run statevector tests" << std::endl;
		std::cout << "2 to run MPS simulator tests" << std::endl;
		std::cout << "3 to run Clifford simulator tests" << std::endl;
		std::cout << "4 to run all tests" << std::endl;
		std::cout << "Command: ";
		std::string dummy;
		getline(std::cin, dummy);
		if (dummy.empty()) continue;
		const char c = dummy.size() == 1 ? dummy[0] : 'x';

		if (c == '0') break;
		else if (c >= '1' && c <= '4') tests(c - '1');
		else
			std::cout << "Sorry, unknown command" << std::endl;
	}
	
	return 0;
}





