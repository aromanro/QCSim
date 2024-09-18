// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iterator>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Tests.h"


int main()
{
	while (true)
	{
		std::cout << "\n\nPick an option:" << std::endl;
		std::cout << "0 to exit" << std::endl;
		std::cout << "1 to run statevector tests" << std::endl;
		std::cout << "2 to run MPS simulator tests" << std::endl;
		std::cout << "3 to run all tests" << std::endl;
		std::cout << "Command: ";
		std::string dummy;
		getline(std::cin, dummy);
		if (dummy.empty()) continue;
		const char c = dummy.size() == 1 ? dummy[0] : 'x';

		if (c == '0') break;
		else if (c >= '1' && c <= '3') tests(c - '1');
		else
			std::cout << "Sorry, unknown command" << std::endl;
	}
	
	return 0;
}





