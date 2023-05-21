#include "Tests.h"
#include "CoinFlipping.h"

#include <iostream>


bool CoinFlippingTests()
{
	std::cout << "\nTesting coin flipping..." << std::endl;

	Games::CoinFlipping<> coinFlipping;

	int wins = 0;
	for (int i = 0; i < 100; ++i)
	{
		const unsigned int state = coinFlipping.Execute();
		if (state)
			++wins;
	}

	std::cout << "Wins: " << wins << (wins == 0 ? " as expected!" : " not expecting that!") << std::endl;

	return wins == 0;
}

bool GamesTests()
{
	return CoinFlippingTests();
}