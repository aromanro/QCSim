#include "Tests.h"
#include "CoinFlipping.h"
#include "MagicSquare.h"

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

bool MagicSquareTests()
{
	std::cout << "\nTesting magic square pseudo-telepathy game..." << std::endl;

	Games::MagicSquare<> magicSquare;
	unsigned int wins = magicSquare.Execute();

	if (wins != magicSquare.getNrPlays())
	{
		std::cout << "Something went wrong! Alice and Bob won " << wins << " times, but they played " << magicSquare.getNrPlays() << " times!" << std::endl;

		return false;
	}

	std::cout << "Alice and Bob won each time!" << std::endl;

	return true;
}

bool GamesTests()
{
	return CoinFlippingTests() && MagicSquareTests();
}