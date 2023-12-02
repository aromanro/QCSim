#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Games {

	// for details check "Quantum Pseudo-Telepathy" by Gilles Brassard, Anne Broadbent, Alain Tapp, 2004
	// https://arxiv.org/abs/quant-ph/0407221


	// the first two qubits are for Alice, the other two belong to Bob
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class MagicSquare :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		MagicSquare(int addseed = 0)
			: BaseClass(4, addseed), distInd(1, 3), nrPlays(100)
		{
			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);

			const double invSqrt2 = 1. / sqrt(2.);
			const std::complex<double> i(0., 1.);

			MatrixClass A1m{
				{ i, 0., 0., 1. },
				{ 0., -i, 1., 0. },
				{ 0., i, 1., 0. },
				{ 1., 0., 0., i }
			};
			A1m *= invSqrt2;
			A1.setOperator(A1m);

			MatrixClass A2m{
				{ i, 1., 1., i },
				{ -i, 1., -1., i },
				{ i, 1., -1., -i },
				{ -i, 1., 1., -i }
			};
			A2m *= 0.5;
			A2.setOperator(A2m);

			MatrixClass A3m{
				{ -1., -1., -1., 1. },
				{ 1., 1., -1., 1. },
				{ 1., -1., 1., 1. },
				{ 1., -1., -1., -1. }
			};
			A3m *= 0.5;
			A3.setOperator(A3m);

			MatrixClass B1m{
				{ i, -i, 1., 1. },
				{ -i, -i, 1., -1. },
				{ 1., 1., -i, i },
				{ -i, i, 1., 1. }
			};
			B1m *= 0.5;
			B1.setOperator(B1m);

			MatrixClass B2m{
				{ -1, i, 1., i },
				{ 1., i, 1., -i },
				{ 1., -i, 1., i },
				{ -1., -i, 1., -i }
			};
			B2m *= 0.5;
			B2.setOperator(B2m);

			MatrixClass B3m{
				{ 1., 0., 0., 1. },
				{ -1., 0., 0., 1. },
				{ 0., 1., 1., 0. },
				{ 0., 1., -1., 0. }
			};
			B3m *= invSqrt2;
			B3.setOperator(B3m);

			BaseClass::reg.setToBasisState(0);
		}

		unsigned int Execute() override
		{
			unsigned int nrWins = 0;

			for (unsigned int i = 0; i < nrPlays; ++i)
			{
				unsigned int resAlice, resBob;
				if (Play(distInd(rng), distInd(rng), resAlice, resBob))
					++nrWins;
			}

			return nrWins;
		}

		// a = the number of the row for Alice, b = the number of the column for Bob
		// resAlice = the row from Alice, bit encoded, resBob = the column from Bob, bit encoded
		// returns true if Alice and Bob won, false otherwise
		bool Play(int a, int b, unsigned int& resAlice, unsigned int& resBob)
		{
			assert(a >= 1 && a <= 3 && b >= 0 && b <= 2);

			Init();

			// the Alice & Bob order does not matter, neither for applying the gates nor for measurements (for example the mesurements could be moved all at the end with the same result)
			if (dist_bool(rng))
			{
				resAlice = AliceOperations(a);
				resBob = BobOperations(b);
			}
			else
			{
				resBob = BobOperations(b);
				resAlice = AliceOperations(a);
			}

			// check the intersection, Alice and Bob need to have the same value for the common matrix element
			unsigned int alice = resAlice;
			alice >>= 3 - b;
			alice &= 1;

			unsigned int bob = resBob; 
			bob >>= 3 - a;
			bob &= 1;

			return alice == bob;
		}

		unsigned int getNrPlays() const
		{
			return nrPlays;
		}

		void setNrPlays(unsigned int n)
		{
			nrPlays = n;
		}

	protected:
		void Init()
		{
			BaseClass::reg.setToBasisState(0);
			// put it in the state 0.5 * (|0011> + |1100> - |1001> - |0110>)
			BaseClass::ApplyGate(h, 0);
			BaseClass::ApplyGate(h, 1);
			BaseClass::ApplyGate(z, 0);
			BaseClass::ApplyGate(z, 1);

			BaseClass::ApplyGate(cnot, 2, 0);
			BaseClass::ApplyGate(cnot, 3, 1);
			BaseClass::ApplyGate(x, 2);
			BaseClass::ApplyGate(x, 3);
		}

		unsigned int AliceOperations(int a)
		{
			switch (a)
			{
			case 1:
				BaseClass::ApplyGate(A1, 0, 1);
				break;
			case 2:
				BaseClass::ApplyGate(A2, 0, 1);
				break;
			case 3:
				BaseClass::ApplyGate(A3, 0, 1);
				break;
			}

			// Alice measures
			unsigned int resAlice = BaseClass::Measure(0, 1) << 1;
			// complete the result, Alice needs to have the row with an even number of ones
			if (resAlice == 2 || resAlice == 4)
				resAlice |= 1;

			return resAlice;
		}

		unsigned int BobOperations(int b)
		{
			switch (b)
			{
			case 1:
				BaseClass::ApplyGate(B1, 2, 3);
				break;
			case 2:
				BaseClass::ApplyGate(B2, 2, 3);
				break;
			case 3:
				BaseClass::ApplyGate(B3, 2, 3);
				break;
			}

			// Bob measures
			unsigned int resBob = BaseClass::Measure(2, 3) << 1;
			// complete the result, Bob needs to have the column with an odd number of ones
			if (resBob == 0 || resBob == 6)
				resBob |= 1;

			return resBob;
		}

		std::mt19937_64 rng;
		std::uniform_int_distribution<> distInd;
		std::bernoulli_distribution dist_bool;

		// gates for creating the entangled state
		QC::Gates::HadamardGate<MatrixClass> h;
		QC::Gates::PauliZGate<MatrixClass> z;
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::CNOTGate<MatrixClass> cnot; // controlled x


		// Alice gates
		QC::Gates::TwoQubitsGate<MatrixClass> A1;
		QC::Gates::TwoQubitsGate<MatrixClass> A2;
		QC::Gates::TwoQubitsGate<MatrixClass> A3;

		// Bob gates
		QC::Gates::TwoQubitsGate<MatrixClass> B1;
		QC::Gates::TwoQubitsGate<MatrixClass> B2;
		QC::Gates::TwoQubitsGate<MatrixClass> B3;

		unsigned int nrPlays;
	};

}


