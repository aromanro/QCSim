#pragma once

#include <iostream>

#include "QuantumFourierTransform.h"

namespace Shor {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorAlgorithm : public QC::QuantumFourierTransform<VectorClass, MatrixClass>
	{
	public:
		ShorAlgorithm(unsigned int C = 15, unsigned int N = 7, unsigned int L = 3, int addseed = 0)
			: QC::QuantumFourierTransform<VectorClass, MatrixClass>(N, 0, L - 1, addseed), Number(C), fRegisterStartQubit(L), A(2)
		{
		}

		unsigned int Execute() override
		{
			Init();

			// apply hadamard over each qubit from the x-register
			// reuse the hadamard gate from the fourier transform base class
			for (unsigned int i = 0; i < fRegisterStartQubit; ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(QC::QuantumFourierTransform<VectorClass, MatrixClass>::hadamard, i);

			// now the f(x)

			// for each l qubit from x (0 - L-1 range)
			// construct a controlled gate by the qubit 

			const unsigned int BasisStatesNo = QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrBasisStates();
			const unsigned int xmask = (1 << fRegisterStartQubit) - 1;
			const unsigned int fmask = ~xmask;


			unsigned long long int An = A;
			for (unsigned int l = 0; l < fRegisterStartQubit; ++l)
			{
				MatrixClass gateOperator = MatrixClass::Zero(BasisStatesNo, BasisStatesNo);

				const unsigned int lbit = 1u << l;

				for (unsigned int k = 0; k < BasisStatesNo; ++k)
				{
					if ((k & xmask) == lbit) //ok, this seems to make it a permutation matrix (unitary) but still doesn't work as expected
					{
						unsigned int f = (k & fmask) >> fRegisterStartQubit;

						if (f >= Number)
							gateOperator(k, k) = 1;
						else
						{
							f = mod(mod(An)*f);
							f <<= fRegisterStartQubit;
						    gateOperator(f | lbit, k) = 1;
						}
					}
					else
						gateOperator(k, k) = 1;
				}

				// apply it
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyOperatorMatrix(gateOperator);

				An *= A;
			}

			// then perform an inverse fourier transform
			QC::QuantumFourierTransform<VectorClass, MatrixClass>::IQFT();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.Measure();
		}

		void setA(unsigned int a)
		{
			A = a;
		}



	protected:
		void Init()
		{
			const unsigned int state = 1 << fRegisterStartQubit;
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.setToBasisState(state);
		}

		static int gcd(int a, int b)
		{
			while (b)
			{
				const int t = b;
				b = a % b;
				a = t;
			}

			return a;
		}

		unsigned int mod(unsigned long long int v)
		{
			return v % Number;
		}

		unsigned int Number;
		unsigned int fRegisterStartQubit;
		unsigned int A;
	};
}
