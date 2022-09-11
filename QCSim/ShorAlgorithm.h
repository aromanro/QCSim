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

			// TODO: FIX IT!!!!!
			// don't forget to uncomment the inverse fourier when period finding works

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

				/*
				// check to see if in every column there is a one (and only one):
				
				for (int z = 0; z < BasisStatesNo; ++z)
				{
					int cnt = 0;
					for (int l = 0; l < BasisStatesNo; ++l)
					{
						if (abs(gateOperator(l, z).real()) > 0.001)
							++cnt;
					}

					if (cnt == 0)
					{
						std::cout << "Ooops, I found a column with all zeros: " << z << std::endl;
						exit(0);
					}
					else if (cnt > 1)
					{
						std::cout << "Ooops, I found a column with more than one non-zero: " << cnt << " in column: " << z << std::endl;

						for (int l = 0; l < BasisStatesNo; ++l)
						{
							if (abs(gateOperator(l, z).real()) > 0.001)
							{
								std::cout << "Line: " << l << " Value: " << gateOperator(l, z) << std::endl;
							}
						}

						exit(0);
					}
				}
				
				// check to see if in every line there is a one (and only one), otherwise (together with the above condition) it's not a permutation matrix:
				
				for (int z = 0; z < BasisStatesNo; ++z)
				{
					int cnt = 0;
					for (int l = 0; l < BasisStatesNo; ++l)
					{
						if (abs(gateOperator(z, l).real()) > 0.001)
							++cnt;
					}

					if (cnt == 0)
					{
						std::cout << "Ooops, I found a line with all zeros: " << z << std::endl;
						exit(0);
					}
					else if (cnt > 1)
					{
						std::cout << "Ooops, I found a line with more than one non-zero: " << cnt << " in line: " << z << std::endl;

						for (int l = 0; l < BasisStatesNo; ++l)
						{
							if (abs(gateOperator(z, l).real()) > 0.001)
							{
								std::cout << "Column: " << l << " Value: " << gateOperator(z, l) << std::endl;
							}
						}

						std::cout << std::endl << gateOperator.block(0, 0, 10, 10) << std::endl;

						exit(0);
					}
				}

				// check unitarity:
				
				MatrixClass mm = gateOperator.adjoint() * gateOperator;
				for (int z = 0; z < BasisStatesNo; ++z)
				{
					if (abs(mm(z,z).real()-1) > 0.00000001)
					{
						std::cout << "i,j: " << z << std::endl;
						exit(1);
					}

					for (int z2 = 0; z2 < BasisStatesNo; ++z2)
					{
						if (z != z2 && abs(mm(z, z2).real()) > 1E-20)
						{
							std::cout << "i: " << z << " j: " << z2 << " Val: " << mm(z, z2) << std::endl;

							std::cout << mm.block(0,0,10,10) << std::endl;

							std::cout << std::endl << gateOperator.block(0, 0, 10, 10) << std::endl;

							exit(2);
						}
					}
				}
				//exit(0);
				*/
				
				
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
