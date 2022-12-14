#pragma once

#include "Function.h"
#include "QuantumFourierTransform.h"

namespace Shor {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class Fx : public QC::Function<VectorClass, MatrixClass>
	{
	public:
		Fx(unsigned int L, unsigned int C)
			: fRegisterStartQubit(L), Number(C), A(2)
		{
		}

		void Apply(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			// for each l qubit from x (0 - L-1 range)
			// construct a controlled gate by the qubit 
			const unsigned int BasisStatesNo = reg.getNrBasisStates();
			const unsigned int xmask = (1 << fRegisterStartQubit) - 1;
			const unsigned int fmask = ~xmask;


			unsigned long long int An = A;
			for (unsigned int l = 0; l < fRegisterStartQubit; ++l)
			{
				MatrixClass gateOperator = MatrixClass::Zero(BasisStatesNo, BasisStatesNo);

				const unsigned int lbit = 1u << l;

				for (unsigned int k = 0; k < BasisStatesNo; ++k)
				{
					const unsigned int xbits = k & xmask;
					if (xbits & lbit)
					{
						unsigned int f = (k & fmask) >> fRegisterStartQubit;

						if (f >= Number)
							gateOperator(k, k) = 1;
						else
						{
							f = mod(mod(An) * f);
							f <<= fRegisterStartQubit;

							gateOperator(f | xbits, k) = 1;
						}
					}
					else
						gateOperator(k, k) = 1;
				}

				// apply it
				reg.ApplyOperatorMatrix(gateOperator);

				An *= An;
			}
		}

		void setParam(unsigned int a)
		{
			A = a;
		}

		unsigned int mod(unsigned long long int v)
		{
			return v % Number;
		}

	protected:
		unsigned int fRegisterStartQubit;
		unsigned int Number;
		unsigned int A;
	};



	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorAlgorithm : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		ShorAlgorithm(unsigned int C = 15, unsigned int N = 7, unsigned int L = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed), 
			fourier(N, 0, L - 1), Number(C), fRegisterStartQubit(L), A(2), fx(L, C)
		{
		}

		unsigned int Execute() override
		{
			Init();

			// apply hadamard over each qubit from the x-register
			// reuse the hadamard gate from the fourier transform base class
			for (unsigned int i = 0; i < fRegisterStartQubit; ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(fourier.hadamard, i);

			// now the f(x)
			fx.Apply(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);

			// it doesn't really matter if you measure the qubits from f and when you do after the above
			// or if you measure them several times in a row
			//QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(fRegisterStartQubit, QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);

			// then perform an inverse fourier transform
			fourier.IQFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);

			// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
			// the f should be measured separately to find out its content

			//return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(0, fRegisterStartQubit - 1);
			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}

		void setA(unsigned int a)
		{
			A = a;
			fx.setParam(A);
		}

		// returns false is Shor algorithm was not used, otherwise true 
		bool factorize(unsigned int& p1, unsigned int& p2, unsigned int numAttempts = 10)
		{
			if (Number > 22)
			{
				// don't allow too large numbers
				p1 = p2 = 0;
				return false;
			}

			if (Number % 2 == 0)
			{
				p1 = 2;
				p2 = Number / 2;
				return false;
			}

			const unsigned int sq = static_cast<unsigned int>(round(sqrt(Number)));
			if (sq * sq == Number)
			{
				p1 = p2 = sq;
				return false;
			}

			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> dist(2, Number - 1);

			// well, this probably should be more optimized, but it seems to work
			for (unsigned int t = 0; t < numAttempts; ++t)
			{
				const int a = dist(gen);

				const unsigned int g = gcd(Number, a);
				if (g > 1)
				{
					//continue;
					p1 = g;
					p2 = Number / g;
					return false;
				}

				setA(a);

				// period finding

				const unsigned int BasisStatesNo = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates();
				const unsigned int xmask = (1 << fRegisterStartQubit) - 1;

				// use a single measurement to guess the period (but not zero)

				unsigned int state = Execute() & xmask;
				while (!state) state = Execute() & xmask;

				const int d = static_cast<int>(pow(2, fRegisterStartQubit));
				const double val = static_cast<double>(state) / d;

				// p is guessed from here using continued fractions
				std::vector<int> nums;
				std::vector<int> divs;
				getFractions(continuedFraction(val, fRegisterStartQubit), nums, divs);

				for (size_t i = 1; i < divs.size(); ++i) // skip first as it's for the integer part
				{
					unsigned int p = divs[i];
					if (p % 2) p *= 2;
					while (p < Number)
					{
						if (fx.mod(static_cast<unsigned int>(pow(A, p))) == 1)
						{
							const unsigned int v = static_cast<unsigned int>(pow(A, p / 2));
							const unsigned int m = fx.mod(v);
							if (m != 1 && m != Number - 1)
							{
								p1 = gcd(m - 1, Number);
								p2 = gcd(m + 1, Number);

								// either both or at least one are factors
								const unsigned int t1 = Number / p1;
								const unsigned int t2 = Number / p2;
								if (p1 * t1 == Number)
									p2 = t1;
								else
									p1 = t2;

								//if (!((p1 == 1 && p2 == Number) || (p1 == Number && p2 == 1)))
								return true;
							}
						}

						p *= 2;
					}
				}
			}

			p1 = 1;
			p2 = Number;

			return false;
		}

	protected:
		void Init()
		{
			const unsigned int state = 1 << fRegisterStartQubit;
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(state);
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

		static std::vector<int> continuedFraction(double val, int limitIter, double limitPrecision = 1E-5)
		{
			std::vector<int> res;

			int intPart = static_cast<int>(val);
			while (limitIter--)
			{
				res.push_back(intPart);
				val -= intPart;
				if (abs(val) < limitPrecision) break;
				val = 1. / val;
				intPart = static_cast<int>(val);
			}

			return res;
		}

		static void getFractions(const std::vector<int> contFrac, std::vector<int>& nums, std::vector<int>& denoms)
		{
			size_t sz = contFrac.size();
			nums.resize(sz);
			denoms.resize(sz);

			int n2 = 0;
			int n1 = 1;
			int d2 = 1;
			int d1 = 0;

			for (size_t i = 0; i < sz; ++i)
			{
				const int n = contFrac[i] * n1 + n2;
				const int d = contFrac[i] * d1 + d2;

				nums[i] = n;
				denoms[i] = d;

				n2 = n1;
				d2 = d1;
				n1 = n;
				d1 = d;
			}
		}

		unsigned int getNBits() const
		{
			return  QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - fRegisterStartQubit;
		}

		QC::QuantumFourierTransform<VectorClass, MatrixClass> fourier;

		unsigned int Number;
		unsigned int fRegisterStartQubit;
		unsigned int A;

		Fx<VectorClass, MatrixClass> fx;
	};
}
