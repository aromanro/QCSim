#pragma once

#include "Function.h"
#include "PhaseEstimation.h"

namespace Shor {

	// this for now follows the implementation described in "Undergraduate computational physics projects on quantum computing" by D. Candela
	// https://pubs.aip.org/aapt/ajp/article-abstract/83/8/688/235341/Undergraduate-computational-physics-projects-on
	// there is a way to make it better, allowing using fewer qubits for a given number, see the description in Lesson 8 from 
	// "Fundamentals In Quantum Algorithms: A Tutorial Series Using Qiskit Continued" by Daniel Koch, Saahil Patel, Laura Wessing, Paul M. Alsing
	// https://arxiv.org/abs/2008.10647


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class FxBase : public QC::Function<VectorClass, MatrixClass>
	{
	public:
		FxBase(unsigned int L, unsigned int M, unsigned int C)
			: fRegisterStartQubit(L), Number(C), A(2), fRegisterNrQubits(M)
		{
		}

		void setParam(unsigned int a)
		{
			A = a;
		}

		unsigned int getParam() const
		{
			return A;
		}

		unsigned int mod(unsigned long long int v)
		{
			return v % Number;
		}

		unsigned int getFunctionStartQubit() const
		{
			return fRegisterStartQubit;
		}

		unsigned int getFunctionNrQubits() const
		{
			return fRegisterNrQubits;
		}

	protected:
		unsigned int fRegisterStartQubit;
		unsigned int Number;
		unsigned int A;
		unsigned int fRegisterNrQubits;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class Fx : public FxBase<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = FxBase<VectorClass, MatrixClass>;

		Fx(unsigned int L, unsigned int M, unsigned int C)
			: BaseClass(L, M, C)
		{
		}

		void Apply(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			// for each l qubit from x (0 - L-1 range)
			// construct a controlled gate by the qubit 
			const unsigned int BasisStatesNo = reg.getNrBasisStates();
			const unsigned int xmask = (1 << BaseClass::fRegisterStartQubit) - 1;
			const unsigned int fmask = ~xmask;


			unsigned long long int An = BaseClass::A;
			unsigned int lbit = 1;
			for (unsigned int l = 0; l < BaseClass::fRegisterStartQubit; ++l)
			{
				MatrixClass gateOperator = MatrixClass::Zero(BasisStatesNo, BasisStatesNo);

				for (unsigned int k = 0; k < BasisStatesNo; ++k)
				{
					const unsigned int xbits = k & xmask;
					if (xbits & lbit)
					{
						unsigned int f = (k & fmask) >> BaseClass::fRegisterStartQubit;

						if (f >= BaseClass::Number)
							gateOperator(k, k) = 1;
						else
						{
							f = BaseClass::mod(BaseClass::mod(An) * f);
							f <<= BaseClass::fRegisterStartQubit;

							gateOperator(f | xbits, k) = 1;
						}
					}
					else
						gateOperator(k, k) = 1;
				}

				assert(checkUnitary(gateOperator));

				// apply it
				reg.ApplyOperatorMatrix(gateOperator);

				An *= An;
				lbit <<= 1;
			}
		}
	};


	template<class Derived, class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class ShorFunction = Fx<VectorClass, MatrixClass>> class ShorAlgorithmBase : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		ShorAlgorithmBase(unsigned int C = 15, unsigned int N = 7, unsigned int L = 3, unsigned int M = 4, int addseed = 0)
			: BaseClass(N, addseed),
			Number(C), fx(L, M, C)
		{
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


			std::random_device rdev;
			std::mt19937 genr(rdev());
			std::uniform_int_distribution<> dist(2, Number - 1);

			// well, this probably should be more optimized, but it seems to work
			for (unsigned int t = 0; t < numAttempts; ++t)
			{
				const int a = dist(genr);

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

				const unsigned int xmask = (1 << fx.getFunctionStartQubit()) - 1;

				// use a single measurement to guess the period (but not zero)

				unsigned int state = static_cast<Derived*>(this)->Execute() & xmask;
				while (!state) state = static_cast<Derived*>(this)->Execute() & xmask;

				const int d = static_cast<int>(pow(2, fx.getFunctionStartQubit()));
				const double val = static_cast<double>(state) / d;

				// p is guessed from here using continued fractions
				std::vector<int> nums;
				std::vector<int> divs;
				getFractions(continuedFraction(val, fx.getFunctionStartQubit()), nums, divs);

				for (size_t i = 1; i < divs.size(); ++i) // skip first as it's for the integer part
				{
					unsigned int p = divs[i];
					if (p % 2) p *= 2;
					while (p < Number)
					{
						if (fx.mod(static_cast<unsigned int>(pow(fx.getParam(), p))) == 1)
						{
							const unsigned int v = static_cast<unsigned int>(pow(fx.getParam(), p / 2));
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

		void setA(unsigned int a)
		{
			fx.setParam(a);
		}

	protected:
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

		unsigned int Number;
		ShorFunction fx;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class ShorFunction = Fx<VectorClass, MatrixClass>> class ShorAlgorithm 
		: public ShorAlgorithmBase<ShorAlgorithm<VectorClass, MatrixClass, ShorFunction>, VectorClass, MatrixClass, ShorFunction >
	{
	public:
		using BaseClass = ShorAlgorithmBase<ShorAlgorithm<VectorClass, MatrixClass, ShorFunction>, VectorClass, MatrixClass>;

		ShorAlgorithm(unsigned int C = 15, unsigned int N = 7, unsigned int L = 3, int addseed = 0)
			: BaseClass(C, N, L, N - L, addseed),
			phaseEstimation(BaseClass::fx, N, L)
		{
		}

		unsigned int Execute() override
		{
			return phaseEstimation.Execute(BaseClass::BaseClass::reg);
		}

		std::map<unsigned int, unsigned int> ExecuteWithMultipleMeasurements(unsigned int nrMeasurements = 10000)
		{
			return phaseEstimation.ExecuteWithMultipleMeasurements(BaseClass::BaseClass::reg, nrMeasurements);
		}

	protected:
		QC::SubAlgo::ShorPhaseEstimation<VectorClass, MatrixClass> phaseEstimation;
	};


	// this is not the best implementation, for a better one see "Circuit for Shor's algorithm using 2n+3 qubits" https://arxiv.org/abs/quant-ph/0205095
	// it's implemented using the same method as described in 
	// "Fundamentals In Quantum Algorithms: A Tutorial Series Using Qiskit Continued" by Daniel Koch, Saahil Patel, Laura Wessing, Paul M. Alsing
	// https://arxiv.org/abs/2008.10647
	
	// I could also have this implemented using the oracle class, but I'm lazy... I would still have to provide a function implementation, so I put all of it here
	// maybe I'll change it later, might need some changes in the oracle class for it

	// there is quite a bit of common code with the implementation using the tensor product, I might more it in base classes later

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class FxWithoutTensorProduct : public FxBase<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = FxBase<VectorClass, MatrixClass>;
		
		FxWithoutTensorProduct(unsigned int L, unsigned int M, unsigned int C)
			: BaseClass(L, M, C)
		{
			std::vector<unsigned int> controlQubits(L);
			std::iota(controlQubits.begin(), controlQubits.end(), 0);

			nControlledNOTs.SetControlQubits(controlQubits);
			nControlledNOTs.SetStartAncillaQubits(L + M);
		}

		void Apply(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			const unsigned int nrStates = 1 << BaseClass::fRegisterStartQubit;
			const unsigned int mask = ((1 << (BaseClass::fRegisterStartQubit + BaseClass::fRegisterNrQubits)) - 1) >> BaseClass::fRegisterStartQubit;

			std::vector<bool> qubits(BaseClass::fRegisterStartQubit, false);

			unsigned int An = BaseClass::mod(BaseClass::A);

			for (unsigned int state = 0; state < nrStates; ++state)
			{
				unsigned int fval = An;
				if (!fval) continue;

				unsigned int v = state;
				for (unsigned int q = 0; q < BaseClass::fRegisterStartQubit; ++q)
				{
					if (qubits[q] != ((v & 1) == 0))
					{
						reg.ApplyGate(x, q);
						qubits[q] = !qubits[q]; // x gate was applied
					}

					v >>= 1;
				}

				nControlledNOTs.ClearGates();

				for (unsigned int q = BaseClass::fRegisterStartQubit; (fval & mask) != 0 && q < BaseClass::fRegisterStartQubit + BaseClass::fRegisterNrQubits; ++q)
				{
					if (fval & 1)
					{
						QC::Gates::AppliedGate<MatrixClass> notGate(x.getRawOperatorMatrix(), q);
						nControlledNOTs.AddGate(notGate);
					}

					fval >>= 1;
				}

				nControlledNOTs.Execute(reg);

				An *= BaseClass::A;
				An = BaseClass::mod(An);
			}

			// undo the x gates
			for (unsigned int q = 0; q < BaseClass::fRegisterStartQubit; ++q)
				if (qubits[q])
					reg.ApplyGate(x, q);
		}

	protected:
		QC::SubAlgo::NControlledGatesWithAncilla<VectorClass, MatrixClass> nControlledNOTs;
		QC::Gates::PauliXGate<MatrixClass> x;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class ShorFunction = FxWithoutTensorProduct<VectorClass, MatrixClass>> class ShorAlgorithmWithoutTensorProduct 
		: public ShorAlgorithmBase<ShorAlgorithmWithoutTensorProduct<VectorClass, MatrixClass, ShorFunction>, VectorClass, MatrixClass, ShorFunction>
	{
	public:
		using BaseClass = ShorAlgorithmBase<ShorAlgorithmWithoutTensorProduct<VectorClass, MatrixClass, ShorFunction>, VectorClass, MatrixClass, ShorFunction>;

		ShorAlgorithmWithoutTensorProduct(unsigned int C = 15, unsigned int L = 3, int M = 4, int addseed = 0)
			: BaseClass(C, 2 * L + M - 1, L, M, addseed),
			fourier(2 * L + M - 1, 0, L - 1)
		{
		}

		unsigned int Execute() override
		{
			ExecuteWithoutMeasurement();

			// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
			// the f should be measured separately to find out its content

			//return BaseClass::BaseClass::Measure(0, BaseClass::getFunctionStartQubit() - 1);
			return BaseClass::BaseClass::Measure();
		}

		std::map<unsigned int, unsigned int> ExecuteWithMultipleMeasurements(unsigned int nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::BaseClass::RepeatedMeasure(nrMeasurements);
		}

	protected:
		void ExecuteWithoutMeasurement()
		{
			Init();

			BaseClass::fx.Apply(BaseClass::reg);

			IQFT();
		}

		void Init()
		{
			BaseClass::BaseClass::setToBasisState(0);

			// apply hadamard over each qubit from the x-register
			ApplyHadamardOnXRegister();
		}

		void ApplyHadamardOnXRegister()
		{
			// apply hadamard over each qubit from the x-register
			// reuse the hadamard gate from the fourier transform base class
			for (unsigned int q = 0; q < BaseClass::fx.getFunctionStartQubit(); ++q)
				BaseClass::BaseClass::ApplyGate(fourier.hadamard, q);
		}

		void IQFT(bool swap = true)
		{
			fourier.IQFT(BaseClass::BaseClass::reg, swap);
		}

		QC::SubAlgo::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};
}
