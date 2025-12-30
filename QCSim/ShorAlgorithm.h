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
		FxBase(size_t L, size_t M, size_t C)
			: fRegisterStartQubit(L), Number(C), A(2), fRegisterNrQubits(M)
		{
		}

		void setParam(size_t a)
		{
			A = a;
		}

		size_t getParam() const
		{
			return A;
		}

		size_t mod(unsigned long long int v)
		{
			return v % Number;
		}

		size_t getFunctionStartQubit() const
		{
			return fRegisterStartQubit;
		}

		size_t getFunctionNrQubits() const
		{
			return fRegisterNrQubits;
		}

	protected:
		size_t fRegisterStartQubit;
		size_t Number;
		size_t A;
		size_t fRegisterNrQubits;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class Fx : public FxBase<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = FxBase<VectorClass, MatrixClass>;

		Fx(size_t L, size_t M, size_t C)
			: BaseClass(L, M, C)
		{
		}

		void Apply(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			// for each l qubit from x (0 - L-1 range)
			// construct a controlled gate by the qubit 
			const size_t BasisStatesNo = reg.getNrBasisStates();
			const size_t xmask = (1ULL << BaseClass::fRegisterStartQubit) - 1;
			const size_t fmask = ~xmask;


			unsigned long long int An = BaseClass::A;
			size_t lbit = 1;
			for (size_t l = 0; l < BaseClass::fRegisterStartQubit; ++l)
			{
				MatrixClass gateOperator = MatrixClass::Zero(BasisStatesNo, BasisStatesNo);

				for (size_t k = 0; k < BasisStatesNo; ++k)
				{
					const size_t xbits = k & xmask;
					if (xbits & lbit)
					{
						size_t f = (k & fmask) >> BaseClass::fRegisterStartQubit;

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

		ShorAlgorithmBase(size_t C = 15, size_t N = 7, size_t L = 3, size_t M = 4, unsigned int addseed = 0)
			: BaseClass(N, addseed),
			Number(C), fx(L, M, C)
		{
		}

		// returns false is Shor algorithm was not used, otherwise true 
		bool factorize(size_t& p1, size_t& p2, size_t numAttempts = 10)
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

			const size_t sq = static_cast<size_t>(round(sqrt(Number)));
			if (sq * sq == Number)
			{
				p1 = p2 = sq;
				return false;
			}


			std::random_device rdev;
			std::mt19937 genr(rdev());
			std::uniform_int_distribution<> dist(2ULL, static_cast<unsigned long long>(Number) - 1ULL);

			// well, this probably should be more optimized, but it seems to work
			for (size_t t = 0; t < numAttempts; ++t)
			{
				const int a = dist(genr);

				const size_t g = gcd(static_cast<int>(Number), a);
				if (g > 1)
				{
					//continue;
					p1 = g;
					p2 = Number / g;
					return false;
				}

				setA(a);

				// period finding

				const size_t xmask = (1ULL << fx.getFunctionStartQubit()) - 1;

				// use a single measurement to guess the period (but not zero)

				size_t state = static_cast<Derived*>(this)->Execute() & xmask;
				while (!state) state = static_cast<Derived*>(this)->Execute() & xmask;

				const int d = static_cast<int>(pow(2, fx.getFunctionStartQubit()));
				const double val = static_cast<double>(state) / d;

				// p is guessed from here using continued fractions
				std::vector<int> nums;
				std::vector<int> divs;
				getFractions(continuedFraction(val, static_cast<int>(fx.getFunctionStartQubit())), nums, divs);

				for (size_t i = 1; i < divs.size(); ++i) // skip first as it's for the integer part
				{
					int p = divs[i];
					if (p % 2) p *= 2;
					while (p < Number)
					{
						if (fx.mod(static_cast<size_t>(pow(fx.getParam(), p))) == 1)
						{
							const int v = static_cast<int>(pow(fx.getParam(), p / 2));
							const int m = static_cast<int>(fx.mod(v));
							if (m != 1 && m != Number - 1)
							{
								p1 = gcd(m - 1, static_cast<int>(Number));
								p2 = gcd(m + 1, static_cast<int>(Number));

								// either both or at least one are factors
								const size_t t1 = Number / p1;
								const size_t t2 = Number / p2;
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

		void setA(size_t a)
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

		size_t Number;
		ShorFunction fx;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class ShorFunction = Fx<VectorClass, MatrixClass>> class ShorAlgorithm 
		: public ShorAlgorithmBase<ShorAlgorithm<VectorClass, MatrixClass, ShorFunction>, VectorClass, MatrixClass, ShorFunction >
	{
	public:
		using BaseClass = ShorAlgorithmBase<ShorAlgorithm<VectorClass, MatrixClass, ShorFunction>, VectorClass, MatrixClass>;

		ShorAlgorithm(size_t C = 15, size_t N = 7, size_t L = 3, unsigned int addseed = 0)
			: BaseClass(C, N, L, N - L, addseed),
			phaseEstimation(BaseClass::fx, N, L)
		{
		}

		size_t Execute() override
		{
			return phaseEstimation.Execute(BaseClass::BaseClass::reg);
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
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
		
		FxWithoutTensorProduct(size_t L, size_t M, size_t C)
			: BaseClass(L, M, C)
		{
			std::vector<size_t> controlQubits(L);
			std::iota(controlQubits.begin(), controlQubits.end(), 0);

			nControlledNOTs.SetControlQubits(controlQubits);
			nControlledNOTs.SetStartAncillaQubits(L + M);
		}

		void Apply(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			const size_t nrStates = 1ULL << BaseClass::fRegisterStartQubit;
			const size_t mask = ((1ULL << (BaseClass::fRegisterStartQubit + BaseClass::fRegisterNrQubits)) - 1) >> BaseClass::fRegisterStartQubit;

			std::vector<bool> qubits(BaseClass::fRegisterStartQubit, false);

			size_t An = BaseClass::mod(BaseClass::A);

			for (size_t state = 0; state < nrStates; ++state)
			{
				size_t fval = An;
				if (!fval) continue;

				size_t v = state;
				for (size_t q = 0; q < BaseClass::fRegisterStartQubit; ++q)
				{
					if (qubits[q] != ((v & 1) == 0))
					{
						reg.ApplyGate(x, q);
						qubits[q] = !qubits[q]; // x gate was applied
					}

					v >>= 1;
				}

				nControlledNOTs.ClearGates();

				for (size_t q = BaseClass::fRegisterStartQubit; (fval & mask) != 0 && q < BaseClass::fRegisterStartQubit + BaseClass::fRegisterNrQubits; ++q)
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
			for (size_t q = 0; q < BaseClass::fRegisterStartQubit; ++q)
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

		ShorAlgorithmWithoutTensorProduct(size_t C = 15, size_t L = 3, int M = 4, unsigned int addseed = 0)
			: BaseClass(C, 2 * L + M - 1, L, M, addseed),
			fourier(2 * L + M - 1, 0, L - 1)
		{
		}

		size_t Execute() override
		{
			ExecuteWithoutMeasurement();

			// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
			// the f should be measured separately to find out its content

			//return BaseClass::BaseClass::Measure(0, BaseClass::getFunctionStartQubit() - 1);
			return BaseClass::BaseClass::Measure();
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
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
			for (size_t q = 0; q < BaseClass::fx.getFunctionStartQubit(); ++q)
				BaseClass::BaseClass::ApplyGate(fourier.hadamard, q);
		}

		void IQFT(bool swap = true)
		{
			fourier.IQFT(BaseClass::BaseClass::reg, swap);
		}

		QC::SubAlgo::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};
}
