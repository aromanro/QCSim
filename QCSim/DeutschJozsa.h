#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Oracle.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace DeutschJozsa {

	// N = 2 reduces to Deutsch's algorithm
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DeutschJozsaAlgorithm(size_t N = 2, int addseed = 0)
			: BaseClass(N, addseed),
			functionType(FunctionType::constantZero)
		{
			assert(N >= 2);
		}

		size_t Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubitsExceptLast();

			return BaseClass::Measure();
		}

		bool WasConstantResult(size_t state) const
		{
			const size_t xmask = (BaseClass::getNrBasisStates() - 1) >> 1;

			return (state & xmask) == 0;
		}

		enum class FunctionType : unsigned char
		{
			constantZero = 0,
			constantOne = 1,
			balanced
		};

		void setFunction(FunctionType ft)
		{
			functionType = ft;

			if (functionType == FunctionType::balanced)
			{
				const size_t funcSize = BaseClass::getNrBasisStates() >> 1;
				fvalues.resize(funcSize);
				const size_t nrOnes = funcSize >> 1;
				for (size_t i = 0; i < nrOnes; ++i)
					fvalues[i] = true;
				for (size_t i = nrOnes; i < funcSize; ++i)
					fvalues[i] = false;

				const unsigned long long int seed = std::chrono::steady_clock::now().time_since_epoch().count();
				auto rng = std::default_random_engine(static_cast<size_t>(seed));
				std::shuffle(fvalues.begin(), fvalues.end(), rng);
			}
			else fvalues.clear();
		}

		FunctionType getFunction() const
		{
			return functionType;
		}

	protected:
		void Init()
		{
			BaseClass::setToQubitState(BaseClass::getNrQubits() - 1);
		}

		void ApplyHadamardOnAllQubits()
		{
			for (size_t q = 0; q < BaseClass::getNrQubits(); ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyHadamardOnAllQubitsExceptLast()
		{
			const size_t limit = BaseClass::getNrQubits() - 1;
			for (size_t q = 0; q < limit; ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyOracle()
		{
			const size_t nrBasisStates = BaseClass::getNrBasisStates();
			MatrixClass U = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			// we have U |y>|x> = |y + f(x)> |x>
			// U is unitary

			// <x1|<y1| U |y2>|x2> = <x1|x2><y1|y2+f(x2)>
			// + is modulo 2 addition

			const size_t mask = nrBasisStates - 1;
			const size_t xmask = mask >> 1;
			const size_t ymask = nrBasisStates >> 1;

			for (size_t stateBra = 0; stateBra < nrBasisStates; ++stateBra)
				for (size_t stateKet = 0; stateKet < nrBasisStates; ++stateKet)
				{
					const size_t xval = (stateBra & xmask);
					const size_t yval = (stateBra & ymask) ? 1 : 0;

					// the operator is sumi sumj |i><j|, so line corresponds to ket and column to bra 
					U(stateKet, stateBra) = ((stateKet & ymask) ? 1U : 0U) == (f(xval) + yval) % 2 && xval == (stateKet & xmask);
				}

			assert(checkUnitary(U));

			BaseClass::ApplyOperatorMatrix(U);
		}

		size_t f(size_t xstate) const
		{
			if (functionType == FunctionType::constantZero) return 0;
			else if (functionType == FunctionType::constantOne) return 1;

			return fvalues[xstate] ? 1 : 0;
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		FunctionType functionType;
		std::vector<bool> fvalues; // for balanced only
	};



	// the same as above, but instead of generating directly the U operator, it's made out of gates
	// needs ancilla qubits for that

	// the 'function' uses NControlledNotWithAncilla gates to implement the function
	// There is one ancilla qubit used by the algorithm (the y + f(x) one)
	// the n controlled not needs max(0, (N - 1) - 2) ancilla qubits (the first -1 is because of the previously mentioned ancilla) 
	// that's how many additional ancilla qubits are needed to reduce the control qubits down to 2 for a ccnot (or 1 for a cnot)

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithmWithGatesOracle :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DeutschJozsaAlgorithmWithGatesOracle(size_t N = 2, int addseed = 0)
			: BaseClass(N == 2 ? 2 : 2 * N - 3, addseed), oracle(N == 2 ? 2 : 2 * N - 3, 0, N - 2, N - 1, N)
			
		{
			assert(N >= 2);
		}

		size_t Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubitsExceptLast();

			return BaseClass::Measure(0, getAlgoQubits() - 2);
		}

		bool WasConstantResult(size_t state) const
		{
			const size_t xmask = (getAlgoNrBasisStates() - 1) >> 1;

			return (state & xmask) == 0;
		}

		enum class FunctionType : unsigned char
		{
			constantZero = 0,
			constantOne = 1,
			balanced
		};

		
		class DeutschJozsaFunction
		{
		public:
			DeutschJozsaFunction() : functionType(FunctionType::constantZero) {}

			bool operator()(size_t state) const
			{
				if (functionType != FunctionType::balanced) return false;
				
				return fvalues[state];
			}

			void setFunction(FunctionType ft, size_t funcSize)
			{
				functionType = ft;

				if (functionType == FunctionType::balanced)
				{
					fvalues.resize(funcSize);
					const size_t nrOnes = funcSize >> 1;
					for (size_t i = 0; i < nrOnes; ++i)
						fvalues[i] = true;
					for (size_t i = nrOnes; i < funcSize; ++i)
						fvalues[i] = false;

					const unsigned long long int seed = std::chrono::steady_clock::now().time_since_epoch().count();
					auto rng = std::default_random_engine(static_cast<size_t>(seed));
					std::shuffle(fvalues.begin(), fvalues.end(), rng);
				}
				else fvalues.clear();
			}

			FunctionType getFunction() const
			{
				return functionType;
			}

		protected:
			FunctionType functionType;
			std::vector<bool> fvalues; // for balanced only
		};

		void setFunction(FunctionType ft)
		{
			func.setFunction(ft, getAlgoNrBasisStates() >> 1);
			oracle.setFunction(func);
		}

		size_t getAlgoQubits() const
		{
			const size_t nrQubits = BaseClass::getNrQubits();
			if (nrQubits == 2) return 2;

			return (nrQubits + 3) / 2;
		}

		size_t getAlgoNrBasisStates() const
		{
			return 1u << getAlgoQubits();
		}

	protected:
		void Init()
		{
			BaseClass::setToQubitState(getAlgoQubits() - 1);
			// or like this:
			//BaseClass::setToBasisState(0);
			//BaseClass::ApplyGate(x, getAlgoQubits() - 1);
		}

		// not on the ancilla ones used for the n controlled not, though
		void ApplyHadamardOnAllQubits()
		{
			for (size_t q = 0; q < getAlgoQubits(); ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyHadamardOnAllQubitsExceptLast()
		{
			const size_t limit = getAlgoQubits() - 1;
			for (size_t q = 0; q < limit; ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyOracle()
		{
			// for constant functions the result is easy, so we can skip all those gates for the function
			if (func.getFunction() == FunctionType::constantZero) return;
			else if (func.getFunction() == FunctionType::constantOne)
				BaseClass::ApplyGate(x, getAlgoQubits() - 1);
			else
			{
				oracle.Execute(BaseClass::reg);
			}
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		DeutschJozsaFunction func;
		QC::Oracles::OracleWithGatesSimpleFunction<DeutschJozsaFunction, VectorClass, MatrixClass> oracle;
		QC::Gates::PauliXGate<MatrixClass> x;
	};

}