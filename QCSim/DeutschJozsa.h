#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace DeutschJozsa {

	// N = 3 reduces to Deutsch's algorithm
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		typedef QC::QuantumAlgorithm<VectorClass, MatrixClass> BaseClass;

		DeutschJozsaAlgorithm(unsigned int N = 3, int addseed = 0)
			: BaseClass(N, addseed),
			functionType(FunctionType::constantZero)
		{
		}

		unsigned int Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubitsExceptLast();

			return BaseClass::Measure();
		}

		bool WasConstantResult(unsigned int state) const
		{
			const unsigned int xmask = (BaseClass::getNrBasisStates() - 1) >> 1;

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
				const unsigned int funcSize = BaseClass::getNrBasisStates() >> 1;
				fvalues.resize(funcSize);
				const unsigned int nrOnes = funcSize >> 1;
				for (unsigned int i = 0; i < nrOnes; ++i)
					fvalues[i] = true;
				for (unsigned int i = nrOnes; i < funcSize; ++i)
					fvalues[i] = false;

				const unsigned long long int seed = std::chrono::steady_clock::now().time_since_epoch().count();
				auto rng = std::default_random_engine(static_cast<unsigned int>(seed));
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
			for (unsigned int i = 0; i < BaseClass::getNrQubits(); ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyHadamardOnAllQubitsExceptLast()
		{
			const unsigned int limit = BaseClass::getNrQubits() - 1;
			for (unsigned int i = 0; i < limit; ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyOracle()
		{
			const unsigned int nrBasisStates = BaseClass::getNrBasisStates();
			MatrixClass U = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			// we have U |y>|x> = |y + f(x)> |x>
			// U is unitary

			// <x1|<y1| U |y2>|x2> = <x1|x2><y1|y2+f(x2)>
			// + is modulo 2 addition

			const unsigned int nrQubits = BaseClass::getNrQubits();
			const unsigned int mask = nrBasisStates - 1;
			const unsigned int xmask = mask >> 1;
			const unsigned int ymask = nrBasisStates >> 1;

			for (unsigned int stateBra = 0; stateBra < nrBasisStates; ++stateBra)
				for (unsigned int stateKet = 0; stateKet < nrBasisStates; ++stateKet)
				{
					const unsigned int xval = (stateBra & xmask);
					const unsigned int yval = (stateBra & ymask) ? 1 : 0;

					// the operator is sumi sumj |i><j|, so line corresponds to ket and column to bra 
					U(stateKet, stateBra) = ((stateKet & ymask) ? 1 : 0) == (f(xval) + yval) % 2 && xval == (stateKet & xmask);
				}
			
			assert(checkUnitary(U));

			BaseClass::ApplyOperatorMatrix(U);
		}

		unsigned int f(unsigned int xstate)
		{
			if (functionType == FunctionType::constantZero) return 0;
			else if (functionType == FunctionType::constantOne) return 1;

			return fvalues[xstate] ? 1 : 0;
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		FunctionType functionType;
		std::vector<bool> fvalues; // for balanced only
	};

}