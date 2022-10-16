#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace DeutschJozsa {

	// N = 3 reduces to Deutsch's algorithm
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		DeutschJozsaAlgorithm(unsigned int N = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed),
			functionType(FunctionType::constantZero)
		{
		}

		unsigned int Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubitsExceptLast();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}

		bool WasConstantResult(unsigned int state) const
		{
			const unsigned int xmask = (QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates() - 1) >> 1;

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
				const unsigned int funcSize = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates() >> 1;
				fvalues.resize(funcSize);
				const unsigned int nrOnes = funcSize >> 1;
				for (unsigned int i = 0; i < nrOnes; ++i)
					fvalues[i] = true;
				for (unsigned int i = nrOnes; i < funcSize; ++i)
					fvalues[i] = false;

				auto rng = std::default_random_engine{};
				std::shuffle(fvalues.begin(), fvalues.end(), rng);
			}
		}

		FunctionType getFunction() const
		{
			return functionType;
		}

	protected:
		void Init()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToQubitState(QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits(); ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, i);
		}

		void ApplyHadamardOnAllQubitsExceptLast()
		{
			const unsigned int limit = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1;
			for (unsigned int i = 0; i < limit; ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, i);
		}

		void ApplyOracle()
		{
			const unsigned int nrBasisStates = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates();
			MatrixClass U = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			// we have U |y>|x> = |y + f(x)> |x>
			// U is unitary

			// <x|<y| U |y>|x> = <x|x><y|y+f(x)> = <y|y + f(x)>
			// + is modulo 2 addition here

			const unsigned int nrQubits = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits();
			const unsigned int mask = nrBasisStates - 1;
			const unsigned int xmask = mask >> 1;
			const unsigned int ymask = nrBasisStates >> 1;

			for (unsigned int stateBra = 0; stateBra < nrBasisStates; ++stateBra)
				for (unsigned int stateKet = 0; stateKet < nrBasisStates; ++stateKet)
				{
					const unsigned int yval = (stateKet & ymask) ? 1 : 0;

					U(stateBra, stateKet) = (stateBra & ymask) == (f(stateKet & xmask) + yval) % 2;
				}

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyOperatorMatrix(U);
		}

		unsigned int f(unsigned int xstate)
		{
			if (functionType == FunctionType::constantZero) return 0;
			else if (functionType == FunctionType::constantOne) return 1;

			return fvalues[xstate] ? 1 : 0;
		}

		QC::HadamardGate<MatrixClass> hadamard;
		FunctionType functionType;
		std::vector<bool> fvalues; // for balanced only
	};

}