#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "NControlledNotWithAncilla.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <numeric>

namespace DeutschJozsa {

	// N = 2 reduces to Deutsch's algorithm
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DeutschJozsaAlgorithm(unsigned int N = 2, int addseed = 0)
			: BaseClass(N, addseed),
			functionType(FunctionType::constantZero)
		{
			assert(N >= 2);
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
			for (unsigned int q = 0; q < BaseClass::getNrQubits(); ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyHadamardOnAllQubitsExceptLast()
		{
			const unsigned int limit = BaseClass::getNrQubits() - 1;
			for (unsigned int q = 0; q < limit; ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyOracle()
		{
			const unsigned int nrBasisStates = BaseClass::getNrBasisStates();
			MatrixClass U = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			// we have U |y>|x> = |y + f(x)> |x>
			// U is unitary

			// <x1|<y1| U |y2>|x2> = <x1|x2><y1|y2+f(x2)>
			// + is modulo 2 addition

			const unsigned int mask = nrBasisStates - 1;
			const unsigned int xmask = mask >> 1;
			const unsigned int ymask = nrBasisStates >> 1;

			for (unsigned int stateBra = 0; stateBra < nrBasisStates; ++stateBra)
				for (unsigned int stateKet = 0; stateKet < nrBasisStates; ++stateKet)
				{
					const unsigned int xval = (stateBra & xmask);
					const unsigned int yval = (stateBra & ymask) ? 1 : 0;

					// the operator is sumi sumj |i><j|, so line corresponds to ket and column to bra 
					U(stateKet, stateBra) = ((stateKet & ymask) ? 1U : 0U) == (f(xval) + yval) % 2 && xval == (stateKet & xmask);
				}

			assert(checkUnitary(U));

			BaseClass::ApplyOperatorMatrix(U);
		}

		unsigned int f(unsigned int xstate) const
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
	// the n controlled not needs max(0, (N - 1) - 1) ancilla qubits (the first -1 is because of the previously mentioned ancilla) 
	// that's how many additional ancilla qubits are needed to reduce the control qubits down to 2 for a ccnot (or 1 for a cnot)

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithmWithGatesOracle :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DeutschJozsaAlgorithmWithGatesOracle(unsigned int N = 2, int addseed = 0)
			: BaseClass(2 * N - 2, addseed),
			functionType(FunctionType::constantZero),
			nControlledNOT(2 * N - 2)
		{
			assert(N >= 2);

			std::vector<unsigned int> controlQubits(N - 1);
			std::iota(controlQubits.begin(), controlQubits.end(), 0);

			nControlledNOT.SetControlQubits(controlQubits);
			nControlledNOT.SetTargetQubit(N - 1);
			nControlledNOT.SetStartAncillaQubits(N);
		}

		unsigned int Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubitsExceptLast();

			return BaseClass::Measure(0, getAlgoQubits() - 2);
		}

		bool WasConstantResult(unsigned int state) const
		{
			const unsigned int xmask = (getAlgoNrBasisStates() - 1) >> 1;

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
				const unsigned int funcSize = getAlgoNrBasisStates() >> 1;

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

		unsigned int getAlgoQubits() const
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();

			return (nrQubits + 2) / 2;
		}

		unsigned int getAlgoNrBasisStates() const
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

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int q = 0; q < getAlgoQubits(); ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyHadamardOnAllQubitsExceptLast()
		{
			const unsigned int limit = getAlgoQubits() - 1;
			for (unsigned int q = 0; q < limit; ++q)
				BaseClass::ApplyGate(hadamard, q);
		}

		void ApplyOracle()
		{
			// for constant functions the result is easy, so we can skip all those gates for the function
			if (functionType == FunctionType::constantZero) return;
			else if (functionType == FunctionType::constantOne)
				BaseClass::ApplyGate(x, getAlgoQubits() - 1);
			else
			{
				const unsigned int nrQubits = getAlgoQubits() - 1;

				/*
				for (unsigned int state = 0; state < fvalues.size(); ++state)
				{
					if (!fvalues[state]) continue;

					unsigned int v = state;

					for (unsigned int q = 0; q < nrQubits; ++q)
					{
						if ((v & 1) == 0)
							BaseClass::ApplyGate(x, q);

						v >>= 1;
					}

					nControlledNOT.Execute(BaseClass::reg);

					// undo the x gates
					v = state;
					for (unsigned int q = 0; q < nrQubits; ++q)
					{
						if ((v & 1) == 0)
							BaseClass::ApplyGate(x, q);

						v >>= 1;
					}
				}
				*/
				
				// the above code is not that good, I'll let it commented in case this one is too hard to understand
				// this one avoids applying x gates twice on the same qubit, between applying n controlled nots
				
				std::vector<bool> qubits(nrQubits, false);

				for (unsigned int state = 0; state < fvalues.size(); ++state)
				{
					if (!fvalues[state]) continue;

					unsigned int v = state;
					for (unsigned int q = 0; q < nrQubits; ++q)
					{
						if (qubits[q] != ((v & 1) == 0))
						{
							BaseClass::ApplyGate(x, q);
							qubits[q] = !qubits[q]; // x gate was applied
						}

						v >>= 1;
					}

					// the last 'uncompute' - that is, at the end of the last n controlled not - is not really necessary
					// unless the ancilla qubits would be used further, so maybe avoit it like this (don't forget to set it back to true):
#define NOT_CLEAR_ANCILLA_AT_THE_END 1
#ifdef NOT_CLEAR_ANCILLA_AT_THE_END
					if (state == fvalues.size() - 1)
						nControlledNOT.SetClearAncillaAtTheEnd(false);
#endif

					nControlledNOT.Execute(BaseClass::reg);
				}

				// undo the x gates
				for (unsigned int q = 0; q < nrQubits; ++q)
					if (qubits[q])
						BaseClass::ApplyGate(x, q);

#ifdef NOT_CLEAR_ANCILLA_AT_THE_END
				nControlledNOT.SetClearAncillaAtTheEnd(true); // don't let it on false, for the next Execute
#endif
			}
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		FunctionType functionType;
		std::vector<bool> fvalues; // for balanced only
		QC::SubAlgo::NControlledNotWithAncilla<VectorClass, MatrixClass> nControlledNOT;
		QC::Gates::PauliXGate<MatrixClass> x;
	};

}