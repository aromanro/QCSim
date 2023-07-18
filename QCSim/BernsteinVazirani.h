#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"
#include "NControlledNotWithAncilla.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <numeric>

namespace BernsteinVazirani {

	template<class MatrixClass = Eigen::MatrixXcd> class Oracle :
		public QC::Gates::QuantumGate<MatrixClass>
	{
	public:
		void setString(unsigned int str)
		{
			stringFunction = str;
		}

		MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const override
		{
			const unsigned int nrBasisStates = 1u << nrQubits;
			MatrixClass extOperatorMat = MatrixClass::Identity(nrBasisStates, nrBasisStates);

			for (unsigned int x = 0; x < nrBasisStates; ++x)
			{
				if (f(x))
					extOperatorMat(x, x) = -1;
			}

			assert(checkUnitary(extOperatorMat));

			return extOperatorMat;
		}

	protected:
		unsigned int f(unsigned int x) const
		{
			// dot product modulo 2 between x and string
			unsigned int prod = x & stringFunction;

			unsigned int accum = 0;
			while (prod)
			{
				accum += prod & 1;
				prod >>= 1;
			}

			return accum % 2;
		}

		unsigned int stringFunction = 0;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BernsteinVaziraniAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		BernsteinVaziraniAlgorithm(unsigned int N = 3, int addseed = 0)
			: BaseClass(N, addseed)
		{
			assert(N >= 1);

			setString(0); // prevent issues if the string is not set before execution
		}

		void setString(unsigned int str)
		{
			Oracle<MatrixClass> o;
			o.setString(str);
			OracleOp = o.getOperatorMatrix(BaseClass::getNrQubits());
		}

		unsigned int Execute() override
		{
			Init();

			BaseClass::ApplyOperatorMatrix(OracleOp);
			ApplyHadamardOnAllQubits();

			return BaseClass::Measure();
		}

	protected:
		void Init()
		{
			//BaseClass::setToBasisState(0);
			//ApplyHadamardOnAllQubits();
			BaseClass::setToEqualSuperposition(); // the same thing as commented above
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < BaseClass::getNrQubits(); ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		MatrixClass OracleOp;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};


	// as above, but using gates instead of the big matrix for the oracle
	// needs ancilla qubits
	// There is one ancilla qubit used by the algorithm
	// the n controlled not needs N - 2 ancilla qubits 
	// that's how many additional ancilla qubits are needed to reduce the control qubits down to 2 for a ccnot (or 1 for a cnot)

	// the algorithm is similar with Deutsch-Jozsa, but the oracle is different


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BernsteinVaziraniAlgorithmWithGatesOracle :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		BernsteinVaziraniAlgorithmWithGatesOracle(unsigned int N = 3, int addseed = 0)
			: BaseClass(2 * N - 1, addseed),
			nControlledNOT(INT_MAX),
			stringFunction(0)
		{
			assert(N >= 1);

			std::vector<unsigned int> controlQubits(N);
			std::iota(controlQubits.begin(), controlQubits.end(), 0);

			nControlledNOT.SetControlQubits(controlQubits);
			nControlledNOT.SetTargetQubit(N);
			nControlledNOT.SetStartAncillaQubits(N + 1);
		}

		void setString(unsigned int str)
		{
			stringFunction = str;
		}

		unsigned int Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubits();

			return BaseClass::Measure(0, getAlgoQubits() - 1);
		}

		unsigned int getAlgoQubits() const
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();

			return (nrQubits + 1) / 2;
		}

	protected:
		void Init()
		{
			BaseClass::setToQubitState(getAlgoQubits());
			// or like this:
			//BaseClass::setToBasisState(0);
			//BaseClass::ApplyGate(x, getAlgoQubits());
		}

		// not on the ancilla ones used for the n controlled not, though
		void ApplyHadamardOnAllQubits()
		{
			// including on the ancilla one
			for (unsigned int q = 0; q <= getAlgoQubits(); ++q)
				BaseClass::ApplyGate(hadamard, q);
		}


		void ApplyOracle()
		{
			const unsigned int nrQubits = getAlgoQubits();
			const unsigned int nrBasisStates = 1u << nrQubits;

			/*
			for (unsigned int state = 0; state < nrBasisStates; ++state)
			{
				unsigned int prod = state & stringFunction;
				unsigned int accum = 0;
				while (prod)
				{
					accum += prod & 1;
					prod >>= 1;
				}
				if (accum % 2 == 0) continue;

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

			for (unsigned int state = 0; state < nrBasisStates; ++state)
			{
				unsigned int prod = state & stringFunction;
				unsigned int accum = 0;
				while (prod)
				{
					accum += prod & 1;
					prod >>= 1;
				}
				if (accum % 2 == 0) continue;
			
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

				nControlledNOT.Execute(BaseClass::reg);
			}

			// undo the x gates
			for (unsigned int q = 0; q < nrQubits; ++q)
				if (qubits[q])
					BaseClass::ApplyGate(x, q);					
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::SubAlgo::NControlledNotWithAncilla<VectorClass, MatrixClass> nControlledNOT;
		QC::Gates::PauliXGate<MatrixClass> x;
		unsigned int stringFunction;
	};
}
