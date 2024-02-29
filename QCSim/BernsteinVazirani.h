#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"
#include "Oracle.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace BernsteinVazirani {

	template<class MatrixClass = Eigen::MatrixXcd> class Oracle :
		public QC::Gates::QuantumGate<MatrixClass>
	{
	public:
		void setString(size_t str)
		{
			stringFunction = str;
		}

		MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
		{
			const size_t nrBasisStates = 1u << nrQubits;
			MatrixClass extOperatorMat = MatrixClass::Identity(nrBasisStates, nrBasisStates);

			for (size_t x = 0; x < nrBasisStates; ++x)
			{
				if (f(x))
					extOperatorMat(x, x) = -1;
			}

			assert(checkUnitary(extOperatorMat));

			return extOperatorMat;
		}

	protected:
		size_t f(size_t x) const
		{
			// dot product modulo 2 between x and string
			size_t prod = x & stringFunction;

			size_t accum = 0;
			while (prod)
			{
				accum += prod & 1;
				prod >>= 1;
			}

			return accum % 2;
		}

		size_t stringFunction = 0;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BernsteinVaziraniAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		BernsteinVaziraniAlgorithm(size_t N = 3, int addseed = 0)
			: BaseClass(N, addseed)
		{
			assert(N >= 1);

			setString(0); // prevent issues if the string is not set before execution
		}

		void setString(size_t str)
		{
			Oracle<MatrixClass> o;
			o.setString(str);
			OracleOp = o.getOperatorMatrix(BaseClass::getNrQubits());
		}

		size_t Execute() override
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
			for (size_t i = 0; i < BaseClass::getNrQubits(); ++i)
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

	class BernsteinVaziraniFunction
	{
	public:
		BernsteinVaziraniFunction(size_t str = 0)
			: stringFunction(str)
		{
		}

		bool operator()(size_t state) const
		{
			size_t prod = state & stringFunction;
			size_t accum = 0;
			while (prod)
			{
				accum += prod & 1;
				prod >>= 1;
			}
			if (accum % 2 == 0)
				return false;

			return true;
		}

	protected:
		size_t stringFunction;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BernsteinVaziraniAlgorithmWithGatesOracle :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		BernsteinVaziraniAlgorithmWithGatesOracle(size_t N = 3, int addseed = 0)
			: BaseClass(2 * N - 1, addseed),
			oracle(2 * N - 1, 0, N - 1, N, N + 1)
		{
			assert(N >= 1);
		}

		void setString(size_t str)
		{
			BernsteinVaziraniFunction func(str);
			oracle.setFunction(func);
		}

		size_t Execute() override
		{
			Init();

			ApplyHadamardOnAllQubits();
			ApplyOracle();
			ApplyHadamardOnAllQubits();

			return BaseClass::Measure(0, getAlgoQubits() - 1);
		}

		size_t getAlgoQubits() const
		{
			const size_t nrQubits = BaseClass::getNrQubits();

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
			for (size_t q = 0; q <= getAlgoQubits(); ++q)
				BaseClass::ApplyGate(hadamard, q);
		}


		void ApplyOracle()
		{
			oracle.Execute(BaseClass::reg);
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Oracles::OracleWithGatesSimpleFunction<BernsteinVaziraniFunction, VectorClass, MatrixClass> oracle;
	};
}
