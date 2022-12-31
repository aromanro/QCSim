#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

#define _USE_MATH_DEFINES
#include <math.h>

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

			return extOperatorMat;
		}

	protected:
		unsigned int f(unsigned int x) const
		{
			// dot product modulo 2 between x and string
			unsigned int prod = (x & stringFunction);

			unsigned int accum = 0;
			while (prod)
			{
				accum += (prod & 1);
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
		BernsteinVaziraniAlgorithm(unsigned int N = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed)
		{
			setString(0); // prevent issues if the string is not set before execution
		}

		void setString(unsigned int str)
		{
			Oracle<MatrixClass> o;
			o.setString(str);
			OracleOp = o.getOperatorMatrix(QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits());
		}

		unsigned int Execute() override
		{
			Init();

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyOperatorMatrix(OracleOp);
			ApplyHadamardOnAllQubits();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}

	protected:
		void Init()
		{
			//reg.setToBasisState(0);
			//ApplyHadamardOnAllQubits();
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToEqualSuperposition(); // the same thing as commented above
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits(); ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, i);
		}

		MatrixClass OracleOp;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};
}
