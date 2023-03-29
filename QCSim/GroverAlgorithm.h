#pragma once

#include "QuantumAlgorithm.h"
#include "Oracle.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace Grover {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class GroverAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		GroverAlgorithm(unsigned int N = 3, int addseed = 0)
			: BaseClass(N, addseed)
		{
			J j;
			JOp = j.getOperatorMatrix(BaseClass::getNrQubits());

			setCorrectQuestionState(0);
		}

		void setCorrectQuestionState(unsigned int state)
		{
			Oracle<MatrixClass> o;
			o.setCorrectQuestionState(state);
			OracleOp = o.getOperatorMatrix(BaseClass::getNrQubits());
		}

		unsigned int Execute() override
		{
			Init();

			const unsigned int repeatNo = static_cast<unsigned int>(round(M_PI / 4. * sqrt(BaseClass::getNrBasisStates())));

			for (unsigned int i = 0; i < repeatNo; ++i)
			{
				BaseClass::ApplyOperatorMatrix(OracleOp);
				ApplyDiffusionOperator();
			}

			return BaseClass::Measure();
		}

	protected:
		void Init()
		{
			//BaseClass::reg.setToBasisState(0);
			//ApplyHadamardOnAllQubits();
			BaseClass::setToEqualSuperposition(); // the same thing as commented above
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < BaseClass::getNrQubits(); ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyDiffusionOperator()
		{
			ApplyHadamardOnAllQubits();
			BaseClass::ApplyOperatorMatrix(JOp);
			ApplyHadamardOnAllQubits();
		}

		MatrixClass OracleOp;
		MatrixClass JOp;

		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};

}
