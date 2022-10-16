#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace DeutschJozsa {

	// TODO: Just a skeleton, implement it!
	// N = 2 reduces to Deutsch's algorithm
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DeutschJozsaAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		DeutschJozsaAlgorithm(unsigned int N = 2, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed)
		{
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

		QC::HadamardGate<MatrixClass> hadamard;
	};

}