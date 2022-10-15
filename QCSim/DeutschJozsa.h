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
	};

}