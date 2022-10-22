#pragma once

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"
#include "QuantumFourierTransform.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QuantumSimulation {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SchrodingerSimulation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		SchrodingerSimulation(unsigned int N = 5, double t = 0, unsigned int nrSteps = 1, bool addSeed = false)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addSeed), simTime(t), steps(nrSteps), fourier(N, 0, N - 1)
		{
		}

		unsigned int Execute() override
		{
			// TODO: implement it
			return 0;
		}

	protected:
		double simTime;
		unsigned int steps;

		QC::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};
}


