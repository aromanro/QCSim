#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

#define _USE_MATH_DEFINES
#include <math.h>

// for details see for example "Experimenting quantum phenomena on NISQ computers using high level quantum programming" by Duc M. Tran and Hung Q. Nguyen
// https://arxiv.org/abs/2111.02896v2

namespace Paradoxes {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class HardyParadox :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		HardyParadox(int addseed = 0)
			: BaseClass(3, addseed), theta0(0.575 * M_PI), theta1(0.575 * M_PI)
		{
		}

		void setTheta0(double t)
		{
			theta0 = t;
		}

		double getTheta0() const
		{
			return theta0;
		}

		void setTheta1(double t)
		{
			theta1 = t;
		}

		double getTheta1() const
		{
			return theta1;
		}

		unsigned int Execute() override
		{
			BaseClass::setToBasisState(0);

			ryGate.SetTheta(theta0);
			BaseClass::ApplyGate(ryGate, 0);
			ryGate.SetTheta(theta1);
			BaseClass::ApplyGate(ryGate, 1);

			BaseClass::ApplyGate(ccnot, 2, 0, 1);

			ryGate.SetTheta(M_PI - theta0);
			BaseClass::ApplyGate(ryGate, 0);
			ryGate.SetTheta(M_PI - theta1);
			BaseClass::ApplyGate(ryGate, 1);

			return BaseClass::Measure();
		}

		double TheoreticalGamma() const
		{
			const double s0 = sin(theta0);
			const double s02 = sin(0.5 * theta0);
			const double s1 = sin(theta1);
			const double c0 = cos(theta0);
			const double c1 = cos(theta1);

			return 0.25 * s0 * s0 * s1 * s1 / (2. * s02 * s02 * c1 + c0 + 3.);
		}

	protected:
		double theta0;
		double theta1;

		QC::Gates::RyGate<MatrixClass> ryGate;
		QC::Gates::ToffoliGate<MatrixClass> ccnot;
	};

}


