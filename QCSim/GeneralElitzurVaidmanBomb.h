#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

#define _USE_MATH_DEFINES
#include <math.h>

// for details see for example "Experimenting quantum phenomena on NISQ computers using high level quantum programming" by Duc M. Tran and Hung Q. Nguyen
// https://arxiv.org/abs/2111.02896v2
// the formula 4 seem to be wrong, though, on the version I looked over

// see also "Interaction-free measurements and counterfactual computation in IBM quantum computers" by J. Alberto Casas and Bryan Zaldivar
// https://arxiv.org/abs/2005.03547

// the one stage Elitzur-Vaidman bomb is basically equivalent with the quantum eraser circuit without the eraser set
// and with the following one with a single stage

namespace Paradoxes {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class GeneralElitzurVaidmanBomb :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		GeneralElitzurVaidmanBomb(unsigned int maxStages, int addseed = 0)
			: BaseClass(maxStages + 1, addseed), stages(maxStages), theta(M_PI / (maxStages + 1.))
		{
			assert((BaseClass::getNrQubits() >= 2));
		}

		void setStages(unsigned int s)
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();
			
			if (s >= nrQubits) s = nrQubits - 1;
			else if (s == 0) s = 1;
			stages = s;
		}

		unsigned int getStages() const
		{
			return stages;
		}

		void setTheta(double t)
		{
			theta = t;
		}

		double getTheta() const
		{
			return theta;
		}

		double getThetai() const
		{
			return (M_PI - theta) / stages;
		}

		unsigned int Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure();
		}

		std::map<unsigned int, unsigned int> ExecuteWithMultipleMeasurements(unsigned int nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(nrMeasurements);
		}

		double TheoreticalEfficiency() const
		{
			// value for a single stage, M_PI / 2 theta (that is, equal beam split)
			// gets 1 / 3 as expected
			const double thetai2 = getThetai();
			const double c = cos(0.5 * thetai2);
			const double c2 = c * c;

			const double s = sin(0.5 * theta);
			const double s2 = s * s;

			const double cl = cos(0.5 * theta);
			const double cl2 = cl * cl;

			const double prod = pow(c2, stages);

			return cl2 * prod / (1. - s2 * prod);
		}

	protected:
		void Init()
		{
			BaseClass::setToBasisState(0);
			// the following has the role of the first beam splitter:
			ryGate.SetTheta(getThetai());
			BaseClass::ApplyGate(ryGate, 0);
		}

		void ExecuteWithoutMeasurement()
		{
			Init();
			// now we're in the state given by the first beam splitter

			for (unsigned int stage = 1; stage < stages; ++stage)
			{
				BaseClass::ApplyGate(cnot, stage);
				BaseClass::ApplyGate(ryGate, 0);
			}

			ryGate.SetTheta(theta);
			BaseClass::ApplyGate(cnot, stages);
			BaseClass::ApplyGate(ryGate, 0);
		}

		unsigned int stages;
		double theta;

		QC::Gates::RyGate<MatrixClass> ryGate;
		QC::Gates::CNOTGate<MatrixClass> cnot;
	};

}

