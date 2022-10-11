#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

namespace BellInequalities {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class CheckCHSHInequality :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		CheckCHSHInequality(int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(2, addseed),
			S((-Q.getOperatorMatrix() - R.getOperatorMatrix()) / sqrt(2.)), // (-Z-X)/sqrt(2)=-H
			T((Q.getOperatorMatrix() - R.getOperatorMatrix()) / sqrt(2.)), // (Z-X)/sqrt(2)
			dist_bool(0, 1)
		{
			uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			timeSeed += addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };

			rng.seed(seed);
		};

		// TODO: Implement it
		// needs to save what the chosen operators were, to be able to compute stats
		// or at least accumulate here the values for statistics
		unsigned int Execute() override
		{
			// start with the Bell state:
			bellState.setBellState11(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);

			// Alice and Bob each get a qubit, they perform measurements on them

			// pick which one to measure at random
			// Alice: 
			const bool aM = dist_bool(gen);
			QC::SingleQubitGate<MatrixClass>& aliceMeasurement = aM ? R : Q;
			

			// Bob:
			const bool bM = dist_bool(gen);
			QC::SingleQubitGate<MatrixClass>& bobMeasurement = bm ? T : S;

			const unsigned int state = QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();

			return state;
		}

	protected:
		BellState<VectorClass, MatrixClass> bellState;

		// Alice measurement basis operators
		QC::PauliZGate<MatrixClass> Q;
		QC::PauliXGate<MatrixClass> R;

		// Bob measurement basis operators
		QC::SingleQubitGate<MatrixClass> S;
		QC::SingleQubitGate<MatrixClass> T;

		std::mt19937_64 rng;
		std::uniform_int_distribution<> dist_bool;
	};

}


