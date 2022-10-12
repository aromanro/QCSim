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
			S((-Q.getRawOperatorMatrix() - R.getRawOperatorMatrix()) / sqrt(2.)), // (-Z-X)/sqrt(2)=-H
			T((Q.getRawOperatorMatrix() - R.getRawOperatorMatrix()) / sqrt(2.)), // (Z-X)/sqrt(2)
			dist_bool(0, 1),
			QSaccum(0), RSaccum(0), RTaccum(0), QTaccum(0), QScount(0), RScount(0), RTcount(0), QTcount(0)
		{
			uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			timeSeed += addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);
		};

		void ResetStatistics()
		{
			QSaccum = RSaccum = RTaccum = QTaccum = 0;
			QScount = RScount = RTcount = QTcount = 0;
		}

		double getValue() const
		{
			return (QScount ? static_cast<double>(QSaccum)/QScount : 0.) + (RScount ? static_cast<double>(RSaccum)/RScount : 0.) + (RTcount ? static_cast<double>(RTaccum) / RTcount : 0.) - (QTcount ? static_cast<double>(QTaccum) / QTcount : 0.);
		}

		unsigned int Execute() override
		{
			// start with the Bell state:
			bellState.setBellState11(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);

			// Alice and Bob each get a qubit, they perform measurements on them

			// pick which one to measure at random
			// Alice: 
			const bool aM = dist_bool(rng) == 1;
			const QC::SingleQubitGate<MatrixClass>& aliceMeasurement = aM ? dynamic_cast<QC::SingleQubitGate<MatrixClass>&>(R) : dynamic_cast<QC::SingleQubitGate<MatrixClass>&>(Q);
			measurementBasis.switchToOperatorBasis(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg, aliceMeasurement.getRawOperatorMatrix(), 0);

			//unsigned int state = QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(0,0);
			//const int res1 = (state & 1) ? 1 : -1;

			// Bob:
			const bool bM = dist_bool(rng) == 1;
			const QC::SingleQubitGate<MatrixClass>& bobMeasurement = bM ? dynamic_cast<QC::SingleQubitGate<MatrixClass>&>(T) : dynamic_cast<QC::SingleQubitGate<MatrixClass>&>(S);
			measurementBasis.switchToOperatorBasis(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg, bobMeasurement.getRawOperatorMatrix(), 1);

			//state = QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, 1);
			//const int res2 = (state & 1) ? 1 : -1;
			
			const unsigned int state = QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
			const int res1 = (state & 1) ? 1 : -1;
			const int res2 = (state & 2) ? 1 : -1;

			const int prod = res1 * res2;

			if (aM && bM) // RT
			{
				RTaccum += prod;
				++RTcount;
			}
			else if (aM) // RS
			{
				RSaccum += prod;
				++RScount;
			}
			else if (bM) // QT
			{
				QTaccum += prod;
				++QTcount;
			}
			else // QS
			{
				QSaccum += prod;
				++QScount;
			}

			return state;
		}

	protected:
		QC::BellState<VectorClass, MatrixClass> bellState;
		QC::MeasurementBasis<VectorClass, MatrixClass> measurementBasis;

		// Alice measurement basis operators
		QC::PauliZGate<MatrixClass> Q;
		QC::PauliXGate<MatrixClass> R;

		// Bob measurement basis operators
		QC::SingleQubitGate<MatrixClass> S;
		QC::SingleQubitGate<MatrixClass> T;

		std::mt19937_64 rng;
		std::uniform_int_distribution<> dist_bool;

		long long int QSaccum;
		long long int RSaccum;
		long long int RTaccum;
		long long int QTaccum;
		unsigned long long int QScount;
		unsigned long long int RScount;
		unsigned long long int RTcount;
		unsigned long long int QTcount;
	};

}


