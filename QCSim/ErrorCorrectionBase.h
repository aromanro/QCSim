#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrectionBase :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrectionBase(unsigned int N, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed), errorQubit(N)
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
		}


		void SetState(std::complex<double> alpha, std::complex<double> beta)
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Clear();
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setRawAmplitude(0, alpha);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setRawAmplitude(1, beta);

			// ensure it's normalized:
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Normalize();
		}

		// set it to bigger than two for 'no error'
		void SetErrorQubit(unsigned int q = 0)
		{
			errorQubit = q;
		}

		unsigned int GetErrorQubit() const
		{
			return errorQubit;
		}

	protected:
		virtual void ApplyError() = 0;

		unsigned int errorQubit;

		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::ToffoliGate<MatrixClass> ccnot;
		QC::Gates::PauliXGate<MatrixClass> x; // error gate (flip or sign)
	};

}

