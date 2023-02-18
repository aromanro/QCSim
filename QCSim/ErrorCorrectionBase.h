#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrectionBase :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		typedef QC::QuantumAlgorithm<VectorClass, MatrixClass> BaseClass;

		ErrorCorrectionBase(unsigned int N, int addseed = 0)
			: BaseClass(N, addseed), errorQubit(N)
		{
			BaseClass::setToBasisState(0);
		}


		void SetState(std::complex<double> alpha, std::complex<double> beta)
		{
			BaseClass::Clear();
			BaseClass::setRawAmplitude(0, alpha);
			BaseClass::setRawAmplitude(1, beta);

			// ensure it's normalized:
			BaseClass::Normalize();
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

