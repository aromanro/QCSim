#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3Qubits :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrection3Qubits(int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(3, addseed), flippedQubit(0)
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

		// set it to bigger than two for 'no flip error'
		void SetFlipQubit(unsigned int q = 0)
		{
			flippedQubit = q;
		}

		unsigned int GetFlipQubit() const
		{
			return flippedQubit;
		}

		unsigned int Execute() override
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 2, 1);

			ApplyError();

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 2, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ccnot, 0, 1, 2);

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError()
		{
			if (flippedQubit > 2) return;

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(x, flippedQubit);
		}

		unsigned int flippedQubit;

		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::ToffoliGate<MatrixClass> ccnot;
		QC::Gates::PauliXGate<MatrixClass> x; // flip error gate
	};

}

