#pragma once

#include "ErrorCorrectionBase.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorCode :
		public ErrorCorrectionBase<VectorClass, MatrixClass>
	{
	public:
		ShorCode(int addseed = 0)
			: ErrorCorrectionBase<VectorClass, MatrixClass>(9, addseed), errorType(Flip)
		{
		}

		enum ErrorType : unsigned int
		{
			Flip,
			Sign,
			Both,
			None // another way of specifying 'no error' (the other one is to set the error qubit larger than the last qubit number)
		};

		void setErrorType(ErrorType t)
		{
			errorType = t;
		}

		ErrorType getErrorType() const
		{
			return errorType;
		}

		unsigned int Execute() override
		{
			Encode();

			ApplyError();

			DetectAndCorrect();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);
		}

	protected:
		void ApplyError() override
		{
			if (ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit >= QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits()) return;

			if (errorType == Flip)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::x, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
			else if (errorType == Sign)
			{
				//QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
				//QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::x, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
				//QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(z, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit); // another way of the above
			}
			else if (errorType == Both)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyOperatorMatrix(std::complex<double>(0., 1.) * y.getOperatorMatrix(QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits(), ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit));
		}

		void Encode()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 3, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 6, 0);

			for (unsigned int qubit = 0; qubit <= 6; qubit += 3)
			{
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, qubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, qubit + 1, qubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, qubit + 2, qubit);
			}
		}

		void DetectAndCorrect()
		{
			for (unsigned int qubit = 0; qubit <= 6; qubit += 3)
			{
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, qubit + 1, qubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, qubit + 2, qubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::ccnot, qubit, qubit + 1, qubit + 2);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, qubit);
			}

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 3, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 6, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::ccnot, 0, 3, 6);
		}

		ErrorType errorType;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliZGate<MatrixClass> z;
		QC::Gates::PauliYGate<MatrixClass> y;
	};

}