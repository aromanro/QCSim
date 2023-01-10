#pragma once

#include "ErrorCorrectionBase.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorCode :
		public ErrorCorrectionBase<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrection3Qubits(int addseed = 0)
			: ErrorCorrectionBase<VectorClass, MatrixClass>(9, addseed), errorType(Flip)
		{
		}

		enum ErrorType : int
		{
			Flip,
			Sign
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
			Prepare();

			ApplyError();

			DetectAndCorrect();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);
		}

	protected:
		void ApplyError() override
		{
			if (ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit >= QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits()) return;

			if(ErrorType == Flip)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::x, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
			else
			{
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::x, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
			}
		}

		void Prepare()
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
	};

}