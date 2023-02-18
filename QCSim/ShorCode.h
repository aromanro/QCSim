#pragma once

#include "ErrorCorrectionBase.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorCode :
		public ErrorCorrectionBase<VectorClass, MatrixClass>
	{
	public:
		typedef ErrorCorrectionBase<VectorClass, MatrixClass> BaseClass;
		typedef QC::QuantumAlgorithm<VectorClass, MatrixClass> AlgorithmClass;

		ShorCode(int addseed = 0)
			: BaseClass(9, addseed), errorType(Flip)
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

			return AlgorithmClass::Measure(1, AlgorithmClass::getNrQubits() - 1);
		}

	protected:
		void ApplyError() override
		{
			if (BaseClass::errorQubit >= AlgorithmClass::getNrQubits()) return;

			if (errorType == Flip)
				AlgorithmClass::ApplyGate(BaseClass::x, BaseClass::errorQubit);
			else if (errorType == Sign)
			{
				//AlgorithmClass::ApplyGate(hadamard, BaseClass::errorQubit);
				//AlgorithmClass::ApplyGate(BaseClass::x,BaseClass::errorQubit);
				//AlgorithmClass::ApplyGate(hadamard, BaseClass::errorQubit);
				AlgorithmClass::ApplyGate(z, BaseClass::errorQubit); // another way of the above
			}
			else if (errorType == Both)
				AlgorithmClass::ApplyOperatorMatrix(std::complex<double>(0., 1.) * y.getOperatorMatrix(AlgorithmClass::getNrQubits(), BaseClass::errorQubit));
		}

		void Encode()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(BaseClass::cnot, 3, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(BaseClass::cnot, 6, 0);

			for (unsigned int qubit = 0; qubit <= 6; qubit += 3)
			{
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, qubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(BaseClass::cnot, qubit + 1, qubit);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(BaseClass::cnot, qubit + 2, qubit);
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