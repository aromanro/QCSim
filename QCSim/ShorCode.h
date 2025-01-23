#pragma once

#include "ErrorCorrectionBase.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorCode :
		public ErrorCorrectionBase<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = ErrorCorrectionBase<VectorClass, MatrixClass>;
		using AlgorithmClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		ShorCode(unsigned int addseed = 0)
			: BaseClass(9, addseed), errorType(Flip), iy(std::complex<double>(0., 1.) * y.getRawOperatorMatrix())
		{
		}

		enum ErrorType : size_t
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

		size_t Execute() override
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
				AlgorithmClass::ApplyGate(iy, BaseClass::errorQubit);
		}

		void Encode()
		{
			AlgorithmClass::ApplyGate(BaseClass::cnot, 3, 0);
			AlgorithmClass::ApplyGate(BaseClass::cnot, 6, 0);

			for (size_t qubit = 0; qubit <= 6; qubit += 3)
			{
				// the following make up a generalized entangling gate for 3 qubits
				// the first two being the two-qubits 'entangling gate'
				// such a gate would create a GHZ state from |000> (that is 1/sqrt(2) * (|000> + |111>))
				AlgorithmClass::ApplyGate(hadamard, qubit);
				AlgorithmClass::ApplyGate(BaseClass::cnot, qubit + 1, qubit);
				AlgorithmClass::ApplyGate(BaseClass::cnot, qubit + 2, qubit);
			}
		}

		void DetectAndCorrect()
		{
			for (size_t qubit = 0; qubit <= 6; qubit += 3)
			{
				AlgorithmClass::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, qubit + 1, qubit);
				AlgorithmClass::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, qubit + 2, qubit);
				AlgorithmClass::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::ccnot, qubit, qubit + 1, qubit + 2);
				AlgorithmClass::ApplyGate(hadamard, qubit);
			}

			AlgorithmClass::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 3, 0);
			AlgorithmClass::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 6, 0);
			AlgorithmClass::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::ccnot, 0, 3, 6);
		}

		ErrorType errorType;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliZGate<MatrixClass> z;
		QC::Gates::PauliYGate<MatrixClass> y;
		QC::Gates::SingleQubitGate<MatrixClass> iy;
	};

}