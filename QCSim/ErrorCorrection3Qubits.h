#pragma once

#include "ErrorCorrectionBase.h"

namespace ErrorCorrection {


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3Qubits :
		public ErrorCorrectionBase<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrection3Qubits(int addseed = 0)
			: ErrorCorrectionBase<VectorClass, MatrixClass>(3, addseed)
		{
		}

	protected:
		void Encode()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 2, 0);
		}

		void DetectAndCorrect()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::cnot, 2, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::ccnot, 0, 1, 2);
		}
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3QubitsFlip :
		public ErrorCorrection3Qubits<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrection3QubitsFlip(int addseed = 0)
			: ErrorCorrection3Qubits<VectorClass, MatrixClass>(addseed)
		{
		}

		unsigned int Execute() override
		{
			ErrorCorrection3Qubits<VectorClass, MatrixClass>::Encode();

			ApplyError();

			ErrorCorrection3Qubits<VectorClass, MatrixClass>::DetectAndCorrect();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError() override
		{
			if (ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit > 2) return;

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::x, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
		}
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3QubitsSign :
		public ErrorCorrection3Qubits<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrection3QubitsSign(int addseed = 0)
			: ErrorCorrection3Qubits<VectorClass, MatrixClass>(addseed)
		{
		}

		unsigned int Execute() override
		{
			ErrorCorrection3Qubits<VectorClass, MatrixClass>::Encode();
			ApplyHadamardOnAllQubits();

			ApplyError();

			ApplyHadamardOnAllQubits(); 
			ErrorCorrection3Qubits<VectorClass, MatrixClass>::DetectAndCorrect();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError() override
		{
			if (ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit > 2) return;

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrectionBase<VectorClass, MatrixClass>::x, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrectionBase<VectorClass, MatrixClass>::errorQubit);
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits(); ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, i);
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};
}

