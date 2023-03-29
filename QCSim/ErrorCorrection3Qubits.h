#pragma once

#include "ErrorCorrectionBase.h"

namespace ErrorCorrection {


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3Qubits :
		public ErrorCorrectionBase<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = ErrorCorrectionBase<VectorClass, MatrixClass>;
		using AlgorithmClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		ErrorCorrection3Qubits(int addseed = 0)
			: BaseClass(3, addseed)
		{
		}

	protected:
		void Encode()
		{
			AlgorithmClass::ApplyGate(BaseClass::cnot, 1, 0);
			AlgorithmClass::ApplyGate(BaseClass::cnot, 2, 0);
		}

		void DetectAndCorrect()
		{
			AlgorithmClass::ApplyGate(BaseClass::cnot, 1, 0);
			AlgorithmClass::ApplyGate(BaseClass::cnot, 2, 0);
			AlgorithmClass::ApplyGate(BaseClass::ccnot, 0, 1, 2);
		}
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3QubitsFlip :
		public ErrorCorrection3Qubits<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = ErrorCorrection3Qubits<VectorClass, MatrixClass>;
		using AlgorithmClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		ErrorCorrection3QubitsFlip(int addseed = 0)
			: BaseClass(addseed)
		{
		}

		unsigned int Execute() override
		{
			BaseClass::Encode();

			ApplyError();

			BaseClass::DetectAndCorrect();

			return AlgorithmClass::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError() override
		{
			if (BaseClass::errorQubit > 2) return;

			AlgorithmClass::ApplyGate(BaseClass::BaseClass::x, BaseClass::BaseClass::errorQubit);
		}
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3QubitsSign :
		public ErrorCorrection3Qubits<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = ErrorCorrection3Qubits<VectorClass, MatrixClass>;
		using AlgorithmClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		ErrorCorrection3QubitsSign(int addseed = 0)
			: BaseClass(addseed)
		{
		}

		unsigned int Execute() override
		{
			BaseClass::Encode();
			ApplyHadamardOnAllQubits();

			ApplyError();

			ApplyHadamardOnAllQubits(); 
			BaseClass::DetectAndCorrect();

			return AlgorithmClass::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError() override
		{
			if (BaseClass::BaseClass::errorQubit > 2) return;

			AlgorithmClass::ApplyGate(hadamard, BaseClass::BaseClass::errorQubit);
			AlgorithmClass::ApplyGate(BaseClass::BaseClass::x, BaseClass::BaseClass::errorQubit);
			AlgorithmClass::ApplyGate(hadamard, BaseClass::BaseClass::errorQubit);
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < AlgorithmClass::getNrQubits(); ++i)
				AlgorithmClass::ApplyGate(hadamard, i);
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};
}

