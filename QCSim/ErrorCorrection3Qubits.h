#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace ErrorCorrection {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ErrorCorrection3Qubits :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		ErrorCorrection3Qubits(int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(3, addseed), errorQubit(3)
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
		

		void Prepare()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 2, 1);
		}

		void DetectAndCorrect()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 2, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ccnot, 0, 1, 2);
		}

		unsigned int errorQubit;

		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::ToffoliGate<MatrixClass> ccnot;
		QC::Gates::PauliXGate<MatrixClass> x; // error gate (flip or sign)
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
			ErrorCorrection3Qubits<VectorClass, MatrixClass>::Prepare();

			ApplyError();

			ErrorCorrection3Qubits<VectorClass, MatrixClass>::DetectAndCorrect();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError() override
		{
			if (ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit > 2) return;

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrection3Qubits<VectorClass, MatrixClass>::x, ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit);
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
			ErrorCorrection3Qubits<VectorClass, MatrixClass>::Prepare();
			ApplyHadamardOnAllQubits();

			ApplyError();

			ApplyHadamardOnAllQubits(); 
			ErrorCorrection3Qubits<VectorClass, MatrixClass>::DetectAndCorrect();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(1, 2); // 11 - qubit 0 was flipped and corrected, 01 - qubit 1 was flipped, 10 - qubit 2 was flipped
		}

	protected:
		void ApplyError() override
		{
			if (ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit > 2) return;

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(ErrorCorrection3Qubits<VectorClass, MatrixClass>::x, ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, ErrorCorrection3Qubits<VectorClass, MatrixClass>::errorQubit);
		}

		void ApplyHadamardOnAllQubits()
		{
			for (unsigned int i = 0; i < QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits(); ++i)
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, i);
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};
}

