#pragma once

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"
#include "QuantumFourierTransform.h"

namespace Adders {

	// see https://arxiv.org/abs/quant-ph/0008033

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DraperAdder : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		DraperAdder(unsigned int N = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(2 * N, addseed), n(N), fourier(2 * N, N)
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
		}

		unsigned int Execute() override
		{
			fourier.IQFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg, false);
			
			// the operations commute with each other, so one can change the ordering

			for (unsigned int qbit = static_cast<int>(n); qbit < static_cast<int>(2 * n); ++qbit)
			{
				double phase = M_PI;
				for (int cbit = static_cast<int>(qbit) - n; cbit >= 0; --cbit)
				{
					fourier.cPhaseShift.SetPhaseShift(phase);
					QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(fourier.cPhaseShift, qbit, cbit);
					phase *= 0.5;
				}
			}
			
			fourier.QFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg, false);

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}


	protected:
		unsigned int n;
		QC::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DraperAdderWithCarry : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		DraperAdderWithCarry(unsigned int N = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(2 * N + 1, addseed), n(N), fourier(2 * N + 1, N)
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
		}

		unsigned int Execute() override
		{
			fourier.IQFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg, false);

			// the operations commute with each other, so one can change the ordering

			// the carry qubit
			double phase = M_PI_2;
			for (int cbit = n - 1; cbit >= 0; --cbit)
			{
				fourier.cPhaseShift.SetPhaseShift(phase);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(fourier.cPhaseShift,2 * n, cbit);
				phase *= 0.5;
			}

			// the other ones
			for (unsigned int qbit = static_cast<int>(n); qbit < static_cast<int>(2 * n); ++qbit)
			{
				double phase = M_PI;
				for (int cbit = static_cast<int>(qbit) - n; cbit >= 0; --cbit)
				{
					fourier.cPhaseShift.SetPhaseShift(phase);
					QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(fourier.cPhaseShift, qbit, cbit);
					phase *= 0.5;
				}
			}

			fourier.QFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg, false);

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}


	protected:
		unsigned int n;
		QC::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};

}


