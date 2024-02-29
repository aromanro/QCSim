#pragma once

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"
#include "QuantumFourierTransform.h"

namespace Adders {

	// see https://arxiv.org/abs/quant-ph/0008033

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DraperAdder : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DraperAdder(size_t N = 3, int addseed = 0)
			: BaseClass(2 * N, addseed), n(N), fourier(2 * N, N)
		{
			BaseClass::setToBasisState(0);
		}

		size_t Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure();
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(nrMeasurements);
		}


	protected:
		void ExecuteWithoutMeasurement()
		{
			fourier.QFT(BaseClass::reg, false);

			// the operations commute with each other, so one can change the ordering

			for (int qbit = static_cast<int>(n); qbit < static_cast<int>(2 * n); ++qbit)
			{
				double phase = M_PI;
				for (int cbit = qbit - n; cbit >= 0; --cbit)
				{
					fourier.cPhaseShift.SetPhaseShift(phase);
					BaseClass::ApplyGate(fourier.cPhaseShift, qbit, cbit);
					phase *= 0.5;
				}
			}

			fourier.IQFT(BaseClass::reg, false);
		}

		size_t n;
		QC::SubAlgo::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DraperAdderWithCarry : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DraperAdderWithCarry(size_t N = 3, int addseed = 0)
			: BaseClass(2 * N + 1, addseed), n(N), fourier(2 * N + 1, N)
		{
			BaseClass::setToBasisState(0);
		}

		size_t Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure();
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(nrMeasurements);
		}

	protected:
		void ExecuteWithoutMeasurement()
		{
			fourier.QFT(BaseClass::reg, false);

			// the operations commute with each other, so one can change the ordering

			// the carry qubit
			double phase = M_PI_2;
			for (int cbit = n - 1; cbit >= 0; --cbit)
			{
				fourier.cPhaseShift.SetPhaseShift(phase);
				BaseClass::ApplyGate(fourier.cPhaseShift, 2 * n, cbit);
				phase *= 0.5;
			}

			// the other ones
			for (int qbit = static_cast<int>(n); qbit < static_cast<int>(2 * n); ++qbit)
			{
				phase = M_PI;
				for (int cbit = qbit - n; cbit >= 0; --cbit)
				{
					fourier.cPhaseShift.SetPhaseShift(phase);
					BaseClass::ApplyGate(fourier.cPhaseShift, qbit, cbit);
					phase *= 0.5;
				}
			}

			fourier.IQFT(BaseClass::reg, false);
		}

		size_t n;
		QC::SubAlgo::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};

}


