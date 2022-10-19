#pragma once

#include "QubitRegister.h"
#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumFourierTransform : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		QuantumFourierTransform(QubitRegister<VectorClass, MatrixClass>& r, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: QuantumSubAlgorithm<VectorClass, MatrixClass>::QuantumSubAlgorithm(r), sQubit(startQubit), eQubit(std::min(r.getNrQubits() - 1, endQubit))
		{
		}

		unsigned int Execute() override
		{
			QFT(sQubit, eQubit);

			return QC::QuantumSubAlgorithm<VectorClass, MatrixClass>::Measure();
		}

		// execute this to avoid measurement
		void QFT()
		{
			QFT(sQubit, eQubit);
			Swap(sQubit, eQubit);
		}

		void IQFT()
		{
			Swap(sQubit, eQubit);
			IQFT(sQubit, eQubit);
		}

		unsigned int getStartQubit() const { return sQubit; };
		unsigned int getEndQubit() const { return eQubit; };

	protected:
		void Swap(unsigned int startQubit, unsigned int endQubit)
		{
			while (startQubit < endQubit)
			{
				QC::QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(swapOp, startQubit, endQubit);
				++startQubit;
				--endQubit;
			}
		}

		void QFT(unsigned int sq, unsigned int eq)
		{
			QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, eq);

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				int div = (int)pow(2, eq - curQubit + 1);
				double phase = M_PI / div;
				for (unsigned int ctrlq = eq; ctrlq >= curQubit; --ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(phaseShift, ctrlq, curQubitm1);

					phase *= 2;
				}

				QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, curQubitm1);
			}
		}

		void IQFT(unsigned int sq, unsigned int eq)
		{
			for (unsigned int curQubit = sq; curQubit < eq; ++curQubit)
			{
				QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, curQubit);

				double phase = M_PI_2;
				for (unsigned int ctrlq = curQubit + 1; ctrlq <= eq; ++ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(phaseShift, curQubit, ctrlq);

					phase /= 2;
				}
			}
			QuantumSubAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, eq);
		}

	public:
		HadamardGate<MatrixClass> hadamard; // public, let others use it

	protected:
		ControlledPhaseShiftGate<MatrixClass> phaseShift;
		SwapGate<MatrixClass> swapOp;

		unsigned int sQubit;
		unsigned int eQubit;
	};

}


