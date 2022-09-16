#pragma once

#include "QubitRegister.h"
#include "QuantumGate.h"


#define _USE_MATH_DEFINES
#include <math.h>

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumFourierTransform : public QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		QuantumFourierTransform(unsigned int N = 3, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX, int addseed = 0)
			: QuantumAlgorithm<VectorClass, MatrixClass>::QuantumAlgorithm(N, addseed), sQubit(startQubit), eQubit(std::min(N - 1, endQubit))
		{	
		}

		unsigned int Execute() override
		{
			QFT(sQubit, eQubit);

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.Measure();
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
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(swapOp, startQubit, endQubit);
				++startQubit;
				--endQubit;
			}
		}

		void QFT(unsigned int sq, unsigned int eq)
		{			
			QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, eq);

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				int div = (int)pow(2, eq - curQubit + 1);
				double phase = M_PI / div;
				for (unsigned int ctrlq = eq; ctrlq >= curQubit; --ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(phaseShift, ctrlq, curQubitm1);

					phase *= 2;
				}

				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, curQubitm1);
			}
		}

		void IQFT(unsigned int sq, unsigned int eq)
		{
			for (unsigned int curQubit = sq; curQubit < eq; ++curQubit)
			{
				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, curQubit);

				double phase = M_PI_2;
				for (unsigned int ctrlq = curQubit + 1; ctrlq <= eq; ++ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(phaseShift, curQubit, ctrlq);

					phase /= 2;
				}
			}
			QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, eq);
		}

		HadamardGate<MatrixClass> hadamard;
		ControlledPhaseShiftGate<MatrixClass> phaseShift;
		SwapGate<MatrixClass> swapOp;

		unsigned int sQubit;
		unsigned int eQubit;
	};

}


