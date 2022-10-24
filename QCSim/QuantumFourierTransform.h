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
		QuantumFourierTransform(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: sQubit(startQubit), eQubit(std::max(startQubit, std::min(N - 1, endQubit)))
		{
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			QFT(reg, sQubit, eQubit);

			return reg.Measure();
		}

		// execute this to avoid measurement
		void QFT(QubitRegister<VectorClass, MatrixClass>& reg) const
		{
			QFT(reg, sQubit, eQubit);
			Swap(reg, sQubit, eQubit);
		}

		void IQFT(QubitRegister<VectorClass, MatrixClass>& reg) const
		{
			Swap(reg, sQubit, eQubit);
			IQFT(reg, sQubit, eQubit);
		}

		unsigned int getStartQubit() const { return sQubit; };
		unsigned int getEndQubit() const { return eQubit; };

	protected:
		void Swap(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int startQubit, unsigned int endQubit) const
		{
			while (startQubit < endQubit)
			{
				reg.ApplyGate(swapOp, startQubit, endQubit);
				++startQubit;
				--endQubit;
			}
		}

		void QFT(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int sq, unsigned int eq) const
		{
			reg.ApplyGate(hadamard, eq);
			
			Gates::ControlledPhaseShiftGate<MatrixClass> phaseShift;

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				int div = (int)pow(2, eq - curQubit + 1);
				double phase = M_PI / div;
				for (unsigned int ctrlq = eq; ctrlq >= curQubit; --ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					reg.ApplyGate(phaseShift, ctrlq, curQubitm1);

					phase *= 2;
				}

				reg.ApplyGate(hadamard, curQubitm1);
			}
		}

		void IQFT(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int sq, unsigned int eq) const
		{
			for (unsigned int curQubit = sq; curQubit < eq; ++curQubit)
			{
				reg.ApplyGate(hadamard, curQubit);

				Gates::ControlledPhaseShiftGate<MatrixClass> phaseShift;

				double phase = M_PI_2;
				for (unsigned int ctrlq = curQubit + 1; ctrlq <= eq; ++ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					reg.ApplyGate(phaseShift, curQubit, ctrlq);

					phase /= 2;
				}
			}
			reg.ApplyGate(hadamard, eq);
		}

	public:
		Gates::HadamardGate<MatrixClass> hadamard; // public, let others use it

	protected:
		Gates::SwapGate<MatrixClass> swapOp;

		unsigned int sQubit;
		unsigned int eQubit;
	};

}


