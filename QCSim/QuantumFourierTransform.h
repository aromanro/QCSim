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
		}

		void IQFT(QubitRegister<VectorClass, MatrixClass>& reg) const
		{
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

		// the sign convention is not as in the Quantum Computation and Quantum Information book
		// also because of the qubits ordering, the circuit is 'mirrored' (they have the binary representation as j1 j2 j3 ... jn, I have it as jn...j1)

		void QFT(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int sq, unsigned int eq) const
		{
			Gates::ControlledPhaseShiftGate<MatrixClass> phaseShift;

			reg.ApplyGate(hadamard, eq);

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				double phase = -M_PI_2;
				for (int ctrlq = curQubitm1; ctrlq >= static_cast<int>(sq); --ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					reg.ApplyGate(phaseShift, curQubit, ctrlq);

					phase *= 0.5;
				}

				reg.ApplyGate(hadamard, curQubitm1);
			}

			Swap(reg, sq, eq);
		}

		// everything should be in reverse here:

		void IQFT(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int sq, unsigned int eq) const
		{
			Gates::ControlledPhaseShiftGate<MatrixClass> phaseShift;

			Swap(reg, sq, eq);

			for (unsigned int curQubit = sq + 1; curQubit <= eq; ++curQubit)
			{
				reg.ApplyGate(hadamard, curQubit - 1);

				int div = (int)pow(2, curQubit - sq);
				double phase = M_PI / div;
				for (unsigned int ctrlq = sq; ctrlq < curQubit; ++ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					reg.ApplyGate(phaseShift, curQubit, ctrlq);

					phase *= 2;
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


