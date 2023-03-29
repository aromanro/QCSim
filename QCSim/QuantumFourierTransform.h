#pragma once

#include "QubitsSwapper.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumFourierTransform : public QubitsSwapper<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QubitsSwapper<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		QuantumFourierTransform(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: BaseClass(N, startQubit, endQubit)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			QFT(reg);

			return reg.Measure();
		}

		// execute this to avoid measurement

		// the sign convention is not as in the Quantum Computation and Quantum Information book
		// also because of the qubits ordering, the circuit is 'mirrored' (they have the binary representation as j1 j2 j3 ... jn, I have it as jn...j1)

		void QFT(RegisterClass& reg, bool doSwap = true)
		{
			const unsigned int sq = BaseClass::BaseClass::getStartQubit();
			const unsigned int eq = BaseClass::BaseClass::getEndQubit();

			reg.ApplyGate(hadamard, eq);

			const double startPhase = M_PI_2;

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				double phase = startPhase; // starts from R2 = phase shift with 2 PI / 2^2 = PI / 2
				for (int ctrlq = curQubitm1; ctrlq >= static_cast<int>(sq); --ctrlq)
				{
					cPhaseShift.SetPhaseShift(phase);
					reg.ApplyGate(cPhaseShift, curQubit, ctrlq);

					phase *= 0.5;
				}

				reg.ApplyGate(hadamard, curQubitm1);
			}

			if (doSwap) BaseClass::Swap(reg);
		}

		void IQFT(RegisterClass& reg, bool doSwap = true)
		{
			const unsigned int sq = BaseClass::BaseClass::getStartQubit();
			const unsigned int eq = BaseClass::BaseClass::getEndQubit();

			if (doSwap) BaseClass::Swap(reg);

			const double startPhase = -M_PI_2;

			for (unsigned int curQubit = sq + 1; curQubit <= eq; ++curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				reg.ApplyGate(hadamard, curQubitm1);

				double phase = startPhase; // starts from R2 = phase shift with 2 PI / 2^2 = PI / 2
				for (int ctrlq = curQubitm1; ctrlq >= static_cast<int>(sq); --ctrlq)
				{
					cPhaseShift.SetPhaseShift(phase);
					reg.ApplyGate(cPhaseShift, curQubit, ctrlq);

					phase *= 0.5;
				}	
			}

			reg.ApplyGate(hadamard, eq);
		}
	
		// public, let others use them
		Gates::HadamardGate<MatrixClass> hadamard; 
		Gates::ControlledPhaseShiftGate<MatrixClass> cPhaseShift;
	};

}


