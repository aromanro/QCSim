#pragma once

#include "QubitRegister.h"
#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QC {

	// maybe it could be used in some other cases
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitsSwapper : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		QubitsSwapper(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: sQubit(startQubit), eQubit(std::max(startQubit, std::min(N - 1, endQubit)))
		{
		}

		unsigned int getStartQubit() const { return sQubit; };
		unsigned int getEndQubit() const { return eQubit; };

		void Swap(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int startQubit, unsigned int endQubit) const
		{
			while (startQubit < endQubit)
			{
				reg.ApplyGate(swapOp, startQubit, endQubit);
				++startQubit;
				--endQubit;
			}
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			Swap(reg, sQubit, eQubit);

			return reg.Measure();
		}

	protected:
		Gates::SwapGate<MatrixClass> swapOp;
		unsigned int sQubit;
		unsigned int eQubit;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumFourierTransform : public QubitsSwapper<VectorClass, MatrixClass>
	{
	public:
		QuantumFourierTransform(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: QubitsSwapper<VectorClass, MatrixClass>(N, startQubit, endQubit)
		{
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			QFT(reg, QubitsSwapper<VectorClass, MatrixClass>::sQubit, QubitsSwapper<VectorClass, MatrixClass>::eQubit);

			return reg.Measure();
		}

		// execute this to avoid measurement
		void QFT(QubitRegister<VectorClass, MatrixClass>& reg) const
		{
			QFT(reg, QubitsSwapper<VectorClass, MatrixClass>::sQubit, QubitsSwapper<VectorClass, MatrixClass>::eQubit);
		}

		void IQFT(QubitRegister<VectorClass, MatrixClass>& reg) const
		{
			QFT(reg, QubitsSwapper<VectorClass, MatrixClass>::sQubit, QubitsSwapper<VectorClass, MatrixClass>::eQubit, true);
		}

		// the sign convention is not as in the Quantum Computation and Quantum Information book
		// also because of the qubits ordering, the circuit is 'mirrored' (they have the binary representation as j1 j2 j3 ... jn, I have it as jn...j1)

		void QFT(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int sq, unsigned int eq, bool inverse = false) const
		{
			Gates::ControlledPhaseShiftGate<MatrixClass> phaseShift;

			reg.ApplyGate(hadamard, eq);

			const double startPhase = (inverse ? M_PI_2 : -M_PI_2);

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				const unsigned int curQubitm1 = curQubit - 1;
				double phase = startPhase; // starts from R2 = phase shift with 2 PI / 2^2 = PI / 2
				for (int ctrlq = curQubitm1; ctrlq >= static_cast<int>(sq); --ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					reg.ApplyGate(phaseShift, curQubit, ctrlq);

					phase *= 0.5;
				}

				reg.ApplyGate(hadamard, curQubitm1);
			}

			QubitsSwapper<VectorClass, MatrixClass>::Swap(reg, sq, eq);
		}
	
		Gates::HadamardGate<MatrixClass> hadamard; // public, let others use it
	};

}


