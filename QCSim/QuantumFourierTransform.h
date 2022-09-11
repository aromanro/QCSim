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
		}

		void IQFT()
		{
			QFT(sQubit, eQubit);
		}

		unsigned int getStartQubit() const { return sQubit; };
		unsigned int getEndQubit() const { return eQubit; };

	protected:

		// TODO: Both QFT and IQFT need after each step to 'flip' the outputs
		// or the below implementation needs rethinking to have it 'upside down'

		void QFT(unsigned int sq, unsigned int eq)
		{
			if (sq == eq)
				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, eq);
			else
			{
				QFT(sq + 1, eq);

				for (unsigned int i = eq; i > sq; --i)
				{
					phaseShift.SetPhaseShift(M_PI / pow(2, sq - i));
					QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(phaseShift, i, sq);
				}

				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, sq);
			}
		}


		void IQFT(unsigned int sq, unsigned int eq)
		{
			if (sq == eq)
				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, eq);
			else
			{
				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, sq);

				double phase = M_PI / 2.;
				for (unsigned int i = sq + 1; i <= eq; ++i)
				{
					phaseShift.SetPhaseShift(phase);
					phase *= 0.5;
					QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(phaseShift, i, sq);
				}

				IQFT(sq + 1, eq);
			}
		}

		HadamardGate<MatrixClass> hadamard;
		ControlledPhaseShiftGate<MatrixClass> phaseShift;

		unsigned int sQubit;
		unsigned int eQubit;
	};

}


