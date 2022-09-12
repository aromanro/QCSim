#pragma once

#include <iostream>

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
			//Swap(sQubit, eQubit);
			IQFT(sQubit, eQubit);
			Swap(sQubit, eQubit);
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
			//std::cout << "Hadamard on " << eq << std::endl;
			//MatrixClass m = hadamard.getOperatorMatrix(QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits(),eq);

			for (unsigned int curQubit = eq; curQubit > sq; --curQubit)
			{
				int div = (int)pow(2, eq - curQubit + 1);
				double phase = M_PI / div;
				for (unsigned int ctrlq = eq; ctrlq >= curQubit; --ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(phaseShift, ctrlq, curQubit - 1);

					//MatrixClass ps = phaseShift.getOperatorMatrix(QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits(), ctrlq, curQubit - 1);
					//m = ps * m;

					//std::cout << curQubit - 1 << " controls pi/" << div << " on " << ctrlq << std::endl;
					phase *= 2;
					//div /= 2;
				}

				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, curQubit - 1);
				
				//MatrixClass m1 = hadamard.getOperatorMatrix(QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits(), curQubit-1);
				//m = m1 * m;
			
				//std::cout << "Hadamard on " << curQubit - 1 << std::endl;
			}

			/*
			QC::SwapGate sg;
			m = sg.getOperatorMatrix(QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits(), 0, 1) * m;
			std::cout << m << std::endl;

			VectorClass v = VectorClass::Ones(QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrBasisStates());
			std::cout << std::endl << m * 0.5 * v << std::endl;
			exit(0);
			*/
		}

		void IQFT(unsigned int sq, unsigned int eq)
		{
			for (unsigned int curQubit = sq; curQubit < eq; ++curQubit)
			{
				QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, curQubit);

				//std::cout << "Hadamard on " << curQubit << std::endl;

				double phase = M_PI_2;
				//int div = 2;
				for (unsigned int ctrlq = curQubit + 1; ctrlq <= eq; ++ctrlq)
				{
					phaseShift.SetPhaseShift(phase);
					QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(phaseShift, curQubit, ctrlq);

					//std::cout << curQubit << " controls pi/" << div << " on " << ctrlq << std::endl;
					phase /= 2;
					//div *= 2;
				}
			}
			QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, eq);

			//std::cout << "Hadamard on " << eq << std::endl;
			//exit(0);
		}

		HadamardGate<MatrixClass> hadamard;
		ControlledPhaseShiftGate<MatrixClass> phaseShift;
		SwapGate<MatrixClass> swapOp;

		unsigned int sQubit;
		unsigned int eQubit;
	};

}


