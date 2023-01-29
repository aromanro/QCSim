#pragma once

#include "QubitRegister.h"
#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace QC {

	// maybe it could be used in some other cases
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitsSwapper : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		QubitsSwapper(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>(N, startQubit, endQubit)
		{
		}


		void Swap(QubitRegister<VectorClass, MatrixClass>& reg) const
		{
			unsigned int startQubit = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>::getStartQubit();
			unsigned int endQubit = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>::getEndQubit();

			while (startQubit < endQubit)
			{
				reg.ApplyGate(swapOp, startQubit, endQubit);
				++startQubit;
				--endQubit;
			}
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			Swap(reg);

			return reg.Measure();
		}

	protected:
		Gates::SwapGate<MatrixClass> swapOp;
	};

}


