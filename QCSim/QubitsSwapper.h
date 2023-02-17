#pragma once

#include "QubitRegister.h"
#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace QC {

	// maybe it could be used in some other cases
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitsSwapper : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		typedef QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass> BaseClass;
		typedef QubitRegister<VectorClass, MatrixClass> RegisterClass;

		QubitsSwapper(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: BaseClass(N, startQubit, endQubit)
		{
		}

		void Swap(RegisterClass& reg) const
		{
			unsigned int startQubit = BaseClass::getStartQubit();
			unsigned int endQubit = BaseClass::getEndQubit();

			while (startQubit < endQubit)
			{
				reg.ApplyGate(swapOp, startQubit, endQubit);
				++startQubit;
				--endQubit;
			}
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			Swap(reg);

			return reg.Measure();
		}

	protected:
		Gates::SwapGate<MatrixClass> swapOp;
	};

}


