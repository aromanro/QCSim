#pragma once

#include "QubitRegister.h"
#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		// maybe it could be used in some other cases
		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitsSwapper : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

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

				return reg.MeasureAll();
			}

		protected:
			Gates::SwapGate<MatrixClass> swapOp;
		};

	}

}


