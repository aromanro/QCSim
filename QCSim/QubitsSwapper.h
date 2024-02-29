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

			QubitsSwapper(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX)
				: BaseClass(N, startQubit, endQubit)
			{
			}

			void Swap(RegisterClass& reg) const
			{
				size_t startQubit = BaseClass::getStartQubit();
				size_t endQubit = BaseClass::getEndQubit();

				while (startQubit < endQubit)
				{
					reg.ApplyGate(swapOp, startQubit, endQubit);
					++startQubit;
					--endQubit;
				}
			}

			size_t Execute(RegisterClass& reg) override
			{
				Swap(reg);

				return reg.MeasureAll();
			}

		protected:
			Gates::SwapGate<MatrixClass> swapOp;
		};

	}

}


