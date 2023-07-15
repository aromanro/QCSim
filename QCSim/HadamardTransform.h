#pragma once

#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class HadamardTransform : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			HadamardTransform(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
				: BaseClass(N, startQubit, endQubit)
			{
			}

			unsigned int Execute(RegisterClass& reg) override
			{
				ApplyHadamardOnAllQubits(reg);

				return 0;
			}

		protected:
			void ApplyHadamardOnAllQubits(RegisterClass& reg) const
			{
				for (unsigned int q = BaseClass::getStartQubit(); q <= BaseClass::getEndQubit(); ++q)
					reg.ApplyGate(hadamard, i);
			}

			QC::Gates::HadamardGate<MatrixClass> hadamard;
		};

	}

}


