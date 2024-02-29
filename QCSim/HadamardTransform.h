#pragma once

#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class GateClass = QC::Gates::HadamardGate<MatrixClass>> class OneQubitGateTransform : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			OneQubitGateTransform(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX)
				: BaseClass(N, startQubit, endQubit)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				ApplyHadamardOnAllQubits(reg);

				return 0;
			}

		protected:
			void ApplyHadamardOnAllQubits(RegisterClass& reg) const
			{
				for (size_t q = BaseClass::getStartQubit(); q <= BaseClass::getEndQubit(); ++q)
					reg.ApplyGate(gate, q);
			}

			GateClass gate;
		};

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class HadamardTransform : public OneQubitGateTransform<VectorClass, MatrixClass, QC::Gates::HadamardGate<MatrixClass>>
		{
		};

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class XTransform : public OneQubitGateTransform<VectorClass, MatrixClass, QC::Gates::PauliXGate<MatrixClass>>
		{
		};

	}

}


