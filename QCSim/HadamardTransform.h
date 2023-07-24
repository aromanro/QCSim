#pragma once

#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class GateClass = QC::Gates::HadamardGate<MatrixClass>> class OneQubitGateTransform : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			OneQubitGateTransform(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
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


