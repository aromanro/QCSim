#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		// See "Generalized GHZ States and Distributed Quantum Computing"
		// https://arxiv.org/abs/quant-ph/0402148v3

		// particular case: two qubits is Bell encoder

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class CatEntangler : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			CatEntangler(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX)
				: BaseClass(N, startQubit, endQubit)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				const size_t startQubit = BaseClass::getStartQubit();
				const size_t endQubit = BaseClass::getEndQubit();

				reg.ApplyGate(hadamard, startQubit);

				for (size_t q = startQubit; q < endQubit; ++q)
					reg.ApplyGate(cnot, q + 1, q);


				return 0; // no measurement, so the return should be ignored
			}

		protected:
			QC::Gates::CNOTGate<MatrixClass> cnot;
			QC::Gates::HadamardGate<MatrixClass> hadamard;
		};

	}

}

