#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		// See "Generalized GHZ States and Distributed Quantum Computing"
		// https://arxiv.org/abs/quant-ph/0402148v3

		// particular case: two qubits is Bell decoder

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class CatDisentangler : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			CatDisentangler(size_t N, size_t restoredQubit = 0, size_t startQubit = 1, size_t endQubit = INT_MAX)
				: BaseClass(N, startQubit, endQubit), rQubit(restoredQubit)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				const size_t startQubit = BaseClass::getStartQubit();
				const size_t endQubit = BaseClass::getEndQubit();

				for (size_t q = startQubit; q <= endQubit; ++q)
					reg.ApplyGate(hadamard, q);

				const size_t measurement = reg.Measure(startQubit, endQubit);

				size_t r = measurement;
				size_t s = 0;
				while (r)
				{
					s += r & 1;
					r >>= 1;
				}

				if (s % 2)
					reg.ApplyGate(z, rQubit);

				return measurement;
			}

		protected:
			size_t rQubit;
			QC::Gates::HadamardGate<MatrixClass> hadamard;
			QC::Gates::PauliZGate<MatrixClass> z;
		};

	}
}


