#pragma once


#include "CatEntangler.h"

namespace QC {

	namespace SubAlgo {

		// See "Generalized GHZ States and Distributed Quantum Computing"
		// https://arxiv.org/abs/quant-ph/0402148v3

		// Also noted with En

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class GeneralizedEntanglingGate : public CatEntangler<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = CatEntangler<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			GeneralizedEntanglingGate(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX)
				: BaseClass(N, startQubit + 1, endQubit)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				const size_t qubit = BaseClass::BaseClass::getStartQubit(); // the start qubit for base class was set to startQubit + 1
				const size_t endQubit = BaseClass::BaseClass::getEndQubit();

				reg.ApplyGate(BaseClass::cnot, qubit, qubit - 1);

				const size_t measQubit = reg.Measure(qubit);

				if (measQubit)
					for (size_t q = qubit + 1; q <= endQubit; ++q)
						reg.ApplyGate(BaseClass::x, q);

				return measQubit;
			}

		protected:
			QC::Gates::PauliXGate<MatrixClass> x;
		};

	}

}


