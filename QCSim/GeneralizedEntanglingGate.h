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

			GeneralizedEntanglingGate(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
				: BaseClass(N, startQubit + 1, endQubit)
			{
			}

			unsigned int Execute(RegisterClass& reg) override
			{
				const unsigned int qubit = BaseClass::BaseClass::getStartQubit(); // the start qubit for base class was set to startQubit + 1
				const unsigned int endQubit = BaseClass::BaseClass::getEndQubit();

				reg.ApplyGate(BaseClass::cnot, qubit, qubit - 1);

				const unsigned int measQubit = reg.Measure(qubit);

				if (measQubit)
					for (unsigned int q = qubit + 1; q <= endQubit; ++q)
						reg.ApplyGate(BaseClass::x, q);

				return measQubit;
			}

		protected:
			QC::Gates::PauliXGate<MatrixClass> x;
		};

	}

}


