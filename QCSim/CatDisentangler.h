#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	// See "Generalized GHZ States and Distributed Quantum Computing"
	// https://arxiv.org/abs/quant-ph/0402148v3

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class CatDisentangler : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		CatDisentangler(unsigned int N, unsigned int restoredQubit = 0, unsigned int startQubit = 1, unsigned int endQubit = INT_MAX)
			: rQubit(restoredQubit), BaseClass(N, startQubit, endQubit)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			const unsigned int startQubit = BaseClass::getStartQubit();
			const unsigned int endQubit = BaseClass::getEndQubit();

			for (unsigned int q = startQubit; q <= endQubit; ++q)
				reg.ApplyGate(hadamard, m);

			const unsigned int measurement = reg.Measure(startQubit, endQubit);

			unsigned int r = measurement;
			unsigned int s = 0;
			while (r)
			{
				if (r & 1)
					++s;
				r >>= 1;
			}

			if (s % 2)
				reg.ApplyGate(z, rQubit);

			return measurement;
		}

	protected:
		unsigned int rQubit;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliZGate<MatrixClass> z;
	};

}


