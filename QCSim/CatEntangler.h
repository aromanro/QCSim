#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	// See "Generalized GHZ States and Distributed Quantum Computing"
	// https://arxiv.org/abs/quant-ph/0402148v3

	// particular case: two qubits is Bell encoder

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class CatEntangler : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		CatEntangler(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: BaseClass(N, startQubit, endQubit)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			const unsigned int startQubit = BaseClass::getStartQubit();
			const unsigned int endQubit = BaseClass::getEndQubit();

			reg.ApplyGate(hadamard, startQubit);
			
			for (unsigned int q = startQubit; q < endQubit; ++q)
				reg.ApplyGate(cnot, q + 1, q);


			return 0; // no measurement, so the return should be ignored
		}

	protected:
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};

}

