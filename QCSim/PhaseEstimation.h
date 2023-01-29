#pragma once

#include "QuantumFourierTransform.h"
#include "Function.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PhaseEstimation : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		PhaseEstimation(QC::Function<VectorClass, MatrixClass>& f, unsigned int N = 7, unsigned int L = 3, int addseed = 0)
			: fRegisterStartQubit(L), fourier(N, 0, L - 1), func(f)
		{
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			// apply hadamard over each qubit from the x-register
			// reuse the hadamard gate from the fourier transform base class
			for (unsigned int i = 0; i < fRegisterStartQubit; ++i)
				reg.ApplyGate(fourier.hadamard, i);

			// now the f(x)
			func.Apply(reg);

			// it doesn't really matter if you measure the qubits from f and when you do after the above
			// or if you measure them several times in a row
			//QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(fRegisterStartQubit, QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);

			// then perform a fourier transform (inverse in the typical quantum literature definition of it)
			fourier.QFT(reg);

			// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
			// the f should be measured separately to find out its content

			//return reg.Measure(0, fRegisterStartQubit - 1);
			return reg.Measure();
		}

		unsigned int getFunctionStartQubit() const
		{
			return fRegisterStartQubit;
		}

	protected:
		unsigned int fRegisterStartQubit;

		QuantumFourierTransform<VectorClass, MatrixClass> fourier;
		QC::Function<VectorClass, MatrixClass>& func;
	};

}

