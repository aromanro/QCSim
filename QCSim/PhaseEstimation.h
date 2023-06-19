#pragma once

#include "QuantumFourierTransform.h"
#include "Function.h"
#include "NQubitsControlledQuantumGate.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorPhaseEstimation : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		ShorPhaseEstimation(QC::Function<VectorClass, MatrixClass>& f, unsigned int N = 7, unsigned int L = 3, int addseed = 0)
			: fRegisterStartQubit(L), fourier(N, 0, L - 1), func(f)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
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

			// then perform an inverse fourier transform
			fourier.IQFT(reg);

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




	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PhaseEstimation : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		PhaseEstimation(const MatrixClass& op, unsigned int N = 7, unsigned int L = 3, int addseed = 0)
			: fRegisterStartQubit(L), nrQubits(N), fourier(N, 0, L - 1), U(op)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			// TODO: check if things are set up all right: size of U, size of reg, etc.
			
			// apply hadamard over each qubit from the x-register
			// reuse the hadamard gate from the fourier transform base class
			for (unsigned int i = 0; i < fRegisterStartQubit; ++i)
				reg.ApplyGate(fourier.hadamard, i);

			
			MatrixClass controlledGate = U;
			const unsigned int lastQubit = fRegisterStartQubit - 1;
			for (unsigned int ctrlQubit = 0; ctrlQubit < lastQubit; ++ctrlQubit)
			{
				NQubitsControlledQuantumGate<VectorClass, MatrixClass> UGate(nrQubits, controlledGate, fRegisterStartQubit, ctrlQubit);
				
				UGate.Execute(reg);

				// power up the controlled gate
				controlledGate *= controlledGate;
			}

			{
				NQubitsControlledQuantumGate<VectorClass, MatrixClass> UGate(nrQubits, controlledGate, fRegisterStartQubit, lastQubit);

				UGate.Execute(reg);
			}

			// it doesn't really matter if you measure the qubits from f and when you do after the above
			// or if you measure them several times in a row
			//QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(fRegisterStartQubit, QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);

			// then perform an inverse fourier transform
			fourier.IQFT(reg);

			// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
			// the f should be measured separately to find out its content

			//return reg.Measure();
			return reg.Measure(0, fRegisterStartQubit - 1);
		}

		unsigned int getFunctionStartQubit() const
		{
			return fRegisterStartQubit;
		}

		unsigned int getNrQubits() const
		{
			return nrQubits;
		}

		const MatrixClass& getU() const
		{
			return U;
		}

	protected:
		unsigned int fRegisterStartQubit;
		unsigned int nrQubits;

		QuantumFourierTransform<VectorClass, MatrixClass> fourier;
		MatrixClass U;
	};


}

