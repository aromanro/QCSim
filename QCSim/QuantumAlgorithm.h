#pragma once

#include "QubitRegister.h"

namespace QC {
	
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumAlgorithm
	{
	public:
		QuantumAlgorithm(unsigned int N = 3, int addseed = 0)
			: reg(N, addseed)
		{
		}

		virtual unsigned int Execute() = 0;

		unsigned int getNrQubits() const { return reg.getNrQubits(); };
		unsigned int getNrBasisStates() const { return reg.getNrBasisStates(); };

		std::complex<double> getBasisStateAmplitude(unsigned int State) const {
			return reg.getBasisStateAmplitude(State);
		}

		void setToBasisState(unsigned int State)
		{
			reg.setToBasisState(State);
		}

		void setToQubitState(unsigned int q)
		{
			reg.setToQubitState(q);
		}

		void setToCatState()
		{
			reg.setToCatState();
		}

		void setToEqualSuperposition()
		{
			reg.setToEqualSuperposition();
		}

		// to be able to set them all, after setting them, call Normalize
		void setRawAmplitude(unsigned int State, std::complex<double> val)
		{
			reg.setRawAmplitude(State, val);
		}

		void Normalize()
		{
			reg.Normalize();
		}

		static unsigned int getQubitState(unsigned int q)
		{
			return 1u << q;
		}

		const QubitRegister<VectorClass, MatrixClass>& getRegister() const
		{
			return reg;
		}

		VectorClass getRegisterStorage() const
		{
			return reg.getRegisterStorage();
		}

		void ApplyGate(const QuantumGate<MatrixClass>& gate, unsigned int qubit, unsigned int controllingQubit = 0)
		{
			reg.ApplyGate(gate, qubit, controllingQubit);
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			reg.ApplyOperatorMatrix(m);
		}

		unsigned int Measure()
		{
			return reg.Measure();
		}

		unsigned int Measure(unsigned int firstQubit, unsigned int secondQubit)
		{
			return reg.Measure(firstQubit, secondQubit);
		}

	protected:
		QubitRegister<VectorClass, MatrixClass> reg;
	};

}
