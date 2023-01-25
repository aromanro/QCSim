#pragma once

#include "QubitRegister.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumAlgorithm
	{
	public:
		QuantumAlgorithm(unsigned int N = 3, int addseed = 0)
			: reg(N, addseed)
		{
			assert(N > 0);
		}

		virtual ~QuantumAlgorithm() {};

		virtual unsigned int Execute() = 0;

		unsigned int getNrQubits() const { return reg.getNrQubits(); };
		unsigned int getNrBasisStates() const { return reg.getNrBasisStates(); };

		std::complex<double> getBasisStateAmplitude(unsigned int State) const
		{
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

		void Clear()
		{
			reg.Clear();
		}

		void Normalize()
		{
			reg.Normalize();
		}

		void AdjustPhaseAndNormalize()
		{
			reg.AdjustPhaseAndNormalize();
		}

		static unsigned int getQubitState(unsigned int q)
		{
			return 1u << q;
		}

		const QubitRegister<VectorClass, MatrixClass>& getRegister() const
		{
			return reg;
		}

		const VectorClass& getRegisterStorage() const
		{
			return reg.getRegisterStorage();
		}

		void setRegisterStorage(const VectorClass& vals)
		{
			reg.setRegisterStorage(vals);
		}

		void ApplyGate(const Gates::QuantumGate<MatrixClass>& gate, unsigned int qubit, unsigned int controllingQubit = 0, unsigned int controllingQubit2 = 0)
		{
			reg.ApplyGate(gate, qubit, controllingQubit, controllingQubit2);
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			reg.ApplyOperatorMatrix(m);
		}

		std::map<unsigned int, unsigned int> RepeatedMeasure(unsigned int nrTimes = 1000)
		{
			return reg.RepeatedMeasure(nrTimes);
		}

		std::pair<unsigned int, unsigned int> RepeatedMeasure(unsigned int qubit, unsigned int nrTimes = 1000)
		{
			return reg.RepeatedMeasure(qubit, nrTimes);
		}

		std::map<unsigned int, unsigned int> RepeatedMeasure(unsigned int firstQubit, unsigned int secondQubit, unsigned int nrTimes = 1000)
		{
			return reg.RepeatedMeasure(firstQubit, secondQubit, nrTimes);
		}

		unsigned int Measure()
		{
			return reg.Measure();
		}

		unsigned int Measure(unsigned int qubit)
		{
			return reg.Measure(qubit);
		}
		
		unsigned int Measure(unsigned int firstQubit, unsigned int secondQubit)
		{
			return reg.Measure(firstQubit, secondQubit);
		}

		// to check how well the computed state matches some 'exact' known one
		double stateFidelity(const VectorClass& state) const
		{
			return reg.stateFidelity(state);
		}

		static void AdjustPhaseAndNormalize(Eigen::VectorXcd& regVals)
		{
			std::complex<double> v = abs(regVals[0]) > 1E-5 ? regVals[0] : std::complex<double>(1, 0);
			std::complex<double> accum(0, 0);
			const unsigned int nrBasisStates = static_cast<unsigned int>(regVals.size());
			for (unsigned int i = 0; i < nrBasisStates; ++i)
			{
				regVals[i] /= v;

				accum += regVals[i] * std::conj(regVals[i]);
			}
			const double norm = 1. / sqrt(accum.real());
			for (unsigned int i = 0; i < nrBasisStates; ++i)
				regVals[i] *= norm;
		}

	protected:
		QubitRegister<VectorClass, MatrixClass> reg;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumSubAlgorithm
	{
	public:
		QuantumSubAlgorithm()
		{
		}
		
		virtual ~QuantumSubAlgorithm() {};

		virtual unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) = 0;
	};

}
