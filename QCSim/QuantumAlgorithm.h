#pragma once

#include "QubitRegisterDebug.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumAlgorithm
	{
	public:
		using RegisterClass = QubitRegisterDebug<VectorClass, MatrixClass>;
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;

		QuantumAlgorithm(size_t N = 3, unsigned int addseed = 0)
			: reg(N, addseed)
		{
			assert(N > 0);
		}

		virtual ~QuantumAlgorithm() {};

		virtual size_t Execute() = 0;

		size_t getNrQubits() const { return reg.getNrQubits(); };
		size_t getNrBasisStates() const { return reg.getNrBasisStates(); };

		std::complex<double> getBasisStateAmplitude(size_t State) const
		{
			return reg.getBasisStateAmplitude(State);
		}

		double getBasisStateProbability(size_t State) const
		{
			return reg.getBasisStateProbability(State);
		}

		void setToBasisState(size_t State)
		{
			reg.setToBasisState(State);
		}

		void setToQubitState(size_t q)
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
		void setRawAmplitude(size_t State, std::complex<double> val)
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

		static size_t getQubitState(size_t q)
		{
			return 1ULL << q;
		}

		const RegisterClass& getRegister() const
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

		void ApplyGate(const GateClass& gate, size_t qubit, size_t controllingQubit = 0, size_t controllingQubit2 = 0)
		{
			reg.ApplyGate(gate, qubit, controllingQubit, controllingQubit2);
		}

		void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate)
		{
			reg.ApplyGate(gate);
		}

		void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates)
		{
			reg.ApplyGates(gates);
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			reg.ApplyOperatorMatrix(m);
		}

		std::map<size_t, size_t> RepeatedMeasure(size_t nrTimes = 1000)
		{
			return reg.RepeatedMeasure(nrTimes);
		}

		std::map<size_t, size_t> RepeatedMeasure(size_t firstQubit, size_t secondQubit, size_t nrTimes = 1000)
		{
			return reg.RepeatedMeasure(firstQubit, secondQubit, nrTimes);
		}

		size_t Measure()
		{
			return reg.MeasureAll();
		}

		size_t Measure(size_t qubit)
		{
			return reg.MeasureQubit(qubit);
		}
		
		size_t Measure(size_t firstQubit, size_t secondQubit)
		{
			return reg.Measure(firstQubit, secondQubit);
		}

		// to check how well the computed state matches some 'exact' known one
		double stateFidelity(const VectorClass& state) const
		{
			return reg.stateFidelity(state);
		}

		static void AdjustPhaseAndNormalize(VectorClass& regVals)
		{
			std::complex<double> v = abs(regVals[0]) > 1E-5 ? regVals[0] : std::complex<double>(1, 0);
			std::complex<double> accum(0, 0);
			const size_t nrBasisStates = static_cast<size_t>(regVals.size());
			for (size_t i = 0; i < nrBasisStates; ++i)
			{
				regVals[i] /= v;

				accum += norm(regVals[i]);
			}
			const double norm = 1. / sqrt(accum.real());
			for (size_t i = 0; i < nrBasisStates; ++i)
				regVals[i] *= norm;
		}

		bool writeToFile(const std::string& name, bool amplitude = true, bool append = false) const
		{
			return reg.writeToFile(name, amplitude, append);
		}

		void displayRegister() const
		{
			reg.displayRegister();
		}

		void ComputeStart()
		{
			reg.ComputeStart();
		}

		void ComputeEnd()
		{
			reg.ComputeEnd();
		}

		void ComputeClear()
		{
			reg.ComputeClear();
		}

		// applies again the recorded gates
		// with this the same operations can be repeated several times
		void Compute()
		{
			reg.Compute();
		}

		// undoes the recorded gates
		// the operations are unitary, so U^-1 = U^t and (U1 * U2)^t = U2^t * U1^t 
		void Uncompute()
		{
			reg.Uncompute();
		}

	protected:
		RegisterClass reg;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumSubAlgorithm
	{
	public:
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		QuantumSubAlgorithm()
		{
		}
		
		virtual ~QuantumSubAlgorithm() {};

		virtual size_t Execute(RegisterClass& reg) = 0;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumSubAlgorithmOnSubregister : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumSubAlgorithm<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		QuantumSubAlgorithmOnSubregister(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX)
			: sQubit(startQubit), eQubit(std::max(startQubit, std::min(N - 1, endQubit)))
		{
		}

		size_t getStartQubit() const { return sQubit; };
		size_t getEndQubit() const { return eQubit; };

	protected:
		size_t sQubit;
		size_t eQubit;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumSubAlgorithmOnSubregisterWithAncilla : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		QuantumSubAlgorithmOnSubregisterWithAncilla(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX, size_t startAncila = INT_MAX)
			: BaseClass(N, startQubit, endQubit), sAncilla(startAncila)
		{
		}

		size_t getStartAncilla() const { return sAncilla; };

	protected:
		size_t sAncilla;
	};
}
