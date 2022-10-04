#pragma once

#include <Eigen/eigen>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	template<class MatrixClass = Eigen::MatrixXcd> class QuantumGate
	{
	public:
		virtual MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const = 0;
	};

	template<class MatrixClass = Eigen::MatrixXcd> class QuantumGateWithOp : public QuantumGate<MatrixClass>
	{
	protected:
		MatrixClass operatorMat;
	};

	template<class MatrixClass = Eigen::MatrixXcd> class SingleQubitGate : public QuantumGateWithOp<MatrixClass>
	{
	public:
		SingleQubitGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat = MatrixClass::Zero(2, 2);
		}

		// controllingQubit is ignored, it will be used for two qubit gates
		MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override
		{
			const unsigned int nrBasisStates = 1u << nrQubits;
			MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			const unsigned int qubitBit = 1u << qubit;

			for (unsigned int i = 0; i < nrBasisStates; ++i)
			{
				const unsigned int ind1 = i | qubitBit;
				for (unsigned int j = 0; j < nrBasisStates; ++j)
					if (ind1 == (j | qubitBit))
						extOperatorMat(i, j) = QuantumGateWithOp<MatrixClass>::operatorMat(i & qubitBit ? 1 : 0, j & qubitBit ? 1 : 0);
			}

			return extOperatorMat;
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class HadamardGate : public SingleQubitGate<MatrixClass>
	{
	public:
		HadamardGate()
		{
			static const double norm = 1. / sqrt(2.);
			QuantumGateWithOp<MatrixClass>::operatorMat(0, 0) = norm;
			QuantumGateWithOp<MatrixClass>::operatorMat(0, 1) = norm;
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 0) = norm;
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 1) = -norm;
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class PhaseShiftGate : public SingleQubitGate<MatrixClass>
	{
	public:
		PhaseShiftGate(double theta = 0)
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(0, 0) = 1;
			SetPhaseShift(theta);
		}

		void SetPhaseShift(double theta)
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 1) = exp(std::complex<double>(0, theta));
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class PauliXGate : public SingleQubitGate<MatrixClass>
	{
	public:
		PauliXGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(0, 1) = 1;
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 0) = 1;
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class PauliYGate : public SingleQubitGate<MatrixClass>
	{
	public:
		PauliYGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(0, 1) = std::complex<double>(0, -1);
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 0) = std::complex<double>(0, 1);
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class PauliZGate : public SingleQubitGate<MatrixClass>
	{
	public:
		PauliZGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(0, 0) = 1;
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 1) = -1;
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class TwoQubitsGate : public QuantumGateWithOp<MatrixClass>
	{
	public:
		TwoQubitsGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat = MatrixClass::Identity(4, 4);
		}

		MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override
		{
			const unsigned int nrBasisStates = 1u << nrQubits;
			MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			const unsigned int qubitBit = 1u << qubit;
			const unsigned int ctrlQubitBit = 1u << controllingQubit;

			for (unsigned int i = 0; i < nrBasisStates; ++i)
			{
				const unsigned int ind1 = i | qubitBit | ctrlQubitBit;
				for (unsigned int j = 0; j < nrBasisStates; ++j)
					if (ind1 == (j | qubitBit | ctrlQubitBit))
						extOperatorMat(i, j) = QuantumGateWithOp<MatrixClass>::operatorMat((i & ctrlQubitBit ? 2 : 0) | (i & qubitBit ? 1 : 0), (j & ctrlQubitBit ? 2 : 0) | (j & qubitBit ? 1 : 0));
			}

			return extOperatorMat;
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class SwapGate : public TwoQubitsGate<MatrixClass>
	{
	public:
		SwapGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 1) = 0;
			QuantumGateWithOp<MatrixClass>::operatorMat(2, 2) = 0;
			QuantumGateWithOp<MatrixClass>::operatorMat(1, 2) = 1;
			QuantumGateWithOp<MatrixClass>::operatorMat(2, 1) = 1;
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class TwoQubitsControlledGate : public TwoQubitsGate<MatrixClass>
	{
	public:
		TwoQubitsControlledGate() {};

		void SetOperation(const MatrixClass& U)
		{
			assert(U.rows() == 2 && U.cols() == 2);
			QuantumGateWithOp<MatrixClass>::operatorMat.block(2, 2, 2, 2) = U;
		}		
	};



	template<class MatrixClass = Eigen::MatrixXcd> class CNOTGate : public TwoQubitsControlledGate<MatrixClass>
	{
	public:
		CNOTGate()
		{
			MatrixClass U = MatrixClass::Zero(2, 2);

			U(0, 1) = 1;
			U(1, 0) = 1;

			TwoQubitsControlledGate<MatrixClass>::SetOperation(U);
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseShiftGate : public TwoQubitsControlledGate<MatrixClass>
	{
	public:
		ControlledPhaseShiftGate(double theta = 0)
		{
			SetPhaseShift(theta);
		}

		void SetPhaseShift(double theta)
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(3, 3) = exp(std::complex<double>(0, theta));
		}
	};

	template<class MatrixClass = Eigen::MatrixXcd> class ControlledZGate : public TwoQubitsControlledGate<MatrixClass>
	{
	public:
		ControlledZGate()
		{
			QuantumGateWithOp<MatrixClass>::operatorMat(3, 3) = -1;
		}
	};
}


