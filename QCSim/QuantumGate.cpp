#include "QuantumGate.h"

namespace QC {
	SingleQubitGate::SingleQubitGate()
	{
		operatorMat = Eigen::MatrixXcd::Zero(2, 2);
	}

	Eigen::MatrixXcd SingleQubitGate::getOperatorMatrix(unsigned int nrQubits, unsigned int qubit, unsigned int controllingQubit) const
	{
		const unsigned int nrBasisStates = 1u << nrQubits;
		Eigen::MatrixXcd extOperatorMat = Eigen::MatrixXcd::Zero(nrBasisStates, nrBasisStates);

		const unsigned int qubitBit = 1u << qubit;

		for (unsigned int i = 0; i < nrBasisStates; ++i)
		{ 
			const unsigned int ind1 = i | qubitBit;
			for (unsigned int j = 0; j < nrBasisStates; ++j)
				if (ind1 == (j | qubitBit))
					extOperatorMat(i, j) = operatorMat(i & qubitBit ? 1 : 0, j & qubitBit ? 1 : 0);
		}

		return extOperatorMat;
	}


	HadamardGate::HadamardGate()
	{
		static const double norm = 1. / sqrt(2.);
		operatorMat(0, 0) = norm;
		operatorMat(0, 1) = norm;
		operatorMat(1, 0) = norm;
		operatorMat(1, 1) = -norm;
	}

	PhaseShiftGate::PhaseShiftGate(double theta)
	{
		operatorMat(0, 0) = 1;
		SetPhaseShift(theta);
	}

	void PhaseShiftGate::SetPhaseShift(double theta)
	{
		operatorMat(1, 1) = exp(std::complex<double>(0, theta));
	}


	TwoQubitsControlledGate::TwoQubitsControlledGate()
	{
		operatorMat = Eigen::MatrixXcd::Identity(4, 4);
	}

	void TwoQubitsControlledGate::SetOperation(const Eigen::MatrixXcd& U)
	{
		assert(U.rows() == 2 && U.cols() == 2);
		operatorMat.block(2, 2, 2, 2) = U;
	}

	Eigen::MatrixXcd TwoQubitsControlledGate::getOperatorMatrix(unsigned int nrQubits, unsigned int qubit, unsigned int controllingQubit) const
	{
		const unsigned int nrBasisStates = 1u << nrQubits;
		Eigen::MatrixXcd extOperatorMat = Eigen::MatrixXcd::Zero(nrBasisStates, nrBasisStates);

		const unsigned int qubitBit = 1u << qubit;
		const unsigned int ctrlQubitBit = 1u << controllingQubit;

		for (unsigned int i = 0; i < nrBasisStates; ++i)
		{
			const unsigned int ind1 = i | qubitBit | ctrlQubitBit;
			for (unsigned int j = 0; j < nrBasisStates; ++j)
				if (ind1 == (j | qubitBit | ctrlQubitBit))
					extOperatorMat(i, j) = operatorMat((i & ctrlQubitBit ? 2 : 0) | (i & qubitBit ? 1 : 0), (j & ctrlQubitBit ? 2 : 0) | (j & qubitBit ? 1 : 0));
		}

		return extOperatorMat;
	}

	CNOTGate::CNOTGate()
	{
		Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(2, 2);

		U(0, 1) = 1;
		U(1, 0) = 1;

		SetOperation(U);
	}

}
