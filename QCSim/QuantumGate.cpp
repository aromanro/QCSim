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
			for (unsigned int j = 0; j < nrBasisStates; ++j)
			{
				// force the qubit bit to 1
				const unsigned int ind1 = i | qubitBit;
				const unsigned int ind2 = j | qubitBit;

				if (ind1 == ind2)
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
}
