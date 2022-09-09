#include "Oracle.h"

namespace Grover {

	Eigen::MatrixXcd Oracle::getOperatorMatrix(unsigned int nrQubits, unsigned int qubit, unsigned int controllingQubit) const
	{
		const unsigned int nrBasisStates = 1u << nrQubits;
		Eigen::MatrixXcd extOperatorMat = Eigen::MatrixXcd::Identity(nrBasisStates, nrBasisStates);

		if (correctQuestionState < nrBasisStates)
			extOperatorMat(correctQuestionState, correctQuestionState) = -1;

		return extOperatorMat;
	}

	Eigen::MatrixXcd J::getOperatorMatrix(unsigned int nrQubits, unsigned int qubit, unsigned int controllingQubit) const
	{
		const unsigned int nrBasisStates = 1u << nrQubits;
		Eigen::MatrixXcd extOperatorMat = Eigen::MatrixXcd::Identity(nrBasisStates, nrBasisStates);

		extOperatorMat(0, 0) = -1;

		return extOperatorMat;
	}
}
