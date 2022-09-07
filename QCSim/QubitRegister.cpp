#include "QubitRegister.h"

namespace QC {
	QubitRegister::QubitRegister(int N)
		: NrQubits(N), NrBasisStates(~(~0u << NrQubits))
	{
		registerStorage = Eigen::VectorXd::Zero(NrBasisStates);
	}

	void QubitRegister::setToBasisState(unsigned int State)
	{
		if (State >= NrBasisStates) return;

		registerStorage.setZero();
		registerStorage(State) = 1;
	}

	void QubitRegister::setToCatState()
	{
		registerStorage.setZero();
		static const double OneOverSqrt2 = 1. / sqrt(2.);

		registerStorage(0) = OneOverSqrt2;
		registerStorage(NrBasisStates - 1) = OneOverSqrt2;
	}
	
	void QubitRegister::setToEqualSuperposition()
	{
		registerStorage.setConstant(1. / sqrt(NrBasisStates));
	}

	void QubitRegister::Normalize()
	{
		const double accum = registerStorage.transpose() * registerStorage;

		registerStorage = 1. / sqrt(accum) * registerStorage;
	}
}
