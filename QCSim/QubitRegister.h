#pragma once

#include <math.h>
#include <Eigen/eigen>

namespace QC {

	class QubitRegister
	{
	public:
		QubitRegister(int N = 3);

		unsigned int getNrQubits() const { return NrQubits; };
		unsigned int getNrBasisStates() const { return NrBasisStates; };

		double getBasisStateAmplitude(unsigned int State) const {
			if (State >= NrBasisStates) return 0;

			return registerStorage(State);
		}

		void setToBasisState(unsigned int State);
		void setToCatState();
		void setToEqualSuperposition();
		void Normalize();

	protected:
		unsigned int NrQubits;
		unsigned int NrBasisStates;

		Eigen::VectorXd registerStorage;
	};

}

