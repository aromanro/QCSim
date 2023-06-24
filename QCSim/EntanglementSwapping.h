#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Teleportation.h"
#include "Utils.h"

namespace Teleportation {


	// Alice and Bob have qubits 0 and 1 and 2 and 3 respectively
	// 4 and 5 are in some other place(s)

	// 0 and 2 are entangled - this entanglement is going to be swapped
	
	// for teleportation:
	// qubit 1 is entangled with 4, which can be in some other place, let's call it Charlie
	// qubit 3 is entangled with 5, which can be in some other place, let's call it Diane

	// by teleporting qubits 0 from Alice to Charlie and 2 from Bob to Diane
	// the entanglement is not shared between Alice and Bob anymore, but between Charlie and Diane

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class EntanglementSwapping :
		public QuantumTeleportation<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumTeleportation<VectorClass, MatrixClass>;
		using AlgorithmClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		EntanglementSwapping(int addseed = 0)
			: BaseClass(6, addseed)
		{
		}

		void SetBellStateToBeSent(bool s1 = false, bool s2 = false)
		{
			QC::BellState bellState;

			bellState.setBellState(BaseClass::BaseClass::reg, 0, 2, s1, s2);
		}

		unsigned int Teleport(bool explicitClassicalTransmission = false)
		{
			unsigned int res = BaseClass::Teleport(0, 1, 4, explicitClassicalTransmission);

			res |= BaseClass::Teleport(2, 3, 5, explicitClassicalTransmission) << 2;

			return res;
		}

		// called only if one wants the teleported qubits to be measured, otherwise just check the register contents
		unsigned int Execute() override
		{
			Teleport();

			return AlgorithmClass::Measure(4, 5);
		}
	};

}


