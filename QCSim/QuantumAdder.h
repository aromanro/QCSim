#pragma once

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace QC {

	// to be used to implement a full quantum adder
	// not tested yet
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TwoQubitsAdder : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		TwoQubitsAdder(unsigned int qubit1, unsigned int qubit2, unsigned int qubitaux)
			: q1(qubit1), q2(qubit2), aux(qubitaux)
		{
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			reg.ApplyGate(toffoli, aux, q1, q2);
			reg.ApplyGate(cnot, q2, q1);

			return 0;
		}

	protected:
		// q1, q2 - the qubits to be added, aux, the auxiliary that should be initialized to |0> and has the 'carry' functionality
		unsigned int q1;
		unsigned int q2;
		unsigned int aux;

		QC::ToffoliGate toffoli;
		QC::CNOTGate cnot;
	};

}
