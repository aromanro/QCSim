#pragma once

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace QC {

	// to be used to implement a full quantum adder
	// not tested yet
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TwoQubitsHalfAdder : public QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		TwoQubitsHalfAdder(unsigned int qubit1, unsigned int qubit2, unsigned int qubitaux)
			: q1(qubit1), q2(qubit2), aux(qubitaux)
		{
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			// TODO: check if the registry is big enough
			reg.ApplyGate(toffoli, aux, q1, q2);
			reg.ApplyGate(cnot, q2, q1);

			return 0;
		}

	protected:
		// q1, q2 - the qubits to be added, aux, the auxiliary that should be initialized to |0> and has the 'carry' functionality
		// 'sum' ends up in q2, carry in aux
		unsigned int q1;
		unsigned int q2;
		unsigned int aux;

		QC::ToffoliGate toffoli;
		QC::CNOTGate cnot;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TwoQubitsFullAdder : public QuantumHalfAdder<VectorClass, MatrixClass>
	{
	public:
		TwoQubitsFullAdderAdder(unsigned int qubit1, unsigned int qubit2, unsigned int cin, unsigned int qubitaux)
			: QuantumHalfAdder<VectorClass, MatrixClass>(qubit1, qubit2, qubitaux), ci(cin),
			halfAdder(qubit2, cin, aux)
		{
		}

		unsigned int Execute(QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			// TODO: check if the registry is big enough
			QuantumHalfAdder<VectorClass, MatrixClass>::Execute(reg);

			halfAdder.Execute(reg);
			reg.ApplyGate(QuantumHalfAdder<VectorClass, MatrixClass>::cnot, qubit2, qubit1);

			return 0;
		}

	protected:
		unsigned int ci;

		QuantumHalfAdder<VectorClass, MatrixClass> halfAdder;
	};
}
