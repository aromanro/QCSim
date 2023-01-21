#pragma once

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace Adders {

	// to be used to implement a full quantum adder
	// not tested yet
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TwoQubitsHalfAdder : public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		TwoQubitsHalfAdder(unsigned int qubit1, unsigned int qubit2, unsigned int qubitaux)
			: q1(qubit1), q2(qubit2), aux(qubitaux)
		{
		}

		unsigned int Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			const unsigned int nrQubits = reg.getNrQubits();
			if (aux >= nrQubits || q2 >= nrQubits || q1 >= nrQubits) return 1; // error code
			else if (q1 == q2 || q1 == aux || q2 == aux) return 2;

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

		QC::Gates::ToffoliGate<MatrixClass> toffoli;
		QC::Gates::CNOTGate<MatrixClass> cnot;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TwoQubitsFullAdder : public TwoQubitsHalfAdder<VectorClass, MatrixClass>
	{
	public:
		TwoQubitsFullAdder(unsigned int qubit1, unsigned int qubit2, unsigned int cin, unsigned int qubitaux)
			: TwoQubitsHalfAdder<VectorClass, MatrixClass>(qubit1, qubit2, qubitaux), ci(cin),
			halfAdder(qubit2, cin, qubitaux)
		{
		}

		unsigned int Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			const unsigned int nrQubits = reg.getNrQubits();
			if (TwoQubitsHalfAdder<VectorClass, MatrixClass>::aux >= nrQubits || TwoQubitsHalfAdder<VectorClass, MatrixClass>::q2 >= nrQubits || TwoQubitsHalfAdder<VectorClass, MatrixClass>::q1 >= nrQubits || ci >= nrQubits) 
				return 1; // error code
			else if (TwoQubitsHalfAdder<VectorClass, MatrixClass>::q1 == TwoQubitsHalfAdder<VectorClass, MatrixClass>::q2 || TwoQubitsHalfAdder<VectorClass, MatrixClass>::q1 == TwoQubitsHalfAdder<VectorClass, MatrixClass>::aux || 
				TwoQubitsHalfAdder<VectorClass, MatrixClass>::q2 == TwoQubitsHalfAdder<VectorClass, MatrixClass>::aux || TwoQubitsHalfAdder<VectorClass, MatrixClass>::aux == ci || TwoQubitsHalfAdder<VectorClass, MatrixClass>::q1 == ci || TwoQubitsHalfAdder<VectorClass, MatrixClass>::q2 == ci) 
				return 2;

			TwoQubitsHalfAdder<VectorClass, MatrixClass>::Execute(reg);

			halfAdder.Execute(reg);
			reg.ApplyGate(TwoQubitsHalfAdder<VectorClass, MatrixClass>::cnot, TwoQubitsHalfAdder<VectorClass, MatrixClass>::q2, TwoQubitsHalfAdder<VectorClass, MatrixClass>::q1);

			return 0;
		}

	protected:
		unsigned int ci;

		TwoQubitsHalfAdder<VectorClass, MatrixClass> halfAdder;
	};

	// could be made more general, with N1 and N2 qubits for the two added numbers and even specifying the qubits location in the registry
	// but this is an example and doing that would make the code too complex without adding much
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NQubitsAdderAlgorithm : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		NQubitsAdderAlgorithm(unsigned int N = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(3*N + 1, addseed) // the qubits for the two N inputs, N outputs and one for carry
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
		}

		unsigned int Execute() override
		{
			const unsigned int nrQubits = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits();
			assert(nrQubits >= 4);

			const unsigned int nrQubitsNumber = (nrQubits - 1) / 3;

			unsigned int n1q = 0;
			unsigned int n2q = nrQubitsNumber;
			unsigned int ci = 2 * nrQubitsNumber;

			while (n1q < nrQubitsNumber)
			{
				TwoQubitsFullAdder fullAdder(n1q, n2q, ci, ci + 1);
				fullAdder.Execute(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);

				++n1q;
				++n2q;
				++ci;
			}

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}
	};

}
