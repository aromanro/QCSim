#pragma once

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

namespace Adders {

	// to be used to implement a full quantum adder
	// not tested yet
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TwoQubitsHalfAdder : public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumSubAlgorithm<VectorClass, MatrixClass>;
		using RegisterClass = QC::QubitRegister<VectorClass, MatrixClass>;


		TwoQubitsHalfAdder(unsigned int qubit1, unsigned int qubit2, unsigned int qubitaux)
			: q1(qubit1), q2(qubit2), aux(qubitaux)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
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
		using BaseClass = TwoQubitsHalfAdder<VectorClass, MatrixClass>;
		using RegisterClass = QC::QubitRegister<VectorClass, MatrixClass>;

		TwoQubitsFullAdder(unsigned int qubit1, unsigned int qubit2, unsigned int cin, unsigned int qubitaux)
			: BaseClass(qubit1, qubit2, qubitaux), ci(cin),
			halfAdder(qubit2, cin, qubitaux)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			const unsigned int nrQubits = reg.getNrQubits();
			if (BaseClass::aux >= nrQubits || BaseClass::q2 >= nrQubits || BaseClass::q1 >= nrQubits || ci >= nrQubits)
				return 1; // error code
			else if (BaseClass::q1 == BaseClass::q2 || BaseClass::q1 == BaseClass::aux ||
				BaseClass::q2 == BaseClass::aux || BaseClass::aux == ci || BaseClass::q1 == ci || BaseClass::q2 == ci)
				return 2;

			BaseClass::Execute(reg);

			halfAdder.Execute(reg);
			reg.ApplyGate(BaseClass::cnot, BaseClass::q2, BaseClass::q1);

			return 0;
		}

	protected:
		unsigned int ci;

		BaseClass halfAdder;
	};

	// could be made more general, with N1 and N2 qubits for the two added numbers and even specifying the qubits location in the registry
	// but this is an example and doing that would make the code too complex without adding much
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NQubitsAdderAlgorithm : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		NQubitsAdderAlgorithm(unsigned int N = 3, int addseed = 0)
			: BaseClass(3*N + 1, addseed) // the qubits for the two N inputs, N outputs and one for carry
		{
			BaseClass::setToBasisState(0);
		}

		unsigned int Execute() override
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();
			assert(nrQubits >= 4);

			const unsigned int nrQubitsNumber = (nrQubits - 1) / 3;

			unsigned int n1q = 0;
			unsigned int n2q = nrQubitsNumber;
			unsigned int ci = 2 * nrQubitsNumber;

			while (n1q < nrQubitsNumber)
			{
				TwoQubitsFullAdder fullAdder(n1q, n2q, ci, ci + 1);
				fullAdder.Execute(BaseClass::reg);

				++n1q;
				++n2q;
				++ci;
			}

			return BaseClass::Measure();
		}
	};

}
