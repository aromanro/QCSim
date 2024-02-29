#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#include "Tests.h"


namespace Distributed {

	// control = qubit 0, controlled = qubit 3
	// qubits 0 and 1 are in one location, qubits 2 and 3 are in another

	// qubits 1 and 2 are cat-entangled, then communication is classical, each party measures its qubit (from 1 and 2) and sends to the other the result

	// for a detailed description (along with more general cases, like a distributed controlled-U) see "Generalized GHZ States and Distributed Quantum Computing"
	// https://arxiv.org/abs/quant-ph/0402148v3


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd, class ControlledUClass = QC::Gates::CNOTGate<MatrixClass>> class DistributedCU : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DistributedCU(size_t nrQubits = 4, size_t ctrlq = 0, size_t targetq = 3, size_t entq1 = 1, size_t entq2 = 2, int addseed = 0)
			: BaseClass(nrQubits, addseed), ctrlQubit(ctrlq), targetQubit(targetq), entQubit1(entq1), entQubit2(entq2)
		{
			BaseClass::setToBasisState(0);
		}


		// could be implemented using subalgorithms, see CatEntangler, CatDisentangler and GeneralizedEntanglingGate
		size_t Execute() override
		{
			// the two gates make up an 'entangling gate', sometimes noted by E or E2
			BaseClass::ApplyGate(hadamard, entQubit1);
			BaseClass::ApplyGate(cnot, entQubit2, entQubit1);

			// teleport the control qubit

			BaseClass::ApplyGate(cnot, entQubit1, ctrlQubit);

			const size_t qubit1measurement = BaseClass::Measure(entQubit1); // measured and sent to the other
			if (qubit1measurement) // the other applies x conditionally on the received measurement
				BaseClass::ApplyGate(x, entQubit2);

			// use teleported qubit to control the controlled-U
			BaseClass::ApplyGate(cU, targetQubit, entQubit2);

			// cat-disentangler follows
			BaseClass::ApplyGate(hadamard, entQubit2);

			const size_t qubit2measurement = BaseClass::Measure(entQubit2); // measured and sent to the other
			if (qubit2measurement) // the other applies z conditionally on the received measurement
				BaseClass::ApplyGate(z, ctrlQubit);

			// return the result of the measurements
			return (qubit2measurement << 1) | qubit1measurement;
		}

	protected:
		size_t ctrlQubit;
		size_t targetQubit;
		size_t entQubit1;
		size_t entQubit2;
		QC::Gates::CNOTGate<MatrixClass> cnot;
		ControlledUClass cU;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::PauliZGate<MatrixClass> z;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DistributedCNOT : public DistributedCU<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = DistributedCU<VectorClass, MatrixClass>;

		DistributedCNOT(int addseed = 0)
			: BaseClass(4, 0, 3, 1, 2, addseed)
		{
		}
	};


}


