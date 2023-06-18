#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Distributed {

	// control = qubit 0, controlled = qubit 3
	// qubits 0 and 1 are in one location, qubits 2 and 3 are in another

	// qubits 1 and 2 are cat-entangled, then communication is classical, each party measures its qubit (from 1 and 2) and sends to the other the result

	// for a detailed description (along with more general cases, like a distributed controlled-U) see "Generalized GHZ States and Distributed Quantum Computing"
	// https://arxiv.org/abs/quant-ph/0402148v3

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class DistributedCNOT : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		DistributedCNOT(int addseed = 0)
			: BaseClass(4, addseed)
		{
			BaseClass::setToBasisState(0);
		}


		// could be implemented using subalgorithms, see CatEntangler, CatDisentangler and GeneralizedEntanglingGate
		unsigned int Execute() override
		{
			// the two gates make up an 'entangling gate', sometimes noted by E or E2
			BaseClass::ApplyGate(hadamard, 1);
			BaseClass::ApplyGate(cnot, 2, 1);

			BaseClass::ApplyGate(cnot, 1, 0);

			const unsigned int qubit1measurement = BaseClass::Measure(1); // measured and sent to the other
			if (qubit1measurement) // the other applies x conditionally on the received measurement
				BaseClass::ApplyGate(x, 2);

			BaseClass::ApplyGate(cnot, 3, 2);

			// cat-disentangler follows
			BaseClass::ApplyGate(hadamard, 2);

			const unsigned int qubit2measurement = BaseClass::Measure(2); // measured and sent to the other
			if (qubit2measurement) // the other applies z conditionally on the received measurement
				BaseClass::ApplyGate(z, 0);

			// return the result of the measurements
			return (qubit2measurement << 1) | qubit1measurement;
		}

	protected:
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::PauliZGate<MatrixClass> z;
	};


}


