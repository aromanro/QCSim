#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Distributed {

	// control qubit = 5, target qubit = 0

	// out 5->3, 0->2

	// see https://en.wikipedia.org/wiki/Quantum_gate_teleportation
	// also "Quantum Teleportation is a Universal Computational Primitive" by Gottesman and Chuang
	// https://arxiv.org/abs/quant-ph/9908010

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class TeleportedCNOT : public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		TeleportedCNOT(int addseed = 0)
			: BaseClass(6, addseed)
		{
			BaseClass::setToBasisState(0);
		}

		size_t Execute() override
		{
			Entangle(1, 2);
			Entangle(3, 4);
			BaseClass::ApplyGate(cnot, 2, 3);

			ApplyTeleportationCircuit(0, 1);
			ApplyTeleportationCircuit(5, 4);

			const int measurements1 = BaseClass::Measure(0, 1);
			const int measurements2 = BaseClass::Measure(4, 5);
			
			if (measurements2 & 1)
			{
				BaseClass::ApplyGate(x, 2);
				BaseClass::ApplyGate(x, 3);
			}
			if (measurements2 & 2) BaseClass::ApplyGate(z, 3);

			if (measurements1 & 2) BaseClass::ApplyGate(x, 2);
			if (measurements1 & 1)
			{
				BaseClass::ApplyGate(z, 2);
				BaseClass::ApplyGate(z, 3);
			}

			return measurements1 | (measurements2 << 4); // leave a gap for the output qubits, use masks to extract them
		}

	protected:
		void Entangle(int q1, int q2)
		{
			// the two gates make up an 'entangling gate', sometimes noted by E or E2
			BaseClass::ApplyGate(hadamard, q1);
			BaseClass::ApplyGate(cnot, q2, q1);
		}

		void ApplyTeleportationCircuit(size_t q1, size_t q2)
		{
			BaseClass::ApplyGate(cnot, q2, q1);
			BaseClass::ApplyGate(hadamard, q1);
		}

		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::PauliZGate<MatrixClass> z;
	};


}

