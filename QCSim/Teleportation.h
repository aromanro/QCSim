#pragma once
#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace Teleportation
{

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumTeleportation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		QuantumTeleportation(unsigned int N = 3, int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed)
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
		}
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumTeleportationRealization : public QuantumTeleportation<VectorClass, MatrixClass>
	{
	public:
		QuantumTeleportationRealization(int addseed = 0)
			: QuantumTeleportation<VectorClass, MatrixClass>(3, addseed)
		{
		}

		void SetState(std::complex<double> alpha, std::complex<double> beta)
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Clear();
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setRawAmplitude(0, alpha);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setRawAmplitude(1, beta);

			// ensure it's normalized:
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Normalize();
		}

		// considers the initial state already set
		// the qubits to be entangled start as |00> and they are entangled to a EPR state using cnot and hadamard gates 
		// the first two qubits are on the Alice side, the third belongs to Bob
		// returns the values of the two measured qubits
		unsigned int Teleport(bool explicitClassicalTransmission = false)
		{
			// make an EPR pair out of the second and third qubit:

			// starting from |00> (default) gets to (|00> + |11>)/sqrt(2)
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, 1);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 2, 1);

			// interacting the sending qubit (0) with the Alice's side of the entangled pair:
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, 0);

			// measurements can be actually done all at the end, but let's pretend
			const unsigned int measuredValues = QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(0, 1);

			// now that they are measured they go to Bob, which uses the values to act on its entangled qubit:

			if (explicitClassicalTransmission)
			{
				const bool firstQubitMeasurement = (measuredValues & 0x1) != 0;
				const bool secondQubitMeasurement = (measuredValues & 0x2) != 0;

				if (secondQubitMeasurement) QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(x, 2);
				if (firstQubitMeasurement) QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(z, 2);
			}
			else
			{
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 2, 1);
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cz, 2, 0);
			}

			return measuredValues;
		}

		// called only if one wants the teleported qubit to be measured, otherwise just check the register contents
		unsigned int Execute() override
		{
			Teleport();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(2, 2); // teleported qubit
		}

	protected:
		// the initial state to be teleported can be set directly or set up from something simple using some gate (hadamard, phase shift)
		// needed to create the Bell states
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::CNOTGate<MatrixClass> cnot; // also needed on Alice side
		// needed on receiving side:
		QC::Gates::ControlledZGate<MatrixClass> cz;

		// needed on Bob side, if the classical transmission is explicit (simple one qubit gates that are applied or not depending on the value of the classical bit):
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::PauliZGate<MatrixClass> z;
	};

}
