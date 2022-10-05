#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace Coding
{
	// there are two qubits, they start in state |00> then they are entangled using a Hadamard and a CNOT gate
	// the first one is kept to encode the classical bits and then sent, the other one is given in advance to the receiver

	// then depending on the classical bits to be sent, the sender (Alice) interacts with her qubit using a not and z gate
	// then the qubit is sent to Bob, where it is decoded, ending with a measurement that gets back the classical bits

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SuperdenseCoding :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		SuperdenseCoding(int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(2, addseed), b1(false), b2(false)
		{
		}

		void Send(bool bit1, bool bit2)
		{
			Init();

			Encode(bit1, bit2);
		}

		unsigned int Receive()
		{
			Decode();

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.Measure();
		}

		// needed in case Execute is used, that one will call Send and Receive one after another
		void SetBits(bool bit1, bool bit2)
		{
			b1 = bit1;
			b2 = bit2;
		}

		unsigned int Execute() override
		{
			Send(b1, b2);

			return Receive();
		}

	protected:
		void Init()
		{
			// prepare the entangled pair
			// starting from |00> (default) gets to (|00> + |11>)/sqrt(2)
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(cnot, 1, 0);
		}

		void Decode()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(cnot, 1, 0);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, 0);
		}

		void Encode(bool bit1, bool bit2)
		{
			// the resemblance with the 'explicit' teleportation code is not a coincidence, teleportation and superdense coding are sort of opposite to each other
			if (bit2) QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(x, 0);
			if (bit1) QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(z, 0);
		}


		// gates used to obtain the entangled qubits and then at the end for decoding
		QC::HadamardGate<MatrixClass> hadamard;
		QC::CNOTGate<MatrixClass> cnot;

		// single qubit gates to act on the transmission qubit (depending on the classical bits to be transmitted)
		QC::PauliXGate<MatrixClass> x;
		QC::PauliZGate<MatrixClass> z;

		// store the bits to send for 'Execute'
		bool b1;
		bool b2;
	};
}

