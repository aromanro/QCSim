#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

namespace Coding
{
	// there are two qubits, they start in state |00> then they are entangled using a Hadamard and a CNOT gate
	// the first one is kept to encode the classical bits and then sent, the other one is given in advance to the receiver

	// then depending on the classical bits to be sent, the sender (Alice) interacts with her qubit using a not/flip/x gate and a z gate
	// then the qubit is sent to Bob, where it is decoded, ending with a measurement that gets back the classical bits

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SuperdenseCoding :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		SuperdenseCoding(int addseed = 0)
			: BaseClass(2, addseed), b1(false), b2(false)
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

			return BaseClass::Measure();
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
			bellState.setBellState00(BaseClass::reg);
			
			// the commented code is equivalent with the above, the two quantum gates making out the entangling gate (noted E or E2)
			//BaseClass::setToBasisState(0);
			//BaseClass::ApplyGate(hadamard, 0);
			//BaseClass::ApplyGate(cnot, 1, 0);
		}

		void Decode()
		{
			// this is the Bell decoding circuit, the inverse of the Bell encoding circuit
			// as the gates are unitary, the inverse is the conjugate transpose
			// so (hadamard * cnot)^dag = cnot^dag * hadamard^dag = cnot * hadamard

			BaseClass::ApplyGate(cnot, 1, 0);
			BaseClass::ApplyGate(hadamard, 0);

			// the results of this will be |00> for the 00 Bell state, |01> for the 01 Bell state, |10> for the 10 Bell state and |11> for the 11 Bell state
			// so measuring the qubits will give the classical bits back
			// see the `BellState` class for more details
		}

		// what this does is to transform the 00 Bell state into one of the other three Bell states (for the 01, 10 or 11 classical bits to be sent) or leave it unchanged for the 00 classical bits
		// the gate(s) need to be applied only on the Alice's qubit
		void Encode(bool bit1, bool bit2)
		{
			// the resemblance with the 'explicit' teleportation code is not a coincidence, teleportation and superdense coding are sort of opposite to each other
			if (bit2) BaseClass::ApplyGate(x, 0);
			if (bit1) BaseClass::ApplyGate(z, 0);
			// the above two gates (for encoding 11) are equivalent with -iY
		}

		QC::BellState<VectorClass, MatrixClass> bellState;

		// gates used at the end for decoding
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::CNOTGate<MatrixClass> cnot;

		// single qubit gates to act on the transmission qubit (depending on the classical bits to be transmitted)
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::PauliZGate<MatrixClass> z;

		// store the bits to send for 'Execute'
		bool b1;
		bool b2;
	};
}

