#pragma once
#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "CatEntangler.h"
#include "CatDisentangler.h"

namespace Teleportation
{

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumTeleportation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		QuantumTeleportation(unsigned int N = 3, int addseed = 0)
			: BaseClass(N, addseed)
		{
			BaseClass::setToBasisState(0);
		}

		void Entangle(unsigned int q1 = 1, unsigned int q2 = 2)
		{
			// make an EPR pair out of the second and third qubit:

			// starting from |00> (default) gets to (|00> + |11>)/sqrt(2)

			// the two gates make up an 'entangling gate', sometimes noted by E or E2
			BaseClass::ApplyGate(hadamard, q1);
			BaseClass::ApplyGate(cnot, q2, q1);
		}

		void ApplyTeleportationCircuit(unsigned int sentQubit = 0, unsigned int q2 = 1)
		{
			// the cnot and hadamard that follow do the inverse of the entangling gate - this way the measurement that follows is a measurement in the Bell basis

			// interacting the sending qubit (0) with the Alice's side of the entangled pair:
			BaseClass::ApplyGate(cnot, q2, sentQubit);

			// in the paper mentioned above, here follows the measurement of the second qubit
			// then based on it the x gate is applied (or not)
			// then the B mark would follow

			// from here it would follow what's after B
			BaseClass::ApplyGate(hadamard, sentQubit);
		}

		void RestoreState(bool firstQubitMeasurement, bool secondQubitMeasurement, unsigned int targetQubit = 2)
		{
			if (secondQubitMeasurement) BaseClass::ApplyGate(x, targetQubit);
			if (firstQubitMeasurement) BaseClass::ApplyGate(z, targetQubit);
		}

		unsigned int Teleport(unsigned int sentQubit = 0, unsigned int q2 = 1, unsigned int targetQubit = 2, bool explicitClassicalTransmission = false)
		{
			// TODO: Using something like the `BellState` class to create the EPR pair, not only the one currently used, but also the other three
			// it works with any of them as long as Alice and Bob agree on which one to use

			Entangle(q2, targetQubit);

			// here we are in the point A marked in fig 7 in the paper mentioned above

			ApplyTeleportationCircuit(sentQubit, q2);

			// measurements can be actually done all at the end, but let's pretend
			// one could do here only the first qubit measurement, the second one being done earlier

			// anyhow, it's basically the same thing
			unsigned int measuredValues = BaseClass::Measure(sentQubit);
			measuredValues |= BaseClass::Measure(q2) << 1;

			// now that they are measured they go to Bob, which uses the values to act on his entangled qubit:

			if (explicitClassicalTransmission)
			{
				const bool firstQubitMeasurement = (measuredValues & 0x1) != 0;
				const bool secondQubitMeasurement = (measuredValues & 0x2) != 0;

				RestoreState(firstQubitMeasurement, secondQubitMeasurement, targetQubit);
			}
			else
			{
				BaseClass::ApplyGate(cnot, targetQubit, q2);
				BaseClass::ApplyGate(cz, targetQubit, sentQubit);
			}

			return measuredValues;
		}

	protected:
		// the initial state to be teleported can be set directly or set up from something simple using some gate (hadamard, phase shift)
		// needed to create the Bell states
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::CNOTGate<MatrixClass> cnot; // also needed on Alice side

		// needed on Bob side, if the classical transmission is explicit (simple one qubit gates that are applied or not depending on the value of the classical bit):
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::PauliZGate<MatrixClass> z;

		// needed on receiving side (with no explicit measurement):
		QC::Gates::ControlledZGate<MatrixClass> cz;
	};

	// teleportation and superdense coding are sort of inverse of each other, see also the `SuperdenseCoding` class
	// both start with an EPR pair, teleportation sends two classical bits to recover a qubit to be sent at the other end, superdense coding sends a qubit to recover the two classical bits at the other end
	
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumTeleportationRealization : public QuantumTeleportation<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumTeleportation<VectorClass, MatrixClass>;
		using AlgorithmClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		QuantumTeleportationRealization(int addseed = 0)
			: BaseClass(3, addseed)
		{
		}

		void SetState(std::complex<double> alpha, std::complex<double> beta)
		{
			AlgorithmClass::Clear();
			AlgorithmClass::setRawAmplitude(0, alpha);
			AlgorithmClass::setRawAmplitude(1, beta);

			// ensure it's normalized:
			AlgorithmClass::Normalize();
		}


		// A more general discussion is in "GENERALIZED GHZ STATES AND DISTRIBUTED QUANTUM COMPUTING"
		// https://arxiv.org/abs/quant-ph/0402148
		// for this see fig 7

		// considers the initial state already set
		// the qubits to be entangled start as |00> and they are entangled to a EPR state using hadamard and cnot gates (together they make the E/E2 entangling gate)
		// the first two qubits are on the Alice side, the third belongs to Bob
		// returns the values of the two measured qubits
		unsigned int Teleport(bool explicitClassicalTransmission = false)
		{
			return BaseClass::Teleport(0, 1, 2, explicitClassicalTransmission);
		}

		// called only if one wants the teleported qubit to be measured, otherwise just check the register contents
		unsigned int Execute() override
		{
			Teleport();

			return AlgorithmClass::Measure(2, 2); // teleported qubit
		}
	};

}

namespace QC {

	namespace SubAlgo {

		// This implementation follows "GENERALIZED GHZ STATES AND DISTRIBUTED QUANTUM COMPUTING"
		// https://arxiv.org/abs/quant-ph/0402148
		// fig 7
		// using the 'cat entangler' and 'cat disentangler' sub-algorithms
		// since teleportation is quite important (very important for distributed quantum computing) it's worth having a sub-algorithm for it, to be used in other algorithms

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class Teleport : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			Teleport(unsigned int N = 3, unsigned int sourceQubit = 0, unsigned int entangledQubit = 1, unsigned int targetQubit = 2)
				: BaseClass(N, entangledQubit, targetQubit), sQubit(sourceQubit), entangler(N, entangledQubit, targetQubit), disentangler(N, targetQubit, sourceQubit, sourceQubit)
			{
			}

			unsigned int Execute(RegisterClass& reg) override
			{
				const unsigned int entQubit = BaseClass::getStartQubit();
				const unsigned int targetQubit = BaseClass::getEndQubit();

				entangler.Execute(reg);

				reg.ApplyGate(cnot, entQubit, sQubit);
				const unsigned int measurement = reg.Measure(entQubit);
				if (measurement) reg.ApplyGate(x, targetQubit);

				const unsigned int measurement2 = disentangler.Execute(reg);

				return measurement2 | (measurement << 1);
			}

		protected:
			unsigned int sQubit;
			Gates::CNOTGate<MatrixClass> cnot;
			Gates::PauliXGate<MatrixClass> x;
			CatEntangler<VectorClass, MatrixClass> entangler;
			CatDisentangler<VectorClass, MatrixClass> disentangler;
		};

	}

}
