#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	// if this looks too complex, check out first NControlledNotWithAncilla
	// this might not be optimal, but it's a good starting point

	// allows only 1 or 2 qubit gates so the controlled gate used by this algorithm is 2 or 3 qubits

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NControlledGateWithAncilla : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		NControlledNotWithAncilla(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: BaseClass(N, startQubit, endQubit), nrQubits(N), startAncillaQubits(1), clearAncillaAtTheEnd(true)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			if (theGate.GetNrQubits() <= 1 || theGate.GetNrQubits() > 3 || controlQubits.empty())
				return 1;

			// 0. If there is only one control qubit, use directly the gate, no ancilla needed
			if (controlQubits.size() == 1)
			{
				if (theGate.GetNrQubits() == 2)
					theGate.setQubit2(controlQubits[0]);
				else
					theGate.setQubit3(controlQubits[0]);
				reg.ApplyGate(theGate);

				return 0;
			}

			if (clearAncillaAtTheEnd) reg.ComputeStart();


			// 1. Walk over the control qubits, use ccnot on pairs of them with an ancilla qubit as target (different ancilla for each pair)
			unsigned int curFreeAncilla = startAncillaQubits;
			for (int i = 0; i < controlQubits.size() - 1; i += 2)
			{
				reg.ApplyGate(ccnot, curFreeAncilla, controlQubits[i], controlQubits[i + 1]);
				++curFreeAncilla;
			}

			// 2. One could remain unpaired, ccnot it with the first used ancilla qubit if that's the case, targeting the next unused ancilla

			unsigned int curAncilla = startAncillaQubits;
			if (controlQubits.size() % 2)
			{
				reg.ApplyGate(ccnot, curFreeAncilla, controlQubits[controlQubits.size() - 1], curAncilla);
				++curFreeAncilla;
				++curAncilla;
			}

			// 3. Start pairing ancilla qubits with ccnot, targeting the next free ancilla qubit, until only one is left
			while (curFreeAncilla - curAncilla > 1)
			{
				reg.ApplyGate(ccnot, curFreeAncilla, curAncilla, curAncilla + 1);
				++curFreeAncilla;
				curAncilla += 2;
			}

			if (clearAncillaAtTheEnd) reg.ComputeEnd();

			// 4. Apply the gate on the last ancilla qubit and the target qubit
			--curFreeAncilla;
			if (theGate.GetNrQubits() == 2)
				theGate.setQubit2(curFreeAncilla);
			else
				theGate.setQubit3(curFreeAncilla);
			reg.ApplyGate(theGate);

			if (clearAncillaAtTheEnd) reg.Uncompute();

			return 0; // no measurement, so the return should be ignored
		}

		const std::vector<unsigned int>& GetControlQubits() const
		{
			return controlQubits;
		}

		void SetControlQubits(const std::vector<unsigned int>& cq)
		{
			controlQubits = cq;
		}

		unsigned int GetStartAncillaQubits() const
		{
			return startAncillaQubits;
		}

		void SetStartAncillaQubits(unsigned int saq)
		{
			startAncillaQubits = saq;
		}

		const Gates::AppliedGate<MatrixClass>& GetGate() const
		{
			return theGate;
		}

		void SetGate(const Gates::AppliedGate<MatrixClass>& gate)
		{
			// make a controlled gate out of it
			MatrixClass m = gate.GetRawOperatorMatrix();
			if (m.rows() != m.cols() || (m.rows() != 2 && m.rows() != 4))
				return;

			if (m.rows() == 2)
			{
				MatrixClass opMat = MatrixClass::Identity(4, 4);
				opMat.block(2, 2, 2, 2) = m;
				theGate.setOperator(opMat);	
			}
			else if (m.rows() == 4)
			{
				MatrixClass opMat = MatrixClass::Identity(8, 8);
				opMat.block(4, 4, 4, 4) = m;
				theGate.setOperator(opMat);
				theGate.setQubit2(gate.getQubit2());
			}
			theGate.setQubit1(gate.getQubit1());
		}

		bool ClearAncillaAtTheEnd() const
		{
			return clearAncillaAtTheEnd;
		}

		void SetClearAncillaAtTheEnd(bool caate)
		{
			clearAncillaAtTheEnd = caate;
		}

	protected:
		Gates::ToffoliGate<MatrixClass> ccnot;
		unsigned int nrQubits;
		std::vector<unsigned int> controlQubits;
		unsigned int startAncillaQubits;

		Gates::AppliedGate<MatrixClass> theGate;
		bool clearAncillaAtTheEnd;
	};

}

