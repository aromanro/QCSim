#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	// this could be simply implemented with NControlledGateWithAncilla, but it's simpler and for helping understanging that one, it's easier to look here first
	// this might not be optimal, but it's a good starting point

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NControlledNotWithAncilla : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		NControlledNotWithAncilla(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX)
			: BaseClass(N, startQubit, endQubit), nrQubits(N), targetQubit(0), startAncillaQubits(1), clearAncillaAtTheEnd(true)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			// TODO: Implement it
			// TODO: Check if there are enough ancilla qubits

			// below is just a sketch of the algorithm, not tested, I have to review it, probably it has issues

			if (controlQubits.empty()) return 1;

			// 0. If the control qubits are <= 2, use directly the ccnot or cnot gate
			if (controlQubits.size() == 2)
			{
				reg.ApplyGate(ccnot, targetQubit, controlQubits[0], controlQubits[1]);
				return 0;
			}
			else if (controlQubits.size() == 1)
			{
				reg.ApplyGate(cnot, targetQubit, controlQubits[0]);
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

			unsigned int curAncillaToCCNOT = startAncillaQubits;
			if (controlQubits.size() % 2)
			{
				reg.ApplyGate(ccnot, curFreeAncilla, controlQubits[controlQubits.size() - 1], curAncillaToCCNOT);
				++curFreeAncilla;
				++curAncillaToCCNOT;
			}

			// 3. Start pairing ancilla qubits with ccnot, targeting the next free ancilla qubit, until either two or one are left
			while (curFreeAncilla - curAncillaToCCNOT > 2)
			{
				reg.ApplyGate(ccnot, curFreeAncilla, curAncillaToCCNOT, curAncillaToCCNOT + 1);
				++curFreeAncilla;
				curAncillaToCCNOT += 2;
			}

			if (clearAncillaAtTheEnd) reg.ComputeEnd();

			// 4. If two are left, ccnot them, otherwise use cnot, targeting the target qubit
			if (curFreeAncilla - curAncillaToCCNOT == 2)
				reg.ApplyGate(ccnot, targetQubit, curAncillaToCCNOT, curAncillaToCCNOT + 1);
			else
				reg.ApplyGate(cnot, targetQubit, curAncillaToCCNOT);

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

		unsigned int GetTargetQubit() const
		{
			return targetQubit;
		}

		void SetTargetQubit(unsigned int tq)
		{
			targetQubit = tq;
		}

		unsigned int GetStartAncillaQubits() const
		{
			return startAncillaQubits;
		}

		void SetStartAncillaQubits(unsigned int saq)
		{
			startAncillaQubits = saq;
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
		Gates::CNOTGate<MatrixClass> cnot;
		Gates::ToffoliGate<MatrixClass> ccnot;
		unsigned int nrQubits;
		std::vector<unsigned int> controlQubits;
		unsigned int targetQubit;
		unsigned int startAncillaQubits;
		bool clearAncillaAtTheEnd;
	};

}


