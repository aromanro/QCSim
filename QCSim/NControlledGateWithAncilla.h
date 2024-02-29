#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		// if this looks too complex, check out first NControlledNotWithAncilla
		// this might not be optimal, but it's a good starting point

		// allows only 1 or 2 qubit gates so the controlled gate used by this algorithm is 2 or 3 qubits

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NControlledGatesWithAncilla : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			NControlledGatesWithAncilla(size_t N = INT_MAX, size_t startQubit = 0, size_t endQubit = INT_MAX)
				: BaseClass(N, startQubit, endQubit), startAncillaQubits(1), clearAncillaAtTheEnd(true)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				if (controlQubits.empty())
					return 1;

				// 0. If there is only one control qubit, use directly the gate, no ancilla needed
				if (controlQubits.size() == 1)
				{
					for (auto& theGate : gates)
					{
						if (theGate.getQubitsNumber() == 2)
							theGate.setQubit2(controlQubits[0]);
						else
							theGate.setQubit3(controlQubits[0]);
						reg.ApplyGate(theGate);
					}

					return 0;
				}

				if (clearAncillaAtTheEnd) reg.ComputeStart();


				// 1. Walk over the control qubits, use ccnot on pairs of them with an ancilla qubit as target (different ancilla for each pair)
				size_t curFreeAncilla = startAncillaQubits;
				for (size_t i = 0; i < controlQubits.size() - 1; i += 2)
				{
					reg.ApplyGate(ccnot, curFreeAncilla, controlQubits[i], controlQubits[i + 1]);
					++curFreeAncilla;
				}

				// 2. One could remain unpaired, ccnot it with the first used ancilla qubit if that's the case, targeting the next unused ancilla

				size_t curAncilla = startAncillaQubits;
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

				for (auto& theGate : gates)
				{
					if (theGate.getQubitsNumber() == 2)
						theGate.setQubit2(curFreeAncilla);
					else
						theGate.setQubit3(curFreeAncilla);

					reg.ApplyGate(theGate);
				}

				if (clearAncillaAtTheEnd) reg.Uncompute();
	
				return 0; // no measurement, so the return should be ignored
			}

			const std::vector<size_t>& GetControlQubits() const
			{
				return controlQubits;
			}

			void SetControlQubits(const std::vector<size_t>& cq)
			{
				controlQubits = cq;
			}

			size_t GetStartAncillaQubits() const
			{
				return startAncillaQubits;
			}

			void SetStartAncillaQubits(size_t saq)
			{
				startAncillaQubits = saq;
			}

			const std::vector<Gates::AppliedGate<MatrixClass>>& GetGates() const
			{
				return gates;
			}

			bool AddGate(const Gates::AppliedGate<MatrixClass>& gate)
			{
				if (gate.getQubitsNumber() < 1 || gate.getQubitsNumber() > 2)
					return false;

				// make a controlled gate out of it
				MatrixClass m = gate.getRawOperatorMatrix();
				if (m.rows() != m.cols() || (m.rows() != 2 && m.rows() != 4))
					return false;

				Gates::AppliedGate<MatrixClass> theGate;
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

				gates.emplace_back(theGate);

				return true;
			}

			void ClearGates()
			{
				gates.clear();
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
			std::vector<size_t> controlQubits;
			size_t startAncillaQubits;
			std::vector<Gates::AppliedGate<MatrixClass>> gates;
			bool clearAncillaAtTheEnd;
		};

	}

}

