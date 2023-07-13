#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	// if this looks too complex, check out first NControlledNotWithAncilla
	// this might not be optimal, but it's a good starting point

	// allows only 1 or 2 qubit gates

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
			if (theGate.GetNrQubits() == 0 || theGate.GetNrQubits() > 2)
				return 1;

			if (clearAncillaAtTheEnd) reg.ComputeStart();

			// TODO: Implement it



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
			theGate = gate;
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
		unsigned int startAncillaQubits;

		Gates::AppliedGate<MatrixClass> theGate;
		bool clearAncillaAtTheEnd;
	};

}

