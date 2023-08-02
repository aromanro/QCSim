#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumFourierTransform.h"

#include "NControlledNotWithAncilla.h"
#include "HadamardTransform.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <set>
#include <numeric>

// the other more complex algorithms that are have functions/oracles have (usually, at this moment Shor doesn't, but I'll probably change that)
// two implementations: one that explicitely contructs a big matrix for the function/oracle, the operator to be applied on the quantum register
// and one where the same is achieved by constructing a circuit out of simple quantum gates, to avoid going through that big matrix construction and multiplication
// the later case was made possible due of the optimization existing now for 1, 2 and 3 qubit gates, otherwise going through the big matrix for each gate would be too slow

// since this is done after the mentioned optimization, I won't bother anymore with the big matrix construction, so I'll only implement the quantum circuit version


namespace QuantumCounting {


	// WARNING: this is just a sketch, I'll have to review it and test it!!!!!!!!!!!!!!
	// I implemented it quite fast and without proper care, so it's probably full of mistakes

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumCountingAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		QuantumCountingAlgorithm(unsigned int PrecisionQubits = 4, unsigned int GroverQubits = 4, int addseed = 0)
			: BaseClass(PrecisionQubits + 2 * GroverQubits, addseed), precisionQubits(PrecisionQubits), groverQubits(GroverQubits), fourier(PrecisionQubits + 2 * GroverQubits, 0, PrecisionQubits)
		{
			assert(PrecisionQubits >= 3);
			assert(GroverQubits >= 1);

			unsigned int firstAncillaQubit = PrecisionQubits + GroverQubits;
			nControlledNOT.SetTargetQubit(firstAncillaQubit);
			nControlledNOT.SetStartAncillaQubits(firstAncillaQubit + 1);
		}

		unsigned int Execute() override
		{
			unsigned int firstAncillaQubit = PrecisionQubits + GroverQubits;

			for (unsigned int q = 0; q < precisionQubits + groverQubits; ++q)
				BaseClass::Apply(hadamard, q);

			BaseClass::Apply(x, firstAncillaQubit);
			BaseClass::Apply(hadamard, firstAncillaQubit);


			for (unsigned int q = 0; q < precisionQubits; ++q)
			{
				unsigned int nrOps = 1 << q;
				for (unsigned int i = 0; i < nrOps; ++i)
					ControlledGrover(q);
			}


			BaseClass::Apply(hadamard, firstAncillaQubit);
			BaseClass::Apply(x, firstAncillaQubit);

			fourier.IQFT(BaseClass::reg);

			return BaseClass::Measure(0, precisionQubits - 1);
		}

		void ClearMarkedStates()
		{
			markedStates.clear();
		}

		// start from zero, don't go over the maximum that is allowed by the GroverQubits
		void AddMarkedState(unsigned int state)
		{
			if (state >= (1 << groverQubits)) return;

			markedStates.insert(state);
		}

	protected:
		void ControlledGrover(unsigned int ctrlQubit)
		{
			for (unsigned int state : markedStates)
				ControlledOracle(crtlQubit, state);

			ControlledDiffusion(ctrlQubit);
		}


		void ControlledOracle(unsigned int ctrlQubit, unsigned int state)
		{
			unsigned int v = state;
			for (unsigned int q = 0; q < groverQubits; ++q)
			{
				if ((v & 1) == 0)
					BaseClass::ApplyGate(cx, precisionQubits + q, ctrlQubit);

				v >>= 1;
			}

			std::vector<unsigned int> controlQubits(groverQubits + 1);
			std::iota(controlQubits.begin(), controlQubits.end() - 1, precisionQubits);
			constrolQubits.back() = ctrlQubit;

			nControlledNOT.SetControlQubits(controlQubits);

			nControlledNOT.Execute(BaseClass::reg);

			unsigned int v = state;
			for (unsigned int q = 0; q < groverQubits; ++q)
			{
				if ((v & 1) == 0)
					BaseClass::ApplyGate(cx, precisionQubits + q, ctrlQubit);

				v >>= 1;
			}
		}


		void ControlledDiffusion(unsigned int ctrlQubit)
		{
			for (unsigned int q = 0; q < qroverQubits; ++q)
				BaseClass::Apply(chadamard, precisionQubits + q, ctrlQubit);

			unsigned int nrGroverStates = 1 << groverQubits;
			for (unsigned int state = 0; state < nrGroverStates; ++state)
				ControlledOracle(ctrlQubit, state);

			for (unsigned int q = 0; q < qroverQubits; ++q)
				BaseClass::Apply(chadamard, precisionQubits + q, ctrlQubit);
		}


		std::set<unsigned int> markedStates;

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::ControlledHadamardGate<MatrixClass> chadamard;

		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::CNOTGate<MatrixClass> cx;

		QC::SubAlgo::NControlledNotWithAncilla<VectorClass, MatrixClass> nControlledNOT;
		QuantumFourierTransform<VectorClass, MatrixClass> fourier;

		unsigned int precisionQubits;
		unsigned int groverQubits;
	};

}