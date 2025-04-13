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
// two implementations: one that explicitely constructs a big matrix for the function/oracle, the operator to be applied on the quantum register
// and one where the same is achieved by constructing a circuit out of simple quantum gates, to avoid going through that big matrix construction and multiplication
// the later case was made possible due of the optimization existing now for 1, 2 and 3 qubit gates (see ApplyGate implementation in register), otherwise going through the big matrix for each gate would be too slow

// since this is done after the mentioned optimization, I won't bother anymore with the big matrix construction, so I'll only implement the quantum circuit version
// the one with the big matrices would be quite slow, anyway


namespace QuantumCounting {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumCountingAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		QuantumCountingAlgorithm(size_t PrecisionQubits = 6, size_t GroverQubits = 4, unsigned int addseed = 0)
			: BaseClass(PrecisionQubits + 2 * GroverQubits, addseed), precisionQubits(PrecisionQubits), groverQubits(GroverQubits), fourier(PrecisionQubits, 0, PrecisionQubits - 1)
		{
			assert(PrecisionQubits >= 3);
			assert(GroverQubits >= 1);

			size_t firstAncillaQubit = PrecisionQubits + GroverQubits;
			nControlledNOT.SetTargetQubit(firstAncillaQubit);
			nControlledNOT.SetStartAncillaQubits(firstAncillaQubit + 1);
		}

		size_t Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure(0, precisionQubits - 1);
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(0, precisionQubits - 1, nrMeasurements);
		}

		void ClearMarkedStates()
		{
			markedStates.clear();
		}

		// start from zero, don't go over the maximum that is allowed by the GroverQubits
		void AddMarkedState(size_t state)
		{
			if (state >= (1ULL << groverQubits)) return;

			markedStates.insert(state);
		}

		void SetMarkedStates(const std::vector<size_t>& states)
		{
			ClearMarkedStates();
			for (size_t state : states)
				AddMarkedState(state);
		}

		size_t GetNumberOfMarkedStates() const
		{
			return static_cast<size_t>(markedStates.size());
		}

		double GetCorrectThetaForMarkedStates() const
		{
			return asin(static_cast<double>(sqrt(static_cast<double>(GetNumberOfMarkedStates())/ static_cast<double>(1ULL << groverQubits)))) * M_1_PI;
		}

		double GetThetaForState(size_t state) const
		{
			double theta = static_cast<double>(state) / static_cast<double>(1ULL << precisionQubits);

			if (theta >= 0.5)
				theta = 1.0 - theta;

			return theta;
		}

		size_t GetCountForState(size_t state) const
		{
			double s = sin(M_PI * GetThetaForState(state));

			return static_cast<size_t>(round(static_cast<double>(1ULL << groverQubits) * s * s));
		}

	private:
		void ExecuteWithoutMeasurement()
		{
			size_t firstAncillaQubit = precisionQubits + groverQubits;

			BaseClass::setToBasisState(0);

			for (size_t q = 0; q < firstAncillaQubit; ++q)
				BaseClass::ApplyGate(hadamard, q);

			BaseClass::ApplyGate(x, firstAncillaQubit);
			BaseClass::ApplyGate(hadamard, firstAncillaQubit);

			for (size_t q = 0; q < precisionQubits; ++q)
			{
				size_t nrOps = 1ULL << q;
				size_t ctrlQubit = precisionQubits - q - 1;
				// or this, if you want to use the IQFT with swapping, but obviously that's slower due of the additional swap gates
				//size_t ctrlQubit = q;
				for (size_t i = 0; i < nrOps; ++i)
					ControlledGrover(ctrlQubit);
			}
			
			BaseClass::ApplyGate(hadamard, firstAncillaQubit);
			BaseClass::ApplyGate(x, firstAncillaQubit);
			
			fourier.IQFT(BaseClass::reg, false);
		}

		void ControlledGrover(size_t ctrlQubit)
		{
			for (size_t state : markedStates)
				ControlledOracle(ctrlQubit, state);

			ControlledDiffusion(ctrlQubit);
		}


		void ControlledOracle(size_t ctrlQubit, size_t state)
		{
			size_t v = state;
			for (size_t q = 0; q < groverQubits; ++q)
			{
				if ((v & 1) == 0)
					BaseClass::ApplyGate(cx, precisionQubits + q, ctrlQubit);

				v >>= 1;
			}
			
			
			std::vector<size_t> controlQubits(groverQubits + 1);
			std::iota(controlQubits.begin(), controlQubits.end(), precisionQubits);
			controlQubits.back() = ctrlQubit;
			
			nControlledNOT.SetControlQubits(controlQubits);
			nControlledNOT.Execute(BaseClass::reg);
			
			v = state;
			for (size_t q = 0; q < groverQubits; ++q)
			{
				if ((v & 1) == 0)
					BaseClass::ApplyGate(cx, precisionQubits + q, ctrlQubit);

				v >>= 1;
			}
		}


		void ControlledDiffusion(size_t ctrlQubit)
		{
			for (size_t q = 0; q < groverQubits; ++q)
				BaseClass::ApplyGate(chadamard, precisionQubits + q, ctrlQubit);

			size_t nrGroverStates = 1ULL << groverQubits;
			for (size_t state = 1; state < nrGroverStates; ++state)
				ControlledOracle(ctrlQubit, state);

			for (size_t q = 0; q < groverQubits; ++q)
				BaseClass::ApplyGate(chadamard, precisionQubits + q, ctrlQubit);
		}


		std::set<size_t> markedStates;

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::ControlledHadamardGate<MatrixClass> chadamard;

		QC::Gates::PauliXGate<MatrixClass> x;
		QC::Gates::CNOTGate<MatrixClass> cx;

		QC::SubAlgo::NControlledNotWithAncilla<VectorClass, MatrixClass> nControlledNOT;
		
		size_t precisionQubits;
		size_t groverQubits;

		QC::SubAlgo::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
	};

}