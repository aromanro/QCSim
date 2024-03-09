#pragma once

#include "QuantumAlgorithm.h"
#include "Oracle.h"

#include "NControlledNotWithAncilla.h"
#include "HadamardTransform.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <numeric>

namespace Grover {

	template<class MatrixClass = Eigen::MatrixXcd> class Oracle :
		public QC::Gates::QuantumGate<MatrixClass>
	{
	public:
		void setCorrectQuestionState(size_t state)
		{
			correctQuestionState = state;
		}

		MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
		{
			const size_t nrBasisStates = 1ULL << nrQubits;
			MatrixClass extOperatorMat = MatrixClass::Identity(nrBasisStates, nrBasisStates);

			if (correctQuestionState < nrBasisStates)
				extOperatorMat(correctQuestionState, correctQuestionState) = -1;

			assert(checkUnitary(extOperatorMat));

			return extOperatorMat;
		}

	protected:
		size_t correctQuestionState = 0;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class GroverAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		GroverAlgorithm(size_t N = 3, int addseed = 0)
			: BaseClass(N, addseed)
		{
			assert(N >= 1);

			Oracle<MatrixClass> o00;
			o00.setCorrectQuestionState(0); // already set to 0 by default, but no harm in setting it explicitly, to show the intention
			OracleOp00 = o00.getOperatorMatrix(BaseClass::getNrQubits());

			setCorrectQuestionState(0);
		}

		void setCorrectQuestionState(size_t state)
		{
			Oracle<MatrixClass> o;
			o.setCorrectQuestionState(state);
			OracleOp = o.getOperatorMatrix(BaseClass::getNrQubits());
		}

		size_t Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure();
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(nrMeasurements);
		}

	protected:
		void Init()
		{
			//BaseClass::reg.setToBasisState(0);
			//ApplyHadamardOnAllQubits();
			BaseClass::setToEqualSuperposition(); // the same thing as commented above
		}

		void ExecuteWithoutMeasurement()
		{
			Init();

			const size_t repeatNo = static_cast<size_t>(round(M_PI / 4. * sqrt(BaseClass::getNrBasisStates())));

			for (size_t i = 0; i < repeatNo; ++i)
			{
				BaseClass::ApplyOperatorMatrix(OracleOp);
				ApplyDiffusionOperator();
			}
		}

		void ApplyHadamardOnAllQubits()
		{
			for (size_t i = 0; i < BaseClass::getNrQubits(); ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyDiffusionOperator()
		{
			ApplyHadamardOnAllQubits();
			BaseClass::ApplyOperatorMatrix(OracleOp00);
			ApplyHadamardOnAllQubits();
		}

		MatrixClass OracleOp;
		MatrixClass OracleOp00;

		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};


	// as above, but using gates instead of the big matrix for the oracle
	// needs ancilla qubits
	// There is one ancilla qubit used by the algorithm
	// the n controlled not needs N - 2 ancilla qubits 
	// that's how many additional ancilla qubits are needed to reduce the control qubits down to 2 for a ccnot (or 1 for a cnot)

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class GroverAlgorithmWithGatesOracle :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		GroverAlgorithmWithGatesOracle(size_t N = 3, int addseed = 0)
			: BaseClass(2 * N - 1, addseed), nControlledNOT(INT_MAX), correctQuestionState(0)
		{
			assert(N >= 1);

			std::vector<size_t> controlQubits(N);
			std::iota(controlQubits.begin(), controlQubits.end(), 0);

			nControlledNOT.SetControlQubits(controlQubits);
			nControlledNOT.SetTargetQubit(N);
			nControlledNOT.SetStartAncillaQubits(N + 1);
		}

		void setCorrectQuestionState(size_t state)
		{
			correctQuestionState = state;
		}

		size_t Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure(0, getAlgoQubits() - 1);
		}

		std::map<size_t, size_t> ExecuteWithMultipleMeasurements(size_t nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(0, getAlgoQubits() - 1, nrMeasurements);
		}

		size_t getAlgoQubits() const
		{
			const size_t nrQubits = BaseClass::getNrQubits();

			return (nrQubits + 1) / 2;
		}

	protected:
		void Init()
		{
			BaseClass::reg.setToBasisState(0);
			ApplyHadamardOnAllQubits();
			BaseClass::ApplyGate(x, getAlgoQubits());
		}

		void ExecuteWithoutMeasurement()
		{
			Init();

			const size_t nrQubits = getAlgoQubits();
			const size_t nrBasisStates = 1ULL << nrQubits;
			const size_t repeatNo = static_cast<size_t>(round(M_PI / 4. * sqrt(nrBasisStates)));

			for (size_t i = 0; i < repeatNo; ++i)
			{
				ApplyOracle(correctQuestionState);
				ApplyDiffusionOperator();
			}
		}

		void ApplyHadamardOnAllQubits()
		{
			for (size_t i = 0; i < getAlgoQubits(); ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyDiffusionOperator()
		{
			ApplyHadamardOnAllQubits();
			ApplyOracle(0);
			ApplyHadamardOnAllQubits();
		}

		void ApplyOracle(size_t state = 0)
		{
			const size_t nrQubits = getAlgoQubits();

			BaseClass::ApplyGate(hadamard, nrQubits);
			
			size_t v = state;
			for (size_t q = 0; q < nrQubits; ++q)
			{
				if ((v & 1) == 0)
					BaseClass::ApplyGate(x, q);

				v >>= 1;
			}

			nControlledNOT.Execute(BaseClass::reg);
	
			v = state;
			for (size_t q = 0; q < nrQubits; ++q)
			{
				if ((v & 1) == 0)
					BaseClass::ApplyGate(x, q);

				v >>= 1;
			}

			BaseClass::ApplyGate(hadamard, nrQubits);
		}

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::SubAlgo::NControlledNotWithAncilla<VectorClass, MatrixClass> nControlledNOT;
		QC::Gates::PauliXGate<MatrixClass> x;
		size_t correctQuestionState;
	};

}
