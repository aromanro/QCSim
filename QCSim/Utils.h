#pragma once

#include "QubitRegister.h"

namespace QC {

	// TODO: implement & use it
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BellState
	{
	public:
		bool setBellState(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1, bool s1 = false, bool s2 = false)
		{
			if (qubit1 == qubit2 || qubit1 >= reg.getNrQubits() || qubit2 >= reg.getNrBasisStates())
				return false;

			const unsigned int state1 = 1u << qubit1;
			const unsigned int state2 = 1u << qubit2;

			reg.setToBasisState((s1 ? state1 : 0) | (s2 ? state2 : 0));

			reg.ApplyGate(hadamard, qubit1);
			reg.ApplyGate(cnot, qubit2, qubit1);

			return true;
		}

		bool setBellState00(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1)
		{
			return setBellState(reg, qubit1, qubit2, false, false);
		}

		bool setBellState01(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1)
		{
			return setBellState(reg, qubit1, qubit2, false, true);
		}

		bool setBellState10(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1)
		{
			return setBellState(reg, qubit1, qubit2, true, false);
		}

		bool setBellState11(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1)
		{
			return setBellState(reg, qubit1, qubit2, true, true);
		}

	protected:
		QC::HadamardGate<MatrixClass> hadamard;
		QC::CNOTGate<MatrixClass> cnot;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class MeasurementBasis
	{
	public:
		bool switchToXBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit = 0)
		{
			if (qubit >= reg.getNrQubits())
				return false;

			reg.ApplyGate(hadamard, qubit);

			return true;
		};

		bool switchToYBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit = 0)
		{
			if (qubit >= reg.getNrQubits())
				return false;

			static const SingleQubitGate<MatrixClass> U(hadamard.getOperatorMatrix() * s.getOperatorMatrix().adjoint());

			reg.ApplyGate(U, qubit);

			return true;
		};

		bool switchToBellBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1)
		{
			if (qubit1 == qubit2 || qubit1 >= reg.getNrQubits() || qubit2 >= reg.getNrBasisStates())
				return false;

			reg.ApplyGate(cnot, qubit2, qubit1);
			reg.ApplyGate(hadamard, qubit1);

			return true;
		};

	protected:
		QC::HadamardGate<MatrixClass> hadamard;
		QC::CNOTGate<MatrixClass> cnot;
		QC::PhaseGate<MatrixClass> s;
	};

}

