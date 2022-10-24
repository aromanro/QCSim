#pragma once

#include <algorithm>

#include "QubitRegister.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BellState
	{
	public:
		bool setBellState(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1, bool s1 = false, bool s2 = false)
		{
			if (qubit1 == qubit2 || qubit1 >= reg.getNrQubits() || qubit2 >= reg.getNrQubits())
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
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::CNOTGate<MatrixClass> cnot;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class MeasurementBasis
	{
	public:
		MeasurementBasis()
		{
			//s.SetPhaseShift(M_PI_2); // no need for this, PhaseShift has already the proper phase
		}

		bool switchToXBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit = 0) const
		{
			if (qubit >= reg.getNrQubits())
				return false;

			reg.ApplyGate(hadamard, qubit);

			return true;
		};

		bool switchToYBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit = 0) const
		{
			if (qubit >= reg.getNrQubits())
				return false;

			//static const SingleQubitGate<MatrixClass> U(hadamard.getRawOperatorMatrix() * s.getRawOperatorMatrix().adjoint());

			reg.ApplyGate(hy, qubit);

			return true;
		};

		bool switchToBellBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1) const
		{
			if (qubit1 == qubit2 || qubit1 >= reg.getNrQubits() || qubit2 >= reg.getNrQubits())
				return false;

			reg.ApplyGate(cnot, qubit2, qubit1);
			reg.ApplyGate(hadamard, qubit1);

			return true;
		};

		// 'computational' is Z
		// if the operator to switch to a basis was U, the one to switch back is U^t
		// if it's a product of two (as in two gates used for changing the basis) then (O1 * O2)^t = O2^t * O1^t

		bool switchToComputationalFromXBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit = 0) const
		{
			if (qubit >= reg.getNrQubits())
				return false;

			reg.ApplyGate(hadamard, qubit);

			return true;
		}

		bool switchToComputationalFromYBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit = 0) const
		{
			if (qubit >= reg.getNrQubits())
				return false;

			//static const QC::SingleQubitGate<MatrixClass> U(s.getRawOperatorMatrix() * hadamard.getRawOperatorMatrix());
			
			reg.ApplyGate(hy, qubit);

			return true;
		}

		bool switchToComputationalFromBellBasis(QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit1 = 0, unsigned int qubit2 = 1) const
		{
			if (qubit1 == qubit2 || qubit1 >= reg.getNrQubits() || qubit2 >= reg.getNrQubits())
				return false;

			reg.ApplyGate(hadamard, qubit1);
			reg.ApplyGate(cnot, qubit2, qubit1);

			return true;
		};

		// the reason of the 'switchBack' flag is for various applications:
		// for example for the cases when one wants to start with some other basis for register initialization then switch to the computational (Z) basis
		// or for the case when one wants to switch to some other basis along the algoritm, measure, then switch back to the computational basis, then continue
		// an application could be in quantum cryptography - see for example BB84 protocol

		bool switchToOperatorBasis(QubitRegister<VectorClass, MatrixClass>& reg, const QC::Gates::SingleQubitGate<MatrixClass>& gate, unsigned int qubit = 0, bool switchBack = false) const
		{
			return switchToOperatorBasis(reg, gate.getRawOperatorMatrix(), qubit, switchBack);
		}

		bool switchToOperatorBasis(QubitRegister<VectorClass, MatrixClass>& reg, const MatrixClass& op, unsigned int qubit = 0, bool switchBack = false) const
		{
			if (qubit >= reg.getNrQubits()) return false;
			else if (op.rows() != 2 || op.cols() != 2) return false;

			const std::complex<double> trace = op(0, 0) + op(1, 1);
			const std::complex<double> det = op(0, 0) * op(1, 1) - op(0, 1) * op(1, 0);
			const std::complex<double> halftrace = 0.5 * trace;
			const std::complex<double> sdelta = std::sqrt(halftrace * halftrace - det);

			const std::complex<double> ev1 = halftrace + sdelta;
			const std::complex<double> ev2 = halftrace - sdelta;

			MatrixClass U = MatrixClass::Identity(2, 2);

			if (abs(op(1, 0)) > 1E-15)
			{
				const std::complex<double> dc = std::conj(op(1, 0));
				U(0, 0) = std::conj(ev1 - op(1, 1));
				U(0, 1) = dc;

				double norm = 1. / sqrt((U(0, 0) * std::conj(U(0, 0)) + U(0, 1) * std::conj(U(0, 1))).real());
				U(0, 0) *= norm;
				U(0, 1) *= norm;

				U(1, 0) = std::conj(ev2 - op(1, 1));
				U(1, 1) = dc;

				norm = 1. / sqrt((U(1, 0) * std::conj(U(1, 0)) + U(1, 1) * std::conj(U(1, 1))).real());
				U(1, 0) *= norm;
				U(1, 1) *= norm;
			}
			else if (abs(op(0, 1)) > 1E-15)
			{
				const std::complex<double> bc = std::conj(op(0, 1));
				U(0, 0) = bc;
				U(0, 1) = std::conj(ev1 - op(0, 0));

				double norm = 1. / sqrt((U(0, 0) * std::conj(U(0, 0)) + U(0, 1) * std::conj(U(0, 1))).real());
				U(0, 0) *= norm;
				U(0, 1) *= norm;

				U(1, 0) = bc;
				U(1, 1) = std::conj(ev2 - op(0, 0));

				norm = 1. / sqrt((U(1, 0) * std::conj(U(1, 0)) + U(1, 1) * std::conj(U(1, 1))).real());
				U(1, 0) *= norm;
				U(1, 1) *= norm;
			}
			else
			{
				// else do almost nothing, the operator is diagonal, probably Z (might also be -Z, in which case swap the values)
				if (op(0, 0).real() < op(1, 1).real())
					U.col(0).swap(U.col(1));
			}

			if (switchBack)
				U.adjointInPlace();

			const QC::Gates::SingleQubitGate<MatrixClass> gate(U);
			reg.ApplyGate(gate, qubit);

			return true;
		}

	protected:
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::HyGate<MatrixClass> hy;
		QC::Gates::CNOTGate<MatrixClass> cnot;
		//QC::Gates::PhaseGate<MatrixClass> s;
	};

}

