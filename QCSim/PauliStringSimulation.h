#pragma once

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QuantumSimulation {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliZStringSimulation :
		public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		PauliZStringSimulation(double t = 0)
		{
			setTimeInterval(t);
		}

		unsigned int Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			const unsigned int nrQubits = reg.getNrQubits();
			const unsigned int lastQubit = nrQubits - 1;

			for (unsigned int qubit = 0; qubit < lastQubit; ++qubit)
				reg.ApplyGate(cnot, qubit + 1, qubit);

			reg.ApplyGate(rz, lastQubit);

			for (unsigned int qubit = lastQubit; qubit > 0; --qubit)
				reg.ApplyGate(cnot, qubit, qubit - 1);

			return 0; // not used here
		}

		void setTimeInterval(double t)
		{
			rz.SetTheta(2. * t);
		}

	protected:
		QC::CNOTGate<MatrixClass> cnot;
		QC::RzGate<MatrixClass> rz;
	};

	// you may decompose a Hamiltonian in a sum of such terms with real coeffcients

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliStringSimulation :
		public PauliZStringSimulation<VectorClass, MatrixClass>
	{
	public:
		PauliStringSimulation(double t = 0, unsigned int N = 3)
			: PauliZStringSimulation<VectorClass, MatrixClass>(t), ops(N, PauliOp::opZ)
		{
		}

		enum class PauliOp : unsigned char
		{
			opZ = 0,
			opX = 1,
			opY = 2
		};

		void setOperatorForQubit(unsigned int qubit, PauliOp op)
		{
			if (qubit >= ops.size()) return;

			ops[qubit] = op;
		}

		PauliOp getOperatorForQubit(unsigned int qubit) const
		{
			if (qubit >= ops.size()) return PauliOp::opZ;

			return ops[qubit];
		}

		unsigned int Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) const override
		{
			const unsigned int nrQubits = reg.getNrQubits();

			static const QC::PauliZGate<MatrixClass> Z;

			// switch to proper basis
			for (unsigned int qubit = 0; qubit < nrQubits; ++qubit)
				if (ops[qubit] == PauliOp::opX)
					measurementBasis.switchToXBasis(reg, qubit);
				else if (ops[qubit] == PauliOp::opY)
					measurementBasis.switchToYBasis(reg, qubit);
				// for Z we don't have to do anything, that one is the computational basis

			PauliZStringSimulation<VectorClass, MatrixClass>::Execute(reg);

			// switch back to computational basis
			for (unsigned int qubit = 0; qubit < nrQubits; ++qubit)
				if (ops[qubit] == PauliOp::opX)
					measurementBasis.switchToComputationalFromXBasis(reg, qubit);
				else if (ops[qubit] == PauliOp::opY)
					measurementBasis.switchToComputationalFromYBasis(reg, qubit);
				// for Z we don't have to do anything, that one is the computational basis

			return 0; // not used here, measurements need to be done explicitely
		}

	protected:
		QC::MeasurementBasis<VectorClass, MatrixClass> measurementBasis;
		std::vector<PauliOp> ops;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliDecomposedHamiltonianSimulation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		PauliDecomposedHamiltonianSimulation(unsigned int N = 3, double t = 0, unsigned int nrSteps = 1, bool addSeed = false)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addSeed), simTime(t), steps(nrSteps)
		{
		}

		void AddTerm(double coeff, const PauliStringSimulation<VectorClass, MatrixClass>& term)
		{
			coeffs.push_back(coeff);
			terms.push_back(term);
		}

		void Clear()
		{
			coeffs.clear();
			terms.clear();
		}

		unsigned int Execute() override
		{
			const double deltat = simTime / steps;

			for (unsigned int i = 0; i < terms.size(); ++i)
				terms[i].setTimeInterval(deltat * coeffs[i]);

			// Suzuki-Trotter expansion
			for (unsigned int step = 0; step < steps; ++step)
			{
				// this is slow and it can be improved by obtaining the whole matrix once and just by 'applying' the operator each step
				// for now I'll leave it as it is but the gains would be quite big here
				for (const PauliStringSimulation<VectorClass, MatrixClass>& term : terms)
					term.Execute(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);
			}

			return 0; // don't measure for this one, must be done explicitely, needed to be able to check the register first
		}

		void SetSimulationTime(double t)
		{
			simTime = t;
		}

		void SetNrSteps(unsigned int nrSteps)
		{
			steps = nrSteps;
		}

		// this is here only for testing purposes
		MatrixClass getEvolutionOperator(double t) const
		{
			static const QC::PauliXGate<MatrixClass> X;
			static const QC::PauliYGate<MatrixClass> Y;
			static const QC::PauliZGate<MatrixClass> Z;

			MatrixClass H = MatrixClass::Zero(QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates(), QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates());

			for (unsigned int i = 0; i < terms.size(); ++i)
			{
				MatrixClass tm;
				for (unsigned int qubit = 0; qubit < QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits(); ++qubit)
				{
					if (terms[i].getOperatorForQubit(qubit) == PauliStringSimulation<>::PauliOp::opX)
					{
						if (qubit == 0) tm = X.getRawOperatorMatrix();
						else tm = Eigen::kroneckerProduct(X.getRawOperatorMatrix(), tm).eval();
					}
					else if (terms[i].getOperatorForQubit(qubit) == PauliStringSimulation<>::PauliOp::opY)
					{
						if (qubit == 0) tm = Y.getRawOperatorMatrix();
						else tm = Eigen::kroneckerProduct(Y.getRawOperatorMatrix(), tm).eval();
					}
					else
					{
						if (qubit == 0) tm = Z.getRawOperatorMatrix();
						else tm = Eigen::kroneckerProduct(Z.getRawOperatorMatrix(), tm).eval();
					}
				}

				H += coeffs[i] * tm;
			}

			return (std::complex<double>(0, -1) * t * H).exp();
		}

	protected:
		double simTime;
		unsigned int steps;

		std::vector<double> coeffs;
		std::vector<PauliStringSimulation<VectorClass, MatrixClass>> terms;
	};
}
