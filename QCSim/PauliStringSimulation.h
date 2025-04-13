#pragma once

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"
#include "PauliString.h"

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

		size_t Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			const size_t nrQubits = reg.getNrQubits();
			const size_t lastQubit = nrQubits - 1;

			for (size_t qubit = 0; qubit < lastQubit; ++qubit)
				reg.ApplyGate(cnot, qubit + 1, qubit);

			reg.ApplyGate(rz, lastQubit);

			for (size_t qubit = lastQubit; qubit > 0; --qubit)
				reg.ApplyGate(cnot, qubit, qubit - 1);

			return 0; // not used here
		}

		void setTimeInterval(double t)
		{
			rz.SetTheta(2. * t);
		}

	protected:
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::RzGate<MatrixClass> rz;
	};

	// you may decompose a Hamiltonian in a sum of such terms with real coefficients
	// for details on fermionic Hamiltonians see Jordan-Wigner transformation in 3.1 in "Simulation of Electronic Structure Hamiltonians Using Quantum Computers" by James D. Whitfield, Jacob Biamonte, and Alan Aspuru-Guzik
	// https://arxiv.org/abs/1001.3855 
	// or this page: https://learn.microsoft.com/en-us/azure/quantum/user-guide/libraries/chemistry/concepts/jordan-wigner

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliStringSimulation :
		public PauliZStringSimulation<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = PauliZStringSimulation<VectorClass, MatrixClass>;

		PauliStringSimulation(double t = 0, size_t N = 3)
			: BaseClass(t), pauliString(N)
		{
		}

		void setOperatorForQubit(size_t qubit, PauliString::PauliString::PauliOp op)
		{
			pauliString.setOperatorForQubit(qubit, op);
		}

		PauliString::PauliString::PauliOp getOperatorForQubit(size_t qubit) const
		{
			return pauliString.getOperatorForQubit(qubit);
		}

		size_t Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			const size_t nrQubits = reg.getNrQubits();

			// switch to proper basis
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
				if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opX)
					measurementBasis.switchToXBasis(reg, qubit);
				else if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opY)
					measurementBasis.switchToYBasis(reg, qubit);
				// for Z we don't have to do anything, that one is the computational basis

			BaseClass::Execute(reg);

			// switch back to computational basis
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
				if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opX)
					measurementBasis.switchToComputationalFromXBasis(reg, qubit);
				else if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opY)
					measurementBasis.switchToComputationalFromYBasis(reg, qubit);
				// for Z we don't have to do anything, that one is the computational basis

			return 0; // not used here, measurements need to be done explicitely
		}

		double getCoefficient() const
		{
			return pauliString.getCoefficient();
		}

		void setCoefficient(double c)
		{
			pauliString.setCoefficient(c);
		}

	protected:
		QC::MeasurementBasis<VectorClass, MatrixClass> measurementBasis;

		PauliString::PauliString pauliString;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliDecomposedHamiltonianSimulation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		PauliDecomposedHamiltonianSimulation(size_t N = 3, double t = 0, size_t nrSteps = 1, bool addSeed = false)
			: BaseClass(N, addSeed), simTime(t), steps(nrSteps)
		{
		}

		void AddTerm(double coeff, const PauliStringSimulation<VectorClass, MatrixClass>& term)
		{
			terms.push_back(term);
			terms.back().setCoefficient(coeff);
		}

		void Clear()
		{
			terms.clear();
		}

		size_t Execute() override
		{
			const double deltat = simTime / steps;

			for (size_t i = 0; i < terms.size(); ++i)
				terms[i].setTimeInterval(deltat * terms[i].getCoefficient());

			// Suzuki-Trotter expansion
			for (size_t step = 0; step < steps; ++step)
			{
				// this is slow and it can be improved by obtaining the whole matrix once and just by 'applying' the operator each step
				// for now I'll leave it as it is but the gains would be quite big here
				for (PauliStringSimulation<VectorClass, MatrixClass>& term : terms)
					term.Execute(BaseClass::reg);
			}

			return 0; // don't measure for this one, must be done explicitely, needed to be able to check the register first
		}

		void setSimulationTime(double t)
		{
			simTime = t;
		}

		double getSimulationTime() const
		{
			return simTime;
		}

		void setNrSteps(size_t nrSteps)
		{
			steps = nrSteps;
		}

		size_t getNrSteps() const
		{
			return steps;
		}

		// this is here only for testing purposes
		MatrixClass getEvolutionOperator(double t) const
		{
			static const QC::Gates::PauliXGate<MatrixClass> X;
			static const QC::Gates::PauliYGate<MatrixClass> Y;
			static const QC::Gates::PauliZGate<MatrixClass> Z;

			MatrixClass H = MatrixClass::Zero(BaseClass::getNrBasisStates(), BaseClass::getNrBasisStates());

			for (size_t i = 0; i < terms.size(); ++i)
			{
				MatrixClass tm;

				if (terms[i].getOperatorForQubit(0) == PauliString::PauliString::PauliOp::opX)
					tm = X.getRawOperatorMatrix();
				else if (terms[i].getOperatorForQubit(0) == PauliString::PauliString::PauliOp::opY)
					tm = Y.getRawOperatorMatrix();
				else
					tm = Z.getRawOperatorMatrix();

				for (size_t qubit = 1; qubit < BaseClass::getNrQubits(); ++qubit)
				{
					if (terms[i].getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opX)
						tm = Eigen::kroneckerProduct(X.getRawOperatorMatrix(), tm).eval();
					else if (terms[i].getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opY)
						tm = Eigen::kroneckerProduct(Y.getRawOperatorMatrix(), tm).eval();
					else
						tm = Eigen::kroneckerProduct(Z.getRawOperatorMatrix(), tm).eval();
				}

				H += terms[i].getCoefficient() * tm;
			}

			return (std::complex<double>(0, -1) * t * H).exp();
		}

	protected:
		double simTime;
		size_t steps;

		std::vector<PauliStringSimulation<VectorClass, MatrixClass>> terms;
	};
}
