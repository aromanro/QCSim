#pragma once

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"
#include "QuantumFourierTransform.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QuantumSimulation {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SchrodingerSimulation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		SchrodingerSimulation(unsigned int N = 5, double t = 0, unsigned int nrSteps = 1, bool addSeed = false)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addSeed), simTime(t), steps(nrSteps), deltax(0.01), fourier(N, 0, N - 1)
		{
			potential.resize(QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates(), 0.);
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0); // useless if the initialization of the start state is done, but better be safe...
		}

		unsigned int Execute() override
		{
			const double deltat = simTime / steps;

			Init(deltat);
			
			// Suzuki-Trotter expansion
			for (unsigned int step = 0; step < steps; ++step)
			{
				ApplyHalfPotentialOperatorEvolution();
				fourier.QFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);
				ApplyKineticOperatorEvolution();
				fourier.IQFT(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg);
				ApplyHalfPotentialOperator();
			}

			return 0;
		}

		// the index is from 0 to nrStates, but the position is going to be from -nrStates/2 * deltax to nrStates/2 * deltax
		double getPotential(unsigned int pos) const
		{
			if (pos >= potential.size()) return 0.;

			return potential[pos];
		}

		void setPotentiak(unsigned int pos, double val)
		{
			if (pos >= potential.size()) return;

			potential[pos] = val;
		}

		double getDeltax() const
		{
			return deltax;
		}

		double setDeltax(double v)
		{
			deltax = v;
		}

		void setSimulationTime(double t)
		{
			simTime = t;
		}

		void setNrSteps(unsigned int nrSteps)
		{
			steps = nrSteps;
		}

	protected:
		void Init(double deltat)
		{
			// ensure that the starting wavefunction is normalized
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Normalize();

			// the potential and kinetic operators can be constructed from controlled phase shift gates
			// I'm not going to do that, at least not yet, because it's going to be too slow
		
			// for how to do that see for example "Quantum simulation of the single-particle Schrodinger equation" by Giuliano Benenti, Giuliano Strini
			// https://arxiv.org/abs/0709.1704

			const unsigned int nrStates = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates();
			kineticOp = MatrixClass::Identity(nrStates, nrStates);
			potentialOp = MatrixClass::Identity(nrStates, nrStates);

			// they are diagonal, the potential is in position basis, the kinetic one is in k-basis, where it's switched using the quantum fourier transform

			const double halfX = deltax * 0.5 * (nrStates - 1);

			for (int i = 0; i < nrStates; ++i)
			{
				const double x = deltax * i - halfX;
				kineticOp(i) = std::exp(std::complex<double>(0, -0.5) * x * x * deltat);
				potentialOp(i) = std::exp(std::complex<double>(0, -0.5) * potential[i] * deltat); // the reason of 0.5 here is that I'm using a Suzuki-Trotter expansion with a better precision than the one used for Pauli strings, see 'Execute'
			}
		}

		void ApplyHalfPotentialOperatorEvolution()
		{
			// with a single operator is simple, it would be quite annoying with a lot of quantum gates
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyOperatorMatrix(potentialOp);
		}

		void ApplyKineticOperatorEvolution()
		{
			// with a single operator is simple, it would be quite annoying with a lot of quantum gates
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyOperatorMatrix(kineticOp);
		}

		double simTime;
		unsigned int steps;
		double deltax;

		QC::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
		std::vector<double> potential;

		MatrixClass kineticOp;
		MatrixClass potentialOp;
	};
}


