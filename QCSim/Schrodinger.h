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
		SchrodingerSimulation(unsigned int N = 8, double t = 1, double dx = 0.1, unsigned int nrSteps = 100, bool addSeed = false)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addSeed), simTime(t), steps(nrSteps), deltax(dx), fourier(N, 0, N - 1)
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
				ApplyHalfPotentialOperatorEvolution();
			}

			return 0;
		}

		// the index is from 0 to nrStates, but the position is going to be from -nrStates/2 * deltax to nrStates/2 * deltax
		double getPotential(unsigned int pos) const
		{
			if (pos >= potential.size()) return 0.;

			return potential[pos];
		}

		void setPotential(unsigned int pos, double val)
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

		// use it to set a square barrier (positive val) or well (negative val) in the middle
		void setConstantPotentialInTheMiddle(double val, unsigned int halfwidth)
		{
			const int nrStates = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates();
			const int halfStates = nrStates / 2;

			for (int i = std::max<int>(1, halfStates - halfwidth); i < std::min<int>(nrStates - 1, halfStates + halfwidth); ++i)
				setPotential(i, val);
		}

		// set a gaussian wavefunction in the register, k makes it 'move'
		void setGaussian(unsigned int pos, double stdev, double k)
		{
			const int nrStates = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates();
			const double a = 1. / (stdev * sqrt(2. * M_PI));

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Clear();
			for (int i = 1; i < nrStates - 1; ++i) // the values at the ends should stay zero
			{
				const double e = (static_cast<double>(i) - static_cast<double>(pos)) / stdev;
				const std::complex<double> val = a * exp(std::complex<double>(-0.5 * e * e, -k * deltax * i));
				QC::QuantumAlgorithm<VectorClass, MatrixClass>::setRawAmplitude(i, val);
			}

			QC::QuantumAlgorithm<VectorClass, MatrixClass>::Normalize();
		}

		//*****************************************************************************************************************************************************
		// Code for finite differences solver - to be used to compare the results from the 'quantum computation'
		//*****************************************************************************************************************************************************

		// for how it could be implemented see: "Computer-Generated Motion Pictures of One-Dimensional Quantum-Mechanical Transmission and Reflection Phenomena" by Abraham Goldberg and Harry M. Schey
		// https://aapt.scitation.org/doi/10.1119/1.1973991

		// the problem with some naive application of finite differences is that the resulting numerical time evolution is not unitary anymore
		// so it needs some special care... there will be some factors of 2 different between that paper and this code, because they consider 2 * m = 1, while I use m = 1 so 1/2 for kinetical term does not dissapear

		void solveWithFiniteDifferences()
		{
			const unsigned int nrStates = QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrBasisStates();
			const double deltat = simTime / steps;
			const double eps2 = deltax * deltax;
			const double lambda = 2. * eps2 / deltat;
			const std::complex<double> ilambda = std::complex<double>(0, 1) * lambda;
			const std::complex<double> oneplusilambda = std::complex<double>(1., 0.) + ilambda;
			const std::complex<double> oneminusilambda = std::complex<double>(1., 0.) - ilambda;

			std::vector<std::complex<double>> e(nrStates, 0);
			std::vector<std::complex<double>> f(nrStates, 0);
			std::vector<std::complex<double>> omegaterm(nrStates, 0);

			// not time dependent so it can be computed once at the beginning
			e[1] = oneminusilambda + eps2 * potential[1];
			for (unsigned int i = 2; i < nrStates - 1; ++i)
			{
				const double eps2pot = eps2 * potential[i];
				e[i] = oneminusilambda + eps2pot - 1. / (2. * e[i - 1]);
				omegaterm[0] = 2. * (oneplusilambda + eps2pot);
			}

			for (unsigned int step = 0; step < steps; ++step)
			{
				// compute the needed values first
				f[1] = omegaterm[1] * QC::QuantumAlgorithm<VectorClass, MatrixClass>::getBasisStateAmplitude(1) - QC::QuantumAlgorithm<VectorClass, MatrixClass>::getBasisStateAmplitude(2);
				for (unsigned int i = 2; i < nrStates - 1; ++i)
				{
					const std::complex<double> omega = omegaterm[i] * QC::QuantumAlgorithm<VectorClass, MatrixClass>::getBasisStateAmplitude(i) - (QC::QuantumAlgorithm<VectorClass, MatrixClass>::getBasisStateAmplitude(i + 1) + QC::QuantumAlgorithm<VectorClass, MatrixClass>::getBasisStateAmplitude(i - 1));
					f[i] = omega + f[i - 1] / e[i - 1];
				}

				// now compute the wavefunction
				// the limits will stay at zero
				for (unsigned int i = nrStates - 2; i > 0; --i)
				{
					const std::complex<double> val = (QC::QuantumAlgorithm<VectorClass, MatrixClass>::getBasisStateAmplitude(i + 1) - f[i]) / e[i];
					QC::QuantumAlgorithm<VectorClass, MatrixClass>::setRawAmplitude(i, val);
				}
			}
		}

		//*****************************************************************************************************************************************************

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

			for (unsigned int i = 0; i < nrStates; ++i)
			{
				const double x = deltax * i - halfX;
				kineticOp(i, i) = std::exp(std::complex<double>(0, -0.5) * x * x * deltat);
				potentialOp(i, i) = std::exp(std::complex<double>(0, -0.5) * potential[i] * deltat); // the reason of 0.5 here is that I'm using a Suzuki-Trotter expansion with a better precision than the one used for Pauli strings, see 'Execute'
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

		//*****************************************************************************************************************************************************
		// Code for finite differences solver - to be used to compare the results from the 'quantum computation'
		//*****************************************************************************************************************************************************




		//*****************************************************************************************************************************************************

		double simTime;
		unsigned int steps;
		double deltax;

		QC::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
		std::vector<double> potential;

		MatrixClass kineticOp;
		MatrixClass potentialOp;
	};
}


