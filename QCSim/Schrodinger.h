#pragma once

#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"
#include "QuantumFourierTransform.h"

#include "FFT.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace QuantumSimulation {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SchrodingerSimulation :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		SchrodingerSimulation(size_t N = 8, double dt = 0.1, double dx = 0.1, size_t nrSteps = 50, bool addSeed = false)
			: BaseClass(N, addSeed), deltat(dt), steps(nrSteps), deltax(dx), fourier(N, 0, N - 1)
		{
			potential.resize(BaseClass::getNrBasisStates(), 0.);
			BaseClass::setToBasisState(0); // useless if the initialization of the start state is done, but better be safe...
		}

		size_t Execute() override
		{
			Init();
			
			// Suzuki-Trotter expansion
			for (size_t step = 0; step < steps; ++step)
			{
				ApplyPotentialOperatorEvolution();
				fourier.IQFT(BaseClass::reg);
				ApplyKineticOperatorEvolution();
				fourier.QFT(BaseClass::reg);
			}

			return 0;
		}

		// the index is from 0 to nrStates, but the position is going to be from -nrStates/2 * deltax to nrStates/2 * deltax
		double getPotential(size_t pos) const
		{
			if (pos >= potential.size()) return 0.;

			return potential[pos];
		}

		void setPotential(size_t pos, double val)
		{
			if (pos >= potential.size()) return;

			potential[pos] = val;
		}

		double getDeltax() const
		{
			return deltax;
		}

		void setDeltax(double v)
		{
			deltax = v;
		}

		void setSimulationTimeInterval(double dt)
		{
			deltat = dt;
		}

		void setNrSteps(size_t nrSteps)
		{
			steps = nrSteps;
		}

		// use it to set a square barrier (positive val) or well (negative val) in the middle
		void setConstantPotentialInTheMiddle(double val, size_t halfwidth)
		{
			const int nrStates = static_cast<int>(BaseClass::getNrBasisStates());
			const int halfStates = nrStates / 2;

			for (int i = std::max<int>(1, halfStates - halfwidth); i < std::min<int>(nrStates - 1, halfStates + halfwidth); ++i)
				setPotential(i, val);

			// also set high values at the limits
			setPotential(0, 1E150);
			setPotential(nrStates - 1, 1E150);
		}

		// use it to set a step potential starting in the middle
		void setConstantPotentialToRight(double val)
		{
			const int nrStates = static_cast<int>(BaseClass::getNrBasisStates());
			const int halfStates = nrStates / 2;

			for (int i = halfStates; i < nrStates; ++i)
				setPotential(i, val);

			// also set high values at the limits
			setPotential(0, 1E15);
			setPotential(nrStates - 1, 1E150);
		}

		// set a gaussian wavefunction in the register, k makes it 'move'
		// set it to the left side
		void setGaussian(size_t pos, double stdev, double k)
		{
			BaseClass::Clear();

			const int nrStates = static_cast<int>(BaseClass::getNrBasisStates());

			stdev *= deltax;
			const double a = 1. / (stdev * sqrt(2. * M_PI));
			const double halfX = 0.5 * deltax * (nrStates - 1.);
			const double m = deltax * pos - halfX;

			for (int i = 1; i < nrStates - 1; ++i) // the values at the ends should stay zero
			{
				const double x = deltax * i - halfX;
				const double e = (x - m) / stdev;
			
				// momentum operator is -i d/dx
				// so it brings down ik to obtain k 
				const std::complex<double> val = a * exp(std::complex<double>(-0.5 * e * e, k * x));
				BaseClass::setRawAmplitude(i, val);
			}

			BaseClass::Normalize();
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
			const size_t nrStates = BaseClass::getNrBasisStates();
			const double eps2 =  2 * deltax * deltax; // the easiest way to add 1/2. for kinetic term, just having multiplied with 2 here
			const double lambda = 2. * eps2 / deltat;
			const std::complex<double> ilambda = std::complex<double>(0, lambda);
			const std::complex<double> twoplusilambda = 2. + ilambda; 
			const std::complex<double> twominusilambda = 2. - ilambda;

			std::vector<std::complex<double>> e(nrStates, 0);
			std::vector<std::complex<double>> f(nrStates, 0);
			std::vector<std::complex<double>> omegaterm(nrStates, 0);

			// not time dependent so it can be computed once at the beginning
			double eps2pot = eps2 * potential[1];
			e[1] = twominusilambda + eps2pot;
			omegaterm[1] = twoplusilambda + eps2pot;
			for (size_t i = 2; i < nrStates - 1; ++i)
			{
				eps2pot = eps2 * potential[i];
				// see eq (17)
				e[i] = twominusilambda + eps2pot - 1. / e[i - 1];
				omegaterm[i] = twoplusilambda + eps2pot; // see the def of omega, above eq (14)
			}

			for (size_t step = 0; step < steps; ++step)
			{
				// compute the needed values first
				f[1] = omegaterm[1] * BaseClass::getBasisStateAmplitude(1) - BaseClass::getBasisStateAmplitude(2);
				for (size_t i = 2; i < nrStates - 1; ++i)
				{
					// see eq (18)
					const std::complex<double> omega = omegaterm[i] * BaseClass::getBasisStateAmplitude(i) - (BaseClass::getBasisStateAmplitude(i + 1) + BaseClass::getBasisStateAmplitude(i - 1));
					f[i] = omega + f[i - 1] / e[i - 1];
				}

				// now compute the wavefunction
				// the limits will stay at zero
				for (size_t i = nrStates - 2; i > 0; --i)
				{
					// see eq (20)
					const std::complex<double> val = (BaseClass::getBasisStateAmplitude(i + 1) - f[i]) / e[i];
					BaseClass::setRawAmplitude(i, val);
				}
			}
		}

		//*****************************************************************************************************************************************************
		// Code for a dumb solver for checks, the above is hard to have values that overlap with the range where the schrodinger solver works
		//*****************************************************************************************************************************************************

		void solveWithFiniteDifferencesSimple()
		{
			InitSimpleFiniteDifferences();

			for (size_t step = 0; step < steps; ++step)
			{
				BaseClass::ApplyOperatorMatrix(evolutionOp);
				ApplyPotentialOperatorEvolution();
				BaseClass::Normalize(); // evolution above is not exactly unitary
			}
		}

		void solveWithFiniteDifferencesSimpler()
		{
			// ensure that the starting wavefunction is normalized
			BaseClass::Normalize();

			const size_t nrStates = BaseClass::getNrBasisStates();
			const size_t nrStatesM1 = nrStates - 1;
			const std::complex<double> j(0, 1);
			const double eps = 0.5 * deltat / (deltax * deltax);
			
			VectorClass newPsi = VectorClass::Zero(nrStates);

			for (size_t step = 0; step < steps; ++step)
			{	
				const VectorClass& oldPsi = BaseClass::getRegisterStorage();

				for (size_t i = 1; i < nrStatesM1; ++i)
					newPsi(i) = oldPsi(i) + j * eps * (oldPsi(i - 1) - 2. * oldPsi(i) + oldPsi(i + 1)) - j * potential[i] * deltat * oldPsi(i);
				
				BaseClass::setRegisterStorage(newPsi); // this also normalizes the wavefunction
			}
		}

		// ok, I gave up, it's hard to set values that are ok for finite differences methods and for the schrodinger simulation, so I went with 'classical' FFT

		void solveWithClassicalFFT()
		{
			Init();

			Fourier::FFT fft;

			// Suzuki-Trotter expansion
			for (size_t step = 0; step < steps; ++step)
			{
				ApplyPotentialOperatorEvolution();

				{
					const VectorClass& oldPsi = BaseClass::getRegisterStorage();

					std::vector<std::complex<double>> vec(oldPsi.size());
					Eigen::Map<VectorClass>(vec.data(), oldPsi.size()) = oldPsi;

					std::vector<std::complex<double>> newvec(oldPsi.size());
					fft.fwd(vec.data(), newvec.data(), static_cast<size_t>(vec.size()));

					const VectorClass newPsi = Eigen::Map<VectorClass>(newvec.data(), newvec.size());
					BaseClass::setRegisterStorage(newPsi);
				}

				ApplyKineticOperatorEvolution();
				
				{
					const VectorClass& oldPsi = BaseClass::getRegisterStorage();

					std::vector<std::complex<double>> vec(oldPsi.size());
					Eigen::Map<VectorClass>(vec.data(), oldPsi.size()) = oldPsi;

					std::vector<std::complex<double>> newvec(oldPsi.size());
					fft.inv(vec.data(), newvec.data(), static_cast<size_t>(vec.size()));

					const VectorClass newPsi = Eigen::Map<VectorClass>(newvec.data(), newvec.size());
					BaseClass::setRegisterStorage(newPsi);
				}
			}
		}

		//*****************************************************************************************************************************************************

	protected:
		void Init()
		{
			// ensure that the starting wavefunction is normalized
			BaseClass::Normalize();

			// the potential and kinetic operators can be constructed from controlled phase shift gates
			// I'm not going to do that, at least not yet, because it's going to be too slow
		
			// for how to do that see for example "Quantum simulation of the single-particle Schrodinger equation" by Giuliano Benenti, Giuliano Strini
			// https://arxiv.org/abs/0709.1704

			const int nrStates = BaseClass::getNrBasisStates();
			const double halfX = 0.5 * deltax * (nrStates - 1.);
	
			kineticOp = MatrixClass::Zero(nrStates, nrStates);
			potentialOp = MatrixClass::Zero(nrStates, nrStates);

			// they are diagonal, the potential is in position basis, the kinetic one is in k-basis, where it's switched using the quantum fourier transform

			for (int i = 0; i < nrStates; ++i)
			{
				const double x = deltax * i - halfX;
	
				kineticOp(i, i) = std::polar(1., -0.5 * x * x * deltat);
				potentialOp(i, i) = std::polar(1., -potential[i] * deltat);
			}
		}

		void InitSimpleFiniteDifferences()
		{
			// ensure that the starting wavefunction is normalized
			BaseClass::Normalize();

			const int nrStates = BaseClass::getNrBasisStates();
			const double eps = 0.5 * deltat / (deltax * deltax);
			const std::complex<double> t = std::complex<double>(1., -2. * eps);
			
			evolutionOp = MatrixClass::Zero(nrStates, nrStates);
			potentialOp = MatrixClass::Zero(nrStates, nrStates);

			for (int i = 0; i < nrStates; ++i)
			{
				evolutionOp(i, i) = t;
				potentialOp(i, i) = std::polar(1., -potential[i] * deltat);
			}

			const int nrStatesMinusOne = nrStates - 1;
			const std::complex<double> v = std::complex(0., eps);
			for (int i = 0; i < nrStatesMinusOne; ++i)
			{
				evolutionOp(i, i + 1) = v;
				evolutionOp(i + 1, i) = v;
			}
			evolutionOp(0, nrStatesMinusOne) = v;
			evolutionOp(nrStatesMinusOne, 0) = v;
		}

		void ApplyPotentialOperatorEvolution()
		{
			// with a single operator is simple, it would be quite annoying with a lot of quantum gates
			BaseClass::ApplyOperatorMatrix(potentialOp);
		}

		void ApplyKineticOperatorEvolution()
		{
			// with a single operator is simple, it would be quite annoying with a lot of quantum gates
			BaseClass::ApplyOperatorMatrix(kineticOp);
		}

		double deltat;
		size_t steps;
		double deltax;

		QC::SubAlgo::QuantumFourierTransform<VectorClass, MatrixClass> fourier;
		std::vector<double> potential;

		MatrixClass kineticOp;
		MatrixClass potentialOp;

		MatrixClass evolutionOp; // not unitary, but should be 'good enough' for tests
	};
}


