#pragma once

#include <math.h>
#include <Eigen/eigen>

#include <random>
#include <complex>

namespace QC {

	class QubitRegister
	{
	public:
		QubitRegister(int N = 3, int addseed = 0);

		unsigned int getNrQubits() const { return NrQubits; };
		unsigned int getNrBasisStates() const { return NrBasisStates; };

		std::complex<double> getBasisStateAmplitude(unsigned int State) const {
			if (State >= NrBasisStates) return std::complex<double>(0);

			return registerStorage(State);
		}

		void setToBasisState(unsigned int State);
		void setToCatState();
		void setToEqualSuperposition();

		// to be able to set them all, after setting them, call Normalize
		void setRawAmplitude(unsigned int State, std::complex<double> val);

		void Normalize();

		unsigned int Measure();

	protected:
		unsigned int NrQubits;
		unsigned int NrBasisStates;

		Eigen::VectorXcd registerStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;
	};

}

