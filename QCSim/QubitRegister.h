#pragma once

#include <math.h>
#include <Eigen/eigen>
#include <Eigen/sparse>


#include <random>
#include <complex>
#include <chrono>

#include "QuantumGate.h"

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegister
	{
	public:
		QubitRegister(int N = 3, int addseed = 0)
			: NrQubits(N), NrBasisStates(1u << NrQubits),
			uniformZeroOne(0, 1)
		{
			registerStorage = VectorClass::Zero(NrBasisStates);

			uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			timeSeed += addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };

			rng.seed(seed);
		}

		unsigned int getNrQubits() const { return NrQubits; };
		unsigned int getNrBasisStates() const { return NrBasisStates; };

		std::complex<double> getBasisStateAmplitude(unsigned int State) const {
			if (State >= NrBasisStates) return std::complex<double>(0);

			return registerStorage(State);
		}

		void setToBasisState(unsigned int State)
		{
			if (State >= NrBasisStates) return;

			registerStorage.setZero();
			registerStorage(State) = 1;
		}

		void setToCatState()
		{
			registerStorage.setZero();
			static const double OneOverSqrt2 = 1. / sqrt(2.);

			registerStorage(0) = OneOverSqrt2;
			registerStorage(NrBasisStates - 1) = OneOverSqrt2;
		}

		void setToEqualSuperposition()
		{
			registerStorage.setConstant(1. / sqrt(NrBasisStates));
		}

		// to be able to set them all, after setting them, call Normalize
		void setRawAmplitude(unsigned int State, std::complex<double> val)
		{
			if (State >= NrBasisStates) return;

			registerStorage(State) = val;
		}

		void Normalize()
		{
			const double accum = (registerStorage.adjoint() * registerStorage)(0).real();
			if (accum < 1E-20) return;

			registerStorage = 1. / sqrt(accum) * registerStorage;
		}

		unsigned int Measure()
		{
			const double prob = uniformZeroOne(rng);
			double accum = 0;
			for (unsigned int i = 0; i < NrBasisStates; ++i)
			{
				accum += (std::conj(registerStorage(i)) * registerStorage(i)).real();
				if (prob < accum)
				{
					setToBasisState(i); // collapse
					return i;
				}
			}
			const unsigned int state = NrBasisStates - 1;
			setToBasisState(state);

			return state;
		}

		void ApplyGate(const QuantumGate<MatrixClass>& gate, unsigned int qubit, unsigned int controllingQubit = 0)
		{
			registerStorage = gate.getOperatorMatrix(NrQubits, qubit, controllingQubit) * registerStorage;
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			registerStorage = m * registerStorage;
		}

	protected:
		unsigned int NrQubits;
		unsigned int NrBasisStates;

		VectorClass registerStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;
	};

}



