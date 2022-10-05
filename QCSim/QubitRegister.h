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
			if (State >= NrBasisStates) return 0;

			return registerStorage(State);
		}

		void setToBasisState(unsigned int State)
		{
			if (State >= NrBasisStates) return;

			Clear();
			registerStorage(State) = 1;
		}

		void setToQubitState(unsigned int q)
		{
			if (q >= NrQubits) return;

			unsigned int state = 1u << q;

			Clear();
			registerStorage(state) = 1;
		}

		void setToCatState()
		{
			Clear();
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

		void Clear()
		{
			registerStorage.setZero();
		}

		void Normalize()
		{
			const double accum = (registerStorage.adjoint()*registerStorage)(0).real();
			if (abs(accum) < 1E-20) return;

			registerStorage *= 1. / sqrt(accum);
		}

		unsigned int Measure()
		{
			const double prob = uniformZeroOne(rng);
			double accum = 0;
			unsigned int state = NrBasisStates - 1;
			for (unsigned int i = 0; i < NrBasisStates; ++i)
			{
				accum += (std::conj(registerStorage(i)) * registerStorage(i)).real();
				if (prob < accum)
				{
					state = i;
					break;
				}
			}

			setToBasisState(state); // collapse

			return state;
		}

		// measure a 'subregister' as a separate register
		// can measure a single qubit, if firstQubit == secondQubit
		// will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on position 0 and so on)
		// so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other measured ones are zero, 2 means that the next one 1 one and all others are zero and so on
		
		unsigned int Measure(unsigned int firstQubit, unsigned int secondQubit)
		{
			const double prob = uniformZeroOne(rng);
			double accum = 0;

			const unsigned int secondQubitp1 = secondQubit + 1;

			const unsigned int firstPartMask = (1u << firstQubit) - 1;
			const unsigned int measuredPartMask = (1u << secondQubitp1) - 1 - firstPartMask;
			const unsigned int secondPartMask = NrBasisStates - 1 - measuredPartMask - firstPartMask;

			const unsigned int secondPartMax = secondPartMask >> secondQubitp1;
			const unsigned int maxMeasuredState = measuredPartMask >> firstQubit;

			unsigned int measuredState = maxMeasuredState;
			
			double norm = 1;
			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;
				double stateProbability = 0;

				for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
				{
					const unsigned int secondPart = secondPartBits << secondQubitp1;
					for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
					{
						const unsigned int wholeState = secondPart | stateRegBits | firstPartBits;
						stateProbability += (std::conj(registerStorage[wholeState]) * registerStorage[wholeState]).real();
					}
				}

				accum += stateProbability;
				if (prob < accum)
				{
					measuredState = state;
					norm = 1. / sqrt(stateProbability);
					break;
				}
			}

			//int cnt = 0;
			// collapse
			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;

				if (state == measuredState)
				{
					for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
					{
						const unsigned int secondPart = secondPartBits << secondQubitp1;
						for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						{
							const unsigned int wholeState = secondPart | stateRegBits | firstPartBits;

							registerStorage[wholeState] *= norm;
							//++cnt;
						}
					}
				}
				else
				{
					for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
					{
						const unsigned int secondPart = secondPartBits << secondQubitp1;
						for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						{
							const unsigned int wholeState = secondPart | stateRegBits | firstPartBits;
							
							registerStorage[wholeState] = 0;
							//++cnt;
						}
					}
				}
			}

			//assert(cnt == NrBasisStates);

			return measuredState;
		}

		void ApplyGate(const QuantumGate<MatrixClass>& gate, unsigned int qubit, unsigned int controllingQubit = 0)
		{
			registerStorage = gate.getOperatorMatrix(NrQubits, qubit, controllingQubit) * registerStorage;
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			registerStorage = m * registerStorage;
		}

		VectorClass getRegisterStorage() const
		{
			return registerStorage;
		}

	protected:
		unsigned int NrQubits;
		unsigned int NrBasisStates;

		VectorClass registerStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;
	};

}



