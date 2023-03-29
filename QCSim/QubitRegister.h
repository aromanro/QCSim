#pragma once

#include <math.h>
#include <Eigen/eigen>

#include <random>
#include <complex>
#include <chrono>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "QuantumGate.h"

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegister
	{
	public:
		using GateClass = Gates::QuantumGate<MatrixClass>;

		QubitRegister(int N = 3, int addseed = 0)
			: NrQubits(N), NrBasisStates(1u << NrQubits),
			uniformZeroOne(0, 1)
		{
			assert(N > 0);

			registerStorage = VectorClass::Zero(NrBasisStates);

			uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			timeSeed += addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };

			rng.seed(seed);
		}

		virtual ~QubitRegister() {}

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

		std::complex<double> getRawAmplitude(unsigned int State) const
		{
			if (State >= NrBasisStates) return 0.;

			return registerStorage(State);
		}

		void Clear()
		{
			registerStorage.setZero();
		}

		void Normalize()
		{
			const double accum = (registerStorage.adjoint() * registerStorage)(0).real();
			if (accum < 1E-20) return;

			registerStorage *= 1. / sqrt(accum);
		}

		// to be able to compare different results
		void AdjustPhaseAndNormalize()
		{
			std::complex<double> v0 = registerStorage[0];
			double av0 = abs(v0);

			if (av0 >= 1E-5)
			{
				for (unsigned int i = 0; i < getNrBasisStates(); ++i)
					registerStorage[i] /= v0;
			}
			else
			{
				v0 = registerStorage[NrBasisStates/2];
				av0 = abs(v0);

				if (av0 >= 1E-5)
				{
					for (unsigned int i = 0; i < getNrBasisStates(); ++i)
						registerStorage[i] /= v0;
				}
				else
				{
					v0 = registerStorage[NrBasisStates - 1];
					av0 = abs(v0);

					if (av0 >= 1E-5)
					{
						for (unsigned int i = 0; i < getNrBasisStates(); ++i)
							registerStorage[i] /= v0;
					}
				}
			}

			Normalize();
		}

		unsigned int Measure()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
			double accum = 0;
			unsigned int state = NrBasisStates - 1;
			for (unsigned int i = 0; i < NrBasisStates; ++i)
			{
				accum += (std::conj(registerStorage(i)) * registerStorage(i)).real();
				if (prob <= accum)
				{
					state = i;
					break;
				}
			}

			setToBasisState(state); // collapse

			return state;
		}


		std::map<unsigned int, unsigned int> RepeatedMeasure(unsigned int nrTimes = 1000)
		{
			std::map<unsigned int, unsigned int> measurements;

			const VectorClass initReg = getRegisterStorage();

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				const unsigned int res = Measure();
				++measurements[res];
				setRegisterStorage(initReg);
			}

			return measurements;
		}

		std::pair<unsigned int, unsigned int> RepeatedMeasure(unsigned int qubit, unsigned int nrTimes = 1000)
		{
			unsigned int val1 = 0;
			unsigned int val2 = 0;

			const VectorClass initReg = getRegisterStorage();

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				const unsigned int res = Measure(qubit);
				if (res) ++val2;
				else val1;
				setRegisterStorage(initReg);
			}

			return std::make_pair(val1, val2);
		}

		std::map<unsigned int, unsigned int> RepeatedMeasure(unsigned int firstQubit, unsigned int secondQubit, unsigned int nrTimes = 1000)
		{
			std::map<unsigned int, unsigned int> measurements;

			const VectorClass initReg = getRegisterStorage();

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				const unsigned int res = Measure(firstQubit, secondQubit);
				++measurements[res];
				setRegisterStorage(initReg);
			}

			return measurements;
		}

		// shortcut for measuring a single qubit
		unsigned int Measure(unsigned int qubit)
		{
			return Measure(qubit, qubit);
		}

		// measure a 'subregister' as a separate register
		// can measure a single qubit, if firstQubit == secondQubit
		// will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on position 0 and so on)
		// so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other measured ones are zero, 2 means that the next one 1 one and all others are zero and so on

		unsigned int Measure(unsigned int firstQubit, unsigned int secondQubit)
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
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
				if (prob <= accum)
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

		// controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for gates with a lower number of qubits
		void ApplyGate(const GateClass& gate, unsigned int qubit, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0)
		{
			registerStorage = gate.getOperatorMatrix(NrQubits, qubit, controllingQubit1, controllingQubit2) * registerStorage;
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			registerStorage = m * registerStorage;
		}

		const VectorClass& getRegisterStorage() const
		{
			return registerStorage;
		}

		void setRegisterStorage(const VectorClass& vals)
		{
			if (registerStorage.size() != vals.size()) return;

			registerStorage = vals;
			Normalize();
		}

		// to check how well the computed state matches some 'exact' known one
		double stateFidelity(const VectorClass& state) const
		{
			if (registerStorage.size() != state.size()) return 0;

			const std::complex<double> p = (registerStorage.adjoint() * state)(0);

			return (conj(p) * p).real();
		}

		// allows saving it into a file to look at the data
		// for example to check the results of a quantum simulation
		bool writeToFile(const std::string& name, bool amplitude = true, bool append = false) const
		{
			try {
				std::ofstream thefile;
				thefile.open(name, std::ios::out | (append ? std::ios::app : std::ios::trunc));

				if (!thefile.is_open()) return false;

				if (append) thefile << std::endl << std::endl;

				for (unsigned int i = 0; i < NrBasisStates; ++i)
				{
					thefile << i << "\t"; 
					if (amplitude) thefile << std::abs(registerStorage(i));
					else thefile << registerStorage(i);
					thefile << std::endl;
				}

				return true;
			}
			catch (...) {};

			return false;
		}

		void displayRegister() const
		{
			const unsigned int nQubits = getNrQubits();
			const unsigned int nStates = getNrBasisStates();

			//std::cout << std::fixed;
			std::cout << std::setprecision(4);

			for (unsigned int state = 0; state < nStates; ++state)
			{
				const std::complex<double> val = getBasisStateAmplitude(state);
				if (abs(real(val)) < 1E-10 && abs(imag(val)) < 1E-10) continue;

				bool r = false;
				if (abs(real(val)) > 1E-10) {
					std::cout << real(val) << " ";
					r = true;
				}
				if (abs(imag(val)) > 1E-10) {
					if (r && imag(val) > 0) std::cout << "+ ";

					if (imag(val) < 0) {
						std::cout << "-";
						if (r) std::cout << " ";
					}
					std::cout << abs(imag(val)) << "i ";
				}

				std::cout << "|";

				unsigned int mask = 1 << (nQubits - 1);
				for (unsigned int qubit = 0; qubit < nQubits; ++qubit)
				{
					std::cout << ((state & mask) ? "1" : "0");
					mask >>= 1;
				}

				std::cout << ">    ";
			}
		}

	protected:
		unsigned int NrQubits;
		unsigned int NrBasisStates;

		VectorClass registerStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;
	};

}



