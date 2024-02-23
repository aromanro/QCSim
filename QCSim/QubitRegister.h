#pragma once

#include "QubitRegisterCalculator.h"

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegister : private QubitRegisterCalculator<VectorClass, MatrixClass>
	{
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
		using BaseClass = QubitRegisterCalculator<VectorClass, MatrixClass>;

		QubitRegister(unsigned int N = 3, int addseed = 0)
			: NrQubits(N), NrBasisStates(1u << NrQubits), uniformZeroOne(0, 1), recordGates(false)
		{
			assert(N > 0);

			registerStorage = VectorClass::Zero(NrBasisStates);
			resultsStorage = VectorClass::Zero(NrBasisStates);

			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);

			registerStorage(0) = 1;
		}

		unsigned int getNrQubits() const { return NrQubits; };
		unsigned int getNrBasisStates() const { return NrBasisStates; };

		std::complex<double> getBasisStateAmplitude(unsigned int State) const {
			if (State >= NrBasisStates) return 0;

			return registerStorage(State);
		}

		double getBasisStateProbability(unsigned int State) const {
			if (State >= NrBasisStates) return 0;

			return std::norm(registerStorage(State));
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

			Clear();
			registerStorage(1u << q) = 1;
		}

		// measurement should give either all 0 or all 1
		void setToCatState()
		{
			Clear();
			static const double OneOverSqrt2 = 1. / sqrt(2.);

			registerStorage(0) = OneOverSqrt2;
			registerStorage(NrBasisStates - 1) = OneOverSqrt2;
		}

		// all states have equal amplitude, so measurement should give any state with equal probability
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
				v0 = registerStorage[NrBasisStates >> 1];
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

		unsigned int MeasureAll()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
			double accum = 0;
			unsigned int state = NrBasisStates - 1;

			for (unsigned int i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					state = i;
					break;
				}
			}

			setToBasisState(state); // collapse

			return state;
		}

		// shortcut for measuring a single qubit
		unsigned int MeasureQubit(unsigned int qubit)
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

			if (firstQubit == secondQubit)
			{
				if (NrBasisStates < 16384 + 8192)
					return BaseClass::MeasureQubit(NrBasisStates, registerStorage, firstQubit, prob);

				return BaseClass::MeasureQubitOmp(NrBasisStates, registerStorage, firstQubit, prob);
			}

			return  BaseClass::Measure(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}


		std::map<unsigned int, unsigned int> RepeatedMeasure(unsigned int nrTimes = 1000)
		{
			std::map<unsigned int, unsigned int> measurements;

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				const unsigned int res = MeasureNoCollapse();
				++measurements[res];
			}

			return measurements;
		}

		std::map<unsigned int, unsigned int> RepeatedMeasure(unsigned int firstQubit, unsigned int secondQubit, unsigned int nrTimes = 1000)
		{
			std::map<unsigned int, unsigned int> measurements;

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				const unsigned int res = MeasureNoCollapse(firstQubit, secondQubit);
				++measurements[res];
			}

			return measurements;
		}



		// controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for gates with a lower number of qubits
		void ApplyGate(const GateClass& gate, unsigned int qubit, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0)
		{
#define OPTIMIZED_TENSOR_PRODUCT 1
#ifdef OPTIMIZED_TENSOR_PRODUCT

			const unsigned int gateQubits = gate.getQubitsNumber();

			assert(gateQubits > 0 && gateQubits <= 3);

			const unsigned int qubitBit = 1u << qubit;

			const MatrixClass& gateMatrix = gate.getRawOperatorMatrix();

			// TODO: perhaps also optimize better for controlled gates
			// TODO: there are ways to optimize further for particular kind of gates, probably I won't bother since it's only a constant factor reduction
			
			if (gateQubits == 1)
			{
				if (NrBasisStates < 16384 + 8192)
					BaseClass::ApplyOneQubitGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, NrBasisStates);
				else
					BaseClass::ApplyOneQubitGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, NrBasisStates);
			}
			else if (gateQubits == 2)
			{
				const unsigned int ctrlQubitBit = 1u << controllingQubit1;

				if (NrBasisStates < 8192 + 4096)
					BaseClass::ApplyTwoQubitsGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
				else
					BaseClass::ApplyTwoQubitsGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
			}
			else
			{
				const unsigned int qubitBit2 = 1u << controllingQubit1;
				const unsigned int ctrlQubitBit = 1u << controllingQubit2;

				if (NrBasisStates < 4096 + 2048)
					BaseClass::ApplyThreeQubitsGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
				else
					BaseClass::ApplyThreeQubitsGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
			}

			registerStorage.swap(resultsStorage);
#else			
			registerStorage = gate.getOperatorMatrix(NrQubits, qubit, controllingQubit1, controllingQubit2) * registerStorage;
#endif

			if (recordGates)
				computeGates.emplace_back(Gates::AppliedGate(gate.getRawOperatorMatrix(), qubit, controllingQubit1, controllingQubit2));
		}

		void ApplyGate(const Gates::AppliedGate<MatrixClass>& gate)
		{
			ApplyGate(gate, gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
		}

		void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>>& gates)
		{
			for (const auto& gate : gates)
				ApplyGate(gate);
		}

		void ApplyOperatorMatrix(const MatrixClass& m)
		{
			registerStorage = m * registerStorage;

			if (recordGates)
				computeGates.emplace_back(Gates::AppliedGate(m));
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

		// warning, you should be sure that the vector is normalized and has the proper size
		void setRegisterStorageFastNoNormalize(VectorClass& vals)
		{
			registerStorage.swap(vals);
		}

		// to check how well the computed state matches some 'exact' known one
		double stateFidelity(const VectorClass& state) const
		{
			if (registerStorage.size() != state.size()) return 0;

			const std::complex<double> p = (registerStorage.adjoint() * state)(0);

			return norm(p);
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


		void displayState(unsigned int state) const
		{
			const unsigned int nQubits = getNrQubits();
			std::cout << "|";

			unsigned int mask = 1 << (nQubits - 1);
			for (unsigned int qubit = 0; qubit < nQubits; ++qubit)
			{
				std::cout << ((state & mask) ? "1" : "0");
				mask >>= 1;
			}

			std::cout << ">    ";
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

				displayState(nQubits);
			}
		}

		void ComputeStart()
		{
			recordGates = true;
			computeGates.clear();
		}

		void ComputeEnd()
		{
			recordGates = false;
		}

		void ComputeClear()
		{
			computeGates.clear();
		}

		// applies again the recorded gates
		// with this the same operations can be repeated several times
		void Compute()
		{
			// avoid recording the gates again if somehow the user forgot to call ComputeEnd
			const bool recordSave = recordGates;
			recordGates = false;

			for (const Gates::AppliedGate<MatrixClass>& gate : computeGates)
			{
				if (gate.getQubitsNumber() > 3)
					ApplyOperatorMatrix(gate.getRawOperatorMatrix());
				else
					ApplyGate(gate);
			}

			recordGates = recordSave;
		}

		// undoes the recorded gates
		// the operations are unitary, so U^-1 = U^t and (U1 * U2)^t = U2^t * U1^t 
		void Uncompute()
		{
			const bool recordSave = recordGates;
			recordGates = false;

			for (auto it = computeGates.crbegin(); it != computeGates.crend(); ++it)
			{
				if (it->getQubitsNumber() > 3)
					ApplyOperatorMatrix(it->getRawOperatorMatrix().adjoint());
				else
				{
					Gates::AppliedGate<MatrixClass> gate(it->getRawOperatorMatrix().adjoint(), it->getQubit1(), it->getQubit2(), it->getQubit3());
					ApplyGate(gate);
				}
			}

			recordGates = recordSave;
		}

	protected:
		// the following ones should be used for 'repeated measurements' that avoid reexecuting the circuit each time
		unsigned int MeasureNoCollapse()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
			double accum = 0;
			unsigned int state = NrBasisStates - 1;
			for (unsigned int i = 0; i < NrBasisStates; ++i)
			{
				accum += norm(registerStorage(i));
				if (prob <= accum)
				{
					state = i;
					break;
				}
			}

			return state;
		}

		// shortcut for measuring a single qubit
		unsigned int MeasureNoCollapse(unsigned int qubit)
		{
			return MeasureNoCollapse(qubit, qubit);
		}

		// measure a 'subregister' as a separate register
		// can measure a single qubit, if firstQubit == secondQubit
		// will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on position 0 and so on)
		// so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other measured ones are zero, 2 means that the next one 1 one and all others are zero and so on

		unsigned int MeasureNoCollapse(unsigned int firstQubit, unsigned int secondQubit)
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 

			if (firstQubit == secondQubit)
			{
				if (NrBasisStates < 16384 + 8192)
					return BaseClass::MeasureQubitNoCollapse(NrBasisStates, registerStorage, firstQubit, prob);

				return BaseClass::MeasureQubitNoCollapseOmp(NrBasisStates, registerStorage, firstQubit, prob);
			}
			
			return BaseClass::MeasureNoCollapse(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}

		unsigned int NrQubits;
		unsigned int NrBasisStates;

		VectorClass registerStorage;
		VectorClass resultsStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;

		std::vector<Gates::AppliedGate<MatrixClass>> computeGates;
		bool recordGates;
	};

}



