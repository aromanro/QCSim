#pragma once

#include "QubitRegisterCalculator.h"

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegister : public QubitRegisterCalculator<VectorClass, MatrixClass>
	{
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
		using BaseClass = QubitRegisterCalculator<VectorClass, MatrixClass>;

		QubitRegister(size_t N = 3, int addseed = 0)
			: NrQubits(N), NrBasisStates(1ULL << NrQubits), 
			registerStorage(VectorClass::Zero(NrBasisStates)),
			uniformZeroOne(0, 1), recordGates(false)
		{
			assert(N > 0);

			resultsStorage.resize(NrBasisStates);

			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);

			registerStorage(0) = 1;
		}

		// this is a special constructor, I need it for something in a derived work (closed source)
		QubitRegister(size_t N, VectorClass& v, int addseed = 0)
			: NrQubits(N), NrBasisStates(1ULL << NrQubits),
			uniformZeroOne(0, 1), recordGates(false)
		{
			assert(N > 0);
			resultsStorage.resize(NrBasisStates);
			registerStorage.swap(v);

			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);
		}

		size_t getNrQubits() const { return NrQubits; };
		size_t getNrBasisStates() const { return NrBasisStates; };

		std::complex<double> getBasisStateAmplitude(size_t State) const {
			if (State >= NrBasisStates) return 0;

			return registerStorage(State);
		}

		double getBasisStateProbability(size_t State) const {
			if (State >= NrBasisStates) return 0;

			return std::norm(registerStorage(State));
		}

		void setToBasisState(size_t State)
		{
			if (State >= NrBasisStates) return;

			Clear();
			registerStorage(State) = 1;
		}

		void setToQubitState(size_t q)
		{
			if (q >= NrQubits) return;

			Clear();
			registerStorage(1ULL << q) = 1;
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
		void setRawAmplitude(size_t State, std::complex<double> val)
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
			const double norm = registerStorage.norm();
			if (norm < 1E-20) return;

			registerStorage *= 1. / norm;
		}

		// to be able to compare different results
		void AdjustPhaseAndNormalize()
		{
			std::complex<double> v0 = registerStorage[0];
			double av0 = abs(v0);

			if (av0 >= 1E-5)
			{
				for (size_t i = 0; i < getNrBasisStates(); ++i)
					registerStorage[i] /= v0;
			}
			else
			{
				v0 = registerStorage[NrBasisStates >> 1];
				av0 = abs(v0);

				if (av0 >= 1E-5)
				{
					for (size_t i = 0; i < getNrBasisStates(); ++i)
						registerStorage[i] /= v0;
				}
				else
				{
					v0 = registerStorage[NrBasisStates - 1];
					av0 = abs(v0);

					if (av0 >= 1E-5)
					{
						for (size_t i = 0; i < getNrBasisStates(); ++i)
							registerStorage[i] /= v0;
					}
				}
			}

			Normalize();
		}

		size_t MeasureAll()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
			double accum = 0;
			size_t state = NrBasisStates - 1;

			for (size_t i = 0; i < NrBasisStates; ++i)
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
		size_t MeasureQubit(size_t qubit)
		{
			return Measure(qubit, qubit);
		}

		// measure a 'subregister' as a separate register
		// can measure a single qubit, if firstQubit == secondQubit
		// will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on position 0 and so on)
		// so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other measured ones are zero, 2 means that the next one 1 one and all others are zero and so on

		size_t Measure(size_t firstQubit, size_t secondQubit)
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 

			if (firstQubit == secondQubit)
			{
				if (NrBasisStates < BaseClass::OneQubitOmpLimit)
					return BaseClass::MeasureQubit(NrBasisStates, registerStorage, firstQubit, prob);

				return BaseClass::MeasureQubitOmp(NrBasisStates, registerStorage, firstQubit, prob);
			}

			return  BaseClass::Measure(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}


		std::map<size_t, size_t> RepeatedMeasure(size_t nrTimes = 1000)
		{
			std::map<size_t, size_t> measurements;

			for (size_t i = 0; i < nrTimes; ++i)
			{
				const size_t res = MeasureNoCollapse();
				++measurements[res];
			}

			return measurements;
		}

		std::map<size_t, size_t> RepeatedMeasure(size_t firstQubit, size_t secondQubit, size_t nrTimes = 1000)
		{
			std::map<size_t, size_t> measurements;

			for (size_t i = 0; i < nrTimes; ++i)
			{
				const size_t res = MeasureNoCollapse(firstQubit, secondQubit);
				++measurements[res];
			}

			return measurements;
		}



		// controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for gates with a lower number of qubits
		void ApplyGate(const GateClass& gate, size_t qubit, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0)
		{
			const size_t gateQubits = gate.getQubitsNumber();

			CheckQubits(gate, qubit, controllingQubit1, controllingQubit2, gateQubits);

#define OPTIMIZED_TENSOR_PRODUCT 1
#ifdef OPTIMIZED_TENSOR_PRODUCT
			assert(gateQubits > 0 && gateQubits <= 3);

			const size_t qubitBit = 1ULL << qubit;

			const MatrixClass& gateMatrix = gate.getRawOperatorMatrix();

			// TODO: perhaps also optimize better for controlled gates
			// TODO: there are ways to optimize further for particular kind of gates, probably I won't bother since it's only a constant factor reduction
			
			bool swapStorage = true;
			if (gateQubits == 1)
			{
				if (NrBasisStates < BaseClass::OneQubitOmpLimit)
					BaseClass::ApplyOneQubitGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, NrBasisStates, swapStorage);
				else
					BaseClass::ApplyOneQubitGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, NrBasisStates, swapStorage);
			}
			else if (gateQubits == 2)
			{
				const size_t ctrlQubitBit = 1ULL << controllingQubit1;

				if (NrBasisStates < BaseClass::TwoQubitOmpLimit)
					BaseClass::ApplyTwoQubitsGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
				else
					BaseClass::ApplyTwoQubitsGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
			}
			else
			{
				const size_t qubitBit2 = 1ULL << controllingQubit1;
				const size_t ctrlQubitBit = 1ULL << controllingQubit2;

				if (NrBasisStates < BaseClass::ThreeQubitOmpLimit)
					BaseClass::ApplyThreeQubitsGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
				else
					BaseClass::ApplyThreeQubitsGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
			}

			if (swapStorage) registerStorage.swap(resultsStorage);
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

				for (size_t i = 0; i < NrBasisStates; ++i)
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


		void displayState(size_t state) const
		{
			const size_t nQubits = getNrQubits();
			std::cout << "|";

			size_t mask = 1ULL << (nQubits - 1);
			for (size_t qubit = 0; qubit < nQubits; ++qubit)
			{
				std::cout << ((state & mask) ? "1" : "0");
				mask >>= 1;
			}

			std::cout << ">    ";
		}

		void displayRegister() const
		{
			const size_t nQubits = getNrQubits();
			const size_t nStates = getNrBasisStates();

			//std::cout << std::fixed;
			std::cout << std::setprecision(4);

			for (size_t state = 0; state < nStates; ++state)
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

		double GetQubitProbability(size_t qubit) const
		{
			if (NrBasisStates < BaseClass::OneQubitOmpLimit)
				return BaseClass::GetQubitProbability(NrBasisStates, registerStorage, qubit);

			return BaseClass::GetQubitProbabilityOmp(NrBasisStates, registerStorage, qubit);
		}

	protected:
		inline void CheckQubits(const GateClass& gate, size_t qubit, size_t controllingQubit1, size_t controllingQubit2, size_t gateQubits) const
		{
			if (NrQubits == 0) throw std::invalid_argument("Qubit number is zero");
		}

		// the following ones should be used for 'repeated measurements' that avoid reexecuting the circuit each time
		size_t MeasureNoCollapse()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
			double accum = 0;
			size_t state = NrBasisStates - 1;
			for (size_t i = 0; i < NrBasisStates; ++i)
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
		size_t MeasureNoCollapse(size_t qubit)
		{
			return MeasureNoCollapse(qubit, qubit);
		}

		// measure a 'subregister' as a separate register
		// can measure a single qubit, if firstQubit == secondQubit
		// will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on position 0 and so on)
		// so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other measured ones are zero, 2 means that the next one 1 one and all others are zero and so on

		size_t MeasureNoCollapse(size_t firstQubit, size_t secondQubit)
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 

			if (firstQubit == secondQubit)
			{
				if (NrBasisStates < BaseClass::OneQubitOmpLimit)
					return BaseClass::MeasureQubitNoCollapse(NrBasisStates, registerStorage, firstQubit, prob);

				return BaseClass::MeasureQubitNoCollapseOmp(NrBasisStates, registerStorage, firstQubit, prob);
			}
			
			return BaseClass::MeasureNoCollapse(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}

		size_t NrQubits;
		size_t NrBasisStates;

		VectorClass registerStorage;
		VectorClass resultsStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;

		std::vector<Gates::AppliedGate<MatrixClass>> computeGates;
		bool recordGates;
	};

}



