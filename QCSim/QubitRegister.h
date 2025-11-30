#pragma once

#include "QubitRegisterCalculator.h"

#include <unordered_map>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegister : public QubitRegisterCalculator<VectorClass, MatrixClass>
	{
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
		using BaseClass = QubitRegisterCalculator<VectorClass, MatrixClass>;

		QubitRegister(size_t N = 3, unsigned int addseed = 0)
			: NrQubits(N), NrBasisStates(1ULL << NrQubits), 
			registerStorage(VectorClass::Zero(NrBasisStates)),
			uniformZeroOne(0, 1), recordGates(false)
		{
			assert(N > 0);

			resultsStorage.resize(NrBasisStates);

			if (addseed == 0)
			{
				std::random_device rd;
				addseed = rd();
			}

			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);

			registerStorage(0) = 1;
		}

		// this is a special constructor, I need it for something in a derived work (closed source)
		QubitRegister(size_t N, VectorClass& v, unsigned int addseed = 0)
			: NrQubits(N), NrBasisStates(1ULL << NrQubits),
			uniformZeroOne(0, 1), recordGates(false)
		{
			assert(N > 0);
			resultsStorage.resize(NrBasisStates);
			registerStorage.swap(v);

			if (addseed == 0)
			{
				std::random_device rd;
				addseed = rd();
			}

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

		void Reset()
		{
			setToBasisState(0);
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
				++i;
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
				if (!BaseClass::GetMultithreading() || NrBasisStates < BaseClass::OneQubitOmpLimit)
					return BaseClass::MeasureQubit(NrBasisStates, registerStorage, firstQubit, prob);

				return BaseClass::MeasureQubitOmp(NrBasisStates, registerStorage, firstQubit, prob);
			}

			if (!BaseClass::GetMultithreading() || NrBasisStates < BaseClass::OneQubitOmpLimit)
				return BaseClass::Measure(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
			
			return  BaseClass::MeasureOmp(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}


		std::map<size_t, size_t> RepeatedMeasure(size_t nrTimes = 1000)
		{
			if (nrTimes == 0) return {}; // nothing to measure

			std::map<size_t, size_t> measurements;

			if (nrTimes == 1) // shortcut for a single measurement
			{
				const size_t meas = MeasureNoCollapse();
				++measurements[meas];
				return measurements;
			}

			// a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N), where N is the number of states), 
			// there is an even O(1) method (see https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

			std::vector<double> probabilities(registerStorage.size());

			double accum = 0;
			for (size_t i = 0; i < static_cast<size_t>(registerStorage.size()); ++i)
			{
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
				++i;
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
			}

			for (size_t shot = 0; shot < nrTimes; ++shot)
			{
				const double prob = 1. - uniformZeroOne(rng);
				const size_t meas = std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

				++measurements[meas];
			}

			return measurements;
		}

		std::unordered_map<size_t, size_t> RepeatedMeasureUnordered(size_t nrTimes = 1000)
		{
			if (nrTimes == 0) return {}; // nothing to measure

			std::unordered_map<size_t, size_t> measurements;

			if (nrTimes == 1) // shortcut for a single measurement
			{
				const size_t meas = MeasureNoCollapse();
				++measurements[meas];
				return measurements;
			}

			// a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N), where N is the number of states), 
			// there is an even O(1) method (see https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

			std::vector<double> probabilities(registerStorage.size());

			double accum = 0;
			for (size_t i = 0; i < registerStorage.size(); ++i)
			{
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
				++i;
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
			}

			for (size_t shot = 0; shot < nrTimes; ++shot)
			{
				const double prob = 1. - uniformZeroOne(rng);
				const size_t meas = std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

				++measurements[meas];
			}

			return measurements;
		}

		std::map<size_t, size_t> RepeatedMeasure(size_t firstQubit, size_t secondQubit, size_t nrTimes = 1000)
		{
			if (nrTimes == 0) return {}; // nothing to measure

			std::map<size_t, size_t> measurements;
			
			const size_t secondQubitp1 = secondQubit + 1;
			const size_t firstPartMask = (1ULL << firstQubit) - 1;
			const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;
			
			if (nrTimes == 1) // shortcut for a single measurement
			{
				const size_t meas = MeasureNoCollapse(firstQubit, secondQubit);
				++measurements[(meas & measuredPartMask) >> firstQubit];
				return measurements;
			}

			// a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N), where N is the number of states), 
			// there is an even O(1) method (see https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

			std::vector<double> probabilities(registerStorage.size());

			double accum = 0;
			for (size_t i = 0; i < static_cast<size_t>(registerStorage.size()); ++i)
			{
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
				++i;
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
			}

			for (size_t shot = 0; shot < nrTimes; ++shot)
			{
				const double prob = 1. - uniformZeroOne(rng);
				const size_t meas = std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

				++measurements[(meas & measuredPartMask) >> firstQubit];
			}

			return measurements;
		}

		std::unordered_map<size_t, size_t> RepeatedMeasureUnordered(size_t firstQubit, size_t secondQubit, size_t nrTimes = 1000)
		{
			if (nrTimes == 0) return {}; // nothing to measure

			std::unordered_map<size_t, size_t> measurements;

			const size_t secondQubitp1 = secondQubit + 1;
			const size_t firstPartMask = (1ULL << firstQubit) - 1;
			const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;

			if (nrTimes == 1) // shortcut for a single measurement
			{
				const size_t meas = MeasureNoCollapse(firstQubit, secondQubit);
				++measurements[(meas & measuredPartMask) >> firstQubit];
				return measurements;
			}

			// a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N), where N is the number of states), 
			// there is an even O(1) method (see https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

			std::vector<double> probabilities(registerStorage.size());

			double accum = 0;
			for (size_t i = 0; i < registerStorage.size(); ++i)
			{
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
				++i;
				accum += std::norm(registerStorage[i]);
				probabilities[i] = accum;
				if (accum > 1.0 - std::numeric_limits<double>::epsilon())
				{
					probabilities.resize(i + 1);
					break;
				}
			}

			for (size_t shot = 0; shot < nrTimes; ++shot)
			{
				const double prob = 1. - uniformZeroOne(rng);
				const size_t meas = std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

				++measurements[(meas & measuredPartMask) >> firstQubit];
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
			// TODO: there are ways to optimize further for particular kind of gates, probably I won't bother
			
			bool swapStorage = true;
			if (gateQubits == 1)
			{
				if (!BaseClass::GetMultithreading() || NrBasisStates < BaseClass::OneQubitOmpLimit)
					BaseClass::ApplyOneQubitGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, NrBasisStates, swapStorage);
				else
					BaseClass::ApplyOneQubitGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, NrBasisStates, swapStorage);
			}
			else if (gateQubits == 2)
			{
				const size_t ctrlQubitBit = 1ULL << controllingQubit1;

				if (!BaseClass::GetMultithreading() || NrBasisStates < BaseClass::TwoQubitOmpLimit)
					BaseClass::ApplyTwoQubitsGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
				else
					BaseClass::ApplyTwoQubitsGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
			}
			else
			{
				const size_t qubitBit2 = 1ULL << controllingQubit1;
				const size_t ctrlQubitBit = 1ULL << controllingQubit2;

				if (!BaseClass::GetMultithreading() || NrBasisStates < BaseClass::ThreeQubitOmpLimit)
					BaseClass::ApplyThreeQubitsGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				else
					BaseClass::ApplyThreeQubitsGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
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
			if (!BaseClass::GetMultithreading() || NrBasisStates < BaseClass::OneQubitOmpLimit)
				return BaseClass::GetQubitProbability(NrBasisStates, registerStorage, qubit);

			return BaseClass::GetQubitProbabilityOmp(NrBasisStates, registerStorage, qubit);
		}

		void SaveState()
		{
			savedSatateStorage = registerStorage;
		}

		void RestoreState()
		{
			registerStorage = savedSatateStorage;
		}

		// the following ones should be used for 'repeated measurements' that avoid reexecuting the circuit each time
		size_t MeasureNoCollapse()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy 
			double accum = 0;
			size_t state = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					state = i;
					break;
				}
				++i;
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					state = i;
					break;
				}
			}

			return state;
		}

		// does not check the gates, that's why it returns a complex number
		// the caller should ensure the hermicity and extract the real part
		std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>>& gates)
		{
			if (gates.empty()) return 1.;

			// TODO: there are faster methods for special cases, like Pauli strings!
			VectorClass savedState = registerStorage;

			ApplyGates(gates);

			const auto res = (savedState.adjoint() * registerStorage)(0);

			registerStorage.swap(savedState); // restore the state

			return res;
		}

		std::unique_ptr<QubitRegister<VectorClass, MatrixClass>> Clone() const
		{
			auto qr = std::make_unique<QubitRegister<VectorClass, MatrixClass>>(1);
			qr->NrQubits = NrQubits;
			qr->NrBasisStates = NrBasisStates;
			qr->registerStorage = registerStorage;
			qr->resultsStorage = resultsStorage;
			qr->savedSatateStorage = savedSatateStorage;
			qr->computeGates = computeGates;
			qr->recordGates = recordGates;
			
			return qr;
		}

	protected:
		inline void CheckQubits(const GateClass& gate, size_t qubit, size_t controllingQubit1, size_t controllingQubit2, size_t gateQubits) const
		{
			if (NrQubits == 0) throw std::invalid_argument("Qubit number is zero");
			else if (NrQubits <= qubit) throw std::invalid_argument("Qubit number is too high");
			else if (gateQubits == 2) {
				if (NrQubits <= controllingQubit1) throw std::invalid_argument("Controlling qubit number is too high");
				else if (qubit == controllingQubit1) throw std::invalid_argument("Qubit and controlling qubit are the same");
			} 
			else if (gateQubits == 3)
			{
				if (NrQubits <= controllingQubit1 || NrQubits <= controllingQubit2) throw std::invalid_argument("Controlling qubit number is too high");
				else if (qubit == controllingQubit1 || qubit == controllingQubit2 || controllingQubit1 == controllingQubit2) throw std::invalid_argument("Qubits must be different");
			}
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
				return BaseClass::MeasureQubitNoCollapse(NrBasisStates, registerStorage, firstQubit, prob);
			
			return BaseClass::MeasureNoCollapse(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
		}

		size_t NrQubits;
		size_t NrBasisStates;

		VectorClass registerStorage;
		VectorClass resultsStorage;

		VectorClass savedSatateStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;

		std::vector<Gates::AppliedGate<MatrixClass>> computeGates;
		bool recordGates;
	};

}



