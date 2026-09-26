#pragma once

#include "QubitRegisterCalculator.h"

#include <unordered_map>
#include <stdexcept>

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


			if (addseed == 0)
			{
				std::random_device rdl;
				addseed = rdl();
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

		void SetSeed(uint64_t theSeed)
		{
			std::seed_seq seed{ uint32_t(theSeed & 0xffffffff), uint32_t(theSeed >> 32) };
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
			const size_t state = BaseClass::SampleBasisState(NrBasisStates, registerStorage, prob, NrBasisStates - 1, UseMultithreading());

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
				if (!UseMultithreading())
					return BaseClass::MeasureQubit(NrBasisStates, registerStorage, firstQubit, prob);

				return BaseClass::MeasureQubitOmp(NrBasisStates, registerStorage, firstQubit, prob);
			}

			if (!UseMultithreading())
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
			for (size_t i = 0; i < (size_t)registerStorage.size(); ++i)
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
			const auto& matrix = gate.getRawOperatorMatrix();
			const auto dimension = matrix.rows();
			if (dimension != matrix.cols() || (dimension != 2 && dimension != 4 && dimension != 8))
				throw std::invalid_argument("Statevector gates must have a 2x2, 4x4 or 8x8 operator matrix");
			const size_t gateQubits = gate.getQubitsNumber();
			const size_t matrixQubits = dimension == 2 ? 1 : (dimension == 4 ? 2 : 3);
			if (gateQubits != matrixQubits)
				throw std::invalid_argument("Gate arity does not match its operator matrix");

			CheckQubits(gate, qubit, controllingQubit1, controllingQubit2, gateQubits);

#define OPTIMIZED_TENSOR_PRODUCT 1
#ifdef OPTIMIZED_TENSOR_PRODUCT
			assert(gateQubits > 0 && gateQubits <= 3);

			const std::array<size_t, 3> bits{
				size_t{1} << qubit,
				gateQubits > 1 ? size_t{1} << controllingQubit1 : 0,
				gateQubits > 2 ? size_t{1} << controllingQubit2 : 0
			};
			BaseClass::ApplyGateInPlace(registerStorage, matrix, gate.getStructure(),
				bits, static_cast<unsigned>(gateQubits), NrBasisStates,
				UseMultithreading());
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
			if (!UseMultithreading())
				return BaseClass::GetQubitProbability(NrBasisStates, registerStorage, qubit);

			return BaseClass::GetQubitProbabilityOmp(NrBasisStates, registerStorage, qubit);
		}

		void SaveState()
		{
			savedStateStorage = registerStorage;
		}

		void RestoreState()
		{
			if (savedStateStorage.size() == 0) return;
			registerStorage = savedStateStorage;
		}

		void RestoreStateDestructive()
		{
			if (savedStateStorage.size() == 0) return;
			registerStorage.swap(savedStateStorage);
			savedStateStorage.resize(0);
		}

		// the following ones should be used for 'repeated measurements' that avoid reexecuting the circuit each time
		size_t MeasureNoCollapse()
		{
			const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy

			return BaseClass::SampleBasisState(NrBasisStates, registerStorage, prob, 0, UseMultithreading());
		}

		// does not check the gates, that's why it returns a complex number
		// the caller should ensure the hermicity and extract the real part
		std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>>& gates)
		{
			if (gates.empty()) return 1.;

			// Pauli strings (the common case) are computed in a single read-only pass, without copying the state
			size_t xMask = 0;
			size_t zMask = 0;
			size_t nrY = 0;
			if (GetPauliStringMasks(gates, xMask, zMask, nrY))
				return BaseClass::PauliExpectationValue(NrBasisStates, registerStorage, xMask, zMask, nrY,
					UseMultithreading());

			VectorClass savedState = registerStorage;

			ApplyGates(gates);

			const auto res = (savedState.adjoint() * registerStorage)(0);

			registerStorage.swap(savedState); // restore the state

			return res;
		}

		std::unique_ptr<QubitRegister<VectorClass, MatrixClass>> Clone() const
		{
			auto qr = std::make_unique<QubitRegister<VectorClass, MatrixClass>>(1);
			qr->SetMultithreading(BaseClass::GetMultithreading());
			qr->NrQubits = NrQubits;
			qr->NrBasisStates = NrBasisStates;
			qr->registerStorage = registerStorage;
			qr->savedStateStorage = savedStateStorage;
			qr->computeGates = computeGates;
			qr->recordGates = recordGates;
			
			return qr;
		}

		MatrixClass getDensityMatrix() const
		{
			return registerStorage * registerStorage.adjoint();
		}

	protected:
		bool UseMultithreading() const
		{
			return BaseClass::GetMultithreading() && NrBasisStates >= BaseClass::GetParallelMinBasisStates();
		}

		inline void CheckQubits(const GateClass& /*gate*/, size_t qubit, size_t controllingQubit1, size_t controllingQubit2, size_t gateQubits) const
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



		// Recognizes a product of single qubit X, Y, Z (or identity) gates, each on a different qubit.
		// Exact comparisons only, anything else (including a repeated qubit) is left to the generic path.
		bool GetPauliStringMasks(const std::vector<Gates::AppliedGate<MatrixClass>>& gates, size_t& xMask, size_t& zMask, size_t& nrY) const
		{
			const std::complex<double> zero(0., 0.);
			const std::complex<double> one(1., 0.);
			const std::complex<double> i(0., 1.);

			size_t usedQubits = 0;
			for (const auto& gate : gates)
			{
				const MatrixClass& m = gate.getRawOperatorMatrix();
				if (m.rows() != 2 || m.cols() != 2) return false;

				const size_t qubit = gate.getQubit1();
				if (qubit >= NrQubits) return false;

				const size_t qubitBit = 1ULL << qubit;
				if (usedQubits & qubitBit) return false;
				usedQubits |= qubitBit;

				if (m(0, 1) == zero && m(1, 0) == zero && m(0, 0) == one)
				{
					if (m(1, 1) == one) continue; // identity
					if (m(1, 1) == -one) // Z
					{
						zMask |= qubitBit;
						continue;
					}
				}
				else if (m(0, 0) == zero && m(1, 1) == zero)
				{
					if (m(0, 1) == one && m(1, 0) == one) // X
					{
						xMask |= qubitBit;
						continue;
					}
					if (m(0, 1) == -i && m(1, 0) == i) // Y
					{
						xMask |= qubitBit;
						zMask |= qubitBit;
						++nrY;
						continue;
					}
				}

				return false;
			}

			return true;
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

		VectorClass savedStateStorage;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;

		std::vector<Gates::AppliedGate<MatrixClass>> computeGates;
		bool recordGates;
	};

}


