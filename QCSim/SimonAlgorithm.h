#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"
#include "Oracle.h"

#include "Tests.h"

#include <unordered_set>


namespace Simon {

	class SimonFunction
	{
	public:
		size_t operator()(size_t x) const
		{
			if (functionTable.empty()) return 0;

			const auto mask = static_cast<size_t>(functionTable.size() - 1);

			return functionTable[x & mask];
		}

		void setString(size_t str, size_t N)
		{
			stringFunction = str;

			const size_t nrBasisStates = 1ULL << N;
			functionTable.resize(nrBasisStates);

			// generate a function compatible with the string
			for (size_t i = 0; i < nrBasisStates; ++i)
				functionTable[i] = i;

			// this is a random way of generating it
			const long long int seed = std::chrono::steady_clock::now().time_since_epoch().count();
			auto rng = std::default_random_engine(static_cast<size_t>(seed));
			std::shuffle(functionTable.begin(), functionTable.end(), rng);

			// add some more randomness, do not look in order to indices when searching for a matching index
			// this is probably unneeded complexity, but...
			std::vector<size_t> indicesTable(functionTable.size());
			for (size_t i = 0; i < nrBasisStates; ++i)
				indicesTable[i] = i;
			std::shuffle(indicesTable.begin(), indicesTable.end(), rng);

			for (size_t x = 0; x < nrBasisStates; ++x)
				for (size_t i = 0; i < nrBasisStates; ++i)
				{
					const size_t y = indicesTable[i];
					if (y == UINT_MAX || x == y) continue;
					// note that is the string function is zero, then the function will be one-to-one
					// because there is always at least a bit that is different between to different values,
					// so the assignment below never takes place
					// if the string is not zero, the function is two-to-one
					if ((x ^ y) == stringFunction)
					{
						functionTable[x] = functionTable[y];
						indicesTable[i] = UINT_MAX;
						break;
					}
				}
		}

		bool checkFunction() const
		{
			if (stringFunction == 0)
			{
				std::unordered_set<size_t> s(functionTable.begin(), functionTable.end());
				if (s.size() != functionTable.size()) return false;
			}
			else
			{
				std::unordered_map<size_t, std::unordered_set<size_t>> m;
				for (size_t i = 0; i < functionTable.size(); ++i)
					m[operator()(i)].insert(i);

				if (m.size() != functionTable.size() / 2) return false;

				for (const auto& s : m)
				{
					if (s.second.size() != 2) return false;

					const size_t a = *s.second.begin();
					const size_t b = *(++s.second.begin());
					const size_t ab = a ^ b;
					if (ab != 0 && ab != stringFunction) return false;
				}
			}

			return true;
		}

	protected:
		size_t stringFunction = 0;
		std::vector<size_t> functionTable;
	};

	template<class MatrixClass = Eigen::MatrixXcd> class Oracle :
		public QC::Gates::QuantumGate<MatrixClass>
	{
	public:
		MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
		{
			const size_t nrBasisStates = 1ULL << nrQubits;
			const size_t N = nrQubits >> 1;

			MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			// we have B |x>|y> = |x> |f(x) + y>
			// B is unitary

			// <y1|<x1| B |x2>|y2> = <x1|x2><y1|y2+f(x2)>
			// + is XOR here

			const size_t xmask = (1u << N) - 1;
			const size_t ymask = xmask << N;

			for (size_t stateBra = 0; stateBra < nrBasisStates; ++stateBra)
				for (size_t stateKet = 0; stateKet < nrBasisStates; ++stateKet)
				{
					const size_t xval = (stateBra & xmask);
					const size_t ypart = (stateBra & ymask);
					//const size_t yval = ypart >> N;

					// the operator is sumi sumj |i><j|, so line corresponds to ket and column to bra 
					extOperatorMat(stateKet, stateBra) = (stateKet & ymask) == ((f(xval) << N) ^ ypart) && xval == (stateKet & xmask);
				}

			assert(checkUnitary(extOperatorMat));

			return extOperatorMat;
		}

		void setString(size_t str, size_t N)
		{
			f.setString(str, N);
		}

	protected:
		SimonFunction f;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SimonAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		SimonAlgorithm(size_t N = 3, int addseed = 0)
			: BaseClass(2 * N, addseed)
		{
			assert(N >= 2);

			setString(0); // prevent issues if the string is not set before execution
		}

		void setString(size_t str)
		{
			const size_t nrQubits = BaseClass::getNrQubits();
			const size_t N = nrQubits >> 1;

			oracle.setString(str, N);
			OracleOp = oracle.getOperatorMatrix(nrQubits);
		}

		size_t Execute() override
		{
			const size_t nrQubits = BaseClass::getNrQubits();
			const size_t N = nrQubits >> 1;
			const size_t mask = (1ULL << N) - 1;
			const size_t nrBasisStates = 1u << N;

			std::unordered_map<size_t, size_t> measurements;
			const size_t nrMeasurements = 300; // make it highly unlikely to fail

			// not exactly how it should be done, but again, I'm lazy

			for (size_t i = 0;i < nrMeasurements; ++i)
			{
				Init();
				BaseClass::ApplyOperatorMatrix(OracleOp);
				ApplyHadamardOnHalfQubits();

				// the measurement is always a value that has b * s = 0 property, where s is the string function
				// unless the string function is zero, in which case all of them could be measured, they have equal probability
				// either measure only the first half of the qubits
				 
				const size_t m = BaseClass::Measure(0, N - 1);
				
				// or all of them but discard the not interesting ones, using a mask:
				//const size_t m = (BaseClass::Measure() & mask);

				++measurements[m];

				if (measurements.size() == nrBasisStates) break;
			}

			if (measurements.size() == nrBasisStates) return 0;

			std::unordered_set<size_t> potential_results;
			for (size_t i = 1; i <= mask; ++i)
				potential_results.insert(i);
			
			// Not exactly the most optimal way, but I'm lazy
			for (auto it = measurements.begin(); it != measurements.end(); ++it)
			{
				if (it->first == 0) continue;

				for (auto pit = potential_results.begin(); pit != potential_results.end();)
				{
					size_t v = ((*pit) & it->first);
					size_t cnt = 0;
					while (v) {
						cnt += v & 1;
						v >>= 1;
					}
					if (cnt % 2) pit = potential_results.erase(pit);
					else ++pit;
				}
			}


			return (potential_results.empty() || potential_results.size() > 1) ? mask + 1 : *potential_results.begin();
		}

	protected:
		void Init()
		{
			BaseClass::setToBasisState(0);
			ApplyHadamardOnHalfQubits();
		}

		void ApplyHadamardOnHalfQubits()
		{
			const size_t halfQubits = BaseClass::getNrQubits() >> 1;
			for (size_t i = 0; i < halfQubits; ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		Oracle<MatrixClass> oracle;
		MatrixClass OracleOp;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
	};



	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SimonAlgorithmWithGatesOracle :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		SimonAlgorithmWithGatesOracle(size_t N = 3, int addseed = 0)
			: BaseClass(3 * N - 1, addseed), // 2 * N qubits for the algorithm, N - 1 for the oracle
			oracle(3 * N - 1, 0, N - 1, N, 2 * N)
		{
			assert(N >= 2);

			setString(0); // prevent issues if the string is not set before execution
		}

		void setString(size_t str)
		{
			const size_t nrQubits = getAlgoQubits();
			const size_t N = nrQubits >> 1;

			func.setString(str, N);
			oracle.setFunction(func);
		}

		size_t Execute() override
		{
			const size_t nrQubits = getAlgoQubits();
			const size_t N = nrQubits >> 1;
			const size_t mask = (1ULL << N) - 1;
			const size_t nrBasisStates = 1ULL << N;

			std::unordered_map<size_t, size_t> measurements;
			const size_t nrMeasurements = 300; // make it highly unlikely to fail

			// not exactly how it should be done, but again, I'm lazy

			for (size_t i = 0; i < nrMeasurements; ++i)
			{
				Init();
				
				ApplyOracle();

				ApplyHadamardOnHalfQubits();

				// the measurement is always a value that has b * s = 0 property, where s is the string function
				// unless the string function is zero, in which case all of them could be measured, they have equal probability
				// either measure only the first half of the qubits

				const size_t m = BaseClass::Measure(0, N - 1);

				// or all of them but discard the not interesting ones, using a mask:
				//const size_t m = (BaseClass::Measure() & mask);

				++measurements[m];

				if (measurements.size() == nrBasisStates) break;
			}

			if (measurements.size() == nrBasisStates) return 0;

			std::unordered_set<size_t> potential_results;
			for (size_t i = 1; i <= mask; ++i)
				potential_results.insert(i);

			// Not exactly the most optimal way, but I'm lazy
			for (auto it = measurements.cbegin(); it != measurements.cend(); ++it)
			{
				if (it->first == 0) continue;

				for (auto pit = potential_results.begin(); pit != potential_results.end();)
				{
					size_t v = ((*pit) & it->first);
					size_t cnt = 0;
					while (v) {
						cnt += v & 1;
						v >>= 1;
					}
					if (cnt % 2) pit = potential_results.erase(pit);
					else ++pit;
				}
			}

			return (potential_results.empty() || potential_results.size() > 1) ? mask + 1 : *potential_results.cbegin();
		}

		size_t getAlgoQubits() const
		{
			const size_t nrQubits = BaseClass::getNrQubits();

			return (nrQubits + 1) / 3 * 2;
		}

	protected:
		void Init()
		{
			BaseClass::setToBasisState(0);
			ApplyHadamardOnHalfQubits();
		}

		void ApplyHadamardOnHalfQubits()
		{
			const size_t halfQubits = getAlgoQubits() >> 1;
			for (size_t i = 0; i < halfQubits; ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyOracle()
		{
			oracle.Execute(BaseClass::reg);
		}

		SimonFunction func;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Oracles::OracleWithGates<SimonFunction, VectorClass, MatrixClass> oracle;
	};
}


