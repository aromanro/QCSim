#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"
#include "NControlledGateWithAncilla.h"

#include "Tests.h"

#include <unordered_set>
#include <numeric>


namespace Simon {

	template<class MatrixClass = Eigen::MatrixXcd> class Oracle :
		public QC::Gates::QuantumGate<MatrixClass>
	{
	public:
		void setString(unsigned int str, unsigned int N)
		{
			stringFunction = str;

			const unsigned int nrBasisStates = 1u << N;
			functionTable.resize(nrBasisStates);

			// generate a function compatible with the string
			for (unsigned int i = 0; i < nrBasisStates; ++i)
				functionTable[i] = i;

			// this is a random way of generating it
			const unsigned long long int seed = std::chrono::steady_clock::now().time_since_epoch().count();
			auto rng = std::default_random_engine(static_cast<unsigned int>(seed));
			std::shuffle(functionTable.begin(), functionTable.end(), rng);

			// add some more randomness, do not look in order to indices when searching for a matching index
			// this is probably unneeded complexity, but...
			std::vector<unsigned int> indicesTable(functionTable.size());
			for (unsigned int i = 0; i < nrBasisStates; ++i)
				indicesTable[i] = i;
			std::shuffle(indicesTable.begin(), indicesTable.end(), rng);
		
			for (unsigned int x = 0; x < nrBasisStates; ++x)
				for (unsigned int i = 0; i < nrBasisStates; ++i)
				{
					const unsigned int y = indicesTable[i];
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

		MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const override
		{
			const unsigned int nrBasisStates = 1u << nrQubits;
			const unsigned int N = nrQubits >> 1;

			MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

			// we have B |x>|y> = |x> |f(x) + y>
			// B is unitary

			// <y1|<x1| B |x2>|y2> = <x1|x2><y1|y2+f(x2)>
			// + is XOR here

			const unsigned int xmask = (1u << N) - 1;
			const unsigned int ymask = xmask << N;

			for (unsigned int stateBra = 0; stateBra < nrBasisStates; ++stateBra)
				for (unsigned int stateKet = 0; stateKet < nrBasisStates; ++stateKet)
				{
					const unsigned int xval = (stateBra & xmask);
					const unsigned int ypart = (stateBra & ymask);
					//const unsigned int yval = ypart >> N;

					// the operator is sumi sumj |i><j|, so line corresponds to ket and column to bra 
					extOperatorMat(stateKet, stateBra) = (stateKet & ymask) == ((f(xval) << N) ^ ypart) && xval == (stateKet & xmask);
				}

			assert(checkUnitary(extOperatorMat));

			return extOperatorMat;
		}

		bool checkFunction() const
		{
			if (stringFunction == 0)
			{
				std::unordered_set<unsigned int> s(functionTable.begin(), functionTable.end());
				if (s.size() != functionTable.size()) return false;
			}
			else
			{
				std::unordered_map<unsigned int, std::unordered_set<unsigned int>> m;
				for (unsigned int i = 0; i < functionTable.size(); ++i)
					m[f(i)].insert(i);

				if (m.size() != functionTable.size() / 2) return false;

				for (const auto& s : m)
				{
					if (s.second.size() != 2) return false;

					const unsigned int a = *s.second.begin();
					const unsigned int b = *(++s.second.begin());
					const unsigned int ab = a ^ b;
					if (ab != 0 && ab != stringFunction) return false;
				}
			}

			return true;
		}

		unsigned int f(unsigned int x) const
		{
			if (functionTable.empty()) return 0;

			const unsigned int mask = static_cast<unsigned int>(functionTable.size() - 1);

			return functionTable[x & mask];
		}

	protected:
		unsigned int stringFunction = 0;
		std::vector<unsigned int> functionTable;
	};

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class SimonAlgorithm :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		SimonAlgorithm(unsigned int N = 3, int addseed = 0)
			: BaseClass(2 * N, addseed)
		{
			assert(N >= 2);

			setString(0); // prevent issues if the string is not set before execution
		}

		void setString(unsigned int str)
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();
			const int unsigned N = nrQubits >> 1;

			oracle.setString(str, N);
			OracleOp = oracle.getOperatorMatrix(nrQubits);
		}

		unsigned int Execute() override
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();
			const int unsigned N = nrQubits >> 1;
			const unsigned int mask = (1u << N) - 1;
			const unsigned int nrBasisStates = 1u << N;

			std::unordered_map<unsigned int, unsigned int> measurements;
			const unsigned int nrMeasurements = 300; // make it highly unlikely to fail

			// not exactly how it should be done, but again, I'm lazy

			for (unsigned int i = 0;i < nrMeasurements; ++i)
			{
				Init();
				BaseClass::ApplyOperatorMatrix(OracleOp);
				ApplyHadamardOnHalfQubits();

				// the measurement is always a value that has b * s = 0 property, where s is the string function
				// unless the string function is zero, in which case all of them could be measured, they have equal probability
				// either measure only the first half of the qubits
				 
				const unsigned int m = BaseClass::Measure(0, N - 1);
				
				// or all of them but discard the not interesting ones, using a mask:
				//const unsigned int m = (BaseClass::Measure() & mask);

				++measurements[m];

				if (measurements.size() == nrBasisStates) break;
			}

			if (measurements.size() == nrBasisStates) return 0;

			std::unordered_set<unsigned int> potential_results;
			for (unsigned int i = 1; i <= mask; ++i)
				potential_results.insert(i);
			
			// Not exactly the most optimal way, but I'm lazy
			for (auto it = measurements.begin(); it != measurements.end(); ++it)
			{
				if (it->first == 0) continue;

				for (auto pit = potential_results.begin(); pit != potential_results.end();)
				{
					unsigned int v = ((*pit) & it->first);
					unsigned int cnt = 0;
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
			const unsigned int halfQubits = BaseClass::getNrQubits() >> 1;
			for (unsigned int i = 0; i < halfQubits; ++i)
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

		SimonAlgorithmWithGatesOracle(unsigned int N = 3, int addseed = 0)
			: BaseClass(3 * N - 1, addseed),
			nControlledNOTs(INT_MAX)
		{
			assert(N >= 2);

			setString(0); // prevent issues if the string is not set before execution

			std::vector<unsigned int> controlQubits(N);
			std::iota(controlQubits.begin(), controlQubits.end(), 0);

			nControlledNOTs.SetControlQubits(controlQubits);
			nControlledNOTs.SetStartAncillaQubits(2 * N);
		}

		void setString(unsigned int str)
		{
			const unsigned int nrQubits = getAlgoQubits();
			const int unsigned N = nrQubits >> 1;

			oracle.setString(str, N);
		}

		unsigned int Execute() override
		{
			const unsigned int nrQubits = getAlgoQubits();
			const int unsigned N = nrQubits >> 1;
			const unsigned int mask = (1u << N) - 1;
			const unsigned int nrBasisStates = 1u << N;

			std::unordered_map<unsigned int, unsigned int> measurements;
			const unsigned int nrMeasurements = 300; // make it highly unlikely to fail

			// not exactly how it should be done, but again, I'm lazy

			for (unsigned int i = 0; i < nrMeasurements; ++i)
			{
				Init();
				
				ApplyOracle();

				ApplyHadamardOnHalfQubits();

				// the measurement is always a value that has b * s = 0 property, where s is the string function
				// unless the string function is zero, in which case all of them could be measured, they have equal probability
				// either measure only the first half of the qubits

				const unsigned int m = BaseClass::Measure(0, N - 1);

				// or all of them but discard the not interesting ones, using a mask:
				//const unsigned int m = (BaseClass::Measure() & mask);

				++measurements[m];

				if (measurements.size() == nrBasisStates) break;
			}

			if (measurements.size() == nrBasisStates) return 0;

			std::unordered_set<unsigned int> potential_results;
			for (unsigned int i = 1; i <= mask; ++i)
				potential_results.insert(i);

			// Not exactly the most optimal way, but I'm lazy
			for (auto it = measurements.begin(); it != measurements.end(); ++it)
			{
				if (it->first == 0) continue;

				for (auto pit = potential_results.begin(); pit != potential_results.end();)
				{
					unsigned int v = ((*pit) & it->first);
					unsigned int cnt = 0;
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

		unsigned int getAlgoQubits() const
		{
			const unsigned int nrQubits = BaseClass::getNrQubits();

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
			const unsigned int halfQubits = getAlgoQubits() >> 1;
			for (unsigned int i = 0; i < halfQubits; ++i)
				BaseClass::ApplyGate(hadamard, i);
		}

		void ApplyOracle()
		{
			const unsigned int nrQubits = getAlgoQubits() >> 1;
			const unsigned int nrBasisStates = 1u << nrQubits;

			/*
			for (unsigned int state = 0; state < nrBasisStates; ++state)
			{
				unsigned int fval = oracle.f(state);
				if (!fval) continue;

				nControlledNOTs.ClearGates();
				unsigned int q = nrQubits;
				while (fval)
				{
					if (fval & 1)
					{
						QC::Gates::AppliedGate<MatrixClass> notGate(x.getRawOperatorMatrix(), q);
						nControlledNOTs.AddGate(notGate);
					}

					fval >>= 1;
					++q;
				}

				unsigned int v = state;
				for (unsigned int q = 0; q < nrQubits; ++q)
				{
					if ((v & 1) == 0)
						BaseClass::ApplyGate(x, q);

					v >>= 1;
				}

				nControlledNOTs.Execute(BaseClass::reg);

				// undo the x gates
				v = state;
				for (unsigned int q = 0; q < nrQubits; ++q)
				{
					if ((v & 1) == 0)
						BaseClass::ApplyGate(x, q);

					v >>= 1;
				}
			}
			*/

			// the above code is not that good, I'll let it commented in case this one is too hard to understand
			// this one avoids applying x gates twice on the same qubit, between applying n controlled nots

			// the reason for not applying NControlledNOTWithAncilla is that it computes/uncomputes for each gate
			// the more general NControlledGatesWithAncilla does a computation for the control ancilla only once, applies all added gates (controlled) then it does the uncomputation

			std::vector<bool> qubits(nrQubits, false);

			for (unsigned int state = 0; state < nrBasisStates; ++state)
			{
				unsigned int fval = oracle.f(state);
				if (!fval) continue;

				nControlledNOTs.ClearGates();
				unsigned int q = nrQubits;
				while (fval)
				{
					if (fval & 1)
					{
						QC::Gates::AppliedGate<MatrixClass> notGate(x.getRawOperatorMatrix(), q);
						nControlledNOTs.AddGate(notGate);
					}

					fval >>= 1;
					++q;
				}

				unsigned int v = state;
				for (unsigned int q = 0; q < nrQubits; ++q)
				{
					if (qubits[q] != ((v & 1) == 0))
					{
						BaseClass::ApplyGate(x, q);
						qubits[q] = !qubits[q]; // x gate was applied
					}

					v >>= 1;
				}

				nControlledNOTs.Execute(BaseClass::reg);
			}

			// undo the x gates
			for (unsigned int q = 0; q < nrQubits; ++q)
				if (qubits[q])
					BaseClass::ApplyGate(x, q);
		}

		Oracle<MatrixClass> oracle;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::PauliXGate<MatrixClass> x;
		QC::SubAlgo::NControlledGatesWithAncilla<VectorClass, MatrixClass> nControlledNOTs;
	};
}


