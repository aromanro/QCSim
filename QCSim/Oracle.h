#pragma once

#include "QuantumAlgorithm.h"
#include "NControlledNotWithAncilla.h"
#include "NControlledGateWithAncilla.h" 

#include <numeric>

namespace QC {
	namespace Oracles {

		template<class Function, class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class OracleWithGatesSimpleFunction :
			public QC::QuantumSubAlgorithmOnSubregisterWithAncilla<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QC::QuantumSubAlgorithmOnSubregisterWithAncilla<VectorClass, MatrixClass>;
			using RegisterClass = typename BaseClass::RegisterClass;

			OracleWithGatesSimpleFunction(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX, size_t targetQubit = 0, size_t startAncilla = INT_MAX)
				: BaseClass(N, startQubit, endQubit), nControlledNOT(INT_MAX), target(targetQubit)
			{
				std::vector<size_t> controlQubits(endQubit - startQubit + 1);
				std::iota(controlQubits.begin(), controlQubits.end(), startQubit);

				nControlledNOT.SetControlQubits(controlQubits);
				nControlledNOT.SetTargetQubit(target);
				nControlledNOT.SetStartAncillaQubits(startAncilla);
			}

			size_t Execute(RegisterClass& reg) override
			{
				size_t stQubit = BaseClass::BaseClass::getStartQubit();
				size_t enQubit = BaseClass::BaseClass::getEndQubit();
				size_t nrQubits = enQubit - stQubit + 1;
				size_t nrBasisStates = 1ULL << nrQubits;

				/*
					for (size_t state = 0; state < nrBasisStates; ++state)
					{
						if (!func(state)) continue;

						size_t realState = state << stQubit;

						size_t v = realState;
						for (size_t q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, stQubit + q);

							v >>= 1;
						}

						nControlledNOT.Execute(reg);

						// undo the x gates
						v = realState;
						for (size_t q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, stQubit + q);

							v >>= 1;
						}
					}
				*/

				// the above code is not that good, I'll let it commented in case this one is too hard to understand
				// this one avoids applying x gates twice on the same qubit, between applying n controlled nots
				std::vector<bool> qubits(nrQubits, false);

				for (size_t state = 0; state < nrBasisStates; ++state)
				{
					if (!func(state)) continue;

					size_t realState = state << stQubit;

					size_t v = realState;
					for (size_t q = 0; q < nrQubits; ++q)
					{
						if (qubits[q] != ((v & 1) == 0))
						{
							reg.ApplyGate(x, stQubit + q);
							qubits[q] = !qubits[q]; // x gate was applied
						}

						v >>= 1;
					}

					nControlledNOT.Execute(reg);
				}

				// undo the x gates
				for (size_t q = 0; q < nrQubits; ++q)
					if (qubits[q])
						reg.ApplyGate(x, stQubit + q);

				return 0;
			}

			void setFunction(const Function& f)
			{
				func = f;
			}

		protected:
			Function func;

			QC::SubAlgo::NControlledNotWithAncilla<VectorClass, MatrixClass> nControlledNOT;
			QC::Gates::PauliXGate<MatrixClass> x;

			size_t target;
		};


		// as above (TODO: maybe use a common base class?), but with a more complex function, not only one that returns a bool

		template<class Function, class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class OracleWithGates :
			public QC::QuantumSubAlgorithmOnSubregisterWithAncilla<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QC::QuantumSubAlgorithmOnSubregisterWithAncilla<VectorClass, MatrixClass>;
			using RegisterClass = typename BaseClass::RegisterClass;

			OracleWithGates(size_t N, size_t startQubit = 0, size_t endQubit = INT_MAX, size_t startTargetQubits = 0, size_t startAncilla = INT_MAX)
				: BaseClass(N, startQubit, endQubit), nControlledNOTs(INT_MAX), target(startTargetQubits)
			{
				std::vector<size_t> controlQubits(endQubit - startQubit + 1);
				std::iota(controlQubits.begin(), controlQubits.end(), startQubit);

				nControlledNOTs.SetControlQubits(controlQubits);
				nControlledNOTs.SetStartAncillaQubits(startAncilla);
			}

			size_t Execute(RegisterClass& reg) override
			{
				size_t stQubit = BaseClass::BaseClass::getStartQubit();
				size_t enQubit = BaseClass::BaseClass::getEndQubit();
				size_t nrQubits = enQubit - stQubit + 1;
				size_t nrBasisStates = 1ULL << nrQubits;

				/*
					for (size_t state = 0; state < nrBasisStates; ++state)
					{
						size_t fval = func(state);
						if (!fval) continue;
						nControlledNOTs.ClearGates();
						
						size_t q = target;
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

						size_t realState = state << stQubit;

						size_t v = realState;
						for (size_t q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, stQubit + q);

							v >>= 1;
						}

						nControlledNOT.Execute(reg);

						// undo the x gates
						v = realState;
						for (size_t q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, stQubit + q);

							v >>= 1;
						}
					}
				*/

				// the above code is not that good, I'll let it commented in case this one is too hard to understand
				// this one avoids applying x gates twice on the same qubit, between applying n controlled nots
				std::vector<bool> qubits(nrQubits, false);

				for (size_t state = 0; state < nrBasisStates; ++state)
				{
					size_t fval = func(state);
					if (!fval) continue;

					nControlledNOTs.ClearGates();
					
					for (size_t q = target; fval; ++q)
					{
						if (fval & 1)
						{
							QC::Gates::AppliedGate<MatrixClass> notGate(x.getRawOperatorMatrix(), q);
							nControlledNOTs.AddGate(notGate);
						}

						fval >>= 1;
					}

					size_t realState = state << stQubit;

					size_t v = realState;
					for (size_t q = 0; q < nrQubits; ++q)
					{
						if (qubits[q] != ((v & 1) == 0))
						{
							reg.ApplyGate(x, stQubit + q);
							qubits[q] = !qubits[q]; // x gate was applied
						}

						v >>= 1;
					}

					nControlledNOTs.Execute(reg);
				}

				// undo the x gates
				for (size_t q = 0; q < nrQubits; ++q)
					if (qubits[q])
						reg.ApplyGate(x, stQubit + q);

				return 0;
			}

			void setFunction(const Function& f)
			{
				func = f;
			}

		private:
			Function func;

			QC::SubAlgo::NControlledGatesWithAncilla<VectorClass, MatrixClass> nControlledNOTs;
			QC::Gates::PauliXGate<MatrixClass> x;

			size_t target;
		};
	}
}
