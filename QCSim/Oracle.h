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

			OracleWithGatesSimpleFunction(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX, unsigned int targetQubit = 0, unsigned int startAncilla = INT_MAX)
				: BaseClass(N, startQubit, endQubit), nControlledNOT(INT_MAX), target(targetQubit)
			{
				std::vector<unsigned int> controlQubits(endQubit - startQubit + 1);
				std::iota(controlQubits.begin(), controlQubits.end(), startQubit);

				nControlledNOT.SetControlQubits(controlQubits);
				nControlledNOT.SetTargetQubit(target);
				nControlledNOT.SetStartAncillaQubits(startAncilla);
			}

			unsigned int Execute(RegisterClass& reg) override
			{
				unsigned int sQubit = BaseClass::BaseClass::getStartQubit();
				unsigned int eQubit = BaseClass::BaseClass::getEndQubit();
				unsigned int nrQubits = eQubit - sQubit + 1;
				unsigned int nrBasisStates = 1 << nrQubits;

				/*
					for (unsigned int state = 0; state < nrBasisStates; ++state)
					{
						if (!func(state)) continue;

						unsigned int realState = state << sQubit;

						unsigned int v = realState;
						for (unsigned int q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, sQubit + q);

							v >>= 1;
						}

						nControlledNOT.Execute(reg);

						// undo the x gates
						v = realState;
						for (unsigned int q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, sQubit + q);

							v >>= 1;
						}
					}
				*/

				// the above code is not that good, I'll let it commented in case this one is too hard to understand
				// this one avoids applying x gates twice on the same qubit, between applying n controlled nots
				std::vector<bool> qubits(nrQubits, false);

				for (unsigned int state = 0; state < nrBasisStates; ++state)
				{
					if (!func(state)) continue;

					unsigned int realState = state << sQubit;

					unsigned int v = realState;
					for (unsigned int q = 0; q < nrQubits; ++q)
					{
						if (qubits[q] != ((v & 1) == 0))
						{
							reg.ApplyGate(x, sQubit + q);
							qubits[q] = !qubits[q]; // x gate was applied
						}

						v >>= 1;
					}

					nControlledNOT.Execute(reg);
				}

				// undo the x gates
				for (unsigned int q = 0; q < nrQubits; ++q)
					if (qubits[q])
						reg.ApplyGate(x, sQubit + q);

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

			unsigned int target;
		};


		// as above (TODO: maybe use a common base class?), but with a more complex function, not only one that returns a bool

		template<class Function, class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class OracleWithGates :
			public QC::QuantumSubAlgorithmOnSubregisterWithAncilla<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QC::QuantumSubAlgorithmOnSubregisterWithAncilla<VectorClass, MatrixClass>;
			using RegisterClass = typename BaseClass::RegisterClass;

			OracleWithGates(unsigned int N, unsigned int startQubit = 0, unsigned int endQubit = INT_MAX, unsigned int startTargetQubits = 0, unsigned int startAncilla = INT_MAX)
				: BaseClass(N, startQubit, endQubit), nControlledNOTs(INT_MAX), target(startTargetQubits)
			{
				std::vector<unsigned int> controlQubits(endQubit - startQubit + 1);
				std::iota(controlQubits.begin(), controlQubits.end(), startQubit);

				nControlledNOTs.SetControlQubits(controlQubits);
				nControlledNOTs.SetStartAncillaQubits(startAncilla);
			}

			unsigned int Execute(RegisterClass& reg) override
			{
				unsigned int sQubit = BaseClass::BaseClass::getStartQubit();
				unsigned int eQubit = BaseClass::BaseClass::getEndQubit();
				unsigned int nrQubits = eQubit - sQubit + 1;
				unsigned int nrBasisStates = 1 << nrQubits;

				/*
					for (unsigned int state = 0; state < nrBasisStates; ++state)
					{
						unsigned int fval = func(state);
						if (!fval) continue;
						nControlledNOTs.ClearGates();
						
						unsigned int q = target;
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

						unsigned int realState = state << sQubit;

						unsigned int v = realState;
						for (unsigned int q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, sQubit + q);

							v >>= 1;
						}

						nControlledNOT.Execute(reg);

						// undo the x gates
						v = realState;
						for (unsigned int q = 0; q < nrQubits; ++q)
						{
							if ((v & 1) == 0)
								reg.ApplyGate(x, sQubit + q);

							v >>= 1;
						}
					}
				*/

				// the above code is not that good, I'll let it commented in case this one is too hard to understand
				// this one avoids applying x gates twice on the same qubit, between applying n controlled nots
				std::vector<bool> qubits(nrQubits, false);

				for (unsigned int state = 0; state < nrBasisStates; ++state)
				{
					unsigned int fval = func(state);
					if (!fval) continue;

					nControlledNOTs.ClearGates();
					
					for (unsigned int q = target; fval; ++q)
					{
						if (fval & 1)
						{
							QC::Gates::AppliedGate<MatrixClass> notGate(x.getRawOperatorMatrix(), q);
							nControlledNOTs.AddGate(notGate);
						}

						fval >>= 1;
					}

					unsigned int realState = state << sQubit;

					unsigned int v = realState;
					for (unsigned int q = 0; q < nrQubits; ++q)
					{
						if (qubits[q] != ((v & 1) == 0))
						{
							reg.ApplyGate(x, sQubit + q);
							qubits[q] = !qubits[q]; // x gate was applied
						}

						v >>= 1;
					}

					nControlledNOTs.Execute(reg);
				}

				// undo the x gates
				for (unsigned int q = 0; q < nrQubits; ++q)
					if (qubits[q])
						reg.ApplyGate(x, sQubit + q);

				return 0;
			}

			void setFunction(const Function& f)
			{
				func = f;
			}

		protected:
			Function func;

			QC::SubAlgo::NControlledGatesWithAncilla<VectorClass, MatrixClass> nControlledNOTs;
			QC::Gates::PauliXGate<MatrixClass> x;

			unsigned int target;
		};
	}
}
