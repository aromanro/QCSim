#pragma once


#ifdef _WIN32
#include <windows.h>
#undef min
#undef max
#endif // _WIN32

#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Eigen>

#include <random>
#include <complex>
#include <chrono>
#include <algorithm>

#include <iostream>
#include <iterator>
#include <iomanip>
#include <fstream>

#include <vector>
#include <thread>

#include "QuantumGate.h"


namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegisterCalculator {
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;

		QubitRegisterCalculator() = default;
		virtual ~QubitRegisterCalculator() = default;

		static inline void ApplyOneQubitGate(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isDiagonal())
			{
				swapStorage = false;
				for (size_t state = 0; state < NrBasisStates; ++state)
					registerStorage(state) *= state & qubitBit ? gateMatrix(1, 1) : gateMatrix(0, 0);
			}
			else
			{
				const size_t notQubitBit = ~qubitBit;

				if (gate.isAntidiagonal())
				{
					for (size_t state = 0; state < NrBasisStates; ++state)
						resultsStorage(state) = state & qubitBit ? gateMatrix(1, 0) * registerStorage(state & notQubitBit) : gateMatrix(0, 1) * registerStorage(state | qubitBit);
				}
				else
				{
					for (size_t state = 0; state < NrBasisStates; ++state)
					{
						const size_t row = state & qubitBit ? 1 : 0;

						resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) + gateMatrix(row, 1) * registerStorage(state | qubitBit);
					}
				}
			}
		}

		static inline void ApplyOneQubitGateOmp(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			//const size_t blockSize = OneQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			if (gate.isDiagonal())
			{
				swapStorage = false;
#pragma omp parallel for 
				//num_threads(processor_count) schedule(static, blockSize)
				for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
					registerStorage(state) *= state & qubitBit ? gateMatrix(1, 1) : gateMatrix(0, 0);
			}
			else
			{
				const size_t notQubitBit = ~qubitBit;

				if (gate.isAntidiagonal())
				{
#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, blockSize)
					for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
						resultsStorage(state) = state & qubitBit ? gateMatrix(1, 0) * registerStorage(state & notQubitBit) : gateMatrix(0, 1) * registerStorage(state | qubitBit);
				}
				else
				{
#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, blockSize)
					for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
					{
						const size_t row = state & qubitBit ? 1 : 0;

						resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) + gateMatrix(row, 1) * registerStorage(state | qubitBit);
					}
				}
			}
		}

		static inline void ApplyTwoQubitsGate(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isSwapGate())
				ApplySwapGate(registerStorage, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
			else if (gate.isControlled())
			{
				if (gate.isDiagonal())
					ApplyDiagonalControlGate(registerStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
				else if (gate.isAntidiagonal())
					ApplyAntidiagonalControlGate(registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
				else
					ApplyGenericControlGate(registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
			}
			else
			{
				const size_t notQubitBit = ~qubitBit;
				const size_t notCtrlQubitBit = ~ctrlQubitBit;
				const size_t orqubits = qubitBit | ctrlQubitBit;

				for (size_t state = 0; state < NrBasisStates; ++state)
				{
					const size_t row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const size_t m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & notCtrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
						gateMatrix(row, 1) * registerStorage((state & notCtrlQubitBit) | qubitBit) +                     // state & ~ctrlQubitBit |  qubitBit      : 01
						gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                         // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                                          // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
		}

		static inline void ApplyTwoQubitsGateOmp(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isSwapGate())
				ApplySwapGateOmp(registerStorage, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
			else if (gate.isControlled())
			{
				if (gate.isDiagonal())
					ApplyDiagonalControlGateOmp(registerStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
				else if (gate.isAntidiagonal())
					ApplyAntidiagonalControlGateOmp(registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
				else
					ApplyGenericControlGateOmp(registerStorage, resultsStorage, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates);
			}
			else
			{
				//const size_t blockSize = TwoQubitOmpLimit / divSchedule;
				//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
				//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

				const size_t notQubitBit = ~qubitBit;
				const size_t notCtrlQubitBit = ~ctrlQubitBit;
				const size_t orqubits = qubitBit | ctrlQubitBit;

#pragma omp parallel for 
				//num_threads(processor_count) schedule(static, blockSize)
				for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
				{
					const size_t row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const size_t m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & notCtrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
						gateMatrix(row, 1) * registerStorage((state & notCtrlQubitBit) | qubitBit) +                     // state & ~ctrlQubitBit |  qubitBit      : 01
						gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                         // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                                          // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
		}

		static inline void ApplyThreeQubitsGate(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isSwapGate())
			{
				swapStorage = false;
				const size_t orqubits = qubitBit | qubitBit2;

				for (size_t state = std::max(ctrlQubitBit, std::min(qubitBit, qubitBit2)); state < NrBasisStates; ++state)
				{
					if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
						continue;

					const bool q1 = state & qubitBit ? 1 : 0;
					if (q1 == 1)
						continue;
					const bool q2 = state & qubitBit2 ? 1 : 0;
					if (q2 == 0)
						continue;

					const size_t swapstate = state ^ orqubits;
					std::swap(registerStorage(state), registerStorage(swapstate));
				}
			}
			else if (gate.isControlled())
				ApplyThreeQubitsControlledGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
			else
				ApplyThreeQubitsGenericGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
		}

		static inline void ApplyThreeQubitsGateOmp(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isSwapGate())
			{
				//const size_t blockSize = ThreeQubitOmpLimit / divSchedule;
				//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
				//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

				swapStorage = false;
				const size_t orqubits = qubitBit | qubitBit2;

				// TODO: is it worth parallelizing the controlled swap gate?
#pragma omp parallel for 
				//num_threads(processor_count) schedule(static, blockSize)
				for (long long int state = std::max(ctrlQubitBit, std::min(qubitBit, qubitBit2)); state < NrBasisStates; ++state)
				{
					if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
						continue;

					const bool q1 = state & qubitBit ? 1 : 0;
					if (q1 == 1)
						continue;
					const bool q2 = state & qubitBit2 ? 1 : 0;
					if (q2 == 0)
						continue;

					const size_t swapstate = state ^ orqubits;
					std::swap(registerStorage(state), registerStorage(swapstate));
				}
			}
			else if (gate.isControlled())
				ApplyThreeQubitsControlledGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
			else
				ApplyThreeQubitsGenericGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
		}

	private:
		static inline void ApplySwapGate(VectorClass& registerStorage, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			const size_t orqubits = qubitBit | ctrlQubitBit;

			swapStorage = false;

			for (size_t state = std::min(qubitBit, ctrlQubitBit); state < NrBasisStates; ++state)
			{
				const bool q1 = state & qubitBit ? 1 : 0;
				if (q1 == 1)
					continue;
				const bool q2 = state & ctrlQubitBit ? 1 : 0;
				if (q2 == 0)
					continue;

				const size_t swapstate = state ^ orqubits;
				std::swap(registerStorage(state), registerStorage(swapstate));
			}
		}

		static inline void ApplySwapGateOmp(VectorClass& registerStorage, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			//const size_t blockSize = TwoQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			const size_t orqubits = qubitBit | ctrlQubitBit;

			swapStorage = false;

			// TODO: is it worth parallelizing the swap gate?
#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = std::min(qubitBit, ctrlQubitBit); state < static_cast<long long int>(NrBasisStates); ++state)
			{
				const bool q1 = state & qubitBit ? 1 : 0;
				if (q1 == 1)
					continue;
				const bool q2 = state & ctrlQubitBit ? 1 : 0;
				if (q2 == 0)
					continue;

				const size_t swapstate = state ^ orqubits;
				std::swap(registerStorage(state), registerStorage(swapstate));
			}
		}

		static inline void ApplyDiagonalControlGate(VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			swapStorage = false;

			for (size_t state = ctrlQubitBit; state < NrBasisStates; ++state)
			{
				const size_t ctrl = (state & ctrlQubitBit);
				if (ctrl == 0)
					continue;

				registerStorage(state) *= state & qubitBit ? gateMatrix(3, 3) : gateMatrix(2, 2);
			}
		}

		static inline void ApplyAntidiagonalControlGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;

			for (size_t state = 0; state < ctrlQubitBit; ++state)
				resultsStorage(state) = registerStorage(state);

			for (size_t state = ctrlQubitBit; state < NrBasisStates; ++state)
			{
				const size_t ctrl = (state & ctrlQubitBit);
				if (ctrl == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				resultsStorage(state) = state & qubitBit ? gateMatrix(3, 2) * registerStorage(state & notQubitBit) : gateMatrix(2, 3) * registerStorage(state | qubitBit);
			}
		}

		static inline void ApplyGenericControlGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t orqubits = qubitBit | ctrlQubitBit;

			for (size_t state = 0; state < ctrlQubitBit; ++state)
				resultsStorage(state) = registerStorage(state);

			for (size_t state = ctrlQubitBit; state < NrBasisStates; ++state)
			{
				const size_t ctrl = (state & ctrlQubitBit);
				if (ctrl == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				const size_t row = 2 | (state & qubitBit ? 1 : 0);
				const size_t m = state & notQubitBit;

				resultsStorage(state) = gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                     // state & ~qubitBit     |  ctrlQubitBit  : 10
					gateMatrix(row, 3) * registerStorage(state | orqubits);                                          // state |  ctrlQubitBit |  qubitBit      : 11
			}
		}


		static inline void ApplyDiagonalControlGateOmp(VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			//const size_t blockSize = TwoQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			swapStorage = false;

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = ctrlQubitBit; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				const size_t ctrl = (state & ctrlQubitBit);
				if (ctrl == 0)
					continue;

				registerStorage(state) *= state & qubitBit ? gateMatrix(3, 3) : gateMatrix(2, 2);
			}
		}

		static inline void ApplyAntidiagonalControlGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			//const size_t blockSize = TwoQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			const size_t notQubitBit = ~qubitBit;

			for (size_t state = 0; state < ctrlQubitBit; ++state)
				resultsStorage(state) = registerStorage(state);

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = ctrlQubitBit; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				const size_t ctrl = (state & ctrlQubitBit);
				if (ctrl == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				resultsStorage(state) = state & qubitBit ? gateMatrix(3, 2) * registerStorage(state & notQubitBit) : gateMatrix(2, 3) * registerStorage(state | qubitBit);
			}
		}

		static inline void ApplyGenericControlGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			//const size_t blockSize = TwoQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			const size_t notQubitBit = ~qubitBit;
			const size_t orqubits = qubitBit | ctrlQubitBit;

			for (size_t state = 0; state < ctrlQubitBit; ++state)
				resultsStorage(state) = registerStorage(state);

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = ctrlQubitBit; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				const size_t ctrl = (state & ctrlQubitBit);
				if (ctrl == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				const size_t row = 2 | (state & qubitBit ? 1 : 0);
				const size_t m = state & notQubitBit;

				resultsStorage(state) = gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                     // state & ~qubitBit     |  ctrlQubitBit  : 10
					gateMatrix(row, 3) * registerStorage(state | orqubits);                                          // state |  ctrlQubitBit |  qubitBit      : 11
			}
		}


		static inline void ApplyThreeQubitsControlledGate(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isControlQubit(1))
			{
				if (gate.isDiagonal())
					ApplyDiagonalCCGate(registerStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				else if (gate.isAntidiagonal())
					ApplyAntidiagonalCCGate(registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				else
					ApplyGenericCCGate(registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
			}
			else
			{
				const size_t notQubitBit = ~qubitBit;
				const size_t notCtrlQubitBit = ~ctrlQubitBit;
				const size_t notQubitBit2 = ~qubitBit2;
				const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
				const size_t orallqubits = qubitBit | ctrlqubits;
				const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;

				for (size_t state = 0; state < ctrlQubitBit; ++state)
					resultsStorage(state) = registerStorage(state);

				for (size_t state = ctrlQubitBit; state < NrBasisStates; ++state)
				{
					const size_t ctrl = (state & ctrlQubitBit);
					if (ctrl == 0)
					{
						resultsStorage(state) = registerStorage(state);
						continue;
					}

					const size_t row = 4 | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const size_t m = state & notQubitBit;
					const size_t m2 = state & notQubitBit2;
					const size_t mnctrlQubitBit = m & notCtrlQubitBit;

					resultsStorage(state) = gateMatrix(row, 4) * registerStorage((m & notQubitBit2) | ctrlQubitBit) +	  // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
						gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				                          // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
						gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								              // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
						gateMatrix(row, 7) * registerStorage(state | orallqubits);                                        // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
				}
			}
		}


		static inline void ApplyDiagonalCCGate(VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			size_t limit = std::max(ctrlQubitBit, qubitBit2);

			swapStorage = false;

			for (size_t state = limit; state < NrBasisStates; ++state)
			{
				if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
					continue;

				registerStorage(state) *= state & qubitBit ? gateMatrix(7, 7) : gateMatrix(6, 6);
			}
		}

		static inline void ApplyAntidiagonalCCGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			size_t limit = std::max(ctrlQubitBit, qubitBit2);

			const size_t notQubitBit = ~qubitBit;

			for (size_t state = 0; state < limit; ++state)
				resultsStorage(state) = registerStorage(state);

			for (size_t state = limit; state < NrBasisStates; ++state)
			{
				if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				resultsStorage(state) = state & qubitBit ? gateMatrix(7, 6) * registerStorage(state & notQubitBit) : gateMatrix(6, 7) * registerStorage(state | qubitBit);
			}
		}

		static inline void ApplyGenericCCGate(VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			size_t limit = std::max(ctrlQubitBit, qubitBit2);

			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;

			for (size_t state = 0; state < limit; ++state)
				resultsStorage(state) = registerStorage(state);

			for (size_t state = limit; state < NrBasisStates; ++state)
			{
				if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				const size_t row = 6 | (state & qubitBit ? 1 : 0);
				const size_t m = state & notQubitBit;
				const size_t m2 = state & notQubitBit2;
				const size_t mnctrlQubitBit = m & notCtrlQubitBit;

				resultsStorage(state) = gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +						  // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
					gateMatrix(row, 7) * registerStorage(state | orallqubits);                                        // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
			}
		}


		static inline void ApplyThreeQubitsGenericGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orqubits = qubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;


			for (size_t state = 0; state < NrBasisStates; ++state)
			{
				const size_t row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
				const size_t m = state & notQubitBit;
				const size_t m2 = state & notQubitBit2;
				const size_t mnctrlQubitBit = m & notCtrlQubitBit;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & notQubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
					gateMatrix(row, 1) * registerStorage((m2 & notCtrlQubitBit) | qubitBit) +				          // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
					gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					              // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
					gateMatrix(row, 3) * registerStorage((state & notCtrlQubitBit) | orqubits) +                      // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
					gateMatrix(row, 4) * registerStorage((m & notQubitBit2) | ctrlQubitBit) +				          // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
					gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				                          // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
					gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								              // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
					gateMatrix(row, 7) * registerStorage(state | orallqubits);                                        // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
			}
		}

		static inline void ApplyThreeQubitsControlledGateOmp(const GateClass& gate, VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			if (gate.isControlQubit(1))
			{
				if (gate.isDiagonal())
					ApplyDiagonalCCGateOmp(registerStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				else if (gate.isAntidiagonal())
					ApplyAntidiagonalCCGateOmp(registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				else
					ApplyGenericCCGateOmp(registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
			}
			else
			{
				//const size_t blockSize = ThreeQubitOmpLimit / divSchedule;
				//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
				//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

				const size_t notQubitBit = ~qubitBit;
				const size_t notCtrlQubitBit = ~ctrlQubitBit;
				const size_t notQubitBit2 = ~qubitBit2;
				const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
				const size_t orallqubits = qubitBit | ctrlqubits;
				const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;

				for (long long int state = 0; state < static_cast<long long int>(ctrlQubitBit); ++state)
					resultsStorage(state) = registerStorage(state);

#pragma omp parallel for 
				//num_threads(processor_count) schedule(static, blockSize)
				for (long long int state = ctrlQubitBit; state < static_cast<long long int>(NrBasisStates); ++state)
				{
					const size_t ctrl = (state & ctrlQubitBit);
					if (ctrl == 0)
					{
						resultsStorage(state) = registerStorage(state);
						continue;
					}

					const size_t row = 4 | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const size_t m = state & notQubitBit;
					const size_t m2 = state & notQubitBit2;
					const size_t mnctrlQubitBit = m & notCtrlQubitBit;

					resultsStorage(state) = gateMatrix(row, 4) * registerStorage((m & notQubitBit2) | ctrlQubitBit) +	  // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
						gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				                          // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
						gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								              // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
						gateMatrix(row, 7) * registerStorage(state | orallqubits);                                        // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
				}
			}
		}

		static inline void ApplyDiagonalCCGateOmp(VectorClass& registerStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			//const size_t blockSize = ThreeQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			long long limit = static_cast<long long>(std::max(ctrlQubitBit, qubitBit2));

			swapStorage = false;

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = limit; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
					continue;

				registerStorage(state) *= state & qubitBit ? gateMatrix(7, 7) : gateMatrix(6, 6);
			}
		}

		static inline void ApplyAntidiagonalCCGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			//const size_t blockSize = ThreeQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			long long limit = static_cast<long long>(std::max(ctrlQubitBit, qubitBit2));

			const size_t notQubitBit = ~qubitBit;

			for (size_t state = 0; state < limit; ++state)
				resultsStorage(state) = registerStorage(state);

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = limit; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				resultsStorage(state) = state & qubitBit ? gateMatrix(7, 6) * registerStorage(state & notQubitBit) : gateMatrix(6, 7) * registerStorage(state | qubitBit);
			}
		}

		static inline void ApplyGenericCCGateOmp(VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates, bool& swapStorage)
		{
			//const size_t blockSize = ThreeQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			long long limit = static_cast<long long>(std::max(ctrlQubitBit, qubitBit2));

			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;

			for (size_t state = 0; state < limit; ++state)
				resultsStorage(state) = registerStorage(state);

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = limit; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & ctrlQubitBit) == 0 || (state & qubitBit2) == 0)
				{
					resultsStorage(state) = registerStorage(state);
					continue;
				}

				const size_t row = 6 | (state & qubitBit ? 1 : 0);
				const size_t m = state & notQubitBit;
				const size_t m2 = state & notQubitBit2;
				const size_t mnctrlQubitBit = m & notCtrlQubitBit;

				resultsStorage(state) = gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +						  // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
					gateMatrix(row, 7) * registerStorage(state | orallqubits);                                        // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
			}
		}

		static inline void ApplyThreeQubitsGenericGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			//const size_t blockSize = ThreeQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orqubits = qubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;

#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				const size_t row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
				const size_t m = state & notQubitBit;
				const size_t m2 = state & notQubitBit2;
				const size_t mnctrlQubitBit = m & notCtrlQubitBit;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & notQubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
					gateMatrix(row, 1) * registerStorage((m2 & notCtrlQubitBit) | qubitBit) +				          // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
					gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					              // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
					gateMatrix(row, 3) * registerStorage((state & notCtrlQubitBit) | orqubits) +                      // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
					gateMatrix(row, 4) * registerStorage((m & notQubitBit2) | ctrlQubitBit) +				          // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
					gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				                          // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
					gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								              // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
					gateMatrix(row, 7) * registerStorage(state | orallqubits);                                        // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
			}
		}

	public:
		//*****************************************************************************************************************************************************************************************
		// 
		//  Measurements
		// 
		//*****************************************************************************************************************************************************************************************

		static inline size_t MeasureQubit(size_t NrBasisStates, VectorClass& registerStorage, size_t qubit, const double prob)
		{
			double accum = 0;

			const size_t measuredQubitMask = 1ULL << qubit;

			size_t pstate = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					pstate = i;
					break;
				}
			}

			size_t measuredState = 0ULL;
			size_t measuredStateMask = 0ULL;
			if ((pstate & measuredQubitMask) != 0)
			{
				measuredState = 1ULL;
				measuredStateMask = measuredQubitMask;
			}

			// find the norm
			accum = 0;
			for (size_t state = measuredStateMask; state < NrBasisStates; ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
			for (size_t state = 0; state < NrBasisStates; ++state)
				registerStorage[state] *= ((state & measuredQubitMask) == measuredStateMask) ? norm : 0;

			return measuredState;
		}

		static inline size_t MeasureQubitOmp(size_t NrBasisStates, VectorClass& registerStorage, size_t qubit, const double prob)
		{
			//const size_t blockSize = OneQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			double accum = 0;

			const size_t measuredQubitMask = 1ULL << qubit;

			size_t pstate = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					pstate = i;
					break;
				}
			}

			size_t measuredState = 0ULL;
			size_t measuredStateMask = 0ULL;
			if ((pstate & measuredQubitMask) != 0)
			{
				measuredState = 1ULL;
				measuredStateMask = measuredQubitMask;
			}

			// find the norm
			accum = 0;

#pragma omp parallel for reduction(+:accum) 
			//num_threads(processor_count) schedule(static, blockSize)		
			for (long long state = measuredStateMask; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
				registerStorage[state] *= ((state & measuredQubitMask) == measuredStateMask) ? norm : 0;

			return measuredState;
		}

		static inline size_t MeasureQubitNoCollapse(size_t NrBasisStates, VectorClass& registerStorage, size_t qubit, const double prob)
		{
			double accum = 0;

			const size_t measuredQubitMask = 1ULL << qubit;

			size_t pstate = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					pstate = i;
					break;
				}
			}

			return (pstate & measuredQubitMask) != 0 ? 1 : 0;
		}

		static inline double GetQubitProbability(size_t NrBasisStates, const VectorClass& registerStorage, size_t qubit)
		{
			double accum = 0;

			const size_t measuredQubitMask = 1ULL << qubit;

			for (size_t state = measuredQubitMask; state < NrBasisStates; ++state)
			{
				if (state & measuredQubitMask)
					accum += std::norm(registerStorage[state]);
			}

			return accum;
		}

		static inline double GetQubitProbabilityOmp(size_t NrBasisStates, const VectorClass& registerStorage, size_t qubit)
		{
			//const size_t blockSize = OneQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			double accum = 0;

			const size_t measuredQubitMask = 1ULL << qubit;

#pragma omp parallel for reduction(+:accum) 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long state = measuredQubitMask; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if (state & measuredQubitMask)
					accum += std::norm(registerStorage[state]);
			}

			return accum;
		}

		static inline size_t Measure(size_t NrBasisStates, VectorClass& registerStorage, size_t firstQubit, size_t secondQubit, const double prob)
		{
			double accum = 0;

			const size_t secondQubitp1 = secondQubit + 1;
			const size_t firstPartMask = (1ULL << firstQubit) - 1;
			const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;

			size_t measuredState = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					measuredState = i;
					break;
				}
			}
			measuredState &= measuredPartMask;

			// find the norm
			accum = 0;
			for (size_t state = measuredState; state < NrBasisStates; ++state)
			{
				if ((state & measuredPartMask) == measuredState)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
			for (size_t state = 0; state < NrBasisStates; ++state)
				registerStorage[state] *= ((state & measuredPartMask) == measuredState) ? norm : 0;

			return measuredState >> firstQubit;
		}

		static inline size_t MeasureOmp(size_t NrBasisStates, VectorClass& registerStorage, size_t firstQubit, size_t secondQubit, const double prob)
		{
			//const size_t blockSize = OneQubitOmpLimit / divSchedule;
			//const size_t NrBlocks = NrBasisStates / blockSize; // it's always a multiple of blockSize, at least two blocks
			//const auto processor_count = std::min<int>(GetNumberOfThreads(), NrBlocks);

			double accum = 0;

			const size_t secondQubitp1 = secondQubit + 1;
			const size_t firstPartMask = (1ULL << firstQubit) - 1;
			const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;

			size_t measuredState = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					measuredState = i;
					break;
				}
			}
			measuredState &= measuredPartMask;

			// find the norm
			accum = 0;

#pragma omp parallel for reduction(+:accum) 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long state = measuredState; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & measuredPartMask) == measuredState)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
#pragma omp parallel for 
			//num_threads(processor_count) schedule(static, blockSize)
			for (long long state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
				registerStorage[state] *= ((state & measuredPartMask) == measuredState) ? norm : 0;

			return measuredState >> firstQubit;
		}

		static inline size_t MeasureNoCollapse(size_t NrBasisStates, VectorClass& registerStorage, size_t firstQubit, size_t secondQubit, const double prob)
		{
			double accum = 0;

			const size_t secondQubitp1 = secondQubit + 1;
			const size_t firstPartMask = (1ULL << firstQubit) - 1;
			const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;

			size_t measuredState = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				accum += std::norm(registerStorage(i));
				if (prob <= accum)
				{
					measuredState = i;
					break;
				}
			}

			return (measuredState & measuredPartMask) >> firstQubit;
		}

		static int GetNumberOfThreads()
		{
			const size_t threads = std::thread::hardware_concurrency();
			return static_cast<int>(threads ? threads : GetCpuInfoNrThreads());
		}


		void SetMultithreading(bool enable = true)
		{
			enableMultithreading = enable;
		}

		bool GetMultithreading() const
		{
			return enableMultithreading;
		}

		//constexpr static auto cone = std::complex<double>(1.0, 0.0);

		//constexpr static int divSchedule = 4;
		constexpr static size_t OneQubitOmpLimit = 32768;
		constexpr static size_t TwoQubitOmpLimit = OneQubitOmpLimit;
		constexpr static size_t ThreeQubitOmpLimit = OneQubitOmpLimit;

	private:
		static size_t GetCpuInfoNrThreads()
		{
#ifdef _WIN32
			SYSTEM_INFO sysinfo;
			GetSystemInfo(&sysinfo);
			return sysinfo.dwNumberOfProcessors;
#else
			std::ifstream cpuinfo("/proc/cpuinfo");

			return std::count(std::istream_iterator<std::string>(cpuinfo), std::istream_iterator<std::string>(), std::string("processor"));
#endif
		}

		bool enableMultithreading = true;
	};

}

