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
			const auto processor_count = GetNumberOfThreads();

			if (gate.isDiagonal())
			{
				swapStorage = false;
#pragma omp parallel for num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
				for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
					registerStorage(state) *= state & qubitBit ? gateMatrix(1, 1) : gateMatrix(0, 0);
			}
			else
			{
				const size_t notQubitBit = ~qubitBit;

				if (gate.isAntidiagonal())
				{
#pragma omp parallel for num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
					for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
						resultsStorage(state) = state & qubitBit ? gateMatrix(1, 0) * registerStorage(state & notQubitBit) : gateMatrix(0, 1) * registerStorage(state | qubitBit);
				}
				else
				{
#pragma omp parallel for num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
					for (long long int state = 0; state < static_cast<long long int>(NrBasisStates); ++state)
					{
						const size_t row = state & qubitBit ? 1 : 0;

						resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) + gateMatrix(row, 1) * registerStorage(state | qubitBit);
					}
				}
			}
		}

		static inline void ApplyTwoQubitsGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t orqubits = qubitBit | ctrlQubitBit;

			if (gate.isControlled())
			{
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
					const size_t m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                     // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                                          // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
			else
			{
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

		static inline void ApplyTwoQubitsGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t orqubits = qubitBit | ctrlQubitBit;
			const auto processor_count = GetNumberOfThreads();

			if (gate.isControlled())
			{
				for (size_t state = 0; state < ctrlQubitBit; ++state)
					resultsStorage(state) = registerStorage(state);

#pragma omp parallel for num_threads(processor_count) schedule(static, TwoQubitOmpLimit / divSchedule)
				for (long long int state = ctrlQubitBit; state < static_cast<long long int>(NrBasisStates); ++state)
				{
					const size_t ctrl = (state & ctrlQubitBit);
					if (ctrl == 0)
					{
						resultsStorage(state) = registerStorage(state);
						continue;
					}

					const size_t row = 2 | (state & qubitBit ? 1 : 0);
					const size_t m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                     // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                                          // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
			else
			{
#pragma omp parallel for num_threads(processor_count) schedule(static, TwoQubitOmpLimit / divSchedule)
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

		static inline void ApplyThreeQubitsGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			if (gate.isControlled())
				ApplyThreeQubitsControlledGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
			else
				ApplyThreeQubitsGenericGate(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
		}

		static inline void ApplyThreeQubitsGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			if (gate.isControlled())
				ApplyThreeQubitsControlledGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
			else
				ApplyThreeQubitsGenericGateOmp(gate, registerStorage, resultsStorage, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates);
		}

	private:
		static inline void ApplyThreeQubitsControlledGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orqubits = qubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;

			size_t limit = ctrlQubitBit;
			if (gate.isControlQubit(1))
				limit = std::max(limit, qubitBit2);

			for (size_t state = 0; state < limit; ++state)
				resultsStorage(state) = registerStorage(state);

			if (gate.isControlQubit(1))
				for (size_t state = limit; state < NrBasisStates; ++state)
				{
					const size_t ctrl = (state & ctrlQubitBit);
					const size_t ctrl2 = (state & qubitBit2);
					if (ctrl == 0 || ctrl2 == 0)
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
			else
				for (size_t state = limit; state < NrBasisStates; ++state)
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

		static inline void ApplyThreeQubitsControlledGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orqubits = qubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;
			const auto processor_count = GetNumberOfThreads();

			long long limit = ctrlQubitBit;
			if (gate.isControlQubit(1))
				limit = std::max(limit, static_cast<long long>(qubitBit2));

			for (long long int state = 0; state < static_cast<long long int>(limit); ++state)
				resultsStorage(state) = registerStorage(state);

			if (gate.isControlQubit(1))
#pragma omp parallel for num_threads(processor_count) schedule(static, ThreeQubitOmpLimit / divSchedule)
				for (long long int state = limit; state < static_cast<long long int>(NrBasisStates); ++state)
				{
					const size_t ctrl = (state & ctrlQubitBit);
					const size_t ctrl2 = (state & qubitBit2);
					if (ctrl == 0 || ctrl2 == 0)
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
			else
#pragma omp parallel for num_threads(processor_count) schedule(static, ThreeQubitOmpLimit / divSchedule)
				for (long long int state = limit; state < static_cast<long long int>(NrBasisStates); ++state)
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

		static inline void ApplyThreeQubitsGenericGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const size_t qubitBit, const size_t qubitBit2, const size_t ctrlQubitBit, const size_t NrBasisStates)
		{
			const size_t notQubitBit = ~qubitBit;
			const size_t notCtrlQubitBit = ~ctrlQubitBit;
			const size_t notQubitBit2 = ~qubitBit2;
			const size_t ctrlqubits = ctrlQubitBit | qubitBit2;
			const size_t orqubits = qubitBit | qubitBit2;
			const size_t orallqubits = qubitBit | ctrlqubits;
			const size_t ctrlorqubit2 = ctrlQubitBit | qubitBit;
			const auto processor_count = GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, ThreeQubitOmpLimit / divSchedule)
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
			double accum = 0;
			const auto processor_count = GetNumberOfThreads();

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

#pragma omp parallel for reduction(+:accum) num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)		
			for (long long state = measuredStateMask; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
#pragma omp parallel for num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
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
			double accum = 0;
			const auto processor_count = GetNumberOfThreads();

			const size_t measuredQubitMask = 1ULL << qubit;

#pragma omp parallel for reduction(+:accum) num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
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
			const auto processor_count = GetNumberOfThreads();
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

#pragma omp parallel for reduction(+:accum) num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
			for (long long state = measuredState; state < static_cast<long long int>(NrBasisStates); ++state)
			{
				if ((state & measuredPartMask) == measuredState)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
#pragma omp parallel for num_threads(processor_count) schedule(static, OneQubitOmpLimit / divSchedule)
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

		constexpr static int divSchedule = 4;
		constexpr static size_t OneQubitOmpLimit = 2048;
		constexpr static size_t TwoQubitOmpLimit = OneQubitOmpLimit / 2;
		constexpr static size_t ThreeQubitOmpLimit = OneQubitOmpLimit / 4;

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

