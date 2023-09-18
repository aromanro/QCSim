#pragma once

#include <math.h>
#include <Eigen/eigen>

#include <random>
#include <complex>
#include <chrono>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <vector>

#include "QuantumGate.h"


namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegisterCalculator {
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;

		QubitRegisterCalculator() = default;
		virtual ~QubitRegisterCalculator() = default;

		static inline void ApplyOneQubitGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;

			if (gate.isDiagonal())
			{
				for (unsigned int state = 0; state < NrBasisStates; ++state)
					resultsStorage(state) = state & qubitBit ? gateMatrix(1, 1) * registerStorage(state | qubitBit) : gateMatrix(0, 0) * registerStorage(state & notQubitBit);
			}
			else if (gate.isAntidiagonal())
			{
				for (unsigned int state = 0; state < NrBasisStates; ++state)
					resultsStorage(state) = state & qubitBit ? gateMatrix(1, 0) * registerStorage(state & notQubitBit) : gateMatrix(0, 1) * registerStorage(state | qubitBit);
			}
			else
			{
				for (unsigned int state = 0; state < NrBasisStates; ++state)
				{
					const unsigned int row = state & qubitBit ? 1 : 0;

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) +
						gateMatrix(row, 1) * registerStorage(state | qubitBit);
				}
			}
		}

		static inline void ApplyOneQubitGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;

			if (gate.isDiagonal())
			{
#pragma omp parallel for
				//schedule(static, 8192)
				for (long long int state = 0; state < NrBasisStates; ++state)
					resultsStorage(state) = state & qubitBit ? gateMatrix(1, 1) * registerStorage(state | qubitBit) : gateMatrix(0, 0) * registerStorage(state & notQubitBit);
			}
			else if (gate.isAntidiagonal())
			{
#pragma omp parallel for 
				//schedule(static, 8192)
				for (long long int state = 0; state < NrBasisStates; ++state)
					resultsStorage(state) = state & qubitBit ? gateMatrix(1, 0) * registerStorage(state & notQubitBit) : gateMatrix(0, 1) * registerStorage(state | qubitBit);
			}
			else
			{
#pragma omp parallel for 
				//schedule(static, 8192)
				for (long long int state = 0; state < NrBasisStates; ++state)
				{
					const unsigned int row = state & qubitBit ? 1 : 0;

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) +
						gateMatrix(row, 1) * registerStorage(state | qubitBit);
				}
			}
		}

		static inline void ApplyTwoQubitsGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int orqubits = qubitBit | ctrlQubitBit;

			if (gate.isControlled())
			{
				for (unsigned int state = 0; state < ctrlQubitBit; ++state)
					resultsStorage(state) = registerStorage(state);

				for (unsigned int state = ctrlQubitBit; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & notCtrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
						gateMatrix(row, 1) * registerStorage(state & notCtrlQubitBit | qubitBit) +                       // state & ~ctrlQubitBit |  qubitBit      : 01
						gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                       // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                         // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
			else
			{
				for (unsigned int state = 0; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & notCtrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
						gateMatrix(row, 1) * registerStorage(state & notCtrlQubitBit | qubitBit) +                       // state & ~ctrlQubitBit |  qubitBit      : 01
						gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                       // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                         // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
		}

		static inline void ApplyTwoQubitsGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int orqubits = qubitBit | ctrlQubitBit;
			
			if (gate.isControlled())
			{
#pragma omp parallel for
				//schedule(static, 4096)
				for (long long int state = 0; state < ctrlQubitBit; ++state)
					resultsStorage(state) = registerStorage(state);

#pragma omp parallel for
				//schedule(static, 4096)
				for (long long int state = ctrlQubitBit; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & notCtrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
						gateMatrix(row, 1) * registerStorage(state & notCtrlQubitBit | qubitBit) +                       // state & ~ctrlQubitBit |  qubitBit      : 01
						gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                       // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                         // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
			else
			{
#pragma omp parallel for
				//schedule(static, 4096) 
				for (long long int state = 0; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit; // ensure it's not computed twice

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & notCtrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
						gateMatrix(row, 1) * registerStorage(state & notCtrlQubitBit | qubitBit) +                       // state & ~ctrlQubitBit |  qubitBit      : 01
						gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                       // state & ~qubitBit     |  ctrlQubitBit  : 10
						gateMatrix(row, 3) * registerStorage(state | orqubits);                         // state |  ctrlQubitBit |  qubitBit      : 11
				}
			}
		}

		static inline void ApplyThreeQubitsGate(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int qubitBit2, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int notQubitBit2 = ~qubitBit2;
			const unsigned int ctrlqubits = ctrlQubitBit | qubitBit2;
			const unsigned int orqubits = qubitBit | qubitBit2;
			const unsigned int orallqubits = qubitBit | ctrlqubits;
			const unsigned int ctrlorqubit2 = ctrlQubitBit | qubitBit;

			// TODO: This can be optimized further by taking advantage of the other control qubit (if it's a cc gate)

			if (gate.isControlled())
			{
				unsigned int limit = ctrlQubitBit;
				if (gate.isControlQubit(1))
					limit = std::max(limit, qubitBit2);

				for (unsigned int state = 0; state < limit; ++state)
					resultsStorage(state) = registerStorage(state);

				for (unsigned int state = limit; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit;
					const unsigned int m2 = state & notQubitBit2;
					const unsigned int mnctrlQubitBit = m & notCtrlQubitBit;

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & notQubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
						gateMatrix(row, 1) * registerStorage(m2 & notCtrlQubitBit | qubitBit) +				            // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
						gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					            // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
						gateMatrix(row, 3) * registerStorage(state & notCtrlQubitBit | orqubits) +            // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
						gateMatrix(row, 4) * registerStorage(m & notQubitBit2 | ctrlQubitBit) +				            // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
						gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				            // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
						gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								            // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
						gateMatrix(row, 7) * registerStorage(state | orallqubits);                            // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
				}
			}
			else
			{
				for (unsigned int state = 0; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit;
					const unsigned int m2 = state & notQubitBit2;
					const unsigned int mnctrlQubitBit = m & notCtrlQubitBit;

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & notQubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
						gateMatrix(row, 1) * registerStorage(m2 & notCtrlQubitBit | qubitBit) +				            // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
						gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					            // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
						gateMatrix(row, 3) * registerStorage(state & notCtrlQubitBit | orqubits) +            // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
						gateMatrix(row, 4) * registerStorage(m & notQubitBit2 | ctrlQubitBit) +				            // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
						gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				            // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
						gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								            // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
						gateMatrix(row, 7) * registerStorage(state | orallqubits);                            // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
				}
			}
		}

		static inline void ApplyThreeQubitsGateOmp(const GateClass& gate, const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int qubitBit2, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int notQubitBit2 = ~qubitBit2;
			const unsigned int ctrlqubits = ctrlQubitBit | qubitBit2;
			const unsigned int orqubits = qubitBit | qubitBit2;
			const unsigned int orallqubits = qubitBit | ctrlqubits;
			const unsigned int ctrlorqubit2 = ctrlQubitBit | qubitBit;

			// TODO: This can be optimized further by taking advantage of the other control qubit (if it's a cc gate)

			if (gate.isControlled())
			{
				long long limit = ctrlQubitBit;
				if (gate.isControlQubit(1))
					limit = std::max(limit, static_cast<long long>(qubitBit2));

#pragma omp parallel for
				//schedule(static, 2048)
				for (long long int state = 0; state < limit; ++state)
					resultsStorage(state) = registerStorage(state);
#pragma omp parallel for
				//schedule(static, 2048)
				for (long long int state = limit; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit;
					const unsigned int m2 = state & notQubitBit2;
					const unsigned int mnctrlQubitBit = m & notCtrlQubitBit;

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & notQubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
						gateMatrix(row, 1) * registerStorage(m2 & notCtrlQubitBit | qubitBit) +				            // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
						gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					            // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
						gateMatrix(row, 3) * registerStorage(state & notCtrlQubitBit | orqubits) +            // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
						gateMatrix(row, 4) * registerStorage(m & notQubitBit2 | ctrlQubitBit) +				            // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
						gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				            // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
						gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								            // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
						gateMatrix(row, 7) * registerStorage(state | orallqubits);                            // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
				}
			}
			else
			{
#pragma omp parallel for
				//schedule(static, 2048)
				for (long long int state = 0; state < NrBasisStates; ++state)
				{
					const unsigned int row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
					const unsigned int m = state & notQubitBit;
					const unsigned int m2 = state & notQubitBit2;
					const unsigned int mnctrlQubitBit = m & notCtrlQubitBit;

					resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & notQubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
						gateMatrix(row, 1) * registerStorage(m2 & notCtrlQubitBit | qubitBit) +				            // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
						gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					            // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
						gateMatrix(row, 3) * registerStorage(state & notCtrlQubitBit | orqubits) +            // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
						gateMatrix(row, 4) * registerStorage(m & notQubitBit2 | ctrlQubitBit) +				            // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
						gateMatrix(row, 5) * registerStorage(m2 | ctrlorqubit2) +				            // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
						gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								            // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
						gateMatrix(row, 7) * registerStorage(state | orallqubits);                            // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
				}
			}
		}

//*****************************************************************************************************************************************************************************************
// 
//  Measurements
// 
//*****************************************************************************************************************************************************************************************

		static inline unsigned int MeasureQubit(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int qubit, const double prob)
		{
			double accum = 0;

			const unsigned int measuredQubitMask = 1u << qubit;
			
			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{	
				if (state & measuredQubitMask)
					accum += std::norm(registerStorage[state]);
			}

			unsigned int measuredState = 0;
			if (prob <= accum)
				measuredState = 1;

			// find the norm
			accum = 0;
			const unsigned int measuredStateMask = measuredState << qubit;
			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
			
			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					registerStorage[state] *= norm;
				else
					registerStorage[state] = 0;
			}

			return measuredState;
		}

		static inline unsigned int MeasureQubitOmp(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int qubit, const double prob)
		{
			double accum = 0;

			const unsigned int measuredQubitMask = 1u << qubit;

#pragma omp parallel for reduction(+:accum) 
			for (long long state = 0; state < NrBasisStates; ++state)
			{
				if (state & measuredQubitMask)
					accum += std::norm(registerStorage[state]);
			}

			unsigned int measuredState = 0;
			if (prob <= accum)
				measuredState = 1;

			// find the norm
			accum = 0;

			const unsigned int measuredStateMask = measuredState << qubit;
#pragma omp parallel for reduction(+:accum) 			
			for (long long state = 0; state < NrBasisStates; ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					accum += std::norm(registerStorage[state]);
			}
			const double norm = 1. / sqrt(accum);

			// collapse
#pragma omp parallel for
			for (long long state = 0; state < NrBasisStates; ++state)
			{
				if ((state & measuredQubitMask) == measuredStateMask)
					registerStorage[state] *= norm;
				else
					registerStorage[state] = 0;
			}

			return measuredState;
		}

		static inline unsigned int MeasureQubitNoCollapse(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int qubit, const double prob)
		{
			double accum = 0;

			const unsigned int measuredQubitMask = 1u << qubit;

			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				if (state & measuredQubitMask)
					accum += std::norm(registerStorage[state]);
			}

			unsigned int measuredState = 0;
			if (prob <= accum)
				measuredState = 1;

			return measuredState;
		}


		// TODO: the following are awful, there might be a better way to optimize this based on cache locality
		// see the above that implement a single qubit measurement
		static inline unsigned int Measure(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int firstQubit, unsigned int secondQubit, const double prob)
		{
			double accum = 0;

			const unsigned int secondQubitp1 = secondQubit + 1;

			const unsigned int firstPartMask = (1u << firstQubit) - 1;
			const unsigned int measuredPartMask = (1u << secondQubitp1) - 1 - firstPartMask;
			const unsigned int secondPartMask = NrBasisStates - 1 - measuredPartMask - firstPartMask;

			const unsigned int secondPartMax = secondPartMask >> secondQubitp1;
			const unsigned int maxMeasuredState = measuredPartMask >> firstQubit;

			unsigned int measuredState = maxMeasuredState;

			double norm = 1;
			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;
				double stateProbability = 0;

				for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
				{
					const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
					for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						stateProbability += std::norm(registerStorage[secondPart | firstPartBits]);
				}

				accum += stateProbability;
				if (prob <= accum)
				{
					measuredState = state;
					norm = 1. / sqrt(stateProbability);
					break;
				}
			}

			//int cnt = 0;
			// collapse

			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;

				if (state == measuredState)
				{
					for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
					{
						const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
						for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						{
							registerStorage[secondPart | firstPartBits] *= norm;
							//++cnt;
						}
					}
				}
				else
				{
					for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
					{
						const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
						for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						{
							registerStorage[secondPart | firstPartBits] = 0;
							//++cnt;
						}
					}
				}
			}

			//assert(cnt == NrBasisStates);

			return measuredState;
		}

		static inline unsigned int MeasureOmp(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int firstQubit, unsigned int secondQubit, const double prob)
		{
			double accum = 0;

			const unsigned int secondQubitp1 = secondQubit + 1;

			const unsigned int firstPartMask = (1u << firstQubit) - 1;
			const unsigned int measuredPartMask = (1u << secondQubitp1) - 1 - firstPartMask;
			const unsigned int secondPartMask = NrBasisStates - 1 - measuredPartMask - firstPartMask;

			const unsigned int secondPartMax = secondPartMask >> secondQubitp1;
			const unsigned int maxMeasuredState = measuredPartMask >> firstQubit;

			unsigned int measuredState = maxMeasuredState;

			double norm = 1;
			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;
				double stateProbability = 0;
#pragma omp parallel for reduction(+:stateProbability) 
				//schedule(static, 8192)
				for (long long secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
				{
					const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
					for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						stateProbability += std::norm(registerStorage[secondPart | firstPartBits]);
				}

				accum += stateProbability;
				if (prob <= accum)
				{
					measuredState = state;
					norm = 1. / sqrt(stateProbability);
					break;
				}
			}

			//int cnt = 0;
			// collapse
//#pragma omp parallel for 
			//schedule(static, 8192)
			for (long long state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;

				if (state == measuredState)
				{
					for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
					{
						const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
						for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						{
							registerStorage[secondPart | firstPartBits] *= norm;
							//++cnt;
						}
					}
				}
				else
				{
					for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
					{
						const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
						for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						{
							registerStorage[secondPart | firstPartBits] = 0;
							//++cnt;
						}
					}
				}
			}

			//assert(cnt == NrBasisStates);

			return measuredState;
		}

		static inline unsigned int MeasureNoCollapse(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int firstQubit, unsigned int secondQubit, const double prob)
		{
			double accum = 0;

			const unsigned int secondQubitp1 = secondQubit + 1;

			const unsigned int firstPartMask = (1u << firstQubit) - 1;
			const unsigned int measuredPartMask = (1u << secondQubitp1) - 1 - firstPartMask;
			const unsigned int secondPartMask = NrBasisStates - 1 - measuredPartMask - firstPartMask;

			const unsigned int secondPartMax = secondPartMask >> secondQubitp1;
			const unsigned int maxMeasuredState = measuredPartMask >> firstQubit;

			unsigned int measuredState = maxMeasuredState;
			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;
				double stateProbability = 0;

				for (unsigned int secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
				{
					const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
					for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						stateProbability += std::norm(registerStorage[secondPart | firstPartBits]);
				}

				accum += stateProbability;
				if (prob <= accum)
				{
					measuredState = state;
					break;
				}
			}

			return measuredState;
		}

		static inline unsigned int MeasureNoCollapseOmp(unsigned int NrBasisStates, VectorClass& registerStorage, unsigned int firstQubit, unsigned int secondQubit, const double prob)
		{
			double accum = 0;

			const unsigned int secondQubitp1 = secondQubit + 1;

			const unsigned int firstPartMask = (1u << firstQubit) - 1;
			const unsigned int measuredPartMask = (1u << secondQubitp1) - 1 - firstPartMask;
			const unsigned int secondPartMask = NrBasisStates - 1 - measuredPartMask - firstPartMask;

			const unsigned int secondPartMax = secondPartMask >> secondQubitp1;
			const unsigned int maxMeasuredState = measuredPartMask >> firstQubit;

			unsigned int measuredState = maxMeasuredState;
			for (unsigned int state = 0; state <= maxMeasuredState; ++state)
			{
				const unsigned int stateRegBits = state << firstQubit;
				double stateProbability = 0;

#pragma omp parallel for reduction(+:stateProbability) 
				//schedule(static, 8192)
				for (long long secondPartBits = 0; secondPartBits <= secondPartMax; ++secondPartBits)
				{
					const unsigned int secondPart = (secondPartBits << secondQubitp1) | stateRegBits;
					for (unsigned int firstPartBits = 0; firstPartBits <= firstPartMask; ++firstPartBits)
						stateProbability += std::norm(registerStorage[secondPart | firstPartBits]);
				}

				accum += stateProbability;
				if (prob <= accum)
				{
					measuredState = state;
					break;
				}
			}

			return measuredState;
		}
	};

}

