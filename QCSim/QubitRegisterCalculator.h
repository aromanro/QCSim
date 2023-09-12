#pragma once

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegisterCalculator {
	public:
		QubitRegisterCalculator() = default;
		virtual ~QubitRegisterCalculator() = default;

		inline void ApplyOneQubitGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int NrBasisStates)
		{
			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = state & qubitBit ? 1 : 0;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & ~qubitBit) +
					gateMatrix(row, 1) * registerStorage(state | qubitBit);
			}
		}

		inline void ApplyOneQubitGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int NrBasisStates)
		{
#pragma omp parallel for schedule(static, 4096)
			for (long long int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = state & qubitBit ? 1 : 0;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & ~qubitBit) +
					gateMatrix(row, 1) * registerStorage(state | qubitBit);
			}
		}

		inline void ApplyTwoQubitsGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
				const unsigned int m = state & ~qubitBit; // ensure it's not computed twice

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & ~ctrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
					gateMatrix(row, 1) * registerStorage(state & ~ctrlQubitBit | qubitBit) +                       // state & ~ctrlQubitBit |  qubitBit      : 01
					gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                       // state & ~qubitBit     |  ctrlQubitBit  : 10
					gateMatrix(row, 3) * registerStorage(state | qubitBit | ctrlQubitBit);                         // state |  ctrlQubitBit |  qubitBit      : 11
			}
		}

		inline void ApplyTwoQubitsGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
#pragma omp parallel for schedule(static, 2048)
			for (long long int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = (state & ctrlQubitBit ? 2 : 0) | (state & qubitBit ? 1 : 0);
				const unsigned int m = state & ~qubitBit; // ensure it's not computed twice

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(m & ~ctrlQubitBit) +                  // state & ~ctrlQubitBit & ~qubitBit      : 00
					gateMatrix(row, 1) * registerStorage(state & ~ctrlQubitBit | qubitBit) +                       // state & ~ctrlQubitBit |  qubitBit      : 01
					gateMatrix(row, 2) * registerStorage(m | ctrlQubitBit) +                                       // state & ~qubitBit     |  ctrlQubitBit  : 10
					gateMatrix(row, 3) * registerStorage(state | qubitBit | ctrlQubitBit);                         // state |  ctrlQubitBit |  qubitBit      : 11
			}
		}

		inline void ApplyThreeQubitsGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int qubitBit2, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
				const unsigned int m = state & ~qubitBit;
				const unsigned int m2 = state & ~qubitBit2;
				const unsigned int ctrlqubits = ctrlQubitBit | qubitBit2;
				const unsigned int mnctrlQubitBit = m & ~ctrlQubitBit;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & ~qubitBit2) +         // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
					gateMatrix(row, 1) * registerStorage(m2 & ~ctrlQubitBit | qubitBit) +				            // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
					gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					            // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
					gateMatrix(row, 3) * registerStorage(state & ~ctrlQubitBit | qubitBit | qubitBit2) +            // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
					gateMatrix(row, 4) * registerStorage(m & ~qubitBit2 | ctrlQubitBit) +				            // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
					gateMatrix(row, 5) * registerStorage(m2 | ctrlQubitBit | qubitBit) +				            // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
					gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								            // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
					gateMatrix(row, 7) * registerStorage(state | qubitBit | ctrlqubits);                            // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
			}
		}

		inline void ApplyThreeQubitsGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int qubitBit2, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
#pragma omp parallel for schedule(static, 1024)
			for (long long int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = (state & ctrlQubitBit ? 4 : 0) | (state & qubitBit2 ? 2 : 0) | (state & qubitBit ? 1 : 0);
				const unsigned int m = state & ~qubitBit;
				const unsigned int m2 = state & ~qubitBit2;
				const unsigned int ctrlqubits = ctrlQubitBit | qubitBit2;
				const unsigned int mnctrlQubitBit = m & ~ctrlQubitBit;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(mnctrlQubitBit & ~qubitBit2) +          // state & ~ctrlQubitBit & ~qubitBit2    & ~qubitBit      : 000
					gateMatrix(row, 1) * registerStorage(m2 & ~ctrlQubitBit | qubitBit) +				             // state & ~ctrlQubitBit & ~qubitBit2    |  qubitBit      : 001
					gateMatrix(row, 2) * registerStorage(mnctrlQubitBit | qubitBit2) +					             // state & ~ctrlQubitBit & ~qubitBit     |  qubitBit2     : 010
					gateMatrix(row, 3) * registerStorage(state & ~ctrlQubitBit | qubitBit | qubitBit2) +             // state & ~ctrlQubitBit |  qubitBit2    |  qubitBit      : 011
					gateMatrix(row, 4) * registerStorage(m & ~qubitBit2 | ctrlQubitBit) +				             // state & ~qubitBit2    & ~qubitBit     |  ctrlQubitBit  : 100
					gateMatrix(row, 5) * registerStorage(m2 | ctrlQubitBit | qubitBit) +				             // state & ~qubitBit2    |  ctrlQubitBit |  qubitBit      : 101
					gateMatrix(row, 6) * registerStorage(m | ctrlqubits) +								             // state & ~qubitBit     |  ctrlQubitBit |  qubitBit2     : 110
					gateMatrix(row, 7) * registerStorage(state | qubitBit | ctrlqubits);                             // state |  ctrlQubitBit |  qubitBit2    |  qubitBit      : 111
			}
		}
	};

}

