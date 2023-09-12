#pragma once

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegisterCalculator {
	public:
		QubitRegisterCalculator() = default;
		virtual ~QubitRegisterCalculator() = default;

		inline void ApplyOneQubitGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;

			for (unsigned int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = state & qubitBit ? 1 : 0;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) +
					gateMatrix(row, 1) * registerStorage(state | qubitBit);
			}
		}

		inline void ApplyOneQubitGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;

#pragma omp parallel for schedule(static, 4096)
			for (long long int state = 0; state < NrBasisStates; ++state)
			{
				const unsigned int row = state & qubitBit ? 1 : 0;

				resultsStorage(state) = gateMatrix(row, 0) * registerStorage(state & notQubitBit) +
					gateMatrix(row, 1) * registerStorage(state | qubitBit);
			}
		}

		inline void ApplyTwoQubitsGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int orqubits = qubitBit | ctrlQubitBit;

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

		inline void ApplyTwoQubitsGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int orqubits = qubitBit | ctrlQubitBit;

#pragma omp parallel for schedule(static, 2048)
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

		inline void ApplyThreeQubitsGate(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int qubitBit2, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int notQubitBit2 = ~qubitBit2;
			const unsigned int ctrlqubits = ctrlQubitBit | qubitBit2;
			const unsigned int orqubits = qubitBit | qubitBit2;
			const unsigned int orallqubits = qubitBit | ctrlqubits;
			const unsigned int ctrlorqubit2 = ctrlQubitBit | qubitBit;

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

		inline void ApplyThreeQubitsGateOmp(const VectorClass& registerStorage, VectorClass& resultsStorage, const MatrixClass& gateMatrix, const unsigned int qubitBit, const unsigned int qubitBit2, const unsigned int ctrlQubitBit, const unsigned int NrBasisStates)
		{
			const unsigned int notQubitBit = ~qubitBit;
			const unsigned int notCtrlQubitBit = ~ctrlQubitBit;
			const unsigned int notQubitBit2 = ~qubitBit2;
			const unsigned int ctrlqubits = ctrlQubitBit | qubitBit2;
			const unsigned int orqubits = qubitBit | qubitBit2;
			const unsigned int orallqubits = qubitBit | ctrlqubits;
			const unsigned int ctrlorqubit2 = ctrlQubitBit | qubitBit;

#pragma omp parallel for schedule(static, 1024)
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
	};

}

