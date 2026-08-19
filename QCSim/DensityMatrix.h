#pragma once

#include "QubitRegisterCalculator.h"
#include "QubitRegister.h"
#include "SimpleGates.h"

#include <random>
#include <vector>
#include <complex>
#include <utility>

// A density-matrix quantum computing simulator.
//
// The state is stored as a rho = 2^n x 2^n complex matrix (MatrixClass), not as a vectorized
// super-operator. Unitary gates are applied as rho' = U rho U^dagger by reusing the exact same
// bit-mask gate kernels from QubitRegisterCalculator, but instead of building a full register and
// copying data back and forth, the kernels are fed directly with the columns and the rows of rho
// (Eigen block expressions used as fake 'registers').
//
// The left multiplication U rho is done by applying U to every column (the ket index); the right
// multiplication by U^dagger is done by applying the (elementwise) conjugated small operator to
// every row (the bra index). Only the small operator matrix is conjugated, never the big rho, and
// the calculators are held as member objects instead of being inherited from.
//
// Qubits are numbered from right to left, starting with zero (same convention as QubitRegister).

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd>
	class DensityMatrix
	{
	public:
		using GateClass = Gates::QuantumGateWithOp<MatrixClass>;

		// whatever rho.col(0) / rho.row(0) return - the Eigen block expressions used as fake 'registers'
		using ColXpr = decltype(std::declval<MatrixClass&>().col(0));
		using RowXpr = decltype(std::declval<MatrixClass&>().row(0));

		// two calculators, one instantiated for the column blocks, one for the row blocks
		using ColCalculator = QubitRegisterCalculator<ColXpr, MatrixClass>;
		using RowCalculator = QubitRegisterCalculator<RowXpr, MatrixClass>;

		DensityMatrix(size_t N = 3, unsigned int addseed = 0)
			: NrQubits(N), NrBasisStates(1ULL << N),
			rho(MatrixClass::Zero(NrBasisStates, NrBasisStates)),
			target(MatrixClass::Zero(NrBasisStates, NrBasisStates)),
			uniformZeroOne(0, 1)
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

			rho(0, 0) = 1.; // |0...0><0...0|
		}

		size_t getNrQubits() const { return NrQubits; }
		size_t getNrBasisStates() const { return NrBasisStates; }

		const MatrixClass& getDensityMatrix() const { return rho; }

		void Clear()
		{
			rho.setZero();
		}

		void setToBasisState(size_t State)
		{
			if (State >= NrBasisStates) return;

			rho.setZero();
			rho(State, State) = 1.;
		}

		void Reset()
		{
			setToBasisState(0);
		}

		// initialize rho = |psi><psi| from a statevector, very convenient for comparing against
		// the statevector simulator
		void setFromStatevector(const VectorClass& psi)
		{
			assert(static_cast<size_t>(psi.size()) == NrBasisStates);

			rho = psi * psi.adjoint();
		}

		// set the whole density matrix directly (the caller is responsible for it being a valid state:
		// Hermitian, positive semidefinite and unit trace)
		void setDensityMatrix(const MatrixClass& newRho)
		{
			assert(static_cast<size_t>(newRho.rows()) == NrBasisStates && static_cast<size_t>(newRho.cols()) == NrBasisStates);

			rho = newRho;
		}

		// set the state to a classical mixture of computational basis states: rho = sum_k p_k |s_k><s_k|
		// the weights are normalized to sum to 1 (negative or zero weights are ignored)
		void setToMixtureOfBasisStates(const std::vector<std::pair<size_t, double>>& mixture)
		{
			rho.setZero();

			double total = 0.;
			for (const auto& [state, weight] : mixture)
				if (weight > 0. && state < NrBasisStates) total += weight;

			if (total <= 0.) return;

			for (const auto& [state, weight] : mixture)
				if (weight > 0. && state < NrBasisStates)
					rho(static_cast<Eigen::Index>(state), static_cast<Eigen::Index>(state)) += weight / total;
		}

		// rho' = U rho U^dagger
		void ApplyGate(const GateClass& gate, size_t qubit, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0)
		{
			const size_t gateQubits = gate.getQubitsNumber();
			const MatrixClass& U = gate.getRawOperatorMatrix();
			const MatrixClass Uconj = U.conjugate(); // small, cheap to conjugate

			// rho <- U rho : apply U to every column (the ket index)
			ApplyGateToColumns(gate, U, gateQubits, qubit, controllingQubit1, controllingQubit2);

			// rho <- (U rho) U^dagger : right multiplication by U^dagger is applying the conjugated
			// small operator to every row (the bra index). Only the small operator was conjugated.
			ApplyGateToRows(gate, Uconj, gateQubits, qubit, controllingQubit1, controllingQubit2);
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

		// Generic completely positive trace preserving channel: rho' = sum_k E_k rho E_k^dagger
		// The Kraus operators are small matrices (2x2 for a single qubit, 4x4 for two qubits) acting
		// on the given qubit(s). This is the fundamental non-unitary operation; a unitary gate is just
		// the special case of a single Kraus operator E_0 = U.
		void ApplyChannel(const std::vector<MatrixClass>& kraus, size_t qubit, size_t controllingQubit1 = 0)
		{
			const MatrixClass original = rho;
			MatrixClass acc = MatrixClass::Zero(NrBasisStates, NrBasisStates);

			for (const auto& E : kraus)
			{
				const Gates::AppliedGate<MatrixClass> gate(E);
				const size_t gateQubits = gate.getQubitsNumber();
				const MatrixClass Econj = E.conjugate();

				rho = original;

				// rho <- E_k rho E_k^dagger (left mult on columns, right mult on rows with the
				// conjugated small operator - the big matrix is never adjointed)
				ApplyGateToColumns(gate, E, gateQubits, qubit, controllingQubit1, 0);
				ApplyGateToRows(gate, Econj, gateQubits, qubit, controllingQubit1, 0);

				acc += rho;
			}

			rho.swap(acc);
		}

		// ---- predefined single qubit noise channels ----

		// bit flip: rho' = (1 - p) rho + p X rho X
		void ApplyBitFlipNoise(size_t qubit, double p)
		{
			ApplyChannel({ std::sqrt(1. - p) * PauliI(), std::sqrt(p) * PauliX() }, qubit);
		}

		// phase flip: rho' = (1 - p) rho + p Z rho Z
		void ApplyPhaseFlipNoise(size_t qubit, double p)
		{
			ApplyChannel({ std::sqrt(1. - p) * PauliI(), std::sqrt(p) * PauliZ() }, qubit);
		}

		// depolarizing: rho' = (1 - p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
		void ApplyDepolarizingNoise(size_t qubit, double p)
		{
			const double s = std::sqrt(p / 3.);
			ApplyChannel({ std::sqrt(1. - p) * PauliI(), s * PauliX(), s * PauliY(), s * PauliZ() }, qubit);
		}

		// amplitude damping (|1> -> |0> relaxation with probability gamma)
		void ApplyAmplitudeDamping(size_t qubit, double gamma)
		{
			MatrixClass E0 = MatrixClass::Zero(2, 2);
			E0(0, 0) = 1.;
			E0(1, 1) = std::sqrt(1. - gamma);

			MatrixClass E1 = MatrixClass::Zero(2, 2);
			E1(0, 1) = std::sqrt(gamma);

			ApplyChannel({ E0, E1 }, qubit);
		}

		// phase damping / dephasing, suppresses the off diagonal coherences by lambda = sqrt(1 - gamma)
		void ApplyPhaseDamping(size_t qubit, double gamma)
		{
			MatrixClass E0 = MatrixClass::Zero(2, 2);
			E0(0, 0) = 1.;
			E0(1, 1) = std::sqrt(1. - gamma);

			MatrixClass E1 = MatrixClass::Zero(2, 2);
			E1(1, 1) = std::sqrt(gamma);

			ApplyChannel({ E0, E1 }, qubit);
		}

		// reset a qubit to |0>: E0 = |0><0|, E1 = |0><1|
		void ApplyReset(size_t qubit)
		{
			MatrixClass E0 = MatrixClass::Zero(2, 2);
			E0(0, 0) = 1.;

			MatrixClass E1 = MatrixClass::Zero(2, 2);
			E1(0, 1) = 1.;

			ApplyChannel({ E0, E1 }, qubit);
		}

		// ---- measurement ----

		double GetQubitProbability(size_t qubit) const
		{
			const size_t mask = 1ULL << qubit;

			double p1 = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
				if (i & mask)
					p1 += rho(i, i).real();

			return p1;
		}

		// sample a computational basis outcome for a single qubit and collapse
		size_t MeasureQubit(size_t qubit)
		{
			const size_t mask = 1ULL << qubit;

			double p0 = 0;
			for (size_t i = 0; i < NrBasisStates; ++i)
				if ((i & mask) == 0)
					p0 += rho(i, i).real();

			const double r = uniformZeroOne(rng);
			const size_t result = (r < p0) ? 0 : 1;
			const double pm = (result == 0) ? p0 : 1. - p0;

			CollapseQubit(qubit, result, pm);

			return result;
		}

		// non selective measurement of a qubit: rho' = P0 rho P0 + P1 rho P1
		// destroys the coherence between the two measurement sectors but keeps the populations
		void DephaseMeasure(size_t qubit)
		{
			const size_t mask = 1ULL << qubit;

			for (size_t i = 0; i < NrBasisStates; ++i)
				for (size_t j = 0; j < NrBasisStates; ++j)
					if (((i & mask) == 0) != ((j & mask) == 0))
						rho(i, j) = 0;
		}

		// ---- diagnostics ----

		std::complex<double> Trace() const
		{
			return rho.trace();
		}

		// Tr(rho^2), equals 1 for a pure state, < 1 for a mixed one
		double Purity() const
		{
			return (rho * rho).trace().real();
		}

		bool IsHermitian(double eps = 1E-10) const
		{
			return (rho - rho.adjoint()).norm() < eps;
		}

		double getBasisStateProbability(size_t State) const
		{
			if (State >= NrBasisStates) return 0;

			return rho(State, State).real();
		}

		// <O> = Tr(rho O), the caller should ensure O is Hermitian and take the real part
		std::complex<double> ExpectationValue(const MatrixClass& O) const
		{
			return (rho * O).trace();
		}

	protected:
		// applies a small gate to every column of rho (the ket index). Each column is a fake 'register'
		// (an Eigen block expression) fed directly to the reused QubitRegisterCalculator kernels; results
		// go in place or into the matching column of the target buffer, which is then swapped in - so no
		// per column temporaries and no copies back and forth.
		void ApplyGateToColumns(const GateClass& gate, const MatrixClass& gateMatrix, size_t gateQubits, size_t qubit, size_t controllingQubit1, size_t controllingQubit2)
		{
			const size_t qubitBit = 1ULL << qubit;

			bool swapStorage = true;
			for (size_t col = 0; col < NrBasisStates; ++col)
			{
				ColXpr src = rho.col(col);
				ColXpr dst = target.col(col);

				swapStorage = true;
				if (gateQubits == 1)
					colCalculator.ApplyOneQubitGate(gate, src, dst, gateMatrix, qubitBit, NrBasisStates, swapStorage);
				else if (gateQubits == 2)
				{
					const size_t ctrlQubitBit = 1ULL << controllingQubit1;
					colCalculator.ApplyTwoQubitsGate(gate, src, dst, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
				}
				else
				{
					const size_t qubitBit2 = 1ULL << controllingQubit1;
					const size_t ctrlQubitBit = 1ULL << controllingQubit2;
					colCalculator.ApplyThreeQubitsGate(gate, src, dst, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				}
			}

			if (swapStorage) rho.swap(target); // the results ended up in the target, make it the new rho
		}

		// applies U^dagger from the right (the bra index) by applying the conjugated small operator to
		// every row of rho. The passed gateMatrix is already the elementwise conjugate, so the kernels
		// that read the matrix (generic / diagonal / antidiagonal / controlled) are directly correct.
		// The hardcoded iSwap / iSwapDag kernels ignore the matrix, so their roles are swapped here to
		// realize the conjugate (iSwap^dagger = iSwapDag); swap gates are self conjugate.
		void ApplyGateToRows(const GateClass& gate, const MatrixClass& gateMatrix, size_t gateQubits, size_t qubit, size_t controllingQubit1, size_t controllingQubit2)
		{
			const size_t qubitBit = 1ULL << qubit;

			const GateClass* dispatchGate = &gate;
			static const Gates::iSwapGate<MatrixClass> iswap;
			static const Gates::iSwapDagGate<MatrixClass> iswapDag;
			if (gateQubits == 2)
			{
				if (gate.IsISwapGate()) dispatchGate = &iswapDag;
				else if (gate.IsISwapDagGate()) dispatchGate = &iswap;
			}

			bool swapStorage = true;
			for (size_t row = 0; row < NrBasisStates; ++row)
			{
				RowXpr src = rho.row(row);
				RowXpr dst = target.row(row);

				swapStorage = true;
				if (gateQubits == 1)
					rowCalculator.ApplyOneQubitGate(*dispatchGate, src, dst, gateMatrix, qubitBit, NrBasisStates, swapStorage);
				else if (gateQubits == 2)
				{
					const size_t ctrlQubitBit = 1ULL << controllingQubit1;
					rowCalculator.ApplyTwoQubitsGate(*dispatchGate, src, dst, gateMatrix, qubitBit, ctrlQubitBit, NrBasisStates, swapStorage);
				}
				else
				{
					const size_t qubitBit2 = 1ULL << controllingQubit1;
					const size_t ctrlQubitBit = 1ULL << controllingQubit2;
					rowCalculator.ApplyThreeQubitsGate(*dispatchGate, src, dst, gateMatrix, qubitBit, qubitBit2, ctrlQubitBit, NrBasisStates, swapStorage);
				}
			}

			if (swapStorage) rho.swap(target);
		}

		void CollapseQubit(size_t qubit, size_t result, double pm)
		{
			const size_t mask = 1ULL << qubit;
			const double invpm = (pm > 1E-20) ? 1. / pm : 0.;

			for (size_t i = 0; i < NrBasisStates; ++i)
			{
				const size_t bi = (i & mask) ? 1 : 0;
				for (size_t j = 0; j < NrBasisStates; ++j)
				{
					const size_t bj = (j & mask) ? 1 : 0;
					if (bi != result || bj != result)
						rho(i, j) = 0;
					else
						rho(i, j) *= invpm;
				}
			}
		}

		static MatrixClass PauliI()
		{
			MatrixClass m = MatrixClass::Identity(2, 2);
			return m;
		}

		static MatrixClass PauliX()
		{
			MatrixClass m = MatrixClass::Zero(2, 2);
			m(0, 1) = 1.;
			m(1, 0) = 1.;
			return m;
		}

		static MatrixClass PauliY()
		{
			MatrixClass m = MatrixClass::Zero(2, 2);
			m(0, 1) = std::complex<double>(0, -1);
			m(1, 0) = std::complex<double>(0, 1);
			return m;
		}

		static MatrixClass PauliZ()
		{
			MatrixClass m = MatrixClass::Zero(2, 2);
			m(0, 0) = 1.;
			m(1, 1) = -1.;
			return m;
		}

		size_t NrQubits;
		size_t NrBasisStates;

		MatrixClass rho;
		MatrixClass target; // reused scratch buffer, swapped in place of rho when a kernel needs a target

		// stateless calculators held as member objects (composition instead of inheritance)
		ColCalculator colCalculator;
		RowCalculator rowCalculator;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne;
	};

}
