#include <iostream>
#include <array>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "Tests.h"

#include "QubitRegister.h"
#include "DensityMatrix.h"

#include <unsupported/Eigen/KroneckerProduct>

#define _USE_MATH_DEFINES
#include <math.h>

#define DM_NR_QUBITS_LIMIT 6

static void DM_FillOneQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::HyGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SDGGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::TGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::TDGGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PhaseShiftGate<>>(0.38 * M_PI));
	gates.emplace_back(std::make_shared<QC::Gates::PauliXGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PauliYGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PauliZGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTDagGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::RxGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::RyGate<>>(M_PI / 7));
	gates.emplace_back(std::make_shared<QC::Gates::RzGate<>>(M_PI / 5));
	gates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 3, M_PI / 5));
}

static void DM_FillTwoQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::SwapGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::iSwapDagGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::DecrementGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledYGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledZGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledHadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledPhaseShiftGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledUGate<>>(M_PI / 3, M_PI / 7));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRxGate<>>(M_PI / 5));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRyGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRzGate<>>(M_PI / 7));
}

static Eigen::MatrixXcd DM_FourierMatrix(Eigen::Index dimension)
{
	Eigen::MatrixXcd fourier(dimension, dimension);
	const double scale = 1. / std::sqrt(static_cast<double>(dimension));
	for (Eigen::Index row = 0; row < dimension; ++row)
		for (Eigen::Index col = 0; col < dimension; ++col)
			fourier(row, col) = scale * std::polar(1., 2. * M_PI * static_cast<double>(row * col) / static_cast<double>(dimension));
	return fourier;
}

static void DM_FillThreeQubitGates(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::ToffoliGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::FredkinGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::CCZGate<>>());

	// A dense complex gate exercises the generic three-qubit kernel rather than a specialized path.
	gates.emplace_back(std::make_shared<QC::Gates::ThreeQubitsGate<>>(DM_FourierMatrix(8)));
}

template<class Callable>
static bool DM_ExpectInvalidArgument(Callable&& callable, const char* description)
{
	try
	{
		callable();
	}
	catch (const std::invalid_argument&)
	{
		return true;
	}
	catch (const std::exception& ex)
	{
		std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
		return false;
	}

	std::cout << description << " did not throw std::invalid_argument" << std::endl;
	return false;
}

template<class Callable>
static bool DM_ExpectDomainError(Callable&& callable, const char* description)
{
	try
	{
		callable();
	}
	catch (const std::domain_error&)
	{
		return true;
	}
	catch (const std::exception& ex)
	{
		std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
		return false;
	}

	std::cout << description << " did not throw std::domain_error" << std::endl;
	return false;
}

// checks that a density matrix equals the outer product |psi><psi| of the statevector simulator
static bool CompareWithStatevector(const QC::DensityMatrix<>& dm, const QC::QubitRegister<>& reg, int nrQubits, double eps = 1E-9)
{
	const Eigen::MatrixXcd expected = reg.getDensityMatrix();
	const Eigen::MatrixXcd& got = dm.getDensityMatrix();

	if (expected.rows() != got.rows() || expected.cols() != got.cols())
	{
		std::cout << "Density matrix dimensions mismatch for " << nrQubits << " qubits" << std::endl;
		return false;
	}

	for (int i = 0; i < expected.rows(); ++i)
		for (int j = 0; j < expected.cols(); ++j)
			if (!approxEqual(expected(i, j), got(i, j), eps))
			{
				std::cout << "Density matrix element (" << i << ", " << j << ") mismatch for " << nrQubits
					<< " qubits: " << expected(i, j) << " (statevector) vs " << got(i, j) << " (density matrix)" << std::endl;
				return false;
			}

	return true;
}

// applies random unitary circuits to both a statevector simulator and the density matrix simulator
// and checks that rho == |psi><psi| after every gate
static bool DensityMatrixUnitaryTest()
{
	std::cout << "\nDensity matrix simulator - unitary circuits compared against the statevector simulator" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		for (int t = 0; t < 10; ++t)
		{
			QC::DensityMatrix<> dm(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);

			for (int i = 0; i < lim; ++i)
			{
				const int g = gateDistr(gen);
				const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
				if (twoQubitsGate && nrQubits < 2) continue;

				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				dm.ApplyGate(*gates[g], qubit1, qubit2);
				reg.ApplyGate(*gates[g], qubit1, qubit2);

				if (!CompareWithStatevector(dm, reg, nrQubits))
					return false;
			}

			// sanity: a pure state must have unit trace and purity 1
			if (!approxEqual(dm.Trace().real(), 1., 1E-9) || !approxEqual(dm.Trace().imag(), 0., 1E-9))
			{
				std::cout << "Trace is not 1 for a pure state: " << dm.Trace() << std::endl;
				return false;
			}

			if (!approxEqual(dm.Purity(), 1., 1E-9))
			{
				std::cout << "Purity is not 1 for a pure state: " << dm.Purity() << std::endl;
				return false;
			}

			if (!dm.IsHermitian(1E-9))
			{
				std::cout << "Density matrix is not Hermitian" << std::endl;
				return false;
			}
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

// Independently embeds a small operator in the full Hilbert space. Local bit zero corresponds to
// qubits[0], local bit one to qubits[1], and so on, matching the simulator's gate convention.
static Eigen::MatrixXcd DM_EmbedOperator(const Eigen::MatrixXcd& op, const std::vector<size_t>& qubits, size_t nrQubits)
{
	const size_t nrBasisStates = 1ULL << nrQubits;
	size_t operatedMask = 0;
	for (const size_t qubit : qubits) operatedMask |= 1ULL << qubit;
	const size_t unaffectedMask = (nrBasisStates - 1) ^ operatedMask;

	Eigen::MatrixXcd embedded = Eigen::MatrixXcd::Zero(static_cast<Eigen::Index>(nrBasisStates), static_cast<Eigen::Index>(nrBasisStates));
	for (size_t row = 0; row < nrBasisStates; ++row)
		for (size_t col = 0; col < nrBasisStates; ++col)
		{
			if ((row & unaffectedMask) != (col & unaffectedMask)) continue;

			size_t localRow = 0;
			size_t localCol = 0;
			for (size_t localQubit = 0; localQubit < qubits.size(); ++localQubit)
			{
				localRow |= ((row >> qubits[localQubit]) & 1ULL) << localQubit;
				localCol |= ((col >> qubits[localQubit]) & 1ULL) << localQubit;
			}
			embedded(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) =
				op(static_cast<Eigen::Index>(localRow), static_cast<Eigen::Index>(localCol));
		}

	return embedded;
}

static Eigen::MatrixXcd DM_DeterministicMixedState(size_t nrQubits)
{
	const Eigen::Index dimension = static_cast<Eigen::Index>(1ULL << nrQubits);
	Eigen::VectorXcd psi1(dimension);
	Eigen::VectorXcd psi2(dimension);
	for (Eigen::Index i = 0; i < dimension; ++i)
	{
		psi1(i) = std::complex<double>(static_cast<double>(i + 1), static_cast<double>(i % 3) - 1.);
		psi2(i) = std::complex<double>(static_cast<double>((2 * i + 1) % 7) - 3., static_cast<double>((3 * i + 2) % 5) - 2.);
	}
	psi1.normalize();
	psi2.normalize();
	return 0.37 * (psi1 * psi1.adjoint()) + 0.63 * (psi2 * psi2.adjoint());
}

// Cross-checks the column/row gate kernels against an independently embedded dense operator on a
// genuinely mixed state. Both specialized gate objects and recorded AppliedGate objects are used.
static bool DensityMatrixDenseGateReferenceTest()
{
	std::cout << "\nDensity matrix simulator - gates compared against independently embedded dense operators" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> oneQubitGates;
	oneQubitGates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
	oneQubitGates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 5., M_PI / 7.));

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> twoQubitGates;
	twoQubitGates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
	twoQubitGates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
	twoQubitGates.emplace_back(std::make_shared<QC::Gates::iSwapDagGate<>>());
	twoQubitGates.emplace_back(std::make_shared<QC::Gates::TwoQubitsGate<>>(DM_FourierMatrix(4)));

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> threeQubitGates;
	DM_FillThreeQubitGates(threeQubitGates);

	constexpr size_t nrQubits = 4;
	const Eigen::MatrixXcd initial = DM_DeterministicMixedState(nrQubits);
	auto checkPlacement = [&](const QC::Gates::QuantumGateWithOp<>& gate, const std::vector<size_t>& qubits) -> bool
	{
		const Eigen::MatrixXcd embedded = DM_EmbedOperator(gate.getRawOperatorMatrix(), qubits, nrQubits);
		const Eigen::MatrixXcd expected = embedded * initial * embedded.adjoint();

		QC::DensityMatrix<> direct(nrQubits);
		direct.setDensityMatrix(initial);
		direct.ApplyGate(gate, qubits[0], qubits.size() > 1 ? qubits[1] : 0, qubits.size() > 2 ? qubits[2] : 0);

		QC::DensityMatrix<> recorded(nrQubits);
		recorded.setDensityMatrix(initial);
		recorded.ApplyGate(QC::Gates::AppliedGate<>(gate.getRawOperatorMatrix(), qubits[0],
			qubits.size() > 1 ? qubits[1] : 0, qubits.size() > 2 ? qubits[2] : 0));

		if ((direct.getDensityMatrix() - expected).norm() > 1E-9 ||
			(recorded.getDensityMatrix() - expected).norm() > 1E-9)
		{
			std::cout << "Dense gate reference mismatch for a " << qubits.size() << "-qubit gate on";
			for (const size_t qubit : qubits) std::cout << " " << qubit;
			std::cout << std::endl;
			return false;
		}
		return true;
	};

	for (const auto& gate : oneQubitGates)
		for (size_t q1 = 0; q1 < nrQubits; ++q1)
			if (!checkPlacement(*gate, { q1 })) return false;

	for (const auto& gate : twoQubitGates)
		for (size_t q1 = 0; q1 < nrQubits; ++q1)
			for (size_t q2 = 0; q2 < nrQubits; ++q2)
				if (q1 != q2 && !checkPlacement(*gate, { q1, q2 })) return false;

	for (const auto& gate : threeQubitGates)
		for (size_t q1 = 0; q1 < nrQubits; ++q1)
			for (size_t q2 = 0; q2 < nrQubits; ++q2)
				for (size_t q3 = 0; q3 < nrQubits; ++q3)
					if (q1 != q2 && q1 != q3 && q2 != q3 && !checkPlacement(*gate, { q1, q2, q3 })) return false;

	std::cout << "Success" << std::endl;
	return true;
}

// Exercises every specialized three-qubit path plus a dense complex 8x8 gate. Every ordered
// placement is checked, including non-adjacent placements in registers larger than three qubits.
static bool DensityMatrixThreeQubitGateTest()
{
	std::cout << "\nDensity matrix simulator - three qubit gates compared against the statevector simulator" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillThreeQubitGates(gates);

	QC::Gates::HadamardGate<> h;
	QC::Gates::RyGate<> ry(M_PI / 7.);
	QC::Gates::CNOTGate<> cnot;

	for (int nrQubits = 3; nrQubits <= 5; ++nrQubits)
	{
		for (const auto& gate : gates)
		{
			for (int q1 = 0; q1 < nrQubits; ++q1)
				for (int q2 = 0; q2 < nrQubits; ++q2)
					for (int q3 = 0; q3 < nrQubits; ++q3)
					{
						if (q1 == q2 || q1 == q3 || q2 == q3) continue;

						QC::DensityMatrix<> dm(nrQubits);
						QC::QubitRegister<> reg(nrQubits);
						for (int q = 0; q < nrQubits; ++q)
						{
							if (q % 2 == 0)
							{
								dm.ApplyGate(h, q);
								reg.ApplyGate(h, q);
							}
							else
							{
								dm.ApplyGate(ry, q);
								reg.ApplyGate(ry, q);
							}
						}
						for (int q = 1; q < nrQubits; ++q)
						{
							dm.ApplyGate(cnot, q, q - 1);
							reg.ApplyGate(cnot, q, q - 1);
						}

						dm.ApplyGate(*gate, q1, q2, q3);
						reg.ApplyGate(*gate, q1, q2, q3);
						if (!CompareWithStatevector(dm, reg, nrQubits)) return false;
					}
		}
	}

	std::cout << "Success" << std::endl;
	return true;
}

// checks the noise channels: trace preservation, that noise turns a pure state into a mixed one,
// and the analytic amplitude damping result
static bool DensityMatrixChannelsTest()
{
	std::cout << "\nDensity matrix simulator - quantum channels" << std::endl;

	// amplitude damping on |1>: rho' = [[gamma, 0], [0, 1 - gamma]]
	{
		const double gamma = 0.3;
		QC::DensityMatrix<> dm(1);
		QC::Gates::PauliXGate<> x;
		dm.ApplyGate(x, 0); // |0> -> |1>
		dm.ApplyAmplitudeDamping(0, gamma);

		const auto& rho = dm.getDensityMatrix();
		if (!approxEqual(rho(0, 0), std::complex<double>(gamma, 0), 1E-9) ||
			!approxEqual(rho(1, 1), std::complex<double>(1. - gamma, 0), 1E-9) ||
			!approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
			!approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
		{
			std::cout << "Amplitude damping result is wrong:\n" << rho << std::endl;
			return false;
		}
	}

	// dephasing of |+> suppresses the off diagonal coherences while keeping the populations
	{
		QC::DensityMatrix<> dm(1);
		QC::Gates::HadamardGate<> h;
		dm.ApplyGate(h, 0); // |+>

		dm.ApplyPhaseDamping(0, 1.0); // full dephasing

		const auto& rho = dm.getDensityMatrix();
		if (!approxEqual(rho(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
			!approxEqual(rho(1, 1), std::complex<double>(0.5, 0), 1E-9) ||
			!approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
			!approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
		{
			std::cout << "Full dephasing result is wrong:\n" << rho << std::endl;
			return false;
		}

		if (approxEqual(dm.Purity(), 1., 1E-9))
		{
			std::cout << "Dephased state should be mixed but purity is 1" << std::endl;
			return false;
		}
	}

	// depolarizing, bit flip and phase flip must all preserve the trace
	{
		std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> oneQubitGates;
		DM_FillOneQubitGates(oneQubitGates);

		std::uniform_int_distribution gateDistr(0, static_cast<int>(oneQubitGates.size()) - 1);

		for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
		{
			std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

			QC::DensityMatrix<> dm(nrQubits);

			for (int i = 0; i < 20; ++i)
			{
				const int q = qubitDistr(gen);
				dm.ApplyGate(*oneQubitGates[gateDistr(gen)], q);
				dm.ApplyDepolarizingNoise(q, 0.1);
				dm.ApplyBitFlipNoise(q, 0.05);
				dm.ApplyPhaseFlipNoise(q, 0.05);

				if (!approxEqual(dm.Trace().real(), 1., 1E-9) || !approxEqual(dm.Trace().imag(), 0., 1E-9))
				{
					std::cout << "A channel did not preserve the trace: " << dm.Trace() << std::endl;
					return false;
				}

				if (!dm.IsHermitian(1E-9))
				{
					std::cout << "A channel broke Hermiticity" << std::endl;
					return false;
				}
			}
		}
	}

	// reset must always bring a qubit to |0>
	{
		QC::DensityMatrix<> dm(1);
		QC::Gates::HadamardGate<> h;
		dm.ApplyGate(h, 0);
		dm.ApplyReset(0);

		if (!approxEqual(dm.GetQubitProbability(0), 0., 1E-9))
		{
			std::cout << "Reset did not bring the qubit to |0>: p1 = " << dm.GetQubitProbability(0) << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

static bool DensityMatrixGenericChannelTest()
{
	std::cout << "\nDensity matrix simulator - generic Kraus channels against explicit dense evolution" << std::endl;

	// One-qubit amplitude damping embedded in the middle of a three-qubit mixed state.
	{
		QC::DensityMatrix<> dm(3);
		const Eigen::MatrixXcd before = DM_DeterministicMixedState(3);
		dm.setDensityMatrix(before);

		constexpr double gamma = 0.29;
		Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(2, 2);
		k0(0, 0) = 1.;
		k0(1, 1) = std::sqrt(1. - gamma);
		Eigen::MatrixXcd k1 = Eigen::MatrixXcd::Zero(2, 2);
		k1(0, 1) = std::sqrt(gamma);

		const Eigen::MatrixXcd fullK0 = DM_EmbedOperator(k0, { 1 }, 3);
		const Eigen::MatrixXcd fullK1 = DM_EmbedOperator(k1, { 1 }, 3);
		const Eigen::MatrixXcd expected = fullK0 * before * fullK0.adjoint() + fullK1 * before * fullK1.adjoint();
		dm.ApplyChannel({ k0, k1 }, 1);
		if ((dm.getDensityMatrix() - expected).norm() > 1E-10)
		{
			std::cout << "Generic one-qubit Kraus evolution does not match the explicit matrix result" << std::endl;
			return false;
		}
	}

	// A correlated channel on non-adjacent qubits zero and two.
	{
		QC::DensityMatrix<> dm(3);
		const Eigen::MatrixXcd before = DM_DeterministicMixedState(3);
		dm.setDensityMatrix(before);
		QC::Gates::CNOTGate<> cnot;
		constexpr double p = 0.37;
		const Eigen::MatrixXcd k0 = std::sqrt(1. - p) * Eigen::MatrixXcd::Identity(4, 4);
		const Eigen::MatrixXcd k1 = std::sqrt(p) * cnot.getRawOperatorMatrix();
		const Eigen::MatrixXcd fullK0 = DM_EmbedOperator(k0, { 0, 2 }, 3);
		const Eigen::MatrixXcd fullK1 = DM_EmbedOperator(k1, { 0, 2 }, 3);
		const Eigen::MatrixXcd expected = fullK0 * before * fullK0.adjoint() + fullK1 * before * fullK1.adjoint();

		dm.ApplyChannel({ k0, k1 }, 0, 2);
		if ((dm.getDensityMatrix() - expected).norm() > 1E-10 ||
			!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-10) || !dm.IsHermitian(1E-10))
		{
			std::cout << "Generic non-adjacent two-qubit Kraus evolution does not match the explicit matrix result" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;
	return true;
}

static bool DensityMatrixInitializationAndValidationTest()
{
	std::cout << "\nDensity matrix simulator - initialization and invalid input handling" << std::endl;

	if (!DM_ExpectInvalidArgument([] { QC::DensityMatrix<> invalid(0); }, "Zero-qubit construction") ||
		!DM_ExpectInvalidArgument([] { QC::DensityMatrix<> invalid(std::numeric_limits<size_t>::digits); }, "Oversized-qubit construction"))
		return false;
	try
	{
		QC::DensityMatrix<> overflow(std::numeric_limits<size_t>::digits / 2);
		std::cout << "Overflowing density-matrix storage construction did not throw" << std::endl;
		return false;
	}
	catch (const std::length_error&)
	{
	}

	QC::DensityMatrix<> dm(2);
	if (!dm.GetMultithreading())
	{
		std::cout << "Default multithreading state should be true" << std::endl;
		return false;
	}
	dm.SetMultithreading(false);
	if (dm.GetMultithreading())
	{
		std::cout << "SetMultithreading(false) failed" << std::endl;
		return false;
	}
	dm.SetMultithreading(true);
	if (!dm.GetMultithreading())
	{
		std::cout << "SetMultithreading(true) failed" << std::endl;
		return false;
	}

	QC::Gates::PauliXGate<> x;
	QC::Gates::CNOTGate<> cnot;
	QC::Gates::ToffoliGate<> toffoli;
	const Eigen::MatrixXcd original = dm.getDensityMatrix();

	Eigen::MatrixXcd malformed = Eigen::MatrixXcd::Identity(3, 3);
	QC::Gates::AppliedGate<> malformedGate(malformed, 0);
	Eigen::MatrixXcd tooLarge = Eigen::MatrixXcd::Identity(16, 16);
	QC::Gates::AppliedGate<> tooLargeGate(tooLarge, 0, 1);
	Eigen::MatrixXcd nonFiniteGateMatrix = Eigen::MatrixXcd::Identity(2, 2);
	nonFiniteGateMatrix(0, 0) = std::numeric_limits<double>::quiet_NaN();
	QC::Gates::AppliedGate<> nonFiniteGate(nonFiniteGateMatrix, 0);

	if (!DM_ExpectInvalidArgument([&] { dm.ApplyGate(x, 2); }, "Out-of-range gate qubit") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyGate(cnot, 0, 0); }, "Duplicate two-qubit gate arguments") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyGate(toffoli, 0, 1, 1); }, "Duplicate three-qubit gate arguments") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyGate(malformedGate); }, "Non-power-of-two gate matrix") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyGate(tooLargeGate); }, "Unsupported four-qubit gate") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyGate(nonFiniteGate); }, "Non-finite gate matrix") ||
		!DM_ExpectInvalidArgument([&] { dm.GetQubitProbability(2); }, "Out-of-range probability qubit") ||
		!DM_ExpectInvalidArgument([&] { dm.MeasureQubit(2); }, "Out-of-range measured qubit") ||
		!DM_ExpectInvalidArgument([&] { dm.DephaseMeasure(2); }, "Out-of-range dephased qubit") ||
		!DM_ExpectInvalidArgument([&] { dm.MeasureNoCollapse(std::set<size_t>{ 0, 2 }); }, "Out-of-range subset qubit") ||
		!DM_ExpectInvalidArgument([&] { dm.RepeatedMeasure(1, 0, 1); }, "Reversed measurement range") ||
		!DM_ExpectInvalidArgument([&] { dm.RepeatedMeasureUnordered(0, 2, 1); }, "Out-of-range measurement range"))
		return false;

	if (!DM_ExpectInvalidArgument([&] { dm.ApplyBitFlipNoise(0, -0.1); }, "Negative noise probability") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyAmplitudeDamping(0, 1.1); }, "Noise probability above one") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyPhaseDamping(0, std::numeric_limits<double>::quiet_NaN()); }, "NaN noise probability"))
		return false;

	const Eigen::MatrixXcd i2 = Eigen::MatrixXcd::Identity(2, 2);
	const Eigen::MatrixXcd i4 = Eigen::MatrixXcd::Identity(4, 4);
	const Eigen::MatrixXcd i8 = Eigen::MatrixXcd::Identity(8, 8);
	Eigen::MatrixXcd nonFiniteKraus = i2;
	nonFiniteKraus(1, 1) = std::numeric_limits<double>::infinity();
	if (!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({}, 0); }, "Empty Kraus channel") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({ 0.5 * i2 }, 0); }, "Non trace-preserving Kraus channel") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({ i2, i4 }, 0); }, "Mismatched Kraus dimensions") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({ i4 }, 0, 0); }, "Duplicate two-qubit channel arguments") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({ i2 }, 2); }, "Out-of-range channel qubit") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({ i8 }, 0); }, "Unsupported three-qubit channel") ||
		!DM_ExpectInvalidArgument([&] { dm.ApplyChannel({ nonFiniteKraus }, 0); }, "Non-finite Kraus operator"))
		return false;

	Eigen::VectorXcd wrongVector = Eigen::VectorXcd::Zero(3);
	Eigen::VectorXcd zeroVector = Eigen::VectorXcd::Zero(4);
	Eigen::VectorXcd nonFiniteVector = Eigen::VectorXcd::Ones(4);
	nonFiniteVector(2) = std::numeric_limits<double>::quiet_NaN();
	Eigen::MatrixXcd wrongMatrix = Eigen::MatrixXcd::Identity(2, 2);
	if (!DM_ExpectInvalidArgument([&] { dm.setFromStatevector(wrongVector); }, "Wrong statevector dimension") ||
		!DM_ExpectInvalidArgument([&] { dm.setFromStatevector(zeroVector); }, "Zero statevector") ||
		!DM_ExpectInvalidArgument([&] { dm.setFromStatevector(nonFiniteVector); }, "Non-finite statevector") ||
		!DM_ExpectInvalidArgument([&] { dm.setDensityMatrix(wrongMatrix); }, "Wrong density-matrix dimension") ||
		!DM_ExpectInvalidArgument([&] { dm.ExpectationValue(wrongMatrix); }, "Wrong observable dimension") ||
		!DM_ExpectInvalidArgument([&] { dm.ExpectationValue(std::string("IA")); }, "Invalid Pauli operator") ||
		!DM_ExpectInvalidArgument([&] { dm.setToBasisState(4); }, "Out-of-range basis state") ||
		!DM_ExpectInvalidArgument([&] { dm.setToMixtureOfBasisStates({}); }, "Empty basis-state mixture") ||
		!DM_ExpectInvalidArgument([&] { dm.setToMixtureOfBasisStates({ { 0, 0. }, { 7, 1. } }); }, "Mixture without a valid positive weight") ||
		!DM_ExpectInvalidArgument([&] { dm.setToMixtureOfBasisStates({ { 0, std::numeric_limits<double>::infinity() } }); }, "Non-finite mixture weight"))
		return false;

	if ((dm.getDensityMatrix() - original).norm() > 1E-12)
	{
		std::cout << "An invalid operation changed the density matrix" << std::endl;
		return false;
	}

	// Non-normalized finite vectors are accepted and normalized into a unit-trace pure state.
	Eigen::VectorXcd psi(4);
	psi << 1., std::complex<double>(0., 2.), -1., 0.5;
	dm.setFromStatevector(psi);
	if (!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-12) || !approxEqual(dm.Purity(), 1., 1E-12))
	{
		std::cout << "setFromStatevector did not normalize a finite non-zero vector" << std::endl;
		return false;
	}

	// Duplicate valid states accumulate; non-positive and out-of-range entries are ignored.
	dm.setToMixtureOfBasisStates({ { 1, 1. }, { 1, 2. }, { 3, 1. }, { 0, -4. }, { 7, 9. } });
	if (!approxEqual(dm.getBasisStateProbability(1), 0.75, 1E-12) ||
		!approxEqual(dm.getBasisStateProbability(3), 0.25, 1E-12) ||
		!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-12) || !approxEqual(dm.Purity(), 0.625, 1E-12))
	{
		std::cout << "Basis-state mixture normalization or duplicate accumulation is wrong" << std::endl;
		return false;
	}

	// All closed-interval endpoints are valid for every predefined noise channel.
	for (const double probability : { 0., 1. })
	{
		QC::DensityMatrix<> boundary(1);
		boundary.ApplyGate(QC::Gates::HadamardGate<>(), 0);
		boundary.ApplyBitFlipNoise(0, probability);
		boundary.ApplyPhaseFlipNoise(0, probability);
		boundary.ApplyDepolarizingNoise(0, probability);
		boundary.ApplyAmplitudeDamping(0, probability);
		boundary.ApplyPhaseDamping(0, probability);
		if (!boundary.getDensityMatrix().allFinite() || !boundary.IsHermitian(1E-10) ||
			!approxEqual(boundary.Trace(), std::complex<double>(1., 0.), 1E-10))
		{
			std::cout << "A noise channel failed at probability endpoint " << probability << std::endl;
			return false;
		}
	}

	QC::DensityMatrix<> cleared(1);
	cleared.Clear();
	if (!DM_ExpectDomainError([&] { cleared.MeasureNoCollapse(); }, "Sampling an empty state")) return false;

	QC::DensityMatrix<> negativePopulation(1);
	Eigen::MatrixXcd invalidRho = Eigen::MatrixXcd::Zero(2, 2);
	invalidRho(0, 0) = 1.1;
	invalidRho(1, 1) = -0.1;
	negativePopulation.setDensityMatrix(invalidRho);
	if (!DM_ExpectDomainError([&] { negativePopulation.RepeatedMeasure(2); }, "Sampling a negative population")) return false;

	QC::DensityMatrix<> complexPopulation(1);
	invalidRho.setZero();
	invalidRho(0, 0) = std::complex<double>(1., 0.01);
	complexPopulation.setDensityMatrix(invalidRho);
	if (!DM_ExpectDomainError([&] { complexPopulation.MeasureQubit(0); }, "Measuring a non-real population")) return false;

	std::cout << "Success" << std::endl;
	return true;
}

// checks that single qubit measurement probabilities match the statevector simulator
static bool DensityMatrixMeasurementTest()
{
	std::cout << "\nDensity matrix simulator - measurement probabilities compared against the statevector simulator" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		for (int t = 0; t < 10; ++t)
		{
			QC::DensityMatrix<> dm(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);
			for (int i = 0; i < lim; ++i)
			{
				const int g = gateDistr(gen);
				const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
				if (twoQubitsGate && nrQubits < 2) continue;

				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				dm.ApplyGate(*gates[g], qubit1, qubit2);
				reg.ApplyGate(*gates[g], qubit1, qubit2);
			}

			for (int q = 0; q < nrQubits; ++q)
				if (!approxEqual(dm.GetQubitProbability(q), reg.GetQubitProbability(q), 1E-9))
				{
					std::cout << "Qubit " << q << " probability mismatch for " << nrQubits << " qubits: "
						<< reg.GetQubitProbability(q) << " (statevector) vs " << dm.GetQubitProbability(q) << " (density matrix)" << std::endl;
					return false;
				}
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

// reference Pauli string expectation value <psi|P|psi> computed directly from the statevector,
// P = (x) P_i where character i acts on qubit i (qubit 0 is the least significant bit)
static std::complex<double> StatevectorPauliExpectation(const QC::QubitRegister<>& reg, const std::string& pauliString)
{
	const size_t nrQubits = reg.getNrQubits();
	const size_t nrBasisStates = reg.getNrBasisStates();

	size_t flipMask = 0;
	size_t signMask = 0;
	size_t yCount = 0;
	for (size_t i = 0; i < nrQubits; ++i)
	{
		const size_t bit = 1ULL << i;
		switch (toupper(static_cast<unsigned char>(pauliString[i])))
		{
		case 'I': break;
		case 'X': flipMask |= bit; break;
		case 'Y': flipMask |= bit; signMask |= bit; ++yCount; break;
		case 'Z': signMask |= bit; break;
		}
	}

	std::complex<double> res = 0.;
	for (size_t k = 0; k < nrBasisStates; ++k)
	{
		size_t s = signMask & k;
		bool negative = false;
		while (s) { negative = !negative; s &= s - 1; }

		// <psi|P|psi> = sum_k conj(psi_{k^flip}) * phase(k) * psi_k
		const std::complex<double> term = std::conj(reg.getBasisStateAmplitude(k ^ flipMask)) * reg.getBasisStateAmplitude(k);
		res += negative ? -term : term;
	}

	static const std::complex<double> iPow[4] = { {1., 0.}, {0., 1.}, {-1., 0.}, {0., -1.} };
	res *= iPow[yCount & 3];

	return res;
}

// builds the full 2^N x 2^N Pauli string matrix (Kronecker product, qubit 0 is the least significant bit)
static Eigen::MatrixXcd DM_BuildPauliMatrix(const std::string& pauliString)
{
	Eigen::Matrix2cd I;
	I << 1., 0., 0., 1.;
	Eigen::Matrix2cd X;
	X << 0., 1., 1., 0.;
	Eigen::Matrix2cd Y;
	Y << 0., std::complex<double>(0., -1.), std::complex<double>(0., 1.), 0.;
	Eigen::Matrix2cd Z;
	Z << 1., 0., 0., -1.;

	auto single = [&](char c) -> Eigen::Matrix2cd {
		switch (toupper(static_cast<unsigned char>(c)))
		{
		case 'X': return X;
		case 'Y': return Y;
		case 'Z': return Z;
		default: return I;
		}
		};

	// qubit 0 is the least significant bit, so it must be the rightmost factor of the Kronecker product
	Eigen::MatrixXcd result = single(pauliString[0]);
	for (size_t i = 1; i < pauliString.size(); ++i)
		result = Eigen::kroneckerProduct(single(pauliString[i]), result).eval();

	return result;
}

static bool DensityMatrixPauliExpectationTest()
{
	std::cout << "\nDensity matrix simulator - Pauli string expectation values compared against the statevector simulator" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	const char paulis[4] = { 'I', 'X', 'Y', 'Z' };
	std::uniform_int_distribution pauliDistr(0, 3);

	for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		for (int t = 0; t < 10; ++t)
		{
			QC::DensityMatrix<> dm(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);
			for (int i = 0; i < lim; ++i)
			{
				const int g = gateDistr(gen);
				const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
				if (twoQubitsGate && nrQubits < 2) continue;

				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				dm.ApplyGate(*gates[g], qubit1, qubit2);
				reg.ApplyGate(*gates[g], qubit1, qubit2);
			}

			// try a few random Pauli strings
			for (int s = 0; s < 5; ++s)
			{
				std::string pauliString(nrQubits, 'I');
				for (int q = 0; q < nrQubits; ++q)
					pauliString[q] = paulis[pauliDistr(gen)];

				const std::complex<double> expected = StatevectorPauliExpectation(reg, pauliString);
				const std::complex<double> got = dm.ExpectationValue(pauliString);

				// also cross check against the full matrix operator overload
				const std::complex<double> viaMatrix = dm.ExpectationValue(DM_BuildPauliMatrix(pauliString));

				if (!approxEqual(got, expected, 1E-9) || !approxEqual(got, viaMatrix, 1E-9))
				{
					std::cout << "Pauli string " << pauliString << " expectation mismatch for " << nrQubits << " qubits: "
						<< expected << " (statevector) vs " << got << " (density matrix) vs " << viaMatrix << " (full matrix)" << std::endl;
					return false;
				}
			}
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

// checks the sampling statistics of MeasureQubit match the analytic GetQubitProbability,
// and that DephaseMeasure only kills the coherences between the two measurement sectors
static bool DensityMatrixMeasurementStatisticsTest()
{
	std::cout << "\nDensity matrix simulator - measurement sampling statistics and non selective measurement" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(15, 25);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	const int nrMeasurements = 20000;

	for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		// build a fixed random circuit so both the reference probabilities and the sampled ones use the same state
		std::vector<QC::Gates::AppliedGate<>> circuit;
		const int lim = nrGatesDistr(gen);
		for (int i = 0; i < lim; ++i)
		{
			const int g = gateDistr(gen);
			const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
			if (twoQubitsGate && nrQubits < 2) continue;

			int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
			int qubit2 = qubit1 + 1;
			if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

			circuit.emplace_back(gates[g]->getRawOperatorMatrix(), qubit1, twoQubitsGate ? qubit2 : 0);
		}

		// analytic single qubit probabilities of the prepared state
		QC::DensityMatrix<> ref(nrQubits);
		for (const auto& gate : circuit)
			ref.ApplyGate(gate);

		std::vector<double> expectedP1(nrQubits);
		for (int q = 0; q < nrQubits; ++q)
			expectedP1[q] = ref.GetQubitProbability(q);

		// sample each qubit independently by re-preparing and measuring just that qubit
		for (int q = 0; q < nrQubits; ++q)
		{
			int ones = 0;
			for (int m = 0; m < nrMeasurements; ++m)
			{
				QC::DensityMatrix<> dm(nrQubits);
				dm.SetRandomSeed(0xD3A51000ULL + static_cast<uint64_t>(nrQubits * nrMeasurements + q * nrMeasurements + m));
				for (const auto& gate : circuit)
					dm.ApplyGate(gate);

				if (dm.MeasureQubit(q) == 1) ++ones;
			}

			const double sampled = static_cast<double>(ones) / nrMeasurements;
			if (std::abs(sampled - expectedP1[q]) > 0.05)
			{
				std::cout << "Measurement statistics mismatch for qubit " << q << " with " << nrQubits << " qubits: "
					<< expectedP1[q] << " (analytic) vs " << sampled << " (sampled)" << std::endl;
				return false;
			}
		}
	}

	// DephaseMeasure keeps the populations but removes the coherences between the measured sectors
	{
		QC::DensityMatrix<> dm(1);
		QC::Gates::HadamardGate<> h;
		dm.ApplyGate(h, 0); // |+>: rho = 1/2 [[1,1],[1,1]]

		const double p1Before = dm.GetQubitProbability(0);
		dm.DephaseMeasure(0);
		const auto& rho = dm.getDensityMatrix();

		if (!approxEqual(dm.GetQubitProbability(0), p1Before, 1E-9))
		{
			std::cout << "DephaseMeasure changed the population: " << p1Before << " -> " << dm.GetQubitProbability(0) << std::endl;
			return false;
		}

		if (!approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) || !approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
		{
			std::cout << "DephaseMeasure did not remove the off diagonal coherences:\n" << rho << std::endl;
			return false;
		}

		if (!approxEqual(rho(0, 0), std::complex<double>(0.5, 0), 1E-9) || !approxEqual(rho(1, 1), std::complex<double>(0.5, 0), 1E-9))
		{
			std::cout << "DephaseMeasure changed the diagonal populations:\n" << rho << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

// Verifies the state update, not just measurement frequencies: rho -> P_m rho P_m / p_m.
static bool DensityMatrixCollapseTest()
{
	std::cout << "\nDensity matrix simulator - projective measurement collapse" << std::endl;

	auto checkProjection = [](QC::DensityMatrix<>& dm, size_t measuredQubit, uint64_t seed,
		const char* description, size_t& result) -> bool
	{
		const Eigen::MatrixXcd before = dm.getDensityMatrix();
		dm.SetRandomSeed(seed);
		result = dm.MeasureQubit(measuredQubit);
		const size_t mask = 1ULL << measuredQubit;

		double probability = 0.;
		for (size_t state = 0; state < dm.getNrBasisStates(); ++state)
			if (((state & mask) != 0) == (result != 0)) probability += before(state, state).real();

		Eigen::MatrixXcd expected = before;
		for (size_t row = 0; row < dm.getNrBasisStates(); ++row)
			for (size_t col = 0; col < dm.getNrBasisStates(); ++col)
				if (((row & mask) != 0) != (result != 0) || ((col & mask) != 0) != (result != 0))
					expected(row, col) = 0.;
				else
					expected(row, col) /= probability;

		if ((dm.getDensityMatrix() - expected).norm() > 1E-10 ||
			!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-10) || !dm.IsHermitian(1E-10))
		{
			std::cout << description << " collapse does not equal P_m rho P_m / p_m" << std::endl;
			return false;
		}
		return true;
	};

	QC::Gates::HadamardGate<> h;
	QC::Gates::RyGate<> ry(M_PI / 5.);
	QC::Gates::CNOTGate<> cnot;
	size_t ignoredResult = 0;

	// Pure separable state.
	QC::DensityMatrix<> pure(2);
	pure.ApplyGate(h, 0);
	pure.ApplyGate(ry, 1);
	if (!checkProjection(pure, 1, 0x5011DULL, "Pure-state", ignoredResult)) return false;

	// Mixed state containing both classical populations and coherences.
	QC::DensityMatrix<> mixed(3);
	mixed.ApplyGate(h, 0);
	mixed.ApplyGate(cnot, 2, 0);
	mixed.ApplyGate(h, 1);
	mixed.ApplyAmplitudeDamping(2, 0.23);
	mixed.ApplyPhaseDamping(1, 0.31);
	if (approxEqual(mixed.Purity(), 1., 1E-10) ||
		!checkProjection(mixed, 1, 0xC011A95EULL, "Mixed-state", ignoredResult)) return false;

	// Entangled Bell state: verify both the projected matrix and the subsequent correlation.
	QC::DensityMatrix<> bell(2);
	bell.ApplyGate(h, 0);
	bell.ApplyGate(cnot, 1, 0);
	size_t q0 = 0;
	if (!checkProjection(bell, 0, 0xBE11ULL, "Entangled-state", q0)) return false;
	const size_t q1 = bell.MeasureQubit(1);
	if (q0 != q1)
	{
		std::cout << "Sequential measurements broke an entangled-qubit correlation" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

// checks the general Hermitian-observable ExpectationValue overload against a direct Tr(rho O)
// computation, verifies a few analytic Pauli expectation edge cases, that setFromStatevector and
// setDensityMatrix reproduce the reference matrix element, and that ExpectationValue throws on a
// malformed Pauli string
static bool DensityMatrixExpectationExtrasTest()
{
	std::cout << "\nDensity matrix simulator - general observable expectation, edge cases and error handling" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		for (int t = 0; t < 10; ++t)
		{
			QC::DensityMatrix<> dm(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);
			for (int i = 0; i < lim; ++i)
			{
				const int g = gateDistr(gen);
				const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
				if (twoQubitsGate && nrQubits < 2) continue;

				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				dm.ApplyGate(*gates[g], qubit1, qubit2);
				reg.ApplyGate(*gates[g], qubit1, qubit2);
			}

			// build a random Hermitian observable O = M + M^dagger and compare the fast overload with Tr(rho O)
			const Eigen::Index dim = static_cast<Eigen::Index>(dm.getNrBasisStates());
			Eigen::MatrixXcd M = Eigen::MatrixXcd::Random(dim, dim);
			const Eigen::MatrixXcd O = M + M.adjoint();

			const std::complex<double> got = dm.ExpectationValue(O);
			std::complex<double> expected = 0.;
			for (Eigen::Index row = 0; row < dim; ++row)
				for (Eigen::Index col = 0; col < dim; ++col)
					expected += dm.getDensityMatrix()(row, col) * O(col, row);

			if (!approxEqual(got, expected, 1E-9))
			{
				std::cout << "General observable expectation mismatch for " << nrQubits << " qubits: "
					<< expected << " (Tr(rho O)) vs " << got << " (ExpectationValue)" << std::endl;
				return false;
			}

			// setFromStatevector must reproduce |psi><psi| of the statevector simulator
			QC::DensityMatrix<> dmFromPsi(nrQubits);
			dmFromPsi.setFromStatevector(reg.getRegisterStorage());
			if ((dmFromPsi.getDensityMatrix() - dm.getDensityMatrix()).norm() > 1E-9)
			{
				std::cout << "setFromStatevector does not reproduce the evolved density matrix for " << nrQubits << " qubits" << std::endl;
				return false;
			}

			// setDensityMatrix round trips the state
			QC::DensityMatrix<> dmCopy(nrQubits);
			dmCopy.setDensityMatrix(dm.getDensityMatrix());
			if ((dmCopy.getDensityMatrix() - dm.getDensityMatrix()).norm() > 1E-12)
			{
				std::cout << "setDensityMatrix does not round trip the state for " << nrQubits << " qubits" << std::endl;
				return false;
			}

			// the all identity Pauli string expectation must equal the trace (1)
			const std::string allI(nrQubits, 'I');
			if (!approxEqual(dm.ExpectationValue(allI), std::complex<double>(1., 0.), 1E-9))
			{
				std::cout << "Identity Pauli string expectation is not the trace for " << nrQubits << " qubits: "
					<< dm.ExpectationValue(allI) << std::endl;
				return false;
			}
		}
	}

	// analytic single qubit edge cases
	{
		// <+|X|+> = 1, <+|Z|+> = 0
		QC::DensityMatrix<> dm(1);
		QC::Gates::HadamardGate<> h;
		dm.ApplyGate(h, 0);
		if (!approxEqual(dm.ExpectationValue(std::string("X")), std::complex<double>(1., 0.), 1E-9) ||
			!approxEqual(dm.ExpectationValue(std::string("Z")), std::complex<double>(0., 0.), 1E-9))
		{
			std::cout << "Single qubit |+> Pauli expectations are wrong: <X>=" << dm.ExpectationValue(std::string("X"))
				<< " <Z>=" << dm.ExpectationValue(std::string("Z")) << std::endl;
			return false;
		}
	}

	{
		// <1|Z|1> = -1
		QC::DensityMatrix<> dm(1);
		QC::Gates::PauliXGate<> x;
		dm.ApplyGate(x, 0);
		if (!approxEqual(dm.ExpectationValue(std::string("Z")), std::complex<double>(-1., 0.), 1E-9))
		{
			std::cout << "Single qubit |1> <Z> is not -1: " << dm.ExpectationValue(std::string("Z")) << std::endl;
			return false;
		}
	}

	// error handling: a Pauli string of the wrong length must throw
	{
		QC::DensityMatrix<> dm(2);
		bool threw = false;
		try
		{
			dm.ExpectationValue(std::string("X")); // only 1 char for 2 qubits
		}
		catch (const std::invalid_argument&)
		{
			threw = true;
		}

		if (!threw)
		{
			std::cout << "ExpectationValue did not throw on a Pauli string of the wrong length" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

// samples full-register and contiguous-subregister distributions using all four RepeatedMeasure
// overloads and compares the empirical histograms against the exact diagonal populations
static bool DensityMatrixSamplingTest()
{
	std::cout << "\nDensity matrix simulator - repeated sampling against the diagonal populations" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(15, 25);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	const int nrSamples = 100000;

	for (int nrQubits = 1; nrQubits < 4; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		// prepare a fixed random state (unitary circuit plus some noise so it is mixed)
		QC::DensityMatrix<> dm(nrQubits);
		const int lim = nrGatesDistr(gen);
		for (int i = 0; i < lim; ++i)
		{
			const int g = gateDistr(gen);
			const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
			if (twoQubitsGate && nrQubits < 2) continue;

			int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
			int qubit2 = qubit1 + 1;
			if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

			dm.ApplyGate(*gates[g], qubit1, qubit2);
		}
		dm.ApplyDepolarizingNoise(qubitDistr(gen), 0.3);
		dm.SetRandomSeed(0x5A6D1100ULL + static_cast<uint64_t>(nrQubits));

		// exact populations and a snapshot used to verify sampling does not change the state
		const size_t nrBasisStates = dm.getNrBasisStates();
		std::vector<double> expected(nrBasisStates);
		for (size_t s = 0; s < nrBasisStates; ++s)
			expected[s] = dm.getBasisStateProbability(s);
		const Eigen::MatrixXcd stateBeforeSampling = dm.getDensityMatrix();

		auto checkHistogram = [nrSamples](const auto& counts, const std::vector<double>& probabilities,
			const std::string& description) -> bool
		{
			size_t total = 0;
			for (const auto& measurement : counts)
				total += measurement.second;
			if (total != static_cast<size_t>(nrSamples))
			{
				std::cout << description << " returned " << total << " samples instead of " << nrSamples << std::endl;
				return false;
			}

			for (size_t s = 0; s < probabilities.size(); ++s)
			{
				const auto it = counts.find(s);
				const size_t count = it == counts.end() ? 0 : it->second;
				const double sampled = static_cast<double>(count) / nrSamples;
				if (std::abs(sampled - probabilities[s]) > 0.02)
				{
					std::cout << description << " mismatch for outcome " << s << ": "
						<< probabilities[s] << " (exact) vs " << sampled << " (sampled)" << std::endl;
					return false;
				}
			}

			return true;
		};

		const auto orderedCounts = dm.RepeatedMeasure(nrSamples);
		if (!checkHistogram(orderedCounts, expected, "Ordered full-register sampling")) return false;

		const auto unorderedCounts = dm.RepeatedMeasureUnordered(nrSamples);
		if (!checkHistogram(unorderedCounts, expected, "Unordered full-register sampling")) return false;

		// sample a contiguous subregister and compare against its exact marginal. Starting at qubit one
		// when possible also checks that the returned outcomes are shifted and packed correctly.
		{
			const size_t firstQubit = nrQubits > 1 ? 1 : 0;
			const size_t secondQubit = static_cast<size_t>(nrQubits - 1);
			const size_t firstPartMask = (1ULL << firstQubit) - 1;
			const size_t measuredPartMask = (1ULL << (secondQubit + 1)) - 1 - firstPartMask;
			const size_t subregisterStates = 1ULL << (secondQubit - firstQubit + 1);
			std::vector<double> expectedSubregister(subregisterStates, 0.);
			for (size_t s = 0; s < nrBasisStates; ++s)
				expectedSubregister[(s & measuredPartMask) >> firstQubit] += expected[s];

			const auto orderedSubregisterCounts = dm.RepeatedMeasure(firstQubit, secondQubit, nrSamples);
			if (!checkHistogram(orderedSubregisterCounts, expectedSubregister, "Ordered subregister sampling")) return false;

			const auto unorderedSubregisterCounts = dm.RepeatedMeasureUnordered(firstQubit, secondQubit, nrSamples);
			if (!checkHistogram(unorderedSubregisterCounts, expectedSubregister, "Unordered subregister sampling")) return false;
		}

		if (!dm.RepeatedMeasure(0).empty() || !dm.RepeatedMeasureUnordered(0).empty() ||
			!dm.RepeatedMeasure(0, 0, 0).empty() || !dm.RepeatedMeasureUnordered(0, 0, 0).empty())
		{
			std::cout << "Zero-shot sampling returned a non-empty histogram" << std::endl;
			return false;
		}

		if ((dm.getDensityMatrix() - stateBeforeSampling).norm() > 1E-12)
		{
			std::cout << "Repeated sampling changed the state for " << nrQubits << " qubits" << std::endl;
			return false;
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

static bool DensityMatrixStateSaveRestoreCloneTest()
{
	std::cout << "\nDensity matrix simulator - state save, restore and clone" << std::endl;

	QC::DensityMatrix<> dm(3);
	QC::Gates::HadamardGate<> h;
	QC::Gates::CNOTGate<> cnot;
	dm.ApplyGate(h, 0);
	dm.ApplyGate(cnot, 2, 0);
	dm.ApplyDepolarizingNoise(1, 0.3);

	const Eigen::MatrixXcd savedState = dm.getDensityMatrix();
	dm.SaveState();

	dm.setToBasisState(7);
	const Eigen::MatrixXcd stateAtClone = dm.getDensityMatrix();
	auto clone = dm.Clone();
	if ((clone->getDensityMatrix() - stateAtClone).norm() > 1E-12)
	{
		std::cout << "Density matrix clone did not copy the current state" << std::endl;
		return false;
	}

	// A regular restore keeps the snapshot available for later restores.
	dm.RestoreState();
	if ((dm.getDensityMatrix() - savedState).norm() > 1E-12)
	{
		std::cout << "Density matrix RestoreState did not restore the snapshot" << std::endl;
		return false;
	}
	dm.setToBasisState(7);
	dm.RestoreState();
	if ((dm.getDensityMatrix() - savedState).norm() > 1E-12)
	{
		std::cout << "Density matrix RestoreState consumed the snapshot" << std::endl;
		return false;
	}

	// A destructive restore consumes the original's snapshot.
	dm.setToBasisState(7);
	dm.RestoreStateDestructive();
	if ((dm.getDensityMatrix() - savedState).norm() > 1E-12)
	{
		std::cout << "Density matrix RestoreStateDestructive did not restore the snapshot" << std::endl;
		return false;
	}
	dm.setToBasisState(0);
	const Eigen::MatrixXcd stateWithoutSnapshot = dm.getDensityMatrix();
	dm.RestoreStateDestructive();
	if ((dm.getDensityMatrix() - stateWithoutSnapshot).norm() > 1E-12)
	{
		std::cout << "Density matrix destructive restore did not consume the snapshot" << std::endl;
		return false;
	}

	// The clone owns an independent copy of both the current state and the saved snapshot.
	clone->RestoreState();
	if ((clone->getDensityMatrix() - savedState).norm() > 1E-12)
	{
		std::cout << "Density matrix clone did not copy the saved snapshot" << std::endl;
		return false;
	}
	clone->setToBasisState(7);
	if ((dm.getDensityMatrix() - stateWithoutSnapshot).norm() > 1E-12)
	{
		std::cout << "Changing a density matrix clone changed the original" << std::endl;
		return false;
	}
	clone->RestoreStateDestructive();
	if ((clone->getDensityMatrix() - savedState).norm() > 1E-12)
	{
		std::cout << "Density matrix clone's destructive restore failed" << std::endl;
		return false;
	}

	// A clone is a complete simulator snapshot, including the random-engine position.
	QC::DensityMatrix<> sampler(1);
	sampler.ApplyGate(h, 0);
	sampler.SetRandomSeed(0xC10EULL);
	(void)sampler.MeasureNoCollapse(); // advance the engine before cloning
	auto samplerClone = sampler.Clone();
	if (sampler.RepeatedMeasure(2000) != samplerClone->RepeatedMeasure(2000))
	{
		std::cout << "Density matrix clone did not preserve the random-engine state" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

// smallest eigenvalue of a Hermitian matrix (used to check positive semidefiniteness)
static double DM_SmallestEigenvalue(const Eigen::MatrixXcd& rho)
{
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(rho);
	return solver.eigenvalues().minCoeff();
}
// unitary + noise evolution: unit trace, hermiticity, positive semidefiniteness, purity bounds,
// plus a few analytic limiting cases (full depolarizing, full amplitude damping, reset idempotence)
static bool DensityMatrixInvariantsTest()
{
	std::cout << "\nDensity matrix simulator - physics invariants (trace, hermiticity, positivity, purity)" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	DM_FillOneQubitGates(gates);
	DM_FillTwoQubitGates(gates);

	std::uniform_int_distribution nrGatesDistr(20, 40);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);
	std::uniform_int_distribution noiseDistr(0, 5);
	std::uniform_real_distribution<double> probDistr(0.05, 0.5);

	for (int nrQubits = 1; nrQubits < DM_NR_QUBITS_LIMIT; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, std::max(0, nrQubits - 2));

		for (int t = 0; t < 10; ++t)
		{
			QC::DensityMatrix<> dm(nrQubits);

			const int lim = nrGatesDistr(gen);
			for (int i = 0; i < lim; ++i)
			{
				const int g = gateDistr(gen);
				const bool twoQubitsGate = gates[g]->getQubitsNumber() == 2;
				if (twoQubitsGate && nrQubits < 2) continue;

				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				dm.ApplyGate(*gates[g], qubit1, qubit2);

				// sprinkle in a random noise channel every few gates
				if (i % 4 == 3)
				{
					const int q = qubitDistr(gen);
					const double p = probDistr(gen);
					switch (noiseDistr(gen))
					{
					case 0: dm.ApplyBitFlipNoise(q, p); break;
					case 1: dm.ApplyPhaseFlipNoise(q, p); break;
					case 2: dm.ApplyDepolarizingNoise(q, p); break;
					case 3: dm.ApplyAmplitudeDamping(q, p); break;
					case 4: dm.ApplyPhaseDamping(q, p); break;
					default: dm.ApplyReset(q); break;
					}
				}
			}

			const auto& rho = dm.getDensityMatrix();

			// trace must remain 1 (real)
			if (!approxEqual(dm.Trace().real(), 1., 1E-9) || !approxEqual(dm.Trace().imag(), 0., 1E-9))
			{
				std::cout << "Trace is not 1 after noisy evolution: " << dm.Trace() << std::endl;
				return false;
			}

			// must stay Hermitian
			if (!dm.IsHermitian(1E-9))
			{
				std::cout << "Density matrix is not Hermitian after noisy evolution" << std::endl;
				return false;
			}

			// must stay positive semidefinite (allow a tiny negative numerical slack)
			const double minEig = DM_SmallestEigenvalue(rho);
			if (minEig < -1E-9)
			{
				std::cout << "Density matrix is not positive semidefinite, smallest eigenvalue: " << minEig << std::endl;
				return false;
			}

			// purity must lie within [1/2^N, 1]
			const double purity = dm.Purity();
			const double minPurity = 1. / static_cast<double>(dm.getNrBasisStates());
			if (purity < minPurity - 1E-9 || purity > 1. + 1E-9)
			{
				std::cout << "Purity " << purity << " out of bounds [" << minPurity << ", 1] for " << nrQubits << " qubits" << std::endl;
				return false;
			}
		}
	}

	// analytic limiting cases on a single qubit
	{
		// depolarizing at p = 3/4 sends any state to the maximally mixed I/2
		// (for the {I, X, Y, Z} Kraus parametrization the fixed point is p = 3/4, not p = 1)
		QC::DensityMatrix<> dm(1);
		QC::Gates::HadamardGate<> h;
		dm.ApplyGate(h, 0);
		dm.ApplyDepolarizingNoise(0, 0.75);

		const auto& rho = dm.getDensityMatrix();
		if (!approxEqual(rho(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
			!approxEqual(rho(1, 1), std::complex<double>(0.5, 0), 1E-9) ||
			!approxEqual(rho(0, 1), std::complex<double>(0, 0), 1E-9) ||
			!approxEqual(rho(1, 0), std::complex<double>(0, 0), 1E-9))
		{
			std::cout << "Depolarizing at p = 3/4 did not produce the maximally mixed state:\n" << rho << std::endl;
			return false;
		}

		if (!approxEqual(dm.Purity(), 0.5, 1E-9))
		{
			std::cout << "Maximally mixed single qubit purity is not 1/2: " << dm.Purity() << std::endl;
			return false;
		}
	}

	{
		// full amplitude damping (gamma = 1) sends |1> to |0>
		QC::DensityMatrix<> dm(1);
		QC::Gates::PauliXGate<> x;
		dm.ApplyGate(x, 0); // |1>
		dm.ApplyAmplitudeDamping(0, 1.);

		const auto& rho = dm.getDensityMatrix();
		if (!approxEqual(rho(0, 0), std::complex<double>(1., 0), 1E-9) ||
			!approxEqual(rho(1, 1), std::complex<double>(0., 0), 1E-9))
		{
			std::cout << "Full amplitude damping did not relax |1> to |0>:\n" << rho << std::endl;
			return false;
		}
	}

	{
		// reset is idempotent: applying it twice equals applying it once
		QC::DensityMatrix<> dm1(2);
		QC::DensityMatrix<> dm2(2);
		QC::Gates::HadamardGate<> h;
		QC::Gates::CNOTGate<> cnot;
		for (auto* dm : { &dm1, &dm2 })
		{
			dm->ApplyGate(h, 0);
			dm->ApplyGate(cnot, 1, 0);
		}

		dm1.ApplyReset(0);
		dm2.ApplyReset(0);
		dm2.ApplyReset(0); // second reset must be a no-op

		const auto& r1 = dm1.getDensityMatrix();
		const auto& r2 = dm2.getDensityMatrix();
		if ((r1 - r2).norm() > 1E-9)
		{
			std::cout << "Reset is not idempotent" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

static bool DensityMatrixNewFeaturesTest()
{
	std::cout << "\nDensity matrix simulator - bitvector initialization, partial trace, overlap and fidelity" << std::endl;

	// 1. Bitvector initialization overloads
	{
		QC::DensityMatrix<> dm(2);
		dm.setToBasisState(std::vector<bool>{ true, false }); // |01> state (qubit 0 = 1, qubit 1 = 0) => state 1
		if (!approxEqual(dm.getBasisStateProbability(1), 1.0, 1E-9))
		{
			std::cout << "setToBasisState(vector<bool>) failed" << std::endl;
			return false;
		}

		dm.setToMixtureOfBasisStates(std::vector<std::pair<std::vector<bool>, double>>{
			{ { false, false }, 0.4 },
			{ { true, true }, 0.6 }
		});
		if (!approxEqual(dm.getBasisStateProbability(0), 0.4, 1E-9) ||
			!approxEqual(dm.getBasisStateProbability(3), 0.6, 1E-9))
		{
			std::cout << "setToMixtureOfBasisStates(vector<bool>) failed" << std::endl;
			return false;
		}
	}

	// 2. Partial trace: Bell state (|00> + |11>)/sqrt(2) -> trace out qubit 1 -> reduced state on qubit 0 is I/2
	{
		QC::DensityMatrix<> dm(2);
		QC::Gates::HadamardGate<> h;
		QC::Gates::CNOTGate<> cnot;
		dm.ApplyGate(h, 0);
		dm.ApplyGate(cnot, 1, 0);

		const Eigen::MatrixXcd rho0 = dm.PartialTrace({ 0 }); // keep qubit 0
		if (!approxEqual(rho0(0, 0), std::complex<double>(0.5, 0), 1E-9) ||
			!approxEqual(rho0(1, 1), std::complex<double>(0.5, 0), 1E-9) ||
			!approxEqual(rho0(0, 1), std::complex<double>(0, 0), 1E-9))
		{
			std::cout << "PartialTrace of Bell state failed:\n" << rho0 << std::endl;
			return false;
		}
	}

	// 3. Overlap and Fidelity
	{
		QC::DensityMatrix<> dm1(1);
		QC::DensityMatrix<> dm2(1);
		QC::Gates::HadamardGate<> h;
		dm1.ApplyGate(h, 0); // |+>
		dm2.ApplyGate(h, 0); // |+>

		if (!approxEqual(dm1.HilbertSchmidtOverlap(dm2).real(), 1.0, 1E-9))
		{
			std::cout << "HilbertSchmidtOverlap of identical pure states is not 1" << std::endl;
			return false;
		}

		QC::QubitRegister<> reg(1);
		reg.ApplyGate(h, 0);
		if (!approxEqual(dm1.FidelityWithStatevector(reg.getRegisterStorage()), 1.0, 1E-9))
		{
			std::cout << "FidelityWithStatevector with identical state is not 1" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;
	return true;
}

bool DensityMatrixTests()
{
	std::cout << "\nDensity matrix simulator tests" << std::endl;

	const std::mt19937 savedGenerator = gen;
	gen.seed(0xD35A17U);
	const bool result = DensityMatrixUnitaryTest() && DensityMatrixDenseGateReferenceTest() && DensityMatrixThreeQubitGateTest() && DensityMatrixChannelsTest() &&
		DensityMatrixGenericChannelTest() && DensityMatrixInitializationAndValidationTest() && DensityMatrixMeasurementTest() &&
		DensityMatrixPauliExpectationTest() && DensityMatrixInvariantsTest() && DensityMatrixMeasurementStatisticsTest() &&
		DensityMatrixCollapseTest() && DensityMatrixExpectationExtrasTest() && DensityMatrixSamplingTest() && DensityMatrixStateSaveRestoreCloneTest() &&
		DensityMatrixNewFeaturesTest();
	gen = savedGenerator;
	return result;
}
