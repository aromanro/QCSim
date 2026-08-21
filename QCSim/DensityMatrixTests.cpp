#include <iostream>
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
			const std::complex<double> expected = (dm.getDensityMatrix() * O).trace();

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

bool DensityMatrixTests()
{
	std::cout << "\nDensity matrix simulator tests" << std::endl;

	return DensityMatrixUnitaryTest() && DensityMatrixChannelsTest() && DensityMatrixMeasurementTest() && DensityMatrixPauliExpectationTest() &&
		DensityMatrixInvariantsTest() && DensityMatrixMeasurementStatisticsTest() && DensityMatrixExpectationExtrasTest() && DensityMatrixSamplingTest();
}
