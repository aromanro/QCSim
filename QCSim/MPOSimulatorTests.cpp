#include <iostream>
#include <iterator>
#include <memory>
#include <vector>
#include <unordered_map>
#include <utility>
#include <array>
#include <cmath>

#include "Tests.h"

#include "QubitRegister.h"
#include "MPOSimulator.h"
#include "DensityMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>


#define NR_QUBITS_LIMIT_MPO 7

// these mirror the helpers from MPSSimulatorTests.cpp, kept local to avoid coupling the two test files
static void FillOneQubitGatesMPO(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
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

static void FillTwoQubitGatesMPO(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::SwapGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledYGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledZGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledHadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledPhaseShiftGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledUGate<>>(M_PI / 3, M_PI / 7));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRxGate<>>(M_PI / 5));
}

// builds the reference density matrix rho = |psi><psi| from the statevector simulator
static Eigen::MatrixXcd ReferenceDensityMatrix(const QC::QubitRegister<>& reg)
{
	return reg.getDensityMatrix();
}

static bool CompareDensityMatrices(const Eigen::MatrixXcd& rhoRef, const Eigen::MatrixXcd& rhoMPO, int nrQubits, double err = 1E-3)
{
	if (rhoRef.rows() != rhoMPO.rows() || rhoRef.cols() != rhoMPO.cols())
	{
		std::cout << "Density matrix dimensions mismatch for the MPO simulator for " << nrQubits << " qubits" << std::endl;
		return false;
	}

	for (Eigen::Index r = 0; r < rhoRef.rows(); ++r)
		for (Eigen::Index c = 0; c < rhoRef.cols(); ++c)
			if (!approxEqual(rhoRef(r, c), rhoMPO(r, c), err))
			{
				std::cout << "Density matrix element (" << r << ", " << c << ") simulation test failed for the MPO simulator for " << nrQubits << " qubits" << std::endl;
				std::cout << "Reference: " << rhoRef(r, c) << " vs MPO: " << rhoMPO(r, c) << std::endl;

				return false;
			}

	return true;
}

// applies the same one and two qubit gates (on adjacent qubits) to both an MPO simulator and a statevector
// register, comparing the resulting density matrix against |psi><psi| after each gate
static bool OneAndTwoQubitGatesTestMPO()
{
	std::cout << "\nMPO simulator state test for both one and two qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);

			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;

				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				mpo.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);
			}

			const Eigen::MatrixXcd rhoRef = ReferenceDensityMatrix(reg);
			const Eigen::MatrixXcd rhoMPO = mpo.getDensityMatrix();

			if (!CompareDensityMatrices(rhoRef, rhoMPO, nrQubits))
				return false;

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// like the above, but using the higher level MPOSimulator that swaps non adjacent qubits together
static bool NonAdjacentGatesTestMPO()
{
	std::cout << "\nMPO simulator state test with non adjacent two qubit gates" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPOSimulator mpo(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);

			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;

				int qubit1 = qubitDistr(gen);
				int qubit2 = qubit1;
				if (twoQubitsGate)
				{
					// pick a different qubit, possibly non adjacent
					while (qubit2 == qubit1)
						qubit2 = qubitDistr(gen);
				}

				mpo.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);
			}

			const Eigen::MatrixXcd rhoRef = ReferenceDensityMatrix(reg);
			const Eigen::MatrixXcd rhoMPO = mpo.getDensityMatrix();

			if (!CompareDensityMatrices(rhoRef, rhoMPO, nrQubits))
				return false;

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// checks the trace stays 1 and the per qubit probabilities match the statevector ones
static bool TraceAndProbabilitiesTestMPO()
{
	std::cout << "\nMPO simulator trace and qubit probabilities test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			const int lim = nrGatesDistr(gen);

			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;

				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				mpo.ApplyGate(*gates[gate], qubit1, qubit2);
				reg.ApplyGate(*gates[gate], qubit1, qubit2);
			}

			const std::complex<double> trace = mpo.Trace();
			if (!approxEqual(trace, std::complex<double>(1., 0.), 1E-3))
			{
				std::cout << "Trace test failed for the MPO simulator for " << nrQubits << " qubits, trace is " << trace << std::endl;
				return false;
			}

			for (int q = 0; q < nrQubits; ++q)
			{
				const double probMPO = mpo.GetProbability(q, false); // probability for |1>
				const double probReg = reg.GetQubitProbability(q);

				if (!approxEqual(probMPO, probReg, 1E-3))
				{
					std::cout << "Qubit " << q << " probability test failed for the MPO simulator for " << nrQubits << " qubits" << std::endl;
					std::cout << "Reference: " << probReg << " vs MPO: " << probMPO << std::endl;
					return false;
				}
			}

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// checks the measurement statistics roughly match the statevector ones
static bool MeasurementsTestMPO()
{
	std::cout << "\nMPO simulator measurements test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	const int nrMeasurements = 10000;

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		// build a fixed random circuit (adjacent gates) to compare statistics on
		std::vector<QC::Gates::AppliedGate<>> circuit;
		const int lim = 25;
		for (int i = 0; i < lim; ++i)
		{
			const int gate = gateDistr(gen);
			const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
			int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
			int qubit2 = qubit1 + 1;
			if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

			QC::Gates::AppliedGate<> appliedGate(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
			circuit.push_back(std::move(appliedGate));
		}

		std::unordered_map<std::vector<bool>, int> measurementsRegMap;
		std::unordered_map<std::vector<bool>, int> measurementsMPOMap;

		for (int t2 = 0; t2 < nrMeasurements; ++t2)
		{
			QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);
			QC::QubitRegister<> reg(nrQubits);

			for (const auto& gate : circuit)
			{
				mpo.ApplyGate(gate);
				reg.ApplyGate(gate);
			}

			std::vector<bool> measurementsReg(nrQubits);
			std::vector<bool> measurementsMPO(nrQubits);

			for (int q = 0; q < nrQubits; ++q)
			{
				measurementsReg[q] = reg.MeasureQubit(q);
				measurementsMPO[q] = mpo.MeasureQubit(q);
			}
			++measurementsRegMap[measurementsReg];
			++measurementsMPOMap[measurementsMPO];
		}

		std::cout << ".";

		for (const auto& [key, value] : measurementsRegMap)
		{
			const double dif = std::abs(static_cast<double>(measurementsMPOMap[key] - value) / nrMeasurements);
			if (dif > 0.05)
			{
				std::cout << "Measurements test failed for the MPO simulator for " << nrQubits << " qubits" << std::endl;
				std::cout << "Might fail due to the randomness of the measurements\n" << std::endl;
				std::cout << "Difference: " << dif << std::endl;
				return false;
			}
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// builds the reference density matrix for a classical mixture sum_k prob_k |state_k><state_k|
// directly as a diagonal matrix (the probabilities are normalized to sum to 1)
static Eigen::MatrixXcd ReferenceMixtureDensityMatrix(const std::vector<std::pair<size_t, double>>& mixture, int nrQubits)
{
	const Eigen::Index dim = static_cast<Eigen::Index>(1ULL << nrQubits);
	Eigen::MatrixXcd rho = Eigen::MatrixXcd::Zero(dim, dim);

	double total = 0.;
	for (const auto& [state, prob] : mixture)
		if (prob > 0.) total += prob;

	if (total <= 0.) return rho;

	for (const auto& [state, prob] : mixture)
		if (prob > 0.)
			rho(static_cast<Eigen::Index>(state), static_cast<Eigen::Index>(state)) += prob / total;

	return rho;
}

// checks that setToMixtureOfBasisStates produces the expected diagonal density matrix with unit trace
static bool MixtureOfBasisStatesTestMPO()
{
	std::cout << "\nMPO simulator mixture of basis states test" << std::endl;

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		const size_t nrBasisStates = 1ULL << nrQubits;
		std::uniform_int_distribution<size_t> stateDistr(0, nrBasisStates - 1);

		for (int t = 0; t < 20; ++t)
		{
			// build a random mixture (possibly with repeated states, which should be merged)
			std::uniform_int_distribution nrTermsDistr(1, static_cast<int>(nrBasisStates));
			const int nrTerms = nrTermsDistr(gen);

			std::vector<std::pair<size_t, double>> mixture;
			for (int k = 0; k < nrTerms; ++k)
				mixture.emplace_back(stateDistr(gen), dist_ampl(gen) + 1.5); // weights in (0.5, 2.5), always positive

			QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);
			mpo.setToMixtureOfBasisStates(mixture);

			const Eigen::MatrixXcd rhoRef = ReferenceMixtureDensityMatrix(mixture, nrQubits);
			const Eigen::MatrixXcd rhoMPO = mpo.getDensityMatrix();

			if (!CompareDensityMatrices(rhoRef, rhoMPO, nrQubits))
				return false;

			const std::complex<double> trace = mpo.Trace();
			if (!approxEqual(trace, std::complex<double>(1., 0.), 1E-3))
			{
				std::cout << "Trace of the mixture is not 1 for the MPO simulator for " << nrQubits << " qubits, it's " << trace << std::endl;
				return false;
			}

			std::cout << ".";
		}
	}

	// also exercise the std::vector<bool> overload on a simple known mixture
	{
		const int nrQubits = 3;
		QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);

		// 0.25 |000><000| + 0.75 |101><101|
		std::vector<std::pair<std::vector<bool>, double>> mixture;
		mixture.emplace_back(std::vector<bool>{false, false, false}, 1.0);
		mixture.emplace_back(std::vector<bool>{true, false, true}, 3.0);
		mpo.setToMixtureOfBasisStates(mixture);

		const Eigen::MatrixXcd rhoMPO = mpo.getDensityMatrix();
		Eigen::MatrixXcd rhoRef = Eigen::MatrixXcd::Zero(8, 8);
		rhoRef(0, 0) = 0.25; // 000
		rhoRef(5, 5) = 0.75; // 101 = bit0 + bit2 = 1 + 4 = 5

		if (!CompareDensityMatrices(rhoRef, rhoMPO, nrQubits))
			return false;

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// checks that evolving a mixture under a circuit equals the probability weighted average of evolving
// each pure component: U (sum_k p_k |s_k><s_k|) U^dagger = sum_k p_k (U|s_k>)(U|s_k>)^dagger
static bool MixtureEvolutionTestMPO()
{
	std::cout << "\nMPO simulator mixture evolution test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		const size_t nrBasisStates = 1ULL << nrQubits;
		std::uniform_int_distribution<size_t> stateDistr(0, nrBasisStates - 1);
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 10; ++t)
		{
			// random normalized mixture of distinct basis states
			std::uniform_int_distribution nrTermsDistr(1, static_cast<int>(nrBasisStates));
			const int nrTerms = nrTermsDistr(gen);

			std::vector<std::pair<size_t, double>> mixture;
			double total = 0.;
			for (int k = 0; k < nrTerms; ++k)
			{
				const double w = dist_ampl(gen) + 1.5;
				mixture.emplace_back(stateDistr(gen), w);
				total += w;
			}

			// build a fixed random circuit (adjacent gates) applied to all simulators
			const int lim = nrGatesDistr(gen);
			std::vector<QC::Gates::AppliedGate<>> circuit;
			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = qubit1 + 1;
				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				circuit.emplace_back(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
			}

			// evolve the mixture in the MPO simulator
			QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);
			mpo.setToMixtureOfBasisStates(mixture);
			for (const auto& gate : circuit)
				mpo.ApplyGate(gate);
			const Eigen::MatrixXcd rhoMPO = mpo.getDensityMatrix();

			// build the reference as the probability weighted average of the evolved pure components
			const Eigen::Index dim = static_cast<Eigen::Index>(nrBasisStates);
			Eigen::MatrixXcd rhoRef = Eigen::MatrixXcd::Zero(dim, dim);
			for (const auto& [state, w] : mixture)
			{
				QC::QubitRegister<> reg(nrQubits);
				reg.setToBasisState(state);
				for (const auto& gate : circuit)
					reg.ApplyGate(gate);

				const auto& psi = reg.getRegisterStorage();
				rhoRef += (w / total) * (psi * psi.adjoint());
			}

			if (!CompareDensityMatrices(rhoRef, rhoMPO, nrQubits))
				return false;

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// checks a non-unitary multi-Kraus channel against dense density matrix evolution
static bool KrausOperatorsTestMPO()
{
	std::cout << "\nMPO simulator Kraus operators test" << std::endl;

	constexpr int nrQubits = 2;
	constexpr double gamma = 0.37;

	Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(2, 2);
	k0(0, 0) = 1.;
	k0(1, 1) = std::sqrt(1. - gamma);

	Eigen::MatrixXcd k1 = Eigen::MatrixXcd::Zero(2, 2);
	k1(0, 1) = std::sqrt(gamma);

	QC::TensorNetworks::MPOSimulator mpo(nrQubits);
	QC::QubitRegister<> reg(nrQubits);

	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;
	mpo.ApplyGate(hGate, 0);
	mpo.ApplyGate(cnotGate, 1, 0);
	reg.ApplyGate(hGate, 0);
	reg.ApplyGate(cnotGate, 1, 0);

	std::vector<Eigen::MatrixXcd> kraus{ k0, k1 };
	mpo.ApplyKrausOperators(kraus, 0);

	const Eigen::MatrixXcd rho = ReferenceDensityMatrix(reg);
	const QC::Gates::SingleQubitGate<> k0Gate(k0);
	const QC::Gates::SingleQubitGate<> k1Gate(k1);
	const Eigen::MatrixXcd k0Full = k0Gate.getOperatorMatrix(nrQubits, 0);
	const Eigen::MatrixXcd k1Full = k1Gate.getOperatorMatrix(nrQubits, 0);
	const Eigen::MatrixXcd rhoRef = k0Full * rho * k0Full.adjoint() + k1Full * rho * k1Full.adjoint();

	if (!CompareDensityMatrices(rhoRef, mpo.getDensityMatrix(), nrQubits))
		return false;

	const std::complex<double> trace = mpo.Trace();
	if (!approxEqual(trace, std::complex<double>(1., 0.), 1E-3))
	{
		std::cout << "Trace after amplitude damping is not 1 for the MPO simulator, it's " << trace << std::endl;
		return false;
	}

	std::cout << ".\nSuccess" << std::endl;

	return true;
}

// builds a random circuit of adjacent one and two qubit gates (two qubit gates act on qubit and qubit+1,
// in either order), returned as a list of AppliedGate so it can be replayed on several simulators
static std::vector<QC::Gates::AppliedGate<>> BuildRandomAdjacentCircuitMPO(const std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates, int nrQubits, int nrGates)
{
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);
	std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
	std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

	std::vector<QC::Gates::AppliedGate<>> circuit;
	circuit.reserve(nrGates);

	for (int i = 0; i < nrGates; ++i)
	{
		const int gate = gateDistr(gen);
		const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
		int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
		int qubit2 = qubit1 + 1;
		if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

		circuit.emplace_back(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
	}

	return circuit;
}

// builds a random circuit where two qubit gates may act on ANY two distinct qubits (not necessarily
// adjacent), in either order. Only usable with simulators that remap/swap qubits internally (the
// MPOSimulator decorator), not with MPOSimulatorImpl which requires adjacent two qubit gates.
static std::vector<QC::Gates::AppliedGate<>> BuildRandomCircuitMPO(const std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates, int nrQubits, int nrGates)
{
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);
	std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

	std::vector<QC::Gates::AppliedGate<>> circuit;
	circuit.reserve(nrGates);

	for (int i = 0; i < nrGates; ++i)
	{
		const int gate = gateDistr(gen);
		const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;

		int qubit1 = qubitDistr(gen);
		int qubit2 = qubit1;
		if (twoQubitsGate)
		{
			do { qubit2 = qubitDistr(gen); } while (qubit2 == qubit1);
		}

		circuit.emplace_back(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
	}

	return circuit;
}

// limits are set so they never actually drop anything (sz == szm). This exercises the limitSize and
// limitEntanglement branches in ApplyTwoQubitGate that the other tests never touch.
static bool CompressionLosslessTestMPO()
{
	std::cout << "\nMPO simulator lossless compression test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);

	// for up to 6 qubits the operator bond dimension is at most 4^3 = 64, so this limit never truncates
	const Eigen::Index largeChi = 1 << 14;

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		for (int t = 0; t < 10; ++t)
		{
			const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, nrGatesDistr(gen));

			QC::QubitRegister<> reg(nrQubits);
			for (const auto& gate : circuit)
				reg.ApplyGate(gate);
			const Eigen::MatrixXcd rhoRef = ReferenceDensityMatrix(reg);

			// (a) a high bond dimension limit: truncation code runs but keeps every singular value
			{
				QC::TensorNetworks::MPOSimulator mpo(nrQubits);
				mpo.setLimitBondDimension(largeChi);
				for (const auto& gate : circuit)
					mpo.ApplyGate(gate);

				if (!CompareDensityMatrices(rhoRef, mpo.getDensityMatrix(), nrQubits))
					return false;

				for (const auto bondDim : mpo.getBondDimensions())
					if (bondDim > largeChi)
					{
						std::cout << "Bond dimension " << bondDim << " exceeds the limit " << largeChi << " for " << nrQubits << " qubits" << std::endl;
						return false;
					}
			}

			// (b) a tiny entanglement threshold: drops only numerically negligible singular values
			{
				QC::TensorNetworks::MPOSimulator mpo(nrQubits);
				mpo.setLimitEntanglement(1E-12);
				for (const auto& gate : circuit)
					mpo.ApplyGate(gate);

				if (!CompareDensityMatrices(rhoRef, mpo.getDensityMatrix(), nrQubits))
					return false;
			}

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// With an aggressive bond dimension limit the represented operator is only an approximation:
// operator-space MPO truncation does not preserve trace/Hermiticity/positivity (see the simulator's
// own documentation), so exact values are NOT asserted here (those are covered by the lossless test).
// Instead this locks the guaranteed invariants: the bond dimension never exceeds the limit and the
// reconstructed density matrix stays finite (no division-by-(near)-zero blow up).
static bool CompressionTruncationTestMPO()
{
	std::cout << "\nMPO simulator bond dimension cap test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	for (int nrQubits = 4; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		for (int t = 0; t < 10; ++t)
		{
			const Eigen::Index chi = 1 + (t % 4); // 1..4, smaller than the untruncated bond (up to 16+)

			QC::TensorNetworks::MPOSimulator mpo(nrQubits);
			mpo.setLimitBondDimension(chi);

			const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, 60);
			for (const auto& gate : circuit)
				mpo.ApplyGate(gate);

			for (const auto bondDim : mpo.getBondDimensions())
				if (bondDim > chi)
				{
					std::cout << "Bond dimension " << bondDim << " exceeds the limit " << chi << " for " << nrQubits << " qubits" << std::endl;
					return false;
				}

			const Eigen::MatrixXcd rho = mpo.getDensityMatrix();
			for (Eigen::Index r = 0; r < rho.rows(); ++r)
				for (Eigen::Index c = 0; c < rho.cols(); ++c)
					if (!std::isfinite(rho(r, c).real()) || !std::isfinite(rho(r, c).imag()))
					{
						std::cout << "Non finite density matrix element after truncation at (" << r << ", " << c << ") for " << nrQubits << " qubits" << std::endl;
						return false;
					}

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// rho -> A rho A^dagger / Tr(A rho A^dagger) for a non-unitary local operator A (amplitude damping K0)
static bool ApplyOperatorAndNormalizeTestMPO()
{
	std::cout << "\nMPO simulator ApplyOperatorAndNormalize test" << std::endl;

	constexpr int nrQubits = 3;
	constexpr double gamma = 0.41;

	Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(2, 2);
	k0(0, 0) = 1.;
	k0(1, 1) = std::sqrt(1. - gamma);

	const QC::Gates::SingleQubitGate<> k0Gate(k0);

	for (int qubit = 0; qubit < nrQubits; ++qubit)
	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		QC::QubitRegister<> reg(nrQubits);

		// some non trivial state with population on every qubit
		QC::Gates::HadamardGate<> hGate;
		QC::Gates::CNOTGate<> cnotGate;
		QC::Gates::RyGate<> ryGate(M_PI / 3);

		std::vector<QC::Gates::AppliedGate<>> circuit;
		circuit.emplace_back(hGate.getRawOperatorMatrix(), 0);
		circuit.emplace_back(cnotGate.getRawOperatorMatrix(), 1, 0);
		circuit.emplace_back(ryGate.getRawOperatorMatrix(), nrQubits - 1);
		circuit.emplace_back(cnotGate.getRawOperatorMatrix(), nrQubits - 1, nrQubits - 2);

		for (const auto& gate : circuit)
		{
			mpo.ApplyGate(gate);
			reg.ApplyGate(gate);
		}

		mpo.ApplyOperatorAndNormalize(k0Gate, qubit);

		// reference: dense A rho A^dagger normalized by its trace
		const Eigen::MatrixXcd rho = ReferenceDensityMatrix(reg);
		const Eigen::MatrixXcd k0Full = k0Gate.getOperatorMatrix(nrQubits, qubit);
		Eigen::MatrixXcd rhoRef = k0Full * rho * k0Full.adjoint();
		rhoRef /= rhoRef.trace();

		if (!CompareDensityMatrices(rhoRef, mpo.getDensityMatrix(), nrQubits))
			return false;

		const std::complex<double> trace = mpo.Trace();
		if (!approxEqual(trace, std::complex<double>(1., 0.), 1E-3))
		{
			std::cout << "Trace after ApplyOperatorAndNormalize on qubit " << qubit << " is not 1, it's " << trace << std::endl;
			return false;
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// a two qubit Kraus channel rho -> (1-p) rho + p U rho U^dagger, with U a (non symmetric) CNOT.
// exercises the two qubit ApplyKrausOperators path, including the decorator swaps for non adjacent qubits
static bool TwoQubitKrausOperatorsTestMPO()
{
	std::cout << "\nMPO simulator two qubit Kraus operators test" << std::endl;

	constexpr double p = 0.3;

	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;
	QC::Gates::RyGate<> ryGate(M_PI / 5);

	// {nrQubits, channelQubit, channelControl}, including a non adjacent case to drive the swaps
	const std::vector<std::array<int, 3>> cases{ {2, 1, 0}, {3, 2, 0} };

	for (const auto& c : cases)
	{
		const int nrQubits = c[0];
		const int channelQubit = c[1];
		const int channelControl = c[2];

		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		QC::QubitRegister<> reg(nrQubits);

		// prepare an entangled, non trivial state
		std::vector<QC::Gates::AppliedGate<>> prep;
		prep.emplace_back(hGate.getRawOperatorMatrix(), 0);
		prep.emplace_back(cnotGate.getRawOperatorMatrix(), 1, 0);
		prep.emplace_back(ryGate.getRawOperatorMatrix(), nrQubits - 1);
		for (const auto& gate : prep)
		{
			mpo.ApplyGate(gate);
			reg.ApplyGate(gate);
		}

		const Eigen::MatrixXcd uRaw = cnotGate.getRawOperatorMatrix();
		const Eigen::MatrixXcd k0 = std::sqrt(1. - p) * Eigen::MatrixXcd::Identity(4, 4);
		const Eigen::MatrixXcd k1 = std::sqrt(p) * uRaw;

		const std::vector<Eigen::MatrixXcd> kraus{ k0, k1 };
		mpo.ApplyKrausOperators(kraus, channelQubit, channelControl);

		// reference via the dense full operators, using the SAME gate-to-qubits convention as the simulator
		const Eigen::MatrixXcd rho = ReferenceDensityMatrix(reg);
		const QC::Gates::TwoQubitsGate<> k0Gate(k0);
		const QC::Gates::TwoQubitsGate<> k1Gate(k1);
		const Eigen::MatrixXcd k0Full = k0Gate.getOperatorMatrix(nrQubits, channelQubit, channelControl);
		const Eigen::MatrixXcd k1Full = k1Gate.getOperatorMatrix(nrQubits, channelQubit, channelControl);
		const Eigen::MatrixXcd rhoRef = k0Full * rho * k0Full.adjoint() + k1Full * rho * k1Full.adjoint();

		if (!CompareDensityMatrices(rhoRef, mpo.getDensityMatrix(), nrQubits))
			return false;

		const std::complex<double> trace = mpo.Trace();
		if (!approxEqual(trace, std::complex<double>(1., 0.), 1E-3))
		{
			std::cout << "Trace after the two qubit Kraus channel is not 1 for " << nrQubits << " qubits, it's " << trace << std::endl;
			return false;
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// getState/setState must snapshot and restore the full state: tensors, bonds and, for the decorator,
// the logical->physical qubit map produced by the swaps of non adjacent gates
static bool StateSaveRestoreTestMPO()
{
	std::cout << "\nMPO simulator state save/restore test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	for (int nrQubits = 3; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

		for (int t = 0; t < 10; ++t)
		{
			QC::TensorNetworks::MPOSimulator mpo(nrQubits);

			// state A: a random circuit including non adjacent gates, so the qubit map gets permuted
			for (int i = 0; i < 20; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = qubitDistr(gen);
				int qubit2 = qubit1;
				if (twoQubitsGate)
					while (qubit2 == qubit1) qubit2 = qubitDistr(gen);

				mpo.ApplyGate(*gates[gate], qubit1, qubit2);
			}

			const Eigen::MatrixXcd rhoA = mpo.getDensityMatrix();
			auto saved = mpo.getState();

			// evolve further into a (almost surely different) state B
			for (int i = 0; i < 20; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = qubitDistr(gen);
				int qubit2 = qubit1;
				if (twoQubitsGate)
					while (qubit2 == qubit1) qubit2 = qubitDistr(gen);

				mpo.ApplyGate(*gates[gate], qubit1, qubit2);
			}

			// restore the snapshot and check we are back at A
			mpo.setState(saved);

			if (!CompareDensityMatrices(rhoA, mpo.getDensityMatrix(), nrQubits))
				return false;

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

static bool TrimTestMPO()
{
	std::cout << "\nMPO simulator Trim test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);
	std::uniform_int_distribution gateDistr(0, static_cast<int>(gates.size()) - 1);

	// a no-op two qubit gate: applying it with truncation enabled on a bond does exactly what Trim
	// does on that bond (contract the two neighbour sites, apply no gate, SVD with truncation)
	const QC::Gates::TwoQubitsGate<> identityTwoQubitGate(Eigen::MatrixXcd::Identity(4, 4));

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		std::uniform_int_distribution qubitDistr2(0, nrQubits - 2);

		for (int t = 0; t < 5; ++t)
		{
			// build a random circuit so it can be replayed identically on several simulators
			std::vector<QC::Gates::AppliedGate<>> circuit;
			const int lim = nrGatesDistr(gen);
			for (int i = 0; i < lim; ++i)
			{
				const int gate = gateDistr(gen);
				const bool twoQubitsGate = gates[gate]->getQubitsNumber() == 2;
				int qubit1 = twoQubitsGate ? qubitDistr2(gen) : qubitDistr(gen);
				int qubit2 = twoQubitsGate ? qubit1 + 1 : qubit1;

				if (twoQubitsGate && dist_bool(gen)) std::swap(qubit1, qubit2);

				circuit.emplace_back(gates[gate]->getRawOperatorMatrix(), qubit1, qubit2);
			}

			const int chi = 1 + (t % 4); // trim down to a bond dimension between 1 and 4

			// reference: apply the circuit without any limit (exact MPO), then lower the limit and
			// truncate with the no-op two qubit gate exactly the bonds Trim would touch, that is
			// only the ones that exceed the (newly lowered) limit. Trimming a bond that is already
			// within the limit is a state preserving re-gauging that would still change the result
			// of subsequent truncations, so the reference must skip those bonds just like Trim does.
			QC::TensorNetworks::MPOSimulatorImpl mpoRef(nrQubits);
			for (const auto& g : circuit) mpoRef.ApplyGate(g);
			mpoRef.setLimitBondDimension(chi);
			const auto refBondDims = mpoRef.getBondDimensions();
			for (int q = 0; q < nrQubits - 1; ++q)
				if (refBondDims[q] > chi)
					mpoRef.ApplyGate(identityTwoQubitGate, q, q + 1);

			// the one under test: apply the same circuit without any limit, then lower the limit and Trim
			QC::TensorNetworks::MPOSimulatorImpl mpoTrim(nrQubits);
			for (const auto& g : circuit) mpoTrim.ApplyGate(g);
			mpoTrim.setLimitBondDimension(chi);
			mpoTrim.Trim();

			// after trimming, every bond dimension must be within the (newly lowered) limit
			for (const auto bd : mpoTrim.getBondDimensions())
				if (bd > chi)
				{
					std::cout << "Trim did not reduce a bond dimension to the limit (" << bd << " > " << chi << ") for " << nrQubits << " qubits" << std::endl;
					return false;
				}

			// Trim must produce the same density matrix as the equivalent no-op two qubit gate truncations
			if (!CompareDensityMatrices(mpoRef.getDensityMatrix(), mpoTrim.getDensityMatrix(), nrQubits))
			{
				std::cout << "Trim density matrix differs from the reference truncation for " << nrQubits << " qubits" << std::endl;
				return false;
			}
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// The MPO simulator and the DensityMatrix simulator both evolve a full density matrix and can handle
// non-unitary operations (Kraus operators, noise channels, classical mixtures). The tests below cross
// check the two dense simulators directly against each other, focusing on the DensityMatrix specific
// functionality rather than only against the pure statevector reference.

// compares the density matrix of both simulators after a random unitary circuit of adjacent gates
static bool UnitaryCircuitVsDensityMatrixMPO()
{
	std::cout << "\nMPO simulator vs DensityMatrix - unitary circuits" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		for (int t = 0; t < 10; ++t)
		{
			const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, nrGatesDistr(gen));

			QC::TensorNetworks::MPOSimulator mpo(nrQubits);
			QC::DensityMatrix<> dm(nrQubits);

			for (const auto& gate : circuit)
			{
				mpo.ApplyGate(gate);
				dm.ApplyGate(gate);
			}

			if (!CompareDensityMatrices(dm.getDensityMatrix(), mpo.getDensityMatrix(), nrQubits))
				return false;

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// applies an amplitude damping single qubit Kraus channel on every qubit in both simulators (after a
// non trivial entangling preparation) and checks the resulting density matrices agree
static bool SingleQubitKrausVsDensityMatrixMPO()
{
	std::cout << "\nMPO simulator vs DensityMatrix - single qubit Kraus channel" << std::endl;

	constexpr double gamma = 0.37;

	Eigen::MatrixXcd k0 = Eigen::MatrixXcd::Zero(2, 2);
	k0(0, 0) = 1.;
	k0(1, 1) = std::sqrt(1. - gamma);

	Eigen::MatrixXcd k1 = Eigen::MatrixXcd::Zero(2, 2);
	k1(0, 1) = std::sqrt(gamma);

	const std::vector<Eigen::MatrixXcd> kraus{ k0, k1 };

	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;
	QC::Gates::RyGate<> ryGate(M_PI / 3);

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		QC::DensityMatrix<> dm(nrQubits);

		// entangling preparation
		std::vector<QC::Gates::AppliedGate<>> prep;
		prep.emplace_back(hGate.getRawOperatorMatrix(), 0);
		for (int q = 1; q < nrQubits; ++q)
			prep.emplace_back(cnotGate.getRawOperatorMatrix(), q, q - 1);
		prep.emplace_back(ryGate.getRawOperatorMatrix(), nrQubits - 1);

		for (const auto& gate : prep)
		{
			mpo.ApplyGate(gate);
			dm.ApplyGate(gate);
		}

		// apply the same channel on every qubit
		for (int q = 0; q < nrQubits; ++q)
		{
			mpo.ApplyKrausOperators(kraus, q);
			dm.ApplyChannel(kraus, q);
		}

		if (!CompareDensityMatrices(dm.getDensityMatrix(), mpo.getDensityMatrix(), nrQubits))
			return false;

		// the trace must stay 1 for a trace preserving channel, for both simulators
		if (!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-3) ||
			!approxEqual(mpo.Trace(), std::complex<double>(1., 0.), 1E-3))
		{
			std::cout << "Trace after the single qubit Kraus channel is not 1 for " << nrQubits << " qubits" << std::endl;
			return false;
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// compares the built in DensityMatrix noise channels against the equivalent Kraus operators applied
// on the MPO simulator (bit flip, phase flip, depolarizing, phase damping)
static bool NoiseChannelsVsDensityMatrixMPO()
{
	std::cout << "\nMPO simulator vs DensityMatrix - noise channels" << std::endl;

	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;
	QC::Gates::RyGate<> ryGate(M_PI / 4);

	constexpr double p = 0.2;
	constexpr double gamma = 0.3;

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		for (int q = 0; q < nrQubits; ++q)
		{
			// bit flip / phase flip / depolarizing / phase damping / amplitude damping / reset,
			// one channel per inner iteration
			for (int channel = 0; channel < 6; ++channel)
			{
				QC::TensorNetworks::MPOSimulator mpo(nrQubits);
				QC::DensityMatrix<> dm(nrQubits);

				// same entangling preparation on both
				std::vector<QC::Gates::AppliedGate<>> prep;
				prep.emplace_back(hGate.getRawOperatorMatrix(), 0);
				for (int c = 1; c < nrQubits; ++c)
					prep.emplace_back(cnotGate.getRawOperatorMatrix(), c, c - 1);
				prep.emplace_back(ryGate.getRawOperatorMatrix(), q);

				for (const auto& gate : prep)
				{
					mpo.ApplyGate(gate);
					dm.ApplyGate(gate);
				}

				// use the equivalent built-in noise helpers on both simulators
				switch (channel)
				{
				case 0: // bit flip
					mpo.ApplyBitFlipNoise(q, p);
					dm.ApplyBitFlipNoise(q, p);
					break;
				case 1: // phase flip
					mpo.ApplyPhaseFlipNoise(q, p);
					dm.ApplyPhaseFlipNoise(q, p);
					break;
				case 2: // depolarizing
					mpo.ApplyDepolarizingNoise(q, p);
					dm.ApplyDepolarizingNoise(q, p);
					break;
				case 3: // phase damping
					mpo.ApplyPhaseDamping(q, gamma);
					dm.ApplyPhaseDamping(q, gamma);
					break;
				case 4: // amplitude damping
					mpo.ApplyAmplitudeDamping(q, gamma);
					dm.ApplyAmplitudeDamping(q, gamma);
					break;
				default: // reset to |0>
					mpo.ApplyReset(q);
					dm.ApplyReset(q);
					break;
				}

				if (!CompareDensityMatrices(dm.getDensityMatrix(), mpo.getDensityMatrix(), nrQubits))
				{
					std::cout << "Mismatch for channel " << channel << " on qubit " << q << " for " << nrQubits << " qubits" << std::endl;
					return false;
				}
			}
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// a two qubit Kraus channel rho -> (1-p) rho + p U rho U^dagger applied on both simulators and compared
static bool TwoQubitKrausVsDensityMatrixMPO()
{
	std::cout << "\nMPO simulator vs DensityMatrix - two qubit Kraus channel" << std::endl;

	constexpr double p = 0.3;

	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;
	QC::Gates::RyGate<> ryGate(M_PI / 5);

	// {nrQubits, channelQubit, channelControl}, including a non adjacent case for the MPO swaps
	const std::vector<std::array<int, 3>> cases{ {2, 1, 0}, {3, 2, 0}, {3, 0, 2} };

	for (const auto& c : cases)
	{
		const int nrQubits = c[0];
		const int channelQubit = c[1];
		const int channelControl = c[2];

		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		QC::DensityMatrix<> dm(nrQubits);

		std::vector<QC::Gates::AppliedGate<>> prep;
		prep.emplace_back(hGate.getRawOperatorMatrix(), 0);
		prep.emplace_back(cnotGate.getRawOperatorMatrix(), 1, 0);
		prep.emplace_back(ryGate.getRawOperatorMatrix(), nrQubits - 1);
		for (const auto& gate : prep)
		{
			mpo.ApplyGate(gate);
			dm.ApplyGate(gate);
		}

		const Eigen::MatrixXcd uRaw = cnotGate.getRawOperatorMatrix();
		const Eigen::MatrixXcd k0 = std::sqrt(1. - p) * Eigen::MatrixXcd::Identity(4, 4);
		const Eigen::MatrixXcd k1 = std::sqrt(p) * uRaw;
		const std::vector<Eigen::MatrixXcd> kraus{ k0, k1 };

		mpo.ApplyKrausOperators(kraus, channelQubit, channelControl);
		dm.ApplyChannel(kraus, channelQubit, channelControl);

		if (!CompareDensityMatrices(dm.getDensityMatrix(), mpo.getDensityMatrix(), nrQubits))
			return false;

		if (!approxEqual(dm.Trace(), std::complex<double>(1., 0.), 1E-3) ||
			!approxEqual(mpo.Trace(), std::complex<double>(1., 0.), 1E-3))
		{
			std::cout << "Trace after the two qubit Kraus channel is not 1 for " << nrQubits << " qubits" << std::endl;
			return false;
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// evolves a classical mixture of basis states under a random circuit in both simulators and compares
// the resulting density matrices, exercising the DensityMatrix setToMixtureOfBasisStates handling
static bool MixtureEvolutionVsDensityMatrixMPO()
{
	std::cout << "\nMPO simulator vs DensityMatrix - mixture evolution" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		const size_t nrBasisStates = 1ULL << nrQubits;
		std::uniform_int_distribution<size_t> stateDistr(0, nrBasisStates - 1);

		for (int t = 0; t < 10; ++t)
		{
			// random normalized mixture of basis states
			std::uniform_int_distribution nrTermsDistr(1, static_cast<int>(nrBasisStates));
			const int nrTerms = nrTermsDistr(gen);

			std::vector<std::pair<size_t, double>> mixture;
			for (int k = 0; k < nrTerms; ++k)
			{
				const double w = dist_ampl(gen) + 1.5;
				mixture.emplace_back(stateDistr(gen), w);
			}

			const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, nrGatesDistr(gen));

			// evolve the mixture in the MPO simulator
			QC::TensorNetworks::MPOSimulator mpo(nrQubits);
			mpo.setToMixtureOfBasisStates(mixture);
			for (const auto& gate : circuit)
				mpo.ApplyGate(gate);

			// evolve the same mixture in the DensityMatrix simulator
			QC::DensityMatrix<> dm(nrQubits);
			dm.setToMixtureOfBasisStates(mixture);
			for (const auto& gate : circuit)
				dm.ApplyGate(gate);

			if (!CompareDensityMatrices(dm.getDensityMatrix(), mpo.getDensityMatrix(), nrQubits))
				return false;

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// applies a random circuit in both simulators, then compares Pauli string expectation values
// computed by the MPO ExpectationValue (chain contraction) against the DensityMatrix ones
static bool PauliExpectationVsDensityMatrixMPO()
{
	std::cout << "\nMPO simulator vs DensityMatrix - Pauli string expectation values" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(15, 30);

	const char paulis[4] = { 'I', 'X', 'Y', 'Z' };
	std::uniform_int_distribution pauliDistr(0, 3);

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		for (int t = 0; t < 10; ++t)
		{
			const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, nrGatesDistr(gen));

			QC::TensorNetworks::MPOSimulator mpo(nrQubits);
			QC::DensityMatrix<> dm(nrQubits);

			for (const auto& gate : circuit)
			{
				mpo.ApplyGate(gate);
				dm.ApplyGate(gate);
			}

			for (int s = 0; s < 5; ++s)
			{
				std::string pauliString(nrQubits, 'I');
				for (int q = 0; q < nrQubits; ++q)
					pauliString[q] = paulis[pauliDistr(gen)];

				const std::complex<double> expected = dm.ExpectationValue(pauliString);
				const std::complex<double> got = mpo.ExpectationValue(pauliString);

				if (!approxEqual(got, expected, 1E-9))
				{
					std::cout << "Pauli string " << pauliString << " expectation mismatch for " << nrQubits << " qubits: "
						<< expected << " (density matrix) vs " << got << " (MPO)" << std::endl;
					return false;
				}
			}

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// smallest eigenvalue of a Hermitian matrix, used to check the MPO density matrix stays positive semidefinite
static double MPO_SmallestEigenvalue(const Eigen::MatrixXcd& rho)
{
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(rho);
	return solver.eigenvalues().minCoeff();
}

// physics invariants for the MPO simulator after a random non-adjacent circuit followed by noise:
// the reconstructed density matrix must be Hermitian, positive semidefinite, unit trace, and its
// purity must lie within the valid bounds; a single Kraus operator equal to a unitary must match the
// plain gate application; ReCanonicalize must not change the represented state
static bool InvariantsTestMPO()
{
	std::cout << "\nMPO simulator - physics invariants (trace, hermiticity, positivity, purity, ReCanonicalize)" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(20, 40);
	std::uniform_int_distribution noiseDistr(0, 5);
	std::uniform_real_distribution<double> probDistr(0.05, 0.5);

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);

		for (int t = 0; t < 5; ++t)
		{
			QC::TensorNetworks::MPOSimulator mpo(nrQubits);

			const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, nrGatesDistr(gen));
			for (const auto& gate : circuit)
				mpo.ApplyGate(gate);

			// a few random noise channels
			for (int n = 0; n < 3; ++n)
			{
				const int q = qubitDistr(gen);
				const double p = probDistr(gen);
				switch (noiseDistr(gen))
				{
				case 0: mpo.ApplyBitFlipNoise(q, p); break;
				case 1: mpo.ApplyPhaseFlipNoise(q, p); break;
				case 2: mpo.ApplyDepolarizingNoise(q, p); break;
				case 3: mpo.ApplyAmplitudeDamping(q, p); break;
				case 4: mpo.ApplyPhaseDamping(q, p); break;
				default: mpo.ApplyReset(q); break;
				}
			}

			// keep a copy of the state before ReCanonicalize
			const Eigen::MatrixXcd rhoBefore = mpo.getDensityMatrix();

			// trace must be 1
			const std::complex<double> trace = mpo.Trace();
			if (!approxEqual(trace, std::complex<double>(1., 0.), 1E-3))
			{
				std::cout << "MPO trace is not 1 after noisy evolution: " << trace << std::endl;
				return false;
			}

			// hermiticity
			if ((rhoBefore - rhoBefore.adjoint()).norm() > 1E-3)
			{
				std::cout << "MPO density matrix is not Hermitian after noisy evolution" << std::endl;
				return false;
			}

			// positive semidefiniteness (allow numerical / truncation slack)
			const double minEig = MPO_SmallestEigenvalue(rhoBefore);
			if (minEig < -1E-3)
			{
				std::cout << "MPO density matrix is not positive semidefinite, smallest eigenvalue: " << minEig << std::endl;
				return false;
			}

			// purity bounds [1/2^N, 1]
			const double purity = (rhoBefore * rhoBefore).trace().real();
			const double minPurity = 1. / static_cast<double>(1ULL << nrQubits);
			if (purity < minPurity - 1E-3 || purity > 1. + 1E-3)
			{
				std::cout << "MPO purity " << purity << " out of bounds for " << nrQubits << " qubits" << std::endl;
				return false;
			}

			// ReCanonicalize must preserve the represented state
			mpo.ReCanonicalize();
			if (!CompareDensityMatrices(rhoBefore, mpo.getDensityMatrix(), nrQubits))
			{
				std::cout << "ReCanonicalize changed the represented MPO state" << std::endl;
				return false;
			}

			std::cout << ".";
		}
	}

	// a single Kraus operator equal to a unitary must equal the plain gate application
	{
		constexpr int nrQubits = 3;
		QC::TensorNetworks::MPOSimulator mpoGate(nrQubits);
		QC::TensorNetworks::MPOSimulator mpoKraus(nrQubits);

		QC::Gates::HadamardGate<> h;
		mpoGate.ApplyGate(h, 0);
		mpoKraus.ApplyGate(h, 0);

		QC::Gates::PauliYGate<> y;
		const Eigen::MatrixXcd yMat = y.getRawOperatorMatrix();

		mpoGate.ApplyGate(y, 1);
		mpoKraus.ApplyKrausOperators(std::vector<Eigen::MatrixXcd>{ yMat }, 1); // single unitary Kraus operator

		if (!CompareDensityMatrices(mpoGate.getDensityMatrix(), mpoKraus.getDensityMatrix(), nrQubits))
		{
			std::cout << "Single unitary Kraus operator does not match the plain gate application" << std::endl;
			return false;
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// compares the decorator (MPOSimulator, which remaps logical qubits and inserts swaps) against the
// adjacent-only implementation (MPOSimulatorImpl) on an equivalent adjacent circuit, ensuring the two
// layers produce the same density matrix, trace and qubit probabilities
static bool DecoratorVsImplTestMPO()
{
	std::cout << "\nMPO simulator - decorator vs implementation on adjacent circuits" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	std::uniform_int_distribution nrGatesDistr(25, 50);

	for (int nrQubits = 2; nrQubits < NR_QUBITS_LIMIT_MPO; ++nrQubits)
	{
		for (int t = 0; t < 10; ++t)
		{
			// adjacent-only circuit so both the decorator and the impl accept it directly
			const auto circuit = BuildRandomAdjacentCircuitMPO(gates, nrQubits, nrGatesDistr(gen));

			QC::TensorNetworks::MPOSimulator mpo(nrQubits);
			QC::TensorNetworks::MPOSimulatorImpl impl(nrQubits);

			for (const auto& gate : circuit)
			{
				mpo.ApplyGate(gate);
				impl.ApplyGate(gate);
			}

			if (!CompareDensityMatrices(impl.getDensityMatrix(), mpo.getDensityMatrix(), nrQubits))
			{
				std::cout << "Decorator and implementation disagree on the density matrix for " << nrQubits << " qubits" << std::endl;
				return false;
			}

			if (!approxEqual(mpo.Trace(), impl.Trace(), 1E-3))
			{
				std::cout << "Decorator and implementation disagree on the trace for " << nrQubits << " qubits" << std::endl;
				return false;
			}

			for (int q = 0; q < nrQubits; ++q)
				if (!approxEqual(mpo.GetProbability(q, false), impl.GetProbability(q, false), 1E-3))
				{
					std::cout << "Decorator and implementation disagree on qubit " << q << " probability for " << nrQubits << " qubits" << std::endl;
					return false;
				}

			std::cout << ".";
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// compares the sampled measurement statistics of the MPO simulator against the DensityMatrix single
// qubit probabilities on the same prepared state, and checks the error handling of ExpectationValue,
// ApplyGate and ApplyKrausOperators for out of range / malformed arguments
static bool MeasurementVsDensityMatrixAndThrowsTestMPO()
{
	std::cout << "\nMPO simulator - measurement statistics vs DensityMatrix and error handling" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	const int nrMeasurements = 20000;

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		// fixed random circuit replayed on both simulators
		const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, 25);

		// analytic single qubit probabilities from the DensityMatrix simulator
		QC::DensityMatrix<> dm(nrQubits);
		for (const auto& gate : circuit)
			dm.ApplyGate(gate);

		std::vector<double> expectedP1(nrQubits);
		for (int q = 0; q < nrQubits; ++q)
			expectedP1[q] = dm.GetQubitProbability(q);

		// sample each qubit by re-preparing the MPO state and measuring just that qubit
		for (int q = 0; q < nrQubits; ++q)
		{
			int ones = 0;
			for (int m = 0; m < nrMeasurements; ++m)
			{
				QC::TensorNetworks::MPOSimulator mpo(nrQubits);
				for (const auto& gate : circuit)
					mpo.ApplyGate(gate);

				if (mpo.MeasureQubit(q)) ++ones;
			}

			const double sampled = static_cast<double>(ones) / nrMeasurements;
			if (std::abs(sampled - expectedP1[q]) > 0.05)
			{
				std::cout << "MPO measurement statistics mismatch for qubit " << q << " with " << nrQubits << " qubits: "
					<< expectedP1[q] << " (DensityMatrix) vs " << sampled << " (MPO sampled)" << std::endl;
				return false;
			}
		}

		std::cout << ".";
	}

	// error handling on the implementation layer (adjacency and index validation)
	{
		QC::TensorNetworks::MPOSimulatorImpl impl(3);

		// ExpectationValue with the wrong Pauli string length must throw
		bool threw = false;
		try { impl.ExpectationValue(std::string("XX")); } // 2 chars for 3 qubits
		catch (const std::invalid_argument&) { threw = true; }
		if (!threw)
		{
			std::cout << "MPO ExpectationValue did not throw on a Pauli string of the wrong length" << std::endl;
			return false;
		}

		// a two qubit gate on non-adjacent qubits must throw on the implementation
		threw = false;
		try
		{
			QC::Gates::CNOTGate<> cnot;
			impl.ApplyGate(cnot, 2, 0); // non adjacent
		}
		catch (const std::invalid_argument&) { threw = true; }
		if (!threw)
		{
			std::cout << "MPO implementation did not throw on a non-adjacent two qubit gate" << std::endl;
			return false;
		}

		// Kraus operators of inconsistent size must throw
		threw = false;
		try
		{
			std::vector<Eigen::MatrixXcd> kraus{ Eigen::MatrixXcd::Identity(2, 2), Eigen::MatrixXcd::Identity(4, 4) };
			impl.ApplyKrausOperators(kraus, 0);
		}
		catch (const std::invalid_argument&) { threw = true; }
		if (!threw)
		{
			std::cout << "MPO implementation did not throw on inconsistent Kraus operator sizes" << std::endl;
			return false;
		}
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

bool MPOSimulatorTests()
{
	std::cout << "\nMPO Simulator Tests" << std::endl;
	return OneAndTwoQubitGatesTestMPO() &&
		NonAdjacentGatesTestMPO() &&
		TraceAndProbabilitiesTestMPO() &&
		MixtureOfBasisStatesTestMPO() &&
		MixtureEvolutionTestMPO() &&
		KrausOperatorsTestMPO() &&
		TwoQubitKrausOperatorsTestMPO() &&
		ApplyOperatorAndNormalizeTestMPO() &&
		CompressionLosslessTestMPO() &&
		CompressionTruncationTestMPO() &&
		StateSaveRestoreTestMPO() &&
		MeasurementsTestMPO() &&
		TrimTestMPO() &&
		UnitaryCircuitVsDensityMatrixMPO() &&
		SingleQubitKrausVsDensityMatrixMPO() &&
		NoiseChannelsVsDensityMatrixMPO() &&
		TwoQubitKrausVsDensityMatrixMPO() &&
		MixtureEvolutionVsDensityMatrixMPO() &&
		PauliExpectationVsDensityMatrixMPO() &&
		InvariantsTestMPO() &&
		DecoratorVsImplTestMPO() &&
		MeasurementVsDensityMatrixAndThrowsTestMPO();
}
