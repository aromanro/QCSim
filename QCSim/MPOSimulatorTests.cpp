#include <iostream>
#include <iterator>
#include <memory>
#include <vector>
#include <unordered_map>
#include <utility>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <set>

#include "Tests.h"

#include "QubitRegister.h"
#include "MPOSimulator.h"
#include "DensityMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>


#define NR_QUBITS_LIMIT_MPO 7

template<class Callable>
static bool MPO_ExpectInvalidArgument(Callable&& callable, const char* description)
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
static bool MPO_ExpectRuntimeError(Callable&& callable, const char* description)
{
	try
	{
		callable();
	}
	catch (const std::runtime_error&)
	{
		return true;
	}
	catch (const std::exception& ex)
	{
		std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
		return false;
	}

	std::cout << description << " did not throw std::runtime_error" << std::endl;
	return false;
}

class MPOWrongDeclaredArityGate final : public QC::Gates::AppliedGate<>
{
public:
	explicit MPOWrongDeclaredArityGate(const Eigen::MatrixXcd& op)
		: QC::Gates::AppliedGate<>(op)
	{
	}

	size_t getQubitsNumber() const override
	{
		return 1;
	}
};

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

static bool NumericalRankStabilityTestMPO()
{
	std::cout << "\nMPO simulator numerical rank stability test" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	struct CircuitStep
	{
		size_t gate;
		Eigen::Index qubit1;
		Eigen::Index qubit2;
	};

	// This circuit creates roundoff-only singular values and then routes a gate
	// across them. Retaining those values makes the Vidal pseudoinverse amplify
	// SVD noise into percent-level density-matrix errors.
	const std::array<CircuitStep, 11> circuit{ {
		{ 2, 0, 0 }, { 24, 2, 4 }, { 22, 0, 1 }, { 18, 0, 3 },
		{ 0, 4, 4 }, { 23, 0, 4 }, { 1, 3, 3 }, { 5, 3, 3 },
		{ 0, 5, 5 }, { 11, 0, 0 }, { 22, 4, 1 }
	} };

	QC::TensorNetworks::MPOSimulator mpo(6);
	QC::QubitRegister<> reg(6);
	for (const auto& step : circuit)
	{
		mpo.ApplyGate(*gates[step.gate], step.qubit1, step.qubit2);
		reg.ApplyGate(*gates[step.gate], step.qubit1, step.qubit2);
	}

	if (!CompareDensityMatrices(ReferenceDensityMatrix(reg), mpo.getDensityMatrix(), 6, 1E-10))
	{
		std::cout << "MPO numerical-rank filtering did not prevent roundoff amplification" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

static bool MeetingPositionCallbackTestMPO()
{
	std::cout << "\nMPO simulator meeting position callback test" << std::endl;

	constexpr int nrQubits = 5;
	QC::TensorNetworks::MPOSimulator mpo(nrQubits);
	QC::QubitRegister<> reg(nrQubits);

	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;
	mpo.ApplyGate(hGate, 0);
	reg.ApplyGate(hGate, 0);

	const auto expectedBondDimensions = mpo.getBondDimensions();
	int callbackCalls = 0;
	bool receivedExpectedBondDimensions = false;
	mpo.SetMeetingPositionCallback(
		[&](const std::vector<Eigen::Index>& bondDimensions)
		{
			++callbackCalls;
			receivedExpectedBondDimensions = bondDimensions == expectedBondDimensions;
			return Eigen::Index{ 2 };
		});

	// Logical qubits 0 and 4 start at the chain ends. Meeting at bond 2
	// moves them to physical positions 2 and 3, respectively.
	mpo.ApplyGate(cnotGate, 4, 0);
	reg.ApplyGate(cnotGate, 4, 0);

	if (callbackCalls != 1 || !receivedExpectedBondDimensions)
	{
		std::cout << "MPO meeting position callback was not called with the current bond dimensions" << std::endl;
		return false;
	}

	const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	const std::vector<Eigen::Index> expectedMap{ 2, 0, 1, 4, 3 };
	if (!state || state->qubitsMap != expectedMap)
	{
		std::cout << "MPO meeting position callback did not route the qubits to the requested bond" << std::endl;
		return false;
	}

	if (!CompareDensityMatrices(ReferenceDensityMatrix(reg), mpo.getDensityMatrix(), nrQubits))
		return false;

	// Once adjacent, applying another two-qubit gate must not invoke the routing callback.
	mpo.ApplyGate(cnotGate, 4, 0);
	if (callbackCalls != 1)
	{
		std::cout << "MPO meeting position callback was called for adjacent qubits" << std::endl;
		return false;
	}

	std::cout << "\nSuccess" << std::endl;
	return true;
}

static bool OptimalMeetingPositionTestMPO()
{
	std::cout << "\nMPO simulator optimal meeting position test" << std::endl;

	constexpr int nrQubits = 5;
	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;

	auto prepareUnevenBonds = [&](QC::TensorNetworks::MPOSimulator& mpo, QC::QubitRegister<>* reg)
	{
		mpo.ApplyGate(hGate, 0);
		mpo.ApplyGate(cnotGate, 1, 0);
		mpo.ApplyGate(hGate, 4);
		mpo.ApplyGate(cnotGate, 3, 4);
		if (reg)
		{
			reg->ApplyGate(hGate, 0);
			reg->ApplyGate(cnotGate, 1, 0);
			reg->ApplyGate(hGate, 4);
			reg->ApplyGate(cnotGate, 3, 4);
		}
	};

	QC::TensorNetworks::MPOSimulator probe(nrQubits);
	prepareUnevenBonds(probe, nullptr);
	const auto bondDimensions = probe.getBondDimensions();
	if (bondDimensions.size() != 4 || bondDimensions[0] <= bondDimensions[1] ||
		bondDimensions[3] <= bondDimensions[1] || bondDimensions[1] != bondDimensions[2])
	{
		std::cout << "MPO optimal meeting position test could not prepare uneven bond dimensions" << std::endl;
		return false;
	}

	const std::vector<Eigen::Index> expectedOptimalMap{ 1, 0, 3, 4, 2 };
	// Product-state heuristic: qubits 0 and 4 straddle the chain middle, so they
	// meet there rather than at the lower qubit. The swap-up case moves qubit 0
	// toward qubit 2 because qubit 2 is closer to the far end.
	const std::vector<Eigen::Index> expectedHeuristicMap{ 1, 0, 3, 4, 2 };
	const std::vector<Eigen::Index> expectedHeuristicSwapUpMap{ 1, 0, 2, 3, 4 };

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		QC::QubitRegister<> reg(nrQubits);
		prepareUnevenBonds(mpo, &reg);

		mpo.ApplyGate(cnotGate, 4, 0);
		reg.ApplyGate(cnotGate, 4, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
		if (!state || state->qubitsMap != expectedOptimalMap)
		{
			std::cout << "MPO default optimal meeting position did not meet at the cheapest bond" << std::endl;
			return false;
		}
		if (!CompareDensityMatrices(ReferenceDensityMatrix(reg), mpo.getDensityMatrix(), nrQubits))
			return false;
	}

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		mpo.SetUseOptimalMeetingPosition(false);
		mpo.ApplyGate(cnotGate, 4, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
		if (!state || state->qubitsMap != expectedHeuristicMap)
		{
			std::cout << "MPO midpoint heuristic was not used after disabling optimal meeting position" << std::endl;
			return false;
		}
	}

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		mpo.SetUseOptimalMeetingPosition(false);
		mpo.ApplyGate(cnotGate, 2, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
		if (!state || state->qubitsMap != expectedHeuristicSwapUpMap)
		{
			std::cout << "MPO heuristic did not move the nearer-end qubit toward the other" << std::endl;
			return false;
		}
	}

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		prepareUnevenBonds(mpo, nullptr);
		mpo.SetMeetingPositionCallback([](const std::vector<Eigen::Index>&) { return Eigen::Index{ 2 }; });
		mpo.ApplyGate(cnotGate, 4, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
		const std::vector<Eigen::Index> expectedCallbackMap{ 2, 0, 1, 4, 3 };
		if (!state || state->qubitsMap != expectedCallbackMap)
		{
			std::cout << "MPO meeting position callback did not take priority over the local optimizer" << std::endl;
			return false;
		}
	}

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		prepareUnevenBonds(mpo, nullptr);
		mpo.SetMeetingPositionCallback([](const std::vector<Eigen::Index>&) { return Eigen::Index{ -1 }; });
		mpo.ApplyGate(cnotGate, 4, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
		if (!state || state->qubitsMap != expectedOptimalMap)
		{
			std::cout << "MPO invalid meeting position did not fall back to the local optimizer" << std::endl;
			return false;
		}
	}

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		mpo.SetUseOptimalMeetingPosition(false);
		mpo.SetMeetingPositionCallback([](const std::vector<Eigen::Index>&) { return Eigen::Index{ -1 }; });
		mpo.ApplyGate(cnotGate, 4, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
		if (!state || state->qubitsMap != expectedHeuristicMap)
		{
			std::cout << "MPO invalid meeting position did not fall back to the routing heuristic" << std::endl;
			return false;
		}
	}

	{
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		mpo.SetUseOptimalMeetingPosition(false);
		auto clone = mpo.Clone();
		clone->ApplyGate(cnotGate, 4, 0);

		const auto state = std::static_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(clone->getState());
		if (!state || state->qubitsMap != expectedHeuristicMap)
		{
			std::cout << "MPO clone did not copy the optimal meeting position flag" << std::endl;
			return false;
		}
	}

	std::cout << "\nSuccess" << std::endl;
	return true;
}

static bool InitialQubitsMapTestMPO()
{
	std::cout << "\nMPO simulator initial qubits map test" << std::endl;

	constexpr int nrQubits = 4;
	QC::TensorNetworks::MPOSimulator mpo(nrQubits);
	QC::QubitRegister<> reg(nrQubits);
	const std::vector<long long int> initialMap{ 2, 0, 3, 1 };
	mpo.SetInitialQubitsMap(initialMap);

	const auto mappedState = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	const std::vector<Eigen::Index> expectedMap{ 2, 0, 3, 1 };
	const std::vector<Eigen::Index> expectedInverse{ 1, 3, 0, 2 };
	if (!mappedState || mappedState->qubitsMap != expectedMap || mappedState->qubitsMapInv != expectedInverse)
	{
		std::cout << "MPO initial qubits map or its inverse was set incorrectly" << std::endl;
		return false;
	}

	QC::Gates::HadamardGate<> h;
	QC::Gates::CNOTGate<> cnot;
	mpo.ApplyGate(h, 0);
	mpo.ApplyGate(cnot, 3, 0);
	reg.ApplyGate(h, 0);
	reg.ApplyGate(cnot, 3, 0);
	if (!CompareDensityMatrices(ReferenceDensityMatrix(reg), mpo.getDensityMatrix(), nrQubits))
		return false;

	const Eigen::MatrixXcd beforeInvalidMap = mpo.getDensityMatrix();
	const auto stateBeforeInvalidMap = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	if (!MPO_ExpectInvalidArgument([&] { mpo.SetInitialQubitsMap({ 0, 1, 2 }); }, "Wrong-size initial MPO qubit map") ||
		!MPO_ExpectInvalidArgument([&] { mpo.SetInitialQubitsMap({ 0, 1, 1, 3 }); }, "Duplicate entry in initial MPO qubit map") ||
		!MPO_ExpectInvalidArgument([&] { mpo.SetInitialQubitsMap({ 0, 1, 2, -1 }); }, "Negative entry in initial MPO qubit map") ||
		!MPO_ExpectInvalidArgument([&] { mpo.SetInitialQubitsMap({ 0, 1, 2, 4 }); }, "Out-of-range entry in initial MPO qubit map"))
		return false;

	const auto stateAfterInvalidMap = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	if (!stateBeforeInvalidMap || !stateAfterInvalidMap ||
		stateAfterInvalidMap->qubitsMap != stateBeforeInvalidMap->qubitsMap ||
		stateAfterInvalidMap->qubitsMapInv != stateBeforeInvalidMap->qubitsMapInv ||
		(mpo.getDensityMatrix() - beforeInvalidMap).norm() > 1E-12)
	{
		std::cout << "A rejected initial MPO qubit map changed the simulator" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

static bool BondDimensionCallbackTestMPO()
{
	std::cout << "\nMPO simulator bond dimension callback test" << std::endl;

	QC::TensorNetworks::MPOSimulator mpo(4);
	QC::Gates::HadamardGate<> hGate;
	QC::Gates::CNOTGate<> cnotGate;

	int callbackCalls = 0;
	std::vector<Eigen::Index> lastBondDimensions;
	mpo.SetBondDimensionCallback(
		[&](const std::vector<Eigen::Index>& bondDimensions)
		{
			++callbackCalls;
			lastBondDimensions = bondDimensions;
		});

	// One-qubit operations cannot change a bond dimension and do not notify.
	mpo.ApplyGate(hGate, 0);
	if (callbackCalls != 0)
	{
		std::cout << "MPO bond dimension callback was called for a one-qubit operation" << std::endl;
		return false;
	}

	// Routing qubits 0 and 3 requires two swaps, followed by the two-qubit operation.
	mpo.ApplyGate(cnotGate, 3, 0);
	if (callbackCalls != 3 || lastBondDimensions != mpo.getBondDimensions())
	{
		std::cout << "MPO bond dimension callback did not report routing and gate updates" << std::endl;
		return false;
	}

	// Collapsing measurements notify as in MPSSimulator; sampling without collapse does not.
	mpo.MeasureQubit(0);
	if (callbackCalls != 4 || lastBondDimensions != mpo.getBondDimensions())
	{
		std::cout << "MPO bond dimension callback did not report a measurement update" << std::endl;
		return false;
	}
	mpo.MeasureNoCollapse();
	if (callbackCalls != 4)
	{
		std::cout << "MPO bond dimension callback was called for sampling without collapse" << std::endl;
		return false;
	}

	// MPO-specific two-qubit mutation paths notify too.
	mpo.ApplyOperatorAndNormalize(cnotGate, 3, 0);
	std::vector<Eigen::MatrixXcd> identityChannel{ Eigen::MatrixXcd::Identity(4, 4) };
	mpo.ApplyKrausOperators(identityChannel, 3, 0);
	if (callbackCalls != 6 || lastBondDimensions != mpo.getBondDimensions())
	{
		std::cout << "MPO bond dimension callback did not report normalized or Kraus updates" << std::endl;
		return false;
	}

	// Clones preserve the callback, and passing nullptr clears it.
	auto clone = mpo.Clone();
	clone->ApplyGate(cnotGate, 3, 0);
	if (callbackCalls != 7 || lastBondDimensions != clone->getBondDimensions())
	{
		std::cout << "MPO clone did not preserve the bond dimension callback" << std::endl;
		return false;
	}
	clone->SetBondDimensionCallback(nullptr);
	clone->ApplyGate(cnotGate, 3, 0);
	if (callbackCalls != 7)
	{
		std::cout << "MPO bond dimension callback was not cleared" << std::endl;
		return false;
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

static bool TruncationModeTestMPO()
{
	using TruncationMode = QC::TensorNetworks::MPOSimulatorInterface::TruncationMode;
	using MPOState = QC::TensorNetworks::MPOSimulatorBaseState;

	std::cout << "\nMPO simulator truncation mode test" << std::endl;

	// The default must be DiscardedWeight (Qiskit Aer's / ITensor's convention) - a deliberate
	// default-behavior change from what every earlier version of this simulator did. Same
	// rationale and API contract as MPSSimulatorInterface::TruncationMode; see MPSSimulatorTests.cpp's
	// TruncationModeTestMPS for the full explanation of the test strategy below.
	{
		QC::TensorNetworks::MPOSimulatorImpl defaultMpo(2);
		if (defaultMpo.getTruncationMode() != TruncationMode::DiscardedWeight)
		{
			std::cout << "Default truncation mode is not DiscardedWeight" << std::endl;
			return false;
		}
	}

	// getter/setter round trip for both modes
	{
		QC::TensorNetworks::MPOSimulatorImpl mpo(2);
		if (!mpo.setTruncationMode(TruncationMode::RelativeToMax) || mpo.getTruncationMode() != TruncationMode::RelativeToMax)
		{
			std::cout << "Truncation mode round trip failed for RelativeToMax" << std::endl;
			return false;
		}

		if (!mpo.setTruncationMode(TruncationMode::DiscardedWeight) || mpo.getTruncationMode() != TruncationMode::DiscardedWeight)
		{
			std::cout << "Truncation mode round trip failed for DiscardedWeight" << std::endl;
			return false;
		}
	}

	// Same construction as the MPS version: two Bell pairs (qubits 0-1 and 2-3) joined by a
	// generic two qubit unitary on qubits 1 and 2. The Bell-pair bonds stay far above any
	// threshold used below in every run, so the joining gate's theta matrix - and therefore its
	// raw SVD spectrum - is identical between the exact reference and every truncation-mode run.
	const QC::Gates::HadamardGate<> hGate;
	const QC::Gates::CNOTGate<> cnotGate;

	Eigen::MatrixXcd seed(4, 4);
	seed << std::complex<double>(0.3, 0.7), std::complex<double>(-0.5, 0.2), std::complex<double>(0.4, -0.6), std::complex<double>(0.1, 0.9),
		std::complex<double>(-0.8, 0.1), std::complex<double>(0.6, 0.5), std::complex<double>(0.2, 0.3), std::complex<double>(-0.4, 0.7),
		std::complex<double>(0.5, -0.3), std::complex<double>(0.9, -0.1), std::complex<double>(-0.6, 0.4), std::complex<double>(0.2, 0.2),
		std::complex<double>(-0.2, 0.6), std::complex<double>(0.3, -0.8), std::complex<double>(0.7, 0.1), std::complex<double>(-0.5, -0.3);
	const Eigen::MatrixXcd joinGate = Eigen::HouseholderQR<Eigen::MatrixXcd>(seed).householderQ();

	auto buildCircuit = [&](QC::TensorNetworks::MPOSimulatorImpl& mpo)
	{
		mpo.ApplyGate(hGate, 0);
		mpo.ApplyGate(cnotGate, 1, 0);
		mpo.ApplyGate(hGate, 2);
		mpo.ApplyGate(cnotGate, 3, 2);
		mpo.ApplyGate(QC::Gates::AppliedGate<>(joinGate, 2, 1));
	};

	constexpr Eigen::Index bondIndex = 1; // the bond between qubit 1 and qubit 2, joined last above

	QC::TensorNetworks::MPOSimulatorImpl mpoRef(4);
	buildCircuit(mpoRef);

	const auto refState = std::static_pointer_cast<MPOState>(mpoRef.getState());
	const Eigen::VectorXd spectrum = refState->lambdas[bondIndex]; // descending-sorted, per Eigen's SVD

	if (spectrum.size() < 3)
	{
		std::cout << "Test setup did not produce a rich enough spectrum (" << spectrum.size() << " singular values)" << std::endl;
		return false;
	}

	const double threshold = 0.3;

	// Independently compute, from the exact spectrum, how many singular values each mode should
	// keep - deliberately not calling MPOSimulatorImpl::ComputeCompressedRank, so this checks the
	// SPECIFICATION (as documented on MPOSimulatorInterface::TruncationMode), not just that the
	// implementation agrees with itself.
	Eigen::Index expectedRelative;
	{
		const double cutoffValue = threshold * spectrum[0];
		Eigen::Index i = spectrum.size() - 1;
		while (i >= 0 && spectrum[i] < cutoffValue) --i;
		expectedRelative = std::max<Eigen::Index>(i + 1, 1);
	}

	Eigen::Index expectedDiscardedWeight;
	{
		const double total = spectrum.squaredNorm();
		double discarded = 0.;
		Eigen::Index keep = spectrum.size();
		for (Eigen::Index i = spectrum.size() - 1; i > 0; --i)
		{
			const double sq = spectrum[i] * spectrum[i];
			if ((discarded + sq) / total >= threshold) break;

			discarded += sq;
			keep = i;
		}
		expectedDiscardedWeight = keep;
	}

	if (expectedRelative == expectedDiscardedWeight)
	{
		std::cout << "Test setup does not discriminate between the two truncation modes at threshold "
			<< threshold << " (both would keep " << expectedRelative << ")" << std::endl;
		return false;
	}

	// RelativeToMax explicitly requested must match the independently-computed expectation - this
	// is the regression test for "RelativeToMax still reproduces the original behavior".
	{
		QC::TensorNetworks::MPOSimulatorImpl mpoRel(4);
		mpoRel.setTruncationMode(TruncationMode::RelativeToMax);
		mpoRel.setLimitEntanglement(threshold);
		buildCircuit(mpoRel);

		const auto bondDims = mpoRel.getBondDimensions();
		if (bondDims[bondIndex] != expectedRelative)
		{
			std::cout << "RelativeToMax kept " << bondDims[bondIndex] << " singular values, expected " << expectedRelative << std::endl;
			return false;
		}
	}

	// DiscardedWeight explicitly requested must match the independently-computed expectation.
	{
		QC::TensorNetworks::MPOSimulatorImpl mpoWeight(4);
		mpoWeight.setTruncationMode(TruncationMode::DiscardedWeight);
		mpoWeight.setLimitEntanglement(threshold);
		buildCircuit(mpoWeight);

		const auto bondDims = mpoWeight.getBondDimensions();
		if (bondDims[bondIndex] != expectedDiscardedWeight)
		{
			std::cout << "DiscardedWeight kept " << bondDims[bondIndex] << " singular values, expected " << expectedDiscardedWeight << std::endl;
			return false;
		}
	}

	// No explicit mode set: must match DiscardedWeight (the default), not RelativeToMax - this is
	// the regression test guarding the default-flip decision.
	{
		QC::TensorNetworks::MPOSimulatorImpl mpoDefault(4);
		mpoDefault.setLimitEntanglement(threshold);
		buildCircuit(mpoDefault);

		const auto bondDims = mpoDefault.getBondDimensions();
		if (bondDims[bondIndex] != expectedDiscardedWeight)
		{
			std::cout << "Default truncation mode did not match DiscardedWeight (kept " << bondDims[bondIndex]
				<< ", expected " << expectedDiscardedWeight << ")" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;

	return true;
}

// Regression test for a real bug found during review: MPOSimulator::Clone() copies
// limitSize/limitEntanglement/chi/singularValueThreshold but, before this test existed,
// silently forgot to copy truncationMode - a clone would always come back at the default
// (DiscardedWeight) even if the original had explicitly opted into RelativeToMax.
static bool CloneTestMPO()
{
	using TruncationMode = QC::TensorNetworks::MPOSimulatorInterface::TruncationMode;

	std::cout << "\nMPO simulator clone preserves truncation mode" << std::endl;

	QC::TensorNetworks::MPOSimulator mpo(2);
	if (!mpo.setTruncationMode(TruncationMode::RelativeToMax))
	{
		std::cout << "Failed to set RelativeToMax on the original simulator" << std::endl;
		return false;
	}

	const auto cloned = mpo.Clone();
	if (cloned->getTruncationMode() != TruncationMode::RelativeToMax)
	{
		std::cout << "Clone did not preserve the RelativeToMax truncation mode (got DiscardedWeight instead)" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;

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

// SaveState/RestoreState and Clone must preserve the full state: tensors, bonds and, for the
// decorator, the logical->physical qubit map produced by swaps of non adjacent gates.
static bool StateSaveRestoreTestMPO()
{
	std::cout << "\nMPO simulator state save, restore and clone test" << std::endl;

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
			mpo.SaveState();

			// Move to a known different state, then clone. The clone must copy this current state while
			// also receiving its own independent copy of the earlier saved snapshot.
			const size_t allOnes = (1ULL << nrQubits) - 1;
			mpo.setToBasisState(allOnes);
			const Eigen::MatrixXcd rhoB = mpo.getDensityMatrix();
			auto clone = mpo.Clone();
			if (!CompareDensityMatrices(rhoB, clone->getDensityMatrix(), nrQubits))
				return false;

			// A regular restore is reusable.
			mpo.RestoreState();
			if (!CompareDensityMatrices(rhoA, mpo.getDensityMatrix(), nrQubits))
				return false;
			mpo.setToBasisState(allOnes);
			mpo.RestoreState();
			if (!CompareDensityMatrices(rhoA, mpo.getDensityMatrix(), nrQubits))
				return false;

			// A destructive restore consumes the original snapshot.
			mpo.setToBasisState(allOnes);
			mpo.RestoreStateDestructive();
			if (!CompareDensityMatrices(rhoA, mpo.getDensityMatrix(), nrQubits))
				return false;
			mpo.setToBasisState(0);
			const Eigen::MatrixXcd rhoWithoutSnapshot = mpo.getDensityMatrix();
			mpo.RestoreStateDestructive();
			if (!CompareDensityMatrices(rhoWithoutSnapshot, mpo.getDensityMatrix(), nrQubits))
				return false;

			// Consuming the original snapshot must not affect the clone's saved snapshot.
			clone->RestoreState();
			if (!CompareDensityMatrices(rhoA, clone->getDensityMatrix(), nrQubits))
				return false;
			clone->setToBasisState(allOnes);
			if (!CompareDensityMatrices(rhoWithoutSnapshot, mpo.getDensityMatrix(), nrQubits))
				return false;
			clone->RestoreStateDestructive();
			if (!CompareDensityMatrices(rhoA, clone->getDensityMatrix(), nrQubits))
				return false;

			std::cout << ".";
		}
	}

	// Clone also preserves simulation limits and the routing callback.
	{
		QC::TensorNetworks::MPOSimulator mpo(4);
		mpo.setLimitBondDimension(1);
		int callbackCalls = 0;
		mpo.SetMeetingPositionCallback([&callbackCalls](const std::vector<Eigen::Index>&) {
			++callbackCalls;
			return Eigen::Index(1);
		});

		auto clone = mpo.Clone();
		QC::Gates::HadamardGate<> h;
		QC::Gates::CNOTGate<> cnot;
		clone->ApplyGate(h, 0);
		clone->ApplyGate(cnot, 0, 3);

		if (callbackCalls != 1)
		{
			std::cout << "MPO clone did not copy the routing callback" << std::endl;
			return false;
		}

		for (const Eigen::Index dimension : clone->getBondDimensions())
			if (dimension > 1)
			{
				std::cout << "MPO clone did not copy the bond dimension limit" << std::endl;
				return false;
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

// samples the full computational basis state distribution using the MPO MeasureNoCollapse (which does
// not collapse the state, so the circuit is executed only once) and compares the empirical histogram
// against the DensityMatrix exact populations for the same prepared state
static bool SamplingNoCollapseTestMPO()
{
	std::cout << "\nMPO simulator - MeasureNoCollapse sampling against the DensityMatrix populations" << std::endl;

	std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>> gates;
	FillOneQubitGatesMPO(gates);
	FillTwoQubitGatesMPO(gates);

	const int nrSamples = 100000;

	for (int nrQubits = 2; nrQubits < 5; ++nrQubits)
	{
		// fixed random circuit prepared on both simulators, plus some noise so the state is mixed
		const auto circuit = BuildRandomCircuitMPO(gates, nrQubits, 25);

		QC::DensityMatrix<> dm(nrQubits);
		QC::TensorNetworks::MPOSimulator mpo(nrQubits);
		for (const auto& gate : circuit)
		{
			dm.ApplyGate(gate);
			mpo.ApplyGate(gate);
		}

		std::uniform_int_distribution qubitDistr(0, nrQubits - 1);
		const int noisyQubit = qubitDistr(gen);
		dm.ApplyDepolarizingNoise(noisyQubit, 0.3);
		mpo.ApplyDepolarizingNoise(noisyQubit, 0.3);

		// exact populations from the DensityMatrix reference
		const size_t nrBasisStates = dm.getNrBasisStates();
		std::vector<double> expected(nrBasisStates);
		for (size_t s = 0; s < nrBasisStates; ++s)
			expected[s] = dm.getBasisStateProbability(s);

		// sample the MPO state without collapsing it (the circuit is not re-run)
		std::vector<int> counts(nrBasisStates, 0);
		for (int m = 0; m < nrSamples; ++m)
		{
			const auto measured = mpo.MeasureNoCollapse();

			size_t state = 0;
			for (const auto& [qubit, val] : measured)
				if (val) state |= (1ULL << qubit);

			++counts[state];
		}

		for (size_t s = 0; s < nrBasisStates; ++s)
		{
			const double sampled = static_cast<double>(counts[s]) / nrSamples;
			if (std::abs(sampled - expected[s]) > 0.02)
			{
				std::cout << "MPO sampling mismatch for state " << s << " with " << nrQubits << " qubits: "
					<< expected[s] << " (DensityMatrix) vs " << sampled << " (MPO sampled)" << std::endl;
				return false;
			}
		}

		// subset sampling: sample only a subset of the qubits and compare against the exact marginal
		{
			std::set<Eigen::Index> subset;
			for (int q = 0; q < nrQubits; q += 2) // qubits 0, 2, 4, ...
				subset.insert(static_cast<Eigen::Index>(q));

			// exact marginal populations over the subset (indexed by the packed subset outcome)
			const size_t subsetStates = 1ULL << subset.size();
			std::vector<double> expectedSubset(subsetStates, 0.);
			for (size_t s = 0; s < nrBasisStates; ++s)
			{
				size_t idx = 0;
				int bit = 0;
				for (const Eigen::Index q : subset)
				{
					if (s & (1ULL << q)) idx |= (1ULL << bit);
					++bit;
				}
				expectedSubset[idx] += expected[s];
			}

			std::vector<int> subsetCounts(subsetStates, 0);
			for (int m = 0; m < nrSamples; ++m)
			{
				const auto measured = mpo.MeasureNoCollapse(subset);
				if (measured.size() != subset.size())
				{
					std::cout << "MPO subset MeasureNoCollapse returned the wrong number of qubits for " << nrQubits << " qubits" << std::endl;
					return false;
				}

				size_t idx = 0;
				int bit = 0;
				for (const Eigen::Index q : subset)
				{
					if (measured.at(q)) idx |= (1ULL << bit);
					++bit;
				}
				++subsetCounts[idx];
			}

			for (size_t s = 0; s < subsetStates; ++s)
			{
				const double sampled = static_cast<double>(subsetCounts[s]) / nrSamples;
				if (std::abs(sampled - expectedSubset[s]) > 0.02)
				{
					std::cout << "MPO subset sampling mismatch for outcome " << s << " with " << nrQubits << " qubits: "
						<< expectedSubset[s] << " (DensityMatrix) vs " << sampled << " (MPO sampled)" << std::endl;
					return false;
				}
			}
		}

		std::cout << ".";
	}

	std::cout << "\nSuccess" << std::endl;

	return true;
}

// Locks the public API validation contract. Invalid calls must throw before changing either the
// represented logical state or the decorator's logical-to-physical qubit permutation.
static bool ValidationAndStateCompatibilityTestMPO()
{
	std::cout << "\nMPO simulator - validation and state compatibility" << std::endl;

	QC::TensorNetworks::MPOSimulator mpo(4);
	QC::Gates::HadamardGate<> h;
	QC::Gates::CNOTGate<> cnot;
	mpo.ApplyGate(h, 0);
	mpo.ApplyGate(cnot, 3, 0); // also make the decorator map non-identity

	const Eigen::MatrixXcd before = mpo.getDensityMatrix();
	const auto beforeState = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	if (!beforeState)
	{
		std::cout << "MPO getState returned an unexpected state type" << std::endl;
		return false;
	}

	Eigen::MatrixXcd empty(0, 0);
	Eigen::MatrixXcd rectangular = Eigen::MatrixXcd::Zero(2, 1);
	Eigen::MatrixXcd three = Eigen::MatrixXcd::Identity(3, 3);
	Eigen::MatrixXcd eight = Eigen::MatrixXcd::Identity(8, 8);
	Eigen::MatrixXcd nonFinite = Eigen::MatrixXcd::Identity(2, 2);
	nonFinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
	const Eigen::MatrixXcd i2 = Eigen::MatrixXcd::Identity(2, 2);
	const Eigen::MatrixXcd i4 = Eigen::MatrixXcd::Identity(4, 4);
	const MPOWrongDeclaredArityGate wrongDeclaredArity(i4);

	if (!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(QC::Gates::AppliedGate<>(empty, 0)); }, "Empty operator matrix") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(QC::Gates::AppliedGate<>(rectangular, 0)); }, "Rectangular operator matrix") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(QC::Gates::AppliedGate<>(three, 0)); }, "Non-power-of-two operator matrix") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(QC::Gates::AppliedGate<>(eight, 3, 0)); }, "Unsupported three-qubit operator") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(QC::Gates::AppliedGate<>(nonFinite, 0)); }, "Non-finite operator matrix") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(wrongDeclaredArity, 3, 0); }, "Operator with mismatched declared arity") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyOperator(QC::Gates::TwoQubitsGate<>(i4), 0, 0); }, "Duplicate operator qubits") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyKrausOperators({ i2, i4 }, 0); }, "Mismatched Kraus dimensions") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyKrausOperators({ i2, nonFinite }, 0); }, "Non-finite Kraus operator") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyKrausOperators({ i4 }, 0, 0); }, "Duplicate Kraus qubits"))
		return false;

	using MPOIndex = QC::TensorNetworks::MPOSimulatorInterface::IndexType;
	if (!MPO_ExpectInvalidArgument([&] { mpo.MeasureQubits(std::set<MPOIndex>{ -1 }); }, "Negative measured qubit") ||
		!MPO_ExpectInvalidArgument([&] { mpo.MeasureQubits(std::set<MPOIndex>{ 0, 4 }); }, "Out-of-range measured qubit") ||
		!MPO_ExpectInvalidArgument([&] { mpo.MeasureNoCollapse(std::set<MPOIndex>{ 4 }); }, "Out-of-range sampled qubit") ||
		!MPO_ExpectInvalidArgument([&] { mpo.MoveAtBeginningOfChain(std::set<MPOIndex>{ -1, 0 }); }, "Invalid moved qubit") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setLimitBondDimension(0); }, "Zero bond dimension limit") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setLimitBondDimension(-2); }, "Negative bond dimension limit") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setLimitEntanglement(-1E-6); }, "Negative singular-value threshold") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setLimitEntanglement(std::numeric_limits<double>::infinity()); }, "Infinite singular-value threshold"))
		return false;

	if (!MPO_ExpectInvalidArgument([&] { mpo.ApplyBitFlipNoise(0, -0.1); }, "Negative bit-flip probability") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyDepolarizingNoise(0, 1.1); }, "Depolarizing probability above one") ||
		!MPO_ExpectInvalidArgument([&] { mpo.ApplyAmplitudeDamping(0, std::numeric_limits<double>::quiet_NaN()); }, "NaN damping probability"))
		return false;

	const auto afterInvalidOperations = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	if (!afterInvalidOperations || beforeState->qubitsMap != afterInvalidOperations->qubitsMap ||
		beforeState->qubitsMapInv != afterInvalidOperations->qubitsMapInv ||
		(mpo.getDensityMatrix() - before).norm() > 1E-12)
	{
		std::cout << "An invalid MPO operation changed the state or qubit map" << std::endl;
		return false;
	}

	// States are implementation-specific: accepting a decorator state in the adjacent implementation
	// would silently discard its qubit map, while downcasting an implementation state is unsafe.
	QC::TensorNetworks::MPOSimulatorImpl adjacent(4);
	adjacent.ApplyGate(h, 0);
	const Eigen::MatrixXcd adjacentBeforeInvalidMeasurement = adjacent.getDensityMatrix();
	if (!MPO_ExpectInvalidArgument([&] { adjacent.MeasureQubits(std::set<MPOIndex>{ 0, 4 }); }, "Implementation subset with an invalid measured qubit") ||
		!MPO_ExpectInvalidArgument([&] { adjacent.MeasureNoCollapse(std::set<MPOIndex>{ 4 }); }, "Implementation subset with an invalid sampled qubit") ||
		!MPO_ExpectInvalidArgument([&] { adjacent.MoveAtBeginningOfChain(std::set<MPOIndex>{ -1 }); }, "Implementation move with an invalid qubit"))
		return false;
	if ((adjacent.getDensityMatrix() - adjacentBeforeInvalidMeasurement).norm() > 1E-12)
	{
		std::cout << "An invalid implementation measurement partially collapsed the MPO" << std::endl;
		return false;
	}

	auto adjacentState = adjacent.getState();
	if (!MPO_ExpectInvalidArgument([&] { mpo.setState(adjacentState); }, "Implementation state restored into decorator"))
		return false;

	auto decoratedState = mpo.getState();
	if (!MPO_ExpectInvalidArgument([&] { adjacent.setState(decoratedState); }, "Decorator state restored into implementation"))
		return false;

	QC::TensorNetworks::MPOSimulator otherSize(3);
	auto otherSizeState = otherSize.getState();
	if (!MPO_ExpectInvalidArgument([&] { mpo.setState(otherSizeState); }, "Wrong-size MPO state"))
		return false;

	auto malformedState = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mpo.getState());
	malformedState->qubitsMap[0] = malformedState->qubitsMap[1];
	if (!MPO_ExpectInvalidArgument([&] { mpo.setState(malformedState); }, "Malformed MPO qubit map"))
		return false;

	if ((mpo.getDensityMatrix() - before).norm() > 1E-12)
	{
		std::cout << "Rejected state restoration changed the MPO" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

// Positive mixture weights are scale-free and finite. Very small or very large common scales must
// therefore produce the same normalized state, while unusable mixtures must be rejected atomically.
static bool MixtureValidationAndScalingTestMPO()
{
	std::cout << "\nMPO simulator - mixture validation and stable normalization" << std::endl;

	for (const double scale : { 1E-300, std::numeric_limits<double>::max() / 4. })
	{
		QC::TensorNetworks::MPOSimulator mpo(2);
		mpo.setToMixtureOfBasisStates({ std::pair<size_t, double>{ 0, scale }, std::pair<size_t, double>{ 1, 2. * scale } });
		if (!approxEqual(mpo.Trace(), std::complex<double>(1., 0.), 1E-12) ||
			!approxEqual(mpo.getBasisStateProbability(0), 1. / 3., 1E-12) ||
			!approxEqual(mpo.getBasisStateProbability(1), 2. / 3., 1E-12))
		{
			std::cout << "MPO mixture normalization is not scale invariant" << std::endl;
			return false;
		}
	}

	QC::TensorNetworks::MPOSimulator mpo(2);
	mpo.setToBasisState(3);
	const Eigen::MatrixXcd before = mpo.getDensityMatrix();
	if (!MPO_ExpectInvalidArgument([&] { mpo.setToMixtureOfBasisStates(std::vector<std::pair<size_t, double>>{}); }, "Empty MPO mixture") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setToMixtureOfBasisStates({ std::pair<size_t, double>{ 0, std::numeric_limits<double>::infinity() } }); }, "Infinite MPO mixture weight") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setToMixtureOfBasisStates({ std::pair<size_t, double>{ 7, 1. }, std::pair<size_t, double>{ 0, 0. } }); }, "MPO mixture without a valid positive term") ||
		!MPO_ExpectInvalidArgument([&] { mpo.setToMixtureOfBasisStates(std::vector<std::pair<std::vector<bool>, double>>{ { std::vector<bool>{ true, false, true }, 1. } }); }, "Oversized vector MPO mixture state"))
		return false;

	if ((mpo.getDensityMatrix() - before).norm() > 1E-12)
	{
		std::cout << "An invalid MPO mixture changed the state" << std::endl;
		return false;
	}

	// Invalid and non-positive entries are ignored when at least one usable term remains.
	mpo.setToMixtureOfBasisStates({ std::pair<size_t, double>{ 1, 1. }, std::pair<size_t, double>{ 7, 9. }, std::pair<size_t, double>{ 2, -4. } });
	if (!approxEqual(mpo.getBasisStateProbability(1), 1., 1E-12) ||
		!approxEqual(mpo.Trace(), std::complex<double>(1., 0.), 1E-12))
	{
		std::cout << "MPO did not ignore unusable mixture entries consistently" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

// The size_t state overload remains useful for wide tensor networks. In particular it must not
// evaluate a width-sized shift when N equals the number of bits in size_t. The vector overload
// must be able to initialize qubits whose indices cannot be represented in a size_t bit mask.
static bool WideBasisInitializationTestMPO()
{
	std::cout << "\nMPO simulator - wide basis-state initialization" << std::endl;

	constexpr size_t digits = std::numeric_limits<size_t>::digits;
	constexpr size_t allBits = std::numeric_limits<size_t>::max();
	for (const size_t nrQubits : { digits, digits + 1 })
	{
		QC::TensorNetworks::MPOSimulatorImpl mpo(nrQubits);
		mpo.setToBasisState(allBits);

		if (!approxEqual(mpo.Trace(), std::complex<double>(1., 0.), 1E-12) ||
			!approxEqual(mpo.GetProbability(0, false), 1., 1E-12) ||
			!approxEqual(mpo.GetProbability(static_cast<Eigen::Index>(digits - 1), false), 1., 1E-12))
		{
			std::cout << "Wide MPO basis initialization lost a representable size_t bit" << std::endl;
			return false;
		}

		if (nrQubits > digits && !approxEqual(mpo.GetProbability(static_cast<Eigen::Index>(digits), true), 1., 1E-12))
		{
			std::cout << "Wide MPO basis initialization did not zero-extend the state" << std::endl;
			return false;
		}
	}

	const size_t vectorQubits = digits + 3;
	std::vector<bool> vectorStateBits(vectorQubits, false);
	vectorStateBits[0] = true;
	vectorStateBits[digits - 1] = true;
	vectorStateBits[digits] = true;
	vectorStateBits[vectorQubits - 1] = true;
	const std::vector<bool> vectorState = std::move(vectorStateBits);

	QC::TensorNetworks::MPOSimulator vectorMpo(vectorQubits);
	QC::TensorNetworks::MPOSimulatorInterface& vectorMpoApi = vectorMpo;
	vectorMpoApi.setToBasisState(vectorState);

	if (!approxEqual(vectorMpo.Trace(), std::complex<double>(1., 0.), 1E-12))
	{
		std::cout << "Vector-initialized wide MPO does not have unit trace" << std::endl;
		return false;
	}

	for (size_t q = 0; q < vectorQubits; ++q)
	{
		const double expectedOneProbability = vectorState[q] ? 1. : 0.;
		if (!approxEqual(vectorMpo.GetProbability(static_cast<Eigen::Index>(q), false), expectedOneProbability, 1E-12))
		{
			std::cout << "Vector MPO basis initialization set qubit " << q << " incorrectly" << std::endl;
			return false;
		}
	}

	const std::vector<bool> shortState{ true, false, true };
	vectorMpoApi.setToBasisState(shortState);
	if (!approxEqual(vectorMpo.GetProbability(0, false), 1., 1E-12) ||
		!approxEqual(vectorMpo.GetProbability(1, false), 0., 1E-12) ||
		!approxEqual(vectorMpo.GetProbability(2, false), 1., 1E-12) ||
		!approxEqual(vectorMpo.GetProbability(static_cast<Eigen::Index>(vectorQubits - 1), false), 0., 1E-12))
	{
		std::cout << "Short MPO basis vector was not zero-extended" << std::endl;
		return false;
	}

	const std::vector<bool> emptyState;
	vectorMpoApi.setToBasisState(emptyState);
	if (!approxEqual(vectorMpo.GetProbability(0, false), 0., 1E-12) ||
		!approxEqual(vectorMpo.GetProbability(static_cast<Eigen::Index>(vectorQubits - 1), false), 0., 1E-12))
	{
		std::cout << "Empty MPO basis vector did not initialize the all-zero state" << std::endl;
		return false;
	}

	QC::TensorNetworks::MPOSimulator mappedMpo(3);
	mappedMpo.setToBasisState(std::vector<bool>{ true, false, true });
	mappedMpo.MoveAtBeginningOfChain({ 2 });
	const auto stateBeforeInvalidCall = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mappedMpo.getState());
	const std::vector<Eigen::Index> identityMap{ 0, 1, 2 };
	if (!stateBeforeInvalidCall || stateBeforeInvalidCall->qubitsMap == identityMap)
	{
		std::cout << "MPO test setup did not create a non-identity logical qubit mapping" << std::endl;
		return false;
	}

	if (!MPO_ExpectInvalidArgument([&] { mappedMpo.setToBasisState(std::vector<bool>(4, false)); }, "Oversized MPO basis vector"))
		return false;

	const auto stateAfterInvalidCall = std::dynamic_pointer_cast<QC::TensorNetworks::MPOSimulatorState>(mappedMpo.getState());
	if (!stateAfterInvalidCall || stateAfterInvalidCall->qubitsMap != stateBeforeInvalidCall->qubitsMap ||
		stateAfterInvalidCall->qubitsMapInv != stateBeforeInvalidCall->qubitsMapInv ||
		!approxEqual(mappedMpo.GetProbability(0, false), 1., 1E-12) ||
		!approxEqual(mappedMpo.GetProbability(1, false), 0., 1E-12) ||
		!approxEqual(mappedMpo.GetProbability(2, false), 1., 1E-12))
	{
		std::cout << "Rejected MPO basis initialization changed the state or logical qubit mapping" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

// Kraus summation must not direct-sum every bond. With compression enabled the configured cap is a
// hard invariant, and a local channel on a product state should not create remote virtual bonds.
static bool KrausBondLimitAndLocalityTestMPO()
{
	std::cout << "\nMPO simulator - Kraus bond limits and locality" << std::endl;

	{
		QC::TensorNetworks::MPOSimulator product(5);
		for (int i = 0; i < 4; ++i)
			product.ApplyDepolarizingNoise(2, 0.2 + 0.1 * i);

		for (const Eigen::Index bond : product.getBondDimensions())
			if (bond != 1)
			{
				std::cout << "A local channel grew a remote/product-state MPO bond to " << bond << std::endl;
				return false;
			}
	}

	QC::TensorNetworks::MPOSimulator mpo(4);
	mpo.setLimitBondDimension(1);
	QC::Gates::HadamardGate<> h;
	QC::Gates::CNOTGate<> cnot;
	mpo.ApplyGate(h, 0);
	mpo.ApplyGate(cnot, 3, 0);

	const Eigen::MatrixXcd cnotMatrix = cnot.getRawOperatorMatrix();
	const Eigen::MatrixXcd k0 = std::sqrt(0.6) * Eigen::MatrixXcd::Identity(4, 4);
	const Eigen::MatrixXcd k1 = std::sqrt(0.4) * cnotMatrix;
	for (int i = 0; i < 3; ++i)
	{
		mpo.ApplyAmplitudeDamping(i, 0.25);
		mpo.ApplyDepolarizingNoise(3 - i, 0.3);
		mpo.ApplyKrausOperators({ k0, k1 }, 3, 0);

		for (const Eigen::Index bond : mpo.getBondDimensions())
			if (bond > 1)
			{
				std::cout << "Kraus application exceeded the configured MPO bond limit: " << bond << std::endl;
				return false;
			}
	}

	const Eigen::MatrixXcd rho = mpo.getDensityMatrix();
	if (!rho.allFinite() || !std::isfinite(mpo.Trace().real()) || !std::isfinite(mpo.Trace().imag()))
	{
		std::cout << "Kraus evolution with a bond cap produced a non-finite MPO" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

// Measurement probabilities are computed relative to the current trace, but a collapsed density
// matrix must itself be normalized. Impossible normalized operations must fail atomically.
static bool CollapseNormalizationAndAtomicFailureTestMPO()
{
	std::cout << "\nMPO simulator - collapse normalization and atomic failure" << std::endl;

	QC::TensorNetworks::MPOSimulator mpo(2);
	QC::Gates::HadamardGate<> h;
	QC::Gates::CNOTGate<> cnot;
	mpo.ApplyGate(h, 0);
	mpo.ApplyGate(cnot, 1, 0);

	const double scale = 0.37;
	mpo.ApplyOperator(QC::Gates::SingleQubitGate<>(std::sqrt(scale) * Eigen::MatrixXcd::Identity(2, 2)), 0);
	if (!approxEqual(mpo.Trace(), std::complex<double>(scale, 0.), 1E-12))
	{
		std::cout << "Non-unitary setup did not produce the expected MPO trace" << std::endl;
		return false;
	}

	const bool outcome = mpo.MeasureQubit(0);
	const Eigen::MatrixXcd collapsed = mpo.getDensityMatrix();
	const size_t expectedState = outcome ? 3 : 0;
	if (!approxEqual(mpo.Trace(), std::complex<double>(1., 0.), 1E-10) ||
		!approxEqual(collapsed(static_cast<Eigen::Index>(expectedState), static_cast<Eigen::Index>(expectedState)), std::complex<double>(1., 0.), 1E-10) ||
		(collapsed - collapsed.adjoint()).norm() > 1E-10)
	{
		std::cout << "MPO measurement did not produce the normalized projected state" << std::endl;
		return false;
	}

	QC::TensorNetworks::MPOSimulatorImpl impossible(1);
	const Eigen::MatrixXcd before = impossible.getDensityMatrix();
	Eigen::MatrixXcd projectOne = Eigen::MatrixXcd::Zero(2, 2);
	projectOne(1, 1) = 1.;
	if (!MPO_ExpectRuntimeError([&] { impossible.ApplyOperatorAndNormalize(QC::Gates::SingleQubitGate<>(projectOne), 0); }, "Impossible normalized MPO operation"))
		return false;
	if ((impossible.getDensityMatrix() - before).norm() > 1E-12)
	{
		std::cout << "An impossible normalized MPO operation corrupted the original state" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}

bool MPOSimulatorTests()
{
	std::cout << "\nMPO Simulator Tests" << std::endl;
	return ValidationAndStateCompatibilityTestMPO() &&
		InitialQubitsMapTestMPO() &&
		MixtureValidationAndScalingTestMPO() &&
		WideBasisInitializationTestMPO() &&
		KrausBondLimitAndLocalityTestMPO() &&
		CollapseNormalizationAndAtomicFailureTestMPO() &&
		OneAndTwoQubitGatesTestMPO() &&
		NonAdjacentGatesTestMPO() &&
		NumericalRankStabilityTestMPO() &&
		MeetingPositionCallbackTestMPO() &&
		OptimalMeetingPositionTestMPO() &&
		BondDimensionCallbackTestMPO() &&
		TraceAndProbabilitiesTestMPO() &&
		MixtureOfBasisStatesTestMPO() &&
		MixtureEvolutionTestMPO() &&
		KrausOperatorsTestMPO() &&
		TwoQubitKrausOperatorsTestMPO() &&
		ApplyOperatorAndNormalizeTestMPO() &&
		CompressionLosslessTestMPO() &&
		CompressionTruncationTestMPO() &&
		TruncationModeTestMPO() &&
		CloneTestMPO() &&
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
		MeasurementVsDensityMatrixAndThrowsTestMPO() &&
		SamplingNoCollapseTestMPO();
}
