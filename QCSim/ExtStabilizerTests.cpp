#include "Tests.h"
#include "QubitRegister.h"
#include "ExtendedStabilizer.h"

#include <limits>
#include <utility>


void ApplyTwoQubitsGate(QC::ExtendedStabilizer& simulator, int code, int qubit1, int qubit2)
{
	switch (code)
	{
	case 9:
		simulator.ApplyCX(qubit1, qubit2);
		break;
	case 10:
		simulator.ApplyCY(qubit1, qubit2);
		break;
	case 11:
		simulator.ApplyCZ(qubit1, qubit2);
		break;
	case 12:
		simulator.ApplySwap(qubit1, qubit2);
		break;
	case 13:
		simulator.ApplyISwap(qubit1, qubit2);
		break;
	case 14:
		simulator.ApplyISwapDag(qubit1, qubit2);
		break;
	}
}


void ApplyGate(QC::ExtendedStabilizer& simulator, int code, int qubit1, int qubit2, double angle = 0.0)
{
	switch (code)
	{
	case 0:
		simulator.ApplyH(qubit1);
		break;
	case 1:
		simulator.ApplyS(qubit1);
		break;
	case 2:
		simulator.ApplySdg(qubit1);
		break;
	case 3:
		simulator.ApplyX(qubit1);
		break;
	case 4:
		simulator.ApplyY(qubit1);
		break;
	case 5:
		simulator.ApplyZ(qubit1);
		break;
	case 6:
		simulator.ApplySx(qubit1);
		break;
	case 7:
		simulator.ApplySxDag(qubit1);
		break;
	case 8:
		simulator.ApplyK(qubit1);
		break;
	case 15:
		simulator.ApplyRx(qubit1, angle);
		break;
	case 16:
		simulator.ApplyRy(qubit1, angle);
		break;
	case 17:
		simulator.ApplyRz(qubit1, angle);
		break;
	default:
		ApplyTwoQubitsGate(simulator, code, qubit1, qubit2);
		break;
	}
}


static void ApplyExtStabilizerTestGate(QC::QubitRegister<>& qubitRegister, QC::ExtendedStabilizer& simulator,
	int code, size_t qubit1, size_t qubit2 = 0, double angle = 0.0)
{
	auto gate = GetGate(code, angle);
	qubitRegister.ApplyGate(*gate, qubit1, qubit2);
	ApplyGate(simulator, code, static_cast<int>(qubit1), static_cast<int>(qubit2), angle);
}


// Comparing every Pauli expectation completely characterizes a small pure state
// up to global phase, and catches relative-phase errors that probabilities miss.
static bool CheckExtStabilizerState(QC::QubitRegister<>& qubitRegister, const QC::ExtendedStabilizer& simulator,
	const std::string& context, double tolerance = 1E-8)
{
	static const QC::Gates::PauliXGate<> xGate;
	static const QC::Gates::PauliYGate<> yGate;
	static const QC::Gates::PauliZGate<> zGate;

	const size_t nrQubits = simulator.GetNrQubits();
	for (size_t q = 0; q < nrQubits; ++q)
	{
		const double expected = qubitRegister.GetQubitProbability(q);
		const double actual = simulator.GetQubitProbability(q);
		if (!approxEqual(expected, actual, tolerance))
		{
			std::cout << "\n" << context << ": probability mismatch for qubit " << q
				<< ", statevector " << expected << ", extended stabilizer " << actual << std::endl;
			return false;
		}
	}

	const size_t nrPauliStrings = 1ULL << (2 * nrQubits);
	for (size_t encoded = 0; encoded < nrPauliStrings; ++encoded)
	{
		std::string pauliString(nrQubits, 'I');
		std::vector<QC::Gates::AppliedGate<>> pauliGates;
		pauliGates.reserve(nrQubits);

		for (size_t q = 0; q < nrQubits; ++q)
		{
			const int pauli = static_cast<int>((encoded >> (2 * q)) & 3ULL);
			switch (pauli)
			{
			case 1:
				pauliString[q] = 'X';
				pauliGates.emplace_back(xGate.getRawOperatorMatrix(), q);
				break;
			case 2:
				pauliString[q] = 'Y';
				pauliGates.emplace_back(yGate.getRawOperatorMatrix(), q);
				break;
			case 3:
				pauliString[q] = 'Z';
				pauliGates.emplace_back(zGate.getRawOperatorMatrix(), q);
				break;
			default:
				break;
			}
		}

		const auto expectedComplex = qubitRegister.ExpectationValue(pauliGates);
		const double expected = expectedComplex.real();
		const double actual = simulator.ExpectationValue(pauliString);
		if (std::abs(expectedComplex.imag()) > tolerance || !approxEqual(expected, actual, tolerance))
		{
			std::cout << "\n" << context << ": expectation mismatch for Pauli string " << pauliString
				<< ", statevector " << expectedComplex << ", extended stabilizer " << actual << std::endl;
			return false;
		}
	}

	return true;
}


static std::complex<double> StatevectorPauliExpectation(QC::QubitRegister<>& qubitRegister,
	const std::string& pauliString)
{
	static const QC::Gates::PauliXGate<> xGate;
	static const QC::Gates::PauliYGate<> yGate;
	static const QC::Gates::PauliZGate<> zGate;

	std::vector<QC::Gates::AppliedGate<>> pauliGates;
	pauliGates.reserve(pauliString.size());
	for (size_t q = 0; q < pauliString.size(); ++q)
		switch (pauliString[q])
		{
		case 'X':
			pauliGates.emplace_back(xGate.getRawOperatorMatrix(), q);
			break;
		case 'Y':
			pauliGates.emplace_back(yGate.getRawOperatorMatrix(), q);
			break;
		case 'Z':
			pauliGates.emplace_back(zGate.getRawOperatorMatrix(), q);
			break;
		default:
			break;
		}

	return qubitRegister.ExpectationValue(pauliGates);
}


// Larger states cannot be checked by enumerating all 4^n Pauli strings. This
// deterministic sample includes every single-qubit Pauli, nearest-neighbour
// correlations, overlapping three-body strings, full-register strings, and a
// reproducible collection of Pauli strings of mixed weight.
static bool CheckExtStabilizerSampledState(QC::QubitRegister<>& qubitRegister,
	const QC::ExtendedStabilizer& simulator, const std::string& context,
	unsigned int sampleSeed, size_t nrRandomStrings, double tolerance = 1E-7)
{
	const size_t nrQubits = simulator.GetNrQubits();
	for (size_t q = 0; q < nrQubits; ++q)
	{
		const double expected = qubitRegister.GetQubitProbability(q);
		const double actual = simulator.GetQubitProbability(q);
		if (!approxEqual(expected, actual, tolerance))
		{
			std::cout << "\n" << context << ": probability mismatch for qubit " << q
				<< ", statevector " << expected << ", extended stabilizer " << actual << std::endl;
			return false;
		}
	}

	std::vector<std::string> pauliStrings;
	pauliStrings.reserve(7 * nrQubits + 3 + nrRandomStrings);
	static constexpr char axes[] = { 'X', 'Y', 'Z' };

	for (size_t q = 0; q < nrQubits; ++q)
		for (const char axis : axes)
		{
			std::string pauli(nrQubits, 'I');
			pauli[q] = axis;
			pauliStrings.push_back(std::move(pauli));
		}

	for (size_t q = 0; q < nrQubits; ++q)
		for (const char axis : axes)
		{
			std::string pauli(nrQubits, 'I');
			pauli[q] = axis;
			pauli[(q + 1) % nrQubits] = axis;
			pauliStrings.push_back(std::move(pauli));
		}

	for (size_t q = 0; q < nrQubits; ++q)
	{
		std::string pauli(nrQubits, 'I');
		pauli[q] = 'X';
		pauli[(q + 1) % nrQubits] = 'Y';
		pauli[(q + 2) % nrQubits] = 'Z';
		pauliStrings.push_back(std::move(pauli));
	}

	for (size_t offset = 0; offset < 3; ++offset)
	{
		std::string pauli(nrQubits, 'I');
		for (size_t q = 0; q < nrQubits; ++q)
			pauli[q] = axes[(q + offset) % 3];
		pauliStrings.push_back(std::move(pauli));
	}

	std::mt19937 pauliGenerator(sampleSeed);
	std::uniform_int_distribution<int> pauliDistribution(0, 3);
	for (size_t sample = 0; sample < nrRandomStrings; ++sample)
	{
		std::string pauli(nrQubits, 'I');
		bool nonIdentity = false;
		for (size_t q = 0; q < nrQubits; ++q)
		{
			const int value = pauliDistribution(pauliGenerator);
			if (value != 0)
			{
				pauli[q] = axes[value - 1];
				nonIdentity = true;
			}
		}
		if (!nonIdentity)
			pauli[sample % nrQubits] = axes[sample % 3];
		pauliStrings.push_back(std::move(pauli));
	}

	for (const auto& pauliString : pauliStrings)
	{
		const auto expected = StatevectorPauliExpectation(qubitRegister, pauliString);
		const double actual = simulator.ExpectationValue(pauliString);
		if (std::abs(expected.imag()) > tolerance || !approxEqual(expected.real(), actual, tolerance))
		{
			std::cout << "\n" << context << ": expectation mismatch for Pauli string " << pauliString
				<< ", statevector " << expected << ", extended stabilizer " << actual << std::endl;
			return false;
		}
	}

	return true;
}


static bool TestExtStabilizerLargerMixedCircuits()
{
	std::cout << "\nExtended Stabilizer larger mixed-circuit tests" << std::endl;

	struct CircuitConfiguration {
		size_t nrQubits;
		size_t nrGates;
		unsigned int seed;
	};
	const std::vector<CircuitConfiguration> configurations = {
		{ 5, 48, 0x51A7C001U },
		{ 6, 56, 0x61A7C002U },
		{ 7, 64, 0x71A7C003U },
		{ 8, 72, 0x81A7C004U }
	};

	for (const auto& configuration : configurations)
	{
		std::mt19937 circuitGenerator(configuration.seed);
		std::uniform_int_distribution<int> qubitDistribution(0,
			static_cast<int>(configuration.nrQubits) - 1);
		std::uniform_int_distribution<int> singleCliffordDistribution(0, 8);
		std::uniform_int_distribution<int> twoQubitCliffordDistribution(9, 14);
		std::uniform_int_distribution<int> secondQubitOffsetDistribution(1,
			static_cast<int>(configuration.nrQubits) - 1);
		std::uniform_real_distribution<double> angleDistribution(-2.4, 2.4);

		QC::QubitRegister<> qubitRegister(configuration.nrQubits);
		QC::ExtendedStabilizer simulator(configuration.nrQubits);
		size_t nrRotations = 0;

		for (size_t gateIndex = 0; gateIndex < configuration.nrGates; ++gateIndex)
		{
			const size_t qubit1 = static_cast<size_t>(qubitDistribution(circuitGenerator));
			const size_t qubit2 = (qubit1
				+ static_cast<size_t>(secondQubitOffsetDistribution(circuitGenerator))) % configuration.nrQubits;

			int code;
			double angle = 0.0;
			if (gateIndex % 4 == 3)
			{
				// Cycle through every supported non-Clifford axis and use a
				// non-special, reproducibly generated angle for each rotation.
				code = 15 + static_cast<int>(nrRotations % 3);
				angle = angleDistribution(circuitGenerator);
				++nrRotations;
			}
			else if (gateIndex % 3 == 1)
				code = twoQubitCliffordDistribution(circuitGenerator);
			else
				code = singleCliffordDistribution(circuitGenerator);

			ApplyExtStabilizerTestGate(qubitRegister, simulator, code, qubit1, qubit2, angle);

			if (gateIndex + 1 == configuration.nrGates / 2)
			{
				const std::string context = "Large mixed circuit midpoint on "
					+ std::to_string(configuration.nrQubits) + " qubits";
				if (!CheckExtStabilizerSampledState(qubitRegister, simulator, context,
					configuration.seed ^ 0x13579BDFU, 12))
					return false;
			}
		}

		const std::string context = "Large mixed circuit final state on "
			+ std::to_string(configuration.nrQubits) + " qubits after "
			+ std::to_string(configuration.nrGates) + " gates and "
			+ std::to_string(nrRotations) + " rotations";
		if (!CheckExtStabilizerSampledState(qubitRegister, simulator, context,
			configuration.seed ^ 0x2468ACE0U, 24))
			return false;

		std::cout << '.';
	}

	std::cout << "\nSuccess" << std::endl;
	return true;
}


static bool TestExtStabilizerClone()
{
	std::cout << "\nExtended Stabilizer clone tests" << std::endl;

	QC::QubitRegister<> qubitRegister(3);
	QC::ExtendedStabilizer simulator(3);

	// Prepare an entangled, coherent state A and retain it as the saved snapshot.
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 1, 0);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 0, 0, 0.37);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 2, 0, -0.41);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 2, 1);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 1, 0, 0.23);
	simulator.SaveState();
	auto savedStatevector = qubitRegister.Clone();

	// Evolve to state B, which must be the clone's current state.
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 2, 0, -0.18);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 11, 0, 2);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 0, 0, 0.31);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 1, 1);
	auto clone = simulator.Clone();

	if (!CheckExtStabilizerState(qubitRegister, simulator,
		"Clone test original state B")
		|| !CheckExtStabilizerState(qubitRegister, *clone,
			"Clone test copied current state B"))
		return false;

	// Restoring the clone must use its independent copy of snapshot A and must
	// not change the original, which remains in state B.
	clone->RestoreState();
	if (!CheckExtStabilizerState(*savedStatevector, *clone,
		"Clone test copied saved state A")
		|| !CheckExtStabilizerState(qubitRegister, simulator,
			"Clone restore changed the original"))
		return false;

	// Mutate only the clone to state C, then restore only the original to A.
	auto cloneBranchStatevector = savedStatevector->Clone();
	ApplyExtStabilizerTestGate(*cloneBranchStatevector, *clone, 3, 2);
	ApplyExtStabilizerTestGate(*cloneBranchStatevector, *clone, 17, 1, 0, -0.27);
	ApplyExtStabilizerTestGate(*cloneBranchStatevector, *clone, 15, 0, 0, 0.19);
	if (!CheckExtStabilizerState(*cloneBranchStatevector, *clone,
		"Clone test independently mutated state C")
		|| !CheckExtStabilizerState(qubitRegister, simulator,
			"Clone mutation changed the original"))
		return false;

	simulator.RestoreState();
	if (!CheckExtStabilizerState(*savedStatevector, simulator,
		"Clone test original saved state A")
		|| !CheckExtStabilizerState(*cloneBranchStatevector, *clone,
			"Original restore changed the clone"))
		return false;

	// A simulator clone also duplicates the random stream, so identical future
	// measurement sequences produce identical outcomes.
	QC::ExtendedStabilizer stochasticSimulator(32);
	for (size_t qubit = 0; qubit < 32; ++qubit)
		stochasticSimulator.ApplyH(qubit);
	auto stochasticClone = stochasticSimulator.Clone();
	for (size_t qubit = 0; qubit < 32; ++qubit)
		if (stochasticSimulator.Measure(qubit)
			!= stochasticClone->Measure(qubit))
		{
			std::cout << "\nClone did not preserve the measurement random stream"
				<< std::endl;
			return false;
		}

	std::cout << "Success" << std::endl;
	return true;
}


static bool TestExtStabilizerCanonicalRotations()
{
	std::cout << "\nExtended Stabilizer Clifford-equivalent rotation tests" << std::endl;

	const double pi = std::acos(-1.0);
	const std::vector<long long> turns{ -8, -5, -4, -3, -2, -1,
		1, 2, 3, 4, 5, 8 };
	for (int axis = 0; axis < 3; ++axis)
		for (const long long quarterTurns : turns)
		{
			QC::QubitRegister<> qubitRegister(3);
			QC::ExtendedStabilizer simulator(3);
			ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
			ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 1, 0);
			ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 2, 0, 0.37);
			ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 0, 0, -0.29);

			const size_t componentsBefore =
				simulator.GetFrames().front().GetFrameSize();
			const int gateCode = 15 + axis;
			const size_t qubit = static_cast<size_t>(axis);
			const double angle = static_cast<double>(quarterTurns) * pi / 2.0;
			ApplyExtStabilizerTestGate(qubitRegister, simulator,
				gateCode, qubit, 0, angle);

			const auto& frame = simulator.GetFrames().front();
			if (frame.GetFrameSize() != componentsBefore
				|| !frame.cliffordBasis.IsConsistent()
				|| !CheckExtStabilizerState(qubitRegister, simulator,
					"Clifford-equivalent rotation"))
			{
				std::cout << "\nCanonical rotation failed for axis " << axis
					<< " and " << quarterTurns << " quarter turns" << std::endl;
				return false;
			}
		}

	// The moving basis uses C_P = exp(i*pi/4) R_P(pi/2), so the exact
	// canonical-rotation global phase is retained in the sole coefficient.
	for (const long long quarterTurns : turns)
	{
		QC::ExtendedStabilizer simulator(1);
		simulator.ApplyRx(0, static_cast<double>(quarterTurns) * pi / 2.0);
		const auto expected = std::polar(1.0,
			-static_cast<double>(quarterTurns) * pi / 4.0);
		const auto& frame = simulator.GetFrames().front();
		if (frame.GetFrameSize() != 1
			|| std::abs(frame.amplitudes.front() - expected) > 1E-12)
		{
			std::cout << "\nCanonical rotation global phase was not retained for "
				<< quarterTurns << " quarter turns" << std::endl;
			return false;
		}
	}

	// An angle merely close to pi/2 remains a genuine non-Clifford rotation.
	QC::ExtendedStabilizer nearbyRotation(1);
	nearbyRotation.ApplyRx(0, pi / 2.0 + 1E-10);
	if (nearbyRotation.GetFrames().front().GetFrameSize() != 2)
	{
		std::cout << "\nA non-canonical nearby angle was incorrectly snapped" << std::endl;
		return false;
	}

	std::cout << "Success" << std::endl;
	return true;
}


static bool TestExtStabilizerApproximationPolicy()
{
	std::cout << "\nExtended Stabilizer approximation-policy tests" << std::endl;

	auto frameNorm = [](const QC::ExtendedFrame& frame)
	{
		double norm = 0.0;
		for (const auto& amplitude : frame.amplitudes)
			norm += std::norm(amplitude);
		return norm;
	};

	// Exact is the default and must not silently remove a nonzero coefficient,
	// even when it is far below the former hard-coded tolerance.
	{
		QC::ExtendedStabilizer simulator(1);
		simulator.ApplyRx(0, 1E-16);
		const auto& policy = simulator.GetApproximationPolicy();
		const auto& statistics = simulator.GetApproximationStatistics();
		if (policy.mode != QC::ExtendedStabilizerApproximationMode::Exact
			|| policy.amplitudeTolerance != 0.0 || policy.maxComponents != 0
			|| simulator.GetFrames().front().GetFrameSize() != 2
			|| statistics.discardedComponents != 0
			|| statistics.cumulativeDiscardedWeight != 0.0
			|| statistics.traceDistanceErrorBound != 0.0)
		{
			std::cout << "\nExact mode silently approximated a tiny rotation" << std::endl;
			return false;
		}
	}

	// Tolerance is expressed in normalized amplitude units. Pruning records the
	// discarded squared norm, exposes its trace-distance contribution, and
	// leaves the retained state normalized.
	{
		const double angle = 1E-4;
		QC::ExtendedStabilizer simulator(1,
			QC::ExtendedStabilizerApproximationPolicy::Approximate(1E-3));
		simulator.ApplyRy(0, angle);
		const auto& frame = simulator.GetFrames().front();
		const auto& statistics = simulator.GetApproximationStatistics();
		const double expectedDiscardedWeight =
			std::pow(std::sin(angle / 2.0), 2);
		if (frame.GetFrameSize() != 1
			|| !approxEqual(frameNorm(frame), 1.0, 1E-12)
			|| !approxEqual(simulator.GetQubitProbability(0), 0.0, 1E-12)
			|| statistics.discardedComponents != 1
			|| statistics.pruningEvents != 1
			|| !approxEqual(statistics.cumulativeDiscardedWeight,
				expectedDiscardedWeight, 1E-16)
			|| !approxEqual(statistics.traceDistanceErrorBound,
				std::sqrt(expectedDiscardedWeight), 1E-12))
		{
			std::cout << "\nTolerance pruning or error accounting failed" << std::endl;
			return false;
		}
	}

	// A hard cap keeps the largest components deterministically and remains
	// normalized. Projector-probability errors must fit the exposed bound.
	{
		QC::QubitRegister<> qubitRegister(2);
		QC::ExtendedStabilizer simulator(2,
			QC::ExtendedStabilizerApproximationPolicy::Approximate(0.0, 2));
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, 0.7);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 1, 0, 0.9);
		const auto& frame = simulator.GetFrames().front();
		const auto& statistics = simulator.GetApproximationStatistics();
		if (frame.GetFrameSize() != 2
			|| !approxEqual(frameNorm(frame), 1.0, 1E-12)
			|| statistics.discardedComponents != 2
			|| statistics.pruningEvents != 1
			|| statistics.traceDistanceErrorBound <= 0.0
			|| std::abs(simulator.GetQubitProbability(0)
				- qubitRegister.GetQubitProbability(0))
				> statistics.traceDistanceErrorBound + 1E-12
			|| std::abs(simulator.GetQubitProbability(1)
				- qubitRegister.GetQubitProbability(1))
				> statistics.traceDistanceErrorBound + 1E-12)
		{
			std::cout << "\nMaximum-component approximation failed" << std::endl;
			return false;
		}
	}

	// Even an intentionally extreme tolerance retains and normalizes the largest
	// component instead of deleting the state.
	{
		QC::ExtendedStabilizer simulator(1,
			QC::ExtendedStabilizerApproximationPolicy::Approximate(2.0));
		simulator.ApplyRx(0, 0.7);
		if (simulator.GetFrames().front().GetFrameSize() != 1
			|| !approxEqual(frameNorm(simulator.GetFrames().front()), 1.0, 1E-12)
			|| simulator.GetApproximationStatistics().discardedComponents != 1)
		{
			std::cout << "\nApproximation removed the final component" << std::endl;
			return false;
		}
	}

	// Policy changes prune the current state. Statistics belong to the saved
	// state, clones copy them, and Reset clears them while retaining the policy.
	{
		QC::ExtendedStabilizer simulator(3);
		simulator.ApplyRx(0, 0.61);
		simulator.ApplyRy(1, 0.73);
		simulator.SetApproximationPolicy(
			QC::ExtendedStabilizerApproximationPolicy::Approximate(0.0, 2));
		if (simulator.GetFrames().front().GetFrameSize() != 2)
			return false;
		simulator.SaveState();
		const auto savedStatistics = simulator.GetApproximationStatistics();
		auto clone = simulator.Clone();
		simulator.ApplyRx(2, 0.47);
		if (simulator.GetApproximationStatistics().pruningEvents
			<= savedStatistics.pruningEvents)
		{
			std::cout << "\nApproximation statistics did not accumulate" << std::endl;
			return false;
		}
		simulator.RestoreState();
		if (simulator.GetApproximationStatistics().pruningEvents
			!= savedStatistics.pruningEvents
			|| clone->GetApproximationStatistics().pruningEvents
				!= savedStatistics.pruningEvents)
		{
			std::cout << "\nSave, restore, or clone lost approximation statistics"
				<< std::endl;
			return false;
		}
		simulator.Reset(2);
		if (simulator.GetApproximationPolicy().mode
				!= QC::ExtendedStabilizerApproximationMode::Approximate
			|| simulator.GetApproximationPolicy().maxComponents != 2
			|| simulator.GetApproximationStatistics().pruningEvents != 0)
		{
			std::cout << "\nReset did not preserve policy and clear statistics"
				<< std::endl;
			return false;
		}
	}

	// A saved state owns its execution policy as well as its frame data. Changing
	// policy after SaveState must not make RestoreState violate that snapshot.
	{
		QC::ExtendedStabilizer simulator(2);
		simulator.ApplyRx(0, 0.61);
		simulator.ApplyRy(1, 0.73);
		simulator.SaveState();
		simulator.SetApproximationPolicy(
			QC::ExtendedStabilizerApproximationPolicy::Approximate(0.0, 2));
		if (simulator.GetFrames().front().GetFrameSize() != 2)
			return false;
		simulator.RestoreState();
		if (simulator.GetApproximationPolicy().mode
				!= QC::ExtendedStabilizerApproximationMode::Exact
			|| simulator.GetFrames().front().GetFrameSize() != 4
			|| simulator.GetApproximationStatistics().pruningEvents != 0)
		{
			std::cout << "\nRestore did not recover the saved approximation policy"
				<< std::endl;
			return false;
		}
	}

	// Repeated off-diagonal measurements exercise the reusable collapse buffers.
	{
		QC::ExtendedStabilizer simulator(3,
			QC::ExtendedStabilizerApproximationPolicy::Approximate(1E-8, 3));
		simulator.ApplyRx(0, 0.37);
		simulator.ApplyRy(1, -0.41);
		simulator.ApplyRx(2, 0.29);
		simulator.ApplyH(0);
		simulator.ApplyCX(1, 0);
		if (simulator.GetApproximationErrorBound() <= 0.0
			|| simulator.GetApproximationErrorBound() >= 1.0)
		{
			std::cout << "\nApproximation did not expose a useful pre-measurement bound"
				<< std::endl;
			return false;
		}
		for (size_t qubit = 0; qubit < 3; ++qubit)
		{
			const bool outcome = simulator.Measure(qubit);
			const auto& frame = simulator.GetFrames().front();
			if (frame.GetFrameSize() > 3
				|| !approxEqual(frameNorm(frame), 1.0, 1E-12)
				|| !approxEqual(simulator.GetQubitProbability(qubit),
					outcome ? 1.0 : 0.0, 1E-12)
				|| simulator.GetApproximationErrorBound() != 1.0
				|| simulator.Measure(qubit) != outcome)
			{
				std::cout << "\nApproximate measurement collapse failed" << std::endl;
				return false;
			}
		}
	}

	// Invalid or meaningless policies are rejected.
	for (const auto& invalidPolicy : {
		QC::ExtendedStabilizerApproximationPolicy{
			QC::ExtendedStabilizerApproximationMode::Exact, 1E-3, 0 },
		QC::ExtendedStabilizerApproximationPolicy::Approximate(0.0, 0),
		QC::ExtendedStabilizerApproximationPolicy::Approximate(-1.0, 1),
		QC::ExtendedStabilizerApproximationPolicy::Approximate(
			std::numeric_limits<double>::quiet_NaN(), 1) })
	{
		bool rejected = false;
		try
		{
			QC::ExtendedStabilizer simulator(1, invalidPolicy);
		}
		catch (const std::invalid_argument&)
		{
			rejected = true;
		}
		if (!rejected)
		{
			std::cout << "\nAn invalid approximation policy was accepted" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;
	return true;
}


static bool TestExtStabilizerLogicalBasisInvariance()
{
	std::cout << "\nExtended Stabilizer logical-basis invariance tests" << std::endl;

	QC::QubitRegister<> qubitRegister(4);
	QC::ExtendedStabilizer simulator(4);

	// Use every rotation axis to create several coherent components in a
	// nontrivial Clifford basis. Component amplitudes and logical labels must
	// then remain exactly unchanged while later Clifford gates rotate that basis.
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 2);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 3, 2);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 1, 1);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, 0.37);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 1, 0, -0.41);
	ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 2, 0, 0.29);

	const auto& preparedFrame = simulator.GetFrames().front();
	if (preparedFrame.GetFrameSize() <= 1)
	{
		std::cout << "\nLogical-basis invariance test did not create multiple components" << std::endl;
		return false;
	}

	const auto amplitudesBeforeCliffords = preparedFrame.amplitudes;
	const auto signsBeforeCliffords = preparedFrame.signs;

	struct CliffordGate {
		int code;
		size_t qubit1;
		size_t qubit2;
	};
	const std::vector<CliffordGate> cliffordGates = {
		{ 0, 3, 0 },  // H
		{ 1, 0, 0 },  // S
		{ 3, 1, 0 },  // X
		{ 4, 2, 0 },  // Y
		{ 5, 3, 0 },  // Z
		{ 9, 2, 0 },  // CX
		{ 11, 1, 3 }, // CZ
		{ 9, 0, 2 },
		{ 0, 1, 0 },
		{ 1, 2, 0 }
	};

	for (size_t gateIndex = 0; gateIndex < cliffordGates.size(); ++gateIndex)
	{
		const auto& gate = cliffordGates[gateIndex];
		ApplyExtStabilizerTestGate(qubitRegister, simulator,
			gate.code, gate.qubit1, gate.qubit2);

		const auto& frame = simulator.GetFrames().front();
		if (frame.amplitudes != amplitudesBeforeCliffords
			|| frame.signs != signsBeforeCliffords)
		{
			std::cout << "\nClifford gate " << gateIndex
				<< " changed logical labels or amplitudes" << std::endl;
			return false;
		}
		if (!frame.cliffordBasis.IsConsistent())
		{
			std::cout << "\nClifford gate " << gateIndex
				<< " corrupted the logical basis map" << std::endl;
			return false;
		}
	}

	if (!CheckExtStabilizerState(qubitRegister, simulator,
		"Logical-basis invariance final state"))
		return false;

	std::cout << "Success" << std::endl;
	return true;
}


static bool CheckExtStabilizerSingleQubitValues(const QC::ExtendedStabilizer& simulator, const std::string& context,
	double expectedProbability, double expectedX, double expectedY, double expectedZ)
{
	const double actualProbability = simulator.GetQubitProbability(0);
	const double actualX = simulator.ExpectationValue("X");
	const double actualY = simulator.ExpectationValue("Y");
	const double actualZ = simulator.ExpectationValue("Z");

	if (!approxEqual(expectedProbability, actualProbability, 1E-9)
		|| !approxEqual(expectedX, actualX, 1E-9)
		|| !approxEqual(expectedY, actualY, 1E-9)
		|| !approxEqual(expectedZ, actualZ, 1E-9))
	{
		std::cout << "\n" << context << ": analytic rotation result mismatch"
			<< "\nExpected P(1), <X>, <Y>, <Z>: " << expectedProbability << ", " << expectedX << ", " << expectedY << ", " << expectedZ
			<< "\nActual P(1), <X>, <Y>, <Z>: " << actualProbability << ", " << actualX << ", " << actualY << ", " << actualZ << std::endl;
		return false;
	}

	return true;
}


static bool CheckExtStabilizerMeasurement(QC::QubitRegister<>& qubitRegister, QC::ExtendedStabilizer& simulator,
	const std::string& context)
{
	const double expectedProbability = qubitRegister.GetQubitProbability(0);
	if (!approxEqual(simulator.GetQubitProbability(0), expectedProbability, 1E-9))
	{
		std::cout << "\n" << context << ": probability differs before measurement, statevector " << expectedProbability
			<< ", extended stabilizer " << simulator.GetQubitProbability(0) << std::endl;
		return false;
	}

	constexpr size_t nrShots = 3000;
	size_t countOne = 0;
	simulator.SaveState();
	for (size_t shot = 0; shot < nrShots; ++shot)
	{
		simulator.RestoreState();
		const bool first = simulator.Measure(0);
		const size_t collapsedComponents =
			simulator.GetFrames().front().GetFrameSize();
		const bool basisIsConsistent =
			simulator.GetFrames().front().cliffordBasis.IsConsistent();
		const bool second = simulator.Measure(0);
		if (first != second || collapsedComponents != 1 || !basisIsConsistent
			|| !approxEqual(simulator.GetQubitProbability(0), first ? 1.0 : 0.0, 1E-9))
		{
			std::cout << "\n" << context
				<< ": measurement did not collapse to one consistent logical component"
				<< std::endl;
			return false;
		}
		if (first) ++countOne;
	}

	simulator.RestoreState();
	if (!approxEqual(simulator.GetQubitProbability(0), expectedProbability, 1E-9))
	{
		std::cout << "\n" << context << ": RestoreState did not restore the measurement probability" << std::endl;
		return false;
	}

	const double measuredProbability = static_cast<double>(countOne) / nrShots;
	if (std::abs(measuredProbability - expectedProbability) > 0.045)
	{
		std::cout << "\n" << context << ": measurement frequency mismatch, statevector " << expectedProbability
			<< ", extended stabilizer " << measuredProbability << std::endl;
		return false;
	}

	return true;
}


static bool TestExtStabilizerNonCliffordRotations()
{
	std::cout << "\nExtended Stabilizer non-Clifford rotation tests" << std::endl;

	const double theta = M_PI / 3.0;
	const double sinTheta = std::sin(theta);
	const double cosTheta = std::cos(theta);
	const double probabilityOne = std::pow(std::sin(theta / 2.0), 2.0);

	// Rx(theta)|0>: the Y expectation checks the phase of the split components.
	{
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, theta);

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Rx analytic test", probabilityOne, 0.0, -sinTheta, cosTheta)
			|| !CheckExtStabilizerState(qubitRegister, simulator, "Rx statevector test"))
			return false;
	}

	// Ry(theta)|0> checks the real relative amplitudes.
	{
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 0, 0, theta);

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Ry analytic test", probabilityOne, sinTheta, 0.0, cosTheta)
			|| !CheckExtStabilizerState(qubitRegister, simulator, "Ry statevector test"))
			return false;
	}

	// Rz(theta)|+> keeps computational-basis probabilities unchanged, so X and Y
	// expectations are needed to detect the relative phase and its sign.
	{
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 0, 0, theta);

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Rz analytic test", 0.5, cosTheta, sinTheta, 0.0)
			|| !CheckExtStabilizerState(qubitRegister, simulator, "Rz statevector test"))
			return false;
	}

	// Repeated rotations around one axis must combine amplitudes belonging to the
	// same stabilizer components exactly as a single rotation by the summed angle.
	{
		constexpr double angle1 = 0.37;
		constexpr double angle2 = -0.82;
		const double totalAngle = angle1 + angle2;

		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, angle1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, angle2);

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Rx composition test",
			std::pow(std::sin(totalAngle / 2.0), 2.0), 0.0, -std::sin(totalAngle), std::cos(totalAngle))
			|| !CheckExtStabilizerState(qubitRegister, simulator, "Rx composition statevector test"))
			return false;
	}

	// Two Rx(pi/2) gates produce |1> up to global phase. This forces exact
	// cancellation when duplicate component sign patterns are combined.
	{
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, M_PI / 2.0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, M_PI / 2.0);

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Rx exact cancellation test", 1.0, 0.0, 0.0, -1.0)
			|| !CheckExtStabilizerState(qubitRegister, simulator, "Rx exact cancellation statevector test"))
			return false;
	}

	// Clifford gates after a non-Clifford split rotate the shared basis while
	// preserving every component's logical label and relative coefficient.
	{
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 0, 0, 0.41);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 1, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 0, 0, -0.29);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 2, 0);

		if (!CheckExtStabilizerState(qubitRegister, simulator, "Rotation followed by Clifford gates test"))
			return false;
	}

	// This state is deliberately expressed in an X stabilizer basis while its Z
	// measurement probability is not 1/2; it exercises cross-component coherence.
	{
		const double coherenceAngle = M_PI / 6.0;
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 0, 0, coherenceAngle);

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Coherent frame test",
			0.5 * (1.0 + std::sin(coherenceAngle)), std::cos(coherenceAngle), 0.0, -std::sin(coherenceAngle))
			|| !CheckExtStabilizerState(qubitRegister, simulator, "Coherent frame statevector test"))
			return false;
	}

	// Entanglement plus rotations about all three axes exercises component signs
	// through subsequent two-qubit Clifford operations.
	{
		QC::QubitRegister<> qubitRegister(2);
		QC::ExtendedStabilizer simulator(2);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 0, 1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 0, 0, 0.37);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 1, 0, -0.61);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 1, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, 0.23);

		if (!CheckExtStabilizerState(qubitRegister, simulator, "Entangled rotation test"))
			return false;
	}

	// Small reproducible mixed circuits provide broader coverage without enabling
	// non-Clifford gates in the much larger random Clifford regression below.
	std::mt19937 fixedGenerator(0xE57AB1E);
	std::uniform_int_distribution<int> singleCliffordDistribution(0, 8);
	std::uniform_int_distribution<int> twoQubitCliffordDistribution(9, 14);
	std::uniform_int_distribution<int> rotationDistribution(15, 17);
	std::uniform_real_distribution<double> angleDistribution(-M_PI, M_PI);

	for (size_t nrQubits = 1; nrQubits <= 4; ++nrQubits)
	{
		std::uniform_int_distribution<int> qubitDistribution(0, static_cast<int>(nrQubits) - 1);
		for (size_t circuitIndex = 0; circuitIndex < 6; ++circuitIndex)
		{
			QC::QubitRegister<> qubitRegister(nrQubits);
			QC::ExtendedStabilizer simulator(nrQubits);

			for (size_t gateIndex = 0; gateIndex < 12; ++gateIndex)
			{
				int code = 0;
				if (gateIndex % 2 == 0)
					code = rotationDistribution(fixedGenerator);
				else if (nrQubits > 1 && gateIndex % 4 == 3)
					code = twoQubitCliffordDistribution(fixedGenerator);
				else
					code = singleCliffordDistribution(fixedGenerator);

				const size_t qubit1 = static_cast<size_t>(qubitDistribution(fixedGenerator));
				const size_t qubit2 = nrQubits > 1 ? (qubit1 + 1) % nrQubits : 0;
				const double angle = code >= 15 ? angleDistribution(fixedGenerator) : 0.0;
				ApplyExtStabilizerTestGate(qubitRegister, simulator, code, qubit1, qubit2, angle);
				const std::string gateContext = "Fixed mixed circuit " + std::to_string(circuitIndex)
					+ " on " + std::to_string(nrQubits) + " qubits after gate "
					+ std::to_string(gateIndex) + " (code " + std::to_string(code) + ")";
				if (!simulator.GetFrames().front().cliffordBasis.IsConsistent())
				{
					std::cout << "\n" << gateContext << ": Clifford basis maps disagree" << std::endl;
					return false;
				}
				if (!CheckExtStabilizerState(qubitRegister, simulator, gateContext))
					return false;
			}

			const std::string context = "Fixed mixed circuit " + std::to_string(circuitIndex)
				+ " on " + std::to_string(nrQubits) + " qubits";
			if (!CheckExtStabilizerState(qubitRegister, simulator, context))
				return false;
		}
	}

	// A stabilizer measurement updates the paired Clifford basis directly, so
	// coherent work can resume without synthesizing a preparation circuit.
	{
		QC::ExtendedStabilizer simulator(1);
		simulator.ApplyH(0);
		const bool outcome = simulator.Measure(0);
		if (!simulator.GetFrames().front().cliffordBasis.IsConsistent())
		{
			std::cout << "\nOne-qubit measurement corrupted the Clifford basis map" << std::endl;
			return false;
		}
		simulator.ApplyH(0);
		simulator.ApplyRz(0, theta);
		const double sign = outcome ? -1.0 : 1.0;

		if (!CheckExtStabilizerSingleQubitValues(simulator, "Measurement followed by rotation test",
			0.5, sign * cosTheta, sign * sinTheta, 0.0))
			return false;
	}

	// Exercise the same direct measurement-basis update on an entangled state.
	// Measuring either half of a Bell pair leaves |00> or |11>; the following
	// H/Rz sequence therefore has a simple outcome-conditioned oracle.
	{
		QC::ExtendedStabilizer simulator(2);
		simulator.ApplyH(0);
		simulator.ApplyCX(1, 0);
		const bool outcome = simulator.Measure(0);
		if (!simulator.GetFrames().front().cliffordBasis.IsConsistent())
		{
			std::cout << "\nEntangled measurement corrupted the Clifford basis map" << std::endl;
			return false;
		}
		simulator.ApplyH(1);
		simulator.ApplyRz(1, theta);
		const double sign = outcome ? -1.0 : 1.0;

		if (!approxEqual(simulator.GetQubitProbability(0), outcome ? 1.0 : 0.0, 1E-9)
			|| !approxEqual(simulator.GetQubitProbability(1), 0.5, 1E-9)
			|| !approxEqual(simulator.ExpectationValue("ZI"), sign, 1E-9)
			|| !approxEqual(simulator.ExpectationValue("IX"), sign * cosTheta, 1E-9)
			|| !approxEqual(simulator.ExpectationValue("IY"), sign * sinTheta, 1E-9))
		{
			std::cout << "\nEntangled measurement followed by rotation test failed" << std::endl;
			return false;
		}
	}

	// Exercise a less structured measurement rebase: condition the statevector
	// on the simulator's random outcome, then compare a four-qubit state after a
	// longer mixture of Clifford gates and rotations about every axis.
	{
		QC::QubitRegister<> qubitRegister(4);
		QC::ExtendedStabilizer simulator(4);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 1, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 2);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 3, 2);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 11, 1, 2);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 1, 3);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 9, 2, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 17, 1, 0, 0.37);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 3, 0, -0.29);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, 0.41);

		const double expectedMeasurementProbability =
			qubitRegister.GetQubitProbability(1);
		if (!approxEqual(simulator.GetQubitProbability(1),
			expectedMeasurementProbability, 1E-9)
			|| expectedMeasurementProbability < 0.05
			|| expectedMeasurementProbability > 0.95)
		{
			std::cout << "\nFour-qubit rotated measurement test did not prepare a valid random measurement" << std::endl;
			return false;
		}

		const bool outcome = simulator.Measure(1);
		auto projectedState = qubitRegister.getRegisterStorage();
		for (size_t basisState = 0;
			basisState < static_cast<size_t>(projectedState.size()); ++basisState)
			if (static_cast<bool>((basisState >> 1) & 1ULL) != outcome)
				projectedState(static_cast<Eigen::Index>(basisState)) = 0.0;
		qubitRegister.setRegisterStorage(projectedState);

		if (!simulator.GetFrames().front().cliffordBasis.IsConsistent()
			|| !CheckExtStabilizerState(qubitRegister, simulator,
				"Four-qubit state immediately after measurement"))
		{
			std::cout << "\nFour-qubit measurement rebase failed" << std::endl;
			return false;
		}

		struct PostMeasurementGate {
			int code;
			size_t qubit1;
			size_t qubit2;
			double angle;
		};
		const std::vector<PostMeasurementGate> postMeasurementGates = {
			{ 0, 3, 0, 0.0 },
			{ 9, 2, 1, 0.0 },
			{ 15, 0, 0, 0.31 },
			{ 16, 2, 0, -0.47 },
			{ 17, 3, 0, 0.28 },
			{ 11, 0, 3, 0.0 },
			{ 6, 1, 0, 0.0 },
			{ 15, 1, 0, -0.19 },
			{ 16, 3, 0, 0.22 },
			{ 17, 0, 0, -0.41 },
			{ 12, 0, 2, 0.0 }
		};

		for (size_t gateIndex = 0; gateIndex < postMeasurementGates.size(); ++gateIndex)
		{
			const auto& gate = postMeasurementGates[gateIndex];
			ApplyExtStabilizerTestGate(qubitRegister, simulator,
				gate.code, gate.qubit1, gate.qubit2, gate.angle);
			const std::string context = "Four-qubit post-measurement circuit after gate "
				+ std::to_string(gateIndex);
			if (!simulator.GetFrames().front().cliffordBasis.IsConsistent()
				|| !CheckExtStabilizerState(qubitRegister, simulator, context))
				return false;
		}
	}

	// Measurement probabilities can contain interference between frame components.
	// Check their distribution, collapse, repeatability, and SaveState/RestoreState.
	{
		const double measurementAngle = M_PI / 6.0;
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 0, 0);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 16, 0, 0, measurementAngle);

		if (!approxEqual(qubitRegister.GetQubitProbability(0), 0.75, 1E-9)
			|| !CheckExtStabilizerMeasurement(qubitRegister, simulator, "Coherent non-Clifford measurement test"))
		{
			std::cout << "\nInvalid result in coherent non-Clifford measurement test" << std::endl;
			return false;
		}
	}

	// Also exercise the simpler biased measurement path produced directly by Rx.
	{
		QC::QubitRegister<> qubitRegister(1);
		QC::ExtendedStabilizer simulator(1);
		ApplyExtStabilizerTestGate(qubitRegister, simulator, 15, 0, 0, M_PI / 3.0);

		if (!approxEqual(qubitRegister.GetQubitProbability(0), 0.25, 1E-9)
			|| !CheckExtStabilizerMeasurement(qubitRegister, simulator, "Rx non-Clifford measurement test"))
		{
			std::cout << "\nInvalid result in Rx non-Clifford measurement test" << std::endl;
			return false;
		}
	}

	// RestoreState without a preceding SaveState is a no-op.
	{
		QC::ExtendedStabilizer simulator(1);
		simulator.ApplyH(0);
		simulator.RestoreState();
		if (!CheckExtStabilizerSingleQubitValues(simulator, "Restore without save test", 0.5, 1.0, 0.0, 0.0))
			return false;
	}

	// Reset starts a new register and must also discard an older snapshot.
	{
		QC::ExtendedStabilizer simulator(1);
		simulator.ApplyX(0);
		simulator.SaveState();
		simulator.Reset(2);
		simulator.RestoreState();
		if (simulator.GetNrQubits() != 2
			|| !approxEqual(simulator.GetQubitProbability(0), 0.0, 1E-12)
			|| !approxEqual(simulator.GetQubitProbability(1), 0.0, 1E-12))
		{
			std::cout << "\nReset did not invalidate an older saved state" << std::endl;
			return false;
		}
	}


	// Exercise packed-word boundaries, including a measurement at the highest
	// qubit and Clifford updates spanning two uint64_t words.
	for (const size_t nrQubits : { size_t(63), size_t(64), size_t(65) })
	{
		// Put a non-Clifford component label on the final bit as well, so the
		// packed Pauli view and component-key packing cross the word boundary.
		QC::ExtendedStabilizer rotatedSimulator(nrQubits);
		const size_t last = nrQubits - 1;
		rotatedSimulator.ApplyRx(last, theta);
		if (rotatedSimulator.GetFrames().front().GetFrameSize() != 2
			|| !approxEqual(rotatedSimulator.GetQubitProbability(last),
				probabilityOne, 1E-12)
			|| !rotatedSimulator.GetFrames().front().cliffordBasis.IsConsistent())
		{
			std::cout << "\nPacked non-Clifford boundary test failed for "
				<< nrQubits << " qubits" << std::endl;
			return false;
		}

		QC::ExtendedStabilizer simulator(nrQubits);
		const size_t middle = nrQubits / 2;
		simulator.ApplyH(last);
		simulator.Measure(last);
		simulator.ApplyH(0);
		simulator.ApplyCX(last, 0);
		simulator.ApplyS(last);
		simulator.ApplyCX(0, last);
		simulator.ApplyY(middle);
		simulator.ApplyCZ(last, middle);
		if (!simulator.GetFrames().front().cliffordBasis.IsConsistent())
		{
			std::cout << "\nPacked Clifford-map boundary test failed for "
				<< nrQubits << " qubits" << std::endl;
			return false;
		}
	}

	// Invalid angles must be rejected before they can corrupt the live state.
	{
		QC::ExtendedStabilizer simulator(1);
		simulator.ApplyRy(0, 0.37);
		const double probabilityBefore = simulator.GetQubitProbability(0);
		const double xBefore = simulator.ExpectationValue("X");
		bool rejected = false;
		try
		{
			simulator.ApplyRz(0, std::numeric_limits<double>::quiet_NaN());
		}
		catch (const std::invalid_argument&)
		{
			rejected = true;
		}

		if (!rejected || !approxEqual(simulator.GetQubitProbability(0), probabilityBefore, 1E-12)
			|| !approxEqual(simulator.ExpectationValue("X"), xBefore, 1E-12))
		{
			std::cout << "\nInvalid rotation angle handling test failed" << std::endl;
			return false;
		}
	}

	std::cout << "Success" << std::endl;
	return true;
}

void ExecuteCircuit(QC::QubitRegister<>& qubitRegister,
	QC::ExtendedStabilizer& extstabSim, const std::vector<int>& gates,
	const std::vector<size_t>& qubits1, const std::vector<size_t>& qubits2)
{
	for (int j = 0; j < static_cast<int>(gates.size()); ++j)
	{
		auto gateptr = GetGate(gates[j]);
		qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
		ApplyGate(extstabSim, gates[j], static_cast<int>(qubits1[j]),
			static_cast<int>(qubits2[j]));
	}
}


static bool TestExtStabilizerMeasurements()
{
	std::cout << "\nExtended Stabilizer measurement tests" << std::endl;

	// Exact probability and conditioned-state checks above carry the correctness
	// burden.  Keep this as a broad but inexpensive statistical smoke test.
	const size_t nrShots = 20000;
	const double errorThreshold = 0.03;
	const size_t nrTests = 4;
	const size_t maxQubits = 8;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(50, 100);

	for (size_t nrQubits = 2; nrQubits < maxQubits; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);

		for (size_t t = 0; t < nrTests; ++t)
		{
			const size_t nrGates = nrGatesDistr(gen);
			std::vector<int> gates(nrGates);
			std::vector<size_t> qubits1(nrGates);
			std::vector<size_t> qubits2(nrGates);

			ConstructCircuit(nrQubits, gates, qubits1, qubits2, gateDistr, qubitDistr);

			// test 1: repeated measurement on the same qubit must be consistent
			{
				QC::ExtendedStabilizer sim(nrQubits);
				QC::QubitRegister<> reg(nrQubits);
				ExecuteCircuit(reg, sim, gates, qubits1, qubits2);

				for (size_t q = 0; q < nrQubits; ++q)
				{
					const bool res1 = sim.Measure(q);
					const bool res2 = sim.Measure(q);
					if (res1 != res2)
					{
						std::cout << std::endl << "Repeated measurement inconsistency for qubit " << q << " with " << nrQubits << " qubits" << std::endl;
						return false;
					}
				}
			}

			// test 2: sampling comparison against statevector
			{
				std::unordered_map<size_t, size_t> stabResults;
				std::unordered_map<size_t, size_t> svResults;

				QC::QubitRegister<> reg(nrQubits);
				QC::ExtendedStabilizer sim(nrQubits);

				ExecuteCircuit(reg, sim, gates, qubits1, qubits2);

				svResults = reg.RepeatedMeasureUnordered(nrShots);

				sim.SaveState();
				for (size_t shot = 0; shot < nrShots; ++shot)
				{
					size_t stabVal = 0;
					for (size_t q = 0; q < nrQubits; ++q)
						if (sim.Measure(q)) stabVal |= 1ULL << q;

					++stabResults[stabVal];

					sim.RestoreState();
				}

				for (const auto& val : svResults)
				{
					const double svFreq = static_cast<double>(val.second) / nrShots;
					const double stabFreq = stabResults.count(val.first) ? static_cast<double>(stabResults[val.first]) / nrShots : 0.0;

					if (std::abs(svFreq - stabFreq) > errorThreshold)
					{
						std::cout << std::endl << "Measurement distribution mismatch for " << nrQubits << " qubits, state " << val.first
							<< ": statevector " << svFreq << ", stabilizer " << stabFreq << std::endl;
						std::cout << "Might fail due to randomness of measurements" << std::endl;
						return false;
					}
				}

				for (const auto& val : stabResults)
				{
					if (svResults.count(val.first)) continue;

					const double stabFreq = static_cast<double>(val.second) / nrShots;
					if (stabFreq > errorThreshold)
					{
						std::cout << std::endl << "Stabilizer produced state " << val.first << " with frequency " << stabFreq
							<< " but statevector never did, for " << nrQubits << " qubits" << std::endl;
						return false;
					}
				}
			}
		}
		std::cout << '.';
	}

	std::cout << "\nSuccess" << std::endl;
	return true;
}

bool TestExtStabilizer()
{
	std::cout << "\nExtended Stabilizer tests" << std::endl;

	if (!TestExtStabilizerNonCliffordRotations())
		return false;
	if (!TestExtStabilizerLargerMixedCircuits())
		return false;
	if (!TestExtStabilizerClone())
		return false;
	if (!TestExtStabilizerCanonicalRotations())
		return false;
	if (!TestExtStabilizerApproximationPolicy())
		return false;
	if (!TestExtStabilizerLogicalBasisInvariance())
		return false;

	std::cout << "\nExtended Stabilizer random Clifford-circuit tests" << std::endl;

	const size_t nrTests = 100;
	const size_t maxQubits = 20;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(50, 100);

	for (size_t nrQubits = 2; nrQubits < maxQubits; ++nrQubits)
	{
		std::uniform_int_distribution qubitDistr(0, static_cast<int>(nrQubits) - 1);

		for (size_t t = 0; t < nrTests; ++t)
		{
			std::unordered_map<size_t, int> results1;
			std::unordered_map<size_t, int> results2;

			// generate random gates, creating a circuit, then apply the random circuits on both simulators
			const size_t nrGates = nrGatesDistr(gen);
			std::vector<int> gates(nrGates);
			std::vector<size_t> qubits1(nrGates);
			std::vector<size_t> qubits2(nrGates);

			ConstructCircuit(nrQubits, gates, qubits1, qubits2, gateDistr, qubitDistr);

			QC::ExtendedStabilizer extstabSim(nrQubits);
			QC::QubitRegister qubitRegister(nrQubits);

			ExecuteCircuit(qubitRegister, extstabSim, gates, qubits1, qubits2);

			for (size_t q = 0; q < nrQubits; ++q)
			{
				double p1 = qubitRegister.GetQubitProbability(q);
				double p2 = extstabSim.GetQubitProbability(q);
				if (!approxEqual(p1, p2, 1E-5))
				{
					std::cout << std::endl << "Probabilities are not equal for statevector and stabilizer simulator for " << nrQubits << " qubits, values: " << p1 << ", " << p2 << std::endl;
					return false;
				}
			}

			std::vector<QC::Gates::AppliedGate<>> expGates;
			expGates.reserve(nrQubits);
			std::string pauliStr;

			ConstructPauliString(nrQubits, pauliStr, expGates);

			const auto exp1 = qubitRegister.ExpectationValue(expGates).real();
			const auto exp2 = extstabSim.ExpectationValue(pauliStr);
			if (!approxEqual(exp1, exp2, 1E-7))
			{
				std::cout << std::endl << "Expectation values are not equal for statevector and stabilizer simulator for " << nrQubits << " qubits, values: " << exp1 << ", " << exp2 << std::endl;

				std::cout << "Pauli string: " << pauliStr << std::endl;

				return false;
			}
		}
		std::cout << '.';
	}

	std::cout << "\nSuccess" << std::endl;

	if (!TestExtStabilizerMeasurements())
		return false;

	return true;
}
