#include "Clifford.h"
#include "Tests.h"

#include "QubitRegister.h"
#include "PauliStringXZCoeff.h"

#include <cstdlib>
#include <future>

std::shared_ptr<QC::Gates::QuantumGateWithOp<>> GetTwoQubitsGate(int code)
{
	switch (code)
	{
	case 9:
		return std::make_shared<QC::Gates::CNOTGate<>>();
	case 10:
		return std::make_shared<QC::Gates::ControlledYGate<>>();
	case 11:
		return std::make_shared<QC::Gates::ControlledZGate<>>();
	case 12:
		return std::make_shared<QC::Gates::SwapGate<>>();
	case 13:
		return std::make_shared<QC::Gates::iSwapGate<>>();
	case 14:
		return std::make_shared<QC::Gates::iSwapDagGate<>>();
	}

	return nullptr;
}


std::shared_ptr<QC::Gates::QuantumGateWithOp<>> GetGate(int code, double param)
{
	switch (code)
	{
	case 0:
		return std::make_shared<QC::Gates::HadamardGate<>>();
	case 1:
		return std::make_shared<QC::Gates::SGate<>>();
	case 2:
		return std::make_shared<QC::Gates::SDGGate<>>();
	case 3:
		return std::make_shared<QC::Gates::PauliXGate<>>();
	case 4:
		return std::make_shared<QC::Gates::PauliYGate<>>();
	case 5:
		return std::make_shared<QC::Gates::PauliZGate<>>();
	case 6:
		return std::make_shared<QC::Gates::SquareRootNOTGate<>>();
	case 7:
		return std::make_shared<QC::Gates::SquareRootNOTDagGate<>>();
	case 8:
		return std::make_shared<QC::Gates::HyGate<>>();
	// for non-clifford
	case 15:
		return std::make_shared<QC::Gates::RxGate<>>(param);
	case 16:
		return std::make_shared<QC::Gates::RyGate<>>(param);
	case 17:
		return std::make_shared<QC::Gates::RzGate<>>(param);
	default:
		return GetTwoQubitsGate(code);
	}

	return nullptr;
}

void ApplyTwoQubitsGate(QC::Clifford::StabilizerSimulator& simulator, int code, int qubit1, int qubit2)
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



void ApplyGate(QC::Clifford::StabilizerSimulator& simulator, int code, int qubit1, int qubit2)
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
	default:
		ApplyTwoQubitsGate(simulator, code, qubit1, qubit2);
		break;
	}
}

void PrintTwoQubitsGate(int code, int qubit1, int qubit2)
{
	switch (code)
	{
	case 9:
		std::cout << "CX " << qubit1 << " " << qubit2 << std::endl;
		break;
	case 10:
		std::cout << "CY " << qubit1 << " " << qubit2 << std::endl;
		break;
	case 11:
		std::cout << "CZ " << qubit1 << " " << qubit2 << std::endl;
		break;
	case 12:
		std::cout << "SWAP " << qubit1 << " " << qubit2 << std::endl;
		break;
	case 13:
		std::cout << "iSWAP " << qubit1 << " " << qubit2 << std::endl;
		break;
	case 14:
		std::cout << "iSWAPDag " << qubit1 << " " << qubit2 << std::endl;
		break;
	}
}

void PrintGate(int code, int qubit1, int qubit2)
{
	switch (code)
	{
	case 0:
		std::cout << "H " << qubit1 << std::endl;
		break;
	case 1:
		std::cout << "S " << qubit1 << std::endl;
		break;
	case 2:
		std::cout << "SDG " << qubit1 << std::endl;
		break;
	case 3:
		std::cout << "X " << qubit1 << std::endl;
		break;
	case 4:
		std::cout << "Y " << qubit1 << std::endl;
		break;
	case 5:
		std::cout << "Z " << qubit1 << std::endl;
		break;
	case 6:
		std::cout << "SX " << qubit1 << std::endl;
		break;
	case 7:
		std::cout << "SXDG " << qubit1 << std::endl;
		break;
	case 8:
		std::cout << "K " << qubit1 << std::endl;
		break;
	default:
		PrintTwoQubitsGate(code, qubit1, qubit2);
	}
}

void ConstructCircuit(size_t nrQubits, std::vector<int>& gates, std::vector<size_t>& qubits1, std::vector<size_t>& qubits2, std::uniform_int_distribution<int>& gateDistr, std::uniform_int_distribution<int>& qubitDistr)
{
	for (int i = 0; i < static_cast<int>(gates.size()); ++i)
	{
		gates[i] = gateDistr(gen);

		qubits1[i] = qubitDistr(gen);
		qubits2[i] = qubitDistr(gen);

		if (qubits2[i] == qubits1[i])
			qubits2[i] = (qubits1[i] + 1) % static_cast<int>(nrQubits);

		if (dist_bool(gen)) std::swap(qubits1[i], qubits2[i]);
	}
}

void ExecuteCircuit(size_t nrShots, size_t nrQubits, const std::vector<int>& gates, const std::vector<size_t>& qubits1, const std::vector<size_t>& qubits2, std::unordered_map<size_t, int>& results1, std::unordered_map<size_t, int>& results2)
{
	size_t remainingCounts = nrShots;
	size_t nrThreads = QC::QubitRegisterCalculator<>::GetNumberOfThreads();
	nrThreads = std::min(nrThreads, std::max<size_t>(remainingCounts, 1ULL));

	std::vector<std::future<void>> tasks(nrThreads);

	const size_t cntPerThread = static_cast<size_t>(ceil(static_cast<double>(remainingCounts) / nrThreads));

	std::mutex resultsMutex;

	for (size_t th = 0; th < nrThreads; ++th)
	{
		const size_t curCnt = std::min(cntPerThread, remainingCounts);
		remainingCounts -= curCnt;

		tasks[th] = std::async(std::launch::async, [&gates, &qubits1, &qubits2, &results1, &results2, curCnt, nrQubits, &resultsMutex]()
			{
				for (int i = 0; i < static_cast<int>(curCnt); ++i)
				{
					QC::QubitRegister qubitRegister(nrQubits);
					QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

					for (int j = 0; j < static_cast<int>(gates.size()); ++j)
					{
						ApplyGate(cliffordSim, gates[j], qubits1[j], qubits2[j]);
						const auto gateptr = GetGate(gates[j]);
						qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
					}

					// now do the measurements
					size_t val1 = 0;
					size_t val2 = 0;
					for (int q = 0; q < static_cast<int>(nrQubits); ++q)
					{
						val1 <<= 1;
						val2 <<= 1;

						if (cliffordSim.MeasureQubit(q)) val1 |= 1;
						if (qubitRegister.MeasureQubit(q)) val2 |= 1;
					}

					{
						const std::lock_guard lock(resultsMutex);
						++results1[val1];
						++results2[val2];
					}
				}
			});
	}

	for (size_t i = 0; i < nrThreads; ++i)
		tasks[i].get();
}

bool CheckProbability(QC::QubitRegister<>& qubitRegister, QC::Clifford::StabilizerSimulator& cliffordSim)
{
	const double probThreshold = 1E-10;
	size_t nrQubits = qubitRegister.getNrQubits();

	for (size_t q = 0; q < nrQubits; ++q)
	{
		const double prob1 = cliffordSim.GetQubitProbability(q);
		const double prob2 = qubitRegister.GetQubitProbability(q);

		if (std::abs(prob1 - prob2) > probThreshold)
		{
			std::cout << "\nFailed qubits probabilities" << std::endl;
			std::cout << "Probability statevector: " << prob2 << ", Probability stabilizer: " << prob1 << std::endl;
			return false;
		}
	}

	return true;
}

bool CheckAllStatesProbability(QC::QubitRegister<>& qubitRegister, QC::Clifford::StabilizerSimulator& cliffordSim)
{
	const double probThreshold = 1E-10;
	size_t nrQubits = qubitRegister.getNrQubits();
	const size_t nrStates = 1ULL << nrQubits;

	for (size_t state = 0; state < nrStates; ++state)
	{
		const double prob1 = cliffordSim.getBasisStateProbability(state);
		const double prob2 = qubitRegister.getBasisStateProbability(state);

		if (std::abs(prob1 - prob2) > probThreshold)
		{
			std::cout << "\nFailed states probabilities" << std::endl;
			std::cout << "Probability statevector: " << prob2 << ", Probability stabilizer: " << prob1 << ", State: " << state << std::endl;
			return false;
		}
	}

	return true;
}

bool CheckMeasurements(QC::Clifford::StabilizerSimulator& cliffordSim)
{
	size_t nrQubits = cliffordSim.getNrQubits();
	for (size_t q = 0; q < nrQubits; ++q)
	{
		const bool res1 = cliffordSim.MeasureQubit(q);
		const bool res2 = cliffordSim.MeasureQubit(q);
		if (res1 != res2)
		{
			std::cout << "\nFailed qubits measurements" << std::endl;
			std::cout << "Measurement 1: " << res2 << ", Measurement 2: " << res1 << ", for qubit: " << q << std::endl;
			return false;
		}
	}
	return true;
}

namespace {
	class CliffordRegressionSimulator : public QC::Clifford::StabilizerSimulator {
	public:
		using StabilizerSimulator::StabilizerSimulator;
		std::mt19937_64 RandomEngine() const { return gen; }
		bool DistributionIsPrepared() const { return (validDistributions & FullDistribution) != 0; }

		bool SameTableau(const CliffordRegressionSimulator& other) const
		{
			if (getNrQubits() != other.getNrQubits()) return false;
			for (size_t q = 0; q < getNrQubits(); ++q)
			{
				if (!(inverseZ[q] == other.inverseZ[q]) ||
					inverseZ[q].PhaseSign != other.inverseZ[q].PhaseSign ||
					!(inverseX[q] == other.inverseX[q]) ||
					inverseX[q].PhaseSign != other.inverseX[q].PhaseSign)
					return false;
			}
			return true;
		}
	};

	bool CliffordCheck(bool condition, const char* message)
	{
		if (!condition) std::cout << "Clifford regression failed: " << message << std::endl;
		return condition;
	}

	template<class Exception, class F> bool CliffordRejects(F&& f)
	{
		try { f(); }
		catch (const Exception&) { return true; }
		catch (...) { return false; }
		return false;
	}

	bool CliffordPackedRegression()
	{
		using Row = QC::Clifford::detail::TableauRow<false>;
		using ConstRow = QC::Clifford::detail::TableauRow<true>;
		static_assert(std::is_copy_constructible_v<Row> && std::is_copy_constructible_v<ConstRow>);
		static_assert(!std::is_copy_assignable_v<Row> && !std::is_move_assignable_v<Row>);
		static_assert(!std::is_copy_assignable_v<ConstRow> && !std::is_move_assignable_v<ConstRow>);
		QC::Clifford::detail::PackedTableau signedRows(2, 65);
		signedRows[0].X[64] = true;
		signedRows[1].CopyFrom(signedRows[0]);
		signedRows[1].PhaseSign = true;
		if (!CliffordCheck(!(signedRows[0] == signedRows[1]), "Row equality includes sign")) return false;
		signedRows[0].CopyFrom(signedRows[1]);
		if (!CliffordCheck(signedRows[0] == signedRows[1], "Explicit row copy includes sign")) return false;
		std::mt19937 rng(9173);
		for (size_t n : { 0, 1, 2, 63, 64, 65, 127, 128, 129, 1025 })
			for (unsigned repeat = 0; repeat < 128; ++repeat)
			{
				QC::Clifford::detail::PackedTableau rows(2, n);
				auto left = rows[0], right = rows[1];
				int phase = 0;
				std::vector<bool> expectedX(n), expectedZ(n);
				for (size_t q = 0; q < n; ++q)
				{
					const int x = rng() & 1, z = rng() & 1, rx = rng() & 1, rz = rng() & 1;
					left.X[q] = x; left.Z[q] = z; right.X[q] = rx; right.Z[q] = rz;
					expectedX[q] = x != rx; expectedZ[q] = z != rz;
					phase += QC::PauliStringXZWithSign::g(x, z, rx, rz);
				}
				left.PhaseSign = bool(rng() & 1); right.PhaseSign = bool(rng() & 1);
				phase += 2 * int(bool(left.PhaseSign) != bool(right.PhaseSign));
				const unsigned extra = (phase & 1) ? 1 : 0;
				phase += extra;
				left.Multiply(right, extra);
				if (!CliffordCheck(bool(left.PhaseSign) == ((phase & 2) != 0), "Packed product phase")) return false;
				for (size_t q = 0; q < n; ++q)
					if (!CliffordCheck(left.X[q] == expectedX[q] && left.Z[q] == expectedZ[q], "Packed product bits")) return false;
				left.Multiply(left);
				if (!CliffordCheck(!left.HasX() && !left.PhaseSign, "Aliased packed square")) return false;
				for (size_t w = 0; w < left.Words(); ++w)
					if (!CliffordCheck(left.Z.words[w] == 0, "Packed square is identity")) return false;
			}
		return true;
	}

	bool CliffordDistributionRegression()
	{
		CliffordRegressionSimulator cached(4);
		QC::QubitRegister<> cachedRef(4);
		cached.ApplyH(0); cached.ApplyCX(1, 0); cached.ApplyH(2);
		cachedRef.ApplyGate(QC::Gates::HadamardGate<>(), 0);
		cachedRef.ApplyGate(QC::Gates::CNOTGate<>(), 1, 0);
		cachedRef.ApplyGate(QC::Gates::HadamardGate<>(), 2);
		for (int gate : {1, 2, 5, 11, 3, 4})
			for (size_t q = 0; q < 4; ++q)
			{
				cached.AllProbabilities();
				ApplyGate(cached, gate, int(q), int((q + 1) % 4));
				cachedRef.ApplyGate(*GetGate(gate), q, (q + 1) % 4);
				if (!CliffordCheck(cached.DistributionIsPrepared() && CheckAllStatesProbability(cachedRef, cached),
					"Diagonal and Pauli gates preserve or patch prepared probabilities")) return false;
				CliffordRegressionSimulator cold(cached);
				cached.SetSeed(781); cold.SetSeed(781);
				if (!CliffordCheck(cached.SampleBasisStates(16) == cold.SampleBasisStates(16),
					"Patched cache matches a fresh distribution with the same seed")) return false;
			}
		std::mt19937 rng(78381);
		for (size_t n = 2; n <= 7; ++n)
		{
			CliffordRegressionSimulator sim(n);
			QC::QubitRegister<> ref(n);
			sim.SetMultithreading(false); ref.SetMultithreading(false);
			for (size_t step = 0; step < 90; ++step)
			{
				const size_t a = rng() % n, b = (a + 1 + rng() % (n - 1)) % n;
				const int gate = rng() % 15;
				ApplyGate(sim, gate, int(a), int(b)); ref.ApplyGate(*GetGate(gate), a, b);
				const auto probabilities = sim.AllProbabilities();
				for (size_t state = 0; state < probabilities.size(); ++state)
					if (!CliffordCheck(approxEqual(probabilities[state], ref.getBasisStateProbability(state)), "Affine distribution vs state vector")) return false;
				const CliffordRegressionSimulator before(sim);
				sim.SetSeed(77);
				const auto shots = sim.SampleBasisStates(32);
				sim.SetSeed(77);
				if (!CliffordCheck(shots == sim.SampleBasisStates(32) && sim.SameTableau(before), "Terminal sampling preserves state and seeded repeatability")) return false;
				for (const auto& shot : shots)
				{
					size_t state = 0;
					for (size_t q = 0; q < n; ++q) if (shot[q]) state |= size_t(1) << q;
					if (!CliffordCheck(probabilities[state] > 0.0, "Samples stay in affine support")) return false;
				}
			}
		}
		CliffordRegressionSimulator bell(65);
		bell.ApplyH(0);
		for (size_t q = 1; q < 65; ++q) bell.ApplyCX(q, 0);
		bell.SetSeed(71);
		size_t ones = 0;
		for (const auto& shot : bell.SampleBasisStates(2048))
		{
			ones += shot[0];
			for (size_t q = 1; q < 65; ++q)
				if (!CliffordCheck(shot[q] == shot[0], "Sampling preserves wide GHZ correlations")) return false;
		}
		return CliffordCheck(ones > 850 && ones < 1200, "Terminal sampler frequencies");
	}

	bool CliffordBoundaryRegression()
	{
		CliffordRegressionSimulator tiny(2);
		tiny.ApplyH(0);
		if (!CliffordCheck(tiny.ContainsBasisState(1) && tiny.Log2BasisStateProbability(1) == -1.0 &&
			!tiny.ContainsBasisState(2) && tiny.Log2BasisStateProbability(2) == -INFINITY &&
			!tiny.ContainsBasisState(4) && tiny.Log2BasisStateProbability(4) == -INFINITY,
			"Support and logarithmic probability bounds")) return false;
		CliffordRegressionSimulator underflow(1076);
		for (size_t q = 0; q < 1074; ++q) underflow.ApplyH(q);
		if (!CliffordCheck(underflow.getBasisStateProbability(0) == std::ldexp(1.0, -1074) &&
			underflow.Log2BasisStateProbability(0) == -1074.0, "Smallest double basis probability")) return false;
		underflow.ApplyH(1074);
		std::vector<bool> wideBits(1076);
		if (!CliffordCheck(underflow.getBasisStateProbability(0) == 0.0 && underflow.ContainsBasisState(0) &&
			underflow.ContainsBasisState(wideBits) && underflow.Log2BasisStateProbability(wideBits) == -1075.0,
			"Supported outcomes remain distinguishable after underflow")) return false;
		wideBits[1075] = true;
		if (!CliffordCheck(!underflow.ContainsBasisState(wideBits) &&
			underflow.Log2BasisStateProbability(wideBits) == -INFINITY &&
			CliffordRejects<std::invalid_argument>([&] { tiny.ContainsBasisState(wideBits); }) &&
			CliffordRejects<std::invalid_argument>([&] { tiny.Log2BasisStateProbability(wideBits); }),
			"Unsupported outcomes and malformed bit vectors")) return false;
		CliffordRegressionSimulator sim(2), zero(2);
		if (!CliffordCheck(
			CliffordRejects<std::out_of_range>([&] { sim.ApplyH(2); }) &&
			CliffordRejects<std::out_of_range>([&] { sim.ApplySwap(2, 2); }) &&
			CliffordRejects<std::out_of_range>([&] { sim.MeasureQubit(2); }) &&
			CliffordRejects<std::out_of_range>([&] { sim.GetQubitProbability(2); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.ApplyCX(0, 0); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.ApplyCY(0, 0); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.ApplyCZ(0, 0); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.ApplyISwap(0, 0); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.ApplyISwapDag(0, 0); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.ExpectationValue("IIZ"); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.getBasisStateProbability(std::vector<bool>{false}); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.getBasisStateProbability(std::vector<bool>(3)); }) &&
			CliffordRejects<std::runtime_error>([&] { sim.ExpectationValue("Q"); }) && sim.SameTableau(zero),
			"Invalid arguments must be rejected before mutation")) return false;
		if (!CliffordCheck(sim.getBasisStateProbability(4) == 0.0 && sim.getBasisStateProbability(0) == 1.0 &&
			sim.ExpectationValue("Z") == 1.0, "Basis index bounds and short Pauli strings")) return false;
		CliffordRegressionSimulator empty(0);
		if (!CliffordCheck(empty.AllProbabilities() == std::vector<double>{1.0} &&
			empty.getBasisStateProbability(1) == 0.0 && empty.ExpectationValue("") == 1.0,
			"Zero-qubit probabilities")) return false;
		if (sizeof(size_t) == 4)
		{
			CliffordRegressionSimulator wide(32);
			if (!CliffordCheck(CliffordRejects<std::length_error>([&] { wide.AllProbabilities(); }),
				"Reject overflowing probability count")) return false;
		}
		if (!CliffordCheck(CliffordRejects<std::length_error>([] {
			QC::Clifford::StabilizerSimulator tooLarge(std::numeric_limits<size_t>::max());
		}), "Tableau allocation overflow guard")) return false;
		CliffordRegressionSimulator seeded(65);
		seeded.SetSeed(984731);
		seeded.ApplyH(0); seeded.SaveState();
		const auto engine = seeded.RandomEngine();
		const CliffordRegressionSimulator snapshot(seeded);
		CliffordRegressionSimulator moved(std::move(seeded));
		if (!CliffordCheck(moved.RandomEngine() == engine && moved.SameTableau(snapshot), "Move constructor preserves RNG and state")) return false;
		CliffordRegressionSimulator assigned(2);
		assigned.getBasisStateProbability(0); // Populate a cache of a different size before the move.
		assigned = std::move(moved);
		if (!CliffordCheck(assigned.RandomEngine() == engine && assigned.SameTableau(snapshot), "Move assignment preserves RNG and state")) return false;
		assigned.Reset(); assigned.RestoreState();
		if (!CliffordCheck(assigned.SameTableau(snapshot), "Move preserves saved state")) return false;
		QC::PauliStringXZWithSign row(65); row.X[64] = true; row.PhaseSign = true;
		QC::PauliStringXZWithSign movedRow(std::move(row));
		return CliffordCheck(movedRow.X[64] && movedRow.PhaseSign && row.X.empty() && row.Z.empty(), "Pauli move transfers buffers");
	}

	bool CliffordCountsRegression()
	{
		std::mt19937 actions(77521);
		for (size_t n : {2, 7, 65, 129})
		{
			CliffordRegressionSimulator sim(n);
			sim.SetMultithreading(false);
			sim.SaveState(); // Sampling must not replace this zero-state snapshot.
			for (size_t step = 0; step < 90; ++step)
			{
				const size_t a = actions() % n, b = (a + 1 + actions() % (n - 1)) % n;
				ApplyGate(sim, int(actions() % 15), int(a), int(b));
				if (step % 9 != 0) continue;
				for (size_t width : {size_t(1), size_t(3), n + 2})
				{
					std::vector<size_t> qubits(width);
					for (size_t q = 0; q < width; ++q) qubits[q] = (n - 1 - q % n);
					const CliffordRegressionSimulator before(sim);
					CliffordRegressionSimulator ref(sim);
					ref.SaveState(); ref.SetSeed(12391); sim.SetSeed(12391);
					std::unordered_map<std::vector<bool>, size_t> expected;
					std::vector<bool> firstExpected;
					for (size_t shot = 0; shot < 31; ++shot)
					{
						std::vector<bool> bits(width);
						for (size_t q = 0; q < width; ++q) bits[q] = ref.MeasureQubit(qubits[q]);
						if (shot == 0) firstExpected = bits;
						++expected[bits]; ref.RestoreState();
					}
					const auto actual = sim.SampleCountsMany(qubits, 31);
					if (!CliffordCheck(actual == expected && sim.SameTableau(before) &&
						sim.RandomEngine() == ref.RandomEngine(), "Marginal samples and RNG match sequential measurement")) return false;
					sim.SetSeed(12391);
					if (!CliffordCheck(sim.SampleCountsMany(qubits, 31) == expected, "Warm marginal cache repeats seeded counts")) return false;
					CliffordRegressionSimulator single(sim);
					single.SetSeed(12391);
					const auto one = single.SampleCountsMany(qubits, 1);
					if (!CliffordCheck(one.size() == 1 && one.begin()->first == firstExpected && one.begin()->second == 1,
						"Cold one-shot paths match sequential measurement")) return false;
					if (width <= std::numeric_limits<size_t>::digits)
					{
						std::unordered_map<size_t, size_t> packed;
						for (const auto& item : expected)
						{
							size_t key = 0;
							for (size_t q = 0; q < width; ++q) if (item.first[q]) key |= size_t(1) << q;
							packed[key] = item.second;
						}
						sim.SetSeed(12391);
						if (!CliffordCheck(sim.SampleCounts(qubits, 31) == packed, "Packed counts preserve requested bit order")) return false;
					}
				}
			}
			sim.RestoreState();
			CliffordRegressionSimulator zero(n);
			if (!CliffordCheck(sim.SameTableau(zero), "Counts preserve the caller's saved state")) return false;
		}

		CliffordRegressionSimulator sim(129);
		sim.ApplyH(0); sim.ApplyCX(128, 0); sim.ApplyX(128);
		std::vector<size_t> all(129);
		for (size_t q = 0; q < all.size(); ++q) all[q] = 128 - q;
		for (int gate : {1, 2, 5, 11, 3, 4, 0, 9})
		{
			sim.SampleCountsMany(all, 8);
			ApplyGate(sim, gate, 0, 128);
			CliffordRegressionSimulator cold(sim);
			sim.SetSeed(711); cold.SetSeed(711);
			if (!CliffordCheck(sim.SampleCountsMany(all, 1) == cold.SampleCountsMany(all, 1) &&
				sim.RandomEngine() == cold.RandomEngine(), "Cold single-shot fallback matches prepared marginal sampling")) return false;
			sim.SetSeed(912); cold.SetSeed(912);
			if (!CliffordCheck(sim.SampleCountsMany(all, 17) == cold.SampleCountsMany(all, 17),
				"Marginal cache updates after diagonal, Pauli and non-diagonal gates")) return false;
		}
		// Non-contiguous selection, a high physical index, and repeated outputs.
		const std::vector<size_t> repeated{128, 0, 128, 64};
		sim.SampleCountsMany(repeated, 8); sim.ApplyX(128);
		CliffordRegressionSimulator cold(sim);
		sim.SetSeed(918); cold.SetSeed(918);
		if (!CliffordCheck(sim.SampleCountsMany(repeated, 64) == cold.SampleCountsMany(repeated, 64),
			"Pauli cache patch handles repeated selected qubits")) return false;
		const auto engine = sim.RandomEngine();
		if (!CliffordCheck(sim.SampleCounts({}, 3).empty() && sim.SampleCountsMany({0}, 0).empty() &&
			CliffordRejects<std::out_of_range>([&] { sim.SampleCounts({129}, 2); }) &&
			CliffordRejects<std::out_of_range>([&] { sim.SampleCountsMany({129}, 2); }) &&
			CliffordRejects<std::invalid_argument>([&] { sim.SampleCounts(all, 2); }) &&
			sim.RandomEngine() == engine, "Counts validate input before consuming randomness")) return false;
		CliffordRegressionSimulator bits(129);
		std::vector<size_t> word(std::numeric_limits<size_t>::digits);
		for (size_t q = 0; q < word.size(); ++q) { word[q] = 128 - q; bits.ApplyX(word[q]); }
		const auto packedWord = bits.SampleCounts(word, 5);
		if (!CliffordCheck(packedWord.size() == 1 && packedWord.at(std::numeric_limits<size_t>::max()) == 5,
			"Packed counts include the highest size_t bit")) return false;
		// A nontrivial zero-preserving Clifford changes the logical encoding of
		// a GHZ state. This exercises Y factors and signs in the rank-one shortcut.
		for (size_t n : {7, 65, 129}) for (int rep = 0; rep < 12; ++rep)
		{
			CliffordRegressionSimulator encoded(n);
			encoded.SetMultithreading(false);
			for (size_t k = 0; k < 16 * n; ++k)
			{
				const size_t a = actions() % n, b = (a + 1 + actions() % (n - 1)) % n;
				switch (actions() % 4)
				{
				case 0: encoded.ApplyS(a); break;
				case 1: encoded.ApplyCZ(a, b); break;
				case 2: encoded.ApplyCX(a, b); break;
				default: encoded.ApplyZ(a); break;
				}
			}
			encoded.ApplyH(0);
			for (size_t q = 1; q < n; ++q) encoded.ApplyCX(q, 0);
			for (size_t q = 0; q < n; ++q) if (actions() & 1) encoded.ApplyX(q);
			CliffordRegressionSimulator reference(encoded);
			encoded.SetSeed(730 + rep); reference.SetSeed(730 + rep);
			std::vector<size_t> selected(n + 2);
			std::vector<bool> expected(n + 2);
			for (size_t q = 0; q < selected.size(); ++q)
			{
				selected[q] = n - 1 - q % n;
				expected[q] = reference.MeasureQubit(selected[q]);
			}
			const auto actual = encoded.SampleCountsMany(selected, 1);
			if (!CliffordCheck(actual.size() == 1 && actual.begin()->first == expected &&
				encoded.RandomEngine() == reference.RandomEngine(), "Rank-one sampling phase and RNG")) return false;
		}
		return true;
	}

	void PrepareCliffordRegressionState(CliffordRegressionSimulator& sim, bool cnotOnly = false)
	{
		std::mt19937 rng(73);
		const size_t n = sim.getNrQubits();
		sim.SetMultithreading(false);
		for (size_t k = 0; k < 32 * n; ++k)
		{
			const size_t a = rng() % n;
			size_t b = rng() % n;
			if (a == b) b = (a + 1) % n;
			const unsigned gate = cnotOnly ? 2 : rng() % 3;
			if (gate == 0) sim.ApplyH(a);
			else if (gate == 1) sim.ApplyS(a);
			else if (n > 1) sim.ApplyCX(a, b);
		}
	}

	bool CliffordResetRegression()
	{
		for (size_t n : { 0, 1, 7, 65, 1025 })
		{
			CliffordRegressionSimulator sim(n), zero(n);
			PrepareCliffordRegressionState(sim);
			if (n != 0) { sim.ApplyH(0); sim.ApplyS(0); sim.ApplyX(0); }
			const CliffordRegressionSimulator before(sim);
			sim.SaveState();
			for (int repeat = 0; repeat < 2; ++repeat)
			{
				sim.Reset();
				if (!CliffordCheck(sim.SameTableau(zero), "Reset must clear every X/Z bit and sign")) return false;
				for (size_t q = 0; q < n; ++q)
					if (!CliffordCheck(sim.GetQubitProbability(q) == 0.0, "Reset must produce |0>")) return false;
			}
			sim.RestoreState();
			if (!CliffordCheck(sim.SameTableau(before), "Reset must preserve the saved state")) return false;
		}
		return true;
	}

	bool CliffordSwapRegression()
	{
		// All two-qubit Paulis, including Y and both signs, must just permute.
		const char labels[] = { 'I', 'Z', 'X', 'Y' };
		for (unsigned pauli = 0; pauli < 16; ++pauli)
			for (bool sign : { false, true })
			{
				QC::PauliStringXZWithSign p(2);
				for (size_t q = 0; q < 2; ++q)
				{
					p.X[q] = (pauli >> (2 * q + 1)) & 1;
					p.Z[q] = (pauli >> (2 * q)) & 1;
				}
				p.PhaseSign = sign;
				QC::PauliStringXZWithCoefficient weighted(2);
				weighted.X = p.X; weighted.Z = p.Z; weighted.Coefficient = -0.375;
				p.ApplySwap(0, 1);
				weighted.ApplySwap(0, 1);
				std::string expected;
				expected += labels[(pauli >> 2) & 3];
				expected += labels[pauli & 3];
				if (!CliffordCheck(p.PauliStringXZ::ToString() == expected && p.PhaseSign == sign &&
					weighted.PauliStringXZ::ToString() == expected && weighted.Coefficient == -0.375,
					"SWAP must preserve signs and coefficients")) return false;
			}

		for (size_t n : { 2, 65, 1025 })
		{
			CliffordRegressionSimulator initial(n);
			PrepareCliffordRegressionState(initial);
			for (bool parallel : { false, true })
			{
				CliffordRegressionSimulator actual(initial), expected(initial);
				actual.SetMultithreading(parallel);
				for (size_t a : { size_t(0), n / 2, n - 1 })
				{
					const size_t b = (a + 1) % n;
					actual.ApplySwap(a, b);
					expected.ApplyCX(a, b); expected.ApplyCX(b, a); expected.ApplyCX(a, b);
					if (!CliffordCheck(actual.SameTableau(expected), "SWAP must equal three CNOTs")) return false;
					actual.ApplySwap(a, a);
					if (!CliffordCheck(actual.SameTableau(expected), "SWAP on one index must be the identity")) return false;
				}
			}
		}
		return true;
	}

	bool CliffordMeasurementRegression()
	{
		CliffordRegressionSimulator single(1);
		single.ApplyS(0); single.ApplyH(0);
		const bool outcome = single.MeasureQubit(0); // Previously asserted on the paired destabilizer.
		if (!CliffordCheck(single.MeasureQubit(0) == outcome && single.GetQubitProbability(0) == double(outcome),
			"Measurement after S,H must collapse consistently")) return false;

		for (size_t n : { 2, 65, 1025 })
		{
			// Two GHZ groups have four equally likely basis outcomes. S,H ensures
			// the probability path also encounters a pivot with an X destabilizer.
			CliffordRegressionSimulator groups(n);
			groups.SetMultithreading(false);
			groups.ApplyS(0); groups.ApplyH(0); groups.ApplyS(1); groups.ApplyH(1);
			for (size_t q = 2; q < n; ++q) groups.ApplyCX(q, q % 2);
			const CliffordRegressionSimulator before(groups);
			for (bool parallel : { false, true })
			{
				groups.SetMultithreading(parallel);
				for (unsigned pattern = 0; pattern < 4; ++pattern)
				{
					std::vector<bool> bits(n);
					for (size_t q = 0; q < n; ++q) bits[q] = (pattern >> (q % 2)) & 1;
					if (!CliffordCheck(groups.getBasisStateProbability(bits) == 0.25 && groups.SameTableau(before),
						"Basis probabilities must preserve the tableau and support")) return false;
					if (n > 2)
					{
						bits[n - 1] = !bits[n - 1];
						if (!CliffordCheck(groups.getBasisStateProbability(bits) == 0.0 && groups.SameTableau(before),
							"Impossible basis states must have zero probability")) return false;
					}
				}
			}
		}

		std::mt19937 rng(75);
		for (size_t n : { 1024, 1025, 1031, 2049 })
			for (int repeat = 0; repeat < 100; ++repeat)
			{
				QC::PauliStringXZWithSign p(n), identity(n);
				for (size_t q = 0; q < n; ++q) { p.X[q] = rng() & 1; p.Z[q] = rng() & 1; }
				p.PhaseSign = (rng() & 1) != 0;
				const auto source = p;
				p.Multiply(source, true);
				if (!CliffordCheck(p == identity && !p.PhaseSign, "P * P must be identity at word boundaries")) return false;
				p = source;
				p.Multiply(p, true);
				if (!CliffordCheck(p == identity && !p.PhaseSign, "Aliased P * P must be identity")) return false;
			}

		for (size_t n : { 63, 64, 65, 511, 512, 513, 1023, 1024, 1025, 1031, 2048, 2049, 4097 })
		{
			CliffordRegressionSimulator initial(n);
			PrepareCliffordRegressionState(initial);
			for (unsigned seed : { 74, 119 })
			{
				CliffordRegressionSimulator serial(initial), parallel(initial);
				serial.SetSeed(seed); parallel.SetSeed(seed);
				parallel.SetMultithreading(true);
				for (size_t q : { size_t(0), n / 2, n - 1, size_t(1) })
				{
					const bool result = serial.MeasureQubit(q);
					if (!CliffordCheck(parallel.MeasureQubit(q) == result && serial.SameTableau(parallel),
						"Parallel measurement must match the complete serial tableau")) return false;
					if (!CliffordCheck(parallel.MeasureQubit(q) == result && serial.SameTableau(parallel),
						"Repeated measurement must be deterministic")) return false;
					serial.ApplyH((q + 1) % n); parallel.ApplyH((q + 1) % n);
				}
			}
		}

		CliffordRegressionSimulator deterministic(1025);
		PrepareCliffordRegressionState(deterministic, true);
		deterministic.SetMultithreading(true);
		for (size_t q : { 0, 512, 1024 })
			if (!CliffordCheck(deterministic.GetQubitProbability(q) == 0.0 && !deterministic.MeasureQubit(q),
				"Dense deterministic measurements must preserve |0>")) return false;
		return true;
	}
}

namespace {
bool CliffordMixedRegression() {
  size_t circuits=0, measurements=0, checks=0, expectations=0;
  std::mt19937 actions(998877);
  for (size_t n=2; n<=8; ++n) for (size_t rep=0; rep<100; ++rep) {
    CliffordRegressionSimulator s(n); s.SetSeed(987654321+n*100+rep); s.SetMultithreading(false);
    QC::QubitRegister<> ref(n); ref.SetMultithreading(false);
    s.SaveState(); ref.SaveState(); bool saved=true;
    for (size_t step=0; step<160; ++step) {
      const auto op=actions()%100;
      const size_t q=actions()%n, other=(q+1+actions()%(n-1))%n;
      if (op<72) {
        const int code=actions()%15;
        ApplyGate(s,code,int(q),int(other));
        ref.ApplyGate(*GetGate(code),q,other);
      } else if (op<87) {
        const bool value=s.MeasureQubit(q);
        auto v=ref.getRegisterStorage();
        for (size_t b=0; b<size_t(v.size()); ++b) if (bool((b>>q)&1)!=value) v[b]=0.;
        if (v.squaredNorm()<1e-10) { std::cout<<"impossible measurement\n"; return false; }
        ref.setRegisterStorage(v); ++measurements;
        if (s.MeasureQubit(q)!=value) return false;
      } else if (op<90) { s.Reset(); ref.Reset(); }
      else if (op<93) { s.SaveState(); ref.SaveState(); saved=true; }
      else if (op<96 && saved) { s.RestoreState(); ref.RestoreState(); }
      else if (saved) { s.RestoreSavedStateDestructive(); ref.RestoreStateDestructive(); saved=false; }

      for (int sample=0; sample<2; ++sample) {
        std::vector<QC::Gates::AppliedGate<>> p; std::string pauli;
        ConstructPauliString(n,pauli,p);
        if (!approxEqual(s.ExpectationValue(pauli),ref.ExpectationValue(p),1e-8)) {
          std::cout<<"expectation mismatch n="<<n<<" rep="<<rep<<" step="<<step<<" op="<<op<<'\n'; return false;
        }
        ++expectations;
      }
      if (step%5==0) {
        const CliffordRegressionSimulator original(s);
        if (!CheckProbability(ref,s) || !CheckAllStatesProbability(ref,s)) {
          std::cout<<"probability mismatch n="<<n<<" rep="<<rep<<" step="<<step<<'\n'; return false;
        }
        if (!s.SameTableau(original)) { std::cout<<"query changed tableau\n"; return false; }
        auto cloned=s.Clone();
        for (size_t b=0; b<(size_t(1)<<n); ++b)
          if (!approxEqual(cloned->getBasisStateProbability(b),ref.getBasisStateProbability(b),1e-8)) return false;
        ++checks;
      }
    }
    ++circuits;
  }
  std::cout<<"PASS "<<circuits<<" mixed circuits, 160 operations each, "<<measurements
           <<" forced-reference measurements, "<<expectations<<" Pauli expectations, "<<checks
           <<" full probability/marginal/clone/tableau checks, with reset/save/restore\n";
  return true;
}
}

bool CliffordRegressionTests()
{
	std::cout << "\nClifford gates, inverse tableau, probability, sampling and boundary regressions" << std::endl;
	if (!CliffordMixedRegression() || !CliffordPackedRegression() || !CliffordDistributionRegression() || !CliffordBoundaryRegression() || !CliffordCountsRegression() ||
		!CliffordResetRegression() || !CliffordSwapRegression() || !CliffordMeasurementRegression()) return false;
	std::cout << "Success" << std::endl;
	return true;
}

bool CliffordSimulatorTests()
{
	if (!CliffordRegressionTests()) return false;

	const size_t nrTests = 10;
	const size_t nrShots = 100000;
	const double errorThreshold = 0.01;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(5, 20);

	std::cout << "\nClifford gates simulator" << std::endl;

	for (size_t nrQubits = 4; nrQubits < 12; ++nrQubits)
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

			QC::QubitRegister qubitRegister(nrQubits);
			QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

			for (int j = 0; j < static_cast<int>(gates.size()); ++j)
			{
				ApplyGate(cliffordSim, gates[j], qubits1[j], qubits2[j]);
				const auto gateptr = GetGate(gates[j]);
				qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
			}

			if (!CheckProbability(qubitRegister, cliffordSim))
				return false;

			// another way of testing is now available
			if (!CheckAllStatesProbability(qubitRegister, cliffordSim))
				return false;

			// applying the measurement again on the same qubit should give the same result
			if (!CheckMeasurements(cliffordSim))
				return false;

			ExecuteCircuit(nrShots, nrQubits, gates, qubits1, qubits2, results1, results2);

			// check to see if the results are close enough
			for (const auto val : results1)
			{
				if (results2.find(val.first) == results2.end()) continue;

				if (std::abs(static_cast<double>(val.second) - results2[val.first]) / nrShots > errorThreshold)
				{
					std::cout << "\nFailed" << std::endl;
					std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
					std::cout << "Result 1: " << static_cast<double>(val.second) / nrShots << ", Result 2: " << static_cast<double>(results2[val.first]) / nrShots << std::endl;
					return false;
				}
			}

			for (const auto val : results2)
			{
				if (results1.find(val.first) == results1.end()) continue;

				if (std::abs(static_cast<double>(val.second) - results1[val.first]) / nrShots > errorThreshold)
				{
					std::cout << "\nFailed" << std::endl;
					std::cout << "Might fail due of the randomness of the measurements\n" << std::endl;
					std::cout << "Result 1: " << static_cast<double>(results1[val.first]) / nrShots << ", Result 2: " << static_cast<double>(val.second) / nrShots << std::endl;
					return false;
				}
			}

			std::cout << '.';
		}
	}

	std::cout << std::endl;

	std::cout << "Success" << std::endl;

	return true;
}

void ConstructPauliString(size_t nrQubits, std::string& pauliStr, std::vector<QC::Gates::AppliedGate<>>& expGates)
{
	static const QC::Gates::PauliXGate xgate;
	static const QC::Gates::PauliYGate ygate;
	static const QC::Gates::PauliZGate zgate;

	std::uniform_int_distribution pauliDistr(0, 3);

	for (int j = 0; j < static_cast<int>(nrQubits); ++j)
	{
		const int p = pauliDistr(gen);
		switch (p)
		{
		case 0:
			pauliStr += 'I';
			break;
		case 1:
			pauliStr += 'X';
			expGates.emplace_back(xgate.getRawOperatorMatrix(), j);
			break;
		case 2:
			pauliStr += 'Y';
			expGates.emplace_back(ygate.getRawOperatorMatrix(), j);
			break;
		case 3:
			pauliStr += 'Z';
			expGates.emplace_back(zgate.getRawOperatorMatrix(), j);
			break;
		}
	}
}


bool CliffordExpectationValuesTests()
{
	const size_t nrTests = 100;

	std::uniform_int_distribution gateDistr(0, 14);
	std::uniform_int_distribution nrGatesDistr(5, 20);

	std::cout << "\nClifford expectation values" << std::endl;

	for (size_t nrQubits = 2; nrQubits < 20; ++nrQubits)
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

			QC::QubitRegister qubitRegister(nrQubits);
			QC::Clifford::StabilizerSimulator cliffordSim(nrQubits);

			for (int j = 0; j < static_cast<int>(gates.size()); ++j)
			{
				ApplyGate(cliffordSim, gates[j], qubits1[j], qubits2[j]);
				const auto gateptr = GetGate(gates[j]);
				qubitRegister.ApplyGate(*gateptr, qubits1[j], qubits2[j]);
			}

			std::vector<QC::Gates::AppliedGate<>> expGates;
			expGates.reserve(nrQubits);
			std::string pauliStr;

			ConstructPauliString(nrQubits, pauliStr, expGates);

			const auto exp1 = qubitRegister.ExpectationValue(expGates);
			const auto exp2 = cliffordSim.ExpectationValue(pauliStr);
			if (!approxEqual(exp1, exp2, 1E-7))
			{
				std::cout << std::endl << "Expectation values are not equal for stabilizer and statevector simulator for " << nrQubits << " qubits, values: " << exp2 << ", " << exp1 << std::endl;

				std::cout << "Pauli string: " << pauliStr << std::endl;

				std::cout << "Circuit:" << std::endl;
				for (int j = 0; j < static_cast<int>(gates.size()); ++j)
					PrintGate(gates[j], qubits1[j], qubits2[j]);

				return false;
			}
		}
		std::cout << '.';
	}

	std::cout << std::endl;

	std::cout << "Success" << std::endl;

	return true;
}
