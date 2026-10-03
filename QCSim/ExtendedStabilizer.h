#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CliffordProbability.h"
#include "Frame.h"
#include "LocalPauliSum.h"

namespace QC {

	enum class ExtendedStabilizerApproximationMode : unsigned char {
		Exact = 0,
		Approximate
	};

	// Approximate pruning first removes components satisfying
	// |amplitude| / ||state|| <= amplitudeTolerance, then keeps at most the
	// largest maxComponents amplitudes. Equal magnitudes retain their earlier
	// component index, at least one component always survives, and the result is
	// renormalized. Switching back to Exact only disables future pruning: it
	// cannot restore discarded components or clear the historical error bound.
	struct ExtendedStabilizerApproximationPolicy {
		ExtendedStabilizerApproximationMode mode =
			ExtendedStabilizerApproximationMode::Exact;
		double amplitudeTolerance = 0.0;
		// Zero means unlimited. An approximate cap must otherwise be at least one.
		size_t maxComponents = 0;

		static ExtendedStabilizerApproximationPolicy Exact() noexcept
		{
			return {};
		}

		static ExtendedStabilizerApproximationPolicy Approximate(
			double tolerance, size_t componentLimit = 0) noexcept
		{
			return { ExtendedStabilizerApproximationMode::Approximate,
				tolerance, componentLimit };
		}
	};

	struct ExtendedStabilizerApproximationStatistics {
		// Sum of the locally discarded squared-norm fractions. It can exceed one
		// across several pruning events and is therefore called a weight, not a
		// probability.
		double cumulativeDiscardedWeight = 0.0;
		// Sum of sqrt(local discarded fractions), capped at one. This bounds the
		// trace distance under unitary evolution. If a later measurement conditions
		// an already-approximated state, it is conservatively promoted to one because
		// postselection can amplify the previous error without a useful finite factor.
		double traceDistanceErrorBound = 0.0;
		size_t discardedComponents = 0;
		size_t pruningEvents = 0;
	};

	// A stabilizer frame represents the state in an orthonormal basis obtained
	// from the computational basis by a Clifford circuit. Clifford gates rotate
	// that basis, while non-Clifford rotations update the coefficients inside it.
	// Multiple ExtendedFrame objects are retained for a future multiframe implementation;
	// the coherent operations below intentionally implement one frame.
	// Simulator instances are not thread-safe, including concurrent const queries,
	// because lookup tables and packed Pauli workspaces are reused internally.
	class ExtendedStabilizer {
		struct CloneTag {};

	public:
		ExtendedStabilizer() = delete;
		ExtendedStabilizer(const ExtendedStabilizer&) = delete;
		ExtendedStabilizer& operator=(const ExtendedStabilizer&) = delete;

		explicit ExtendedStabilizer(size_t nrQubits)
			: ExtendedStabilizer(nrQubits,
				ExtendedStabilizerApproximationPolicy::Exact())
		{
		}

		ExtendedStabilizer(size_t nrQubits,
			const ExtendedStabilizerApproximationPolicy& policy)
			: approximationPolicy(policy), gen(std::random_device{}()),
			dist(0.0, 1.0)
		{
			ValidateApproximationPolicy(approximationPolicy);
			frames.emplace_back(nrQubits);
		}

		void SetSeed(uint64_t theSeed)
		{
			std::seed_seq seed{ uint32_t(theSeed & 0xffffffff), uint32_t(theSeed >> 32) };
			gen.seed(seed);
		}

		size_t GetNrQubits() const
		{
			if (frames.empty()) return 0;
			return frames.front().GetNrQubits();
		}

		void Reset(size_t nrQubits)
		{
			frames.clear();
			savedFrames.clear();
			approximationStatistics = {};
			savedApproximationStatistics = {};
			savedApproximationPolicy = approximationPolicy;
			frames.emplace_back(nrQubits);
		}

		void SetMultithreading(bool enable = true)
		{
			enableMultithreading = enable;
		}

		bool GetMultithreading() const
		{
			return enableMultithreading;
		}

		void setToBasisState(size_t State)
		{
			if (!IsRepresentable(State))
				throw std::invalid_argument("Basis state is outside the register");
			setToBasisState(BasisStateBits(State));
		}

		void setToBasisState(const std::vector<bool>& state)
		{
			const size_t nrQubits = GetNrQubits();
			if (state.size() > nrQubits)
				throw std::invalid_argument("Basis state bitvector exceeds qubit count");

			frames.clear();
			savedFrames.clear();
			approximationStatistics = {};
			savedApproximationStatistics = {};

			ExtendedFrame frame(nrQubits);
			frame.amplitudes = { {1.0, 0.0} };
			frame.signs.clear();

			const size_t nrWords = frame.signs.GetNrWords();
			std::vector<ExtendedFrame::Word> words(nrWords, 0);
			for (size_t i = 0; i < state.size(); ++i)
				if (state[i])
					SetPackedBit(words.data(), i);

			frame.signs.Append(words.data());
			frames.push_back(std::move(frame));
		}

		// Outcomes with bits beyond the register are impossible. Qubits above
		// the width of size_t are zero in this overload.
		double getBasisStateProbability(size_t State) const
		{
			if (!IsRepresentable(State)) return 0.0;
			return getBasisStateProbability(BasisStateBits(State));
		}

		double getBasisStateProbability(const std::vector<bool>& state) const
		{
			const size_t nrQubits = GetNrQubits();
			if (state.size() != nrQubits)
				throw std::invalid_argument("State size does not match qubit count");

			if (IsStabilizerState())
			{
				// One basis state U|b> is uniform on an affine support.
				const auto& frame = frames.front();
				Clifford::detail::BasisDistribution distribution;
				distribution.BuildFromInverse(frame.cliffordBasis.ZImages(),
					frame.signs.LabelWords(0));
				return distribution.Probability(state);
			}

			auto cloneSim = Clone();
			double totalProb = 1.0;

			for (size_t q = 0; q < nrQubits; ++q)
			{
				const bool targetOutcome = state[q];
				const double p1 = cloneSim->GetQubitProbability(q);
				const double pTarget = targetOutcome ? p1 : (1.0 - p1);

				if (pTarget <= 1E-15 || !std::isfinite(pTarget))
					return 0.0;

				totalProb *= pTarget;
				const bool outcome = cloneSim->MeasureConditioned(q, targetOutcome);
				if (outcome != targetOutcome)
					return 0.0;
			}

			return ClampProbability(totalProb);
		}

		// Computational-basis samples of the selected qubits: bit i of each key
		// is qubits[i], and a repeated index repeats the same value. The live
		// and saved states are preserved; the RNG advances exactly as for the
		// same shots of sequential Measure calls from a restored state.
		std::unordered_map<size_t, size_t> SampleCounts(
			const std::vector<size_t>& qubits, size_t shots)
		{
			if (qubits.size() > std::numeric_limits<size_t>::digits)
				throw std::invalid_argument(
					"Use SampleCountsMany for outcomes wider than size_t");
			for (size_t qubit : qubits) ValidateQubit(qubit);
			std::unordered_map<size_t, size_t> counts;
			if (shots == 0 || qubits.empty()) return counts;
			ForEachSample(qubits, shots, [&](const std::vector<ExtendedFrame::Word>& bits) {
				++counts[static_cast<size_t>(bits[0])];
			});
			return counts;
		}

		std::unordered_map<std::vector<bool>, size_t> SampleCountsMany(
			const std::vector<size_t>& qubits, size_t shots)
		{
			for (size_t qubit : qubits) ValidateQubit(qubit);
			std::unordered_map<std::vector<bool>, size_t> counts;
			if (shots == 0 || qubits.empty()) return counts;
			std::vector<bool> result(qubits.size());
			ForEachSample(qubits, shots, [&](const std::vector<ExtendedFrame::Word>& bits) {
				for (size_t bit = 0; bit < result.size(); ++bit)
					result[bit] = (bits[bit / 64] >> (bit % 64)) & 1U;
				++counts[result];
			});
			return counts;
		}

		// Primarily useful for reproducible measurement runs and regression tests.
		void SetRandomSeed(std::mt19937::result_type seed)
		{
			gen.seed(seed);
			dist.reset();
		}

		// Clifford gates rotate every frame's basis with the direct inverse-map
		// updates shared with the stabilizer simulator; amplitudes are unchanged.
		void ApplyH(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyH(qubit);
		}

		void ApplyS(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyS(qubit);
		}

		void ApplySdg(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplySdg(qubit);
		}

		void ApplyX(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyX(qubit);
		}

		void ApplyY(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyY(qubit);
		}

		void ApplyZ(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyZ(qubit);
		}

		void ApplySx(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplySx(qubit);
		}

		void ApplySxDag(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplySxDag(qubit);
		}

		void ApplyK(size_t qubit)
		{
			ValidateQubit(qubit);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyK(qubit);
		}

		void ApplyCX(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyCX(target, control);
		}

		void ApplyCY(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyCY(target, control);
		}

		void ApplyCZ(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyCZ(target, control);
		}

		void ApplySwap(size_t qubit1, size_t qubit2)
		{
			ValidateTwoQubits(qubit1, qubit2);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplySwap(qubit1, qubit2);
		}

		void ApplyISwap(size_t qubit1, size_t qubit2)
		{
			ValidateTwoQubits(qubit1, qubit2);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyISwap(qubit1, qubit2);
		}

		void ApplyISwapDag(size_t qubit1, size_t qubit2)
		{
			ValidateTwoQubits(qubit1, qubit2);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyISwapDag(qubit1, qubit2);
		}

		void ApplyRx(size_t qubit, double angle)
		{
			ValidateQubit(qubit);
			ApplyAxisRotation(qubit, angle, 'X');
		}

		void ApplyRy(size_t qubit, double angle)
		{
			ValidateQubit(qubit);
			ApplyAxisRotation(qubit, angle, 'Y');
		}

		void ApplyRz(size_t qubit, double angle)
		{
			ValidateQubit(qubit);
			ApplyAxisRotation(qubit, angle, 'Z');
		}

		// The gates below are not decomposed into rotations: each is compiled
		// into one Pauli sum on its qubits followed by Clifford basis updates,
		// and applied in one pass over every frame. At angles where a gate is a
		// Clifford, it only updates the basis. Two-qubit gates take (target,
		// control), as ApplyCX does. Global phase is not tracked.

		// U(theta, phi, lambda) = Rz(phi) Ry(theta) Rz(lambda), up to global phase.
		void ApplyU(size_t qubit, double theta, double phi, double lambda)
		{
			ValidateQubit(qubit);
			RequireFiniteAngles({ theta, phi, lambda });
			auto& gate = BeginGate(1);
			AddU(gate, 0, theta, phi, lambda);
			const size_t qubits[] = { qubit };
			ApplyLocalGate(gate, qubits);
		}

		// Controlled e^(i gamma) U(theta, phi, lambda).
		void ApplyCU(size_t target, size_t control, double theta, double phi,
			double lambda, double gamma = 0.0)
		{
			ValidateTwoQubits(target, control);
			RequireFiniteAngles({ theta, phi, lambda, gamma });
			// Local qubit 0 is the control and 1 the target.
			auto& gate = BeginGate(2);
			gate.Rotate(0, 'Z', gamma)
				.Rotate(1, 'Z', 0.5 * (lambda - phi))
				.Rotate(0, 'Z', 0.5 * (lambda + phi))
				.CX(1, 0);
			AddU(gate, 1, -0.5 * theta, 0.0, -0.5 * (phi + lambda));
			gate.CX(1, 0);
			AddU(gate, 1, 0.5 * theta, phi, 0.0);
			ApplyLocalGate(gate, target, control);
		}

		void ApplyCRx(size_t target, size_t control, double angle)
		{
			ValidateTwoQubits(target, control);
			RequireFiniteAngles({ angle });
			auto& gate = BeginGate(2);
			gate.H(1).CX(1, 0).Rotate(1, 'Z', -0.5 * angle).CX(1, 0)
				.Rotate(1, 'Z', 0.5 * angle).H(1);
			ApplyLocalGate(gate, target, control);
		}

		void ApplyCRy(size_t target, size_t control, double angle)
		{
			ValidateTwoQubits(target, control);
			RequireFiniteAngles({ angle });
			auto& gate = BeginGate(2);
			gate.Rotate(1, 'Y', 0.5 * angle).CX(1, 0)
				.Rotate(1, 'Y', -0.5 * angle).CX(1, 0);
			ApplyLocalGate(gate, target, control);
		}

		void ApplyCRz(size_t target, size_t control, double angle)
		{
			ValidateTwoQubits(target, control);
			RequireFiniteAngles({ angle });
			auto& gate = BeginGate(2);
			gate.Rotate(1, 'Z', 0.5 * angle).CX(1, 0)
				.Rotate(1, 'Z', -0.5 * angle).CX(1, 0);
			ApplyLocalGate(gate, target, control);
		}

		// Controlled Hadamard: Ry(pi/4) CZ Ry(-pi/4) on the target.
		void ApplyCH(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			const double eighthTurn = 0.25 * std::acos(-1.0);
			auto& gate = BeginGate(2);
			gate.Rotate(1, 'Y', -eighthTurn).CZ(1, 0).Rotate(1, 'Y', eighthTurn);
			ApplyLocalGate(gate, target, control);
		}

		// Controlled phase diag(1, 1, 1, e^(i lambda)).
		void ApplyCP(size_t target, size_t control, double lambda)
		{
			ValidateTwoQubits(target, control);
			RequireFiniteAngles({ lambda });
			// Phases at multiples of pi/2 are exact; at odd multiples of pi the
			// gate is CZ.
			std::complex<double> phase;
			long long quarterTurns = 0;
			if (detail::TryGetQuarterTurns(lambda, quarterTurns))
			{
				const unsigned turns = unsigned(((quarterTurns % 4) + 4) % 4);
				if (turns == 0) return;
				if (turns == 2)
				{
					ApplyCZ(target, control);
					return;
				}
				phase = detail::LocalPauliSum::IPower(turns);
			}
			else
				phase = std::polar(1.0, lambda);
			auto& gate = BeginGate(2);
			gate.Multiply(Controlled(0, detail::LocalPauliSum::Projector(1, false)
				+ detail::LocalPauliSum::Projector(1, true) * detail::LocalPauliSum::Identity(phase)));
			ApplyLocalGate(gate, target, control);
		}

		void ApplyCS(size_t target, size_t control)
		{
			ApplyCP(target, control, 0.5 * std::acos(-1.0));
		}

		void ApplyCSdg(size_t target, size_t control)
		{
			ApplyCP(target, control, -0.5 * std::acos(-1.0));
		}

		// Controlled square root of X, ((1 + i) I + (1 - i) X) / 2, and its inverse.
		void ApplyCSx(size_t target, size_t control)
		{
			ApplyControlledSquareRootX(target, control, false);
		}

		void ApplyCSxDag(size_t target, size_t control)
		{
			ApplyControlledSquareRootX(target, control, true);
		}

		void ApplyCCX(size_t target, size_t control1, size_t control2)
		{
			ValidateThreeQubits(target, control1, control2);
			// Local qubits 0 and 1 are the controls, 2 the target.
			auto& gate = BeginGate(3);
			gate.Multiply(Controlled(0, Controlled(1, detail::LocalPauliSum::Pauli(2, 'X'))));
			const size_t qubits[] = { control1, control2, target };
			ApplyLocalGate(gate, qubits);
		}

		void ApplyCSwap(size_t target1, size_t target2, size_t control)
		{
			ValidateThreeQubits(target1, target2, control);
			// Local qubit 0 is the control; SWAP = (II + XX + YY + ZZ) / 2.
			using Sum = detail::LocalPauliSum;
			Sum swap = Sum::Identity(0.5);
			for (const char pauli : { 'X', 'Y', 'Z' })
				swap = swap + Sum::Pauli(1, pauli, 0.5) * Sum::Pauli(2, pauli);
			auto& gate = BeginGate(3);
			gate.Multiply(Controlled(0, swap));
			const size_t qubits[] = { control, target1, target2 };
			ApplyLocalGate(gate, qubits);
		}

		bool Measure(size_t qubit, const bool* forcedOutcome = nullptr)
		{
			ValidateQubit(qubit);
			AccountForMeasurementConditioning();
			auto& frame = frames.front();
			if (frame.GetFrameSize() == 1)
				return MeasureSingleStabilizer(qubit, forcedOutcome);

			const auto observable = frame.cliffordBasis.ImageZ(qubit);
			CompilePauliAction(frame, observable, pauliActionWorkspace);
			const auto& action = pauliActionWorkspace;
			if (!HasPauliFlip(action))
				return MeasureDiagonalFrame(frame, action, forcedOutcome);
			return MeasureOffDiagonalFrame(frame, qubit, action, forcedOutcome);
		}

		bool MeasureConditioned(size_t qubit, bool forcedOutcome)
		{
			return Measure(qubit, &forcedOutcome);
		}

		double GetQubitProbability(size_t qubit) const
		{
			ValidateQubit(qubit);
			const auto& frame = frames.front();
			const auto observable = frame.cliffordBasis.ImageZ(qubit);
			return ClampProbability(0.5 * (1.0 - PauliExpectation(frames.front(), observable)));
		}

		double ExpectationValue(const std::string& pauliString) const
		{
			if (pauliString.size() > GetNrQubits())
				throw std::invalid_argument("Pauli string is longer than the register");

			// The logical images of the one-qubit factors are multiplied directly
			// into packed scratch storage; no physical Pauli string is built.
			// For a single basis state |b>, XOR accumulation decides most cases
			// without phase arithmetic: an off-diagonal image has expectation 0,
			// and a product of diagonal rows (-1)^s Z^z has (-1)^(s + z.b).
			const auto& frame = frames.front();
			const auto& basis = frame.cliffordBasis;
			const bool singleComponent = frame.GetFrameSize() == 1;
			bool isIdentity = true, allDiagonal = true;
			basis.ClearImage();
			for (size_t qubit = 0; qubit < pauliString.size(); ++qubit)
			{
				const char pauli = ToPauli(pauliString[qubit]);
				if (pauli == 'I') continue;
				isIdentity = false;
				if (singleComponent && basis.XorImage(qubit, pauli)) allDiagonal = false;
			}

			if (isIdentity) return 1.0;
			if (singleComponent)
			{
				const auto image = basis.Image();
				if (image.HasX()) return 0.0;
				if (allDiagonal)
				{
					const auto* label = frame.signs.LabelWords(0);
					ExtendedFrame::Word parity = 0;
					for (size_t word = 0; word < image.Words(); ++word)
						parity ^= image.Z.words[word] & label[word];
					const bool negative = bool(image.PhaseSign)
						!= ((Clifford::detail::Popcount(parity) & 1U) != 0);
					return negative ? -1.0 : 1.0;
				}
			}
			basis.ClearImage();
			for (size_t qubit = 0; qubit < pauliString.size(); ++qubit)
			{
				const char pauli = ToPauli(pauliString[qubit]);
				if (pauli != 'I') basis.MultiplyImage(qubit, pauli);
			}
			return PauliExpectation(frame, basis.Image());
		}

		void SaveState()
		{
			savedFrames = frames;
			savedApproximationPolicy = approximationPolicy;
			savedApproximationStatistics = approximationStatistics;
		}

		void RestoreState()
		{
			if (savedFrames.empty()) return;

			frames = savedFrames;
			approximationPolicy = savedApproximationPolicy;
			approximationStatistics = savedApproximationStatistics;
		}

		std::unique_ptr<ExtendedStabilizer> Clone() const
		{
			return std::unique_ptr<ExtendedStabilizer>(
				new ExtendedStabilizer(*this, CloneTag{}));
		}

		const std::vector<ExtendedFrame>& GetFrames() const noexcept
		{
			return frames;
		}

		const ExtendedStabilizerApproximationPolicy&
			GetApproximationPolicy() const noexcept
		{
			return approximationPolicy;
		}

		const ExtendedStabilizerApproximationStatistics&
			GetApproximationStatistics() const noexcept
		{
			return approximationStatistics;
		}

		double GetApproximationErrorBound() const noexcept
		{
			return approximationStatistics.traceDistanceErrorBound;
		}

		void SetApproximationPolicy(
			const ExtendedStabilizerApproximationPolicy& policy)
		{
			ValidateApproximationPolicy(policy);
			approximationPolicy = policy;
			if (policy.mode == ExtendedStabilizerApproximationMode::Approximate)
				for (auto& frame : frames)
					PruneComponents(frame);
		}

	private:
		struct ScaledNormAccumulator {
			void Add(double magnitude) noexcept
			{
				if (magnitude == 0.0) return;
				if (scale < magnitude)
				{
					const double ratio = scale / magnitude;
					sumSquares = 1.0 + sumSquares * ratio * ratio;
					scale = magnitude;
				}
				else
				{
					const double ratio = magnitude / scale;
					sumSquares += ratio * ratio;
				}
			}

			double Value() const noexcept
			{
				return scale == 0.0 ? 0.0 : scale * std::sqrt(sumSquares);
			}

			double scale = 0.0;
			double sumSquares = 1.0;
		};

		ExtendedStabilizer(const ExtendedStabilizer& other, CloneTag)
			: frames(other.frames), savedFrames(other.savedFrames),
			approximationPolicy(other.approximationPolicy),
			savedApproximationPolicy(other.savedApproximationPolicy),
			approximationStatistics(other.approximationStatistics),
			savedApproximationStatistics(other.savedApproximationStatistics),
			gen(other.gen), dist(other.dist),
			enableMultithreading(other.enableMultithreading)
		{
		}

		bool IsRepresentable(size_t state) const noexcept
		{
			const size_t nrQubits = GetNrQubits();
			return nrQubits >= static_cast<size_t>(std::numeric_limits<size_t>::digits)
				|| (state >> nrQubits) == 0;
		}

		std::vector<bool> BasisStateBits(size_t state) const
		{
			const size_t nrQubits = GetNrQubits();
			std::vector<bool> bits(nrQubits, false);
			const size_t width = std::min(nrQubits,
				static_cast<size_t>(std::numeric_limits<size_t>::digits));
			for (size_t qubit = 0; qubit < width; ++qubit)
				bits[qubit] = ((state >> qubit) & 1U) != 0;
			return bits;
		}

		// A single component is an ordinary stabilizer state U|b>.
		bool IsStabilizerState() const noexcept
		{
			return frames.size() == 1 && frames.front().GetFrameSize() == 1;
		}

		template<class Consumer>
		void ForEachSample(const std::vector<size_t>& qubits, size_t shots, Consumer consume)
		{
			std::vector<ExtendedFrame::Word> bits((qubits.size() + 63) / 64);
			if (IsStabilizerState())
			{
				// Sample the affine support directly. It consumes one fair draw per
				// random outcome, in the order of sequential measurements.
				const auto& frame = frames.front();
				Clifford::detail::BasisDistribution distribution;
				distribution.BuildFromInverse(frame.cliffordBasis.ZImages(), qubits,
					frame.signs.LabelWords(0));
				auto coin = [this](std::mt19937_64& engine) { return dist(engine) < 0.5; };
				for (size_t shot = 0; shot < shots; ++shot)
				{
					distribution.SampleInto(bits, gen, coin);
					consume(bits);
				}
				return;
			}

			// Superpositions of basis states measure a private copy per shot.
			auto sampler = Clone();
			sampler->SaveState();
			for (size_t shot = 0; shot < shots; ++shot)
			{
				sampler->RestoreState();
				std::fill(bits.begin(), bits.end(), ExtendedFrame::Word(0));
				for (size_t bit = 0; bit < qubits.size(); ++bit)
					if (sampler->Measure(qubits[bit]))
						bits[bit / 64] |= ExtendedFrame::Word(1) << (bit % 64);
				consume(bits);
			}
			gen = sampler->gen;
			dist = sampler->dist;
		}

		void AccountForMeasurementConditioning() noexcept
		{
			if (approximationStatistics.traceDistanceErrorBound > 0.0)
				approximationStatistics.traceDistanceErrorBound = 1.0;
		}

		static void ValidateApproximationPolicy(
			const ExtendedStabilizerApproximationPolicy& policy)
		{
			if (policy.mode != ExtendedStabilizerApproximationMode::Exact
				&& policy.mode != ExtendedStabilizerApproximationMode::Approximate)
				throw std::invalid_argument("Unknown approximation mode");
			if (!std::isfinite(policy.amplitudeTolerance)
				|| policy.amplitudeTolerance < 0.0)
				throw std::invalid_argument(
					"Component amplitude tolerance must be finite and nonnegative");
			if (policy.mode == ExtendedStabilizerApproximationMode::Exact)
			{
				if (policy.amplitudeTolerance != 0.0
					|| policy.maxComponents != 0)
					throw std::invalid_argument(
						"Exact mode cannot specify pruning parameters");
				return;
			}
			if (policy.amplitudeTolerance == 0.0
				&& policy.maxComponents == 0)
				throw std::invalid_argument(
					"Approximate mode needs a tolerance or component limit");
		}

		static std::complex<double> QuarterTurnGlobalPhase(long long quarterTurns)
		{
			int phase = static_cast<int>(quarterTurns % 8);
			if (phase < 0) phase += 8;
			const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);
			switch (phase)
			{
			case 0: return { 1.0, 0.0 };
			case 1: return { inverseSqrtTwo, -inverseSqrtTwo };
			case 2: return { 0.0, -1.0 };
			case 3: return { -inverseSqrtTwo, -inverseSqrtTwo };
			case 4: return { -1.0, 0.0 };
			case 5: return { -inverseSqrtTwo, inverseSqrtTwo };
			case 6: return { 0.0, 1.0 };
			default: return { inverseSqrtTwo, inverseSqrtTwo };
			}
		}

		void ApplyCliffordRotation(size_t physicalQubit, char axis,
			long long quarterTurns)
		{
			int mapTurns = static_cast<int>(quarterTurns % 4);
			if (mapTurns < 0) mapTurns += 4;
			const auto globalPhase = QuarterTurnGlobalPhase(quarterTurns);
			for (auto& frame : frames)
			{
				if (mapTurns == 1)
					frame.cliffordBasis.ApplyQuarterTurn(physicalQubit, axis);
				else if (mapTurns == 2)
				{
					frame.cliffordBasis.ApplyQuarterTurn(physicalQubit, axis);
					frame.cliffordBasis.ApplyQuarterTurn(physicalQubit, axis);
				}
				else if (mapTurns == 3)
					frame.cliffordBasis.ApplyQuarterTurn(physicalQubit, axis, true);

				for (auto& amplitude : frame.amplitudes)
					amplitude *= globalPhase;
			}
		}

		void ValidateQubit(size_t qubit) const
		{
			if (qubit >= GetNrQubits())
				throw std::out_of_range("Qubit index is outside the register");
		}

		void ValidateTwoQubits(size_t qubit1, size_t qubit2) const
		{
			ValidateQubit(qubit1);
			ValidateQubit(qubit2);
			if (qubit1 == qubit2)
				throw std::invalid_argument("A two-qubit gate needs two distinct qubits");
		}

		static char ToPauli(char pauli)
		{
			switch (pauli)
			{
			case 'I': case 'i': return 'I';
			case 'X': case 'x': return 'X';
			case 'Y': case 'y': return 'Y';
			case 'Z': case 'z': return 'Z';
			default: throw std::runtime_error("Invalid operator in the Pauli string");
			}
		}

		static double ClampProbability(double probability)
		{
			return std::max(0.0, std::min(1.0, probability));
		}

		struct PauliAction {
			void ReferenceMasks(const ExtendedFrame::Word* flip,
				const ExtendedFrame::Word* phase, size_t words) noexcept
			{
				flipMask = flip;
				phaseMask = phase;
				nrWords = words;
			}

			const ExtendedFrame::Word* flipMask = nullptr;
			const ExtendedFrame::Word* phaseMask = nullptr;
			size_t nrWords = 0;
			std::complex<double> basePhase{ 1.0, 0.0 };
		};

		static void SetPackedBit(ExtendedFrame::Word* words,
			size_t bit) noexcept
		{
			words[bit / PackedComponentLabels::BitsPerWord]
				|= ExtendedFrame::Word(1)
				<< (bit % PackedComponentLabels::BitsPerWord);
		}

		static void CompilePauliAction(const ExtendedFrame& frame,
			const CliffordBasisMap::Row& pauli,
			PauliAction& action)
		{
			const size_t nrWords = frame.signs.GetNrWords();
			action.ReferenceMasks(pauli.X.words, pauli.Z.words, nrWords);
			action.basePhase = { 1.0, 0.0 };

			size_t nrY = 0;
			for (size_t word = 0; word < nrWords; ++word)
				nrY += Clifford::detail::Popcount(action.flipMask[word] & action.phaseMask[word]);

			switch (nrY % 4)
			{
			case 0:
				action.basePhase = { 1.0, 0.0 };
				break;
			case 1:
				action.basePhase = { 0.0, 1.0 };
				break;
			case 2:
				action.basePhase = { -1.0, 0.0 };
				break;
			default:
				action.basePhase = { 0.0, -1.0 };
				break;
			}

			if (pauli.PhaseSign)
				action.basePhase = -action.basePhase;
		}

		static std::complex<double> PauliPhase(
			const ExtendedFrame::Word* sourceSigns,
			const PauliAction& action)
		{
			unsigned parity = 0;
			for (size_t word = 0; word < action.nrWords; ++word)
				parity ^= Clifford::detail::Popcount(sourceSigns[word] & action.phaseMask[word]) & 1U;
			const bool negate = parity != 0;
			return negate ? -action.basePhase : action.basePhase;
		}

		static bool HasPauliFlip(const PauliAction& action)
		{
			for (size_t word = 0; word < action.nrWords; ++word)
				if (action.flipMask[word] != 0) return true;
			return false;
		}

		static bool DiagonalPauliOutcome(
			const ExtendedFrame::Word* logicalLabel,
			const PauliAction& action)
		{
			// A diagonal logical Pauli has no Y factors, hence basePhase is
			// exactly +1 or -1.  Its eigenvalue on |b> is
			// (-1)^(baseMinus + z.b).
			bool outcome = action.basePhase.real() < 0.0;
			for (size_t word = 0; word < action.nrWords; ++word)
				if ((Clifford::detail::Popcount(action.phaseMask[word] & logicalLabel[word]) & 1U) != 0)
					outcome = !outcome;
			return outcome;
		}

		bool MeasureDiagonalFrame(ExtendedFrame& frame,
			const PauliAction& action, const bool* forcedOutcome = nullptr)
		{
			double probabilityZero = 0.0;
			double probabilityOne = 0.0;
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
			{
				double& probability = DiagonalPauliOutcome(
					frame.signs.LabelWords(component), action)
					? probabilityOne : probabilityZero;
				probability += std::norm(frame.amplitudes[component]);
			}

			const double totalProbability = probabilityZero + probabilityOne;
			if (totalProbability <= 0.0)
				throw std::runtime_error(
					"Cannot measure a frame with zero total probability");
			const double normalizedProbabilityOne = ClampProbability(
				probabilityOne / totalProbability);
			const bool outcome = forcedOutcome ? *forcedOutcome
				: (normalizedProbabilityOne >= 1.0
					|| (normalizedProbabilityOne > 0.0
						&& dist(gen) < normalizedProbabilityOne));

			const double outcomeProbability = outcome
				? probabilityOne : probabilityZero;
			if (outcomeProbability <= 0.0)
				throw std::runtime_error(
					"Cannot normalize a zero-probability measurement outcome");
			const double inverseNorm = 1.0 / std::sqrt(outcomeProbability);

			size_t write = 0;
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
				if (DiagonalPauliOutcome(frame.signs.LabelWords(component), action) == outcome)
				{
					if (write != component)
					{
						frame.amplitudes[write] = frame.amplitudes[component];
						frame.signs.CopyLabel(write, component);
					}
					frame.amplitudes[write] *= inverseNorm;
					++write;
				}

			frame.amplitudes.resize(write);
			frame.signs.resize(write);
			frame.InvalidateComponentIndex();
			if (write == 0)
				throw std::runtime_error(
					"Measurement removed every frame component");
			return outcome;
		}

		bool MeasureOffDiagonalFrame(ExtendedFrame& frame,
			size_t physicalQubit, const PauliAction& action, const bool* forcedOutcome = nullptr)
		{
			// The action's flip mask is the X part of U^dagger Z U.
			size_t pivot = 0;
			Clifford::detail::FindPivot(frame.cliffordBasis.ImageZ(physicalQubit), pivot);

			frame.EnsureComponentIndex(frame.GetFrameSize());
			const size_t notFound = PackedComponentIndex::NotFound;
			const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);

			// Since x[pivot] is one, every logical label belongs to a unique pair
			// {r, r xor x} with r[pivot] == 0. Cache the pairs in frame-owned
			// storage: repeated measurements allocate nothing at steady state and
			// the probability/collapse passes share one set of hash lookups.
			auto& pairs = frame.measurementPairs;
			pairs.clear();
			pairs.reserve(frame.GetFrameSize());
			for (size_t component = 0;
				component < frame.GetFrameSize(); ++component)
			{
				const auto* componentLabel = frame.signs.LabelWords(component);
				const size_t partner = frame.FindXorComponent(componentLabel,
					action.flipMask);
				size_t root = component;
				size_t target = partner;
				if (frame.signs.Get(component, pivot))
				{
					if (partner != notFound) continue;
					root = notFound;
					target = component;
				}
				pairs.emplace_back(root, target);
			}

			auto pairAmplitudes = [&](size_t root, size_t target)
			{
				const std::complex<double> rootAmplitude = root == notFound
					? std::complex<double>(0.0, 0.0)
					: frame.amplitudes[root];
				const std::complex<double> targetAmplitude = target == notFound
					? std::complex<double>(0.0, 0.0)
					: frame.amplitudes[target];

				// If Q|b> = phase(b)|b xor x>, projection onto outcome m
				// gives the new-basis coefficient
				// (a_r + (-1)^m phase(r xor x) a_(r xor x))/sqrt(2).
				const auto mappedTargetAmplitude =
					target == notFound ? std::complex<double>(0.0, 0.0)
					: PauliPhase(frame.signs.LabelWords(target), action)
						* targetAmplitude;
				return std::pair<std::complex<double>, std::complex<double>>{
					inverseSqrtTwo * (rootAmplitude + mappedTargetAmplitude),
					inverseSqrtTwo * (rootAmplitude - mappedTargetAmplitude) };
			};

			double probabilityZero = 0.0;
			double probabilityOne = 0.0;
			for (const auto& pair : pairs)
			{
				const size_t root = pair.first;
				const size_t target = pair.second;
				const auto amplitudes = pairAmplitudes(root, target);
				const auto& outcomeZeroAmplitude = amplitudes.first;
				const auto& outcomeOneAmplitude = amplitudes.second;
				probabilityZero += std::norm(outcomeZeroAmplitude);
				probabilityOne += std::norm(outcomeOneAmplitude);
			}

			const double totalProbability = probabilityZero + probabilityOne;
			if (totalProbability <= 0.0)
				throw std::runtime_error(
					"Cannot measure a frame with zero total probability");
			const double normalizedProbabilityOne = ClampProbability(
				probabilityOne / totalProbability);
			const bool outcome = forcedOutcome ? *forcedOutcome
				: (normalizedProbabilityOne >= 1.0
					|| (normalizedProbabilityOne > 0.0
						&& dist(gen) < normalizedProbabilityOne));

			const double outcomeProbability = outcome
				? probabilityOne : probabilityZero;
			if (outcomeProbability <= 0.0)
				throw std::runtime_error(
					"Cannot normalize a zero-probability measurement outcome");
			const double inverseOutcomeNorm = 1.0 / std::sqrt(outcomeProbability);
			auto& collapsedAmplitudes = frame.nextAmplitudes;
			auto& collapsedSigns = frame.nextSigns;
			collapsedAmplitudes.clear();
			collapsedSigns.clear();
			collapsedAmplitudes.reserve(frame.GetFrameSize());
			collapsedSigns.reserve(frame.GetFrameSize());
			for (const auto& pair : pairs)
			{
				const size_t root = pair.first;
				const size_t target = pair.second;
				const auto amplitudes = pairAmplitudes(root, target);
				const auto amplitude = (outcome ? amplitudes.second : amplitudes.first)
					* inverseOutcomeNorm;
				if (amplitude.real() == 0.0 && amplitude.imag() == 0.0)
					continue;
				collapsedAmplitudes.push_back(amplitude);
				if (root != notFound)
					collapsedSigns.Append(frame.signs.LabelWords(root));
				else
					collapsedSigns.AppendXor(frame.signs.LabelWords(target),
						action.flipMask);
				collapsedSigns.Set(collapsedSigns.size() - 1, pivot, outcome);
			}

			if (collapsedAmplitudes.empty())
				throw std::runtime_error(
					"Measurement removed every frame component");

			frame.cliffordBasis.RebaseZ(physicalQubit, pivot, enableMultithreading);
			frame.amplitudes.swap(collapsedAmplitudes);
			frame.signs.swap(collapsedSigns);
			frame.InvalidateComponentIndex();
			PruneComponents(frame);
			return outcome;
		}

		void PruneComponents(ExtendedFrame& frame)
		{
			const bool approximate = approximationPolicy.mode
				== ExtendedStabilizerApproximationMode::Approximate;
			const size_t originalSize = frame.GetFrameSize();

			// Exact mode removes only coefficients whose real and imaginary parts are
			// both zero. Squaring a tiny but representable coefficient can underflow,
			// so std::norm must not decide exact liveness.
			if (!approximate)
			{
				size_t write = 0;
				for (size_t component = 0; component < originalSize; ++component)
				{
					const auto amplitude = frame.amplitudes[component];
					if (amplitude.real() == 0.0 && amplitude.imag() == 0.0)
						continue;
					if (write != component)
					{
						frame.amplitudes[write] = amplitude;
						frame.signs.CopyLabel(write, component);
					}
					++write;
				}
				if (write == 0)
					throw std::runtime_error(
						"A frame operation cancelled every component");
				if (write != originalSize)
				{
					frame.amplitudes.resize(write);
					frame.signs.resize(write);
					frame.InvalidateComponentIndex();
				}
				return;
			}

			auto& retained = frame.componentOrderWorkspace;
			retained.clear();
			retained.reserve(originalSize);
			auto& magnitudes = frame.componentMagnitudeWorkspace;
			magnitudes.resize(originalSize);

			ScaledNormAccumulator totalNormAccumulator;
			double largestMagnitude = -1.0;
			size_t largestComponent = 0;
			size_t positiveComponents = 0;
			for (size_t component = 0; component < originalSize; ++component)
			{
				const double magnitude = std::abs(frame.amplitudes[component]);
				magnitudes[component] = magnitude;
				if (!std::isfinite(magnitude))
					throw std::runtime_error("A frame contains a non-finite amplitude");
				totalNormAccumulator.Add(magnitude);
				if (magnitude == 0.0) continue;
				++positiveComponents;
				if (magnitude > largestMagnitude)
				{
					largestMagnitude = magnitude;
					largestComponent = component;
				}
			}

			const double totalNorm = totalNormAccumulator.Value();
			if (totalNorm <= 0.0 || positiveComponents == 0)
				throw std::runtime_error("A frame operation cancelled every component");
			const double amplitudeCutoff =
				approximationPolicy.amplitudeTolerance * totalNorm;
			for (size_t component = 0; component < originalSize; ++component)
			{
				const double magnitude = magnitudes[component];
				if (magnitude > amplitudeCutoff)
					retained.push_back(component);
			}
			// Approximation is never allowed to erase the complete state.
			if (retained.empty()) retained.push_back(largestComponent);

			const size_t componentLimit = approximate
				? approximationPolicy.maxComponents : 0;
			if (componentLimit != 0 && retained.size() > componentLimit)
			{
				auto heavier = [&magnitudes](size_t left, size_t right)
				{
					const double leftMagnitude = magnitudes[left];
					const double rightMagnitude = magnitudes[right];
					return leftMagnitude != rightMagnitude
						? leftMagnitude > rightMagnitude : left < right;
				};
				std::nth_element(retained.begin(),
					retained.begin() + componentLimit, retained.end(), heavier);
				retained.resize(componentLimit);
				std::sort(retained.begin(), retained.end());
			}

			// Sum retained and discarded norms independently. Subtracting two totals
			// loses small discarded branches to catastrophic cancellation and can turn
			// a real approximation error into a reported zero bound.
			ScaledNormAccumulator retainedNormAccumulator;
			ScaledNormAccumulator discardedNormAccumulator;
			size_t retainedPosition = 0;
			for (size_t component = 0; component < originalSize; ++component)
			{
				const double magnitude = magnitudes[component];
				if (retainedPosition < retained.size()
					&& retained[retainedPosition] == component)
				{
					retainedNormAccumulator.Add(magnitude);
					++retainedPosition;
				}
				else if (magnitude != 0.0)
					discardedNormAccumulator.Add(magnitude);
			}
			const double retainedNorm = retainedNormAccumulator.Value();
			const double discardedNorm = discardedNormAccumulator.Value();
			if (retainedNorm <= 0.0)
				throw std::runtime_error("Approximation retained a zero-norm state");

			if (retained.size() != originalSize)
			{
				for (size_t write = 0; write < retained.size(); ++write)
				{
					const size_t source = retained[write];
					if (write == source) continue;
					frame.amplitudes[write] = frame.amplitudes[source];
					frame.signs.CopyLabel(write, source);
				}
				frame.amplitudes.resize(retained.size());
				frame.signs.resize(retained.size());
				frame.InvalidateComponentIndex();
			}

			const size_t approximateDiscardedComponents =
				positiveComponents - retained.size();
			if (approximateDiscardedComponents == 0) return;

			const double inverseNorm = 1.0 / retainedNorm;
			for (auto& amplitude : frame.amplitudes)
				amplitude *= inverseNorm;
			const double localTraceDistance = std::min(1.0,
				std::nextafter(discardedNorm / totalNorm,
					std::numeric_limits<double>::infinity()));
			const double localDiscardedFraction =
				localTraceDistance * localTraceDistance;
			approximationStatistics.cumulativeDiscardedWeight +=
				localDiscardedFraction;
			approximationStatistics.traceDistanceErrorBound = std::min(1.0,
				std::nextafter(
					approximationStatistics.traceDistanceErrorBound
						+ localTraceDistance,
					std::numeric_limits<double>::infinity()));
			approximationStatistics.discardedComponents +=
				approximateDiscardedComponents;
			++approximationStatistics.pruningEvents;
		}

		static void AddU(detail::LocalGate& gate, size_t qubit, double theta,
			double phi, double lambda)
		{
			gate.Rotate(qubit, 'Z', lambda).Rotate(qubit, 'Y', theta).Rotate(qubit, 'Z', phi);
		}

		// |0><0| on the control plus |1><1| on the control times op.
		static detail::LocalPauliSum Controlled(size_t control,
			const detail::LocalPauliSum& op)
		{
			return detail::LocalPauliSum::Projector(control, false)
				+ detail::LocalPauliSum::Projector(control, true) * op;
		}

		void ApplyControlledSquareRootX(size_t target, size_t control, bool inverse)
		{
			ValidateTwoQubits(target, control);
			const std::complex<double> plus(0.5, 0.5), minus(0.5, -0.5);
			auto& gate = BeginGate(2);
			gate.Multiply(Controlled(0, detail::LocalPauliSum::Identity(inverse ? minus : plus)
				+ detail::LocalPauliSum::Pauli(1, 'X', inverse ? plus : minus)));
			ApplyLocalGate(gate, target, control);
		}

		static void RequireFiniteAngles(std::initializer_list<double> angles)
		{
			for (const double angle : angles)
				if (!std::isfinite(angle))
					throw std::invalid_argument("Rotation angle must be finite");
		}

		void ValidateThreeQubits(size_t qubit1, size_t qubit2, size_t qubit3) const
		{
			ValidateTwoQubits(qubit1, qubit2);
			ValidateQubit(qubit3);
			if (qubit3 == qubit1 || qubit3 == qubit2)
				throw std::invalid_argument("A three-qubit gate needs three distinct qubits");
		}

		// The reusable gate being compiled; gates do not nest.
		detail::LocalGate& BeginGate(size_t qubits)
		{
			localGate.Reset(qubits);
			return localGate;
		}

		// Local qubit 0 is the control and 1 the target.
		void ApplyLocalGate(const detail::LocalGate& gate, size_t target, size_t control)
		{
			const size_t qubits[] = { control, target };
			ApplyLocalGate(gate, qubits);
		}

		// qubits maps the gate's local qubits to physical qubits.
		void ApplyLocalGate(const detail::LocalGate& gate, const size_t* qubits)
		{
			const bool hasSum = !gate.Sum().IsScalar();
			for (auto& frame : frames)
			{
				if (hasSum)
					ApplyPauliSum(frame, gate.Sum(), qubits, gate.GetNrQubits());
				gate.ReplayCliffords(frame.cliffordBasis, qubits);
			}
		}

		struct PauliSumWorkspace {
			static constexpr size_t MaxOffsets = detail::LocalPauliSum::Keys;
			static constexpr size_t MaxGenerators = 2 * detail::LocalPauli::MaxQubits;

			struct Term {
				std::complex<double> coefficient;  // includes the image's phase
				size_t offset;                     // coset offset index of its X part
				bool hasZ;
			};

			// Term images then reduced generators, and the coset offsets; both
			// grow to the largest sum seen.
			Clifford::detail::PackedTableau rows;
			Clifford::detail::PackedTableau offsets;
			std::vector<Term> terms;
			std::vector<uint8_t> offsetParity;  // term-major parity(z . offset)
			std::vector<uint8_t> visited;
			std::vector<ExtendedFrame::Word> label;
		};

		static void EnsureRows(Clifford::detail::PackedTableau& rows, size_t count, size_t qubits)
		{
			if (rows.GetNrQubits() == qubits && rows.size() >= count) return;
			rows = Clifford::detail::PackedTableau(
				rows.GetNrQubits() == qubits ? std::max(count, rows.size()) : count, qubits);
		}

		static bool OddParity(const ExtendedFrame::Word* left,
			const ExtendedFrame::Word* right, size_t nrWords) noexcept
		{
			ExtendedFrame::Word parity = 0;
			for (size_t word = 0; word < nrWords; ++word) parity ^= left[word] & right[word];
			return (Clifford::detail::Popcount(parity) & 1U) != 0;
		}

		static double L1(const std::complex<double>& value) noexcept
		{
			return std::abs(value.real()) + std::abs(value.imag());
		}

		// Apply a Pauli sum on a few physical qubits to one frame in one pass.
		// Each term's logical image P maps |b> to phase(b) |b xor x(P)>. The X
		// parts span 2^m label offsets, so the sum acts independently on each
		// coset of labels; each coset is visited once, gathered, multiplied by
		// the sum's 2^m x 2^m action, and written back, appending new labels.
		// A result within the rounding error of its own contributions is an
		// exact zero, so cancellations leave no residue; genuinely small
		// amplitudes, coming from small contributions, are kept.
		void ApplyPauliSum(ExtendedFrame& frame, const detail::LocalPauliSum& sum,
			const size_t* qubits, size_t nrLocalQubits)
		{
			auto& w = sumWorkspace;
			const auto& basis = frame.cliffordBasis;
			const size_t nrWords = frame.signs.GetNrWords();
			const size_t nrTerms = sum.Terms();
			const size_t generatorBase = nrTerms;
			EnsureRows(w.rows, nrTerms + std::min(nrTerms, PauliSumWorkspace::MaxGenerators),
				frame.GetNrQubits());

			// Logical images of the terms, with i^|x&z| and their signs folded
			// into the coefficients.
			w.terms.clear();
			sum.ForEachTerm([&](detail::LocalPauli pauli, std::complex<double> coefficient) {
				auto image = w.rows[w.terms.size()];
				image.Clear();
				for (size_t qubit = 0; qubit < nrLocalQubits; ++qubit)
				{
					const char factor = pauli.Factor(qubit);
					if (factor != 'I') basis.MultiplyImageInto(image, qubits[qubit], factor);
				}
				unsigned nrY = 0;
				bool hasZ = false;
				for (size_t word = 0; word < nrWords; ++word)
				{
					nrY += Clifford::detail::Popcount(image.X.words[word] & image.Z.words[word]);
					hasZ |= image.Z.words[word] != 0;
				}
				const auto phase = detail::LocalPauliSum::IPower(nrY);
				w.terms.push_back({ (image.PhaseSign ? -coefficient : coefficient) * phase, 0, hasZ });
			});

			// Generators of the X parts in reduced echelon form: each reduced row
			// keeps its pivot alone among them and records which generators
			// (original term X parts) it combines.
			std::array<size_t, PauliSumWorkspace::MaxGenerators> pivots{}, generatorTerms{};
			std::array<unsigned, PauliSumWorkspace::MaxGenerators> combinations{};
			size_t nrGenerators = 0;
			for (size_t term = 0; term < nrTerms; ++term)
			{
				auto reduced = w.rows[generatorBase + nrGenerators];
				const auto image = w.rows[term];
				std::copy_n(image.X.words, nrWords, reduced.X.words);
				unsigned combination = 0;
				for (size_t g = 0; g < nrGenerators; ++g)
					if (reduced.X[pivots[g]])
					{
						const auto other = w.rows[generatorBase + g];
						for (size_t word = 0; word < nrWords; ++word) reduced.X.words[word] ^= other.X.words[word];
						combination ^= combinations[g];
					}
				size_t pivot = 0;
				if (!Clifford::detail::FindPivot(reduced, pivot))
				{
					w.terms[term].offset = combination;
					continue;
				}
				combination ^= 1U << nrGenerators;
				for (size_t g = 0; g < nrGenerators; ++g)
				{
					auto other = w.rows[generatorBase + g];
					if (!other.X[pivot]) continue;
					for (size_t word = 0; word < nrWords; ++word) other.X.words[word] ^= reduced.X.words[word];
					combinations[g] ^= combination;
				}
				pivots[nrGenerators] = pivot;
				combinations[nrGenerators] = combination;
				generatorTerms[nrGenerators] = term;
				w.terms[term].offset = size_t(1) << nrGenerators;
				++nrGenerators;
			}

			// Offset u flips the labels by the generators selected by its bits.
			const size_t nrOffsets = size_t(1) << nrGenerators;
			EnsureRows(w.offsets, nrOffsets, frame.GetNrQubits());
			for (size_t u = 0; u < nrOffsets; ++u)
			{
				auto offset = w.offsets[u];
				std::fill_n(offset.X.words, nrWords, ExtendedFrame::Word(0));
				for (size_t g = 0; g < nrGenerators; ++g)
					if ((u >> g) & 1U)
					{
						const auto generator = w.rows[generatorTerms[g]];
						for (size_t word = 0; word < nrWords; ++word) offset.X.words[word] ^= generator.X.words[word];
					}
			}
			w.offsetParity.assign(nrTerms * nrOffsets, 0);
			for (size_t term = 0; term < nrTerms; ++term)
				if (w.terms[term].hasZ)
					for (size_t u = 1; u < nrOffsets; ++u)
						w.offsetParity[term * nrOffsets + u] = OddParity(w.rows[term].Z.words,
							w.offsets[u].X.words, nrWords);

			const double tolerance = 8.0 * double(nrTerms + 1) * std::numeric_limits<double>::epsilon();
			const size_t originalSize = frame.GetFrameSize();
			if (nrGenerators == 0)
			{
				// A diagonal sum rescales every amplitude in place.
				for (size_t component = 0; component < originalSize; ++component)
				{
					const auto* label = frame.signs.LabelWords(component);
					std::complex<double> factor(0.0, 0.0);
					double bound = 0.0;
					for (size_t term = 0; term < nrTerms; ++term)
					{
						const auto& t = w.terms[term];
						const bool odd = t.hasZ && OddParity(w.rows[term].Z.words, label, nrWords);
						factor += odd ? -t.coefficient : t.coefficient;
						bound += L1(t.coefficient);
					}
					frame.amplitudes[component] = L1(factor) <= tolerance * bound
						? std::complex<double>(0.0, 0.0) : frame.amplitudes[component] * factor;
				}
				PruneComponents(frame);
				return;
			}

			// The coset loop reads per-term data from local arrays.
			constexpr size_t maxTerms = detail::LocalPauliSum::Keys;
			std::array<std::complex<double>, maxTerms> coefficients, negatedCoefficients;
			std::array<double, maxTerms> magnitudes;
			std::array<size_t, maxTerms> termOffsets;
			for (size_t term = 0; term < nrTerms; ++term)
			{
				coefficients[term] = w.terms[term].coefficient;
				negatedCoefficients[term] = -coefficients[term];
				magnitudes[term] = L1(coefficients[term]);
				termOffsets[term] = w.terms[term].offset;
			}
			std::array<size_t, PauliSumWorkspace::MaxOffsets> members;
			std::array<std::complex<double>, PauliSumWorkspace::MaxOffsets> amplitudes, results;
			std::array<double, PauliSumWorkspace::MaxOffsets> amplitudeMagnitudes;
			std::array<bool, maxTerms> labelParity;

			const size_t notFound = PackedComponentIndex::NotFound;
			frame.EnsureComponentIndex(2 * originalSize);
			w.visited.assign(originalSize, 0);
			w.label.resize(nrWords);
			bool appended = false;
			for (size_t component = 0; component < originalSize; ++component)
			{
				if (w.visited[component]) continue;
				// The coset is the labels L xor offset; L is copied because
				// appending may reallocate the label storage.
				std::copy_n(frame.signs.LabelWords(component), nrWords, w.label.data());
				for (size_t u = 0; u < nrOffsets; ++u)
				{
					const size_t member = u == 0 ? component
						: frame.FindXorComponent(w.label.data(), w.offsets[u].X.words);
					members[u] = member;
					if (member == notFound)
						amplitudes[u] = 0.0;
					else
					{
						amplitudes[u] = frame.amplitudes[member];
						w.visited[member] = 1;
					}
					amplitudeMagnitudes[u] = L1(amplitudes[u]);
				}

				// out[v] = sum_j c_j (-1)^(z_j . (L xor offset_u)) in[u], u = v xor offset_j.
				// The bound sums the contributions' magnitudes for the rounding test.
				for (size_t term = 0; term < nrTerms; ++term)
					labelParity[term] = w.terms[term].hasZ
						&& OddParity(w.rows[term].Z.words, w.label.data(), nrWords);
				for (size_t v = 0; v < nrOffsets; ++v)
				{
					std::complex<double> result(0.0, 0.0);
					double bound = 0.0;
					for (size_t term = 0; term < nrTerms; ++term)
					{
						const size_t u = v ^ termOffsets[term];
						const bool odd = labelParity[term] != bool(w.offsetParity[term * nrOffsets + u]);
						result += (odd ? negatedCoefficients[term] : coefficients[term]) * amplitudes[u];
						bound += magnitudes[term] * amplitudeMagnitudes[u];
					}
					results[v] = L1(result) <= tolerance * bound
						? std::complex<double>(0.0, 0.0) : result;
				}

				for (size_t v = 0; v < nrOffsets; ++v)
				{
					if (members[v] != notFound)
						frame.amplitudes[members[v]] = results[v];
					else if (results[v] != std::complex<double>(0.0, 0.0))
					{
						frame.signs.AppendXor(w.label.data(), w.offsets[v].X.words);
						frame.amplitudes.push_back(results[v]);
						appended = true;
					}
				}
			}

			// Appended labels are not indexed.
			if (appended) frame.InvalidateComponentIndex();
			PruneComponents(frame);
		}

		static double PauliExpectation(const ExtendedFrame& frame,
			const PauliAction& action)
		{
			if (!HasPauliFlip(action))
			{
				double expectation = 0.0;
				for (size_t component = 0; component < frame.GetFrameSize(); ++component)
					expectation += std::norm(frame.amplitudes[component])
						* PauliPhase(frame.signs.LabelWords(component), action).real();
				return expectation;
			}

			// An off-diagonal Pauli maps the only occupied label to an orthogonal
			// label, so a one-component expectation is exactly zero.
			if (frame.GetFrameSize() == 1)
				return 0.0;

			frame.EnsureComponentIndex(frame.GetFrameSize());
			const size_t notFound = PackedComponentIndex::NotFound;

			std::complex<double> expectation(0.0, 0.0);
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
			{
				const auto* sourceLabel = frame.signs.LabelWords(component);
				const size_t target = frame.FindXorComponent(sourceLabel,
					action.flipMask);
				if (target == notFound) continue;

				expectation += std::conj(frame.amplitudes[target])
					* PauliPhase(sourceLabel, action)
					* frame.amplitudes[component];
			}

			return expectation.real();
		}

		double PauliExpectation(const ExtendedFrame& frame,
			const CliffordBasisMap::Row& pauli) const
		{
			CompilePauliAction(frame, pauli, pauliActionWorkspace);
			return PauliExpectation(frame, pauliActionWorkspace);
		}

		bool MeasureSingleStabilizer(size_t qubit, const bool* forcedOutcome = nullptr)
		{
			auto& frame = frames.front();
			const auto observable = frame.cliffordBasis.ImageZ(qubit);
			const auto* label = frame.signs.LabelWords(0);
			size_t pivot = 0;
			if (!Clifford::detail::FindPivot(observable, pivot))
			{
				// Q = (-1)^sign Z^z has eigenvalue (-1)^(sign + z.b) on |b>.
				unsigned parity = 0;
				for (size_t word = 0; word < observable.Words(); ++word)
					parity ^= Clifford::detail::Popcount(observable.Z.words[word] & label[word]);
				return bool(observable.PhaseSign) != bool(parity & 1U);
			}

			const bool outcome = forcedOutcome ? *forcedOutcome : (dist(gen) < 0.5);
			if (frame.signs.Get(0, pivot))
			{
				// Let Q = U^dagger Z_q U and
				// Q|b> = phase(b)|b xor x>.  When b contains the pivot,
				// its normalized coefficient in the measured basis is
				// (-1)^outcome phase(b) a_b.  Preserve that phase explicitly;
				// it is global for one component, but becomes relative if frames
				// are joined later.  Read the packed Pauli directly so this fast
				// path does not compile a PauliAction.
				size_t nrY = 0;
				unsigned labelParity = 0;
				for (size_t word = 0; word < observable.Words(); ++word)
				{
					const auto xBits = observable.X.words[word];
					const auto zBits = observable.Z.words[word];
					nrY += Clifford::detail::Popcount(xBits & zBits);
					labelParity ^= Clifford::detail::Popcount(zBits & label[word]) & 1U;
				}
				bool negate = bool(observable.PhaseSign) != outcome;
				if (labelParity != 0) negate = !negate;

				std::complex<double> phase;
				switch (nrY % 4)
				{
				case 0: phase = { 1.0, 0.0 }; break;
				case 1: phase = { 0.0, 1.0 }; break;
				case 2: phase = { -1.0, 0.0 }; break;
				default: phase = { 0.0, -1.0 }; break;
				}
				frame.amplitudes.front() *= negate ? -phase : phase;
			}
			frame.cliffordBasis.MeasureZ(qubit, pivot, outcome,
				frame.signs.LabelWords(0), enableMultithreading);
			frame.InvalidateComponentIndex();
			return outcome;
		}

		void ApplyAxisRotation(size_t physicalQubit, double angle, char axis)
		{
			if (!std::isfinite(angle))
				throw std::invalid_argument("Rotation angle must be finite");
			if (angle == 0.0)
				return;
			long long quarterTurns = 0;
			if (detail::TryGetQuarterTurns(angle, quarterTurns))
			{
				ApplyCliffordRotation(physicalQubit, axis, quarterTurns);
				return;
			}

			const double halfAngle = 0.5 * angle;
			auto rotation = detail::LocalPauliSum::Identity(std::cos(halfAngle));
			rotation.Add(detail::LocalPauli::Single(0, axis), { 0.0, -std::sin(halfAngle) });
			const size_t qubits[] = { physicalQubit };
			for (auto& frame : frames)
				ApplyPauliSum(frame, rotation, qubits, 1);
		}

		std::vector<ExtendedFrame> frames;
		std::vector<ExtendedFrame> savedFrames;
		ExtendedStabilizerApproximationPolicy approximationPolicy;
		ExtendedStabilizerApproximationPolicy savedApproximationPolicy;
		ExtendedStabilizerApproximationStatistics approximationStatistics;
		ExtendedStabilizerApproximationStatistics savedApproximationStatistics;
		mutable PauliAction pauliActionWorkspace;
		PauliSumWorkspace sumWorkspace;
		detail::LocalGate localGate{ 1 };

		std::mt19937_64 gen;
		std::uniform_real_distribution<double> dist;
		bool enableMultithreading = true;
	};

}
