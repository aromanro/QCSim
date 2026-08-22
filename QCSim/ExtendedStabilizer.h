#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <memory>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Frame.h"

namespace QC {

	enum class ExtendedStabilizerApproximationMode : unsigned char {
		Exact = 0,
		Approximate
	};

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
	class ExtendedStabilizer {
	public:
		ExtendedStabilizer() = delete;

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

		void ApplySdg(size_t qubit)
		{
			ValidateQubit(qubit);
			ApplyZ(qubit);
			ApplyS(qubit);
		}

		void ApplyK(size_t qubit)
		{
			ValidateQubit(qubit);
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		void ApplySx(size_t qubit)
		{
			ValidateQubit(qubit);
			ApplyZ(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyZ(qubit);
			ApplyS(qubit);
		}

		void ApplySxDag(size_t qubit)
		{
			ValidateQubit(qubit);
			ApplyS(qubit);
			ApplyH(qubit);
			ApplyS(qubit);
		}

		void ApplyCY(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			ApplyZ(target);
			ApplyS(target);
			ApplyCX(target, control);
			ApplyS(target);
		}

		void ApplyCZ(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			ApplyH(target);
			ApplyCX(target, control);
			ApplyH(target);
		}

		void ApplySwap(size_t qubit1, size_t qubit2)
		{
			ValidateTwoQubits(qubit1, qubit2);
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
		}

		void ApplyISwap(size_t qubit1, size_t qubit2)
		{
			ValidateTwoQubits(qubit1, qubit2);
			ApplyS(qubit1);
			ApplyH(qubit1);
			ApplyS(qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyCX(qubit1, qubit2);
			ApplyH(qubit2);
		}

		void ApplyISwapDag(size_t qubit1, size_t qubit2)
		{
			ValidateTwoQubits(qubit1, qubit2);
			ApplyH(qubit2);
			ApplyCX(qubit1, qubit2);
			ApplyCX(qubit2, qubit1);
			ApplyZ(qubit2);
			ApplyS(qubit2);
			ApplyH(qubit1);
			ApplyZ(qubit1);
			ApplyS(qubit1);
		}

		void ApplyCX(size_t target, size_t control)
		{
			ValidateTwoQubits(target, control);
			for (auto& frame : frames)
				frame.cliffordBasis.ApplyCX(target, control);
		}

		void ApplyRx(size_t qubit, double angle)
		{
			ValidateQubit(qubit);
			ApplyPackedAxisRotation(qubit, angle, true);
		}

		void ApplyRy(size_t qubit, double angle)
		{
			PauliStringXZ pauli(GetNrQubits());
			ValidateQubit(qubit);
			pauli.SetY(qubit);
			ApplyPauliRotation(pauli, angle);
		}

		void ApplyRz(size_t qubit, double angle)
		{
			ValidateQubit(qubit);
			ApplyPackedAxisRotation(qubit, angle, false);
		}

		bool Measure(size_t qubit)
		{
			ValidateQubit(qubit);
			AccountForMeasurementConditioning();
			auto& frame = frames.front();
			if (frame.GetFrameSize() == 1)
				return MeasureSingleStabilizer(qubit);

			const auto observable = BasisPauliForQubit(frame, qubit, false);
			CompilePauliAction(frame, observable, pauliActionWorkspace);
			const auto& action = pauliActionWorkspace;
			if (!HasPauliFlip(action))
				return MeasureDiagonalFrame(frame, action);
			return MeasureOffDiagonalFrame(frame, qubit, action);
		}

		double GetQubitProbability(size_t qubit) const
		{
			ValidateQubit(qubit);
			const auto& frame = frames.front();
			const auto observable = BasisPauliForQubit(frame, qubit, false);
			return ClampProbability(0.5 * (1.0 - PauliExpectation(frames.front(), observable)));
		}

		double ExpectationValue(const std::string& pauliString) const
		{
			if (pauliString.size() > GetNrQubits())
				throw std::invalid_argument("Pauli string is longer than the register");

			PauliStringXZ physicalPauli(GetNrQubits());
			bool isIdentity = true;
			for (size_t qubit = 0; qubit < pauliString.size(); ++qubit)
			{
				const char pauli = static_cast<char>(std::toupper(static_cast<unsigned char>(pauliString[qubit])));
				switch (pauli)
				{
				case 'I':
					break;
				case 'X':
					physicalPauli.SetX(qubit);
					isIdentity = false;
					break;
				case 'Y':
					physicalPauli.SetY(qubit);
					isIdentity = false;
					break;
				case 'Z':
					physicalPauli.SetZ(qubit);
					isIdentity = false;
					break;
				default:
					throw std::runtime_error("Invalid operator in the Pauli string");
				}
			}

			if (isIdentity) return 1.0;
			const auto& frame = frames.front();
			return PauliExpectation(frame, frame.cliffordBasis.TransformToBasis(physicalPauli));
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
			auto simulator = std::make_unique<ExtendedStabilizer>(
				GetNrQubits(), approximationPolicy);
			simulator->frames = frames;
			simulator->savedFrames = savedFrames;
			simulator->savedApproximationPolicy = savedApproximationPolicy;
			simulator->approximationStatistics = approximationStatistics;
			simulator->savedApproximationStatistics =
				savedApproximationStatistics;
			simulator->gen = gen;
			simulator->dist = dist;
			return simulator;
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
		void AccountForMeasurementConditioning() noexcept
		{
			if (approximationStatistics.traceDistanceErrorBound > 0.0)
				approximationStatistics.traceDistanceErrorBound = 1.0;
		}

		static void ValidateApproximationPolicy(
			const ExtendedStabilizerApproximationPolicy& policy)
		{
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

		static bool TryGetQuarterTurns(double angle, long long& quarterTurns)
		{
			const double pi = std::acos(-1.0);
			const double quarterTurnAngle = 0.5 * pi;
			const double turns = angle / quarterTurnAngle;
			if (std::abs(turns) >= static_cast<double>(
				std::numeric_limits<long long>::max()))
				return false;

			const auto roundedTurns = static_cast<long long>(std::llround(turns));
			const double reconstructed = static_cast<double>(roundedTurns)
				* quarterTurnAngle;
			if (angle != reconstructed)
				return false;

			quarterTurns = roundedTurns;
			return true;
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

		void ApplyCliffordRotation(const PauliStringXZ& physicalPauli,
			long long quarterTurns)
		{
			int mapTurns = static_cast<int>(quarterTurns % 4);
			if (mapTurns < 0) mapTurns += 4;
			const auto globalPhase = QuarterTurnGlobalPhase(quarterTurns);
			for (auto& frame : frames)
			{
				if (mapTurns == 1)
					frame.cliffordBasis.ApplyPauliQuarterTurn(physicalPauli);
				else if (mapTurns == 2)
				{
					frame.cliffordBasis.ApplyPauliQuarterTurn(physicalPauli);
					frame.cliffordBasis.ApplyPauliQuarterTurn(physicalPauli);
				}
				else if (mapTurns == 3)
					frame.cliffordBasis.ApplyPauliQuarterTurn(physicalPauli, true);

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

		static double ClampProbability(double probability)
		{
			return std::max(0.0, std::min(1.0, probability));
		}

		struct PauliAction {
			void OwnMasks(size_t words)
			{
				ownedMasks.resize(2 * words);
				std::fill(ownedMasks.begin(), ownedMasks.end(),
					ExtendedFrame::Word(0));
				flipMask = words == 0 ? nullptr : ownedMasks.data();
				phaseMask = words == 0 ? nullptr
					: ownedMasks.data() + words;
				nrWords = words;
			}

			void ReferenceMasks(const ExtendedFrame::Word* flip,
				const ExtendedFrame::Word* phase, size_t words) noexcept
			{
				flipMask = flip;
				phaseMask = phase;
				nrWords = words;
			}

			std::vector<ExtendedFrame::Word> ownedMasks;
			const ExtendedFrame::Word* flipMask = nullptr;
			const ExtendedFrame::Word* phaseMask = nullptr;
			size_t nrWords = 0;
			std::complex<double> basePhase{ 1.0, 0.0 };
		};

		static unsigned PopCount(ExtendedFrame::Word value) noexcept
		{
#ifdef _MSC_VER
			return static_cast<unsigned>(__popcnt64(value));
#else
			return static_cast<unsigned>(__builtin_popcountll(value));
#endif
		}

		static bool GetPackedBit(const ExtendedFrame::Word* words,
			size_t bit) noexcept
		{
			return (words[bit / PackedComponentLabels::BitsPerWord]
				& (ExtendedFrame::Word(1)
					<< (bit % PackedComponentLabels::BitsPerWord))) != 0;
		}

		static void SetPackedBit(ExtendedFrame::Word* words,
			size_t bit) noexcept
		{
			words[bit / PackedComponentLabels::BitsPerWord]
				|= ExtendedFrame::Word(1)
				<< (bit % PackedComponentLabels::BitsPerWord);
		}

		static CliffordBasisMap::PackedPauliView BasisPauliForQubit(
			const ExtendedFrame& frame, size_t qubit, bool useX)
		{
			if (useX)
				return frame.cliffordBasis.TransformXToBasisPacked(qubit);
			return frame.cliffordBasis.TransformZToBasisPacked(qubit);
		}

		static void CompilePauliAction(const ExtendedFrame& frame,
			const PauliStringXZWithSign& pauli, PauliAction& action)
		{
			const size_t nrQubits = frame.GetNrQubits();
			const size_t nrWords = frame.signs.GetNrWords();
			action.OwnMasks(nrWords);
			action.basePhase = { 1.0, 0.0 };

			// Component signs are logical labels, so the compiled action is exactly
			// the X/Z masks of U^dagger P U.
			for (size_t logical = 0; logical < nrQubits; ++logical)
			{
				if (pauli.X[logical])
					SetPackedBit(action.ownedMasks.data(), logical);
				if (pauli.Z[logical])
					SetPackedBit(action.ownedMasks.data() + nrWords, logical);
			}

			size_t nrY = 0;
			for (size_t logical = 0; logical < nrQubits; ++logical)
				if (pauli.X[logical] && pauli.Z[logical])
					++nrY;
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

		static void CompilePauliAction(const ExtendedFrame& frame,
			const CliffordBasisMap::PackedPauliView& pauli,
			PauliAction& action)
		{
			const size_t nrWords = frame.signs.GetNrWords();
			action.ReferenceMasks(pauli.GetXWords(), pauli.GetZWords(), nrWords);
			action.basePhase = { 1.0, 0.0 };

			size_t nrY = 0;
			for (size_t word = 0; word < nrWords; ++word)
				nrY += PopCount(action.flipMask[word] & action.phaseMask[word]);

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

			if (pauli.GetPhaseSign())
				action.basePhase = -action.basePhase;
		}

		static std::complex<double> PauliPhase(
			const ExtendedFrame::Word* sourceSigns,
			const PauliAction& action)
		{
			unsigned parity = 0;
			for (size_t word = 0; word < action.nrWords; ++word)
				parity ^= PopCount(sourceSigns[word] & action.phaseMask[word]) & 1U;
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
				if ((PopCount(action.phaseMask[word] & logicalLabel[word]) & 1U) != 0)
					outcome = !outcome;
			return outcome;
		}

		bool MeasureDiagonalFrame(ExtendedFrame& frame,
			const PauliAction& action)
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
			const bool outcome = normalizedProbabilityOne >= 1.0
				|| (normalizedProbabilityOne > 0.0
					&& dist(gen) < normalizedProbabilityOne);

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
			size_t physicalQubit, const PauliAction& action)
		{
			size_t pivot = 0;
			while (!GetPackedBit(action.flipMask, pivot))
				++pivot;

			frame.EnsureComponentIndex(frame.GetFrameSize());
			const size_t notFound = PackedComponentIndex::NotFound;
			const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);

			// Since x[pivot] is one, every logical label belongs to a unique
			// pair {r, r xor x} with r[pivot] == 0.  Process each such pair once,
			// including pairs whose zero-pivot member is absent from the frame.
			auto visitPairs = [&](const auto& visitor)
			{
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
					visitor(root, target);
				}
			};

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
			visitPairs([&](size_t root, size_t target)
			{
				const auto amplitudes = pairAmplitudes(root, target);
				const auto& outcomeZeroAmplitude = amplitudes.first;
				const auto& outcomeOneAmplitude = amplitudes.second;
				probabilityZero += std::norm(outcomeZeroAmplitude);
				probabilityOne += std::norm(outcomeOneAmplitude);
			});

			const double totalProbability = probabilityZero + probabilityOne;
			if (totalProbability <= 0.0)
				throw std::runtime_error(
					"Cannot measure a frame with zero total probability");
			const double normalizedProbabilityOne = ClampProbability(
				probabilityOne / totalProbability);
			const bool outcome = normalizedProbabilityOne >= 1.0
				|| (normalizedProbabilityOne > 0.0
					&& dist(gen) < normalizedProbabilityOne);

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
			visitPairs([&](size_t root, size_t target)
			{
				const auto amplitudes = pairAmplitudes(root, target);
				const auto amplitude = (outcome ? amplitudes.second : amplitudes.first)
					* inverseOutcomeNorm;
				if (std::norm(amplitude) == 0.0) return;
				collapsedAmplitudes.push_back(amplitude);
				if (root != notFound)
					collapsedSigns.Append(frame.signs.LabelWords(root));
				else
					collapsedSigns.AppendXor(frame.signs.LabelWords(target),
						action.flipMask);
				collapsedSigns.Set(collapsedSigns.size() - 1, pivot, outcome);
			});

			if (collapsedAmplitudes.empty())
				throw std::runtime_error(
					"Measurement removed every frame component");

			frame.cliffordBasis.RebaseZ(physicalQubit);
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
			const double toleranceSquared = approximate
				? approximationPolicy.amplitudeTolerance
					* approximationPolicy.amplitudeTolerance
				: 0.0;
			const size_t originalSize = frame.GetFrameSize();
			auto& retained = frame.componentOrderWorkspace;
			retained.clear();
			retained.reserve(originalSize);

			double totalWeight = 0.0;
			double largestWeight = -1.0;
			size_t largestComponent = 0;
			size_t positiveComponents = 0;
			for (size_t component = 0; component < originalSize; ++component)
			{
				const double weight = std::norm(frame.amplitudes[component]);
				totalWeight += weight;
				if (weight <= 0.0) continue;
				++positiveComponents;
				if (weight > largestWeight)
				{
					largestWeight = weight;
					largestComponent = component;
				}
			}

			if (totalWeight <= 0.0 || positiveComponents == 0)
				throw std::runtime_error("A frame operation cancelled every component");
			const double normalizedCutoff = toleranceSquared * totalWeight;
			for (size_t component = 0; component < originalSize; ++component)
			{
				const double weight = std::norm(frame.amplitudes[component]);
				if (weight > 0.0
					&& (!approximate || weight > normalizedCutoff))
					retained.push_back(component);
			}
			// Approximation is never allowed to erase the complete state.
			if (retained.empty()) retained.push_back(largestComponent);

			const size_t componentLimit = approximate
				? approximationPolicy.maxComponents : 0;
			if (componentLimit != 0 && retained.size() > componentLimit)
			{
				auto heavier = [&frame](size_t left, size_t right)
				{
					const double leftWeight = std::norm(frame.amplitudes[left]);
					const double rightWeight = std::norm(frame.amplitudes[right]);
					return leftWeight != rightWeight
						? leftWeight > rightWeight : left < right;
				};
				std::nth_element(retained.begin(),
					retained.begin() + componentLimit, retained.end(), heavier);
				retained.resize(componentLimit);
				std::sort(retained.begin(), retained.end());
			}

			double retainedWeight = 0.0;
			for (const size_t component : retained)
				retainedWeight += std::norm(frame.amplitudes[component]);
			if (retainedWeight <= 0.0)
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

			const size_t approximateDiscardedComponents = approximate
				? positiveComponents - retained.size() : 0;
			if (approximateDiscardedComponents == 0) return;

			const double inverseNorm = 1.0 / std::sqrt(retainedWeight);
			for (auto& amplitude : frame.amplitudes)
				amplitude *= inverseNorm;
			const double discardedWeight = std::max(0.0,
				totalWeight - retainedWeight);
			const double localDiscardedFraction = std::min(1.0,
				discardedWeight / totalWeight);
			approximationStatistics.cumulativeDiscardedWeight +=
				localDiscardedFraction;
			approximationStatistics.traceDistanceErrorBound = std::min(1.0,
				approximationStatistics.traceDistanceErrorBound
					+ std::sqrt(localDiscardedFraction));
			approximationStatistics.discardedComponents +=
				approximateDiscardedComponents;
			++approximationStatistics.pruningEvents;
		}

		void ApplyPauliCombination(ExtendedFrame& frame,
			const PauliAction& action,
			const std::complex<double>& identityCoefficient, const std::complex<double>& pauliCoefficient)
		{
			// A diagonal logical Pauli does not create new basis labels.  Updating
			// coefficients in place avoids all component copies and hash lookups for
			// common rotations such as Rz in the computational basis.
			if (!HasPauliFlip(action))
			{
				for (size_t component = 0; component < frame.GetFrameSize(); ++component)
					frame.amplitudes[component] *= identityCoefficient
						+ pauliCoefficient
							* PauliPhase(frame.signs.LabelWords(component), action);
				PruneComponents(frame);
				return;
			}

			// A Pauli permutes labels in disjoint two-element orbits. Process each
			// occupied pair once; only an orbit with one missing endpoint grows the
			// frame. The reusable flat index is updated as endpoints are appended.
			const size_t originalSize = frame.GetFrameSize();
			frame.amplitudes.reserve(2 * originalSize);
			frame.signs.reserve(2 * originalSize);
			frame.EnsureComponentIndex(2 * originalSize);
			const size_t notFound = PackedComponentIndex::NotFound;
			for (size_t component = 0; component < originalSize; ++component)
			{
				const auto* sourceLabel = frame.signs.LabelWords(component);
				const size_t target = frame.FindXorComponent(sourceLabel,
					action.flipMask);
				if (target != notFound)
				{
					if (target < component) continue;
					const auto sourceAmplitude = frame.amplitudes[component];
					const auto targetAmplitude = frame.amplitudes[target];
					frame.amplitudes[component] = identityCoefficient * sourceAmplitude
						+ pauliCoefficient
							* PauliPhase(frame.signs.LabelWords(target), action)
							* targetAmplitude;
					frame.amplitudes[target] = identityCoefficient * targetAmplitude
						+ pauliCoefficient * PauliPhase(sourceLabel, action)
							* sourceAmplitude;
					continue;
				}

				const auto sourceAmplitude = frame.amplitudes[component];
				const auto targetAmplitude = pauliCoefficient
					* PauliPhase(sourceLabel, action) * sourceAmplitude;
				frame.amplitudes[component] = identityCoefficient * sourceAmplitude;
				frame.signs.AppendXor(sourceLabel, action.flipMask);
				frame.amplitudes.push_back(targetAmplitude);
				frame.InsertIndexedComponent(frame.GetFrameSize() - 1);
			}

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
			const PauliStringXZWithSign& pauli) const
		{
			CompilePauliAction(frame, pauli, pauliActionWorkspace);
			return PauliExpectation(frame, pauliActionWorkspace);
		}

		double PauliExpectation(const ExtendedFrame& frame,
			const CliffordBasisMap::PackedPauliView& pauli) const
		{
			CompilePauliAction(frame, pauli, pauliActionWorkspace);
			return PauliExpectation(frame, pauliActionWorkspace);
		}

		bool MeasureSingleStabilizer(size_t qubit)
		{
			auto& frame = frames.front();
			const auto observable = BasisPauliForQubit(frame, qubit, false);
			const size_t nrQubits = frame.GetNrQubits();
			size_t pivot = nrQubits;
			for (size_t word = 0; word < observable.GetNrWords(); ++word)
			{
				auto xBits = observable.GetXWords()[word];
				if (xBits == 0) continue;

				size_t bit = 0;
				while ((xBits & 1ULL) == 0)
				{
					xBits >>= 1;
					++bit;
				}
				pivot = 8 * sizeof(xBits) * word + bit;
				break;
			}

			if (pivot == nrQubits)
			{
				bool outcome = observable.GetPhaseSign();
				for (size_t logical = 0; logical < nrQubits; ++logical)
					if (observable.Z(logical) && frame.signs.Get(0, logical))
						outcome = !outcome;
				return outcome;
			}

			const bool outcome = dist(gen) < 0.5;
			if (frame.signs.Get(0, pivot))
			{
				// Let Q = U^dagger Z_q U and
				// Q|b> = phase(b)|b xor x>.  When b contains the pivot,
				// its normalized coefficient in the measured basis is
				// (-1)^outcome phase(b) a_b.  Preserve that phase explicitly;
				// it is global for one component, but becomes relative if frames
				// are joined later.  Read the packed Pauli directly so this fast
				// path does not allocate a PauliAction and its two masks.
				size_t nrY = 0;
				bool negate = observable.GetPhaseSign() != outcome;
				for (size_t logical = 0; logical < nrQubits; ++logical)
				{
					const bool hasX = observable.X(logical);
					const bool hasZ = observable.Z(logical);
					if (hasX && hasZ) ++nrY;
					if (hasZ && frame.signs.Get(0, logical)) negate = !negate;
				}

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
			frame.cliffordBasis.MeasureZ(qubit, outcome,
				frame.signs.LabelWords(0), frame.signs.GetNrWords());
			frame.InvalidateComponentIndex();
			return outcome;
		}

		void ApplyPauliRotation(const PauliStringXZ& physicalPauli, double angle)
		{
			if (!std::isfinite(angle))
				throw std::invalid_argument("Rotation angle must be finite");
			if (angle == 0.0)
				return;
			long long quarterTurns = 0;
			if (TryGetQuarterTurns(angle, quarterTurns))
			{
				ApplyCliffordRotation(physicalPauli, quarterTurns);
				return;
			}

			const double halfAngle = 0.5 * angle;
			const std::complex<double> pauliCoefficient(0.0, -std::sin(halfAngle));
			for (auto& frame : frames)
			{
				const auto pauli = frame.cliffordBasis.TransformToBasis(physicalPauli);
				CompilePauliAction(frame, pauli, pauliActionWorkspace);
				ApplyPauliCombination(frame, pauliActionWorkspace,
					std::cos(halfAngle), pauliCoefficient);
			}
		}

		void ApplyPackedAxisRotation(size_t physicalQubit, double angle, bool useX)
		{
			if (!std::isfinite(angle))
				throw std::invalid_argument("Rotation angle must be finite");
			if (angle == 0.0)
				return;
			long long quarterTurns = 0;
			if (TryGetQuarterTurns(angle, quarterTurns))
			{
				PauliStringXZ physicalPauli(GetNrQubits());
				if (useX) physicalPauli.SetX(physicalQubit);
				else physicalPauli.SetZ(physicalQubit);
				ApplyCliffordRotation(physicalPauli, quarterTurns);
				return;
			}

			const double halfAngle = 0.5 * angle;
			const std::complex<double> pauliCoefficient(0.0, -std::sin(halfAngle));
			for (auto& frame : frames)
			{
				const auto pauli = useX
					? frame.cliffordBasis.TransformXToBasisPacked(physicalQubit)
					: frame.cliffordBasis.TransformZToBasisPacked(physicalQubit);
				CompilePauliAction(frame, pauli, pauliActionWorkspace);
				ApplyPauliCombination(frame, pauliActionWorkspace,
					std::cos(halfAngle), pauliCoefficient);
			}
		}

		std::vector<ExtendedFrame> frames;
		std::vector<ExtendedFrame> savedFrames;
		ExtendedStabilizerApproximationPolicy approximationPolicy;
		ExtendedStabilizerApproximationPolicy savedApproximationPolicy;
		ExtendedStabilizerApproximationStatistics approximationStatistics;
		ExtendedStabilizerApproximationStatistics savedApproximationStatistics;
		mutable PauliAction pauliActionWorkspace;

		std::mt19937 gen;
		std::uniform_real_distribution<double> dist;
	};

}
