#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Frame.h"

namespace QC {

	// A stabilizer frame represents the state in an orthonormal basis obtained
	// from the computational basis by a Clifford circuit. Clifford gates rotate
	// that basis, while non-Clifford rotations update the coefficients inside it.
	// Multiple ExtendedFrame objects are retained for a future multiframe implementation;
	// the coherent operations below intentionally implement one frame.
	class ExtendedStabilizer {
	public:
		ExtendedStabilizer() = delete;

		explicit ExtendedStabilizer(size_t nrQubits)
			: gen(std::random_device{}()), dist(0.0, 1.0)
		{
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
			auto& frame = frames.front();
			if (frame.GetFrameSize() == 1)
				return MeasureSingleStabilizer(qubit);

			const auto observable = BasisPauliForQubit(frame, qubit, false);
			const auto action = CompilePauliAction(frame, observable);
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
		}

		void RestoreState()
		{
			if (savedFrames.empty()) return;

			frames = savedFrames;
		}

		std::unique_ptr<ExtendedStabilizer> Clone() const
		{
			auto simulator = std::make_unique<ExtendedStabilizer>(GetNrQubits());
			simulator->frames = frames;
			simulator->savedFrames = savedFrames;
			return simulator;
		}

		const std::vector<ExtendedFrame>& GetFrames() const noexcept
		{
			return frames;
		}

	private:
		static constexpr double ComponentTolerance = 1E-14;

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
			ExtendedFrame::Signs flipMask;
			ExtendedFrame::Signs phaseMask;
			std::complex<double> basePhase;
		};

		static CliffordBasisMap::PackedPauliView BasisPauliForQubit(
			const ExtendedFrame& frame, size_t qubit, bool useX)
		{
			if (useX)
				return frame.cliffordBasis.TransformXToBasisPacked(qubit);
			return frame.cliffordBasis.TransformZToBasisPacked(qubit);
		}

		static PauliAction CompilePauliAction(const ExtendedFrame& frame,
			const PauliStringXZWithSign& pauli)
		{
			const size_t nrQubits = frame.GetNrQubits();
			PauliAction action{
				ExtendedFrame::Signs(nrQubits, false),
				ExtendedFrame::Signs(nrQubits, false),
				std::complex<double>(1.0, 0.0)
			};

			// Component signs are logical labels, so the compiled action is exactly
			// the X/Z masks of U^dagger P U.
			for (size_t logical = 0; logical < nrQubits; ++logical)
			{
				action.flipMask[logical] = pauli.X[logical];
				action.phaseMask[logical] = pauli.Z[logical];
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
			return action;
		}

		static PauliAction CompilePauliAction(const ExtendedFrame& frame,
			const CliffordBasisMap::PackedPauliView& pauli)
		{
			const size_t nrQubits = frame.GetNrQubits();
			PauliAction action{
				ExtendedFrame::Signs(nrQubits, false),
				ExtendedFrame::Signs(nrQubits, false),
				std::complex<double>(1.0, 0.0)
			};

			size_t nrY = 0;
			for (size_t logical = 0; logical < nrQubits; ++logical)
			{
				action.flipMask[logical] = pauli.X(logical);
				action.phaseMask[logical] = pauli.Z(logical);
				if (pauli.X(logical) && pauli.Z(logical))
					++nrY;
			}

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
			return action;
		}

		static ExtendedFrame::Signs PauliTargetSigns(
			const ExtendedFrame::Signs& sourceSigns,
			const PauliAction& action)
		{
			ExtendedFrame::Signs targetSigns(sourceSigns);
			for (size_t row = 0; row < targetSigns.size(); ++row)
				if (action.flipMask[row])
					targetSigns[row] = !targetSigns[row];
			return targetSigns;
		}

		static std::complex<double> PauliPhase(
			const ExtendedFrame::Signs& sourceSigns,
			const PauliAction& action)
		{
			bool negate = false;
			for (size_t row = 0; row < sourceSigns.size(); ++row)
				if (action.phaseMask[row] && sourceSigns[row])
					negate = !negate;
			return negate ? -action.basePhase : action.basePhase;
		}

		static std::string SignsKey(const ExtendedFrame::Signs& signs)
		{
			std::string key((signs.size() + 7) / 8, '\0');
			for (size_t bit = 0; bit < signs.size(); ++bit)
				if (signs[bit])
				{
					const size_t byte = bit / 8;
					const auto value = static_cast<unsigned char>(key[byte]) | static_cast<unsigned char>(1U << (bit % 8));
					key[byte] = static_cast<char>(value);
				}
			return key;
		}

		static bool HasPauliFlip(const PauliAction& action)
		{
			return std::find(action.flipMask.begin(), action.flipMask.end(), true)
				!= action.flipMask.end();
		}

		static bool DiagonalPauliOutcome(
			const ExtendedFrame::Signs& logicalLabel,
			const PauliAction& action)
		{
			// A diagonal logical Pauli has no Y factors, hence basePhase is
			// exactly +1 or -1.  Its eigenvalue on |b> is
			// (-1)^(baseMinus + z.b).
			bool outcome = action.basePhase.real() < 0.0;
			for (size_t logical = 0; logical < logicalLabel.size(); ++logical)
				if (action.phaseMask[logical] && logicalLabel[logical])
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
					frame.signs[component], action)
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
				if (DiagonalPauliOutcome(frame.signs[component], action) == outcome)
				{
					if (write != component)
					{
						frame.amplitudes[write] = frame.amplitudes[component];
						frame.signs[write] = std::move(frame.signs[component]);
					}
					frame.amplitudes[write] *= inverseNorm;
					++write;
				}

			frame.amplitudes.resize(write);
			frame.signs.resize(write);
			if (write == 0)
				throw std::runtime_error(
					"Measurement removed every frame component");
			return outcome;
		}

		bool MeasureOffDiagonalFrame(ExtendedFrame& frame,
			size_t physicalQubit, const PauliAction& action)
		{
			const size_t nrQubits = frame.GetNrQubits();
			size_t pivot = 0;
			while (!action.flipMask[pivot])
				++pivot;

			std::unordered_map<std::string, size_t> indices;
			indices.reserve(frame.GetFrameSize());
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
				indices.emplace(SignsKey(frame.signs[component]), component);

			struct MeasurementPair {
				ExtendedFrame::Signs root;
				std::complex<double> outcomeZeroAmplitude;
				std::complex<double> outcomeOneAmplitude;
			};
			std::vector<MeasurementPair> pairs;
			pairs.reserve(frame.GetFrameSize());
			double probabilityZero = 0.0;
			double probabilityOne = 0.0;
			const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);

			// Since x[pivot] is one, every logical label belongs to a unique
			// pair {r, r xor x} with r[pivot] == 0.  Process each such pair once,
			// including pairs whose zero-pivot member is absent from the frame.
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
			{
				ExtendedFrame::Signs root(frame.signs[component]);
				if (root[pivot])
				{
					root = PauliTargetSigns(root, action);
					if (indices.find(SignsKey(root)) != indices.end())
						continue;
				}

				const auto rootEntry = indices.find(SignsKey(root));
				const std::complex<double> rootAmplitude = rootEntry == indices.end()
					? std::complex<double>(0.0, 0.0)
					: frame.amplitudes[rootEntry->second];
				const auto target = PauliTargetSigns(root, action);
				const auto targetEntry = indices.find(SignsKey(target));
				const std::complex<double> targetAmplitude = targetEntry == indices.end()
					? std::complex<double>(0.0, 0.0)
					: frame.amplitudes[targetEntry->second];

				// If Q|b> = phase(b)|b xor x>, projection onto outcome m
				// gives the new-basis coefficient
				// (a_r + (-1)^m phase(r xor x) a_(r xor x))/sqrt(2).
				const auto mappedTargetAmplitude =
					PauliPhase(target, action) * targetAmplitude;
				const auto outcomeZeroAmplitude = inverseSqrtTwo
					* (rootAmplitude + mappedTargetAmplitude);
				const auto outcomeOneAmplitude = inverseSqrtTwo
					* (rootAmplitude - mappedTargetAmplitude);
				probabilityZero += std::norm(outcomeZeroAmplitude);
				probabilityOne += std::norm(outcomeOneAmplitude);
				pairs.push_back({ std::move(root), outcomeZeroAmplitude,
					outcomeOneAmplitude });
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

			std::vector<std::complex<double>> collapsedAmplitudes;
			std::vector<ExtendedFrame::Signs> collapsedSigns;
			collapsedAmplitudes.reserve(pairs.size());
			collapsedSigns.reserve(pairs.size());
			double retainedProbability = 0.0;
			for (auto& pair : pairs)
			{
				const auto amplitude = outcome
					? pair.outcomeOneAmplitude : pair.outcomeZeroAmplitude;
				if (std::abs(amplitude) <= ComponentTolerance)
					continue;
				retainedProbability += std::norm(amplitude);
				pair.root[pivot] = outcome;
				collapsedAmplitudes.push_back(amplitude);
				collapsedSigns.push_back(std::move(pair.root));
			}

			if (retainedProbability <= 0.0)
				throw std::runtime_error(
					"Measurement removed every frame component");
			const double inverseNorm = 1.0 / std::sqrt(retainedProbability);
			for (auto& amplitude : collapsedAmplitudes)
				amplitude *= inverseNorm;

			// The map update is independent of the supplied component label;
			// that argument only receives the corresponding new logical label.
			// With a zero-pivot representative the returned non-pivot bits are
			// unchanged, exactly matching the labels constructed above.
			ExtendedFrame::Signs dummyLabel(nrQubits, false);
			frame.cliffordBasis.MeasureZ(physicalQubit, outcome, dummyLabel);
			frame.amplitudes = std::move(collapsedAmplitudes);
			frame.signs = std::move(collapsedSigns);
			return outcome;
		}

		static void AccumulateComponent(std::vector<std::complex<double>>& amplitudes,
			std::vector<ExtendedFrame::Signs>& signs,
			std::unordered_map<std::string, size_t>& indices,
			const ExtendedFrame::Signs& componentSigns,
			const std::complex<double>& amplitude)
		{
			const std::string key = SignsKey(componentSigns);
			const auto found = indices.find(key);
			if (found != indices.end())
			{
				amplitudes[found->second] += amplitude;
				return;
			}

			indices.emplace(key, amplitudes.size());
			amplitudes.push_back(amplitude);
			signs.push_back(componentSigns);
		}

		static void RemoveZeroComponents(ExtendedFrame& frame)
		{
			constexpr double componentToleranceSquared =
				ComponentTolerance * ComponentTolerance;
			size_t write = 0;
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
				if (std::norm(frame.amplitudes[component]) > componentToleranceSquared)
				{
					if (write != component)
					{
						frame.amplitudes[write] = frame.amplitudes[component];
						frame.signs[write] = std::move(frame.signs[component]);
					}
					++write;
				}

			frame.amplitudes.resize(write);
			frame.signs.resize(write);
			if (write == 0)
				throw std::runtime_error("A frame operation cancelled every component");
		}

		static void ApplyPauliCombination(ExtendedFrame& frame,
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
						+ pauliCoefficient * PauliPhase(frame.signs[component], action);
				RemoveZeroComponents(frame);
				return;
			}

			std::vector<std::complex<double>> amplitudes;
			std::vector<ExtendedFrame::Signs> signs;
			std::unordered_map<std::string, size_t> indices;
			amplitudes.reserve(2 * frame.GetFrameSize());
			signs.reserve(2 * frame.GetFrameSize());
			indices.reserve(2 * frame.GetFrameSize());

			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
			{
				const auto& sourceSigns = frame.signs[component];
				const auto& sourceAmplitude = frame.amplitudes[component];
				AccumulateComponent(amplitudes, signs, indices, sourceSigns,
					identityCoefficient * sourceAmplitude);

				const auto targetSigns = PauliTargetSigns(sourceSigns, action);
				AccumulateComponent(amplitudes, signs, indices, targetSigns,
					pauliCoefficient * PauliPhase(sourceSigns, action) * sourceAmplitude);
			}

			frame.amplitudes = std::move(amplitudes);
			frame.signs = std::move(signs);
			RemoveZeroComponents(frame);
		}

		static double PauliExpectation(const ExtendedFrame& frame,
			const PauliAction& action)
		{
			if (!HasPauliFlip(action))
			{
				double expectation = 0.0;
				for (size_t component = 0; component < frame.GetFrameSize(); ++component)
					expectation += std::norm(frame.amplitudes[component])
						* PauliPhase(frame.signs[component], action).real();
				return expectation;
			}

			// An off-diagonal Pauli maps the only occupied label to an orthogonal
			// label, so a one-component expectation is exactly zero.
			if (frame.GetFrameSize() == 1)
				return 0.0;

			std::unordered_map<std::string, size_t> indices;
			indices.reserve(frame.GetFrameSize());
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
				indices.emplace(SignsKey(frame.signs[component]), component);

			std::complex<double> expectation(0.0, 0.0);
			for (size_t component = 0; component < frame.GetFrameSize(); ++component)
			{
				const auto targetSigns = PauliTargetSigns(frame.signs[component], action);
				const auto target = indices.find(SignsKey(targetSigns));
				if (target == indices.end()) continue;

				expectation += std::conj(frame.amplitudes[target->second])
					* PauliPhase(frame.signs[component], action)
					* frame.amplitudes[component];
			}

			return expectation.real();
		}

		static double PauliExpectation(const ExtendedFrame& frame,
			const PauliStringXZWithSign& pauli)
		{
			return PauliExpectation(frame, CompilePauliAction(frame, pauli));
		}

		static double PauliExpectation(const ExtendedFrame& frame,
			const CliffordBasisMap::PackedPauliView& pauli)
		{
			return PauliExpectation(frame, CompilePauliAction(frame, pauli));
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
					if (observable.Z(logical) && frame.signs.front()[logical])
						outcome = !outcome;
				return outcome;
			}

			const bool outcome = dist(gen) < 0.5;
			const auto& oldLabel = frame.signs.front();
			if (oldLabel[pivot])
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
					if (hasZ && oldLabel[logical]) negate = !negate;
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
			frame.cliffordBasis.MeasureZ(qubit, outcome, frame.signs.front());
			return outcome;
		}

		void ApplyPauliRotation(const PauliStringXZ& physicalPauli, double angle)
		{
			if (!std::isfinite(angle))
				throw std::invalid_argument("Rotation angle must be finite");
			if (angle == 0.0)
				return;

			const double halfAngle = 0.5 * angle;
			const std::complex<double> pauliCoefficient(0.0, -std::sin(halfAngle));
			for (auto& frame : frames)
			{
				const auto pauli = frame.cliffordBasis.TransformToBasis(physicalPauli);
				const auto action = CompilePauliAction(frame, pauli);
				ApplyPauliCombination(frame, action, std::cos(halfAngle), pauliCoefficient);
			}
		}

		void ApplyPackedAxisRotation(size_t physicalQubit, double angle, bool useX)
		{
			if (!std::isfinite(angle))
				throw std::invalid_argument("Rotation angle must be finite");
			if (angle == 0.0)
				return;

			const double halfAngle = 0.5 * angle;
			const std::complex<double> pauliCoefficient(0.0, -std::sin(halfAngle));
			for (auto& frame : frames)
			{
				const auto pauli = useX
					? frame.cliffordBasis.TransformXToBasisPacked(physicalQubit)
					: frame.cliffordBasis.TransformZToBasisPacked(physicalQubit);
				const auto action = CompilePauliAction(frame, pauli);
				ApplyPauliCombination(frame, action,
					std::cos(halfAngle), pauliCoefficient);
			}
		}

		std::vector<ExtendedFrame> frames;
		std::vector<ExtendedFrame> savedFrames;

		std::mt19937 gen;
		std::uniform_real_distribution<double> dist;
	};

}
