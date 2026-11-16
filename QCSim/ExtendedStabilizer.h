#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CliffordProbability.h"
#include "Frame.h"

// The matrix path uses vector instructions for its rounding screen and for
// the norm where the build enables AVX2, under the statevector kernels'
// switches; without them it runs the same tests in scalar code.
#if !defined(QCSIM_DISABLE_SIMD) && !defined(QCSIM_DISABLE_AVX2) && (defined(__AVX2__) || defined(_M_AVX2))
#include <immintrin.h>
#define QCSIM_AVX2 1
#endif
#include "LocalCliffordRecognition.h"
#include "LocalGateExtensions.h"
#include "LocalPauliSum.h"

namespace QC
{

enum class ExtendedStabilizerApproximationMode : unsigned char
{
    Exact = 0,
    Approximate
};

// Approximate pruning first removes components satisfying
// |amplitude| / ||state|| <= amplitudeTolerance, then keeps at most the
// largest maxComponents amplitudes. Equal magnitudes retain their earlier
// component index, at least one component always survives, and the result is
// renormalized. Switching back to Exact only disables future pruning: it
// cannot restore discarded components or clear the historical error bound.
struct ExtendedStabilizerApproximationPolicy
{
    ExtendedStabilizerApproximationMode mode = ExtendedStabilizerApproximationMode::Exact;
    double amplitudeTolerance = 0.0;
    // Zero means unlimited. An approximate cap must otherwise be at least one.
    size_t maxComponents = 0;

    static ExtendedStabilizerApproximationPolicy Exact() noexcept
    {
        return {};
    }

    static ExtendedStabilizerApproximationPolicy Approximate(double tolerance, size_t componentLimit = 0) noexcept
    {
        return {ExtendedStabilizerApproximationMode::Approximate, tolerance, componentLimit};
    }
};

struct ExtendedStabilizerApproximationStatistics
{
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
class ExtendedStabilizer
{
    struct CloneTag
    {
    };

  public:
    ExtendedStabilizer() = delete;
    ExtendedStabilizer(const ExtendedStabilizer &) = delete;
    ExtendedStabilizer &operator=(const ExtendedStabilizer &) = delete;

    explicit ExtendedStabilizer(size_t nrQubits)
        : ExtendedStabilizer(nrQubits, ExtendedStabilizerApproximationPolicy::Exact())
    {
    }

    ExtendedStabilizer(size_t nrQubits, const ExtendedStabilizerApproximationPolicy &policy)
        : approximationPolicy(policy), gen(std::random_device{}()), dist(0.0, 1.0)
    {
        ValidateApproximationPolicy(approximationPolicy);
        frames.emplace_back(nrQubits);
        unguardedQubits = nrQubits;
    }

    // Seeds the engine with the value itself. A std::seed_seq costs about ten
    // microseconds per call, too much for callers that reseed for every shot.
    void SetSeed(uint64_t theSeed)
    {
        gen.seed(theSeed);
        dist.reset();
    }

    size_t GetNrQubits() const
    {
        if (frames.empty())
            return 0;
        return frames.front().GetNrQubits();
    }

    void Reset(size_t nrQubits)
    {
        ClearUnitaryHistory();
        frames.clear();
        savedFrames.clear();
        approximationStatistics = {};
        savedApproximationStatistics = {};
        savedApproximationPolicy = approximationPolicy;
        frames.emplace_back(nrQubits);
        SetHistoryQubits(unitaryHistory.qubitMask);
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

    void setToBasisState(const std::vector<bool> &state)
    {
        const size_t nrQubits = GetNrQubits();
        if (state.size() > nrQubits)
            throw std::invalid_argument("Basis state bitvector exceeds qubit count");

        ClearUnitaryHistory();
        frames.clear();
        savedFrames.clear();
        approximationStatistics = {};
        savedApproximationStatistics = {};

        ExtendedFrame frame(nrQubits);
        frame.amplitudes = {{1.0, 0.0}};
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
        if (!IsRepresentable(State))
            return 0.0;
        return getBasisStateProbability(BasisStateBits(State));
    }

    double getBasisStateProbability(const std::vector<bool> &state) const
    {
        const size_t nrQubits = GetNrQubits();
        if (state.size() != nrQubits)
            throw std::invalid_argument("State size does not match qubit count");

        if (IsStabilizerState())
        {
            // One basis state U|b> is uniform on an affine support.
            const auto &frame = frames.front();
            Clifford::detail::BasisDistribution distribution;
            distribution.BuildFromInverse(frame.cliffordBasis.ZImages(), frame.signs.LabelWords(0));
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
    std::unordered_map<size_t, size_t> SampleCounts(const std::vector<size_t> &qubits, size_t shots)
    {
        if (qubits.size() > std::numeric_limits<size_t>::digits)
            throw std::invalid_argument("Use SampleCountsMany for outcomes wider than size_t");
        for (size_t qubit : qubits)
            ValidateQubit(qubit);
        std::unordered_map<size_t, size_t> counts;
        if (shots == 0 || qubits.empty())
            return counts;
        ForEachSample(qubits, shots,
                      [&](const std::vector<ExtendedFrame::Word> &bits) { ++counts[static_cast<size_t>(bits[0])]; });
        return counts;
    }

    std::unordered_map<std::vector<bool>, size_t> SampleCountsMany(const std::vector<size_t> &qubits, size_t shots)
    {
        for (size_t qubit : qubits)
            ValidateQubit(qubit);
        std::unordered_map<std::vector<bool>, size_t> counts;
        if (shots == 0 || qubits.empty())
            return counts;
        std::vector<bool> result(qubits.size());
        ForEachSample(qubits, shots, [&](const std::vector<ExtendedFrame::Word> &bits) {
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
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyH(qubit);
    }

    void ApplyS(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyS(qubit);
    }

    void ApplySdg(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplySdg(qubit);
    }

    void ApplyX(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyX(qubit);
    }

    void ApplyY(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyY(qubit);
    }

    void ApplyZ(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyZ(qubit);
    }

    void ApplySx(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplySx(qubit);
    }

    void ApplySxDag(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplySxDag(qubit);
    }

    void ApplyK(size_t qubit)
    {
        BeginNativeGate(qubit);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyK(qubit);
    }

    void ApplyCX(size_t target, size_t control)
    {
        BeginNativeGate(target, control);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyCX(target, control);
    }

    void ApplyCY(size_t target, size_t control)
    {
        BeginNativeGate(target, control);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyCY(target, control);
    }

    void ApplyCZ(size_t target, size_t control)
    {
        BeginNativeGate(target, control);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyCZ(target, control);
    }

    void ApplySwap(size_t qubit1, size_t qubit2)
    {
        BeginNativeGate(qubit1, qubit2);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplySwap(qubit1, qubit2);
    }

    void ApplyISwap(size_t qubit1, size_t qubit2)
    {
        BeginNativeGate(qubit1, qubit2);
        for (auto &frame : frames)
            frame.cliffordBasis.ApplyISwap(qubit1, qubit2);
    }

    void ApplyISwapDag(size_t qubit1, size_t qubit2)
    {
        BeginNativeGate(qubit1, qubit2);
        for (auto &frame : frames)
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
        RequireFiniteAngles({theta, phi, lambda});
        auto &gate = BeginGate(1);
        AddU(gate, 0, theta, phi, lambda);
        const size_t qubits[] = {qubit};
        ApplyLocalGate(gate, qubits);
    }

    // Controlled e^(i gamma) U(theta, phi, lambda).
    void ApplyCU(size_t target, size_t control, double theta, double phi, double lambda, double gamma = 0.0)
    {
        ValidateTwoQubits(target, control);
        RequireFiniteAngles({theta, phi, lambda, gamma});
        // Local qubit 0 is the control and 1 the target.
        auto &gate = BeginGate(2);
        gate.Rotate(0, 'Z', gamma).Rotate(1, 'Z', 0.5 * (lambda - phi)).Rotate(0, 'Z', 0.5 * (lambda + phi)).CX(1, 0);
        AddU(gate, 1, -0.5 * theta, 0.0, -0.5 * (phi + lambda));
        gate.CX(1, 0);
        AddU(gate, 1, 0.5 * theta, phi, 0.0);
        ApplyLocalGate(gate, target, control);
    }

    void ApplyCRx(size_t target, size_t control, double angle)
    {
        ValidateTwoQubits(target, control);
        RequireFiniteAngles({angle});
        auto &gate = BeginGate(2);
        gate.H(1).CX(1, 0).Rotate(1, 'Z', -0.5 * angle).CX(1, 0).Rotate(1, 'Z', 0.5 * angle).H(1);
        ApplyLocalGate(gate, target, control);
    }

    void ApplyCRy(size_t target, size_t control, double angle)
    {
        ValidateTwoQubits(target, control);
        RequireFiniteAngles({angle});
        auto &gate = BeginGate(2);
        gate.Rotate(1, 'Y', 0.5 * angle).CX(1, 0).Rotate(1, 'Y', -0.5 * angle).CX(1, 0);
        ApplyLocalGate(gate, target, control);
    }

    void ApplyCRz(size_t target, size_t control, double angle)
    {
        ValidateTwoQubits(target, control);
        RequireFiniteAngles({angle});
        auto &gate = BeginGate(2);
        gate.Rotate(1, 'Z', 0.5 * angle).CX(1, 0).Rotate(1, 'Z', -0.5 * angle).CX(1, 0);
        ApplyLocalGate(gate, target, control);
    }

    // Controlled Hadamard: Ry(pi/4) CZ Ry(-pi/4) on the target.
    void ApplyCH(size_t target, size_t control)
    {
        ValidateTwoQubits(target, control);
        const double eighthTurn = 0.25 * std::acos(-1.0);
        auto &gate = BeginGate(2);
        gate.Rotate(1, 'Y', -eighthTurn).CZ(1, 0).Rotate(1, 'Y', eighthTurn);
        ApplyLocalGate(gate, target, control);
    }

    // Controlled phase diag(1, 1, 1, e^(i lambda)).
    void ApplyCP(size_t target, size_t control, double lambda)
    {
        ValidateTwoQubits(target, control);
        RequireFiniteAngles({lambda});
        // Phases at multiples of pi/2 are exact; at odd multiples of pi the
        // gate is CZ.
        std::complex<double> phase;
        long long quarterTurns = 0;
        if (detail::TryGetQuarterTurns(lambda, quarterTurns))
        {
            const unsigned turns = unsigned(((quarterTurns % 4) + 4) % 4);
            if (turns == 0)
                return;
            if (turns == 2)
            {
                ApplyCZ(target, control);
                return;
            }
            phase = detail::LocalPauliSum::IPower(turns);
        }
        else
            phase = std::polar(1.0, lambda);
        auto &gate = BeginGate(2);
        gate.Multiply(
            Controlled(0, detail::LocalPauliSum::Projector(1, false) +
                              detail::LocalPauliSum::Projector(1, true) * detail::LocalPauliSum::Identity(phase)));
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
        auto &gate = BeginGate(3);
        gate.Multiply(Controlled(0, Controlled(1, detail::LocalPauliSum::Pauli(2, 'X'))));
        const size_t qubits[] = {control1, control2, target};
        ApplyLocalGate(gate, qubits);
    }

    void ApplyCSwap(size_t target1, size_t target2, size_t control)
    {
        ValidateThreeQubits(target1, target2, control);
        // Local qubit 0 is the control; SWAP = (II + XX + YY + ZZ) / 2.
        using Sum = detail::LocalPauliSum;
        Sum swap = Sum::Identity(0.5);
        for (const char pauli : {'X', 'Y', 'Z'})
            swap = swap + Sum::Pauli(1, pauli, 0.5) * Sum::Pauli(2, pauli);
        auto &gate = BeginGate(3);
        gate.Multiply(Controlled(0, swap));
        const size_t qubits[] = {control, target1, target2};
        ApplyLocalGate(gate, qubits);
    }

    // Matrix-gate arguments follow low-to-high local bits, as QubitRegister
    // does. The generic local expansion retains every nonzero coefficient.
    // A compact gate on at most three qubits acts through its dense matrix on
    // its own qubits: the target first, then the controls; a diagonal in the
    // order of its table bits.
    template <class AppliedGate> void ApplyGate(const AppliedGate &source)
    {
        if (const Gates::StructuredGate *compact = CompactOperation(source, 0))
        {
            const Gates::UnitaryLocalMatrix matrix = compact->DenseMatrix<Gates::UnitaryLocalMatrix>(3);
            size_t qubits[3] = {0, 0, 0};
            unsigned arity = 0;
            if (compact->kind != Gates::StructuredGate::Kind::Diagonal)
                qubits[arity++] = compact->target;
            for (const size_t qubit : compact->qubits)
                qubits[arity++] = qubit;
            ApplyMatrixGate(matrix, arity, qubits, [] { return Gates::UnitaryGateDescriptor{}; });
            return;
        }
        ApplyGate(source, source.getQubit1(), source.getQubit2(), source.getQubit3());
    }

    template <class Gate> void ApplyGate(const Gate &source, size_t q0, size_t q1 = 0, size_t q2 = 0)
    {
        const unsigned arity = static_cast<unsigned>(source.getQubitsNumber());
        const auto &matrix = source.getRawOperatorMatrix();
        if (!arity || arity > 3 || matrix.rows() != (1U << arity))
            throw std::invalid_argument("Gate requires a finite square matrix on one to three qubits");
        const size_t qubits[] = {q0, q1, q2};
        ApplyMatrixGate(matrix, arity, qubits, [&source] { return Gates::GetUnitaryGateDescriptor(source); });
    }

    template <class Gate> void ApplyControlledGate(const Gate &source, size_t target, size_t control)
    {
        ValidateTwoQubits(target, control);
        const Gates::ControlledGate<> controlled(source);
        ApplyGate(controlled, target, control);
    }

    template <class Gate>
    void ApplyDoublyControlledGate(const Gate &source, size_t target, size_t control1, size_t control2)
    {
        ValidateThreeQubits(target, control1, control2);
        const Gates::DoublyControlledGate<> controlled(source);
        ApplyGate(controlled, target, control1, control2);
    }

    void ApplyRXX(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'X', 'X', angle);
    }

    void ApplyRYY(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'Y', 'Y', angle);
    }

    void ApplyRZZ(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'Z', 'Z', angle);
    }

    void ApplyRZX(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'Z', 'X', angle);
    }

    void ApplyRXY(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'X', 'Y', angle);
    }

    void ApplyXXPlusYY(size_t a, size_t b, double theta, double beta = 0.)
    {
        ApplyXXYY(a, b, theta, beta, true);
    }

    void ApplyXXMinusYY(size_t a, size_t b, double theta, double beta = 0.)
    {
        ApplyXXYY(a, b, theta, beta, false);
    }

    void ApplyDCX(size_t a, size_t b)
    {
        ValidateTwoQubits(a, b);
        ApplyCX(b, a);
        ApplyCX(a, b);
    }

    void ApplyECR(size_t a, size_t b)
    {
        ValidateTwoQubits(a, b);
        ApplyS(a);
        ApplySx(b);
        ApplyCX(b, a);
        ApplyX(a);
    }

    bool Measure(size_t qubit, const bool *forcedOutcome = nullptr)
    {
        ValidateQubit(qubit);
        ClearUnitaryHistory();
        AccountForMeasurementConditioning();
        auto &frame = frames.front();
        if (frame.GetFrameSize() == 1)
            return MeasureSingleStabilizer(qubit, forcedOutcome);

        const auto observable = frame.cliffordBasis.ImageZ(qubit);
        CompilePauliAction(frame, observable, pauliActionWorkspace);
        const auto &action = pauliActionWorkspace;
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
        const auto &frame = frames.front();
        const auto observable = frame.cliffordBasis.ImageZ(qubit);
        return ClampProbability(0.5 * (1.0 - PauliExpectation(frames.front(), observable)));
    }

    double ExpectationValue(const std::string &pauliString) const
    {
        if (pauliString.size() > GetNrQubits())
            throw std::invalid_argument("Pauli string is longer than the register");

        // The logical images of the one-qubit factors are multiplied directly
        // into packed scratch storage; no physical Pauli string is built.
        // For a single basis state |b>, XOR accumulation decides most cases
        // without phase arithmetic: an off-diagonal image has expectation 0,
        // and a product of diagonal rows (-1)^s Z^z has (-1)^(s + z.b).
        const auto &frame = frames.front();
        const auto &basis = frame.cliffordBasis;
        const bool singleComponent = frame.GetFrameSize() == 1;
        bool isIdentity = true, allDiagonal = true;
        basis.ClearImage();
        for (size_t qubit = 0; qubit < pauliString.size(); ++qubit)
        {
            const char pauli = ToPauli(pauliString[qubit]);
            if (pauli == 'I')
                continue;
            isIdentity = false;
            if (singleComponent && basis.XorImage(qubit, pauli))
                allDiagonal = false;
        }

        if (isIdentity)
            return 1.0;
        if (singleComponent)
        {
            const auto image = basis.Image();
            if (image.HasX())
                return 0.0;
            if (allDiagonal)
            {
                const auto *label = frame.signs.LabelWords(0);
                ExtendedFrame::Word parity = 0;
                for (size_t word = 0; word < image.Words(); ++word)
                    parity ^= image.Z.words[word] & label[word];
                const bool negative = bool(image.PhaseSign) != ((Clifford::detail::Popcount(parity) & 1U) != 0);
                return negative ? -1.0 : 1.0;
            }
        }
        basis.ClearImage();
        for (size_t qubit = 0; qubit < pauliString.size(); ++qubit)
        {
            const char pauli = ToPauli(pauliString[qubit]);
            if (pauli != 'I')
                basis.MultiplyImage(qubit, pauli);
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
        if (savedFrames.empty())
            return;

        ClearUnitaryHistory();
        frames = savedFrames;
        approximationPolicy = savedApproximationPolicy;
        approximationStatistics = savedApproximationStatistics;
    }

    std::unique_ptr<ExtendedStabilizer> Clone() const
    {
        return std::unique_ptr<ExtendedStabilizer>(new ExtendedStabilizer(*this, CloneTag{}));
    }

    const std::vector<ExtendedFrame> &GetFrames() const noexcept
    {
        return frames;
    }

    const ExtendedStabilizerApproximationPolicy &GetApproximationPolicy() const noexcept
    {
        return approximationPolicy;
    }

    const ExtendedStabilizerApproximationStatistics &GetApproximationStatistics() const noexcept
    {
        return approximationStatistics;
    }

    double GetApproximationErrorBound() const noexcept
    {
        return approximationStatistics.traceDistanceErrorBound;
    }

    void SetApproximationPolicy(const ExtendedStabilizerApproximationPolicy &policy)
    {
        ValidateApproximationPolicy(policy);
        ClearUnitaryHistory();
        approximationPolicy = policy;
        if (policy.mode == ExtendedStabilizerApproximationMode::Approximate)
            for (auto &frame : frames)
                PruneComponents(frame);
    }

  private:
    struct ScaledNormAccumulator
    {
        void Add(double magnitude) noexcept
        {
            if (magnitude == 0.0)
                return;
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

    ExtendedStabilizer(const ExtendedStabilizer &other, CloneTag)
        : frames(other.frames), savedFrames(other.savedFrames), approximationPolicy(other.approximationPolicy),
          savedApproximationPolicy(other.savedApproximationPolicy),
          approximationStatistics(other.approximationStatistics),
          savedApproximationStatistics(other.savedApproximationStatistics), gen(other.gen), dist(other.dist),
          enableMultithreading(other.enableMultithreading)
    {
        // The copy starts without a unitary history.
        unguardedQubits = GetNrQubits();
    }

    bool IsRepresentable(size_t state) const noexcept
    {
        const size_t nrQubits = GetNrQubits();
        return nrQubits >= static_cast<size_t>(std::numeric_limits<size_t>::digits) || (state >> nrQubits) == 0;
    }

    std::vector<bool> BasisStateBits(size_t state) const
    {
        const size_t nrQubits = GetNrQubits();
        std::vector<bool> bits(nrQubits, false);
        const size_t width = std::min(nrQubits, static_cast<size_t>(std::numeric_limits<size_t>::digits));
        for (size_t qubit = 0; qubit < width; ++qubit)
            bits[qubit] = ((state >> qubit) & 1U) != 0;
        return bits;
    }

    // A single component is an ordinary stabilizer state U|b>.
    bool IsStabilizerState() const noexcept
    {
        return frames.size() == 1 && frames.front().GetFrameSize() == 1;
    }

    template <class Consumer> void ForEachSample(const std::vector<size_t> &qubits, size_t shots, Consumer consume)
    {
        std::vector<ExtendedFrame::Word> bits((qubits.size() + 63) / 64);
        if (IsStabilizerState())
        {
            // Sample the affine support directly. It consumes one fair draw per
            // random outcome, in the order of sequential measurements.
            const auto &frame = frames.front();
            Clifford::detail::BasisDistribution distribution;
            distribution.BuildFromInverse(frame.cliffordBasis.ZImages(), qubits, frame.signs.LabelWords(0));
            auto coin = [this](std::mt19937_64 &engine) { return dist(engine) < 0.5; };
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

    static void ValidateApproximationPolicy(const ExtendedStabilizerApproximationPolicy &policy)
    {
        if (policy.mode != ExtendedStabilizerApproximationMode::Exact &&
            policy.mode != ExtendedStabilizerApproximationMode::Approximate)
            throw std::invalid_argument("Unknown approximation mode");
        if (!std::isfinite(policy.amplitudeTolerance) || policy.amplitudeTolerance < 0.0)
            throw std::invalid_argument("Component amplitude tolerance must be finite and nonnegative");
        if (policy.mode == ExtendedStabilizerApproximationMode::Exact)
        {
            if (policy.amplitudeTolerance != 0.0 || policy.maxComponents != 0)
                throw std::invalid_argument("Exact mode cannot specify pruning parameters");
            return;
        }
        if (policy.amplitudeTolerance == 0.0 && policy.maxComponents == 0)
            throw std::invalid_argument("Approximate mode needs a tolerance or component limit");
    }

    static std::complex<double> QuarterTurnGlobalPhase(long long quarterTurns)
    {
        int phase = static_cast<int>(quarterTurns % 8);
        if (phase < 0)
            phase += 8;
        const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);
        switch (phase)
        {
        case 0:
            return {1.0, 0.0};
        case 1:
            return {inverseSqrtTwo, -inverseSqrtTwo};
        case 2:
            return {0.0, -1.0};
        case 3:
            return {-inverseSqrtTwo, -inverseSqrtTwo};
        case 4:
            return {-1.0, 0.0};
        case 5:
            return {-inverseSqrtTwo, inverseSqrtTwo};
        case 6:
            return {0.0, 1.0};
        default:
            return {inverseSqrtTwo, inverseSqrtTwo};
        }
    }

    void ApplyCliffordRotation(size_t physicalQubit, char axis, long long quarterTurns)
    {
        int mapTurns = static_cast<int>(quarterTurns % 4);
        if (mapTurns < 0)
            mapTurns += 4;
        const auto globalPhase = QuarterTurnGlobalPhase(quarterTurns);
        for (auto &frame : frames)
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

            for (auto &amplitude : frame.amplitudes)
                amplitude *= globalPhase;
        }
    }

    void ValidateQubit(size_t qubit) const
    {
        if (qubit >= GetNrQubits())
            throw std::out_of_range("Qubit index is outside the register");
    }

    // The entry of a native Clifford gate: its qubits are validated, and a
    // unitary history that covers them ends. While there is no history one
    // comparison per qubit settles both (see unguardedQubits); otherwise, and
    // for a qubit that is not valid, the gate takes the checks out of line.
    void BeginNativeGate(size_t qubit)
    {
        if (qubit >= unguardedQubits)
            GuardNativeGate(qubit);
    }

    void BeginNativeGate(size_t qubit1, size_t qubit2)
    {
        if (qubit1 >= unguardedQubits || qubit2 >= unguardedQubits || qubit1 == qubit2)
            GuardNativeGate(qubit1, qubit2);
    }

    EIGEN_DONT_INLINE void GuardNativeGate(size_t qubit)
    {
        ValidateQubit(qubit);
        InvalidateOverlappingUnitaryHistory(&qubit, 1);
    }

    EIGEN_DONT_INLINE void GuardNativeGate(size_t qubit1, size_t qubit2)
    {
        ValidateTwoQubits(qubit1, qubit2);
        const size_t qubits[] = {qubit1, qubit2};
        InvalidateOverlappingUnitaryHistory(qubits, 2);
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
        case 'I':
        case 'i':
            return 'I';
        case 'X':
        case 'x':
            return 'X';
        case 'Y':
        case 'y':
            return 'Y';
        case 'Z':
        case 'z':
            return 'Z';
        default:
            throw std::runtime_error("Invalid operator in the Pauli string");
        }
    }

    static double ClampProbability(double probability)
    {
        return std::max(0.0, std::min(1.0, probability));
    }

    struct PauliAction
    {
        void ReferenceMasks(const ExtendedFrame::Word *flip, const ExtendedFrame::Word *phase, size_t words) noexcept
        {
            flipMask = flip;
            phaseMask = phase;
            nrWords = words;
        }

        const ExtendedFrame::Word *flipMask = nullptr;
        const ExtendedFrame::Word *phaseMask = nullptr;
        size_t nrWords = 0;
        std::complex<double> basePhase{1.0, 0.0};
    };

    static void SetPackedBit(ExtendedFrame::Word *words, size_t bit) noexcept
    {
        words[bit / PackedComponentLabels::BitsPerWord] |= ExtendedFrame::Word(1)
                                                           << (bit % PackedComponentLabels::BitsPerWord);
    }

    static void CompilePauliAction(const ExtendedFrame &frame, const CliffordBasisMap::Row &pauli, PauliAction &action)
    {
        const size_t nrWords = frame.signs.GetNrWords();
        action.ReferenceMasks(pauli.X.words, pauli.Z.words, nrWords);
        action.basePhase = {1.0, 0.0};

        size_t nrY = 0;
        for (size_t word = 0; word < nrWords; ++word)
            nrY += Clifford::detail::Popcount(action.flipMask[word] & action.phaseMask[word]);

        switch (nrY % 4)
        {
        case 0:
            action.basePhase = {1.0, 0.0};
            break;
        case 1:
            action.basePhase = {0.0, 1.0};
            break;
        case 2:
            action.basePhase = {-1.0, 0.0};
            break;
        default:
            action.basePhase = {0.0, -1.0};
            break;
        }

        if (pauli.PhaseSign)
            action.basePhase = -action.basePhase;
    }

    static std::complex<double> PauliPhase(const ExtendedFrame::Word *sourceSigns, const PauliAction &action)
    {
        unsigned parity = 0;
        for (size_t word = 0; word < action.nrWords; ++word)
            parity ^= Clifford::detail::Popcount(sourceSigns[word] & action.phaseMask[word]) & 1U;
        const bool negate = parity != 0;
        return negate ? -action.basePhase : action.basePhase;
    }

    static bool HasPauliFlip(const PauliAction &action)
    {
        for (size_t word = 0; word < action.nrWords; ++word)
            if (action.flipMask[word] != 0)
                return true;
        return false;
    }

    static bool DiagonalPauliOutcome(const ExtendedFrame::Word *logicalLabel, const PauliAction &action)
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

    bool MeasureDiagonalFrame(ExtendedFrame &frame, const PauliAction &action, const bool *forcedOutcome = nullptr)
    {
        double probabilityZero = 0.0;
        double probabilityOne = 0.0;
        for (size_t component = 0; component < frame.GetFrameSize(); ++component)
        {
            double &probability =
                DiagonalPauliOutcome(frame.signs.LabelWords(component), action) ? probabilityOne : probabilityZero;
            probability += std::norm(frame.amplitudes[component]);
        }

        const double totalProbability = probabilityZero + probabilityOne;
        if (totalProbability <= 0.0)
            throw std::runtime_error("Cannot measure a frame with zero total probability");
        const double normalizedProbabilityOne = ClampProbability(probabilityOne / totalProbability);
        const bool outcome = forcedOutcome ? *forcedOutcome
                                           : (normalizedProbabilityOne >= 1.0 ||
                                              (normalizedProbabilityOne > 0.0 && dist(gen) < normalizedProbabilityOne));

        const double outcomeProbability = outcome ? probabilityOne : probabilityZero;
        if (outcomeProbability <= 0.0)
            throw std::runtime_error("Cannot normalize a zero-probability measurement outcome");
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
            throw std::runtime_error("Measurement removed every frame component");
        return outcome;
    }

    bool MeasureOffDiagonalFrame(ExtendedFrame &frame, size_t physicalQubit, const PauliAction &action,
                                 const bool *forcedOutcome = nullptr)
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
        auto &pairs = frame.measurementPairs;
        pairs.clear();
        pairs.reserve(frame.GetFrameSize());
        for (size_t component = 0; component < frame.GetFrameSize(); ++component)
        {
            const auto *componentLabel = frame.signs.LabelWords(component);
            const size_t partner = frame.FindXorComponent(componentLabel, action.flipMask);
            size_t root = component;
            size_t target = partner;
            if (frame.signs.Get(component, pivot))
            {
                if (partner != notFound)
                    continue;
                root = notFound;
                target = component;
            }
            pairs.emplace_back(root, target);
        }

        auto pairAmplitudes = [&](size_t root, size_t target) {
            const std::complex<double> rootAmplitude =
                root == notFound ? std::complex<double>(0.0, 0.0) : frame.amplitudes[root];
            const std::complex<double> targetAmplitude =
                target == notFound ? std::complex<double>(0.0, 0.0) : frame.amplitudes[target];

            // If Q|b> = phase(b)|b xor x>, projection onto outcome m
            // gives the new-basis coefficient
            // (a_r + (-1)^m phase(r xor x) a_(r xor x))/sqrt(2).
            const auto mappedTargetAmplitude =
                target == notFound ? std::complex<double>(0.0, 0.0)
                                   : PauliPhase(frame.signs.LabelWords(target), action) * targetAmplitude;
            return std::pair<std::complex<double>, std::complex<double>>{
                inverseSqrtTwo * (rootAmplitude + mappedTargetAmplitude),
                inverseSqrtTwo * (rootAmplitude - mappedTargetAmplitude)};
        };

        double probabilityZero = 0.0;
        double probabilityOne = 0.0;
        for (const auto &pair : pairs)
        {
            const size_t root = pair.first;
            const size_t target = pair.second;
            const auto amplitudes = pairAmplitudes(root, target);
            const auto &outcomeZeroAmplitude = amplitudes.first;
            const auto &outcomeOneAmplitude = amplitudes.second;
            probabilityZero += std::norm(outcomeZeroAmplitude);
            probabilityOne += std::norm(outcomeOneAmplitude);
        }

        const double totalProbability = probabilityZero + probabilityOne;
        if (totalProbability <= 0.0)
            throw std::runtime_error("Cannot measure a frame with zero total probability");
        const double normalizedProbabilityOne = ClampProbability(probabilityOne / totalProbability);
        const bool outcome = forcedOutcome ? *forcedOutcome
                                           : (normalizedProbabilityOne >= 1.0 ||
                                              (normalizedProbabilityOne > 0.0 && dist(gen) < normalizedProbabilityOne));

        const double outcomeProbability = outcome ? probabilityOne : probabilityZero;
        if (outcomeProbability <= 0.0)
            throw std::runtime_error("Cannot normalize a zero-probability measurement outcome");
        const double inverseOutcomeNorm = 1.0 / std::sqrt(outcomeProbability);
        auto &collapsedAmplitudes = frame.nextAmplitudes;
        auto &collapsedSigns = frame.nextSigns;
        collapsedAmplitudes.clear();
        collapsedSigns.clear();
        collapsedAmplitudes.reserve(frame.GetFrameSize());
        collapsedSigns.reserve(frame.GetFrameSize());
        for (const auto &pair : pairs)
        {
            const size_t root = pair.first;
            const size_t target = pair.second;
            const auto amplitudes = pairAmplitudes(root, target);
            const auto amplitude = (outcome ? amplitudes.second : amplitudes.first) * inverseOutcomeNorm;
            if (amplitude.real() == 0.0 && amplitude.imag() == 0.0)
                continue;
            collapsedAmplitudes.push_back(amplitude);
            if (root != notFound)
                collapsedSigns.Append(frame.signs.LabelWords(root));
            else
                collapsedSigns.AppendXor(frame.signs.LabelWords(target), action.flipMask);
            collapsedSigns.Set(collapsedSigns.size() - 1, pivot, outcome);
        }

        if (collapsedAmplitudes.empty())
            throw std::runtime_error("Measurement removed every frame component");

        frame.cliffordBasis.RebaseZ(physicalQubit, pivot, enableMultithreading);
        frame.amplitudes.swap(collapsedAmplitudes);
        frame.signs.swap(collapsedSigns);
        frame.InvalidateComponentIndex();
        PruneComponents(frame);
        return outcome;
    }

    void PruneComponents(ExtendedFrame &frame)
    {
        const bool approximate = approximationPolicy.mode == ExtendedStabilizerApproximationMode::Approximate;
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
                throw std::runtime_error("A frame operation cancelled every component");
            if (write != originalSize)
            {
                frame.amplitudes.resize(write);
                frame.signs.resize(write);
                frame.InvalidateComponentIndex();
            }
            return;
        }

        auto &retained = frame.componentOrderWorkspace;
        retained.clear();
        retained.reserve(originalSize);
        auto &magnitudes = frame.componentMagnitudeWorkspace;
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
            if (magnitude == 0.0)
                continue;
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
        const double amplitudeCutoff = approximationPolicy.amplitudeTolerance * totalNorm;
        for (size_t component = 0; component < originalSize; ++component)
        {
            const double magnitude = magnitudes[component];
            if (magnitude > amplitudeCutoff)
                retained.push_back(component);
        }
        // Approximation is never allowed to erase the complete state.
        if (retained.empty())
            retained.push_back(largestComponent);

        const size_t componentLimit = approximate ? approximationPolicy.maxComponents : 0;
        if (componentLimit != 0 && retained.size() > componentLimit)
        {
            auto heavier = [&magnitudes](size_t left, size_t right) {
                const double leftMagnitude = magnitudes[left];
                const double rightMagnitude = magnitudes[right];
                return leftMagnitude != rightMagnitude ? leftMagnitude > rightMagnitude : left < right;
            };
            std::nth_element(retained.begin(), retained.begin() + componentLimit, retained.end(), heavier);
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
            if (retainedPosition < retained.size() && retained[retainedPosition] == component)
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
                if (write == source)
                    continue;
                frame.amplitudes[write] = frame.amplitudes[source];
                frame.signs.CopyLabel(write, source);
            }
            frame.amplitudes.resize(retained.size());
            frame.signs.resize(retained.size());
            frame.InvalidateComponentIndex();
        }

        const size_t approximateDiscardedComponents = positiveComponents - retained.size();
        if (approximateDiscardedComponents == 0)
            return;

        const double inverseNorm = 1.0 / retainedNorm;
        for (auto &amplitude : frame.amplitudes)
            amplitude *= inverseNorm;
        const double localTraceDistance =
            std::min(1.0, std::nextafter(discardedNorm / totalNorm, std::numeric_limits<double>::infinity()));
        const double localDiscardedFraction = localTraceDistance * localTraceDistance;
        approximationStatistics.cumulativeDiscardedWeight += localDiscardedFraction;
        approximationStatistics.traceDistanceErrorBound =
            std::min(1.0, std::nextafter(approximationStatistics.traceDistanceErrorBound + localTraceDistance,
                                         std::numeric_limits<double>::infinity()));
        approximationStatistics.discardedComponents += approximateDiscardedComponents;
        ++approximationStatistics.pruningEvents;
    }

    static void AddU(detail::LocalGate &gate, size_t qubit, double theta, double phi, double lambda)
    {
        gate.Rotate(qubit, 'Z', lambda).Rotate(qubit, 'Y', theta).Rotate(qubit, 'Z', phi);
    }

    // |0><0| on the control plus |1><1| on the control times op.
    static detail::LocalPauliSum Controlled(size_t control, const detail::LocalPauliSum &op)
    {
        return detail::LocalPauliSum::Projector(control, false) + detail::LocalPauliSum::Projector(control, true) * op;
    }

    void ApplyControlledSquareRootX(size_t target, size_t control, bool inverse)
    {
        ValidateTwoQubits(target, control);
        const std::complex<double> plus(0.5, 0.5), minus(0.5, -0.5);
        auto &gate = BeginGate(2);
        gate.Multiply(Controlled(0, detail::LocalPauliSum::Identity(inverse ? minus : plus) +
                                        detail::LocalPauliSum::Pauli(1, 'X', inverse ? plus : minus)));
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
    detail::LocalGate &BeginGate(size_t qubits)
    {
        localGate.Reset(qubits);
        return localGate;
    }

    // Local qubit 0 is the control and 1 the target.
    void ValidateLocalQubits(const size_t *qubits, unsigned arity) const
    {
        for (unsigned q = 0; q < arity; ++q)
        {
            ValidateQubit(qubits[q]);
            for (unsigned previous = 0; previous < q; ++previous)
                if (qubits[previous] == qubits[q])
                    throw std::invalid_argument("Repeated gate qubit");
        }
    }

    void ApplyPairRotation(size_t a, size_t b, char first, char second, double angle)
    {
        ValidateTwoQubits(a, b);
        RequireFiniteAngles({angle});
        auto &gate = BeginGate(2);
        detail::AppendPairRotation(gate, first, second, angle);
        const size_t qubits[] = {a, b};
        ApplyLocalGate(gate, qubits);
    }

    void ApplyXXYY(size_t a, size_t b, double theta, double beta, bool plus)
    {
        ValidateTwoQubits(a, b);
        RequireFiniteAngles({theta, beta});
        // A local compiled circuit keeps exact Clifford factors separate from
        // the non-Clifford sum, including the iSWAP points of XX+YY.
        auto &gate = BeginGate(2);
        const size_t phaseQubit = plus ? 0 : 1;
        gate.Rotate(phaseQubit, 'Z', plus ? beta : -beta);
        detail::AppendPairRotation(gate, 'X', 'X', .5 * theta);
        detail::AppendPairRotation(gate, 'Y', 'Y', (plus ? .5 : -.5) * theta);
        gate.Rotate(phaseQubit, 'Z', plus ? -beta : beta);
        const size_t qubits[] = {a, b};
        ApplyLocalGate(gate, qubits);
    }

    void ApplyLocalGate(const detail::LocalGate &gate, size_t target, size_t control)
    {
        const size_t qubits[] = {control, target};
        ApplyLocalGate(gate, qubits);
    }

    // qubits maps the gate's local qubits to physical qubits.
    void ApplyLocalGate(const detail::LocalGate &gate, const size_t *qubits)
    {
        const bool hasSum = !gate.Sum().IsScalar();
        if (hasSum)
            ClearUnitaryHistory();
        else
            InvalidateOverlappingUnitaryHistory(qubits, gate.GetNrQubits());
        for (auto &frame : frames)
        {
            if (hasSum)
                ApplyPauliSum(frame, gate.Sum(), qubits, gate.GetNrQubits());
            gate.ReplayCliffords(frame.cliffordBasis, qubits);
        }
    }

    // A compact gate carries its own qubits; other gate types have none.
    template <class Gate> static auto CompactOperation(const Gate &gate, int) -> decltype(gate.getStructuredGate())
    {
        return gate.getStructuredGate();
    }

    template <class Gate> static const Gates::StructuredGate *CompactOperation(const Gate &, long)
    {
        return nullptr;
    }

    // The components of a frame, detached from its Clifford basis.
    struct FrameComponents
    {
        std::vector<std::complex<double>> amplitudes;
        PackedComponentLabels signs;
    };

    // A compiled matrix gate. It is immutable once cached, and the entries of
    // the unitary history share it.
    struct MatrixAction
    {
        std::shared_ptr<const detail::LocalGate> clifford;
        // The operator applied where every control selects; the complete
        // operator when the matrix shows no control.
        detail::LocalPauliSum sum;
        // Fingerprints of the matrix and, for a verified constructor, of its
        // adjoint: they only reject candidates when matching inverses.
        uint64_t matrixFingerprint = 0, inverseFingerprint = 0;
        // Local qubit and selecting value of each control visible in the matrix.
        std::array<uint8_t, 2> controlQubits{};
        std::array<bool, 2> controlValues{};
        unsigned nrControls = 0;
        // Factor of the part no control selects: one over the removed scale.
        double inactiveWeight = 1.;
        bool cannotAnnihilate = false;
        bool unitary = false;
        // A verified constructor, or a raw matrix that is unitary within
        // 1e-12: not a state filter for the trace-distance bound. This only
        // decides that bookkeeping; no coefficient depends on it.
        bool preservesBound = false;
    };

    struct CompiledAction
    {
        Gates::UnitaryGateDescriptor descriptor;
        std::shared_ptr<const MatrixAction> action;
    };

    // One cache per thread for every matrix type.
    static detail::LocalMatrixCache<CompiledAction> &MatrixActionCache()
    {
        thread_local detail::LocalMatrixCache<CompiledAction> cache;
        return cache;
    }

    // A control is visible in a matrix when the rows and columns where one
    // local qubit has the other value hold exactly the identity, and no entry
    // couples the two values of that qubit. Nothing is rounded: the test is
    // exact equality with zero and one.
    static bool FindVisibleControl(const Gates::UnitaryLocalMatrix &block, unsigned width, unsigned &position,
                                   bool &value)
    {
        const unsigned dimension = 1U << width;
        for (unsigned qubit = width; qubit-- > 0;)
            for (const bool selecting : {true, false})
            {
                bool visible = true;
                for (unsigned row = 0; row < dimension && visible; ++row)
                    for (unsigned col = 0; col < dimension && visible; ++col)
                    {
                        if ((bool((row >> qubit) & 1U) == selecting) && (bool((col >> qubit) & 1U) == selecting))
                            continue;
                        visible = block(row, col) == std::complex<double>(row == col ? 1. : 0.);
                    }
                if (visible)
                {
                    position = qubit;
                    value = selecting;
                    return true;
                }
            }
        return false;
    }

    template <class Matrix>
    static const std::shared_ptr<const MatrixAction> &CompileMatrixAction(
        const Matrix &matrix, unsigned arity, const Gates::UnitaryGateDescriptor &descriptor)
    {
        using LocalMatrix = Gates::UnitaryLocalMatrix;
        const auto key = detail::MatrixKey(matrix);
        auto &cache = MatrixActionCache();
        const unsigned set = cache.Set(key.hash, UINT64_C(0xD6E8FEB86659FD93));
        if (auto *slot = cache.Find(key, set, [&descriptor](const CompiledAction &compiled) {
                return (!compiled.descriptor && !descriptor) || compiled.descriptor == descriptor;
            }))
            return slot->value.action;
        auto compiled = std::make_shared<MatrixAction>();
        MatrixAction &action = *compiled;
        action.clifford = detail::FindCanonicalCliffordOwned(matrix, key);
        if (!action.clifford)
        {
            action.unitary = bool(descriptor) && Gates::VerifyUnitaryGateDescriptor(matrix, descriptor);
            // Compute once per compilation, not on each history comparison.
            // A fingerprint only rejects candidates; restoration still
            // requires an exact descriptor inverse or full matrix equality.
            action.matrixFingerprint = detail::MatrixKeyFingerprint(key, false);
            if (action.unitary)
                action.inverseFingerprint = detail::MatrixKeyFingerprint(key, true);
            // Generic finite matrices act as normalized state filters.
            // Remove a positive global scale before compiling their sums.
            double scale = 0.;
            for (Eigen::Index row = 0; row < matrix.rows(); ++row)
                for (Eigen::Index col = 0; col < matrix.cols(); ++col)
                {
                    const std::complex<double> value = matrix(row, col);
                    scale = std::max(scale, std::max(std::abs(value.real()), std::abs(value.imag())));
                }
            if (scale == 0.)
                throw std::domain_error("A zero operator annihilates the stabilizer state");
            const auto scaled = [scale](const LocalMatrix &source) {
                LocalMatrix result(source.rows(), source.cols());
                for (Eigen::Index row = 0; row < source.rows(); ++row)
                    for (Eigen::Index col = 0; col < source.cols(); ++col)
                    {
                        const std::complex<double> value = source(row, col);
                        result(row, col) = {value.real() / scale, value.imag() / scale};
                    }
                return result;
            };
            LocalMatrix block(matrix.rows(), matrix.cols());
            for (Eigen::Index row = 0; row < matrix.rows(); ++row)
                for (Eigen::Index col = 0; col < matrix.cols(); ++col)
                    block(row, col) = matrix(row, col);
            const LocalMatrix normalized = scaled(block);

            // Peel off the controls the matrix shows exactly. A controlled
            // operator is P0 (x) I + P1 (x) V: compiling V alone and applying
            // it only to the selected part keeps the unselected part
            // bit-identical, where the flat Pauli expansion would have to
            // cancel pairs of rounded products.
            std::array<uint8_t, 3> remaining{{0, 1, 2}};
            unsigned width = arity;
            while (width > 1 && action.nrControls < 2)
            {
                unsigned position = 0;
                bool value = false;
                if (!FindVisibleControl(block, width, position, value))
                    break;
                action.controlQubits[action.nrControls] = remaining[position];
                action.controlValues[action.nrControls++] = value;
                const unsigned reduced = 1U << (width - 1);
                const auto expand = [position, value](unsigned index) {
                    return (index & ((1U << position) - 1U)) | (unsigned(value) << position) |
                           ((index >> position) << (position + 1));
                };
                LocalMatrix selected(reduced, reduced);
                for (unsigned row = 0; row < reduced; ++row)
                    for (unsigned col = 0; col < reduced; ++col)
                        selected(row, col) = block(expand(row), expand(col));
                block = selected;
                for (unsigned qubit = position; qubit + 1 < width; ++qubit)
                    remaining[qubit] = remaining[qubit + 1];
                --width;
            }
            if (action.nrControls == 0)
                action.sum = detail::MatrixPauliSum(normalized, arity);
            else
            {
                // The selected block's sum, on the gate's own local qubits.
                detail::MatrixPauliSum(scaled(block), width)
                    .ForEachTerm([&](detail::LocalPauli pauli, std::complex<double> coefficient) {
                        detail::LocalPauli mapped;
                        for (unsigned qubit = 0; qubit < width; ++qubit)
                        {
                            mapped.x |= uint8_t(((pauli.x >> qubit) & 1U) << remaining[qubit]);
                            mapped.z |= uint8_t(((pauli.z >> qubit) & 1U) << remaining[qubit]);
                        }
                        action.sum.Add(mapped, coefficient);
                    });
                action.inactiveWeight = 1. / scale;
            }
            // This certificate controls rollback storage only. It neither
            // declares the matrix unitary nor changes any coefficient. A
            // comfortably positive Gershgorin bound on V†V proves that no
            // nonzero input can be annihilated; near-singular filters retain
            // the transactional path. The wide margin covers rounding in
            // the small (at most 8x8) Gram matrix.
            if (action.unitary)
                action.cannotAnnihilate = action.preservesBound = true;
            else
            {
                const LocalMatrix gram = normalized.adjoint() * normalized;
                const double squaredScale = scale * scale;
                double largestDiagonal = 0., lowerBound = std::numeric_limits<double>::infinity();
                bool nearUnitary = true;
                for (Eigen::Index row = 0; row < gram.rows(); ++row)
                {
                    const double diagonal = gram(row, row).real();
                    double offDiagonal = 0.;
                    for (Eigen::Index col = 0; col < gram.cols(); ++col)
                    {
                        if (col != row)
                            offDiagonal += std::abs(gram(row, col));
                        // The Gram matrix of the matrix as given, against the identity.
                        if (!(std::abs(gram(row, col) * squaredScale - (col == row ? 1. : 0.)) <= 1e-12))
                            nearUnitary = false;
                    }
                    largestDiagonal = std::max(largestDiagonal, diagonal);
                    lowerBound = std::min(lowerBound, diagonal - offDiagonal);
                }
                action.cannotAnnihilate = lowerBound > .25 * largestDiagonal;
                action.preservesBound = nearUnitary;
            }
        }
        auto &slot = cache.Insert(key, set);
        slot.value.descriptor = descriptor;
        slot.value.action = std::move(compiled);
        return slot.value.action;
    }

    template <class Matrix, class DescriptorOf>
    void ApplyMatrixGate(const Matrix &matrix, unsigned arity, const size_t *qubits, DescriptorOf descriptorOf)
    {
        ValidateLocalQubits(qubits, arity);
        // Common adjacent exact adjoints need no matrix hashing, compilation,
        // or replay. The saved matrix was certified when the anchor was made.
        if (unitaryHistory.size == 1 && RestoreAdjacentUnitaryInverse(matrix, qubits, arity))
            return;
        const Gates::UnitaryGateDescriptor descriptor = descriptorOf();
        const auto &compiled = CompileMatrixAction(matrix, arity, descriptor);
        const MatrixAction &action = *compiled;
        if (action.clifford)
        {
            ApplyLocalGate(*action.clifford, qubits);
            return;
        }
        if (unitaryHistory.size != 0 && RestoreUnitaryInverse(matrix, descriptor, action, qubits, arity))
            return;
        const bool remember = action.unitary && approximationPolicy.mode == ExtendedStabilizerApproximationMode::Exact;
        if (remember)
        {
            // One reusable anchor for a bounded sequence, rather than a full
            // component copy at every controlled gate. No amplitudes are
            // discarded and public queries always see the materialized state.
            if (unitaryHistory.size == UnitaryHistory::Capacity)
                ClearUnitaryHistory();
            if (unitaryHistory.size == 0)
            {
                if (unitaryHistory.entries.empty())
                    unitaryHistory.entries.resize(UnitaryHistory::Capacity);
                for (auto &frame : frames)
                {
                    frame.nextAmplitudes = frame.amplitudes;
                    frame.nextSigns = frame.signs;
                }
                unitaryHistory.statistics = approximationStatistics;
                unitaryHistory.firstAdjoint = matrix.adjoint();
                SetHistoryQubits(0);
                unitaryHistory.idle = 0;
            }
        }
        else
            ClearUnitaryHistory();
        bool normalized = false;
        try
        {
            normalized = ApplyMatrixSum(action, qubits, arity);
        }
        catch (...)
        {
            if (remember)
                ReplayUnitaryHistory();
            throw;
        }
        if (remember)
        {
            auto &entry = unitaryHistory.entries[unitaryHistory.size++];
            entry.descriptor = descriptor;
            entry.action = compiled;
            entry.arity = arity;
            entry.qubits = {qubits[0], qubits[1], qubits[2]};
            CanonicalizeUnitaryQubits(entry.qubits.data(), descriptor);
            SetHistoryQubits(unitaryHistory.qubitMask | QubitMask(qubits, arity));
            // The frames have unit norm now, so no amplitude has an L1 norm
            // above the square root of two. The history ends before any
            // other operation changes an amplitude, and the bound with it.
            amplitudeBound = normalized ? 1.5 : -1.0;
        }
    }

    // Returns whether the frames were left with unit norm (a scalar phase does not look).
    bool ApplyMatrixSum(const MatrixAction &action, const size_t *qubits, unsigned arity)
    {
        const auto &sum = action.sum;
        if (action.nrControls == 0 && sum.IsScalar())
        {
            std::complex<double> factor = 0.;
            sum.ForEachTerm([&](detail::LocalPauli, std::complex<double> value) { factor = value; });
            if (factor == std::complex<double>(0.))
                throw std::domain_error("A zero operator annihilates the stabilizer state");
            const auto phase = factor / std::abs(factor);
            for (auto &frame : frames)
                for (auto &amplitude : frame.amplitudes)
                    amplitude *= phase;
            return false;
        }
        // Only component labels and amplitudes change here. Back up those
        // into reusable workspaces, not the quadratic-size Clifford tableau.
        // Finish every allocation before changing any physical component.
        if (!action.cannotAnnihilate)
            for (auto &frame : frames)
            {
                frame.nextAmplitudes = frame.amplitudes;
                frame.nextSigns = frame.signs;
            }
        const auto previousStatistics = approximationStatistics;
        try
        {
            for (auto &frame : frames)
            {
                if (action.nrControls == 0)
                    ApplyPauliSum<true>(frame, sum, qubits, arity);
                else
                    ApplyMatrixTerms(frame, action, qubits, arity);
                bool vacant = false;
                const double squared = SquaredNorm(frame.amplitudes.data(), frame.amplitudes.size(), vacant);
                double length;
                if (squared >= std::numeric_limits<double>::min() && std::isfinite(squared))
                    length = std::sqrt(squared);
                else
                {
                    // Retain subnormal filtered states rather than mistaking
                    // underflowed squared magnitudes for annihilation.
                    ScaledNormAccumulator norm;
                    for (const auto &amplitude : frame.amplitudes)
                        norm.Add(std::abs(amplitude));
                    length = norm.Value();
                }
                if (!(length > 0.) || !std::isfinite(length))
                    throw std::domain_error("The operator annihilates the stabilizer state");
                // std::complex<double> has the interleaved real/imaginary
                // layout required here. Packet division preserves tiny
                // lengths without forming an overflowing reciprocal. A length
                // that is one to rounding is left alone: dividing by it would
                // only add a rounding to every amplitude, also to those of
                // components the gate did not touch.
                if (!(std::abs(length - 1.) <= 4. * std::numeric_limits<double>::epsilon()))
                    vacant = DivideAmplitudes(frame.amplitudes.data(), frame.amplitudes.size(), length);
                // Exact pruning only removes amplitudes that are exactly
                // zero, and the pass that last wrote them saw whether there
                // is one.
                if (vacant || approximationPolicy.mode != ExtendedStabilizerApproximationMode::Exact)
                    PruneComponents(frame);
            }
        }
        catch (...)
        {
            if (!action.cannotAnnihilate)
                for (auto &frame : frames)
                {
                    frame.amplitudes.swap(frame.nextAmplitudes);
                    frame.signs.swap(frame.nextSigns);
                    frame.InvalidateComponentIndex();
                }
            approximationStatistics = previousStatistics;
            throw;
        }
        // Conditioning a previously approximated state can amplify the
        // discarded component arbitrarily. Verified unitary constructors and
        // raw matrices that are unitary to rounding preserve the existing
        // bound; a state filter retains the same conservative bound used
        // after measurement/postselection.
        if (!action.preservesBound && previousStatistics.traceDistanceErrorBound > 0.)
            approximationStatistics.traceDistanceErrorBound = 1.;
        return true;
    }

    // Scratch amplitudes for one gate. A std::array of std::complex is zero
    // filled wherever it is declared, which costs more than a small gate on a
    // small frame; this storage is not, and every element is written before
    // it is read.
    template <size_t Size> struct ScratchAmplitudes
    {
        std::complex<double> *data() noexcept
        {
            return reinterpret_cast<std::complex<double> *>(storage);
        }

        const std::complex<double> *data() const noexcept
        {
            return reinterpret_cast<const std::complex<double> *>(storage);
        }

        std::complex<double> &operator[](size_t index) noexcept
        {
            return data()[index];
        }

        const std::complex<double> &operator[](size_t index) const noexcept
        {
            return data()[index];
        }

        // Stores an amplitude in one piece. A copy made in two halves cannot
        // hand its value to a read of the whole amplitude that follows: the
        // read then waits for both halves to reach the cache.
        void Set(size_t index, const std::complex<double> &value) noexcept
        {
#ifdef QCSIM_AVX2
            _mm_storeu_pd(reinterpret_cast<double *>(storage) + 2 * index, _mm_set_pd(value.imag(), value.real()));
#else
            data()[index] = value;
#endif
        }

        alignas(std::complex<double>) unsigned char storage[Size * sizeof(std::complex<double>)];
    };

    // The squared norm of a frame's amplitudes and whether one of them is
    // exactly zero, in one pass: the squares of the real and of the
    // imaginary parts are summed apart, in two partial sums over pairs of
    // amplitudes.
    static double SquaredNorm(const std::complex<double> *amplitudes, size_t size, bool &hasZero) noexcept
    {
        double real = 0.0, imaginary = 0.0;
        size_t index = 0;
        bool zero = false;
#ifdef QCSIM_AVX2
        if (size >= 2)
        {
            const double *values = reinterpret_cast<const double *>(amplitudes);
            const __m256d none = _mm256_setzero_pd();
            __m256d zeros = none;
            const auto squares = [values, none, &zeros](size_t at) {
                const __m256d pair = _mm256_loadu_pd(values + 2 * at);
                // An amplitude is zero where both of its parts are.
                const __m256d equal = _mm256_cmp_pd(pair, none, _CMP_EQ_OQ);
                zeros = _mm256_or_pd(zeros, _mm256_and_pd(equal, _mm256_permute_pd(equal, 5)));
                return _mm256_mul_pd(pair, pair);
            };
            __m256d first = squares(0);
            index = 2;
            if (size >= 4)
            {
                __m256d second = squares(2);
                for (index = 4; index + 4 <= size; index += 4)
                {
                    first = _mm256_add_pd(first, squares(index));
                    second = _mm256_add_pd(second, squares(index + 2));
                }
                first = _mm256_add_pd(first, second);
                if (index + 2 <= size)
                {
                    first = _mm256_add_pd(first, squares(index));
                    index += 2;
                }
            }
            const __m128d sums = _mm_add_pd(_mm256_castpd256_pd128(first), _mm256_extractf128_pd(first, 1));
            real = _mm_cvtsd_f64(sums);
            imaginary = _mm_cvtsd_f64(_mm_unpackhi_pd(sums, sums));
            zero = _mm256_movemask_pd(zeros) != 0;
        }
#endif
        for (; index < size; ++index)
        {
            const double realPart = amplitudes[index].real(), imaginaryPart = amplitudes[index].imag();
            real += realPart * realPart;
            imaginary += imaginaryPart * imaginaryPart;
            zero |= (realPart == 0.0) & (imaginaryPart == 0.0);
        }
        hasZero = zero;
        return real + imaginary;
    }

    // Divides the amplitudes by the norm, part by part, and returns whether
    // one of them is exactly zero afterwards: one that already was, or a
    // tiny one that underflowed.
    static bool DivideAmplitudes(std::complex<double> *amplitudes, size_t size, double length) noexcept
    {
        double *values = reinterpret_cast<double *>(amplitudes);
        size_t index = 0;
        bool zero = false;
#ifdef QCSIM_AVX2
        const __m256d divisor = _mm256_set1_pd(length), none = _mm256_setzero_pd();
        __m256d zeros = none;
        for (; index + 2 <= size; index += 2)
        {
            const __m256d pair = _mm256_div_pd(_mm256_loadu_pd(values + 2 * index), divisor);
            _mm256_storeu_pd(values + 2 * index, pair);
            const __m256d equal = _mm256_cmp_pd(pair, none, _CMP_EQ_OQ);
            zeros = _mm256_or_pd(zeros, _mm256_and_pd(equal, _mm256_permute_pd(equal, 5)));
        }
        zero = _mm256_movemask_pd(zeros) != 0;
#endif
        for (; index < size; ++index)
        {
            const double real = values[2 * index] / length, imaginary = values[2 * index + 1] / length;
            values[2 * index] = real;
            values[2 * index + 1] = imaginary;
            zero |= (real == 0.0) & (imaginary == 0.0);
        }
        return zero;
    }

    // An upper bound of the L1 norm of every amplitude of a frame: twice its
    // largest component. One pass over the frame finds it; while the frames
    // are known to be normalized (see amplitudeBound) there is none.
    double LargestAmplitude(const ExtendedFrame &frame) const noexcept
    {
        if (amplitudeBound > 0.0)
            return amplitudeBound;
        return 2.0 * Eigen::Map<const Eigen::ArrayXd>(reinterpret_cast<const double *>(frame.amplitudes.data()),
                                                      2 * frame.amplitudes.size())
                         .abs()
                         .maxCoeff();
    }

    // Collects, over the results of a coset, whether some real or imaginary
    // part is nonzero and at most the floor: only such a part can be the
    // rounding of its own sum. It reads each result as it is computed, not
    // from the array it was just stored to.
    struct RoundingScreen
    {
#ifdef QCSIM_AVX2
        explicit RoundingScreen(double floor) noexcept : limit(_mm_set1_pd(floor)), found(_mm_setzero_pd())
        {
        }

        void Add(const std::complex<double> &result) noexcept
        {
            const __m128d magnitude = _mm_andnot_pd(_mm_set1_pd(-0.0), _mm_set_pd(result.imag(), result.real()));
            found =
                _mm_or_pd(found, _mm_and_pd(_mm_cmple_pd(magnitude, limit), _mm_cmpgt_pd(magnitude, _mm_setzero_pd())));
        }

        bool Found() const noexcept
        {
            return _mm_movemask_pd(found) != 0;
        }

        __m128d limit, found;
#else
        explicit RoundingScreen(double floor) noexcept : limit(floor)
        {
        }

        void Add(const std::complex<double> &result) noexcept
        {
            const double real = std::abs(result.real()), imaginary = std::abs(result.imag());
            found |= ((real != 0.0) & (real <= limit)) | ((imaginary != 0.0) & (imaginary <= limit));
        }

        bool Found() const noexcept
        {
            return found;
        }

        double limit;
        bool found = false;
#endif
    };

    // Clear the parts of a coset's results that are the rounding of their own
    // sums: each real and imaginary part against the magnitudes of the real
    // products that entered it (and of the part that passed the controls,
    // when there is one). A factor of one or zero keeps the test free of
    // branches; whether a part cancels is not predictable. The passes only
    // come here for a coset whose screen found a candidate.
    //
    // Each product of magnitudes goes through an absolute value, which does
    // not change it, before it is summed: the compiler then cannot fuse a
    // product into the sum, and the bound is the same number in every build
    // and in every form of this test.
    EIGEN_STRONG_INLINE static void ClearRoundedParts(std::complex<double> *results,
                                                      const std::complex<double> *sources,
                                                      const std::complex<double> *rests, size_t nrOffsets,
                                                      const size_t *termOffsets, const double *realMagnitudes,
                                                      const double *imaginaryMagnitudes, size_t firstTerm,
                                                      size_t nrTerms, double rounding) noexcept
    {
#ifdef QCSIM_AVX2
        const __m128d sign = _mm_set1_pd(-0.0), one = _mm_set1_pd(1.0), scale = _mm_set1_pd(rounding);
        for (size_t v = 0; v < nrOffsets; ++v)
        {
            // The real and the imaginary sum, side by side.
            __m128d sums = rests ? _mm_andnot_pd(sign, _mm_loadu_pd(reinterpret_cast<const double *>(rests + v)))
                                 : _mm_setzero_pd();
            for (size_t term = firstTerm; term < nrTerms; ++term)
            {
                const __m128d parts = _mm_andnot_pd(
                    sign, _mm_loadu_pd(reinterpret_cast<const double *>(sources + (v ^ termOffsets[term]))));
                sums = _mm_add_pd(sums,
                                  _mm_add_pd(_mm_andnot_pd(sign, _mm_mul_pd(_mm_set1_pd(realMagnitudes[term]), parts)),
                                             _mm_andnot_pd(sign, _mm_mul_pd(_mm_set1_pd(imaginaryMagnitudes[term]),
                                                                            _mm_shuffle_pd(parts, parts, 1)))));
            }
            double *value = reinterpret_cast<double *>(results + v);
            const __m128d result = _mm_loadu_pd(value);
            const __m128d keep = _mm_and_pd(_mm_cmpgt_pd(_mm_andnot_pd(sign, result), _mm_mul_pd(scale, sums)), one);
            _mm_storeu_pd(value, _mm_mul_pd(result, keep));
        }
#else
        for (size_t v = 0; v < nrOffsets; ++v)
        {
            double realTerms = rests ? std::abs(rests[v].real()) : 0.0;
            double imaginaryTerms = rests ? std::abs(rests[v].imag()) : 0.0;
            for (size_t term = firstTerm; term < nrTerms; ++term)
            {
                const auto &source = sources[v ^ termOffsets[term]];
                const double realPart = std::abs(source.real()), imaginaryPart = std::abs(source.imag());
                realTerms +=
                    std::abs(realMagnitudes[term] * realPart) + std::abs(imaginaryMagnitudes[term] * imaginaryPart);
                imaginaryTerms +=
                    std::abs(realMagnitudes[term] * imaginaryPart) + std::abs(imaginaryMagnitudes[term] * realPart);
            }
            const double real = results[v].real(), imaginary = results[v].imag();
            results[v] = {real * double(std::abs(real) > rounding * realTerms),
                          imaginary * double(std::abs(imaginary) > rounding * imaginaryTerms)};
        }
#endif
    }

    // The same test for the two results of a pair without controls, out of
    // the pair loop. With vector instructions both are tested at once, the
    // first in the low half and the second in the high half: a term with
    // offset zero reads the two amplitudes in place, one with offset one
    // exchanged. Each amplitude is read as it was stored, in one piece.
    EIGEN_DONT_INLINE static void ClearRoundedPair(std::complex<double> *results, const std::complex<double> *pair,
                                                   const size_t *termOffsets, const double *realMagnitudes,
                                                   const double *imaginaryMagnitudes, size_t nrTerms,
                                                   double rounding) noexcept
    {
#ifdef QCSIM_AVX2
        const auto both = [](const std::complex<double> *values) {
            return _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadu_pd(reinterpret_cast<const double *>(values))),
                                        _mm_loadu_pd(reinterpret_cast<const double *>(values + 1)), 1);
        };
        const __m256d sign = _mm256_set1_pd(-0.0);
        const __m256d inPlace = _mm256_andnot_pd(sign, both(pair));
        const __m256d exchanged = _mm256_permute2f128_pd(inPlace, inPlace, 1);
        __m256d sums = _mm256_setzero_pd();
        for (size_t term = 0; term < nrTerms; ++term)
        {
            const __m256d parts = termOffsets[term] ? exchanged : inPlace;
            sums = _mm256_add_pd(
                sums, _mm256_add_pd(_mm256_andnot_pd(sign, _mm256_mul_pd(_mm256_set1_pd(realMagnitudes[term]), parts)),
                                    _mm256_andnot_pd(sign, _mm256_mul_pd(_mm256_set1_pd(imaginaryMagnitudes[term]),
                                                                         _mm256_permute_pd(parts, 5)))));
        }
        const __m256d result = both(results);
        const __m256d keep = _mm256_and_pd(
            _mm256_cmp_pd(_mm256_andnot_pd(sign, result), _mm256_mul_pd(_mm256_set1_pd(rounding), sums), _CMP_GT_OQ),
            _mm256_set1_pd(1.0));
        _mm256_storeu_pd(reinterpret_cast<double *>(results), _mm256_mul_pd(result, keep));
#else
        ClearRoundedParts(results, pair, nullptr, 2, termOffsets, realMagnitudes, imaginaryMagnitudes, 0, nrTerms,
                          rounding);
#endif
    }

    // The test out of line, for the general loops: a coset of any size keeps
    // its code out of them.
    EIGEN_DONT_INLINE static void ClearRoundedSums(std::complex<double> *results, const std::complex<double> *sources,
                                                   const std::complex<double> *rests, size_t nrOffsets,
                                                   const size_t *termOffsets, const double *realMagnitudes,
                                                   const double *imaginaryMagnitudes, size_t firstTerm, size_t nrTerms,
                                                   double rounding) noexcept
    {
        ClearRoundedParts(results, sources, rests, nrOffsets, termOffsets, realMagnitudes, imaginaryMagnitudes,
                          firstTerm, nrTerms, rounding);
    }

    // Apply a compiled matrix with visible controls to one frame in one pass
    // over the cosets of its labels, as ApplyPauliSum does for a flat sum.
    // With b the part of the state every control selects and V the compiled
    // sum, the result is w (a - b) + V b.
    //
    // A control's Z has the logical image P; the selected part is (b + s P b)/2
    // with s = -1 for a control value of one. Its two terms are exact (a
    // unit phase and a halving), so a state the control does not select gives
    // an exact zero and passes through bit for bit, and a selected one is
    // taken whole. No pair of rounded products has to cancel.
    //
    // The parts of a result are tested against the rounding of their own
    // sums as in the matrix variant of ApplyPauliSum, with the part that
    // passed the controls as one more term.
    void ApplyMatrixTerms(ExtendedFrame &frame, const MatrixAction &action, const size_t *qubits, size_t nrLocalQubits)
    {
        constexpr size_t maxTerms = detail::LocalPauliSum::Keys + 2;
        auto &w = sumWorkspace;
        const auto &basis = frame.cliffordBasis;
        const size_t nrWords = frame.signs.GetNrWords();
        const size_t nrControls = action.nrControls;
        const size_t nrTerms = nrControls + action.sum.Terms();
        const size_t generatorBase = nrTerms;
        EnsureRows(w.rows, nrTerms + std::min(nrTerms, PauliSumWorkspace::MaxGenerators), frame.GetNrQubits());

        // The controls' Z first, then the terms of the sum.
        std::array<detail::LocalPauli, maxTerms> paulis;
        ScratchAmplitudes<maxTerms> coefficients, negatedCoefficients;
        for (size_t control = 0; control < nrControls; ++control)
        {
            paulis[control] = detail::LocalPauli::Single(action.controlQubits[control], 'Z');
            coefficients[control] = {1.0, 0.0};
        }
        {
            size_t count = nrControls;
            action.sum.ForEachTerm(
                [&paulis, &coefficients, &count](detail::LocalPauli pauli, std::complex<double> value) {
                    paulis[count] = pauli;
                    coefficients[count++] = value;
                });
        }

        // Logical images of the terms, with i^|x&z| and their signs folded
        // into the coefficients.
        w.terms.clear();
        for (size_t term = 0; term < nrTerms; ++term)
        {
            auto image = w.rows[term];
            image.Clear();
            for (size_t qubit = 0; qubit < nrLocalQubits; ++qubit)
            {
                const char factor = paulis[term].Factor(qubit);
                if (factor != 'I')
                    basis.MultiplyImageInto(image, qubits[qubit], factor);
            }
            unsigned nrY = 0;
            bool hasZ = false;
            for (size_t word = 0; word < nrWords; ++word)
            {
                nrY += Clifford::detail::Popcount(image.X.words[word] & image.Z.words[word]);
                hasZ |= image.Z.words[word] != 0;
            }
            const auto phase = detail::LocalPauliSum::IPower(nrY);
            coefficients[term] = (image.PhaseSign ? -coefficients[term] : coefficients[term]) * phase;
            w.terms.push_back({coefficients[term], 0, hasZ});
        }

        // Generators of the X parts in reduced echelon form, as in ApplyPauliSum.
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
                    for (size_t word = 0; word < nrWords; ++word)
                        reduced.X.words[word] ^= other.X.words[word];
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
                if (!other.X[pivot])
                    continue;
                for (size_t word = 0; word < nrWords; ++word)
                    other.X.words[word] ^= reduced.X.words[word];
                combinations[g] ^= combination;
            }
            pivots[nrGenerators] = pivot;
            combinations[nrGenerators] = combination;
            generatorTerms[nrGenerators] = term;
            w.terms[term].offset = size_t(1) << nrGenerators;
            ++nrGenerators;
        }

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
                    for (size_t word = 0; word < nrWords; ++word)
                        offset.X.words[word] ^= generator.X.words[word];
                }
        }
        w.offsetParity.assign(nrTerms * nrOffsets, 0);
        for (size_t term = 0; term < nrTerms; ++term)
            if (w.terms[term].hasZ)
                for (size_t u = 1; u < nrOffsets; ++u)
                    w.offsetParity[term * nrOffsets + u] =
                        OddParity(w.rows[term].Z.words, w.offsets[u].X.words, nrWords);

        const double unit = 0.5 * std::numeric_limits<double>::epsilon();
        const size_t nrSumTerms = nrTerms - nrControls;
        const double weight = action.inactiveWeight;
        const size_t originalSize = frame.GetFrameSize();
        if (nrGenerators == 0)
        {
            // Every image is diagonal: a control keeps or passes a whole
            // component, and the sum rescales the selected ones in place.
            // The factor is a sum of exact coefficients: n - 1 additions.
            double realBound = 0.0, imaginaryBound = 0.0;
            for (size_t term = nrControls; term < nrTerms; ++term)
            {
                realBound += std::abs(coefficients[term].real());
                imaginaryBound += std::abs(coefficients[term].imag());
            }
            const double rounding = nrSumTerms > 1 ? double(nrSumTerms - 1) * unit : 0.0;
            realBound *= rounding;
            imaginaryBound *= rounding;
            for (size_t component = 0; component < originalSize; ++component)
            {
                const auto *label = frame.signs.LabelWords(component);
                bool selected = true;
                for (size_t control = 0; control < nrControls && selected; ++control)
                {
                    // The image is (+/-) Z^z: the control qubit is one where it is -1.
                    const bool one =
                        (coefficients[control].real() < 0.0) != OddParity(w.rows[control].Z.words, label, nrWords);
                    selected = one == action.controlValues[control];
                }
                if (!selected)
                {
                    if (weight != 1.0)
                        frame.amplitudes[component] *= weight;
                    continue;
                }
                double real = 0.0, imaginary = 0.0;
                for (size_t term = nrControls; term < nrTerms; ++term)
                {
                    const auto &t = w.terms[term];
                    const bool odd = t.hasZ && OddParity(w.rows[term].Z.words, label, nrWords);
                    real += odd ? -t.coefficient.real() : t.coefficient.real();
                    imaginary += odd ? -t.coefficient.imag() : t.coefficient.imag();
                }
                if (std::abs(real) <= realBound)
                    real = 0.0;
                if (std::abs(imaginary) <= imaginaryBound)
                    imaginary = 0.0;
                frame.amplitudes[component] *= std::complex<double>(real, imaginary);
            }
            return;
        }

        // The coset loop reads per-term data from local arrays.
        std::array<double, maxTerms> realMagnitudes, imaginaryMagnitudes;
        std::array<size_t, maxTerms> termOffsets;
        double coefficientBound = 0.0;
        for (size_t term = 0; term < nrTerms; ++term)
        {
            negatedCoefficients[term] = -coefficients[term];
            realMagnitudes[term] = std::abs(coefficients[term].real());
            imaginaryMagnitudes[term] = std::abs(coefficients[term].imag());
            termOffsets[term] = w.terms[term].offset;
            if (term >= nrControls)
                coefficientBound += realMagnitudes[term] + imaginaryMagnitudes[term];
        }
        // A control enters as (b + s P b)/2, with s = -1 where it selects one.
        for (size_t control = 0; control < nrControls; ++control)
            if (action.controlValues[control])
                std::swap(coefficients[control], negatedCoefficients[control]);
        // A result sums the n products and the part that passed the controls.
        const double rounding = double(nrSumTerms + 2) * unit;
        // No part of a result above this floor can be within the rounding of
        // its sum. Twice the largest amplitude component of the frame bounds
        // the L1 norm of every amplitude and of every selected part (an
        // average of two amplitudes up to unit phases), hence the products;
        // the part that passed the controls is a difference of two of them.
        const double largest = LargestAmplitude(frame);
        const double floor = rounding * (coefficientBound + 2.0 * weight) * largest;
        // Each product of the sum takes its operands from tables, as in the
        // matrix variant of ApplyPauliSum.
        ScratchAmplitudes<2 * maxTerms> signedCoefficients;
        w.readSources.resize(nrOffsets * nrSumTerms);
        w.readSigns.resize(nrOffsets * nrSumTerms);
        for (size_t term = nrControls; term < nrTerms; ++term)
        {
            signedCoefficients[2 * term] = coefficients[term];
            signedCoefficients[2 * term + 1] = negatedCoefficients[term];
        }
        for (size_t v = 0, read = 0; v < nrOffsets; ++v)
            for (size_t term = nrControls; term < nrTerms; ++term, ++read)
            {
                const size_t u = v ^ termOffsets[term];
                w.readSources[read] = uint8_t(u);
                w.readSigns[read] = uint8_t(2 * term + w.offsetParity[term * nrOffsets + u]);
            }
        const uint8_t *const readSources = w.readSources.data(), *const readSigns = w.readSigns.data();

        std::array<size_t, PauliSumWorkspace::MaxOffsets> members;
        ScratchAmplitudes<PauliSumWorkspace::MaxOffsets> amplitudes, results, rests;
        std::array<ScratchAmplitudes<PauliSumWorkspace::MaxOffsets>, 2> selections;
        std::array<bool, maxTerms> labelParity;

        const size_t notFound = PackedComponentIndex::NotFound;
        const std::complex<double> none(0.0, 0.0);
        frame.EnsureComponentIndex(2 * originalSize);
        w.visited.assign(originalSize, 0);
        w.label.resize(nrWords);
        bool appended = false;
        for (size_t component = 0; component < originalSize; ++component)
        {
            if (w.visited[component])
                continue;
            // The coset is the labels L xor offset; L is copied because
            // appending may reallocate the label storage.
            std::copy_n(frame.signs.LabelWords(component), nrWords, w.label.data());
            for (size_t u = 0; u < nrOffsets; ++u)
            {
                const size_t member = u == 0 ? component : frame.FindXorComponent(w.label.data(), w.offsets[u].X.words);
                members[u] = member;
                if (member == notFound)
                    amplitudes.Set(u, none);
                else
                {
                    amplitudes.Set(u, frame.amplitudes[member]);
                    w.visited[member] = 1;
                }
            }
            for (size_t term = 0; term < nrTerms; ++term)
                labelParity[term] = w.terms[term].hasZ && OddParity(w.rows[term].Z.words, w.label.data(), nrWords);

            // The part every control selects.
            const std::complex<double> *selected = amplitudes.data();
            for (size_t control = 0; control < nrControls; ++control)
            {
                auto &next = selections[control];
                const size_t offset = termOffsets[control];
                const uint8_t *parity = &w.offsetParity[control * nrOffsets];
                if (offset == 0)
                {
                    // Diagonal on this coset: each label is kept or dropped whole.
                    const bool positive = coefficients[control].real() > 0.0;
                    for (size_t v = 0; v < nrOffsets; ++v)
                        next.Set(v, positive != (labelParity[control] != bool(parity[v])) ? selected[v] : none);
                }
                else
                    for (size_t v = 0; v < nrOffsets; ++v)
                    {
                        const size_t u = v ^ offset;
                        const bool odd = labelParity[control] != bool(parity[u]);
                        next.Set(v, 0.5 * (selected[v] +
                                           (odd ? negatedCoefficients[control] : coefficients[control]) * selected[u]));
                    }
                selected = next.data();
            }

            // out[v] = w (a[v] - b[v]) + sum_j c_j (-1)^(z_j . (L xor offset_u)) b[u], u = v xor offset_j.
            RoundingScreen screen(floor);
            for (size_t v = 0, read = 0; v < nrOffsets; ++v)
            {
                std::complex<double> result(0.0, 0.0);
                for (size_t term = nrControls; term < nrTerms; ++term, ++read)
                    result +=
                        signedCoefficients[readSigns[read] ^ unsigned(labelParity[term])] * selected[readSources[read]];
                std::complex<double> rest = amplitudes[v] - selected[v];
                if (weight != 1.0)
                    rest *= weight;
                result += rest;
                rests.Set(v, rest);
                results.Set(v, result);
                screen.Add(result);
            }
            if (screen.Found())
                ClearRoundedSums(results.data(), selected, rests.data(), nrOffsets, termOffsets.data(),
                                 realMagnitudes.data(), imaginaryMagnitudes.data(), nrControls, nrTerms, rounding);

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
        if (appended)
            frame.InvalidateComponentIndex();
    }

    static size_t QubitMask(const size_t *qubits, unsigned arity) noexcept
    {
        size_t mask = 0;
        for (unsigned i = 0; i < arity; ++i)
            mask |= size_t{1} << (qubits[i] % std::numeric_limits<size_t>::digits);
        return mask;
    }

    static bool QubitsOverlap(const size_t *a, unsigned aSize, const size_t *b, unsigned bSize) noexcept
    {
        for (unsigned i = 0; i < aSize; ++i)
            for (unsigned j = 0; j < bSize; ++j)
                if (a[i] == b[j])
                    return true;
        return false;
    }

    static void CanonicalizeUnitaryQubits(size_t *qubits, const Gates::UnitaryGateDescriptor &descriptor) noexcept
    {
        // CC gates exchange their controls; Fredkin exchanges its targets.
        if (descriptor.controls == 2 && qubits[1] > qubits[2])
            std::swap(qubits[1], qubits[2]);
        else if (descriptor.kind == Gates::UnitaryGateDescriptor::Kind::Fredkin && qubits[0] > qubits[1])
            std::swap(qubits[0], qubits[1]);
    }

    // The mask is nonzero while entries are retained, and also while an
    // unused anchor still waits to be released, so the common case of no
    // history at all costs one test.
    void ClearUnitaryHistory() noexcept
    {
        if (unitaryHistory.qubitMask)
            ResetUnitaryHistory();
    }

    // The qubits the history covers (all of them while only its anchor's
    // storage is left to release), and with them the limit below which a
    // native gate does not look at the history.
    void SetHistoryQubits(size_t mask) noexcept
    {
        unitaryHistory.qubitMask = mask;
        unguardedQubits = mask ? 0 : GetNrQubits();
    }

    // Drop the entries and the later states kept for matching inverses. The
    // anchor's storage outlives a few more state-changing operations, since
    // the next matrix gate usually reopens the history; a large one that
    // stays unused is then released instead of doubling the memory of a
    // long native computation.
    EIGEN_DONT_INLINE void ResetUnitaryHistory() noexcept
    {
        auto &history = unitaryHistory;
        amplitudeBound = -1.0;
        if (history.size != 0)
        {
            history.size = 0;
            history.idle = 0;
            std::vector<UnitaryHistory::Checkpoint>().swap(history.checkpoints);
            std::vector<FrameComponents>().swap(history.held);
            size_t bytes = 0;
            for (const auto &frame : frames)
                bytes += frame.nextAmplitudes.capacity() *
                         (sizeof(std::complex<double>) + frame.signs.GetNrWords() * sizeof(ExtendedFrame::Word));
            SetHistoryQubits(bytes >= UnitaryHistory::ReleasableBytes ? ~size_t(0) : 0);
            return;
        }
        if (++history.idle < UnitaryHistory::IdleOperations)
            return;
        for (auto &frame : frames)
        {
            std::vector<std::complex<double>>().swap(frame.nextAmplitudes);
            PackedComponentLabels(frame.GetNrQubits(), 0).swap(frame.nextSigns);
        }
        SetHistoryQubits(0);
    }

    static bool SameUnitaryQubits(const size_t *a, const size_t *b, unsigned arity,
                                  const Gates::UnitaryGateDescriptor &descriptor) noexcept
    {
        if (descriptor.kind == Gates::UnitaryGateDescriptor::Kind::Fredkin)
            return a[2] == b[2] && ((a[0] == b[0] && a[1] == b[1]) || (a[0] == b[1] && a[1] == b[0]));
        if (a[0] != b[0])
            return false;
        if (descriptor.controls == 2)
            return (a[1] == b[1] && a[2] == b[2]) || (a[1] == b[2] && a[2] == b[1]);
        return (arity < 2 || a[1] == b[1]) && (arity < 3 || a[2] == b[2]);
    }

    template <class Matrix>
    bool RestoreAdjacentUnitaryInverse(const Matrix &matrix, const size_t *qubits, unsigned arity)
    {
        const auto &entry = unitaryHistory.entries[0];
        if (arity != entry.arity)
            return false;
        if (!SameUnitaryQubits(qubits, entry.qubits.data(), arity, entry.descriptor) ||
            matrix.rows() != unitaryHistory.firstAdjoint.rows() ||
            matrix.cols() != unitaryHistory.firstAdjoint.cols() ||
            !(matrix.array() == unitaryHistory.firstAdjoint.array()).all())
            return false;
        RestoreUnitaryAnchor();
        return true;
    }

    void RestoreUnitaryAnchor() noexcept
    {
        for (auto &frame : frames)
        {
            frame.amplitudes.swap(frame.nextAmplitudes);
            frame.signs.swap(frame.nextSigns);
            frame.InvalidateComponentIndex();
        }
        approximationStatistics = unitaryHistory.statistics;
        ClearUnitaryHistory();
    }

    // Called by every native gate. Without a history the mask is zero, and
    // this test is all a native gate pays for the history.
    void InvalidateOverlappingUnitaryHistory(const size_t *qubits, unsigned arity) noexcept
    {
        if (unitaryHistory.qubitMask)
            InvalidateUnitaryHistoryOn(qubits, arity);
    }

    EIGEN_DONT_INLINE void InvalidateUnitaryHistoryOn(const size_t *qubits, unsigned arity) noexcept
    {
        if (!(unitaryHistory.qubitMask & QubitMask(qubits, arity)))
            return;
        if (unitaryHistory.size == 0)
        {
            // Only an unused anchor is left: count this operation against it.
            ResetUnitaryHistory();
            return;
        }
        for (size_t i = 0; i < unitaryHistory.size; ++i)
        {
            const auto &entry = unitaryHistory.entries[i];
            if (QubitsOverlap(qubits, arity, entry.qubits.data(), entry.arity))
            {
                ResetUnitaryHistory();
                return;
            }
        }
    }

    // Rebuild the live components from the anchor and every retained entry.
    // Entries own their compiled actions, so nothing is rebuilt or looked up.
    void ReplayUnitaryHistory()
    {
        amplitudeBound = -1.0;
        for (auto &frame : frames)
        {
            frame.amplitudes = frame.nextAmplitudes;
            frame.signs = frame.nextSigns;
            frame.InvalidateComponentIndex();
        }
        approximationStatistics = unitaryHistory.statistics;
        for (size_t i = 0; i < unitaryHistory.size; ++i)
        {
            const auto &entry = unitaryHistory.entries[i];
            ApplyMatrixSum(*entry.action, entry.qubits.data(), entry.arity);
        }
    }

    template <class Matrix>
    bool RestoreUnitaryInverse(const Matrix &matrix, const Gates::UnitaryGateDescriptor &descriptor,
                               const MatrixAction &action, const size_t *qubits, unsigned arity)
    {
        if (unitaryHistory.size == 0 || !(unitaryHistory.qubitMask & QubitMask(qubits, arity)))
            return false;
        for (size_t i = unitaryHistory.size; i-- > 0;)
        {
            const auto &entry = unitaryHistory.entries[i];
            if (!QubitsOverlap(qubits, arity, entry.qubits.data(), entry.arity))
                continue;
            // A different overlapping gate is a barrier, even if an earlier
            // entry has the inverse descriptor. Disjoint later gates commute.
            if (arity != entry.arity)
                return false;
            const bool semanticInverse =
                action.unitary && Gates::AreUnitaryGateDescriptorsInverse(entry.descriptor, descriptor);
            if (!semanticInverse && entry.action->inverseFingerprint != action.matrixFingerprint)
                return false;
            if (!SameUnitaryQubits(qubits, entry.qubits.data(), arity, entry.descriptor))
                return false;
            if (!semanticInverse && !MatchesRawUnitaryInverse(matrix, entry.descriptor))
                return false;
            if (unitaryHistory.size == 1)
            {
                RestoreUnitaryAnchor();
                return true;
            }

            RestoreMatchedUnitaryInverse(i);
            return true;
        }
        return false;
    }

    template <class Matrix>
    EIGEN_DONT_INLINE static bool MatchesRawUnitaryInverse(const Matrix &matrix,
                                                           const Gates::UnitaryGateDescriptor &descriptor)
    {
        const auto original = Gates::UnitaryGateMatrix(descriptor);
        return (matrix.array() == original.adjoint().array()).all();
    }

    // Remove one retained entry: the live components become those of the
    // anchor followed by every other entry in order, bit for bit. The replay
    // starts from the nearest state kept at or before the removed entry, and
    // leaves copies of a few of the states it passes (at distances 1, 2, 4,
    // ... below the removed entry), so that undoing a nest of gates in
    // reverse order replays about one gate per inverse instead of the whole
    // history. Each kept state is a plain copy of what the replay computed:
    // exactness does not depend on which states are kept. Their number is
    // bounded, and so is their total size.
    EIGEN_DONT_INLINE void RestoreMatchedUnitaryInverse(size_t removed)
    {
        auto &history = unitaryHistory;
        auto &checkpoints = history.checkpoints;
        // States after the removed entry include it.
        while (!checkpoints.empty() && checkpoints.back().position > removed)
            checkpoints.pop_back();
        const bool last = removed + 1 == history.size;
        const size_t start = checkpoints.empty() ? 0 : checkpoints.back().position;
        amplitudeBound = -1.0;

        // States worth keeping on the way, nearest first, within the budget.
        size_t bytes = 0;
        for (const auto &frame : frames)
            bytes += frame.amplitudes.size() *
                     (sizeof(std::complex<double>) + frame.signs.GetNrWords() * sizeof(ExtendedFrame::Word));
        const size_t allowed =
            std::min(UnitaryHistory::MaxCheckpoints, UnitaryHistory::CheckpointBytes / std::max<size_t>(bytes, 1));
        std::array<size_t, UnitaryHistory::MaxCheckpoints> wanted{};
        size_t nrWanted = 0;
        for (size_t distance = 1; checkpoints.size() + nrWanted < allowed && distance < removed - start; distance *= 2)
            wanted[nrWanted++] = removed - distance;

        // Park the live components: they return untouched if the replay
        // fails, and their storage is reused by the next replay otherwise.
        auto &held = history.held;
        if (held.size() != frames.size())
        {
            held.clear();
            for (const auto &frame : frames)
                held.push_back({{}, PackedComponentLabels(frame.GetNrQubits(), 0)});
        }
        for (size_t f = 0; f < frames.size(); ++f)
        {
            held[f].amplitudes.swap(frames[f].amplitudes);
            held[f].signs.swap(frames[f].signs);
        }
        const auto statistics = approximationStatistics;
        try
        {
            if (start == removed && last)
                // The wanted state itself was kept: take it.
                for (size_t f = 0; f < frames.size(); ++f)
                {
                    frames[f].amplitudes.swap(checkpoints.back().frames[f].amplitudes);
                    frames[f].signs.swap(checkpoints.back().frames[f].signs);
                }
            else
                for (size_t f = 0; f < frames.size(); ++f)
                {
                    frames[f].amplitudes = start ? checkpoints.back().frames[f].amplitudes : frames[f].nextAmplitudes;
                    frames[f].signs = start ? checkpoints.back().frames[f].signs : frames[f].nextSigns;
                }
            for (auto &frame : frames)
                frame.InvalidateComponentIndex();
            approximationStatistics = history.statistics;
            for (size_t i = start; i < history.size; ++i)
            {
                if (i == removed)
                    continue;
                if (i > start && std::find(wanted.begin(), wanted.begin() + nrWanted, i) != wanted.begin() + nrWanted)
                {
                    UnitaryHistory::Checkpoint checkpoint;
                    checkpoint.position = i;
                    checkpoint.frames.reserve(frames.size());
                    for (const auto &frame : frames)
                        checkpoint.frames.push_back({frame.amplitudes, frame.signs});
                    checkpoints.push_back(std::move(checkpoint));
                }
                const auto &entry = history.entries[i];
                ApplyMatrixSum(*entry.action, entry.qubits.data(), entry.arity);
            }
        }
        catch (...)
        {
            // Nothing is applied after taking a kept state whole, so the
            // parked components always return to an unchanged history.
            for (size_t f = 0; f < frames.size(); ++f)
            {
                frames[f].amplitudes.swap(held[f].amplitudes);
                frames[f].signs.swap(held[f].signs);
                frames[f].InvalidateComponentIndex();
            }
            approximationStatistics = statistics;
            throw;
        }
        for (size_t j = removed + 1; j < history.size; ++j)
            history.entries[j - 1] = std::move(history.entries[j]);
        --history.size;
        // A kept state at the new end is the live state itself.
        while (!checkpoints.empty() && checkpoints.back().position >= history.size)
            checkpoints.pop_back();
        if (history.size == 1)
            history.firstAdjoint = Gates::UnitaryGateMatrix(history.entries[0].descriptor).adjoint();
        size_t covered = 0;
        for (size_t j = 0; j < history.size; ++j)
        {
            const auto &remaining = history.entries[j];
            covered |= QubitMask(remaining.qubits.data(), remaining.arity);
        }
        SetHistoryQubits(covered);
    }

    struct PauliSumWorkspace
    {
        static constexpr size_t MaxOffsets = detail::LocalPauliSum::Keys;
        static constexpr size_t MaxGenerators = 2 * detail::LocalPauli::MaxQubits;

        struct Term
        {
            std::complex<double> coefficient; // includes the image's phase
            size_t offset;                    // coset offset index of its X part
            bool hasZ;
        };

        // Term images then reduced generators, and the coset offsets; both
        // grow to the largest sum seen.
        Clifford::detail::PackedTableau rows;
        Clifford::detail::PackedTableau offsets;
        std::vector<Term> terms;
        std::vector<uint8_t> offsetParity; // term-major parity(z . offset)
        // The reads of the matrix passes in the order they happen, result by
        // result and term by term: the coset member a term reads, and its
        // entry in the table of signed coefficients (twice the term, plus
        // one where the offset flips the sign).
        std::vector<uint8_t> readSources, readSigns;
        std::vector<uint8_t> visited;
        std::vector<ExtendedFrame::Word> label;
    };

    static void EnsureRows(Clifford::detail::PackedTableau &rows, size_t count, size_t qubits)
    {
        if (rows.GetNrQubits() == qubits && rows.size() >= count)
            return;
        rows = Clifford::detail::PackedTableau(rows.GetNrQubits() == qubits ? std::max(count, rows.size()) : count,
                                               qubits);
    }

    static bool OddParity(const ExtendedFrame::Word *left, const ExtendedFrame::Word *right, size_t nrWords) noexcept
    {
        ExtendedFrame::Word parity = 0;
        for (size_t word = 0; word < nrWords; ++word)
            parity ^= left[word] & right[word];
        return (Clifford::detail::Popcount(parity) & 1U) != 0;
    }

    static double L1(const std::complex<double> &value) noexcept
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
    //
    // MatrixPath is the variant for compiled matrices, whose coefficients
    // must not be dropped by a tolerance. A result is the sum of n rounded
    // products; its rounding error is at most (n + 1) u times the sum of the
    // magnitudes of the real products that enter each of its two parts
    // (u = epsilon/2: each product is the difference of two rounded real
    // products, 2u, and each of the n - 1 additions adds at most u of the
    // partial sum). A real or imaginary part within that bound of zero is the
    // rounding of its own cancellation and is cleared; a small part with no
    // cancellation behind it is far above its bound and stays, and the parts
    // are tested separately, so a genuine tiny imaginary part survives next
    // to a cancelled real one. The bound is only evaluated for a coset where
    // some part is below a floor no cancellation can exceed, so a pass in
    // which nothing is that small pays two comparisons per result. The
    // caller normalizes and prunes.
    template <bool MatrixPath = false>
    void ApplyPauliSum(ExtendedFrame &frame, const detail::LocalPauliSum &sum, const size_t *qubits,
                       size_t nrLocalQubits)
    {
        auto &w = sumWorkspace;
        const auto &basis = frame.cliffordBasis;
        const size_t nrWords = frame.signs.GetNrWords();
        const size_t nrTerms = sum.Terms();
        const size_t generatorBase = nrTerms;
        EnsureRows(w.rows, nrTerms + std::min(nrTerms, PauliSumWorkspace::MaxGenerators), frame.GetNrQubits());

        // Logical images of the terms, with i^|x&z| and their signs folded
        // into the coefficients.
        w.terms.clear();
        sum.ForEachTerm([&](detail::LocalPauli pauli, std::complex<double> coefficient) {
            auto image = w.rows[w.terms.size()];
            image.Clear();
            for (size_t qubit = 0; qubit < nrLocalQubits; ++qubit)
            {
                const char factor = pauli.Factor(qubit);
                if (factor != 'I')
                    basis.MultiplyImageInto(image, qubits[qubit], factor);
            }
            unsigned nrY = 0;
            bool hasZ = false;
            for (size_t word = 0; word < nrWords; ++word)
            {
                nrY += Clifford::detail::Popcount(image.X.words[word] & image.Z.words[word]);
                hasZ |= image.Z.words[word] != 0;
            }
            const auto phase = detail::LocalPauliSum::IPower(nrY);
            w.terms.push_back({(image.PhaseSign ? -coefficient : coefficient) * phase, 0, hasZ});
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
                    for (size_t word = 0; word < nrWords; ++word)
                        reduced.X.words[word] ^= other.X.words[word];
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
                if (!other.X[pivot])
                    continue;
                for (size_t word = 0; word < nrWords; ++word)
                    other.X.words[word] ^= reduced.X.words[word];
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
                    for (size_t word = 0; word < nrWords; ++word)
                        offset.X.words[word] ^= generator.X.words[word];
                }
        }
        w.offsetParity.assign(nrTerms * nrOffsets, 0);
        for (size_t term = 0; term < nrTerms; ++term)
            if (w.terms[term].hasZ)
                for (size_t u = 1; u < nrOffsets; ++u)
                    w.offsetParity[term * nrOffsets + u] =
                        OddParity(w.rows[term].Z.words, w.offsets[u].X.words, nrWords);

        const double tolerance = 8.0 * double(nrTerms + 1) * std::numeric_limits<double>::epsilon();
        const size_t originalSize = frame.GetFrameSize();
        // The matrix variant's rounding unit, and twice the largest amplitude
        // component of the frame: no term reads an amplitude with a larger
        // L1 norm, so one pass over the frame replaces a bound per coset.
        const double unit = 0.5 * std::numeric_limits<double>::epsilon();
        double largest = 0.0;
        if constexpr (MatrixPath)
            largest = LargestAmplitude(frame);
        if (nrGenerators == 0)
        {
            // A diagonal sum rescales every amplitude in place. In the
            // matrix variant the factor is a sum of exact coefficients,
            // n - 1 additions, tested part by part.
            double realFloor = 0.0, imaginaryFloor = 0.0;
            if constexpr (MatrixPath)
            {
                for (size_t term = 0; term < nrTerms; ++term)
                {
                    realFloor += std::abs(w.terms[term].coefficient.real());
                    imaginaryFloor += std::abs(w.terms[term].coefficient.imag());
                }
                const double additions = nrTerms > 1 ? double(nrTerms - 1) * unit : 0.0;
                realFloor *= additions;
                imaginaryFloor *= additions;
            }
            if constexpr (MatrixPath)
                if (nrTerms <= 4)
                {
                    // A few terms give a few factors: one per pattern of
                    // their signs, summed in the order of the general loop
                    // and tested once instead of once per component.
                    ScratchAmplitudes<16> factors;
                    for (size_t pattern = 0; pattern < (size_t(1) << nrTerms); ++pattern)
                    {
                        std::complex<double> factor(0.0, 0.0);
                        for (size_t term = 0; term < nrTerms; ++term)
                            factor += ((pattern >> term) & 1U) ? -w.terms[term].coefficient : w.terms[term].coefficient;
                        factors[pattern] = {std::abs(factor.real()) <= realFloor ? 0.0 : factor.real(),
                                            std::abs(factor.imag()) <= imaginaryFloor ? 0.0 : factor.imag()};
                    }
                    for (size_t component = 0; component < originalSize; ++component)
                    {
                        const auto *label = frame.signs.LabelWords(component);
                        size_t pattern = 0;
                        for (size_t term = 0; term < nrTerms; ++term)
                            pattern |= size_t(w.terms[term].hasZ && OddParity(w.rows[term].Z.words, label, nrWords))
                                       << term;
                        frame.amplitudes[component] *= factors[pattern];
                    }
                    return;
                }
            for (size_t component = 0; component < originalSize; ++component)
            {
                const auto *label = frame.signs.LabelWords(component);
                std::complex<double> factor(0.0, 0.0);
                double bound = 0.0;
                for (size_t term = 0; term < nrTerms; ++term)
                {
                    const auto &t = w.terms[term];
                    const bool odd = t.hasZ && OddParity(w.rows[term].Z.words, label, nrWords);
                    factor += odd ? -t.coefficient : t.coefficient;
                    if constexpr (!MatrixPath)
                        bound += L1(t.coefficient);
                }
                if constexpr (MatrixPath)
                    frame.amplitudes[component] *=
                        std::complex<double>(std::abs(factor.real()) <= realFloor ? 0.0 : factor.real(),
                                             std::abs(factor.imag()) <= imaginaryFloor ? 0.0 : factor.imag());
                else
                    frame.amplitudes[component] = L1(factor) <= tolerance * bound
                                                      ? std::complex<double>(0.0, 0.0)
                                                      : frame.amplitudes[component] * factor;
            }
            if constexpr (!MatrixPath)
                PruneComponents(frame);
            return;
        }

        // The coset loop reads per-term data from local arrays.
        constexpr size_t maxTerms = detail::LocalPauliSum::Keys;
        std::conditional_t<MatrixPath, ScratchAmplitudes<maxTerms>, std::array<std::complex<double>, maxTerms>>
            coefficients, negatedCoefficients;
        std::array<double, maxTerms> magnitudes, realMagnitudes, imaginaryMagnitudes;
        std::array<size_t, maxTerms> termOffsets;
        double floor = 0.0;
        for (size_t term = 0; term < nrTerms; ++term)
        {
            coefficients[term] = w.terms[term].coefficient;
            negatedCoefficients[term] = -coefficients[term];
            if constexpr (MatrixPath)
            {
                realMagnitudes[term] = std::abs(coefficients[term].real());
                imaginaryMagnitudes[term] = std::abs(coefficients[term].imag());
                floor += realMagnitudes[term] + imaginaryMagnitudes[term];
            }
            else
                magnitudes[term] = L1(coefficients[term]);
            termOffsets[term] = w.terms[term].offset;
        }
        // No part above this floor can be within the rounding of its sum:
        // the products that enter it are bounded by the coefficients' L1
        // norms times the largest amplitude.
        const double rounding = double(nrTerms + 1) * unit;
        floor *= rounding * largest;
        // The matrix variant takes each product's operands from tables: the
        // coefficient with its sign by index, with no branch on a parity
        // that changes from coset to coset.
        ScratchAmplitudes<2 * maxTerms> signedCoefficients;
        if constexpr (MatrixPath)
        {
            w.readSources.resize(nrOffsets * nrTerms);
            w.readSigns.resize(nrOffsets * nrTerms);
            for (size_t term = 0; term < nrTerms; ++term)
            {
                signedCoefficients[2 * term] = coefficients[term];
                signedCoefficients[2 * term + 1] = negatedCoefficients[term];
            }
            for (size_t v = 0, read = 0; v < nrOffsets; ++v)
                for (size_t term = 0; term < nrTerms; ++term, ++read)
                {
                    const size_t u = v ^ termOffsets[term];
                    w.readSources[read] = uint8_t(u);
                    w.readSigns[read] = uint8_t(2 * term + w.offsetParity[term * nrOffsets + u]);
                }
        }
        const uint8_t *const readSources = w.readSources.data(), *const readSigns = w.readSigns.data();
        std::array<size_t, PauliSumWorkspace::MaxOffsets> members;
        std::conditional_t<MatrixPath, ScratchAmplitudes<PauliSumWorkspace::MaxOffsets>,
                           std::array<std::complex<double>, PauliSumWorkspace::MaxOffsets>>
            amplitudes, results;
        std::array<double, PauliSumWorkspace::MaxOffsets> amplitudeMagnitudes;
        std::array<bool, maxTerms> labelParity;

        const size_t notFound = PackedComponentIndex::NotFound;
        frame.EnsureComponentIndex(2 * originalSize);
        w.visited.assign(originalSize, 0);
        w.label.resize(nrWords);
        bool appended = false;
        if constexpr (MatrixPath)
            if (nrOffsets == 2)
            {
                // One generator: the cosets are pairs, a label and its
                // partner, and the two results are summed side by side, in
                // the order of the general loop below.
                const ExtendedFrame::Word *const partnerOffset = w.offsets[1].X.words;
                const std::complex<double> none(0.0, 0.0);
                ScratchAmplitudes<2> pair, sums;
                bool anyZ = false;
                for (size_t term = 0; term < nrTerms; ++term)
                {
                    labelParity[term] = false;
                    anyZ |= w.terms[term].hasZ;
                }
                for (size_t component = 0; component < originalSize; ++component)
                {
                    if (w.visited[component])
                        continue;
                    std::copy_n(frame.signs.LabelWords(component), nrWords, w.label.data());
                    const size_t partner = frame.FindXorComponent(w.label.data(), partnerOffset);
                    pair.Set(0, frame.amplitudes[component]);
                    pair.Set(1, partner == notFound ? none : frame.amplitudes[partner]);
                    if (partner != notFound)
                        w.visited[partner] = 1;
                    if (anyZ)
                        for (size_t term = 0; term < nrTerms; ++term)
                            labelParity[term] =
                                w.terms[term].hasZ && OddParity(w.rows[term].Z.words, w.label.data(), nrWords);
                    std::complex<double> first(0.0, 0.0), second(0.0, 0.0);
                    for (size_t term = 0; term < nrTerms; ++term)
                    {
                        const unsigned parity = labelParity[term];
                        first += signedCoefficients[readSigns[term] ^ parity] * pair[readSources[term]];
                        second +=
                            signedCoefficients[readSigns[nrTerms + term] ^ parity] * pair[readSources[nrTerms + term]];
                    }
                    sums.Set(0, first);
                    sums.Set(1, second);
                    RoundingScreen screen(floor);
                    screen.Add(first);
                    screen.Add(second);
                    if (screen.Found())
                        ClearRoundedPair(sums.data(), pair.data(), termOffsets.data(), realMagnitudes.data(),
                                         imaginaryMagnitudes.data(), nrTerms, rounding);
                    frame.amplitudes[component] = sums[0];
                    if (partner != notFound)
                        frame.amplitudes[partner] = sums[1];
                    else if (sums[1] != none)
                    {
                        frame.signs.AppendXor(w.label.data(), partnerOffset);
                        frame.amplitudes.push_back(sums[1]);
                        appended = true;
                    }
                }
                if (appended)
                    frame.InvalidateComponentIndex();
                return;
            }
        for (size_t component = 0; component < originalSize; ++component)
        {
            if (w.visited[component])
                continue;
            // The coset is the labels L xor offset; L is copied because
            // appending may reallocate the label storage.
            std::copy_n(frame.signs.LabelWords(component), nrWords, w.label.data());
            for (size_t u = 0; u < nrOffsets; ++u)
            {
                const size_t member = u == 0 ? component : frame.FindXorComponent(w.label.data(), w.offsets[u].X.words);
                members[u] = member;
                if (member == notFound)
                    amplitudes[u] = 0.0;
                else
                {
                    amplitudes[u] = frame.amplitudes[member];
                    w.visited[member] = 1;
                }
                if constexpr (!MatrixPath)
                    amplitudeMagnitudes[u] = L1(amplitudes[u]);
            }

            // out[v] = sum_j c_j (-1)^(z_j . (L xor offset_u)) in[u], u = v xor offset_j.
            // The bound sums the contributions' magnitudes for the rounding test.
            for (size_t term = 0; term < nrTerms; ++term)
                labelParity[term] = w.terms[term].hasZ && OddParity(w.rows[term].Z.words, w.label.data(), nrWords);
            RoundingScreen screen(floor);
            for (size_t v = 0, read = 0; v < nrOffsets; ++v)
            {
                std::complex<double> result(0.0, 0.0);
                if constexpr (MatrixPath)
                {
                    for (size_t term = 0; term < nrTerms; ++term, ++read)
                        result += signedCoefficients[readSigns[read] ^ unsigned(labelParity[term])] *
                                  amplitudes[readSources[read]];
                    results.Set(v, result);
                    screen.Add(result);
                }
                else
                {
                    double bound = 0.0;
                    for (size_t term = 0; term < nrTerms; ++term)
                    {
                        const size_t u = v ^ termOffsets[term];
                        const bool odd = labelParity[term] != bool(w.offsetParity[term * nrOffsets + u]);
                        result += (odd ? negatedCoefficients[term] : coefficients[term]) * amplitudes[u];
                        bound += magnitudes[term] * amplitudeMagnitudes[u];
                    }
                    results[v] = L1(result) <= tolerance * bound ? std::complex<double>(0.0, 0.0) : result;
                }
            }
            if constexpr (MatrixPath)
                if (screen.Found())
                    ClearRoundedSums(results.data(), amplitudes.data(), nullptr, nrOffsets, termOffsets.data(),
                                     realMagnitudes.data(), imaginaryMagnitudes.data(), 0, nrTerms, rounding);

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
        if (appended)
            frame.InvalidateComponentIndex();
        if constexpr (!MatrixPath)
            PruneComponents(frame);
    }

    static double PauliExpectation(const ExtendedFrame &frame, const PauliAction &action)
    {
        if (!HasPauliFlip(action))
        {
            double expectation = 0.0;
            for (size_t component = 0; component < frame.GetFrameSize(); ++component)
                expectation += std::norm(frame.amplitudes[component]) *
                               PauliPhase(frame.signs.LabelWords(component), action).real();
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
            const auto *sourceLabel = frame.signs.LabelWords(component);
            const size_t target = frame.FindXorComponent(sourceLabel, action.flipMask);
            if (target == notFound)
                continue;

            expectation +=
                std::conj(frame.amplitudes[target]) * PauliPhase(sourceLabel, action) * frame.amplitudes[component];
        }

        return expectation.real();
    }

    double PauliExpectation(const ExtendedFrame &frame, const CliffordBasisMap::Row &pauli) const
    {
        CompilePauliAction(frame, pauli, pauliActionWorkspace);
        return PauliExpectation(frame, pauliActionWorkspace);
    }

    bool MeasureSingleStabilizer(size_t qubit, const bool *forcedOutcome = nullptr)
    {
        auto &frame = frames.front();
        const auto observable = frame.cliffordBasis.ImageZ(qubit);
        const auto *label = frame.signs.LabelWords(0);
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
            if (labelParity != 0)
                negate = !negate;

            std::complex<double> phase;
            switch (nrY % 4)
            {
            case 0:
                phase = {1.0, 0.0};
                break;
            case 1:
                phase = {0.0, 1.0};
                break;
            case 2:
                phase = {-1.0, 0.0};
                break;
            default:
                phase = {0.0, -1.0};
                break;
            }
            frame.amplitudes.front() *= negate ? -phase : phase;
        }
        frame.cliffordBasis.MeasureZ(qubit, pivot, outcome, frame.signs.LabelWords(0), enableMultithreading);
        frame.InvalidateComponentIndex();
        return outcome;
    }

    void ApplyAxisRotation(size_t physicalQubit, double angle, char axis)
    {
        if (!std::isfinite(angle))
            throw std::invalid_argument("Rotation angle must be finite");
        if (angle == 0.0)
            return;
        ClearUnitaryHistory();
        long long quarterTurns = 0;
        if (detail::TryGetQuarterTurns(angle, quarterTurns))
        {
            ApplyCliffordRotation(physicalQubit, axis, quarterTurns);
            return;
        }

        const double halfAngle = 0.5 * angle;
        auto rotation = detail::LocalPauliSum::Identity(std::cos(halfAngle));
        rotation.Add(detail::LocalPauli::Single(0, axis), {0.0, -std::sin(halfAngle)});
        const size_t qubits[] = {physicalQubit};
        for (auto &frame : frames)
            ApplyPauliSum(frame, rotation, qubits, 1);
    }

    std::vector<ExtendedFrame> frames;
    // A native gate on qubits below this index does not look at the unitary
    // history: the register size while there is none, zero while there is
    // one. It sits beside the frames every gate reads.
    size_t unguardedQubits = 0;
    std::vector<ExtendedFrame> savedFrames;
    ExtendedStabilizerApproximationPolicy approximationPolicy;
    ExtendedStabilizerApproximationPolicy savedApproximationPolicy;
    ExtendedStabilizerApproximationStatistics approximationStatistics;
    ExtendedStabilizerApproximationStatistics savedApproximationStatistics;
    mutable PauliAction pauliActionWorkspace;
    PauliSumWorkspace sumWorkspace;
    detail::LocalGate localGate{1};
    // Positive while every frame is known to have unit norm: an upper bound
    // of the L1 norm of each amplitude. A remembered matrix gate sets it; it
    // lasts as long as the unitary history, which every other operation that
    // changes an amplitude clears.
    double amplitudeBound = -1.0;

    struct UnitaryHistory
    {
        static constexpr size_t Capacity = 32;
        // Later states a replay may keep besides the anchor, and their total size.
        static constexpr size_t MaxCheckpoints = 5;
        static constexpr size_t CheckpointBytes = size_t(64) << 20;
        // An unused anchor of at least this size is released after this many
        // further state-changing operations.
        static constexpr size_t ReleasableBytes = size_t(1) << 20;
        static constexpr unsigned IdleOperations = 16;

        struct Entry
        {
            Gates::UnitaryGateDescriptor descriptor;
            // Compiled once and shared with the cache: a replay looks nothing up.
            std::shared_ptr<const MatrixAction> action;
            std::array<size_t, 3> qubits{};
            unsigned arity = 0;
        };

        // The components after the first `position` entries.
        struct Checkpoint
        {
            size_t position = 0;
            std::vector<FrameComponents> frames;
        };

        size_t size = 0;
        size_t qubitMask = 0;
        unsigned idle = 0;
        ExtendedStabilizerApproximationStatistics statistics;
        Gates::UnitaryLocalMatrix firstAdjoint;
        // Allocated when the first matrix gate is remembered: a clone made
        // only to sample or to evaluate a probability carries no history.
        std::vector<Entry> entries;
        std::vector<Checkpoint> checkpoints; // ascending positions inside (0, size)
        std::vector<FrameComponents> held;
    } unitaryHistory;

    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist;
    bool enableMultithreading = true;
};

} // namespace QC
