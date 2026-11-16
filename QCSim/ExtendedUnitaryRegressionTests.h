#pragma once

#include "ExtendedStabilizer.h"
#include "QuantumGate.h"
#include <array>
#include <functional>
#include <iostream>
#include <utility>

namespace QC
{
namespace ExtendedUnitaryRegression
{
inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline void Prepare(ExtendedStabilizer &state, bool tiny)
{
    state.ApplyH(0);
    state.ApplyH(1);
    state.ApplyH(2);
    if (tiny)
    {
        state.ApplyRy(0, 1e-200);
        state.ApplyRy(1, .317);
        state.ApplyCX(2, 1);
    }
}

template <class Gate> inline void RoundTrip(const Gate &gate, bool tiny)
{
    ExtendedStabilizer state(3);
    Prepare(state, tiny);
    const auto amplitudes = state.GetFrames().front().amplitudes;
    const auto signs = state.GetFrames().front().signs;
    state.ApplyGate(gate, 2, 0, 1);
    state.SaveState(); // A public checkpoint must remain independent.
    const auto savedAmplitudes = state.GetFrames().front().amplitudes;
    const Gates::AppliedGate<> inverse(gate.getRawOperatorMatrix().adjoint().eval(), 2, 0, 1);
    state.ApplyGate(inverse);
    Require(state.GetFrames().front().amplitudes == amplitudes && state.GetFrames().front().signs == signs,
            "A certified controlled inverse pair failed to restore its exact input components");
    state.RestoreState();
    Require(state.GetFrames().front().amplitudes == savedAmplitudes,
            "Controlled inverse restoration overwrote the user's saved state");
    // Restore/clone must not inherit a stale internal undo checkpoint.
    auto clone = state.Clone();
    clone->ApplyGate(inverse);
    if (!tiny)
        for (unsigned index = 0; index < 8; ++index)
            Require(std::abs(clone->getBasisStateProbability(index) - .125) < 2e-13,
                    "Cloning a controlled gate state retained a stale internal checkpoint");
}

inline void Compare(const ExtendedStabilizer &actual, const ExtendedStabilizer &reference)
{
    const size_t qubits = actual.GetNrQubits();
    Require(qubits == reference.GetNrQubits(), "History reference has a different register size");
    for (size_t code = 0; code < (size_t(1) << (2 * qubits)); ++code)
    {
        std::string observable(qubits, 'I');
        for (size_t q = 0; q < qubits; ++q)
            observable[q] = "IXYZ"[(code >> (2 * q)) & 3U];
        Require(std::abs(actual.ExpectationValue(observable) - reference.ExpectationValue(observable)) < 3e-12,
                "An intervening mutation failed to invalidate controlled inverse restoration");
    }
}

inline void ExactComponents(const ExtendedStabilizer &state, const ExtendedFrame &expected, const char *message)
{
    Require(state.GetFrames().size() == 1 && state.GetFrames().front().amplitudes == expected.amplitudes &&
                state.GetFrames().front().signs == expected.signs,
            message);
}

inline void RequireTinyComponent(const ExtendedFrame &frame)
{
    bool found = false;
    for (const auto amplitude : frame.amplitudes)
    {
        const double magnitude = std::max(std::abs(amplitude.real()), std::abs(amplitude.imag()));
        found = found || (magnitude > 0. && magnitude < 1e-190);
    }
    Require(found, "History precision setup lacks the intended tiny nonzero coefficient");
}

inline void ApplyRaw(ExtendedStabilizer &state, const Gates::AppliedGate<> &gate)
{
    state.ApplyGate(
        Gates::AppliedGate<>(gate.getRawOperatorMatrix(), gate.getQubit1(), gate.getQubit2(), gate.getQubit3()));
}

inline void NestedHistory(bool tiny)
{
    ExtendedStabilizer state(3), reference(3);
    Prepare(state, tiny);
    Prepare(reference, tiny);
    const ExtendedFrame initial = state.GetFrames().front();
    if (tiny)
        RequireTinyComponent(initial);
    const std::array<Gates::AppliedGate<>, 3> gates{
        Gates::AppliedGate<>(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), 2, 0, 1),
        Gates::AppliedGate<>(Gates::XXPlusYYGate<>(.37, -.29), 0, 2), Gates::AppliedGate<>(Gates::TGate<>(), 1)};
    std::vector<ExtendedFrame> snapshots{initial};
    for (const auto &gate : gates)
    {
        state.ApplyGate(gate);
        ApplyRaw(reference, gate);
        // State queries and GetFrames must already expose each forward gate.
        Compare(state, reference);
        snapshots.push_back(state.GetFrames().front());
    }
    for (size_t count = gates.size(); count != 0; --count)
    {
        state.ApplyGate(gates[count - 1].Adjoint());
        ExactComponents(state, snapshots[count - 1], "Nested inverse failed to replay the retained history exactly");
    }
    ExactComponents(state, initial, "Nested typed uncompute lost the exact original coefficients");

    // Provenance on an AppliedGate adjoint must also work when it comes first.
    for (const auto &gate : gates)
    {
        state.ApplyGate(gate.Adjoint());
        state.ApplyGate(gate);
        ExactComponents(state, initial, "Inverse-first typed AppliedGate did not cancel exactly");
    }
}

template <class Forward, class Inverse> inline void SemanticInversePair(const Forward &forward, const Inverse &inverse)
{
    ExtendedStabilizer state(3);
    Prepare(state, true);
    const ExtendedFrame initial = state.GetFrames().front();
    RequireTinyComponent(initial);
    state.ApplyGate(forward, 0, 1, 2);
    state.ApplyGate(inverse, 0, 1, 2);
    ExactComponents(state, initial, "An exact semantic inverse failed to retain tiny input coefficients");
}

inline void SemanticInverses()
{
    SemanticInversePair(Gates::TGate<>(), Gates::TDGGate<>());
    SemanticInversePair(Gates::TDGGate<>(), Gates::TGate<>());
    SemanticInversePair(Gates::XXPlusYYGate<>(.43, -.29), Gates::XXPlusYYGate<>(-.43, -.29));
    SemanticInversePair(Gates::XXMinusYYGate<>(.43, -.29), Gates::XXMinusYYGate<>(-.43, -.29));
    SemanticInversePair(Gates::UGate<>(.43, 2.54, 2.707, 2.863), Gates::UGate<>(-.43, -2.707, -2.54, -2.863));
    SemanticInversePair(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), Gates::CCUGate<>(-.43, -2.707, -2.54, -2.863));
    ExtendedStabilizer controls(3);
    Prepare(controls, true);
    const ExtendedFrame initial = controls.GetFrames().front();
    controls.ApplyGate(Gates::CCHGate<>(), 2, 0, 1);
    controls.ApplyGate(Gates::CCHGate<>(), 2, 1, 0);
    ExactComponents(controls, initial, "Swapping symmetric controls prevented exact inverse cancellation");
}

template <class Gate> inline Gates::AppliedGate<> WithControls(const Gate &gate, unsigned controls)
{
    if (controls == 1)
        return Gates::AppliedGate<>(Gates::ControlledGate<>(gate), 0, 1);
    if (controls == 2)
        return Gates::AppliedGate<>(Gates::DoublyControlledGate<>(gate), 0, 1, 2);
    return Gates::AppliedGate<>(gate, 0);
}

inline void CrossFamilySequence(const Gates::AppliedGate<> &first, const Gates::AppliedGate<> &inverse, bool tiny,
                                bool rawInverse)
{
    ExtendedStabilizer state(4), reference(4);
    Prepare(state, tiny);
    Prepare(reference, tiny);
    state.ApplyCX(3, 0);
    reference.ApplyCX(3, 0);
    const ExtendedFrame initial = state.GetFrames().front();
    if (tiny)
        RequireTinyComponent(initial);
    const Gates::AppliedGate<> middle(Gates::TGate<>(), 3);
    state.ApplyGate(first);
    state.ApplyGate(middle);
    reference.ApplyGate(middle);
    if (rawInverse)
        ApplyRaw(state, inverse);
    else
        state.ApplyGate(inverse);
    ExactComponents(state, reference.GetFrames().front(),
                    "A cross-family inverse failed to replay the disjoint middle gate exactly");
    Compare(state, reference);
    state.ApplyGate(middle.Adjoint());
    ExactComponents(state, initial, "Cross-family cancellation lost the original tiny components");
}

inline void CrossFamilyHistory()
{
    const double quarterPi = std::acos(-1.) / 4.;
    for (unsigned controls = 0; controls <= 2; ++controls)
    {
        const std::array<std::pair<Gates::AppliedGate<>, Gates::AppliedGate<>>, 2> pairs{
            {{WithControls(Gates::PhaseShiftGate<>(quarterPi), controls), WithControls(Gates::TDGGate<>(), controls)},
             {WithControls(Gates::UGate<>(.43, 0., 0., 0.), controls), WithControls(Gates::RyGate<>(-.43), controls)}}};
        for (const auto &pair : pairs)
            for (bool tiny : {false, true})
                for (bool rawInverse : {false, true})
                {
                    CrossFamilySequence(pair.first, pair.second, tiny, rawInverse);
                    CrossFamilySequence(pair.second, pair.first, tiny, rawInverse);
                    CrossFamilySequence(pair.first.Adjoint(), pair.second.Adjoint(), tiny, rawInverse);
                    CrossFamilySequence(pair.second.Adjoint(), pair.first.Adjoint(), tiny, rawInverse);
                }
    }

    // A nonzero Euler phase must not disappear when U is compared with Ry,
    // even when that phase is too small for its squared norm to be representable.
    // U(-t, phi, 0) Ry(t) is diag(1, e^(i phi)): the phase shows only on a
    // state with a |1> part, so the qubit starts in |+>. On |0> the product is
    // the identity, and its only trace would be the rounding of a sum that
    // cancels exactly, which the matrix path removes.
    ExtendedStabilizer state(4), reference(4);
    state.ApplyH(0);
    reference.ApplyH(0);
    const Gates::AppliedGate<> first(Gates::RyGate<>(.43), 0), nearby(Gates::UGate<>(-.43, 1e-200, 0., 0.), 0),
        middle(Gates::TGate<>(), 3);
    const ExtendedFrame initial = state.GetFrames().front();
    for (const auto &gate : {first, middle, nearby})
    {
        state.ApplyGate(gate);
        ApplyRaw(reference, gate);
    }
    ExactComponents(state, reference.GetFrames().front(), "A tiny cross-family phase was changed by history matching");
    Require(state.GetFrames().front().amplitudes != initial.amplitudes ||
                !(state.GetFrames().front().signs == initial.signs),
            "A genuine tiny cross-family phase was incorrectly cancelled");
    bool tinyPhase = false;
    for (const auto amplitude : state.GetFrames().front().amplitudes)
        tinyPhase = tinyPhase || (std::abs(amplitude.imag()) > 0. && std::abs(amplitude.imag()) < 1e-190);
    Require(tinyPhase, "Cross-family matching discarded the genuine tiny imaginary coefficient");
}

inline void PrepareMatrixInverse(ExtendedStabilizer &state)
{
    for (size_t qubit = 0; qubit < 4; ++qubit)
        state.ApplyH(qubit);
    // Put the tiny coefficient in the same quadrature as an ordinary Rz
    // contribution: forward arithmetic alone cannot retain it in that sum.
    state.ApplyRz(0, 1e-200);
    state.ApplyRy(1, .317);
    state.ApplyCX(2, 1);
    RequireTinyComponent(state.GetFrames().front());
}

inline void MatrixInverseSequence(const Gates::AppliedGate<> &first, const Gates::AppliedGate<> &inverse,
                                  bool rawInverse)
{
    Require((inverse.getRawOperatorMatrix().array() == first.getRawOperatorMatrix().adjoint().array()).all(),
            "Matrix inverse regression requires exactly equal stored adjoints");
    ExtendedStabilizer state(4), reference(4);
    PrepareMatrixInverse(state);
    PrepareMatrixInverse(reference);
    const ExtendedFrame initial = state.GetFrames().front();
    const Gates::AppliedGate<> middle(Gates::TGate<>(), 3);
    state.ApplyGate(first);
    state.ApplyGate(middle);
    reference.ApplyGate(middle);
    if (rawInverse)
        ApplyRaw(state, inverse);
    else
        state.ApplyGate(inverse);
    ExactComponents(state, reference.GetFrames().front(), "An exact matrix inverse missed retained history");
    RequireTinyComponent(state.GetFrames().front());
    state.ApplyGate(middle.Adjoint());
    ExactComponents(state, initial, "Exact matrix inverse restoration lost tiny original coefficients");
}

inline void MatrixInverseHistory()
{
    constexpr double angle = 1e-16;
    const double pi = std::acos(-1.);
    for (unsigned controls = 0; controls <= 2; ++controls)
    {
        const auto rz = WithControls(Gates::RzGate<>(angle), controls);
        const auto reflection = WithControls(Gates::UGate<>(.37, 0., pi), controls);
        const std::array<std::pair<Gates::AppliedGate<>, Gates::AppliedGate<>>, 3> pairs{
            {{rz, WithControls(Gates::UGate<>(0., 0., -angle, .5 * angle), controls)},
             {WithControls(Gates::UGate<>(0., .3, 0., 0.), controls),
              WithControls(Gates::UGate<>(0., -.3, 0., 0.), controls)},
             {reflection, reflection}}};
        for (const auto &pair : pairs)
            for (bool rawInverse : {false, true})
            {
                MatrixInverseSequence(pair.first, pair.second, rawInverse);
                MatrixInverseSequence(pair.second, pair.first, rawInverse);
                MatrixInverseSequence(pair.first.Adjoint(), pair.second.Adjoint(), rawInverse);
            }

        const std::array<Gates::AppliedGate<>, 3> noninverses{
            WithControls(Gates::UGate<>(0., 0., angle, -.5 * angle), controls),
            WithControls(Gates::UGate<>(0., 0., -std::nextafter(angle, 1.), .5 * angle), controls),
            WithControls(Gates::UGate<>(0., 0., -angle, pi), controls)};
        for (const auto &noninverse : noninverses)
        {
            Require(!(noninverse.getRawOperatorMatrix().array() == rz.getRawOperatorMatrix().adjoint().array()).all(),
                    "Negative inverse regression accidentally constructed an exact adjoint");
            ExtendedStabilizer state(4), rawReference(4);
            PrepareMatrixInverse(state);
            PrepareMatrixInverse(rawReference);
            for (const auto &gate : {rz, Gates::AppliedGate<>(Gates::TGate<>(), 3), noninverse})
            {
                state.ApplyGate(gate);
                ApplyRaw(rawReference, gate);
            }
            ExactComponents(state, rawReference.GetFrames().front(),
                            "Matrix inverse matching changed a forward, neighboring or phase-shifted gate");
        }
    }
}

inline void PrepareFredkin(ExtendedStabilizer &state)
{
    Prepare(state, true);
    state.ApplyRy(2, .213);
    state.ApplyCX(3, 0);
    RequireTinyComponent(state.GetFrames().front());
}

inline void FredkinSymmetryHistory()
{
    const Gates::AppliedGate<> forward(Gates::FredkinGate<>(), 0, 1, 2), swapped(Gates::FredkinGate<>(), 1, 0, 2),
        middle(Gates::TGate<>(), 3), overlapping(Gates::UGate<>(.41, .23, -.17), 0);
    ExtendedStabilizer state(4), reference(4);
    PrepareFredkin(state);
    PrepareFredkin(reference);
    const ExtendedFrame initial = state.GetFrames().front();

    for (bool rawInverse : {false, true})
    {
        state.ApplyGate(forward);
        if (rawInverse)
            ApplyRaw(state, swapped);
        else
            state.ApplyGate(swapped);
        ExactComponents(state, initial, "Swapped Fredkin targets failed adjacent exact restoration");
        state.ApplyGate(forward.Adjoint());
        state.ApplyGate(swapped);
        ExactComponents(state, initial, "Inverse-first Fredkin target symmetry lost tiny components");
    }

    state.ApplyGate(forward);
    const ExtendedFrame afterForward = state.GetFrames().front();
    state.ApplyGate(overlapping);
    state.ApplyGate(overlapping.Adjoint());
    ExactComponents(state, afterForward, "Nested Fredkin history failed to restore the forward state");
    state.ApplyGate(swapped);
    ExactComponents(state, initial, "Nested Fredkin target symmetry failed exact restoration");

    state.ApplyGate(forward);
    state.ApplyGate(middle);
    reference.ApplyGate(middle);
    state.ApplyGate(swapped);
    ExactComponents(state, reference.GetFrames().front(), "Fredkin cancellation lost a disjoint retained gate");
    Compare(state, reference);
    state.ApplyGate(middle.Adjoint());
    ExactComponents(state, initial, "Fredkin cancellation consumed the reusable anchor");

    state.ApplyGate(forward);
    state.ApplyH(3);
    reference.ApplyGate(middle.Adjoint());
    reference.ApplyH(3);
    state.ApplyGate(swapped);
    ExactComponents(state, reference.GetFrames().front(), "A disjoint Clifford invalidated Fredkin target symmetry");
    Compare(state, reference);
}

inline void FredkinHistoryBoundaries()
{
    const Gates::AppliedGate<> forward(Gates::FredkinGate<>(), 0, 1, 2), swapped(Gates::FredkinGate<>(), 1, 0, 2),
        changedControl(Gates::FredkinGate<>(), 0, 2, 1), middle(Gates::TGate<>(), 3);
    ExtendedStabilizer state(4), reference(4), original(4);
    PrepareFredkin(state);
    PrepareFredkin(reference);
    PrepareFredkin(original);
    const ExtendedFrame initial = state.GetFrames().front();
    for (const auto &gate : {forward, changedControl})
    {
        state.ApplyGate(gate);
        ApplyRaw(reference, gate);
    }
    Compare(state, reference);
    double changedProbability = 0.;
    for (size_t index = 0; index < 16; ++index)
        changedProbability = std::max(changedProbability, std::abs(state.getBasisStateProbability(index) -
                                                                   original.getBasisStateProbability(index)));
    Require(changedProbability > 1e-6, "Fredkin with a different control was incorrectly cancelled");

    ExtendedStabilizer savedState(4);
    PrepareFredkin(savedState);
    savedState.ApplyGate(forward);
    savedState.ApplyGate(middle);
    const ExtendedFrame saved = savedState.GetFrames().front();
    savedState.SaveState();
    auto clone = savedState.Clone();
    savedState.ApplyGate(middle.Adjoint());
    savedState.ApplyGate(swapped);
    ExactComponents(savedState, initial, "Fredkin history failed after saving the public state");
    savedState.RestoreState();
    ExactComponents(savedState, saved, "Fredkin cancellation modified the user's saved state");
    for (const auto &gate : {middle.Adjoint(), swapped})
    {
        savedState.ApplyGate(gate);
        clone->ApplyGate(gate);
        ExactComponents(savedState, clone->GetFrames().front(), "Restore/clone inherited stale Fredkin history");
        Compare(savedState, *clone);
    }

    // A raw perturbation breaks the target symmetry certificate. Start with
    // sparse support so its 1e-200 contribution is visible as a new component.
    ExtendedStabilizer filtered(4), rawReference(4);
    filtered.ApplyRy(0, 1e-200);
    rawReference.ApplyRy(0, 1e-200);
    const ExtendedFrame unfiltered = filtered.GetFrames().front();
    Gates::FredkinGate<> mutated;
    Eigen::MatrixXcd filter = mutated.getRawOperatorMatrix();
    filter(1, 0) = 1e-200;
    mutated.setOperator(filter);
    const Gates::AppliedGate<> perturbed(mutated, 1, 0, 2);
    for (const auto &gate : {forward, middle, perturbed})
    {
        filtered.ApplyGate(gate);
        ApplyRaw(rawReference, gate);
    }
    ExactComponents(filtered, rawReference.GetFrames().front(), "Fredkin symmetry changed a perturbed raw filter");
    Require(filtered.GetFrames().front().amplitudes != unfiltered.amplitudes ||
                !(filtered.GetFrames().front().signs == unfiltered.signs),
            "Fredkin history cancelled a genuine tiny matrix perturbation");
    RequireTinyComponent(filtered.GetFrames().front());
}

inline void DisjointHistory()
{
    const std::vector<std::function<void(ExtendedStabilizer &)>> cliffords{
        [](auto &state) { state.ApplyH(3); }, [](auto &state) { state.ApplyS(4); },
        [](auto &state) { state.ApplyCX(4, 3); }, [](auto &state) { state.ApplySwap(3, 4); }};
    for (const auto &clifford : cliffords)
    {
        ExtendedStabilizer state(5), reference(5);
        Prepare(state, true);
        Prepare(reference, true);
        // Entanglement across the two supports prevents a product-state-only
        // test from hiding a wrong relationship between anchor and basis.
        state.ApplyCX(3, 0);
        reference.ApplyCX(3, 0);
        state.ApplyH(4);
        reference.ApplyH(4);
        state.ApplyGate(Gates::CCHGate<>(), 2, 0, 1);
        clifford(state);
        clifford(reference);
        state.ApplyGate(Gates::CCHGate<>(), 2, 0, 1);
        ExactComponents(state, reference.GetFrames().front(),
                        "Disjoint Clifford changed the history anchor coefficients");
        Compare(state, reference);
    }
    ExtendedStabilizer state(4), reference(4);
    Prepare(state, true);
    Prepare(reference, true);
    state.ApplyCX(3, 0);
    reference.ApplyCX(3, 0);
    const ExtendedFrame initial = state.GetFrames().front();
    state.ApplyGate(Gates::CCHGate<>(), 2, 0, 1);
    state.ApplyGate(Gates::TGate<>(), 3);
    reference.ApplyGate(Gates::TGate<>(), 3);
    state.ApplyGate(Gates::CCHGate<>(), 2, 1, 0);
    ExactComponents(state, reference.GetFrames().front(), "Earlier inverse did not replay the disjoint remaining gate");
    Compare(state, reference);
    state.ApplyGate(Gates::TDGGate<>(), 3);
    ExactComponents(state, initial, "Disjoint commuting cancellation consumed the reusable anchor");
}

inline void NoncommutingHistory()
{
    ExtendedStabilizer state(3), reference(3);
    Prepare(state, true);
    Prepare(reference, true);
    const Gates::AppliedGate<> first(Gates::RxGate<>(.51), 0), middle(Gates::RzGate<>(.73), 0);
    for (const auto &gate : {first, middle, first.Adjoint(), middle.Adjoint()})
    {
        state.ApplyGate(gate);
        ApplyRaw(reference, gate);
        Compare(state, reference);
    }
    ExtendedStabilizer nearby(3);
    Prepare(nearby, true);
    const ExtendedFrame initial = nearby.GetFrames().front();
    nearby.ApplyGate(Gates::CCPhaseShiftGate<>(.3), 0, 1, 2);
    nearby.ApplyGate(Gates::CCPhaseShiftGate<>(-std::nextafter(.3, 1.)), 0, 1, 2);
    Require(nearby.GetFrames().front().amplitudes != initial.amplitudes ||
                !(nearby.GetFrames().front().signs == initial.signs),
            "A neighboring semantic inverse was incorrectly folded to identity");
}

inline void HistoryCapacity()
{
    for (size_t count : {size_t(32), size_t(33)})
    {
        ExtendedStabilizer state(3), reference(3);
        Prepare(state, true);
        Prepare(reference, true);
        const ExtendedFrame initial = state.GetFrames().front();
        std::vector<Gates::AppliedGate<>> gates;
        std::vector<ExtendedFrame> snapshots{initial};
        for (size_t index = 0; index < count; ++index)
        {
            gates.emplace_back(Gates::CCUGate<>(.11 + .007 * index, -.37 + .013 * index, .23 - .019 * index, .17), 0, 1,
                               2);
            state.ApplyGate(gates.back());
            ApplyRaw(reference, gates.back());
            snapshots.push_back(state.GetFrames().front());
        }
        Compare(state, reference);
        for (size_t remaining = count; remaining != 0; --remaining)
        {
            const auto inverse = gates[remaining - 1].Adjoint();
            state.ApplyGate(inverse);
            ApplyRaw(reference, inverse);
            Compare(state, reference);
            if (count == 32 || remaining == count)
                ExactComponents(state, snapshots[remaining - 1], "History capacity boundary lost an immediate inverse");
        }
        if (count == 32)
            ExactComponents(state, initial, "A full bounded history failed to restore its exact input");
        // Beyond the retained window, ordinary immediate execution remains
        // correct even though exact restoration of the old anchor is not promised.
        ExtendedStabilizer original(3);
        Prepare(original, true);
        Compare(state, original);
    }
}

inline void HistoryLifecycle()
{
    ExtendedStabilizer state(3);
    Prepare(state, true);
    const ExtendedFrame initial = state.GetFrames().front();
    const Gates::AppliedGate<> first(Gates::CCRxGate<>(.31), 0, 1, 2), second(Gates::CCPhaseShiftGate<>(.47), 0, 1, 2);
    state.ApplyGate(first);
    state.ApplyGate(second);
    const ExtendedFrame saved = state.GetFrames().front();
    state.SaveState();
    auto clone = state.Clone();
    bool rejected = false;
    try
    {
        state.ApplyGate(Gates::CCHGate<>(), 0, 0, 2);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    Require(rejected, "Invalid gate indices were accepted while history was active");
    ExactComponents(state, saved, "An invalid gate changed the materialized history state");
    state.ApplyGate(second.Adjoint());
    state.ApplyGate(first.Adjoint());
    ExactComponents(state, initial, "Invalid input or SaveState corrupted a multi-entry anchor");
    state.RestoreState();
    ExactComponents(state, saved, "History cancellation changed the user's saved frame");
    for (const auto &inverse : {second.Adjoint(), first.Adjoint()})
    {
        state.ApplyGate(inverse);
        clone->ApplyGate(inverse);
        ExactComponents(state, clone->GetFrames().front(), "Restore/clone inherited stale unitary history");
        Compare(state, *clone);
    }

    ExtendedStabilizer approximate(3, ExtendedStabilizerApproximationPolicy::Approximate(0., 1));
    approximate.ApplyRy(0, .317);
    approximate.SetApproximationPolicy(ExtendedStabilizerApproximationPolicy::Exact());
    approximate.ApplyX(1);
    approximate.ApplyX(2);
    const auto statistics = approximate.GetApproximationStatistics();
    Require(statistics.traceDistanceErrorBound > 0., "History statistics setup lacks prior approximation error");
    approximate.ApplyGate(first);
    approximate.ApplyGate(second);
    approximate.ApplyGate(second.Adjoint());
    approximate.ApplyGate(first.Adjoint());
    const auto after = approximate.GetApproximationStatistics();
    Require(after.traceDistanceErrorBound == statistics.traceDistanceErrorBound &&
                after.cumulativeDiscardedWeight == statistics.cumulativeDiscardedWeight &&
                after.discardedComponents == statistics.discardedComponents &&
                after.pruningEvents == statistics.pruningEvents,
            "Exact history restoration erased or inflated prior approximation statistics");

    // The quick overlap mask aliases qubits modulo the machine word width;
    // the following exact index check must distinguish physical qubits 0/64.
    ExtendedStabilizer wide(65);
    Prepare(wide, true);
    const ExtendedFrame wideInput = wide.GetFrames().front();
    wide.ApplyGate(Gates::TGate<>(), 0);
    wide.ApplyH(64);
    wide.ApplyGate(Gates::TDGGate<>(), 0);
    ExactComponents(wide, wideInput, "A modulo-mask collision incorrectly discarded a disjoint history anchor");
    Require(wide.GetQubitProbability(64) == .5, "Disjoint wide Clifford was lost while restoring the anchor");
}

inline void InterveningMutations()
{
    const Gates::CCHGate<> gate;
    const Gates::AppliedGate<> raw(gate.getRawOperatorMatrix(), 0, 1, 2);
    const std::vector<std::function<void(ExtendedStabilizer &)>> mutations = {
        [](auto &s) { s.ApplyH(0); },
        [](auto &s) { s.ApplyS(0); },
        [](auto &s) { s.ApplySdg(0); },
        [](auto &s) { s.ApplyX(0); },
        [](auto &s) { s.ApplyY(0); },
        [](auto &s) { s.ApplyZ(0); },
        [](auto &s) { s.ApplySx(0); },
        [](auto &s) { s.ApplySxDag(0); },
        [](auto &s) { s.ApplyK(0); },
        [](auto &s) { s.ApplyCX(0, 1); },
        [](auto &s) { s.ApplyCY(0, 1); },
        [](auto &s) { s.ApplyCZ(0, 1); },
        [](auto &s) { s.ApplySwap(0, 1); },
        [](auto &s) { s.ApplyISwap(0, 1); },
        [](auto &s) { s.ApplyISwapDag(0, 1); },
        [](auto &s) { s.ApplyRx(0, .31); },
        [](auto &s) { s.ApplyRy(0, .31); },
        [](auto &s) { s.ApplyRz(0, .31); },
        [](auto &s) { s.ApplyRy(0, M_PI_2); },
        [](auto &s) { s.ApplyCH(0, 1); },
        [](auto &s) { s.ApplyRXX(0, 1, .31); },
        [](auto &s) { s.ApplyGate(Gates::HadamardGate<>(), 0); },
        [](auto &s) { s.ApplyGate(Gates::UGate<>(.31, .17, -.52), 0); },
        [](auto &s) { s.MeasureConditioned(0, false); },
        [](auto &s) { s.setToBasisState(3); },
        [](auto &s) { s.Reset(3); },
        [](auto &s) { s.SetApproximationPolicy(ExtendedStabilizerApproximationPolicy::Approximate(0., 2)); }};
    for (const auto &mutate : mutations)
    {
        ExtendedStabilizer actual(3), reference(3);
        Prepare(actual, true);
        Prepare(reference, true);
        actual.ApplyGate(gate, 0, 1, 2);
        reference.ApplyGate(raw);
        mutate(actual);
        mutate(reference);
        actual.ApplyGate(gate, 0, 1, 2);
        reference.ApplyGate(raw);
        Compare(actual, reference);
    }
    ExtendedStabilizer actual(3), reference(3);
    Prepare(actual, true);
    Prepare(reference, true);
    actual.ApplyGate(gate, 0, 1, 2);
    reference.ApplyGate(raw);
    actual.ApplyGate(gate, 1, 0, 2);
    reference.ApplyGate(raw, 1, 0, 2);
    Compare(actual, reference);
}

inline void FiltersAndNearbyInverse()
{
    ExtendedStabilizer actual(3), reference(3);
    Prepare(actual, true);
    Prepare(reference, true);
    const Gates::CCRxGate<> first(.3), nearby(-.3 + 1e-12);
    actual.ApplyGate(first, 0, 1, 2);
    reference.ApplyGate(Gates::AppliedGate<>(first.getRawOperatorMatrix(), 0, 1, 2));
    actual.ApplyGate(nearby, 0, 1, 2);
    reference.ApplyGate(Gates::AppliedGate<>(nearby.getRawOperatorMatrix(), 0, 1, 2));
    Compare(actual, reference);

    Gates::CCHGate<> mutated;
    Eigen::MatrixXcd filter = mutated.getRawOperatorMatrix();
    filter(0, 0) = 2.;
    mutated.setOperator(filter);
    actual.ApplyGate(mutated, 0, 1, 2);
    reference.ApplyGate(Gates::AppliedGate<>(filter, 0, 1, 2));
    actual.ApplyGate(Gates::AppliedGate<>(filter.adjoint().eval(), 0, 1, 2));
    reference.ApplyGate(Gates::AppliedGate<>(filter.adjoint().eval(), 0, 1, 2));
    Compare(actual, reference);

    ExtendedStabilizer approximate(3, ExtendedStabilizerApproximationPolicy::Approximate(0., 1));
    approximate.ApplyRy(0, .317);
    const double before = approximate.GetApproximationErrorBound();
    Require(before > 0. && before < 1., "Approximation-bound regression setup failed");
    approximate.SetApproximationPolicy(ExtendedStabilizerApproximationPolicy::Exact());
    approximate.ApplyGate(first, 0, 1, 2);
    Require(approximate.GetApproximationErrorBound() == before,
            "A verified unitary was treated as arbitrary postselection");
    approximate.ApplyGate(Gates::AppliedGate<>(filter, 0, 1, 2));
    Require(approximate.GetApproximationErrorBound() == 1.,
            "A genuine state filter understated its prior approximation error");
}

// A matrix that shows a control exactly (the identity where the control does
// not select, nothing coupling the two values) acts only on the selected
// part: an unselected state must come back bit for bit, with no component
// left by cancelling products.
inline void ExactControlStructure()
{
    // The control of every gate is zero: |0...0> stays one component.
    ExtendedStabilizer chain(32);
    double angle = -2.9;
    const auto next = [&angle] { return angle = std::fmod(angle + 1.7183, 6.) - 3.; };
    for (size_t q = 0; q + 1 < 32; ++q)
    {
        const double theta = next(), phi = next(), lambda = next(), gamma = next();
        chain.ApplyGate(Gates::ControlledUGate<>(theta, phi, lambda, gamma), q + 1, q);
    }
    Require(chain.GetFrames().front().GetFrameSize() == 1 &&
                chain.GetFrames().front().amplitudes.front() == std::complex<double>(1.),
            "Controlled gates with an unset control changed the state or left components");

    ExtendedStabilizer single(2), doubly(3), euler(3);
    single.ApplyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 1, 0);
    doubly.ApplyDoublyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 2, 0, 1);
    euler.ApplyGate(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), 2, 0, 1);
    euler.ApplyX(0);
    euler.ApplyGate(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), 2, 0, 1);
    Require(single.GetFrames().front().GetFrameSize() == 1 && doubly.GetFrames().front().GetFrameSize() == 1 &&
                euler.GetFrames().front().GetFrameSize() == 1,
            "A controlled wrapper left a component for an unset control");

    // An unselected tiny component is not touched, even next to a large one.
    ExtendedStabilizer tiny(3);
    tiny.ApplyRy(0, 1e-200);
    tiny.ApplyRy(2, .317);
    const ExtendedFrame before = tiny.GetFrames().front();
    RequireTinyComponent(before);
    tiny.ApplyGate(Gates::ControlledUGate<>(.4, .2, -.3, .1), 0, 1);
    tiny.ApplyGate(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), 0, 1, 2);
    ExactComponents(tiny, before, "A controlled gate changed components its control does not select");

    // Raw matrices: the control on the low local qubit, a control that
    // selects zero, and two controls of different values. Each must act as
    // the native gate between X gates, and leave an unselected state alone.
    const Eigen::MatrixXcd u = Gates::UGate<>(.4, .2, -.3, .1).getRawOperatorMatrix();
    const auto controlled = [&u](unsigned width, unsigned target, unsigned selecting) {
        // U on local qubit `target` where the other local qubits equal the bits of `selecting`
        const unsigned dimension = 1U << width;
        Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Identity(dimension, dimension);
        for (unsigned row = 0; row < dimension; ++row)
            for (unsigned col = 0; col < dimension; ++col)
            {
                const unsigned others = ~(1U << target) & (dimension - 1);
                if ((row & others) != (selecting & others) || (col & others) != (selecting & others))
                    continue;
                matrix(row, col) = u((row >> target) & 1U, (col >> target) & 1U);
            }
        return matrix;
    };
    struct Case
    {
        unsigned width, target, selecting;
    };
    for (const Case &c : {Case{2, 1, 1}, Case{2, 1, 0}, Case{2, 0, 0}, Case{3, 0, 2}, Case{3, 1, 4}, Case{3, 2, 0}})
    {
        const Eigen::MatrixXcd matrix = controlled(c.width, c.target, c.selecting);
        ExtendedStabilizer actual(3), reference(3);
        Prepare(actual, false);
        Prepare(reference, false);
        actual.ApplyRy(1, .417);
        reference.ApplyRy(1, .417);
        const size_t qubits[] = {2, 0, 1};
        actual.ApplyGate(Gates::AppliedGate<>(matrix, qubits[0], qubits[1], qubits[2]));
        // the same operator from X gates and the wrappers on the selected value one
        for (unsigned bit = 0; bit < c.width; ++bit)
            if (bit != c.target && !((c.selecting >> bit) & 1U))
                reference.ApplyX(qubits[bit]);
        if (c.width == 2)
            reference.ApplyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), qubits[c.target], qubits[1 - c.target]);
        else
            reference.ApplyDoublyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), qubits[c.target],
                                                qubits[(c.target + 1) % 3], qubits[(c.target + 2) % 3]);
        for (unsigned bit = 0; bit < c.width; ++bit)
            if (bit != c.target && !((c.selecting >> bit) & 1U))
                reference.ApplyX(qubits[bit]);
        Compare(actual, reference);

        // every control on its selecting value except the first one after the target
        ExtendedStabilizer unselected(3);
        const unsigned wrong = (c.target + 1) % c.width;
        for (unsigned bit = 0; bit < c.width; ++bit)
            if (bit != c.target && bool((c.selecting >> bit) & 1U) != (bit == wrong))
                unselected.ApplyX(qubits[bit]);
        unselected.ApplyRy(qubits[c.target], .3);
        const ExtendedFrame untouched = unselected.GetFrames().front();
        unselected.ApplyGate(Gates::AppliedGate<>(matrix, qubits[0], qubits[1], qubits[2]));
        ExactComponents(unselected, untouched, "A raw controlled matrix changed a state its control does not select");
    }
}

// A result that is the rounding of its own cancelling sum leaves no
// component; a small result with no cancellation behind it stays.
inline void RoundedCancellations()
{
    ExtendedStabilizer phases(1);
    phases.ApplyH(0);
    phases.ApplyGate(Gates::TGate<>(), 0);
    phases.ApplyS(0); // an overlapping Clifford: the history cannot cancel the pair
    phases.ApplyGate(Gates::TDGGate<>(), 0);
    Require(phases.GetFrames().front().GetFrameSize() == 1, "T, S, Tdg left a cancelled component behind");

    // Compute and uncompute through the matrix path with the history cleared in between.
    ExtendedStabilizer pair(3);
    Prepare(pair, false);
    pair.ApplyGate(Gates::CCUGate<>(.43, 2.54, 2.707, 2.863), 2, 0, 1);
    Require(pair.GetFrames().front().GetFrameSize() > 1, "Cancellation regression setup did not branch");
    pair.ApplyZ(2);
    pair.ApplyZ(2);
    pair.ApplyGate(Gates::CCUGate<>(-.43, -2.707, -2.54, -2.863), 2, 0, 1);
    Require(pair.GetFrames().front().GetFrameSize() == 1, "An uncomputed controlled gate left components behind");
    for (unsigned basis = 0; basis < 8; ++basis)
        Require(std::abs(pair.getBasisStateProbability(basis) - .125) < 2e-13,
                "Clearing rounded cancellations changed the state");

    // The residual of a neighbouring inverse is genuine: about 5e-16, several
    // times the rounding bound of its two products. It must survive in the
    // flat path and behind a control.
    ExtendedStabilizer flat(1), behind(2);
    flat.ApplyGate(Gates::RyGate<>(.731), 0);
    flat.ApplyZ(0);
    flat.ApplyZ(0);
    flat.ApplyGate(Gates::RyGate<>(-.731 + 1e-15), 0);
    Require(flat.GetFrames().front().GetFrameSize() == 2 && std::abs(flat.ExpectationValue("X")) > 2e-16,
            "The matrix path discarded a genuine near-cancellation");
    behind.ApplyX(1);
    behind.ApplyGate(Gates::ControlledRyGate<>(.731), 0, 1);
    behind.ApplyZ(0);
    behind.ApplyZ(0);
    behind.ApplyGate(Gates::ControlledRyGate<>(-.731 + 1e-15), 0, 1);
    Require(behind.GetFrames().front().GetFrameSize() == 2 && std::abs(behind.ExpectationValue("XI")) > 2e-16,
            "A controlled matrix discarded a genuine near-cancellation");

    // No pruned residue: the approximation statistics stay untouched, also
    // after a measurement, which would promote any positive bound to one.
    ExtendedStabilizer approximate(3, ExtendedStabilizerApproximationPolicy::Approximate(1e-9));
    approximate.ApplyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 1, 0);
    approximate.ApplyDoublyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 2, 0, 1);
    approximate.ApplyH(2);
    approximate.SetSeed(1);
    approximate.Measure(2);
    const auto &statistics = approximate.GetApproximationStatistics();
    Require(statistics.traceDistanceErrorBound == 0. && statistics.pruningEvents == 0 &&
                statistics.discardedComponents == 0 && approximate.GetFrames().front().GetFrameSize() == 1,
            "A controlled wrapper left a residue for the approximation to prune");
}

// Inverses that arrive in an arbitrary admissible order: the components must
// always be those of a simulator that only applied the surviving gates.
inline void HistoryReplayOrders()
{
    struct Step
    {
        double angle;
        size_t qubit;
        bool x;
    };
    std::vector<Step> surviving;
    const auto apply = [](ExtendedStabilizer &s, const Step &step) {
        if (step.x)
            s.ApplyGate(Gates::RxGate<>(step.angle), step.qubit);
        else
            s.ApplyGate(Gates::RzGate<>(step.angle), step.qubit);
    };
    const auto check = [&](const ExtendedStabilizer &state, const char *message) {
        // A Z pair after each gate keeps the reference from retaining any history.
        ExtendedStabilizer reference(4);
        Prepare(reference, true);
        reference.ApplyCX(3, 0);
        for (const Step &step : surviving)
        {
            apply(reference, step);
            reference.ApplyZ(step.qubit);
            reference.ApplyZ(step.qubit);
        }
        ExactComponents(state, reference.GetFrames().front(), message);
    };
    ExtendedStabilizer state(4);
    Prepare(state, true);
    state.ApplyCX(3, 0);
    for (size_t k = 0; k < 28; ++k)
    {
        surviving.push_back({.1 + .013 * double(k), k % 4, (k / 4) % 2 == 1});
        apply(state, surviving.back());
    }
    check(state, "A retained sequence of matrix gates differs from its forward application");
    // Qubit by qubit, last gate first: each inverse matches the last entry on
    // its qubit, which is rarely the last entry of the history.
    for (size_t qubit : {size_t(2), size_t(0), size_t(3), size_t(1)})
        for (size_t k = surviving.size(); k-- > 0;)
        {
            if (surviving[k].qubit != qubit)
                continue;
            Step inverse = surviving[k];
            inverse.angle = -inverse.angle;
            surviving.erase(surviving.begin() + k);
            apply(state, inverse);
            check(state, "Removing a retained gate did not reproduce the remaining gates exactly");
        }
    Require(surviving.empty(), "History order regression did not remove every gate");
}

// A raw matrix that is unitary to rounding is not a state filter for the
// error bound; this decides the bookkeeping only.
inline void RawUnitaryBound()
{
    const Eigen::MatrixXcd unitary = Gates::CCUGate<>(.43, 2.54, 2.707, 2.863).getRawOperatorMatrix();
    for (double scale : {1., 1. + 1e-13, 1. + 1e-9})
    {
        ExtendedStabilizer approximate(3, ExtendedStabilizerApproximationPolicy::Approximate(0., 1)), exact(3);
        approximate.ApplyRy(0, .317);
        exact.ApplyRy(0, .317);
        const double before = approximate.GetApproximationErrorBound();
        Require(before > 0. && before < 1., "Raw unitary bound regression setup failed");
        approximate.SetApproximationPolicy(ExtendedStabilizerApproximationPolicy::Exact());
        const Eigen::MatrixXcd raw = scale * unitary;
        approximate.ApplyGate(Gates::AppliedGate<>(raw, 0, 1, 2));
        Require(approximate.GetApproximationErrorBound() == (scale < 1. + 1e-10 ? before : 1.),
                "The error bound mistook a raw unitary for a filter, or a filter for a unitary");
        // The tolerance never reaches the amplitudes: the same matrix gives
        // the same components whatever the bound says.
        ExtendedStabilizer first(3), second(3, ExtendedStabilizerApproximationPolicy::Approximate(0., 64));
        first.ApplyRy(0, .317);
        second.ApplyRy(0, .317);
        first.ApplyGate(Gates::AppliedGate<>(raw, 0, 1, 2));
        second.ApplyGate(Gates::AppliedGate<>(raw, 0, 1, 2));
        ExactComponents(first, second.GetFrames().front(), "The unitarity tolerance changed an amplitude");
    }
}

// Compact gates on at most three qubits act through their dense matrices.
inline void CompactGates()
{
    using Applied = Gates::AppliedGate<>;
    const Eigen::MatrixXcd u = Gates::UGate<>(.4, .2, -.3, .1).getRawOperatorMatrix();
    const std::complex<double> phase = std::polar(1., .3);
    const auto compare = [](const Applied &compact, const std::function<void(ExtendedStabilizer &)> &native) {
        ExtendedStabilizer actual(4), reference(4);
        Prepare(actual, false);
        Prepare(reference, false);
        actual.ApplyRy(1, .417);
        reference.ApplyRy(1, .417);
        actual.ApplyGate(compact);
        native(reference);
        Compare(actual, reference);
    };
    compare(Applied::MultiControlled(u, 2, {0, 1}),
            [](auto &s) { s.ApplyDoublyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 2, 0, 1); });
    compare(Applied::MultiControlled(u, 2, {0}),
            [](auto &s) { s.ApplyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 2, 0); });
    compare(Applied::MultiControlled(u, 1, {}), [](auto &s) { s.ApplyGate(Gates::UGate<>(.4, .2, -.3, .1), 1); });
    compare(Applied::MultiControlled(u, 2, {0, 1}, 1), [](auto &s) {
        // selected where qubit 0 is one and qubit 1 is zero
        s.ApplyX(1);
        s.ApplyDoublyControlledGate(Gates::UGate<>(.4, .2, -.3, .1), 2, 0, 1);
        s.ApplyX(1);
    });
    compare(Applied::Diagonal({1., 1., 1., phase}, {0, 2}), [](auto &s) { s.ApplyCP(2, 0, .3); });
    compare(Applied::MultiplexedRotation(Gates::RotationAxis::Y, {.21, .58}, 1, {0}), [](auto &s) {
        s.ApplyRy(1, .21);
        s.ApplyCRy(1, 0, .58 - .21);
    });
    compare(Applied::MultiControlled(Gates::PauliXGate<>().getRawOperatorMatrix(), 3, {1}),
            [](auto &s) { s.ApplyCX(3, 1); });

    // A compact gate with an unset control leaves one component, as its dense form does.
    ExtendedStabilizer unset(3);
    unset.ApplyGate(Applied::MultiControlled(u, 2, {0, 1}));
    unset.ApplyGate(Applied::MultiControlled(Gates::PauliXGate<>().getRawOperatorMatrix(), 2, {0, 1}));
    Require(unset.GetFrames().front().GetFrameSize() == 1, "A compact controlled gate left components behind");

    // More than three qubits are rejected before any change.
    ExtendedStabilizer wide(5);
    Prepare(wide, true);
    const ExtendedFrame before = wide.GetFrames().front();
    bool rejected = false;
    try
    {
        wide.ApplyGate(Applied::MultiControlled(u, 4, {0, 1, 2}));
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    Require(rejected, "A compact gate on four qubits was accepted");
    ExactComponents(wide, before, "A rejected compact gate changed the state");
    rejected = false;
    try
    {
        wide.ApplyGate(Applied::Diagonal(std::vector<std::complex<double>>(16, phase), {0, 1, 2, 3}));
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    Require(rejected, "A compact diagonal on four qubits was accepted");
    ExactComponents(wide, before, "A rejected compact diagonal changed the state");
}

inline void Run()
{
    for (bool tiny : {false, true})
    {
        RoundTrip(Gates::CCHGate<>(), tiny);
        RoundTrip(Gates::ControlledHadamardGate<>(), tiny);
        RoundTrip(Gates::CCRxGate<>(.3), tiny);
        RoundTrip(Gates::CCUGate<>(.43, -.71, .19, .23), tiny);
        RoundTrip(Gates::CCPhaseShiftGate<>(.3), tiny);
    }
    ExtendedStabilizer wrapped(3);
    Prepare(wrapped, false);
    wrapped.ApplyDoublyControlledGate(Gates::HadamardGate<>(), 0, 1, 2);
    wrapped.ApplyDoublyControlledGate(Gates::HadamardGate<>(), 0, 1, 2);
    Require(wrapped.GetFrames().front().GetFrameSize() == 1,
            "Controlled wrappers discarded the verified constructor descriptor");
    InterveningMutations();
    FiltersAndNearbyInverse();
    NestedHistory(false);
    NestedHistory(true);
    SemanticInverses();
    CrossFamilyHistory();
    MatrixInverseHistory();
    FredkinSymmetryHistory();
    FredkinHistoryBoundaries();
    DisjointHistory();
    NoncommutingHistory();
    HistoryCapacity();
    HistoryLifecycle();
    ExactControlStructure();
    RoundedCancellations();
    HistoryReplayOrders();
    RawUnitaryBound();
    CompactGates();
    std::cout << "Extended bounded unitary history, mutation and precision regressions passed.\n";
}
} // namespace ExtendedUnitaryRegression
} // namespace QC
