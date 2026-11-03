#pragma once

#include "PauliPropagationEngine.h"
#include <cstring>
#include <iostream>
#include <map>
#include <random>
#include <string>

namespace QC
{
namespace PauliRuntimeEngineRegression
{
inline void Require(bool value, const char *message)
{
    if (!value)
        throw std::runtime_error(message);
}

template <class T> std::vector<T> ReferenceDedup(const std::vector<T> &input, const PauliDetail::Settings &s)
{
    using Key = std::pair<std::vector<uint64_t>, std::vector<uint64_t>>;
    std::map<Key, size_t> slots;
    std::vector<T> result;
    for (const auto &term : input)
    {
        if (s.weight < s.qubits && PauliDetail::Weight(term) > s.weight)
            continue;
        Key key{{term.X.begin(), term.X.end()}, {term.Z.begin(), term.Z.end()}};
        const auto inserted = slots.emplace(std::move(key), result.size());
        if (inserted.second)
            result.push_back(term);
        else
            result[inserted.first->second].Coefficient += term.Coefficient;
    }
    result.erase(std::remove_if(result.begin(), result.end(),
                                [&](const T &term) {
                                    return !(std::abs(term.Coefficient) > s.coefficient &&
                                             (s.weight >= s.qubits || PauliDetail::Weight(term) <= s.weight));
                                }),
                 result.end());
    return result;
}

template <class T> void Same(const std::vector<T> &actual, const std::vector<T> &expected)
{
    Require(actual.size() == expected.size(), "runtime engine output size mismatch");
    for (size_t i = 0; i < actual.size(); ++i)
    {
        Require(actual[i].X == expected[i].X && actual[i].Z == expected[i].Z,
                "runtime engine changed first-occurrence order");
        Require(std::memcmp(&actual[i].Coefficient, &expected[i].Coefficient, sizeof(double)) == 0,
                "runtime engine changed coefficient addition order");
    }
}

template <class T> std::vector<T> Input(size_t width, size_t count, size_t unique)
{
    std::vector<T> result;
    result.reserve(count);
    for (size_t i = 0; i < count; ++i)
    {
        T term(width);
        const size_t label = (i * 251 + 19) % unique;
        term.X[0] = (label & 0xff) << 1;
        term.Z[0] = (label >> 8) | 1;
        if (width > 64)
        {
            term.X[(width - 1) / 64] |= uint64_t(label & 1) << ((width - 1) % 64);
            term.Z[(width - 1) / 64] |= uint64_t((label >> 1) & 1) << ((width - 1) % 64);
        }
        // Cancellation and unequal magnitudes expose reassociation mistakes.
        term.Coefficient = i % 7 == 0 ? 1e15 : i % 7 == 1 ? -1e15 : double(int(i % 11) - 5) * .125;
        result.push_back(std::move(term));
    }
    return result;
}

template <class T> void DedupChecks(size_t width)
{
    PauliDetail::Workspace<T> ws;
    PauliDetail::Settings settings;
    settings.qubits = width;
    for (size_t count : {size_t(0), size_t(7), size_t(513), size_t(8192), size_t(9), size_t(2048)})
        for (size_t unique : {size_t(1), size_t(13), size_t(257), size_t(8192)})
            for (double cutoff : {0., .3})
                for (size_t weight : {size_t(4), std::numeric_limits<size_t>::max()})
                {
                    settings.coefficient = cutoff;
                    settings.weight = weight;
                    ws.terms = Input<T>(width, count, unique);
                    const auto expected = ReferenceDedup(ws.terms, settings);
                    PauliDetail::Deduplicate(ws, settings);
                    Same(ws.terms, expected);
                }
    // A fresh duplicate-heavy table should not reserve a full-frontier index.
    PauliDetail::Workspace<T> duplicates;
    duplicates.terms = Input<T>(width, 16384, 1);
    settings.coefficient = 0.;
    settings.weight = std::numeric_limits<size_t>::max();
    PauliDetail::Deduplicate(duplicates, settings);
    Require(duplicates.hashSlots.size() <= 128, "duplicate-heavy dedup capacity hint overallocated");
}

template <class T> void AlternatingDedupChecks(size_t width)
{
    PauliDetail::Workspace<T> ws;
    PauliDetail::Settings settings;
    settings.qubits = width;
    const std::array<std::pair<size_t, size_t>, 3> shapes{{{16384, 16384}, {2048, 2048}, {16384, 64}}};
    for (unsigned cycle = 0; cycle < 4; ++cycle)
        for (const auto &shape : shapes)
        {
            ws.terms = Input<T>(width, shape.first, shape.second);
            const auto expected = ReferenceDedup(ws.terms, settings);
            const size_t rehashes = ws.dedupRehashCount;
            PauliDetail::Deduplicate(ws, settings);
            Same(ws.terms, expected);
            if (cycle && shape.first == 16384 && shape.second == 16384)
                Require(ws.dedupRehashCount == rehashes,
                        "retained dedup table repeated intermediate rehashes after shrinking");
        }
    Require(ws.hashSlots.capacity() >= 32768, "alternating dedup discarded retained storage");
    Require(ws.hashSlots.size() < ws.hashSlots.capacity(), "alternating dedup did not exercise logical shrinking");
}

template <class T> void SmallRetainedDedupChecks(size_t width)
{
    auto large = Input<T>(width, 16384, 16384);
    const size_t first = static_cast<size_t>(PauliDetail::Mix(UINT64_C(0x9e3779b97f4a7c15))) % large.size();
    const size_t second = static_cast<size_t>(PauliDetail::Mix(UINT64_C(0x9e3779b97f4a7c16))) % large.size();
    large[second].X = large[first].X;
    large[second].Z = large[first].Z;
    std::array<size_t, 64> hashes;
    size_t collisions = 0;
    for (size_t sample = 0; sample < hashes.size(); ++sample)
    {
        const size_t index =
            static_cast<size_t>(PauliDetail::Mix(sample + UINT64_C(0x9e3779b97f4a7c15))) % large.size();
        hashes[sample] = PauliDetail::Hash(large[index]);
        for (size_t previous = 0; previous < sample; ++previous)
            collisions += hashes[previous] == hashes[sample];
    }
    Require(collisions == 1, "small-retained regression must have exactly one sampled duplicate");
    PauliDetail::Settings settings;
    settings.qubits = width;
    PauliDetail::Workspace<T> ws;
    ws.terms = large;
    PauliDetail::Deduplicate(ws, settings);
    ws.terms = Input<T>(width, 7, 7);
    PauliDetail::Deduplicate(ws, settings);
    Require(ws.hashSlots.size() == 64 && ws.hashSlots.capacity() >= 32768,
            "small-retained regression did not reach the intended logical-size branch");
    ws.terms = large;
    const auto expected = ReferenceDedup(ws.terms, settings);
    const size_t rehashes = ws.dedupRehashCount;
    PauliDetail::Deduplicate(ws, settings);
    Same(ws.terms, expected);
    Require(ws.dedupRehashCount == rehashes, "small logical table bypassed conservative retained-capacity planning");
}

class Revive : public OperatorX
{
  public:
    Revive() : OperatorX(0)
    {
    }

    void Apply(PauliStringXZWithCoefficient &term, PauliStringStorage &) const override
    {
        term.Coefficient += .125;
    }

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<Revive>(*this);
    }
};

template <class T> void ScalarBlockChecks(size_t width)
{
    using Operation = PauliDetail::Operation;
    using Kind = OperationType;
    auto input = Input<T>(width, 16384, 4096);
    for (size_t i = 0; i < input.size(); i += 127)
        input[i].Coefficient = std::numeric_limits<double>::denorm_min();
    ThreadPool<> pool(4);
    for (int variant = 0; variant < 5; ++variant)
    {
        std::vector<Operation> operations;
        for (int q = 1; q <= 9; ++q)
            operations.emplace_back(Kind::H, q);
        if (width > 64)
            operations.emplace_back(Kind::H, static_cast<int>(width - 1));
        for (size_t gate = 0; gate < 8; ++gate)
            operations.push_back(Operation::Rotation(Kind::RY, 0, .13 + .01 * gate));
        if (variant == 1)
        {
            operations[operations.size() - 3] = Operation(Kind::H, 0);
            operations[operations.size() - 6] = Operation::Rotation(Kind::RZ, 0, .57);
        }
        if (variant == 2)
        {
            auto projector = Operation(Kind::PROJ, 0);
            projector.coefficient = .5;
            operations[operations.size() - 4] = projector;
        }
        if (variant == 3)
            operations[operations.size() - 4] = Operation::Import(std::make_unique<Revive>());
        if (variant == 4)
        {
            operations[operations.size() - 2] = Operation::Rotation(Kind::RZ, 0, .57);
            operations[operations.size() - 3] = Operation::Rotation(Kind::RX, 0, -.31);
            operations[operations.size() - 5] = Operation(Kind::H, 0);
        }
        for (bool expansion : {false, true})
        {
            PauliDetail::Settings settings;
            settings.qubits = width;
            settings.dedup = 4;
            settings.trims = 3;
            settings.parallelThreshold = 4096;
            settings.batch = 2048;
            if (variant == 1)
                settings.coefficient = .05;
            PauliDetail::Workspace<T> reference, actual;
            reference.terms = input;
            actual.terms = input;
            if (!expansion && (variant == 0 || variant == 2 || variant == 4))
                for (size_t i = 0; i < input.size(); ++i)
                {
                    // Scalar block compilation reassociates exact-mode sums.
                    // Test a conditioned observable expansion against the
                    // original full operator propagation.
                    const double coefficient =
                        i % 127 == 0 ? std::numeric_limits<double>::denorm_min() : double(int(i % 11) - 5) * 1e-5;
                    reference.terms[i].Coefficient = actual.terms[i].Coefficient = coefficient;
                }
            const auto scalarInput = actual.terms;
            const double want = PauliDetail::Execute(reference, operations, settings, true);
            settings.pool = &pool;
            const double got = PauliDetail::Execute(actual, operations, settings, expansion);
            if (expansion || variant == 1 || variant == 3)
                Same(actual.terms, reference.terms);
            Require(std::abs(got - want) <= 1e-12 * std::max(1., std::abs(want)),
                    "scalar block fusion changed scalar result");
            if (!expansion && variant == 0)
            {
                ThreadPool<> oneWorker(1);
                auto oneSettings = settings;
                oneSettings.pool = &oneWorker;
                PauliDetail::Workspace<T> one;
                one.terms = scalarInput;
                const double deterministic = PauliDetail::Execute(one, operations, oneSettings, false);
                Same(one.terms, actual.terms);
                Require(std::memcmp(&got, &deterministic, sizeof(double)) == 0,
                        "scalar block fusion depends on worker count");
                if (width == 64)
                {
                    const auto concurrent = [&] {
                        PauliDetail::Workspace<T> local;
                        local.terms = scalarInput;
                        return PauliDetail::Execute(local, operations, settings, false);
                    };
                    auto first = std::async(std::launch::async, concurrent);
                    auto second = std::async(std::launch::async, concurrent);
                    Require(first.get() == got && second.get() == got,
                            "concurrent scalar block queries changed the result");
                }
            }
        }
    }
    const auto table = PauliDetail::ParameterizedTransfer(Kind::CRY, .713);
    const auto operation = Operation::Local(Kind::CRY, 0, 1, 0, table);
    auto large = Input<T>(width, 65536, 16384);
    auto expected = large;
    PauliDetail::Settings settings;
    settings.qubits = width;
    PauliDetail::ExecuteLocal(expected, operation, settings);
    settings.pool = &pool;
    PauliDetail::ExecuteLocal(large, operation, settings);
    Same(large, expected);
}

inline void FusionChecks()
{
    using Operation = PauliDetail::Operation;
    using Kind = OperationType;
    std::vector<Operation> operations{Operation::Rotation(Kind::RX, 0, .31),
                                      Operation(Kind::H, 0),
                                      Operation::Rotation(Kind::RZ, 0, -.17),
                                      Operation(Kind::S, 0),
                                      Operation::Rotation(Kind::RY, 0, -.73),
                                      Operation(Kind::K, 0),
                                      Operation(Kind::SXDG, 0),
                                      Operation::Rotation(Kind::RX, 0, std::acos(-1.) / 2),
                                      Operation(Kind::SDG, 0)};
    const auto table = PauliDetail::CompileScalarBlock(operations, 0, operations.size() - 1);
    for (unsigned input = 0; input < 4; ++input)
    {
        PauliStringXZWithCoefficient initial(1);
        initial.X[0] = input & 1;
        initial.Z[0] = (input >> 1) & 1;
        PauliStringStorage terms{initial};
        for (size_t next = operations.size(); next;)
        {
            auto operation = operations[--next].Legacy();
            PauliStringStorage branches;
            for (auto &term : terms)
                operation->Apply(term, branches);
            for (auto &term : branches)
                terms.push_back(std::move(term));
        }
        std::array<double, 4> expected{}, actual{};
        for (const auto &term : terms)
            expected[unsigned(term.X[0]) | (unsigned(term.Z[0]) << 1)] += term.Coefficient;
        for (size_t entry = table->offsets[input]; entry < table->offsets[input + 1]; ++entry)
            actual[table->entries[entry].pauli] += table->entries[entry].coefficient;
        for (unsigned output = 0; output < 4; ++output)
            Require(std::abs(expected[output] - actual[output]) < 2e-12,
                    "fused scalar block disagrees with independent virtual operators");
    }
}

inline void NonfiniteFusionChecks()
{
    using Term = PauliDetail::Term<1>;
    using Operation = PauliDetail::Operation;
    using Kind = OperationType;
    const std::vector<Operation> operations{Operation(Kind::H, 0), Operation::Rotation(Kind::RY, 0, .23),
                                            Operation::Rotation(Kind::RY, 0, .23)};
    ThreadPool<> pool(4);
    PauliDetail::Settings settings;
    settings.qubits = 1;
    settings.dedup = 4;
    settings.parallelThreshold = 4096;
    settings.batch = 2048;
    for (double coefficient : {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()})
        for (bool isolated : {false, true})
        {
            PauliDetail::Workspace<Term> reference, actual;
            reference.terms.assign(16384, Term{});
            for (size_t i = 0; i < reference.terms.size(); ++i)
            {
                reference.terms[i].X[0] = 1;
                reference.terms[i].Coefficient =
                    isolated && i + 1 < reference.terms.size() ? 1. / reference.terms.size() : coefficient;
            }
            actual.terms = reference.terms;
            settings.pool = nullptr;
            const double expected = PauliDetail::Execute(reference, operations, settings, true);
            settings.pool = &pool;
            const double result = PauliDetail::Execute(actual, operations, settings, false);
            const bool same = std::isnan(expected) ? std::isnan(result)
                              : std::isinf(expected)
                                  ? std::isinf(result) && std::signbit(expected) == std::signbit(result)
                                  : expected == result;
            Require(same, "scalar fusion changed nonfinite arithmetic classification");
            // The two RY expansions produce Inf-Inf in the eventual Z label;
            // legacy trimming removes that NaN and the surviving X contributes0.
            Require(expected == 0., "nonfinite legacy regression no longer exercises cancellation");
        }
}

inline void NonfiniteAngleFusionChecks()
{
    using Term = PauliDetail::Term<1>;
    using Operation = PauliDetail::Operation;
    using Kind = OperationType;
    size_t checks = 0;
    PauliDetail::Settings settings;
    settings.qubits = 15;
    for (Kind axis : {Kind::RX, Kind::RY, Kind::RZ})
        for (double angle : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()})
            for (unsigned separation : {0u, 1u, 2u})
            {
                std::vector<Operation> operations;
                for (int qubit = 1; qubit < 15; ++qubit)
                    operations.emplace_back(Kind::H, qubit);
                if (axis == Kind::RZ)
                    operations.emplace_back(Kind::H, 0);
                operations.push_back(Operation::Rotation(axis, 0, angle));
                // Finite inverse rotations may be fused before a nonfinite
                // rotation, even outside the 32-operation fusion window.
                if (separation == 1)
                    operations.emplace_back(Kind::H, 1);
                if (separation == 2)
                    for (unsigned index = 0; index < 34; ++index)
                        operations.emplace_back(Kind::H, 0);
                const auto inverseAxis = axis == Kind::RX ? Kind::RY : Kind::RX;
                operations.push_back(Operation::Rotation(inverseAxis, 0, -.3));
                operations.push_back(Operation::Rotation(inverseAxis, 0, .3));
                PauliDetail::Workspace<Term> reference, actual;
                reference.terms.resize(16384);
                for (size_t index = 0; index < reference.terms.size(); ++index)
                {
                    auto &term = reference.terms[index];
                    term.X[0] = axis == Kind::RZ ? 0 : 1;
                    term.Z[0] = (uint64_t(index) << 1) | (axis == Kind::RX ? 0 : 1);
                    term.Coefficient = 1. / reference.terms.size();
                }
                actual.terms = reference.terms;
                const double expected = PauliDetail::Execute(reference, operations, settings, true);
                const double result = PauliDetail::Execute(actual, operations, settings, false);
                Require(std::isnan(expected), "Nonfinite angle regression no longer exercises NaN propagation");
                Require(std::isnan(result), "Scalar fusion hid nonfinite rotation arithmetic");
                ++checks;
            }
    std::cout << "PASS " << checks << " nonfinite RX/RY/RZ fusion comparisons across block boundaries\n";
}

template <class Term> inline void ContractionOverflowChecks(size_t width)
{
    using Operation = PauliDetail::Operation;
    using Kind = OperationType;
    PauliDetail::Settings settings;
    settings.qubits = width;
    size_t checks = 0;
    for (bool scheduledDedup : {false, true})
        for (size_t count : {size_t(63), size_t(64), size_t(65), size_t(128)})
            for (double coefficient : {1e308, -1e308, std::numeric_limits<double>::max()})
                for (bool cancelling : {false, true})
                {
                    settings.dedup = scheduledDedup ? 2 : std::numeric_limits<int>::max();
                    std::vector<Operation> operations{Operation::Rotation(Kind::RY, 0, 1e-308)};
                    if (scheduledDedup)
                    {
                        operations.emplace_back(Kind::H, 0);
                        operations.emplace_back(Kind::H, 0);
                    }
                    PauliDetail::Workspace<Term> reference, actual;
                    for (size_t index = 0; index < count; ++index)
                    {
                        Term term(width);
                        term.X[0] = 1;
                        term.Z[0] = uint64_t(index) << 1;
                        // Exercise the packed/dynamic high words as finished labels.
                        if (width > 64)
                            term.Z[(width - 1) / 64] |= uint64_t(index & 1) << ((width - 1) % 64);
                        term.Coefficient = cancelling && (index & 1) ? -coefficient : coefficient;
                        reference.terms.push_back(std::move(term));
                    }
                    actual.terms = reference.terms;
                    const double expected = PauliDetail::Execute(reference, operations, settings, true);
                    const double result = PauliDetail::Execute(actual, operations, settings, false);
                    Require(std::isfinite(expected) && std::isfinite(result),
                            "Contraction overflowed finite arithmetic");
                    Require(std::abs(expected - result) <= 2e-12 * std::max(1., std::abs(expected)),
                            "Contraction changed extreme finite coefficients");
                    if (!cancelling)
                        Require(std::abs(expected - count * (coefficient * std::sin(1e-308))) < 2e-12,
                                "Contraction overflow reference disagrees with analytic expectation");
                    ++checks;
                }
    std::cout << "PASS " << checks << " contraction overflow/boundary comparisons at width " << width << '\n';
}

inline void NonfiniteContractionChecks()
{
    using Term = PauliDetail::Term<1>;
    using Operation = PauliDetail::Operation;
    using Kind = OperationType;
    PauliDetail::Settings settings;
    settings.qubits = 7;
    size_t checks = 0;
    for (double angle : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()})
        for (size_t count : {size_t(63), size_t(64), size_t(65), size_t(128)})
            for (bool delayed : {false, true})
            {
                std::vector<Operation> operations{Operation::Rotation(Kind::RY, 0, angle)};
                // A finite scalar contraction can precede a nonfinite rotation
                // beyond the fusion window, so inspect the entire query.
                if (delayed)
                    for (unsigned index = 0; index < 34; ++index)
                        operations.emplace_back(Kind::H, 0);
                PauliDetail::Workspace<Term> reference, actual;
                for (size_t index = 0; index < count; ++index)
                {
                    Term term;
                    term.X[0] = 1;
                    term.Z[0] = uint64_t(index % 64) << 1;
                    term.Coefficient = index & 1 ? -1. : 1.;
                    reference.terms.push_back(term);
                }
                actual.terms = reference.terms;
                const double expected = PauliDetail::Execute(reference, operations, settings, true);
                const double result = PauliDetail::Execute(actual, operations, settings, false);
                Require(std::isnan(expected) && std::isnan(result), "Contraction hid nonfinite rotation arithmetic");
                ++checks;
            }
    std::cout << "PASS " << checks << " nonfinite contraction comparisons below the fusion threshold\n";
}

inline void Run()
{
    FusionChecks();
    NonfiniteFusionChecks();
    NonfiniteAngleFusionChecks();
    ContractionOverflowChecks<PauliDetail::Term<1>>(8);
    ContractionOverflowChecks<PauliDetail::Term<8>>(257);
    ContractionOverflowChecks<PauliDetail::Term<0>>(513);
    NonfiniteContractionChecks();
    DedupChecks<PauliDetail::Term<1>>(64);
    DedupChecks<PauliDetail::Term<8>>(257);
    DedupChecks<PauliDetail::Term<0>>(513);
    AlternatingDedupChecks<PauliDetail::Term<1>>(64);
    AlternatingDedupChecks<PauliDetail::Term<8>>(257);
    AlternatingDedupChecks<PauliDetail::Term<0>>(513);
    SmallRetainedDedupChecks<PauliDetail::Term<1>>(64);
    SmallRetainedDedupChecks<PauliDetail::Term<8>>(257);
    SmallRetainedDedupChecks<PauliDetail::Term<0>>(513);
    std::cout
        << "PASS runtime dedup reference, bitwise accumulation order, wide keys, scratch reuse and capacity hints\n";
    ScalarBlockChecks<PauliDetail::Term<1>>(64);
    ScalarBlockChecks<PauliDetail::Term<8>>(257);
    ScalarBlockChecks<PauliDetail::Term<0>>(513);
    std::cout << "PASS scalar block adjoints/references, worker determinism/concurrency, expansion/cutoff/custom "
                 "preservation and native counts\n";
}
} // namespace PauliRuntimeEngineRegression
} // namespace QC
