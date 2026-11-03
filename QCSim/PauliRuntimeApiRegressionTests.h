#pragma once

#include "PauliPropagator.h"
#include <atomic>
#include <future>
#include <iostream>
#include <stdexcept>

namespace QC
{
namespace PauliRuntimeApiRegression
{
inline void Require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

inline void Close(double actual, double expected, const char *message)
{
    Require(std::isfinite(actual) && std::isfinite(expected) && std::abs(actual - expected) < 3e-10, message);
}

inline void ContractionChecks()
{
    size_t checks = 0;
    for (double coefficient : {1e308, -1e308, std::numeric_limits<double>::max()})
    {
        PauliPropagator sim;
        sim.SetNrQubits(7);
        sim.ApplyRY(0, 1e-308);
        PauliStringStorage input;
        for (size_t index = 0; index < 64; ++index)
        {
            PauliStringXZWithCoefficient term(7);
            term.X[0] = true;
            for (size_t qubit = 1; qubit < 7; ++qubit)
                term.Z[qubit] = ((index >> (qubit - 1)) & 1) != 0;
            term.Coefficient = coefficient;
            input.push_back(std::move(term));
        }
        double expected = 0.;
        for (const auto &term : input)
            expected += sim.ExpectationValue(term);
        Close(sim.ExpectationValue(input), expected, "Public batched expectation overflowed during contraction");
        ++checks;
    }
    for (double angle : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()})
    {
        PauliPropagator sim;
        sim.SetNrQubits(7);
        sim.ApplyRY(0, angle);
        for (int qubit = 1; qubit < 7; ++qubit)
        {
            sim.ApplyRY(qubit, -.3);
            sim.ApplyRY(qubit, .3);
        }
        PauliStringStorage output;
        const double reference = sim.ExpectationValue("XXXXXXX", output);
        Require(std::isnan(reference) && std::isnan(sim.ExpectationValue("XXXXXXX")),
                "Public scalar contraction hid a nonfinite angle");
        ++checks;
    }
    std::cout << "PASS " << checks << " public API extreme/nonfinite scalar contraction comparisons\n";
}

inline void ReferenceFrom(const PauliPropagator &source, PauliPropagator &reference)
{
    reference.SetNrQubits(source.GetNrQubits());
    reference.SetOperations(source.GetOperations());
    reference.SetCoefficientThreshold(source.GetCoefficientThreshold());
    reference.SetPauliWeightThreshold(source.GetPauliWeightThreshold());
    reference.SetStepsBetweenTrims(source.StepsBetweenTrims());
    reference.SetStepsBetweenDeduplication(source.StepsBetweenDeduplication());
    reference.SetSamplingCacheMaxNodes(0);
}

inline void SamplingMatches(PauliPropagator &source, const std::vector<int> &qubits)
{
    PauliPropagator reference;
    ReferenceFrom(source, reference);
    source.SetSeed(239);
    reference.SetSeed(239);
    for (size_t repeat = 0; repeat < 3; ++repeat)
        Require(source.SampleCounts(qubits, 73) == reference.SampleCounts(qubits, 73),
                "Persistent sampling cache changed seeded results");
    source.SetSeed(912);
    reference.SetSeed(912);
    for (size_t repeat = 0; repeat < 20; ++repeat)
        Require(source.Sample(qubits) == reference.Sample(qubits), "Persistent single-shot cache changed results");
}

inline void ProbabilitiesMatch(PauliPropagator &source)
{
    PauliPropagator reference;
    ReferenceFrom(source, reference);
    const auto count = source.GetOperations().size();
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
        double total = 0.;
        for (size_t outcome = 0; outcome < (size_t{1} << source.GetNrQubits()); ++outcome)
        {
            const double actual = source.Probability(outcome);
            Close(actual, reference.Probability(outcome), "Moment probability/reference mismatch");
            total += actual;
        }
        if (source.GetCoefficientThreshold() == 0. &&
            source.GetPauliWeightThreshold() >= static_cast<size_t>(source.GetNrQubits()))
            Close(total, 1., "Exact cached probabilities lost normalization");
    }
    Require(count == source.GetOperations().size(), "Probability cache mutated the recorded circuit");
}

class NestedScale : public OperatorX
{
  public:
    explicit NestedScale(const PauliPropagator &nested) : OperatorX(0), nested(&nested)
    {
    }

    void Apply(PauliStringXZWithCoefficient &term, PauliStringStorage &) const override
    {
        term.Coefficient *= nested->ExpectationValue("Z");
    }

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<NestedScale>(*this);
    }

  private:
    const PauliPropagator *nested;
};

class ThrowOnce : public OperatorX
{
  public:
    ThrowOnce() : OperatorX(0)
    {
    }

    void Apply(PauliStringXZWithCoefficient &, PauliStringStorage &) const override
    {
        if (first.exchange(false))
            throw std::runtime_error("Expected custom query failure");
    }

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<ThrowOnce>();
    }

  private:
    mutable std::atomic<bool> first{true};
};

class MutableEstimate : public OperatorX
{
  public:
    explicit MutableEstimate(const double &value) : OperatorX(0), value(&value)
    {
    }

    void Apply(PauliStringXZWithCoefficient &term, PauliStringStorage &) const override
    {
        term.X[0] = term.Z[0] = false;
        term.Coefficient *= *value;
    }

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<MutableEstimate>(*this);
    }

  private:
    const double *value;
};

inline void ValidationBoundaries()
{
    constexpr uint64_t seed = 91;
    for (size_t count :
         {size_t(0), size_t(1), size_t(6), size_t(8), size_t(9), size_t(13), size_t(16), size_t(17), size_t(64)})
    {
        PauliPropagator p;
        p.SetNrQubits(64);
        p.SetSeed(seed);
        std::seed_seq sequence{uint32_t(seed), uint32_t(seed >> 32)};
        std::mt19937_64 reference(sequence);
        std::uniform_real_distribution<double> uniform(0., 1.);
        std::vector<int> qubits(count);
        std::iota(qubits.rbegin(), qubits.rend(), 0);
        Require(p.SampleCounts(qubits, 0).empty(), "Valid zero-shot subset was rejected");
        if (count && count <= 13)
        {
            const auto counts = p.SampleCounts(qubits, 1);
            Require(counts.size() == 1 && counts.begin()->second == 1 &&
                        std::none_of(counts.begin()->first.begin(), counts.begin()->first.end(),
                                     [](bool bit) { return bit; }),
                    "Valid deterministic sampling changed at validation boundary");
            for (size_t draw = 0; draw < count; ++draw)
                (void)uniform(reference);
        }
        if (count > 1)
        {
            auto invalid = qubits;
            invalid[1] = invalid[0];
            invalid.back() = -1;
            bool duplicate = false;
            try
            {
                p.SampleCounts(invalid, 0);
            }
            catch (const std::invalid_argument &)
            {
                duplicate = true;
            }
            Require(duplicate, "Duplicate validation did not precede later range error");
            invalid = qubits;
            invalid[0] = -1;
            invalid.back() = invalid[1];
            bool range = false;
            try
            {
                p.SampleCounts(invalid, 0);
            }
            catch (const std::out_of_range &)
            {
                range = true;
            }
            Require(range, "Range validation did not precede later duplicate");
        }
        if (count)
        {
            auto invalid = qubits;
            invalid.back() = 64;
            bool range = false;
            try
            {
                p.SampleCounts(invalid, 0);
            }
            catch (const std::out_of_range &)
            {
                range = true;
            }
            Require(range, "Upper-bound readout qubit was accepted");
        }
        // Validation, including failures and zero shots, must consume no RNG.
        // A later H measurement independently reveals the expected next draws.
        p.ApplyH(0);
        for (size_t draw = 0; draw < 32; ++draw)
            Require(p.Sample({0})[0] == (uniform(reference) < .5), "Validation altered sampling RNG progression");
    }
    std::cout << "PASS <=8/9/13/16/17/large sampling validation, exception order and RNG progression\n";
}

inline void Run()
{
    ContractionChecks();
    ValidationBoundaries();
    for (int width : {4, 65, 129, 256, 257, 511, 512, 513})
    {
        PauliPropagator p;
        p.SetNrQubits(width);
        p.SetStepsBetweenDeduplication(3);
        const std::array<int, 3> active{0, width / 2, width - 1};
        const std::array<double, 3> angles{.21, .35, .56};
        PauliStringXZWithCoefficient observable(width);
        double expected = 1.;
        for (size_t i = 0; i < active.size(); ++i)
        {
            p.ApplyRY(active[i], angles[i]);
            observable.X[active[i]] = true;
            expected *= std::sin(angles[i]);
        }
        for (size_t repeat = 0; repeat < 10; ++repeat)
        {
            Close(p.ExpectationValue(observable), expected, "Reused packed-width workspace changed expectation");
            PauliStringStorage expansion;
            Close(p.ExpectationValue(observable, expansion), expected, "Reused output workspace changed expectation");
            Require(!expansion.empty(), "Output expansion was discarded");
        }
        SamplingMatches(p, {active[2], active[0], active[1]});
        p.SetCoefficientThreshold(.015);
        p.SetPauliWeightThreshold(2);
        SamplingMatches(p, {active[2], active[0], active[1]});
    }

    PauliPropagator p;
    p.SetNrQubits(5);
    p.SetStepsBetweenDeduplication(3);
    for (int q = 0; q < 5; ++q)
    {
        p.ApplyRY(q, .13 + .07 * q);
        if (q)
            p.ApplyCX(q - 1, q);
    }
    ProbabilitiesMatch(p);
    {
        PauliPropagator reference;
        ReferenceFrom(p, reference);
        p.SetSamplingCacheMaxNodes(0);
        p.SetSamplingCacheMaxNodes(8192);
        p.SetSeed(713);
        reference.SetSeed(713);
        // Cover both initial normalized-projector evaluation and promotion to
        // full moments without consuming any of the later sampling draws.
        for (size_t outcome = 0; outcome < 20; ++outcome)
            p.Probability(outcome);
        Require(p.SampleCounts({0, 1, 2, 3, 4}, 127) == reference.SampleCounts({0, 1, 2, 3, 4}, 127),
                "Probability evaluation consumed sampling RNG draws");
    }
    SamplingMatches(p, {4, 0, 3, 1});
    SamplingMatches(p, {0, 1, 3, 4});
    p.SaveState();
    p.ApplyH(2);
    ProbabilitiesMatch(p);
    p.RestoreState();
    ProbabilitiesMatch(p);
    p.SetCoefficientThreshold(.012);
    ProbabilitiesMatch(p);
    p.SetCoefficientThreshold(0.);
    p.SetPauliWeightThreshold(2);
    SamplingMatches(p, {4, 0, 3, 1});
    p.SetPauliWeightThreshold(size_t(-1));
    p.SetStepsBetweenTrims(2);
    p.SetStepsBetweenDeduplication(2);
    SamplingMatches(p, {4, 0, 3, 1});
    PauliPropagator clone;
    clone.SetNrQubits(5);
    clone.ShareOperationsFrom(p);
    ProbabilitiesMatch(clone);
    p.ApplyX(0);
    ProbabilitiesMatch(p);
    ProbabilitiesMatch(clone);
    p.ClearOperations();
    Close(p.Probability(0), 1., "Cleared program reused stale moment probabilities");
    p.ApplyX(4);
    Close(p.Probability(16), 1., "Appended program reused stale moment probabilities");

    PauliPropagator nested;
    nested.SetNrQubits(1);
    nested.ApplyRY(0, .43);
    PauliPropagator outer;
    outer.SetNrQubits(1);
    std::vector<std::unique_ptr<Operator>> operators;
    operators.push_back(std::make_unique<NestedScale>(nested));
    outer.SetOperations(std::move(operators));
    for (size_t repeat = 0; repeat < 10; ++repeat)
        Close(outer.ExpectationValue("Z"), std::cos(.43), "Reentrant query reused an active workspace");
    operators.clear();
    operators.push_back(std::make_unique<ThrowOnce>());
    outer.SetOperations(std::move(operators));
    bool threw = false;
    try
    {
        outer.ExpectationValue("Z");
    }
    catch (const std::runtime_error &)
    {
        threw = true;
    }
    Require(threw, "Expected custom failure was lost");
    Close(outer.ExpectationValue("Z"), 1., "Exception left a leased workspace busy or stale");

    double estimate = 1.;
    operators.clear();
    operators.push_back(std::make_unique<MutableEstimate>(estimate));
    outer.SetOperations(std::move(operators));
    outer.SetSeed(77);
    Require(outer.SampleCounts({0}, 11).at({false}) == 11, "Custom positive estimate");
    estimate = -1.;
    outer.SetSeed(77);
    Require(outer.SampleCounts({0}, 11).at({true}) == 11, "Custom state was cached across sampling calls");

    std::vector<std::future<double>> queries;
    for (size_t worker = 0; worker < 4; ++worker)
        queries.push_back(std::async(std::launch::async, [&] {
            double value = 0.;
            for (size_t repeat = 0; repeat < 50; ++repeat)
                value += nested.ExpectationValue("Z");
            return value;
        }));
    for (auto &query : queries)
        Close(query.get(), 50 * std::cos(.43), "Thread-local workspaces leaked across concurrent queries");
    std::cout << "PASS runtime API workspace leases, reentrancy/exceptions, 256/257/511/512/513 boundaries, "
                 "moment probabilities and persistent sampling invalidation\n";
}
} // namespace PauliRuntimeApiRegression
} // namespace QC
