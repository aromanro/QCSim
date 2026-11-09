#pragma once

#include "PathIntegralLightconeRegressionTests.h"
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

namespace QC
{
namespace PathIntegralRngRegression
{
using namespace PathIntegralLightconeRegression;
using Map = PathIntegral::AmplitudeMap;

inline Map State(size_t count = 8)
{
    Map state;
    const size_t width = count > 8 ? 8 : 3;
    for (size_t index = 0; index < count; ++index)
        state[PathIntegral::FastVectorBool(Bits(index, width))] = double(1 + index % 9);
    return state;
}

inline void SameSamples(Simulator &left, Simulator &right, size_t &checks)
{
    for (size_t sample = 0; sample < 128; ++sample)
    {
        Require(left.MeasureNoCollapse() == right.MeasureNoCollapse(), "Path measurement RNG sequence changed");
        ++checks;
    }
}

inline PathIntegral::FastVectorBool ExpectedState(const Map &state, double target)
{
    double cumulative = 0.;
    PathIntegral::FastVectorBool last;
    for (const auto &entry : state)
    {
        last = PathIntegral::FastVectorBool(entry.first);
        cumulative += std::norm(entry.second);
        if (target < cumulative)
            return last;
    }
    return last;
}

inline void SeededReference(size_t &checks)
{
    // Integer amplitudes make both the small and block CDF sums exact. Compare
    // the public draw sequence with the same standard engine and distribution
    // used before lazy initialization, including seeds with high bits set.
    for (uint64_t seed : {uint64_t(0), uint64_t(1), uint64_t(1) << 32, ~uint64_t(0)})
        for (size_t count : {size_t(8), size_t(129)})
        {
            const auto state = State(count);
            Simulator simulator;
            simulator.GetAmplitudes() = state;
            simulator.SetSeed(seed);
            std::seed_seq sequence{uint32_t(seed & 0xffffffff), uint32_t(seed >> 32)};
            std::mt19937_64 engine(sequence);
            std::uniform_real_distribution<double> uniform(0., 1.);
            double total = 0.;
            for (const auto &entry : state)
                total += std::norm(entry.second);
            for (size_t sample = 0; sample < 128; ++sample)
            {
                const double raw = uniform(engine) * total;
                const double target = raw < total ? raw : std::nextafter(total, 0.);
                Require(simulator.MeasureNoCollapse() == ExpectedState(state, target),
                        "Seeded path measurement differs from the original engine sequence");
                ++checks;
            }

            // Reseeding an already used engine must restart the same stream.
            Simulator fresh;
            fresh.GetAmplitudes() = state;
            fresh.SetSeed(seed);
            simulator.SetSeed(seed);
            SameSamples(simulator, fresh, checks);
        }

    const auto state = State();
    Simulator simulator;
    simulator.SetSeed(0x123456789abcdef0ULL);
    std::seed_seq sequence{uint32_t(0x9abcdef0), uint32_t(0x12345678)};
    std::mt19937_64 engine(sequence);
    std::uniform_real_distribution<double> uniform(0., 1.);
    for (size_t sample = 0; sample < 128; ++sample)
    {
        simulator.GetAmplitudes() = state;
        const size_t qubit = sample % 3;
        double total = 0., zero = 0.;
        for (const auto &entry : state)
        {
            const double mass = std::norm(entry.second);
            total += mass;
            if (!entry.first.get(qubit))
                zero += mass;
        }
        const bool expected = uniform(engine) * total >= zero;
        Require(simulator.MeasureQubit(qubit) == expected, "Seeded collapsing path measurement changed");
        for (const auto &entry : simulator.GetAmplitudes())
            Require(entry.first.get(qubit) == expected, "Path measurement retained the wrong sector");
        ++checks;
    }
}

inline void Lifecycle(size_t &checks)
{
    static_assert(std::is_nothrow_move_constructible<Simulator>::value, "Path move construction must remain noexcept");
    static_assert(std::is_nothrow_move_assignable<Simulator>::value, "Path move assignment must remain noexcept");
    for (bool initialized : {false, true})
    {
        Simulator source;
        source.GetAmplitudes() = State();
        if (initialized)
            source.MeasureNoCollapse();
        auto copied = source;
        Simulator assigned;
        assigned.SetSeed(17);
        assigned = source;
        auto control = source;
        SameSamples(copied, control, checks);
        control = source;
        SameSamples(assigned, control, checks);

        control = source;
        auto moved = std::move(source);
        SameSamples(moved, control, checks);
        control = moved;
        Simulator moveAssigned;
        moveAssigned.SetSeed(19);
        moveAssigned = std::move(moved);
        SameSamples(moveAssigned, control, checks);

        control = moveAssigned;
        moveAssigned.Reset();
        moveAssigned.GetAmplitudes() = State();
        SameSamples(moveAssigned, control, checks);
    }

    Simulator source;
    source.SetTrimValue(0.);
    const std::vector<Gate> circuit{{MatrixFor(0, 0.), 0}, {MatrixFor(0, 0.), 1}};
    source.PropagateAll(circuit, Bits(0, 3));
    source.SetSeed(123);
    source.MeasureNoCollapse();
    auto cloned = source.Clone();
    Require(cloned.GetForwardGateCount() == source.GetForwardGateCount(), "Path clone lost the recorded circuit");
    Require(cloned.GetAmplitudes().size() == source.GetAmplitudes().size(), "Path clone lost the stored state");
    auto cloneCopy = cloned;
    SameSamples(cloned, cloneCopy, checks);
    cloned.SetSeed(456);
    source.SetSeed(456);
    SameSamples(cloned, source, checks);
}

template <class Exception, class Action> inline void Throws(Action action, size_t &checks)
{
    bool caught = false;
    try
    {
        action();
    }
    catch (const Exception &)
    {
        caught = true;
    }
    Require(caught, "Invalid path measurement did not preserve its exception");
    ++checks;
}

inline void InvalidMeasurements(size_t &checks)
{
    const auto state = State();
    for (bool initialized : {false, true})
    {
        Simulator simulator;
        simulator.GetAmplitudes() = state;
        if (initialized)
        {
            simulator.SetSeed(1234);
            simulator.MeasureNoCollapse();
        }
        auto control = simulator;
        simulator.GetAmplitudes().clear();
        Throws<std::domain_error>([&] { simulator.MeasureNoCollapse(); }, checks);
        Throws<std::domain_error>([&] { simulator.MeasureQubit(0); }, checks);
        simulator.GetAmplitudes() = state;
        Throws<std::out_of_range>([&] { simulator.MeasureQubit(3); }, checks);
        for (double value : {0., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        {
            for (auto entry : simulator.GetAmplitudes())
                entry.second = value;
            Throws<std::domain_error>([&] { simulator.MeasureNoCollapse(); }, checks);
            Throws<std::domain_error>([&] { simulator.MeasureQubit(0); }, checks);
        }
        simulator.GetAmplitudes() = state;
        SameSamples(simulator, control, checks);
    }
}

// States whose total is below 2^-900, or whose squares underflow or overflow,
// are weighed from amplitudes scaled by an exact power of two. Both
// measurements must follow the law of the amplitude ratios, draw one uniform
// per call, and leave a collapsed state normalized.
inline void ExtremeTotals(size_t &checks)
{
    const double tiny = std::sqrt(3. * std::numeric_limits<double>::denorm_min());
    struct Case
    {
        std::complex<double> zero, one;
        double ones;
    };
    const Case cases[] = {
        // Squares below the smallest double.
        {{std::sqrt(.7) * 0x1p-540, 0.}, {0., std::sqrt(.3) * 0x1p-540}, .3},
        // Squares above the largest double.
        {{std::sqrt(.7) * 0x1p600, 0.}, {0., std::sqrt(.3) * 0x1p600}, .3},
        // Equal subnormal weights: an unscaled target rounds to the subnormal
        // grid and picks the first row with probability 5/12.
        {{tiny, 0.}, {tiny, 0.}, .5},
    };
    for (const auto &c : cases)
    {
        Map state;
        state[PathIntegral::FastVectorBool(Bits(0, 1))] = c.zero;
        state[PathIntegral::FastVectorBool(Bits(1, 1))] = c.one;
        Simulator simulator, control;
        simulator.SetSeed(2024);
        control.SetSeed(2024);
        simulator.GetAmplitudes() = state;
        control.GetAmplitudes() = State();
        constexpr size_t draws = 20000, collapses = 4000;
        size_t ones = 0;
        for (size_t draw = 0; draw < draws; ++draw)
        {
            ones += simulator.MeasureNoCollapse().get(0);
            control.MeasureNoCollapse();
        }
        Require(std::abs(double(ones) - c.ones * draws) < 6. * std::sqrt(c.ones * (1. - c.ones) * draws),
                "Path measurement of an extreme total lost its law");
        ++checks;
        ones = 0;
        for (size_t draw = 0; draw < collapses; ++draw)
        {
            simulator.GetAmplitudes() = state;
            control.GetAmplitudes() = State();
            const bool one = simulator.MeasureQubit(0);
            control.MeasureQubit(0);
            ones += one;
            double mass = 0.;
            for (const auto &entry : simulator.GetAmplitudes())
            {
                Require(entry.first.get(0) == one, "Collapsing path measurement retained the wrong sector");
                mass += std::norm(entry.second);
            }
            Require(std::abs(mass - 1.) < 1e-12, "Collapsing path measurement of an extreme total is not normalized");
        }
        Require(std::abs(double(ones) - c.ones * collapses) < 6. * std::sqrt(c.ones * (1. - c.ones) * collapses),
                "Collapsing path measurement of an extreme total lost its law");
        ++checks;
        // One uniform per call: both engines must still agree.
        simulator.GetAmplitudes() = State();
        control.GetAmplitudes() = State();
        SameSamples(simulator, control, checks);
    }
}

inline void Run()
{
    size_t checks = 0;
    SeededReference(checks);
    Lifecycle(checks);
    InvalidMeasurements(checks);
    ExtremeTotals(checks);
    std::cout << "PASS " << checks
              << " path RNG reference, copy/move/reset/clone, invalid-measurement and extreme-total comparisons\n";
}
} // namespace PathIntegralRngRegression
} // namespace QC
