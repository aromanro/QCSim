#pragma once

#include "PathIntegral.h"
#include <iostream>
#include <limits>
#include <random>

namespace QC
{
namespace PathExpectationRegression
{
using Simulator = PathIntegral::PathIntegralSimulator;
using Map = PathIntegral::AmplitudeMap;
using Bits = PathIntegral::FastVectorBool;
using Complex = std::complex<double>;

inline double Reference(Simulator &sim, const std::string &pauli)
{
    const auto &input = sim.GetAmplitudes();
    auto transformed = input;
    for (size_t qubit = 0; qubit < pauli.size(); ++qubit)
    {
        const char c = pauli[qubit];
        if (c == 'X' || c == 'x')
            sim.PropagateStep(Gates::AppliedGate<>(Gates::PauliXGate<>().getRawOperatorMatrix(), qubit), transformed);
        else if (c == 'Y' || c == 'y')
            sim.PropagateStep(Gates::AppliedGate<>(Gates::PauliYGate<>().getRawOperatorMatrix(), qubit), transformed);
        else if (c == 'Z' || c == 'z')
            sim.PropagateStep(Gates::AppliedGate<>(Gates::PauliZGate<>().getRawOperatorMatrix(), qubit), transformed);
    }
    double result = 0.;
    for (const auto &entry : transformed)
    {
        const auto original = input.find(entry.first);
        if (original != input.end())
            result += std::real(std::conj(original->second) * entry.second);
    }
    return result;
}

inline void Close(double actual, double expected)
{
    if ((std::isnan(actual) && std::isnan(expected)) || actual == expected)
        return;
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual - expected) > 2e-11)
        throw std::runtime_error("Path integral direct expectation mismatch");
}

inline std::vector<size_t> Active(size_t width)
{
    std::vector<size_t> sites{0, 1, width / 2, width / 2 + 1, width - 2, width - 1};
    std::sort(sites.begin(), sites.end());
    sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
    return sites;
}

inline Map Initial(size_t width)
{
    Map result;
    const auto active = Active(width);
    for (size_t state = 0; state < (size_t{1} << active.size()); ++state)
    {
        Bits bits(width);
        for (size_t bit = 0; bit < active.size(); ++bit)
            bits.set(active[bit], ((state >> bit) & 1) != 0);
        if (state % 7 != 6)
            result[bits] = Complex((int(state % 7) - 3) * .13, (int(state % 5) - 2) * .17);
    }
    return result;
}

inline double DenseReference(const Map &input, const std::string &pauli, double epsilon)
{
    const size_t width = input.QubitCount();
    std::vector<Complex> amplitudes(size_t{1} << width);
    for (const auto &entry : input)
    {
        size_t state = 0;
        for (size_t q = 0; q < width; ++q)
            state |= size_t(entry.first.get(q)) << q;
        amplitudes[state] = entry.second;
    }
    bool nonidentity = false;
    for (char c : pauli)
        nonidentity = nonidentity || c == 'X' || c == 'Y' || c == 'Z' || c == 'x' || c == 'y' || c == 'z';
    if (nonidentity && !(1. > epsilon))
        return 0.;
    double result = 0.;
    for (size_t state = 0; state < amplitudes.size(); ++state)
    {
        if (nonidentity && std::norm(amplitudes[state]) < epsilon)
            continue;
        size_t destination = state;
        Complex phase = 1.;
        for (size_t q = 0; q < pauli.size(); ++q)
        {
            const char c = pauli[q];
            const bool one = ((state >> q) & 1) != 0;
            if (c == 'X' || c == 'x' || c == 'Y' || c == 'y')
                destination ^= size_t{1} << q;
            if (c == 'Y' || c == 'y')
                phase *= one ? Complex(0., -1.) : Complex(0., 1.);
            else if (c == 'Z' || c == 'z')
                phase *= one ? -1. : 1.;
        }
        result += std::real(std::conj(amplitudes[destination]) * phase * amplitudes[state]);
    }
    return result;
}

inline void Run()
{
    size_t checks = 0;
    std::mt19937 random(816317);
    for (size_t width : {size_t(3), size_t(6), size_t(65), size_t(129), size_t(513), size_t(1024)})
    {
        Simulator sim;
        sim.GetAmplitudes() = Initial(width);
        const auto original = sim.GetAmplitudes();
        const auto active = Active(width);
        for (double epsilon : {0., 1e-15, .04, .3, 1.})
        {
            sim.SetTrimValue(epsilon);
            for (size_t trial = 0; trial < 64; ++trial)
            {
                std::string pauli(width, 'I');
                for (size_t site : active)
                    pauli[site] = "IXYZixyz?"[random() % 9];
                const auto expected = Reference(sim, pauli);
                Close(sim.ExpectationValue(pauli), expected);
                if (width <= 6)
                    Close(sim.ExpectationValue(pauli), DenseReference(original, pauli, epsilon));
                ++checks;
            }
            Close(sim.ExpectationValue(""), Reference(sim, ""));
            Close(sim.ExpectationValue(std::string(width, 'I')), Reference(sim, std::string(width, 'I')));
        }
        if (sim.GetAmplitudes().size() != original.size())
            throw std::runtime_error("Direct expectation changed sparse support");
        for (const auto &entry : original)
            if (sim.GetAmplitudes().at(entry.first) != entry.second)
                throw std::runtime_error("Direct expectation changed amplitudes");
    }
    Simulator unusual;
    unusual.GetAmplitudes() = Initial(3);
    for (double epsilon :
         {-1e-6, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), .9})
    {
        unusual.SetTrimValue(epsilon);
        for (const std::string p : {"III", "XYZ", "yyy", "XIIIX"})
        {
            bool referenceThrows = false, actualThrows = false;
            double reference = 0., actual = 0.;
            try
            {
                reference = Reference(unusual, p);
            }
            catch (const std::out_of_range &)
            {
                referenceThrows = true;
            }
            try
            {
                actual = unusual.ExpectationValue(p);
            }
            catch (const std::out_of_range &)
            {
                actualThrows = true;
            }
            if (referenceThrows != actualThrows)
                throw std::runtime_error("Direct expectation changed invalid-qubit behavior");
            if (!referenceThrows)
                Close(actual, reference);
            ++checks;
        }
    }
    for (double nonfinite : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
        unusual.SetTrimValue(0.);
        unusual.GetAmplitudes().begin()->second = Complex(nonfinite, 0.);
        for (const std::string p : {"III", "XXX", "YZZ"})
        {
            Close(unusual.ExpectationValue(p), Reference(unusual, p));
            ++checks;
        }
    }
    unusual.Reset();
    Close(unusual.ExpectationValue("XYZ"), 0.);
    std::cout << "PASS " << checks << " direct path expectation reference comparisons\n";
}
} // namespace PathExpectationRegression
} // namespace QC
