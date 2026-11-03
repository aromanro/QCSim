#pragma once

#include "PathIntegralStorage.h"
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

#define QC_PATH_INTEGRAL_DIRECT_EXPECTATION 1

namespace QC
{
namespace PathIntegral
{
namespace Detail
{
// Preserve the old gate-by-gate semantics for unusual thresholds, nonfinite
// amplitudes, and overlong strings. Normal finite requests use the mask kernel.
inline double PauliExpectationLegacy(const AmplitudeMap &input, double epsilon, const std::string &pauli)
{
    auto current = input;
    const std::complex<double> imaginary(0., 1.);
    for (size_t qubit = 0; qubit < pauli.size(); ++qubit)
    {
        const char symbol = pauli[qubit];
        const bool x = symbol == 'X' || symbol == 'x';
        const bool y = symbol == 'Y' || symbol == 'y';
        const bool z = symbol == 'Z' || symbol == 'z';
        if (!x && !y && !z)
            continue;
        if (current.empty())
            break;
        if (qubit >= current.QubitCount())
            throw std::out_of_range("Path integral gate qubit is outside the register");
        AmplitudeMap next;
        for (const auto &entry : current)
        {
            if (std::norm(entry.second) < epsilon)
                continue;
            const unsigned column = unsigned(entry.first.get(qubit));
            for (unsigned row = 0; row < 2; ++row)
            {
                const std::complex<double> coefficient =
                    z ? (row == column ? (column ? -1. : 1.) : 0.)
                      : (row != column ? (y ? (column ? -imaginary : imaginary) : 1.) : 0.);
                if (!(std::norm(coefficient) > epsilon))
                    continue;
                auto state = static_cast<FastVectorBool>(entry.first);
                state.set(qubit, row != 0);
                next[state] += coefficient * entry.second;
            }
        }
        current = std::move(next);
    }
    double result = 0.;
    for (const auto &entry : current)
    {
        const auto original = input.find(entry.first);
        if (original != input.end())
            result += std::real(std::conj(original->second) * entry.second);
    }
    return result;
}

inline unsigned PauliParity(uint64_t word) noexcept
{
    word ^= word >> 32;
    word ^= word >> 16;
    word ^= word >> 8;
    word ^= word >> 4;
    return unsigned((0x6996u >> (word & 15)) & 1u);
}

inline double PauliExpectation(const AmplitudeMap &input, double epsilon, const std::string &pauli)
{
    if (input.empty())
        return 0.;
    const size_t width = input.QubitCount();
    if (!std::isfinite(epsilon) || epsilon < 0. || pauli.size() > width)
        return PauliExpectationLegacy(input, epsilon, pauli);
    std::array<uint64_t, FastVectorBool::MaxWords> flip{}, phase{};
    unsigned yCount = 0;
    bool nonidentity = false, flips = false;
    for (size_t qubit = 0; qubit < pauli.size(); ++qubit)
    {
        const char symbol = pauli[qubit];
        const bool x = symbol == 'X' || symbol == 'x';
        const bool y = symbol == 'Y' || symbol == 'y';
        const bool z = symbol == 'Z' || symbol == 'z';
        const uint64_t bit = uint64_t{1} << (qubit % 64);
        if (x || y)
        {
            flip[qubit / 64] |= bit;
            flips = true;
        }
        if (z || y)
            phase[qubit / 64] |= bit;
        yCount += y;
        nonidentity = nonidentity || x || y || z;
    }
    if (nonidentity && !(1. > epsilon))
        return 0.;
    const std::array<std::complex<double>, 4> yPhase{1., {0., 1.}, -1., {0., -1.}};
    const auto factor = yPhase[yCount % 4];
    const size_t words = (width + 63) / 64;
    std::array<uint64_t, FastVectorBool::MaxWords> partnerWords;
    double result = 0.;
    for (const auto &entry : input)
    {
        if (!std::isfinite(entry.second.real()) || !std::isfinite(entry.second.imag()))
            return PauliExpectationLegacy(input, epsilon, pauli);
        const double mass = std::norm(entry.second);
        if (nonidentity && mass < epsilon)
            continue;
        uint64_t parityBits = 0;
        const auto *key = entry.first.getWords();
        for (size_t word = 0; word < words; ++word)
            parityBits ^= key[word] & phase[word];
        const double sign = PauliParity(parityBits) ? -1. : 1.;
        if (!flips)
        {
            result += sign * mass;
            continue;
        }
        for (size_t word = 0; word < words; ++word)
            partnerWords[word] = key[word] ^ flip[word];
        const auto partner = input.find(StateView(partnerWords.data(), width));
        if (partner != input.end())
            result += sign * std::real(std::conj(partner->second) * (factor * entry.second));
    }
    return result;
}
} // namespace Detail
} // namespace PathIntegral
} // namespace QC
