#pragma once

#include <cmath>
#include <limits>

namespace QC
{
namespace detail
{
// Recognize only the exact floating-point encoding of k*pi/2 for |k|<=2^20.
// At huge magnitudes floating-point spacing can make arbitrary angles match
// a rounded product; those inputs retain ordinary trigonometric values.
// A neighboring representable angle is distinct and must never be snapped.
inline bool TryGetQuarterTurns(double angle, long long &quarterTurns)
{
    constexpr double quarterTurnAngle = 1.57079632679489661923132169163975144;
    // The only encoded quarter turn in this open interval is zero. Most
    // variational angles take this exact rejection path, avoiding division
    // and integer rounding without changing any near-quarter decision.
    if (std::abs(angle) < quarterTurnAngle)
    {
        quarterTurns = 0;
        return angle == 0.;
    }
    // The candidate count needs no division and no library rounding call: the product with
    // 2/pi is within a few units in the last place of the quotient, far from a half for every
    // accepted count, and the exact comparison below is what decides.
    constexpr double quarterTurnsPerRadian = 0.636619772367581343075535053490057448;
    const double turns = std::rint(angle * quarterTurnsPerRadian);
    if (!(std::abs(turns) <= 1048576.))
        return false;
    if (angle != turns * quarterTurnAngle)
        return false;
    quarterTurns = static_cast<long long>(turns);
    return true;
}

inline bool QuarterTurnCosSin(double angle, double &cosine, double &sine)
{
    if (angle == 0.)
    {
        cosine = 1.;
        sine = 0.;
        return true;
    }
    long long turns;
    if (TryGetQuarterTurns(angle, turns))
    {
        static constexpr double values[4][2] = {{1., 0.}, {0., 1.}, {-1., 0.}, {0., -1.}};
        const auto &value = values[((turns % 4) + 4) % 4];
        cosine = value[0];
        sine = value[1];
        return true;
    }
    cosine = std::cos(angle);
    sine = std::sin(angle);
    return false;
}

// Half-angle coefficients of an encoded quarter-turn rotation. Odd turns
// use the same rounded sqrt(1/2) for both magnitudes, preserving Clifford
// sparsity without classifying almost-Clifford rotations as exact ones.
inline bool RotationHalfAngleCosSin(double angle, double &cosine, double &sine)
{
    if (angle == 0.)
    {
        cosine = 1.;
        sine = 0.;
        return true;
    }
    long long turns;
    if (TryGetQuarterTurns(angle, turns))
    {
        const double h = std::sqrt(.5);
        const double values[8][2] = {{1., 0.}, {h, h}, {0., 1.}, {-h, h}, {-1., 0.}, {-h, -h}, {0., -1.}, {h, -h}};
        const auto &value = values[((turns % 8) + 8) % 8];
        cosine = value[0];
        sine = value[1];
        return true;
    }
    cosine = std::cos(.5 * angle);
    sine = std::sin(.5 * angle);
    return false;
}
} // namespace detail
} // namespace QC
