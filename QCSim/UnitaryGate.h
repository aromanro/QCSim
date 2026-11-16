#pragma once

#include "GateAngles.h"
#include <Eigen/Eigen>
#include <array>
#include <complex>
#include <limits>

namespace QC
{
namespace Gates
{
// A candidate records a constructor's analytic meaning, not a unitary promise.
// A bridge must verify the current matrix before using it. Ordinary dense
// execution does not inspect this metadata or allocate a matrix snapshot.
struct UnitaryGateDescriptor
{
    enum class Kind : unsigned char
    {
        None,
        Identity,
        H,
        Hy,
        S,
        Sdg,
        T,
        Tdg,
        X,
        Y,
        Z,
        Sx,
        Sxdg,
        Splitter,
        Phase,
        Rx,
        Ry,
        Rz,
        U,
        PairRotation,
        XXPlusYY,
        XXMinusYY,
        Fredkin
    };
    Kind kind = Kind::None;
    std::array<double, 4> angles{};
    unsigned controls = 0;
    char first = 0, second = 0;
    bool adjoint = false;
    // The parameter-dependent coefficients exactly as the constructor stored
    // them in the matrix (see UnitaryGateMatrix for the layout of each kind).
    // Rebuilding the matrix copies these instead of evaluating trigonometry or
    // complex products again: the same expression can round differently at
    // another call site (fused multiply-add contraction, compile-time folding),
    // which must not make a constructor's own matrix fail its verification.
    std::array<std::complex<double>, 4> built{};
    bool hasBuilt = false;

    explicit operator bool() const
    {
        return kind != Kind::None;
    }

    bool operator==(const UnitaryGateDescriptor &other) const
    {
        return kind == other.kind && angles == other.angles && controls == other.controls && first == other.first &&
               second == other.second && adjoint == other.adjoint && hasBuilt == other.hasBuilt &&
               (!hasBuilt || built == other.built);
    }

    // A candidate together with the coefficients its constructor stored.
    static UnitaryGateDescriptor Built(Kind kind, const std::array<double, 4> &angles, unsigned controls,
                                       const std::array<std::complex<double>, 4> &built, char first = 0,
                                       char second = 0)
    {
        UnitaryGateDescriptor descriptor;
        descriptor.kind = kind;
        descriptor.angles = angles;
        descriptor.controls = controls;
        descriptor.first = first;
        descriptor.second = second;
        descriptor.built = built;
        descriptor.hasBuilt = true;
        return descriptor;
    }
};

// Compare candidates only after their current matrices have been verified.
// This is an exact semantic relation, not a tolerance-based matrix inverse
// test: tiny parameter changes and nearby angles must remain observable.
// A template only so that the definition needs no `inline`, which GCC reports
// next to the attribute that keeps this cold path out of its callers.
template <class Descriptor>
EIGEN_DONT_INLINE bool AreUnitaryGateDescriptorAliasesInverse(const Descriptor &a, const Descriptor &b);

inline bool AreUnitaryGateDescriptorsInverse(const UnitaryGateDescriptor &a, const UnitaryGateDescriptor &b)
{
    using Kind = UnitaryGateDescriptor::Kind;
    if (a.controls != b.controls || a.controls > 2)
        return false;
    const auto supported = [](const UnitaryGateDescriptor &d) {
        switch (d.kind)
        {
        case Kind::None:
            return false;
        case Kind::PairRotation:
            return d.controls == 0 && (d.first == 'X' || d.first == 'Y' || d.first == 'Z') &&
                   (d.second == 'X' || d.second == 'Y' || d.second == 'Z') && std::isfinite(d.angles[0]);
        case Kind::XXPlusYY:
        case Kind::XXMinusYY:
            return d.controls == 0 && std::isfinite(d.angles[0]) && std::isfinite(d.angles[1]);
        case Kind::Fredkin:
            return d.controls == 0;
        case Kind::U:
            return std::isfinite(d.angles[0]) && std::isfinite(d.angles[1]) && std::isfinite(d.angles[2]) &&
                   std::isfinite(d.angles[3]);
        case Kind::Phase:
        case Kind::Rx:
        case Kind::Ry:
        case Kind::Rz:
            return std::isfinite(d.angles[0]);
        case Kind::Identity:
        case Kind::H:
        case Kind::Hy:
        case Kind::S:
        case Kind::Sdg:
        case Kind::T:
        case Kind::Tdg:
        case Kind::X:
        case Kind::Y:
        case Kind::Z:
        case Kind::Sx:
        case Kind::Sxdg:
        case Kind::Splitter:
            return true;
        default:
            return false;
        }
    };
    const auto matches = [&]() {
        if (a.adjoint != b.adjoint)
        {
            if (a.kind != b.kind)
                return false;
            switch (a.kind)
            {
            case Kind::U:
                return a.angles == b.angles;
            case Kind::PairRotation:
                return a.first == b.first && a.second == b.second && a.angles[0] == b.angles[0];
            case Kind::XXPlusYY:
            case Kind::XXMinusYY:
                return a.angles[0] == b.angles[0] && a.angles[1] == b.angles[1];
            case Kind::Phase:
            case Kind::Rx:
            case Kind::Ry:
            case Kind::Rz:
                return a.angles[0] == b.angles[0];
            default:
                return true;
            }
        }
        switch (a.kind)
        {
        case Kind::S:
            return b.kind == Kind::Sdg;
        case Kind::Sdg:
            return b.kind == Kind::S;
        case Kind::T:
            return b.kind == Kind::Tdg;
        case Kind::Tdg:
            return b.kind == Kind::T;
        case Kind::Sx:
            return b.kind == Kind::Sxdg;
        case Kind::Sxdg:
            return b.kind == Kind::Sx;
        case Kind::Identity:
        case Kind::H:
        case Kind::Hy:
        case Kind::X:
        case Kind::Y:
        case Kind::Z:
        case Kind::Fredkin:
            return a.kind == b.kind;
        case Kind::Phase:
        case Kind::Rx:
        case Kind::Ry:
        case Kind::Rz:
            return a.kind == b.kind && a.angles[0] == -b.angles[0];
        case Kind::PairRotation:
            return a.kind == b.kind && a.first == b.first && a.second == b.second && a.angles[0] == -b.angles[0];
        case Kind::XXPlusYY:
        case Kind::XXMinusYY:
            return a.kind == b.kind && a.angles[0] == -b.angles[0] && a.angles[1] == b.angles[1];
        case Kind::U:
            return b.kind == Kind::U && a.angles[0] == -b.angles[0] && a.angles[1] == -b.angles[2] &&
                   a.angles[2] == -b.angles[1] && a.angles[3] == -b.angles[3];
        default:
            return false;
        }
    };
    // Most forward gates fail the first exact angle comparison. Validate
    // parameter domains only for candidates that actually match an inverse.
    if (matches())
        return supported(a) && supported(b);
    // Keep ordinary forward gates on the same-kind rejection path. Only
    // cross-family candidates need the cold exact-alias normalization.
    if (a.kind == b.kind)
        return false;
    const auto fixedPhase = [](Kind kind) {
        return kind == Kind::Identity || kind == Kind::Z || (kind >= Kind::S && kind <= Kind::Tdg);
    };
    const auto reducibleU = [](const UnitaryGateDescriptor &d) {
        return d.angles[3] == 0. && ((d.angles[1] == 0. && d.angles[2] == 0.) ||
                                     (d.angles[0] == 0. && (d.angles[1] == 0. || d.angles[2] == 0.)));
    };
    const auto couldAlias = [&]() {
        if (a.kind == Kind::U)
            return reducibleU(a) && (b.kind == Kind::Ry || b.kind == Kind::Phase || fixedPhase(b.kind));
        if (b.kind == Kind::U)
            return reducibleU(b) && (a.kind == Kind::Ry || a.kind == Kind::Phase || fixedPhase(a.kind));
        if (a.kind == Kind::Phase)
            return fixedPhase(b.kind) || (b.kind == Kind::Ry && b.angles[0] == 0.);
        if (b.kind == Kind::Phase)
            return fixedPhase(a.kind) || (a.kind == Kind::Ry && a.angles[0] == 0.);
        return (a.kind == Kind::Identity && b.kind == Kind::Ry && b.angles[0] == 0.) ||
               (b.kind == Kind::Identity && a.kind == Kind::Ry && a.angles[0] == 0.);
    };
    return couldAlias() && AreUnitaryGateDescriptorAliasesInverse(a, b) && supported(a) && supported(b);
}

template <class Descriptor>
EIGEN_DONT_INLINE bool AreUnitaryGateDescriptorAliasesInverse(const Descriptor &a, const Descriptor &b)
{
    using Kind = UnitaryGateDescriptor::Kind;
    const auto normalize = [](UnitaryGateDescriptor d) {
        // The normalized candidate is compared by kind and angles only; the
        // stored coefficients of the original kind no longer describe it.
        d.hasBuilt = false;
        d.built = {};
        if (d.kind == Kind::U && d.angles[3] == 0.)
        {
            double angle;
            if (d.angles[1] == 0. && d.angles[2] == 0.)
            {
                angle = d.angles[0];
                d.kind = Kind::Ry;
            }
            else if (d.angles[0] == 0. && (d.angles[1] == 0. || d.angles[2] == 0.))
            {
                angle = d.angles[d.angles[1] == 0. ? 2 : 1];
                d.kind = Kind::Phase;
            }
            else
                return d;
            d.angles = {angle};
        }
        if (d.kind == Kind::Ry && d.angles[0] == 0.)
        {
            d.kind = Kind::Identity;
            d.angles = {};
        }
        if (d.kind == Kind::Phase)
        {
            constexpr double quarter = 1.57079632679489661923132169163975144;
            const double angle = d.angles[0];
            if (angle == 0.)
                d.kind = Kind::Identity;
            else if (angle == .5 * quarter)
                d.kind = Kind::T;
            else if (angle == -.5 * quarter)
                d.kind = Kind::Tdg;
            else if (angle == quarter)
                d.kind = Kind::S;
            else if (angle == -quarter)
                d.kind = Kind::Sdg;
            else if (angle == 2. * quarter || angle == -2. * quarter)
                d.kind = Kind::Z;
            if (d.kind != Kind::Phase)
                d.angles = {};
        }
        return d;
    };
    const auto first = normalize(a), second = normalize(b);
    if (first.kind == a.kind && second.kind == b.kind)
        return false;
    // Normalization only removes exact zero parameters or substitutes exact
    // fixed phases. The caller still validates the original descriptors.
    return AreUnitaryGateDescriptorsInverse(first, second);
}

EIGEN_STRONG_INLINE std::complex<double> UnitaryPhase(double angle)
{
    // Keep ordinary constructor phases as direct polar calls so compilers
    // can fold constant Euler parameters. No nonzero quarter turn lies
    // inside this open interval; zero is already exact in std::polar.
    if (std::abs(angle) < 1.57079632679489661923132169163975144)
        return std::polar(1., angle);
    double c, s;
    detail::QuarterTurnCosSin(angle, c, s);
    return {c, s};
}

using UnitaryLocalMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 8, 8>;

// The product of two unit phases. A phase times its exact conjugate must come out real:
// the two cross terms cancel only if they round alike, which a compiler that fuses one of
// them into a multiply-add breaks. That case is therefore decided on the operands, not left
// to the arithmetic; everywhere else the products are ordinary.
inline std::complex<double> MultiplyUnitaryPhases(const std::complex<double> &a, const std::complex<double> &b)
{
    const double real = a.real() * b.real() - a.imag() * b.imag();
    if (a.real() == b.real() && a.imag() == -b.imag())
        return {real, 0.};
    return {real, a.real() * b.imag() + a.imag() * b.real()};
}

EIGEN_STRONG_INLINE std::array<std::complex<double>, 4> UnitaryEulerCoefficients(double theta, double phi,
                                                                                 double lambda, double gamma)
{
    double c, s;
    detail::RotationHalfAngleCosSin(theta, c, s);
    const auto g = UnitaryPhase(gamma), p = UnitaryPhase(phi), l = UnitaryPhase(lambda);
    // Compose phases as unit complex factors, never by adding their angles.
    // Group p*l before g so cancelling large phases do not erase a tiny
    // common phase (see MultiplyUnitaryPhases for the exact cancellation).
    return {g * c, -MultiplyUnitaryPhases(g, l) * s, MultiplyUnitaryPhases(g, p) * s,
            MultiplyUnitaryPhases(g, MultiplyUnitaryPhases(p, l)) * c};
}

inline UnitaryLocalMatrix UnitaryGateMatrix(const UnitaryGateDescriptor &descriptor)
{
    using Kind = UnitaryGateDescriptor::Kind;
    const bool pair = descriptor.kind == Kind::PairRotation || descriptor.kind == Kind::XXPlusYY ||
                      descriptor.kind == Kind::XXMinusYY;
    if (!descriptor || descriptor.controls > 2 || ((pair || descriptor.kind == Kind::Fredkin) && descriptor.controls))
        return {};
    if (descriptor.kind == Kind::PairRotation)
    {
        const auto validAxis = [](char axis) { return axis == 'X' || axis == 'Y' || axis == 'Z'; };
        if (!validAxis(descriptor.first) || !validAxis(descriptor.second))
            return {};
    }
    const unsigned localDimension = descriptor.kind == Kind::Fredkin ? 8 : pair ? 4 : 2;
    const unsigned dimension = localDimension << descriptor.controls;
    for (double angle : descriptor.angles)
        if (!std::isfinite(angle))
            return {};
    // Stored constructor coefficients, when present, are copied below instead
    // of being evaluated again. Layout: T, Tdg and Phase keep the phase entry;
    // Rx, Ry, Rz and pair rotations keep (cos, sin) of the half angle; U keeps
    // its four entries in row order; XX+-YY keeps the cosine and the two
    // off-diagonal entries of its active sector.
    const bool stored = descriptor.hasBuilt;
    const auto &built = descriptor.built;
    if (stored)
        for (const auto &value : built)
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                return {};
    UnitaryLocalMatrix result = UnitaryLocalMatrix::Identity(dimension, dimension);
    const unsigned offset = dimension - localDimension;
    auto block = result.block(offset, offset, localDimension, localDimension);
    block.setZero();
    const double h = 1. / std::sqrt(2.);
    const auto theta = descriptor.angles[0];
    double c = 0., s = 0.;
    switch (descriptor.kind)
    {
    case Kind::Fredkin:
        block.setIdentity();
        block(5, 5) = block(6, 6) = 0.;
        block(5, 6) = block(6, 5) = 1.;
        break;
    case Kind::Identity:
        block.setIdentity();
        break;
    case Kind::H:
        block << h, h, h, -h;
        break;
    case Kind::Hy:
        block << h, std::complex<double>(0., -h), std::complex<double>(0., h), -h;
        break;
    case Kind::S:
        block(0, 0) = 1.;
        block(1, 1) = {0., 1.};
        break;
    case Kind::Sdg:
        block(0, 0) = 1.;
        block(1, 1) = {0., -1.};
        break;
    case Kind::T:
        block(0, 0) = 1.;
        block(1, 1) = stored ? built[0] : std::polar(1., std::acos(-1.) / 4.);
        break;
    case Kind::Tdg:
        block(0, 0) = 1.;
        block(1, 1) = stored ? built[0] : std::polar(1., -std::acos(-1.) / 4.);
        break;
    case Kind::X:
        block(0, 1) = block(1, 0) = 1.;
        break;
    case Kind::Y:
        block(0, 1) = {0., -1.};
        block(1, 0) = {0., 1.};
        break;
    case Kind::Z:
        block(0, 0) = 1.;
        block(1, 1) = -1.;
        break;
    case Kind::Sx:
    case Kind::Sxdg: {
        const double sign = descriptor.kind == Kind::Sx ? 1. : -1.;
        block(0, 0) = block(1, 1) = {.5, sign * .5};
        block(0, 1) = block(1, 0) = {.5, -sign * .5};
        break;
    }
    case Kind::Splitter:
        block(0, 0) = block(1, 1) = h;
        block(0, 1) = block(1, 0) = {0., h};
        break;
    case Kind::Phase:
        block(0, 0) = 1.;
        block(1, 1) = stored ? built[0] : UnitaryPhase(theta);
        break;
    case Kind::Rx:
    case Kind::Ry:
    case Kind::Rz:
        if (stored)
        {
            c = built[0].real();
            s = built[0].imag();
        }
        else
            detail::RotationHalfAngleCosSin(theta, c, s);
        if (descriptor.kind == Kind::Rz)
        {
            block(0, 0) = {c, -s};
            block(1, 1) = {c, s};
        }
        else
        {
            block(0, 0) = block(1, 1) = c;
            block(0, 1) = descriptor.kind == Kind::Rx ? std::complex<double>(0., -s) : std::complex<double>(-s);
            block(1, 0) = descriptor.kind == Kind::Rx ? std::complex<double>(0., -s) : std::complex<double>(s);
        }
        break;
    case Kind::U: {
        const auto values =
            stored ? built
                   : UnitaryEulerCoefficients(theta, descriptor.angles[1], descriptor.angles[2], descriptor.angles[3]);
        block(0, 0) = values[0];
        block(0, 1) = values[1];
        block(1, 0) = values[2];
        block(1, 1) = values[3];
        break;
    }
    case Kind::PairRotation: {
        if (stored)
        {
            c = built[0].real();
            s = built[0].imag();
        }
        else
            detail::RotationHalfAngleCosSin(theta, c, s);
        const char axes[] = {descriptor.first, descriptor.second};
        for (unsigned col = 0; col < 4; ++col)
        {
            unsigned row = col;
            std::complex<double> phase(0., -s);
            for (unsigned q = 0; q < 2; ++q)
            {
                const bool bit = (col >> q) & 1U;
                if (axes[q] != 'Z')
                    row ^= 1U << q;
                if (axes[q] == 'Y')
                    phase *= std::complex<double>(0., bit ? -1. : 1.);
                else if (axes[q] == 'Z' && bit)
                    phase = -phase;
            }
            block(col, col) = c;
            block(row, col) += phase;
        }
        break;
    }
    case Kind::XXPlusYY:
    case Kind::XXMinusYY: {
        block.setIdentity();
        const unsigned a = descriptor.kind == Kind::XXPlusYY ? 1 : 0;
        const unsigned b = descriptor.kind == Kind::XXPlusYY ? 2 : 3;
        if (stored)
        {
            block(a, a) = block(b, b) = built[0].real();
            block(a, b) = built[1];
            block(b, a) = built[2];
            break;
        }
        detail::RotationHalfAngleCosSin(theta, c, s);
        double cb, sb;
        detail::QuarterTurnCosSin(descriptor.angles[1], cb, sb);
        block(a, a) = block(b, b) = c;
        block(a, b) = std::complex<double>(0., -s) * std::complex<double>(cb, -sb);
        block(b, a) = std::complex<double>(0., -s) * std::complex<double>(cb, sb);
        break;
    }
    default:
        return {};
    }
    if (descriptor.adjoint)
        result.adjointInPlace();
    return result;
}

template <class Matrix> bool VerifyUnitaryGateDescriptor(const Matrix &matrix, const UnitaryGateDescriptor &descriptor)
{
    const auto expected = UnitaryGateMatrix(descriptor);
    if (!expected.size() || matrix.rows() != expected.rows() || matrix.cols() != expected.cols())
        return false;
    for (Eigen::Index col = 0; col < matrix.cols(); ++col)
        for (Eigen::Index row = 0; row < matrix.rows(); ++row)
            if (matrix(row, col) != expected(row, col))
                return false;
    return true;
}

template <class Gate> auto GetUnitaryGateDescriptorImpl(const Gate &gate, int) -> decltype(gate.getUnitaryDescriptor())
{
    return gate.getUnitaryDescriptor();
}

template <class Gate> UnitaryGateDescriptor GetUnitaryGateDescriptorImpl(const Gate &, long)
{
    return {};
}

// This cheap candidate query deliberately does not re-evaluate trigonometry.
// Cache by matrix contents AND this descriptor; verify on a compilation miss.
template <class Gate> UnitaryGateDescriptor GetUnitaryGateDescriptor(const Gate &gate)
{
    return GetUnitaryGateDescriptorImpl(gate, 0);
}
} // namespace Gates
} // namespace QC
