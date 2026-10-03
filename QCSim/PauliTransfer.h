#pragma once

#include "PauliPropTypes.h"
#include "LocalPauli.h"
#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

namespace QC { namespace PauliDetail {

// Sparse columns of P -> U^dagger P U. Local labels interleave x,z bits;
// q0 is the target for controlled gates. Tables are immutable once recorded.
// Bit p of unchanged is set when the gate maps local Pauli p to itself; the
// kernel skips those terms, which are most terms of sparse strings.
struct LocalTransfer {
    struct Entry { double coefficient; unsigned char pauli; };
    std::array<unsigned short, 65> offsets{};
    std::vector<Entry> entries;
    uint64_t unchanged = 0;
    unsigned char qubits = 1, maxOutputs = 1;
    bool clifford = true;

    explicit LocalTransfer(int width) : qubits(static_cast<unsigned char>(width)) {
        entries.reserve(width == 3 ? 232 : width == 2 ? 104 : 10);
    }
    template<size_t N> void Column(unsigned input, const std::array<double, N>& values) {
        offsets[input] = static_cast<unsigned short>(entries.size());
        for (unsigned output = 0; output < (1U << (2 * qubits)); ++output)
            if (values[output] != 0.) entries.push_back({values[output], static_cast<unsigned char>(output)});
        offsets[input + 1] = static_cast<unsigned short>(entries.size());
        const auto count = offsets[input + 1] - offsets[input];
        maxOutputs = std::max(maxOutputs, static_cast<unsigned char>(count));
        clifford = clifford && count == 1 && std::abs(entries[offsets[input]].coefficient) == 1.;
        if (count == 1 && entries[offsets[input]].pauli == input && entries[offsets[input]].coefficient == 1.)
            unchanged |= uint64_t(1) << input;
    }
};

struct Trig {
    double c, s;
    bool quarter;
    explicit Trig(double angle) {
        if (!std::isfinite(angle)) throw std::invalid_argument("Gate angle must be finite");
        quarter = detail::QuarterTurnCosSin(angle, c, s);
    }
};

using Rotation3 = std::array<std::array<double, 3>, 3>;
inline Rotation3 Multiply(const Rotation3& a, const Rotation3& b) {
    Rotation3 result{};
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 3; ++k) result[i][j] += a[i][k] * b[k][j];
    return result;
}
inline Rotation3 AxisRotation(int axis, const Trig& t) {
    Rotation3 r{};
    r[axis][axis] = 1.;
    const int a = (axis + 1) % 3, b = (axis + 2) % 3;
    r[a][a] = r[b][b] = t.c;
    r[a][b] = t.s; r[b][a] = -t.s;
    return r;
}
inline Rotation3 URotation(double theta, double phi, double lambda) {
    return Multiply(Multiply(AxisRotation(2, Trig(lambda)), AxisRotation(1, Trig(theta))), AxisRotation(2, Trig(phi)));
}
inline unsigned AxisCode(size_t axis) { return axis == 0 ? 1U : axis == 1 ? 3U : 2U; }
inline size_t CodeAxis(unsigned code) { return code == 1 ? 0 : code == 3 ? 1 : 2; }

inline detail::LocalPauli DecodeLocal(unsigned code, unsigned qubits) {
    detail::LocalPauli p;
    for (unsigned q = 0; q < qubits; ++q) {
        p.x |= ((code >> (2*q)) & 1U) << q;
        p.z |= ((code >> (2*q+1)) & 1U) << q;
    }
    return p;
}
inline unsigned EncodeLocal(detail::LocalPauli p, unsigned qubits) {
    unsigned code = 0;
    for (unsigned q = 0; q < qubits; ++q)
        code |= ((p.x >> q) & 1U) << (2*q) | ((p.z >> q) & 1U) << (2*q+1);
    return code;
}

inline std::shared_ptr<const LocalTransfer> CompileU(double theta, double phi, double lambda) {
    auto result = std::make_shared<LocalTransfer>(1);
    const auto rotation = URotation(theta, phi, lambda);
    for (unsigned p = 0; p < 4; ++p) {
        std::array<double, 4> values{};
        if (!p) values[0] = 1.;
        else for (size_t q = 0; q < 3; ++q) values[AxisCode(q)] = rotation[q][CodeAxis(p)];
        result->Column(p, values);
    }
    return result;
}

// V's Pauli coefficients are ordered I,X,Y,Z. The rotation describes V^dagger P V.
// diagonalDifference[j] is 1-R[j][j], evaluated without subtracting near-equal
// numbers. This preserves O(angle^2) branches even when cos(angle) rounds to 1.
inline std::shared_ptr<const LocalTransfer> CompileControlled(
    const std::array<std::complex<double>, 4>& v, const Rotation3& rotation,
    const std::array<double, 3>& diagonalDifference) {
    auto result = std::make_shared<LocalTransfer>(2);
    for (unsigned input = 0; input < 16; ++input) {
        std::array<double, 16> values{};
        const unsigned control = input >> 2, target = input & 3;
        if (control == 0 || control == 2) {
            if (target == 0) values[input] = 1.;
            else {
                const size_t j = CodeAxis(target);
                for (size_t i = 0; i < 3; ++i) {
                    const unsigned p = AxisCode(i);
                    const double same = i == j ? 1. - 0.5 * diagonalDifference[j] : 0.5 * rotation[i][j];
                    const double other = i == j ? 0.5 * diagonalDifference[j] : -0.5 * rotation[i][j];
                    values[p | (control << 2)] = same;
                    values[p | ((control ^ 2) << 2)] = other;
                }
            }
        } else {
            const auto left = DecodeLocal(target, 1);
            for (size_t i = 0; i < 4; ++i) {
                const auto right = DecodeLocal(i == 0 ? 0 : AxisCode(i-1), 1);
                const unsigned p = EncodeLocal({uint8_t(left.x ^ right.x), uint8_t(left.z ^ right.z)}, 1);
                const auto a = v[i] * detail::LocalPauliSum::IPower(detail::LocalPauliSum::ProductPhase(left, right));
                values[p | (1U << 2)] = control == 1 ? a.real() : a.imag();
                values[p | (3U << 2)] = control == 1 ? -a.imag() : a.real();
            }
        }
        result->Column(input, values);
    }
    return result;
}

inline std::shared_ptr<const LocalTransfer> CompileCR(int axis, double angle) {
    const Trig full(angle), half(0.5 * angle);
    std::array<std::complex<double>, 4> v{};
    v[0] = half.c; v[axis + 1] = {0., -half.s};
    std::array<double, 3> difference{};
    for (int j = 0; j < 3; ++j)
        if (j != axis) difference[j] = full.quarter ? 1. - full.c : 2. * half.s * half.s;
    return CompileControlled(v, AxisRotation(axis, full), difference);
}

inline std::shared_ptr<const LocalTransfer> CompileCP(double angle) {
    const Trig full(angle), half(0.5 * angle);
    const double d = full.quarter ? 1. - full.c : 2. * half.s * half.s;
    return CompileControlled({std::complex<double>(1. - 0.5*d, 0.5*full.s), 0., 0.,
        std::complex<double>(0.5*d, -0.5*full.s)}, AxisRotation(2, full), {d, d, 0.});
}

inline std::shared_ptr<const LocalTransfer> CompileCU(double theta, double phi, double lambda, double gamma) {
    const Trig t(theta), p(phi), l(lambda), g(gamma), half(0.5*theta);
    const Trig plus(0.5*phi + 0.5*lambda), minus(0.5*phi - 0.5*lambda);
    const std::array<double, 3> b{{-half.s * minus.s, half.s * minus.c, half.c * plus.s}};
    const auto phase = std::complex<double>(plus.c, plus.s) * std::complex<double>(g.c, g.s);
    std::array<std::complex<double>, 4> v{{phase * (half.c * plus.c),
        phase * std::complex<double>(0., -b[0]), phase * std::complex<double>(0., -b[1]),
        phase * std::complex<double>(0., -b[2])}};
    const auto r = URotation(theta, phi, lambda);
    std::array<double, 3> difference{};
    for (size_t j = 0; j < 3; ++j)
        difference[j] = t.quarter && p.quarter && l.quarter ? 1. - r[j][j]
            : 2. * (b[(j+1)%3]*b[(j+1)%3] + b[(j+2)%3]*b[(j+2)%3]);
    return CompileControlled(v, r, difference);
}

// Compile exact dyadic three-qubit unitaries with the local Pauli algebra
// shared with the extended stabilizer. All pair products are combined here,
// never in the hot loop.
inline std::shared_ptr<const LocalTransfer> CompileExact(const detail::LocalPauliSum& u, int qubits) {
    auto result = std::make_shared<LocalTransfer>(qubits);
    for (unsigned input = 0; input < (1U << (2*qubits)); ++input) {
        std::array<double, 64> values{};
        const auto p = DecodeLocal(input, qubits);
        u.ForEachTerm([&](detail::LocalPauli a, std::complex<double> ca) {
            const detail::LocalPauli ap{uint8_t(a.x ^ p.x), uint8_t(a.z ^ p.z)};
            u.ForEachTerm([&](detail::LocalPauli b, std::complex<double> cb) {
                const unsigned phase = detail::LocalPauliSum::ProductPhase(a, p) + detail::LocalPauliSum::ProductPhase(ap, b);
                const auto coefficient = std::conj(ca) * cb * detail::LocalPauliSum::IPower(phase);
                values[EncodeLocal({uint8_t(ap.x ^ b.x), uint8_t(ap.z ^ b.z)}, qubits)] += coefficient.real();
            });
        });
        result->Column(input, values);
    }
    return result;
}

inline std::shared_ptr<const LocalTransfer> FixedTransfer(OperationType type) {
    constexpr double halfPi = 1.57079632679489661923132169163975144;
    using Sum = detail::LocalPauliSum;
    switch (type) {
    case OperationType::CS: { static const auto map = CompileCP(halfPi); return map; }
    case OperationType::CSDAG: { static const auto map = CompileCP(-halfPi); return map; }
    case OperationType::CSX: case OperationType::CSXDAG: {
        const auto build = [halfPi](bool inverse) {
            const double sign = inverse ? -1. : 1.;
            return CompileControlled({std::complex<double>(.5, .5*sign),
                std::complex<double>(.5, -.5*sign), 0., 0.}, AxisRotation(0, Trig(sign*halfPi)), {0., 1., 1.});
        };
        static const auto forward = build(false), inverse = build(true);
        return type == OperationType::CSX ? forward : inverse;
    }
    case OperationType::CH: {
        static const auto map = CompileControlled({0., std::sqrt(.5), 0., std::sqrt(.5)},
            Rotation3{{{{0.,0.,1.}},{{0.,-1.,0.}},{{1.,0.,0.}}}}, {1.,2.,1.});
        return map;
    }
    case OperationType::CCX: {
        static const auto map = CompileExact(Sum::Identity() + Sum::Projector(1, true) * Sum::Projector(2, true)
            * (Sum::Pauli(0, 'X') + Sum::Identity(-1.)), 3);
        return map;
    }
    case OperationType::CSWAP: {
        static const auto map = [] {
            auto swap = Sum::Identity(.5);
            for (char axis : {'X','Y','Z'}) swap = swap + Sum::Pauli(0, axis, .5) * Sum::Pauli(1, axis);
            return CompileExact(Sum::Projector(2, false) + Sum::Projector(2, true) * swap, 3);
        }();
        return map;
    }
    default: throw std::invalid_argument("Not a fixed local Pauli gate");
    }
}

// The table of a parameterized gate; angles not used by the gate must be zero.
// Circuits often repeat angles (QAOA layers, Trotter steps), so each thread
// keeps the last tables it compiled in a small direct-mapped cache, keyed by
// the exact angle bits. A hit shares the immutable table.
inline std::shared_ptr<const LocalTransfer> ParameterizedTransfer(OperationType type,
    double a, double b = 0., double c = 0., double d = 0.) {
    struct Slot {
        OperationType type = OperationType::X;
        std::array<uint64_t, 4> bits{};
        std::shared_ptr<const LocalTransfer> table;
    };
    thread_local std::array<Slot, 16> cache;
    const double angles[] = {a, b, c, d};
    std::array<uint64_t, 4> bits;
    std::memcpy(bits.data(), angles, sizeof(angles));
    uint64_t hash = static_cast<uint64_t>(type);
    for (const auto word : bits) hash = (hash ^ word) * UINT64_C(0x9E3779B97F4A7C15);
    auto& slot = cache[hash >> 60];
    if (slot.table && slot.type == type && slot.bits == bits) return slot.table;
    std::shared_ptr<const LocalTransfer> table;
    switch (type) {
    case OperationType::U: table = CompileU(a, b, c); break;
    case OperationType::CU: table = CompileCU(a, b, c, d); break;
    case OperationType::CRX: table = CompileCR(0, a); break;
    case OperationType::CRY: table = CompileCR(1, a); break;
    case OperationType::CRZ: table = CompileCR(2, a); break;
    case OperationType::CP: table = CompileCP(a); break;
    default: throw std::invalid_argument("Not a parameterized local Pauli gate");
    }
    slot.type = type;
    slot.bits = bits;
    slot.table = table;
    return table;
}

}}
