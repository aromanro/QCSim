#pragma once

#include "LocalPauli.h"
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace QC
{
namespace detail
{

// The generic bridges cache immutable compilations by matrix value bits,
// not by a gate object's address or dynamic type. Mutating a gate therefore
// cannot reuse an obsolete transfer. At most eight rows are supported.
struct LocalMatrixKey
{
    unsigned dimension = 0;
    // Only the dimension's active prefix is initialized, compared or copied.
    // A one-qubit lookup must not clear/copy the unused 8x8 capacity.
    std::array<uint64_t, 128> bits;
    uint64_t hash = UINT64_C(1469598103934665603);

    LocalMatrixKey() = default;

    LocalMatrixKey(const LocalMatrixKey &other) : dimension(other.dimension), hash(other.hash)
    {
        std::copy_n(other.bits.begin(), 2 * dimension * dimension, bits.begin());
    }

    LocalMatrixKey &operator=(const LocalMatrixKey &other)
    {
        if (this == &other)
            return *this;
        dimension = other.dimension;
        hash = other.hash;
        std::copy_n(other.bits.begin(), 2 * dimension * dimension, bits.begin());
        return *this;
    }

    bool operator==(const LocalMatrixKey &other) const
    {
        return dimension == other.dimension &&
               std::equal(bits.begin(), bits.begin() + 2 * dimension * dimension, other.bits.begin());
    }
};

template <class Matrix> LocalMatrixKey MatrixKey(const Matrix &matrix)
{
    const auto dimension = matrix.rows();
    if ((dimension != 2 && dimension != 4 && dimension != 8) || matrix.cols() != dimension || !matrix.allFinite())
        throw std::invalid_argument("Gate requires a finite square matrix on one to three qubits");
    LocalMatrixKey key;
    key.dimension = unsigned(dimension);
    size_t index = 0;
    for (unsigned row = 0; row < key.dimension; ++row)
        for (unsigned col = 0; col < key.dimension; ++col)
        {
            const std::complex<double> value = matrix(row, col);
            const double parts[] = {value.real(), value.imag()};
            for (double part : parts)
            {
                // Conjugation changes the sign of zero without changing the
                // operator. Use one encoding in both the key and its hash.
                if (part == 0.)
                    part = 0.;
                std::memcpy(&key.bits[index], &part, sizeof(part));
                key.hash = (key.hash ^ key.bits[index++]) * UINT64_C(1099511628211);
            }
        }
    key.hash = (key.hash ^ key.dimension) * UINT64_C(1099511628211);
    return key;
}

// The key has already been validated by MatrixKey. Transpose its logical
// row-major bit pairs and conjugate nonzero imaginary parts; canonical zeros
// must keep their positive-zero encoding. This hashes no unused capacity.
inline uint64_t AdjointMatrixKeyHash(const LocalMatrixKey &key) noexcept
{
    static_assert(sizeof(double) == sizeof(uint64_t) && std::numeric_limits<double>::is_iec559,
                  "Matrix-key conjugation requires IEEE 754 binary64");
    uint64_t hash = UINT64_C(1469598103934665603);
    for (unsigned row = 0; row < key.dimension; ++row)
        for (unsigned col = 0; col < key.dimension; ++col)
        {
            const unsigned index = 2 * (col * key.dimension + row);
            const uint64_t imaginary = key.bits[index + 1];
            hash = (hash ^ key.bits[index]) * UINT64_C(1099511628211);
            hash = (hash ^ (imaginary == 0 ? 0 : imaginary ^ UINT64_C(0x8000000000000000))) * UINT64_C(1099511628211);
        }
    return (hash ^ key.dimension) * UINT64_C(1099511628211);
}

// A second hash of a key, or of its adjoint, that rotates between words. FNV
// over whole 64-bit words carries a flipped sign bit only into the top bit of
// the result, so the plain hashes of a matrix and of its adjoint coincide
// whenever conjugation flips an even number of signs, as for every diagonal
// rotation. Compilations store this one to reject inverse candidates; it is
// computed once per compilation, not per gate.
inline uint64_t MatrixKeyFingerprint(const LocalMatrixKey &key, bool adjoint) noexcept
{
    uint64_t hash = UINT64_C(0x9E3779B97F4A7C15) ^ key.dimension;
    const auto mix = [&hash](uint64_t word) {
        hash = (hash ^ word) * UINT64_C(0xBF58476D1CE4E5B9);
        hash = (hash << 27) | (hash >> 37);
    };
    for (unsigned row = 0; row < key.dimension; ++row)
        for (unsigned col = 0; col < key.dimension; ++col)
        {
            const unsigned index = 2 * (adjoint ? col * key.dimension + row : row * key.dimension + col);
            const uint64_t imaginary = key.bits[index + 1];
            mix(key.bits[index]);
            mix(adjoint && imaginary != 0 ? imaginary ^ UINT64_C(0x8000000000000000) : imaginary);
        }
    return hash;
}

// A small set-associative cache of immutable compilations keyed by matrix
// contents. A hit compares the complete key; the tags (the key's hash, never
// zero) only reject the other entries of a set quickly. A set fills in order
// and is then replaced round-robin, so two matrices that share a set no
// longer evict each other.
template <class Value> struct LocalMatrixCache
{
    static constexpr unsigned Sets = 16, Ways = 4;

    struct Slot
    {
        LocalMatrixKey key;
        Value value;
    };

    // Caches of different compilations take different multipliers, so a pair
    // of matrices that shares a set in one cache rarely shares one in another.
    static unsigned Set(uint64_t hash, uint64_t multiplier) noexcept
    {
        return unsigned(((hash ^ (hash >> 29)) * multiplier) >> 60);
    }

    template <class Accept> Slot *Find(const LocalMatrixKey &key, unsigned set, Accept accept)
    {
        const uint64_t tag = key.hash | 1U;
        for (unsigned index = set * Ways; index < (set + 1) * Ways; ++index)
            if (tags[index] == tag && slots[index].key == key && accept(slots[index].value))
                return &slots[index];
        return nullptr;
    }

    // The slot that receives a new compilation of this key: a free one of
    // the set, or the next one in turn. The caller assigns the value.
    Slot &Insert(const LocalMatrixKey &key, unsigned set)
    {
        unsigned index = set * Ways;
        while (index < (set + 1) * Ways && tags[index] != 0)
            ++index;
        if (index == (set + 1) * Ways)
        {
            index = set * Ways + next[set];
            next[set] = uint8_t((next[set] + 1) % Ways);
        }
        slots[index].key = key;
        tags[index] = key.hash | 1U;
        return slots[index];
    }

    std::array<uint64_t, Sets * Ways> tags{};
    std::array<uint8_t, Sets> next{};
    std::array<Slot, Sets * Ways> slots;
};

inline LocalPauli PairPauli(char first, char second)
{
    const auto a = LocalPauli::Single(0, first), b = LocalPauli::Single(1, second);
    return {uint8_t(a.x | b.x), uint8_t(a.z | b.z)};
}

// Exact finite Pauli expansion, with no tolerance pruning of small coefficients.
// The generic bridge is deliberately limited to the local algebra's three bits.
template <class Matrix> LocalPauliSum MatrixPauliSum(const Matrix &matrix, unsigned qubits)
{
    if (!qubits || qubits > 3 || matrix.rows() != (1U << qubits) || matrix.cols() != matrix.rows() ||
        !matrix.allFinite())
        throw std::invalid_argument("Gate requires a finite square matrix on one to three qubits");
    const unsigned dimension = 1U << qubits;
    LocalPauliSum sum;
    for (unsigned x = 0; x < dimension; ++x)
        for (unsigned z = 0; z < dimension; ++z)
        {
            std::complex<double> coefficient = 0.;
            const auto phase = LocalPauliSum::IPower(Clifford::detail::Popcount(x & z));
            for (unsigned col = 0; col < dimension; ++col)
                coefficient += ((Clifford::detail::Popcount(col & z) & 1U) ? -phase : phase) * matrix(col, col ^ x);
            sum.Add({uint8_t(x), uint8_t(z)}, coefficient / double(dimension));
        }
    return sum;
}

template <class Gate> void AppendPairRotation(Gate &gate, char first, char second, double angle)
{
    const char axes[] = {first, second};
    for (unsigned q = 0; q < 2; ++q)
    {
        if (axes[q] == 'Y')
            gate.Sdg(q);
        if (axes[q] != 'Z')
            gate.H(q);
    }
    gate.CX(1, 0).Rotate(1, 'Z', angle).CX(1, 0);
    for (unsigned q = 0; q < 2; ++q)
    {
        if (axes[q] != 'Z')
            gate.H(q);
        if (axes[q] == 'Y')
            gate.S(q);
    }
}

inline LocalPauliSum XXYYPauliSum(bool plus, double theta, double beta)
{
    if (!std::isfinite(theta) || !std::isfinite(beta))
        throw std::invalid_argument("Gate angles must be finite");
    double c, s, cb, sb;
    QuarterTurnCosSin(theta * .5, c, s);
    QuarterTurnCosSin(beta, cb, sb);
    // cos(theta/2)-1 is evaluated through sin(theta/4), preserving tiny angles.
    double delta = -2. * std::sin(theta * .25) * std::sin(theta * .25);
    long long turns;
    if (TryGetQuarterTurns(theta * .5, turns))
        delta = c - 1.;
    const double sign = plus ? 1. : -1.;
    auto sum = LocalPauliSum::Identity(1. + .5 * delta);
    sum.Add(PairPauli('Z', 'Z'), -.5 * sign * delta);
    sum.Add(PairPauli('X', 'X'), {0., -.5 * s * cb});
    sum.Add(PairPauli('Y', 'Y'), {0., -.5 * sign * s * cb});
    sum.Add(PairPauli('X', 'Y'), {0., -.5 * s * sb});
    sum.Add(PairPauli('Y', 'X'), {0., .5 * sign * s * sb});
    return sum;
}

} // namespace detail
} // namespace QC
