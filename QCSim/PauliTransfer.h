#pragma once

#include "LocalCliffordRecognition.h"
#include "LocalGateExtensions.h"
#include "LocalPauli.h"
#include "PauliPropTypes.h"
#include <algorithm>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace QC
{
namespace PauliDetail
{

// Sparse columns of P -> U^dagger P U. Local labels interleave x,z bits;
// q0 is the target for controlled gates. Tables are immutable once recorded.
// Bit p of unchanged is set when the gate maps local Pauli p to itself; the
// kernel skips those terms, which are most terms of sparse strings.
struct LocalTransfer
{
    struct Entry
    {
        double coefficient;
        unsigned char pauli;
    };

    std::array<unsigned short, 65> offsets{};
    std::vector<Entry> entries;
    uint64_t unchanged = 0;
    unsigned char qubits = 1, maxOutputs = 1;
    bool clifford = true;
    // The members below describe the table to the recorder; the kernels read
    // only the ones above, which keep their places. The two flags use what
    // was padding and the angle adds one double, so a shared table still
    // takes the allocation size it took without them.
    // When rotationAxis is 0, 1 or 2 this one-qubit table is the action of the
    // native rotation about X, Y or Z by rotationAngle, so a recorder may use
    // that kernel instead of the table (see RotationCosineSine).
    signed char rotationAxis = -1;
    // The table of a matrix M with M^dagger M away from the identity by more
    // than 1e-12: it does not preserve the trace, so outcome probabilities,
    // which take <I> = 1, are not defined after it. This only tells whether
    // probabilities may be asked for; it never changes a coefficient.
    bool nonunitary = false;
    double rotationAngle = 0.;

    explicit LocalTransfer(int width) : qubits(static_cast<unsigned char>(width))
    {
        entries.reserve(width == 3 ? 232 : width == 2 ? 104 : 10);
    }

    template <size_t N> void Column(unsigned input, const std::array<double, N> &values)
    {
        offsets[input] = static_cast<unsigned short>(entries.size());
        for (unsigned output = 0; output < (1U << (2 * qubits)); ++output)
            if (values[output] != 0.)
                entries.push_back({values[output], static_cast<unsigned char>(output)});
        Close(input);
    }

    // The same column from the outputs listed in a bit mask: a compiler that
    // writes a few of the 64 outputs need not scan the others. Bits are
    // visited in increasing order, so the entries are those of Column. The
    // values of the listed outputs are set back to zero for the next column.
    template <size_t N> void TakeColumn(unsigned input, std::array<double, N> &values, uint64_t outputs)
    {
        offsets[input] = static_cast<unsigned short>(entries.size());
        for (; outputs; outputs &= outputs - 1)
        {
            const unsigned output = static_cast<unsigned>(Clifford::detail::TrailingZero(outputs));
            if (values[output] != 0.)
            {
                entries.push_back({values[output], static_cast<unsigned char>(output)});
                values[output] = 0.;
            }
        }
        Close(input);
    }

    // The identity input has the identity as its only output, with any
    // coefficient: the gate rescales terms it does not act on but leaves their
    // support alone.
    bool KeepsIdentitySupport() const
    {
        return offsets[1] == 1 && entries[0].pauli == 0;
    }

  private:
    void Close(unsigned input)
    {
        offsets[input + 1] = static_cast<unsigned short>(entries.size());
        const auto count = offsets[input + 1] - offsets[input];
        maxOutputs = std::max(maxOutputs, static_cast<unsigned char>(count));
        clifford = clifford && count == 1 && std::abs(entries[offsets[input]].coefficient) == 1.;
        if (count == 1 && entries[offsets[input]].pauli == input && entries[offsets[input]].coefficient == 1.)
            unchanged |= uint64_t(1) << input;
    }
};

struct Trig
{
    // The cosine and the sine are written one after the other. The flag
    // stands between them so that they cannot be read back as one vector,
    // which would have to wait for both writes to leave the store buffer.
    double c;
    bool quarter;
    double s;

    explicit Trig(double angle)
    {
        if (!std::isfinite(angle))
            throw std::invalid_argument("Gate angle must be finite");
        quarter = detail::QuarterTurnCosSin(angle, c, s);
    }
};

struct EulerHalfPhaseFactors
{
    double plusC, plusS, minusC, minusS;

    EulerHalfPhaseFactors(double phi, double lambda)
    {
        long long p, l;
        if (detail::TryGetQuarterTurns(phi, p) && detail::TryGetQuarterTurns(lambda, l))
        {
            constexpr double quarter = 1.57079632679489661923132169163975144;
            detail::RotationHalfAngleCosSin(double((p + l) % 8) * quarter, plusC, plusS);
            detail::RotationHalfAngleCosSin(double((p - l) % 8) * quarter, minusC, minusS);
            return;
        }
        double pc, ps, lc, ls;
        detail::RotationHalfAngleCosSin(phi, pc, ps);
        detail::RotationHalfAngleCosSin(lambda, lc, ls);
        // Product identities avoid overflow and lost small phases in phi +/- lambda.
        plusC = pc * lc - ps * ls;
        plusS = ps * lc + pc * ls;
        minusC = pc * lc + ps * ls;
        minusS = ps * lc - pc * ls;
    }
};

using Rotation3 = std::array<std::array<double, 3>, 3>;

inline Rotation3 Multiply(const Rotation3 &a, const Rotation3 &b)
{
    Rotation3 result{};
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 3; ++k)
                result[i][j] += a[i][k] * b[k][j];
    return result;
}

inline Rotation3 AxisRotation(int axis, const Trig &t)
{
    Rotation3 r{};
    r[axis][axis] = 1.;
    const int a = (axis + 1) % 3, b = (axis + 2) % 3;
    r[a][a] = r[b][b] = t.c;
    r[a][b] = t.s;
    r[b][a] = -t.s;
    return r;
}

inline Rotation3 URotation(double theta, double phi, double lambda)
{
    return Multiply(Multiply(AxisRotation(2, Trig(lambda)), AxisRotation(1, Trig(theta))), AxisRotation(2, Trig(phi)));
}

// The same product written out, for the table of a single-qubit gate, which
// is compiled for every distinct gate of a circuit. Every entry has the value
// the two matrix products give it, from the same operations in the same order
// (a sum starts from zero and takes its products one by one); the products
// with the zeros and ones of the three factors, which change no value, are
// left out. Nothing goes through a matrix in memory on the way, so no entry
// is read back in a wider piece than it was written in, which would have to
// wait for the writes to leave the store buffer.
inline Rotation3 URotationWrittenOut(double theta, double phi, double lambda)
{
    const Trig l(lambda), t(theta), p(phi);
    // The rotation about Z by lambda times the rotation about Y by theta.
    const double m00 = l.c * t.c, m01 = l.s, m02 = l.c * -t.s;
    const double m10 = -l.s * t.c, m11 = l.c, m12 = -l.s * -t.s;
    const double m20 = t.s, m22 = t.c;
    // An entry of that product times the rotation about Z by phi.
    const auto pair = [](double first, double firstFactor, double second, double secondFactor) {
        double sum = 0.;
        sum += first * firstFactor;
        sum += second * secondFactor;
        return sum;
    };
    Rotation3 r;
    r[0][0] = pair(m00, p.c, m01, -p.s);
    r[0][1] = pair(m00, p.s, m01, p.c);
    r[0][2] = m02;
    r[1][0] = pair(m10, p.c, m11, -p.s);
    r[1][1] = pair(m10, p.s, m11, p.c);
    r[1][2] = m12;
    r[2][0] = m20 * p.c;
    r[2][1] = m20 * p.s;
    r[2][2] = m22;
    return r;
}

inline unsigned AxisCode(size_t axis)
{
    return axis == 0 ? 1U : axis == 1 ? 3U : 2U;
}

inline size_t CodeAxis(unsigned code)
{
    return code == 1 ? 0 : code == 3 ? 1 : 2;
}

inline detail::LocalPauli DecodeLocal(unsigned code, unsigned qubits)
{
    detail::LocalPauli p;
    for (unsigned q = 0; q < qubits; ++q)
    {
        p.x |= ((code >> (2 * q)) & 1U) << q;
        p.z |= ((code >> (2 * q + 1)) & 1U) << q;
    }
    return p;
}

inline unsigned EncodeLocal(detail::LocalPauli p, unsigned qubits)
{
    unsigned code = 0;
    for (unsigned q = 0; q < qubits; ++q)
        code |= ((p.x >> q) & 1U) << (2 * q) | ((p.z >> q) & 1U) << (2 * q + 1);
    return code;
}

// The cosine and sine of a one-qubit table that is a rotation about one axis
// (rotationAxis), read back from its entries: the Pauli after next along the
// axis order goes to itself with the cosine and to the next one with the
// sine. An entry that is absent is an exact zero.
inline void RotationCosineSine(const LocalTransfer &table, double &cosine, double &sine)
{
    const unsigned last = AxisCode((size_t(table.rotationAxis) + 2) % 3);
    cosine = sine = 0.;
    for (unsigned i = table.offsets[last]; i < table.offsets[last + 1]; ++i)
        (table.entries[i].pauli == last ? cosine : sine) = table.entries[i].coefficient;
}

inline std::shared_ptr<const LocalTransfer> CompileU(double theta, double phi, double lambda)
{
    auto result = std::make_shared<LocalTransfer>(1);
    const auto rotation = URotationWrittenOut(theta, phi, lambda);
    for (unsigned p = 0; p < 4; ++p)
    {
        std::array<double, 4> values{};
        if (!p)
            values[0] = 1.;
        else
            for (size_t q = 0; q < 3; ++q)
                values[AxisCode(q)] = rotation[q][CodeAxis(p)];
        result->Column(p, values);
    }
    return result;
}

// V's Pauli coefficients are ordered I,X,Y,Z. The rotation describes V^dagger P V.
// diagonalDifference[j] is 1-R[j][j], evaluated without subtracting near-equal
// numbers. This preserves O(angle^2) branches even when cos(angle) rounds to 1.
inline std::shared_ptr<const LocalTransfer> CompileControlled(const std::array<std::complex<double>, 4> &v,
                                                              const Rotation3 &rotation,
                                                              const std::array<double, 3> &diagonalDifference)
{
    auto result = std::make_shared<LocalTransfer>(2);
    for (unsigned input = 0; input < 16; ++input)
    {
        std::array<double, 16> values{};
        const unsigned control = input >> 2, target = input & 3;
        if (control == 0 || control == 2)
        {
            if (target == 0)
                values[input] = 1.;
            else
            {
                const size_t j = CodeAxis(target);
                for (size_t i = 0; i < 3; ++i)
                {
                    const unsigned p = AxisCode(i);
                    const double same = i == j ? 1. - 0.5 * diagonalDifference[j] : 0.5 * rotation[i][j];
                    const double other = i == j ? 0.5 * diagonalDifference[j] : -0.5 * rotation[i][j];
                    values[p | (control << 2)] = same;
                    values[p | ((control ^ 2) << 2)] = other;
                }
            }
        }
        else
        {
            const auto left = DecodeLocal(target, 1);
            for (size_t i = 0; i < 4; ++i)
            {
                const auto right = DecodeLocal(i == 0 ? 0 : AxisCode(i - 1), 1);
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

inline std::shared_ptr<const LocalTransfer> CompileCR(int axis, double angle)
{
    const Trig full(angle), half(0.5 * angle);
    std::array<std::complex<double>, 4> v{};
    v[0] = half.c;
    v[axis + 1] = {0., -half.s};
    std::array<double, 3> difference{};
    for (int j = 0; j < 3; ++j)
        if (j != axis)
            difference[j] = full.quarter ? 1. - full.c : 2. * half.s * half.s;
    return CompileControlled(v, AxisRotation(axis, full), difference);
}

inline std::shared_ptr<const LocalTransfer> CompileCP(double angle)
{
    const Trig full(angle), half(0.5 * angle);
    const double d = full.quarter ? 1. - full.c : 2. * half.s * half.s;
    return CompileControlled(
        {std::complex<double>(1. - 0.5 * d, 0.5 * full.s), 0., 0., std::complex<double>(0.5 * d, -0.5 * full.s)},
        AxisRotation(2, full), {d, d, 0.});
}

inline std::shared_ptr<const LocalTransfer> CompileCU(double theta, double phi, double lambda, double gamma)
{
    const Trig t(theta), p(phi), l(lambda), g(gamma);
    double halfC, halfS;
    detail::RotationHalfAngleCosSin(theta, halfC, halfS);
    const EulerHalfPhaseFactors phases(phi, lambda);
    const std::array<double, 3> b{{-halfS * phases.minusS, halfS * phases.minusC, halfC * phases.plusS}};
    const auto phase = std::complex<double>(phases.plusC, phases.plusS) * std::complex<double>(g.c, g.s);
    std::array<std::complex<double>, 4> v{{phase * (halfC * phases.plusC), phase * std::complex<double>(0., -b[0]),
                                           phase * std::complex<double>(0., -b[1]),
                                           phase * std::complex<double>(0., -b[2])}};
    const auto r = URotation(theta, phi, lambda);
    std::array<double, 3> difference{};
    for (size_t j = 0; j < 3; ++j)
        difference[j] = t.quarter && p.quarter && l.quarter
                            ? 1. - r[j][j]
                            : 2. * (b[(j + 1) % 3] * b[(j + 1) % 3] + b[(j + 2) % 3] * b[(j + 2) % 3]);
    return CompileControlled(v, r, difference);
}

// Compile exact dyadic three-qubit unitaries with the local Pauli algebra
// shared with the extended stabilizer. All pair products are combined here,
// never in the hot loop.
inline std::shared_ptr<const LocalTransfer> CompileExact(const detail::LocalPauliSum &u, int qubits, double scale = 1.)
{
    auto result = std::make_shared<LocalTransfer>(qubits);
    for (unsigned input = 0; input < (1U << (2 * qubits)); ++input)
    {
        std::array<double, 64> values{};
        const auto p = DecodeLocal(input, qubits);
        u.ForEachTerm([&](detail::LocalPauli a, std::complex<double> ca) {
            const detail::LocalPauli ap{uint8_t(a.x ^ p.x), uint8_t(a.z ^ p.z)};
            u.ForEachTerm([&](detail::LocalPauli b, std::complex<double> cb) {
                const unsigned phase =
                    detail::LocalPauliSum::ProductPhase(a, p) + detail::LocalPauliSum::ProductPhase(ap, b);
                const auto coefficient = std::conj(ca) * cb * detail::LocalPauliSum::IPower(phase);
                values[EncodeLocal({uint8_t(ap.x ^ b.x), uint8_t(ap.z ^ b.z)}, qubits)] += coefficient.real();
            });
        });
        if (scale != 1.)
            for (auto &value : values)
                value *= scale;
        result->Column(input, values);
    }
    return result;
}

inline std::shared_ptr<const LocalTransfer> CompileCanonicalClifford(const detail::LocalGate &gate, unsigned qubits)
{
    auto result = std::make_shared<LocalTransfer>(qubits);
    for (unsigned input = 0; input < (1U << (2 * qubits)); ++input)
    {
        std::array<double, 64> values{};
        bool negative;
        const auto output = gate.ConjugatePauli(DecodeLocal(input, qubits), negative);
        values[EncodeLocal(output, qubits)] = negative ? -1. : 1.;
        result->Column(input, values);
    }
    return result;
}

// A bounded per-thread cache of compiled tables: thirty-two entries in sets
// of four. The hash of a key names one entry, its home, and a key enters the
// cache there. A repeated gate therefore finds its table with one comparison
// at an address that depends on the hash alone, as in a cache with a single
// entry per hash value. The entry a new key finds at its home is not lost: it
// moves to the place of the oldest entry of its set, so a circuit that
// alternates between a handful of gates keeps all of their tables, where a
// single entry per hash value would compile again at every collision. Entries
// away from home are found through one word per set, a byte of the key hash
// for each entry, which tells which entries can hold a key at all: a miss
// compares no further keys, and a circuit of all-distinct gates records as
// fast as without the sets. Slots hold the table they describe in a member
// named table, which is empty until the slot is used, and move by assignment.
template <class Slot, unsigned SetBits = 3> struct TransferCache
{
    static constexpr unsigned Ways = 4, Sets = 1U << SetBits;

    std::array<Slot, Sets * Ways> slots;
    std::array<uint32_t, Sets> tags{};
    // The age of every entry of each set, a byte each: 0 for the oldest, 3
    // for the newest, stored relative to the ages 0, 1, 2, 3 of an unused set.
    std::array<uint32_t, Sets> ages{};

    // matches(slot) tells whether a used slot holds the wanted key in full.
    // The entry at the home of the key, if it holds the key.
    template <class Matches> Slot *Home(uint64_t hash, const Matches &matches)
    {
        Slot *const home = &slots[HomeOf(hash)];
        return home->table && matches(*home) ? home : nullptr;
    }

    // An entry away from the home of the key that holds the key.
    template <class Matches> Slot *Moved(uint64_t hash, const Matches &matches)
    {
        const size_t home = HomeOf(hash), set = home / Ways;
        // A byte of the difference is zero where an entry has the wanted tag.
        // The marks find every such byte, and at worst a byte above one as
        // well: the keys decide. The home itself is not looked at again.
        const uint32_t difference = tags[set] ^ (TagOf(hash) * UINT32_C(0x01010101));
        uint32_t marks = (difference - UINT32_C(0x01010101)) & ~difference & UINT32_C(0x80808080) &
                         ~(UINT32_C(0x80) << (8 * (home % Ways)));
        for (; marks; marks &= marks - 1)
        {
            Slot &slot =
                slots[set * Ways + (static_cast<size_t>(Clifford::detail::TrailingZero(uint64_t(marks))) >> 3)];
            if (slot.table && matches(slot))
                return &slot;
        }
        return nullptr;
    }

    // Any entry that holds the key.
    template <class Matches> Slot *Find(uint64_t hash, const Matches &matches)
    {
        Slot *const home = Home(hash, matches);
        return home ? home : Moved(hash, matches);
    }

    // The slot a new key of this hash takes: its home. The entry found there
    // keeps its age and takes the place of the oldest entry of the set, which
    // leaves the cache; unused entries are the oldest. When the home holds
    // the oldest entry itself the move changes nothing. No step depends on
    // which of these happens: a circuit of all-distinct gates, where it is a
    // matter of chance, has nothing to predict here.
    Slot &Replace(uint64_t hash)
    {
        const size_t home = HomeOf(hash), set = home / Ways;
        const unsigned way = static_cast<unsigned>(home % Ways);
        Slot *const first = &slots[set * Ways];
        const uint32_t before = ages[set] ^ UINT32_C(0x03020100);
        // Exactly one age is zero.
        const unsigned oldest = static_cast<unsigned>(Clifford::detail::TrailingZero(
                                    uint64_t((before - UINT32_C(0x01010101)) & ~before & UINT32_C(0x80808080)))) >>
                                3;
        first[oldest] = std::move(first[way]);
        // Every entry is one step older; the moved one takes its age along,
        // and the home holds the newest.
        uint32_t after = ((before | UINT32_C(0x80808080)) - UINT32_C(0x01010101)) & UINT32_C(0x03030303);
        after = WithByte(after, oldest, (((before >> (8 * way)) & UINT32_C(0xFF)) - 1U) & 3U);
        ages[set] = WithByte(after, way, 3U) ^ UINT32_C(0x03020100);
        const uint32_t tag = tags[set];
        tags[set] = WithByte(WithByte(tag, oldest, (tag >> (8 * way)) & UINT32_C(0xFF)), way, TagOf(hash));
        return first[way];
    }

  private:
    // The word with one of its bytes replaced.
    static uint32_t WithByte(uint32_t word, unsigned index, uint32_t value)
    {
        return (word & ~(UINT32_C(0xFF) << (8 * index))) | (value << (8 * index));
    }

    // The highest bits of the hash choose the set and the way in it.
    static size_t HomeOf(uint64_t hash)
    {
        return static_cast<size_t>(hash >> (62 - SetBits));
    }

    // The eight bits of the hash below those.
    static uint32_t TagOf(uint64_t hash)
    {
        return static_cast<uint32_t>(hash >> (54 - SetBits)) & UINT32_C(0xFF);
    }
};

// A sum of products of two matrix entries, together with what bounds the
// rounding it carries: the absolute products and the absolute partial sums.
struct RoundedSum
{
    double value = 0., products = 0., partials = 0.;

    void Add(double x, double y)
    {
        const double product = x * y;
        if (product == 0.)
            return;
        const bool first = products == 0.;
        value += product;
        products += std::abs(product);
        // Adding the first product to zero is exact.
        if (!first)
            partials += std::abs(value);
    }

    // What the value can owe to rounding alone, in units of half an epsilon.
    // Each product is rounded once and its two factors are stored numbers,
    // each within half a unit in the last place of the value it stands for;
    // each addition rounds its partial sum once.
    double Rounding() const
    {
        return 3. * products + partials;
    }
};

// What the sums of CompileMatrix need of the Paulis on a given number of
// qubits, by the bits (x, z) of the Pauli at index x * dimension + z.
struct MatrixPauliTables
{
    // The label of the Pauli and the power of i it carries.
    std::array<unsigned char, 64> labels, powers;
    // For every x, the z whose Pauli carries an even power of i, then the
    // others: the first take the real part of a matrix entry and the others
    // the imaginary part.
    std::array<std::array<unsigned char, 8>, 8> order;
    std::array<unsigned char, 8> evens;

    explicit MatrixPauliTables(unsigned qubits) : labels{}, powers{}, order{}, evens{}
    {
        const unsigned dimension = 1U << qubits;
        for (unsigned x = 0; x < dimension; ++x)
        {
            unsigned odd = dimension;
            for (unsigned z = 0; z < dimension; ++z)
            {
                labels[x * dimension + z] = static_cast<unsigned char>(EncodeLocal({uint8_t(x), uint8_t(z)}, qubits));
                powers[x * dimension + z] = static_cast<unsigned char>(Clifford::detail::Popcount(x & z) & 3U);
                if (powers[x * dimension + z] & 1U)
                    order[x][--odd] = static_cast<unsigned char>(z);
                else
                    order[x][evens[x]++] = static_cast<unsigned char>(z);
            }
        }
    }

    // The tables of one, two or three qubits, built once.
    static const MatrixPauliTables &Of(unsigned qubits)
    {
        static const MatrixPauliTables tables[3] = {MatrixPauliTables(1), MatrixPauliTables(2), MatrixPauliTables(3)};
        return tables[qubits - 1];
    }
};

// The transfer P -> M^dagger P M of any matrix, straight from its entries:
//   R[p][q] = Re sum_r phase_q(r) sum_j conj(M(j, r)) phase_p(j^x_p) M(j^x_p, r^x_q) / 2^n.
// A coefficient that vanishes in exact arithmetic comes out of this sum as
// the rounding of the terms that cancelled, a few 1e-17 for a unitary matrix,
// and as a table entry it would make every term branch, including the terms
// the gate does not act on. A coefficient is therefore dropped when, and only
// when, it does not exceed the rounding its own sum can carry (RoundedSum).
// Nothing else is changed: no entry is rescaled or pinned, so a matrix that
// is not unitary keeps its actual action, and a coefficient made of products
// that do not cancel (an entry of 1e-200, say) always survives, because the
// bound is then a few epsilons of that coefficient itself.
//
// A qubit whose basis value no entry of the matrix changes (a control, or
// every qubit of a diagonal matrix) costs no sums for the inputs with Z or Y
// on it. With Z_s the Z factors of an input on such qubits, the column of
// Z_s P is the column of P with Z_s multiplied into every output, and a sign
// for the outputs whose own Z factor meets a Y of the input: the sums of the
// two columns have the same terms in the same order up to that sign, so the
// copy is the column the sums would give, bit for bit.
template <class Matrix> std::shared_ptr<const LocalTransfer> CompileMatrix(const Matrix &matrix, unsigned qubits)
{
    const unsigned dimension = 1U << qubits, count = dimension * dimension;
    std::array<double, 64> real, imaginary;
    // The nonzero columns of every row: sparse matrices, controlled ones
    // above all, cost what their entries cost.
    std::array<std::array<unsigned char, 8>, 8> columns;
    std::array<unsigned char, 8> width{};
    // The basis bits that some entry changes between its column and its row.
    unsigned moved = 0;
    for (unsigned row = 0; row < dimension; ++row)
        for (unsigned col = 0; col < dimension; ++col)
        {
            const std::complex<double> value = matrix(row, col);
            real[row * dimension + col] = value.real();
            imaginary[row * dimension + col] = value.imag();
            if (value != std::complex<double>(0.))
            {
                columns[row][width[row]++] = static_cast<unsigned char>(col);
                moved |= row ^ col;
            }
        }
    // The label bits of Z on the qubits that no entry changes.
    unsigned fixedZ = 0;
    for (unsigned q = 0; q < qubits; ++q)
        if (!((moved >> q) & 1U))
            fixedZ |= 2U << (2 * q);
    const MatrixPauliTables &paulis = MatrixPauliTables::Of(qubits);
    const auto &labels = paulis.labels;
    const auto &powers = paulis.powers;
    const auto &order = paulis.order;
    const auto &evens = paulis.evens;
    constexpr double half = .5 * std::numeric_limits<double>::epsilon();
    auto result = std::make_shared<LocalTransfer>(qubits);
    // Real and imaginary parts of the entries of M^dagger P M, row major.
    std::array<RoundedSum, 64> re{}, im{};
    std::array<double, 64> values{}, sums{}, rounding{};
    for (unsigned input = 0; input < count; ++input)
    {
        if (const unsigned shift = input & fixedZ)
        {
            // The column of the input without these Z factors is already in
            // the table: its label is smaller. The sign is that of the
            // outputs with Z where the input has Y, among the fixed qubits.
            const unsigned source = input ^ shift, negated = shift & (input << 1);
            uint64_t outputs = 0;
            for (unsigned i = result->offsets[source]; i < result->offsets[source + 1]; ++i)
            {
                const LocalTransfer::Entry entry = result->entries[i];
                const unsigned output = entry.pauli ^ shift;
                values[output] =
                    (Clifford::detail::Popcount(entry.pauli & negated) & 1U) ? -entry.coefficient : entry.coefficient;
                outputs |= uint64_t(1) << output;
            }
            result->TakeColumn(input, values, outputs);
            continue;
        }
        const auto p = DecodeLocal(input, qubits);
        uint64_t entries = 0;
        for (unsigned j = 0; j < dimension; ++j)
        {
            // P|source> = i^turn |j>.
            const unsigned source = j ^ p.x;
            const unsigned turn =
                (powers[p.x * dimension + p.z] + 2U * (Clifford::detail::Popcount(source & p.z) & 1U)) & 3U;
            for (unsigned ri = 0; ri < width[j]; ++ri)
            {
                const unsigned r = columns[j][ri];
                const double ar = real[j * dimension + r], ai = imaginary[j * dimension + r];
                for (unsigned ci = 0; ci < width[source]; ++ci)
                {
                    const unsigned c = columns[source][ci];
                    const double br = real[source * dimension + c], bi = imaginary[source * dimension + c];
                    // conj(a) b = (ar br + ai bi) + i (ar bi - ai br), times i^turn.
                    RoundedSum &x = re[r * dimension + c], &y = im[r * dimension + c];
                    entries |= uint64_t(1) << (r * dimension + c);
                    switch (turn)
                    {
                    case 0:
                        x.Add(ar, br);
                        x.Add(ai, bi);
                        y.Add(ar, bi);
                        y.Add(-ai, br);
                        break;
                    case 1:
                        x.Add(-ar, bi);
                        x.Add(ai, br);
                        y.Add(ar, br);
                        y.Add(ai, bi);
                        break;
                    case 2:
                        x.Add(-ar, br);
                        x.Add(-ai, bi);
                        y.Add(-ar, bi);
                        y.Add(ai, br);
                        break;
                    default:
                        x.Add(ar, bi);
                        x.Add(-ai, br);
                        y.Add(-ar, br);
                        y.Add(-ai, bi);
                        break;
                    }
                }
            }
        }
        // The output with bits (x, z) takes the real part of i^turn times the
        // entry (r, r ^ x), for every row r. Entries are visited row by row,
        // so each output adds its rows in increasing order. A part of an
        // entry without products is skipped with all the outputs it feeds.
        uint64_t outputs = 0;
        for (uint64_t rest = entries; rest; rest &= rest - 1)
        {
            const unsigned index = static_cast<unsigned>(Clifford::detail::TrailingZero(rest));
            const unsigned r = index / dimension, x = r ^ (index % dimension);
            for (unsigned odd = 0; odd < 2; ++odd)
            {
                RoundedSum &term = odd ? im[index] : re[index];
                if (term.products == 0.)
                    continue;
                const double value = term.value, bound = term.Rounding();
                term = RoundedSum();
                for (unsigned k = odd ? evens[x] : 0U, end = odd ? dimension : evens[x]; k < end; ++k)
                {
                    const unsigned z = order[x][k];
                    const unsigned turn =
                        (powers[x * dimension + z] + 2U * (Clifford::detail::Popcount(r & z) & 1U)) & 3U;
                    const unsigned output = labels[x * dimension + z];
                    const uint64_t bit = uint64_t(1) << output;
                    sums[output] += (turn == 0 || turn == 3) ? value : -value;
                    rounding[output] += bound;
                    // Adding the first row to zero is exact.
                    if (outputs & bit)
                        rounding[output] += std::abs(sums[output]);
                    outputs |= bit;
                }
            }
        }
        // The identity input gives the Pauli coefficients of M^dagger M.
        if (!input)
            for (unsigned output = 0; output < count; ++output)
                if (!(std::abs(sums[output] / dimension - (output ? 0. : 1.)) <= 1e-12))
                    result->nonunitary = true;
        for (uint64_t rest = outputs; rest; rest &= rest - 1)
        {
            const unsigned output = static_cast<unsigned>(Clifford::detail::TrailingZero(rest));
            if (std::abs(sums[output]) > half * rounding[output])
                values[output] = sums[output] / dimension;
            sums[output] = rounding[output] = 0.;
        }
        result->TakeColumn(input, values, outputs);
    }
    return result;
}

// Compile I + |s><s| (V-I) without expanding the control projectors into
// quadratic Pauli sums. The only changed control blocks are the selected
// diagonal block, or two off-diagonal blocks. v is V, the controls are the
// bits of mask, their selected values those of selected, and target is the
// one-qubit transfer of V. Without differences, 1 - target[p][p] is formed
// from the table; a verified unitary supplies it without cancellation.
inline std::shared_ptr<const LocalTransfer> CompileControlledTarget(const Eigen::Matrix2cd &v, unsigned qubits,
                                                                    unsigned mask, unsigned selected,
                                                                    const LocalTransfer &targetTransfer,
                                                                    const std::array<double, 4> *diagonalDifferences)
{
    unsigned target = 0;
    while ((mask >> target) & 1U)
        ++target;
    // The Pauli coefficients of V by local label: I, X, Z, Y.
    const std::array<std::complex<double>, 4> coefficients{{(v(0, 0) + v(1, 1)) / 2., (v(0, 1) + v(1, 0)) / 2.,
                                                            (v(0, 0) - v(1, 1)) / 2.,
                                                            std::complex<double>(0., 1.) * (v(0, 1) - v(1, 0)) / 2.}};
    const unsigned controls = qubits - 1, count = 1U << controls;
    const auto compress = [&](unsigned bits) {
        return (bits & ((1U << target) - 1)) | ((bits >> (target + 1)) << target);
    };
    const auto expand = [&](unsigned bits) {
        return (bits & ((1U << target) - 1)) | ((bits >> target) << (target + 1));
    };
    const unsigned active = compress(selected);
    // The power of i in the phase of P(x, z) on a control basis state.
    const auto power = [](unsigned x, unsigned z, unsigned column) {
        return Clifford::detail::Popcount(x & z) + 2U * (Clifford::detail::Popcount(z & column) & 1U);
    };
    // The real part of i^turn times value: multiplying by a power of i only
    // exchanges and negates the components.
    const auto turned = [](unsigned turn, std::complex<double> value) {
        switch (turn & 3U)
        {
        case 0:
            return value.real();
        case 1:
            return -value.imag();
        case 2:
            return -value.real();
        default:
            return value.imag();
        }
    };
    // The output label of every (control X pattern, control Z pattern, target
    // Pauli) triple.
    std::array<unsigned char, 64> labels;
    for (unsigned x = 0; x < count; ++x)
        for (unsigned controlZ = 0; controlZ < count; ++controlZ)
            for (unsigned q = 0; q < 4; ++q)
            {
                const detail::LocalPauli output{uint8_t(expand(x) | ((q & 1U) << target)),
                                                uint8_t(expand(controlZ) | ((q >> 1) << target))};
                labels[(x * count + controlZ) * 4 + q] = static_cast<unsigned char>(EncodeLocal(output, qubits));
            }
    // (V^dagger P - P) and (P V - P) by target Pauli P and output Pauli: the
    // target factors of the two off-diagonal control blocks.
    std::array<std::array<std::complex<double>, 4>, 4> left{}, right{};
    for (unsigned local = 0; local < 4; ++local)
    {
        const auto t = DecodeLocal(local, 1);
        for (unsigned q = 0; q < 4; ++q)
        {
            const auto a = DecodeLocal(q, 1);
            const unsigned output = EncodeLocal({uint8_t(t.x ^ a.x), uint8_t(t.z ^ a.z)}, 1);
            left[local][output] =
                std::conj(coefficients[q]) * detail::LocalPauliSum::IPower(detail::LocalPauliSum::ProductPhase(a, t));
            right[local][output] =
                coefficients[q] * detail::LocalPauliSum::IPower(detail::LocalPauliSum::ProductPhase(t, a));
        }
        left[local][local] -= 1.;
        right[local][local] -= 1.;
    }
    // (V^dagger P V - P) by target Pauli and output Pauli, for the selected block.
    std::array<std::array<double, 4>, 4> diagonal{};
    for (unsigned local = 0; local < 4; ++local)
    {
        for (unsigned i = targetTransfer.offsets[local]; i < targetTransfer.offsets[local + 1]; ++i)
            diagonal[local][targetTransfer.entries[i].pauli] = targetTransfer.entries[i].coefficient;
        if (diagonalDifferences)
            diagonal[local][local] = -(*diagonalDifferences)[local];
        else
            diagonal[local][local] -= 1.;
    }
    auto result = std::make_shared<LocalTransfer>(qubits);
    std::array<double, 64> values{};
    for (unsigned input = 0; input < (1U << (2 * qubits)); ++input)
    {
        values[input] = 1.;
        uint64_t outputs = uint64_t(1) << input;
        const auto p = DecodeLocal(input, qubits);
        const unsigned x = compress(p.x), z = compress(p.z);
        const unsigned local = ((p.x >> target) & 1U) | (((p.z >> target) & 1U) << 1);
        const unsigned char *const label = &labels[x * count * 4];
        if (!x)
        {
            for (unsigned controlZ = 0; controlZ < count; ++controlZ)
            {
                const bool negative = Clifford::detail::Popcount(active & (z ^ controlZ)) & 1U;
                for (unsigned q = 0; q < 4; ++q)
                {
                    const double correction = (negative ? -diagonal[local][q] : diagonal[local][q]) / count;
                    values[label[4 * controlZ + q]] += correction;
                    outputs |= uint64_t(1) << label[4 * controlZ + q];
                }
            }
        }
        else
        {
            const unsigned leftBase = power(x, z, active ^ x), rightBase = power(x, z, active);
            for (unsigned controlZ = 0; controlZ < count; ++controlZ)
            {
                const unsigned leftTurn = leftBase + power(x, controlZ, active);
                const unsigned rightTurn = rightBase + power(x, controlZ, active ^ x);
                for (unsigned q = 0; q < 4; ++q)
                {
                    const double correction =
                        (turned(leftTurn, left[local][q]) + turned(rightTurn, right[local][q])) / count;
                    values[label[4 * controlZ + q]] += correction;
                    outputs |= uint64_t(1) << label[4 * controlZ + q];
                }
            }
        }
        result->TakeColumn(input, values, outputs);
    }
    return result;
}

// The transfer of a matrix that is the identity except for one single-qubit
// block V selected by the other bits (structure). An exact Clifford V, or
// the analytic target of a verified unitary, keeps the identity transfer
// exact; any other matrix is compiled from its entries.
template <class Matrix>
std::shared_ptr<const LocalTransfer> CompileControlledBlock(
    const Matrix &matrix, unsigned qubits, const Gates::GateStructure &structure,
    const std::shared_ptr<const LocalTransfer> &analyticTarget = {},
    const std::array<double, 4> *diagonalDifferences = nullptr)
{
    const unsigned mask = structure.controlMask;
    unsigned target = 0;
    while ((mask >> target) & 1U)
        ++target;
    const unsigned selected = structure.controlValue;
    Eigen::Matrix2cd v;
    for (unsigned row = 0; row < 2; ++row)
        for (unsigned col = 0; col < 2; ++col)
            v(row, col) = matrix(selected | (row << target), selected | (col << target));
    if (analyticTarget)
        return CompileControlledTarget(v, qubits, mask, selected, *analyticTarget, diagonalDifferences);
    if (const auto *canonical = detail::FindCanonicalClifford(v))
        return CompileControlledTarget(v, qubits, mask, selected, *CompileCanonicalClifford(*canonical, 1), nullptr);
    return CompileMatrix(matrix, qubits);
}

// A Pauli-product rotation has at most two outputs per input Pauli. Computing
// sin(theta) directly retains small branches and makes exact quarter turns
// genuine signed permutations, without numerical coefficient pruning.
inline std::shared_ptr<const LocalTransfer> CompilePairRotation(char first, char second, double angle)
{
    const Trig trig(angle);
    const auto axis = detail::PairPauli(first, second);
    auto result = std::make_shared<LocalTransfer>(2);
    for (unsigned input = 0; input < 16; ++input)
    {
        std::array<double, 16> values{};
        const auto p = DecodeLocal(input, 2);
        if ((Clifford::detail::Popcount((axis.x & p.z) ^ (axis.z & p.x)) & 1U) == 0)
            values[input] = 1.;
        else
        {
            values[input] = trig.c;
            const auto phase = detail::LocalPauliSum::IPower(1U + detail::LocalPauliSum::ProductPhase(axis, p));
            values[EncodeLocal({uint8_t(axis.x ^ p.x), uint8_t(axis.z ^ p.z)}, 2)] = phase.real() * trig.s;
        }
        result->Column(input, values);
    }
    return result;
}

// Closed-form XX+YY adjoint columns. Algebraic cancellations and the identity
// column are resolved before evaluating floating-point coefficients. Squared
// sine factors retain genuine second-order branches near zero. XX-YY is the
// conjugation of XX+YY by X on local bit zero.
inline std::shared_ptr<const LocalTransfer> CompileXXYY(bool plus, double theta, double beta)
{
    const Trig full(theta), phase(beta);
    double c, s;
    detail::RotationHalfAngleCosSin(theta, c, s);
    const double sb = phase.s, cb = phase.c, s2 = s * s;
    const double a = s * sb, b = s * cb, d = .5 * full.s;
    const double e = d * sb, f = d * cb, p = s2 * sb * sb, q = s2 * cb * cb, r = s2 * cb * sb;
    const double notP = p > .5 ? c * c + q : 1. - p;
    const double notQ = q > .5 ? c * c + p : 1. - q;
    auto result = std::make_shared<LocalTransfer>(2);
    for (unsigned input = 0; input < 16; ++input)
    {
        std::array<double, 16> v{};
        switch (input)
        {
        case 0:
        case 10:
            v[input] = 1.;
            break;
        case 1:
            v[1] = c;
            v[6] = -a;
            v[14] = b;
            break;
        case 2:
            v[2] = c * c;
            v[5] = e;
            v[7] = f;
            v[8] = s2;
            v[13] = -f;
            v[15] = e;
            break;
        case 3:
            v[3] = c;
            v[6] = -b;
            v[14] = -a;
            break;
        case 4:
            v[4] = c;
            v[9] = a;
            v[11] = b;
            break;
        case 5:
            v[2] = -e;
            v[5] = notP;
            v[7] = -r;
            v[8] = e;
            v[13] = r;
            v[15] = -p;
            break;
        case 6:
            v[1] = a;
            v[3] = b;
            v[6] = c;
            break;
        case 7:
            v[2] = -f;
            v[5] = -r;
            v[7] = notQ;
            v[8] = f;
            v[13] = q;
            v[15] = -r;
            break;
        case 8:
            v[2] = s2;
            v[5] = -e;
            v[7] = -f;
            v[8] = c * c;
            v[13] = f;
            v[15] = -e;
            break;
        case 9:
            v[4] = -a;
            v[9] = c;
            v[12] = b;
            break;
        case 11:
            v[4] = -b;
            v[11] = c;
            v[12] = -a;
            break;
        case 12:
            v[9] = -b;
            v[11] = a;
            v[12] = c;
            break;
        case 13:
            v[2] = f;
            v[5] = r;
            v[7] = q;
            v[8] = -f;
            v[13] = notQ;
            v[15] = r;
            break;
        case 14:
            v[1] = -b;
            v[3] = a;
            v[14] = c;
            break;
        case 15:
            v[2] = -e;
            v[5] = -p;
            v[7] = -r;
            v[8] = e;
            v[13] = r;
            v[15] = notP;
            break;
        }
        if (!plus)
            for (unsigned output = 0; output < 16; ++output)
                if (((input ^ output) >> 1) & 1U)
                    v[output] = -v[output];
        result->Column(input, v);
    }
    return result;
}

// The real Pauli transfer of a unitary inverse is the transpose. Gathering the
// sparse columns neither evaluates nor negates any trigonometry again.
inline std::shared_ptr<const LocalTransfer> TransposeTransfer(const LocalTransfer &forward)
{
    auto inverse = std::make_shared<LocalTransfer>(forward.qubits);
    const unsigned count = 1U << (2 * forward.qubits);
    for (unsigned input = 0; input < count; ++input)
    {
        std::array<double, 64> values{};
        for (unsigned output = 0; output < count; ++output)
            for (unsigned entry = forward.offsets[output]; entry < forward.offsets[output + 1]; ++entry)
                if (forward.entries[entry].pauli == input)
                    values[output] = forward.entries[entry].coefficient;
        inverse->Column(input, values);
    }
    // cos is even and sin is odd, bit for bit: the transpose is the rotation
    // by the opposite angle.
    if (forward.rotationAxis >= 0)
    {
        inverse->rotationAxis = forward.rotationAxis;
        inverse->rotationAngle = -forward.rotationAngle;
    }
    return inverse;
}

// The one-qubit transfers of the fixed Clifford constructors, compiled once
// from their constructor matrices. Empty for any other kind.
inline std::shared_ptr<const LocalTransfer> FixedCliffordTarget(Gates::UnitaryGateDescriptor::Kind kind)
{
    using Kind = Gates::UnitaryGateDescriptor::Kind;
    static const auto tables = [] {
        std::array<std::shared_ptr<const LocalTransfer>, static_cast<size_t>(Kind::Splitter) + 1> result;
        for (Kind fixed : {Kind::Identity, Kind::H, Kind::Hy, Kind::S, Kind::Sdg, Kind::X, Kind::Y, Kind::Z, Kind::Sx,
                           Kind::Sxdg, Kind::Splitter})
        {
            Gates::UnitaryGateDescriptor descriptor;
            descriptor.kind = fixed;
            if (const auto *canonical = detail::FindCanonicalClifford(Gates::UnitaryGateMatrix(descriptor)))
                result[static_cast<size_t>(fixed)] = CompileCanonicalClifford(*canonical, 1);
        }
        return result;
    }();
    return static_cast<size_t>(kind) < tables.size() ? tables[static_cast<size_t>(kind)] : nullptr;
}

// Analytic target actions preserve exact identity and small diagonal
// differences. v is the verified single-qubit matrix of the forward candidate
// (descriptor.adjoint is not looked at here) and controls the number of
// controls above it; the target is local bit 0. Empty for other candidates.
inline std::shared_ptr<const LocalTransfer> CompileVerifiedSingle(const Eigen::Matrix2cd &v, unsigned controls,
                                                                  const Gates::UnitaryGateDescriptor &descriptor)
{
    using Kind = Gates::UnitaryGateDescriptor::Kind;
    std::shared_ptr<const LocalTransfer> target;
    std::array<double, 4> difference{};
    switch (descriptor.kind)
    {
    case Kind::Identity:
    case Kind::H:
    case Kind::Hy:
    case Kind::S:
    case Kind::Sdg:
    case Kind::X:
    case Kind::Y:
    case Kind::Z:
    case Kind::Sx:
    case Kind::Sxdg:
    case Kind::Splitter:
        target = FixedCliffordTarget(descriptor.kind);
        if (!target)
            return {};
        break;
    case Kind::U:
        if (const auto *canonical = detail::FindCanonicalClifford(v))
            target = CompileCanonicalClifford(*canonical, 1);
        else
        {
            const auto &a = descriptor.angles;
            target = CompileU(a[0], a[1], a[2]);
            double c, s;
            detail::RotationHalfAngleCosSin(a[0], c, s);
            const EulerHalfPhaseFactors phases(a[1], a[2]);
            const std::array<double, 3> b{{-s * phases.minusS, s * phases.minusC, c * phases.plusS}};
            for (unsigned j = 0; j < 3; ++j)
                difference[AxisCode(j)] = 2. * (b[(j + 1) % 3] * b[(j + 1) % 3] + b[(j + 2) % 3] * b[(j + 2) % 3]);
        }
        break;
    case Kind::Rx:
    case Kind::Ry:
    case Kind::Rz:
    case Kind::Phase:
    case Kind::T:
    case Kind::Tdg: {
        // An encoded quarter turn gets exact cosines and sines, which makes the
        // rotation below the signed permutation of that Clifford: the matrix
        // need not be searched for among the Cliffords.
        const int axis = descriptor.kind == Kind::Rx ? 0 : descriptor.kind == Kind::Ry ? 1 : 2;
        const double angle = descriptor.kind == Kind::T     ? std::acos(-1.) / 4.
                             : descriptor.kind == Kind::Tdg ? -std::acos(-1.) / 4.
                                                            : descriptor.angles[0];
        const Trig trig(angle);
        const auto rotation = AxisRotation(axis, trig);
        auto result = std::make_shared<LocalTransfer>(1);
        for (unsigned input = 0; input < 4; ++input)
        {
            std::array<double, 4> values{};
            if (!input)
                values[0] = 1.;
            else
                for (unsigned row = 0; row < 3; ++row)
                    values[AxisCode(row)] = rotation[row][CodeAxis(input)];
            result->Column(input, values);
        }
        result->rotationAxis = static_cast<signed char>(axis);
        result->rotationAngle = angle;
        target = result;
        double c, s;
        detail::RotationHalfAngleCosSin(angle, c, s);
        for (unsigned j = 0; j < 3; ++j)
            if (j != unsigned(axis))
                difference[AxisCode(j)] = trig.quarter ? 1. - trig.c : 2. * s * s;
        break;
    }
    default:
        return {};
    }
    if (!controls)
        return target;
    if (target->clifford)
        for (unsigned input = 1; input < 4; ++input)
        {
            difference[input] = 1.;
            for (unsigned entry = target->offsets[input]; entry < target->offsets[input + 1]; ++entry)
                if (target->entries[entry].pauli == input)
                    difference[input] -= target->entries[entry].coefficient;
        }
    const unsigned mask = ((1U << controls) - 1U) << 1;
    return CompileControlledTarget(v, controls + 1, mask, mask, *target, &difference);
}

// The same for a candidate that may be an adjoint: v is then the matrix of
// the adjoint, and the forward transfer is transposed.
inline std::shared_ptr<const LocalTransfer> CompileVerifiedControlled(const Eigen::Matrix2cd &v, unsigned controls,
                                                                      const Gates::UnitaryGateDescriptor &descriptor)
{
    if (!descriptor.adjoint)
        return CompileVerifiedSingle(v, controls, descriptor);
    const Eigen::Matrix2cd forwardMatrix = v.adjoint();
    const auto forward = CompileVerifiedSingle(forwardMatrix, controls, descriptor);
    return forward ? TransposeTransfer(*forward) : nullptr;
}

// The caller has checked the complete current gate matrix against its
// constructor descriptor; arbitrary/mutated matrices never enter.
template <class Matrix>
std::shared_ptr<const LocalTransfer> CompileVerifiedUnitary(const Matrix &matrix, unsigned qubits,
                                                            const Gates::UnitaryGateDescriptor &descriptor)
{
    using Kind = Gates::UnitaryGateDescriptor::Kind;
    if (descriptor.adjoint)
    {
        auto forwardDescriptor = descriptor;
        forwardDescriptor.adjoint = false;
        const Gates::UnitaryLocalMatrix forwardMatrix = matrix.adjoint();
        const auto forward = CompileVerifiedUnitary(forwardMatrix, qubits, forwardDescriptor);
        return forward ? TransposeTransfer(*forward) : nullptr;
    }
    if (descriptor.kind == Kind::Fredkin)
        return CompileExact(detail::MatrixPauliSum(matrix, qubits), qubits);
    if (descriptor.kind == Kind::PairRotation && !descriptor.controls)
        return CompilePairRotation(descriptor.first, descriptor.second, descriptor.angles[0]);
    if ((descriptor.kind == Kind::XXPlusYY || descriptor.kind == Kind::XXMinusYY) && !descriptor.controls)
        return CompileXXYY(descriptor.kind == Kind::XXPlusYY, descriptor.angles[0], descriptor.angles[1]);
    // A verified single-qubit candidate is the identity outside its last two
    // rows and columns, which hold the single-qubit matrix.
    const unsigned dimension = 1U << qubits;
    if (descriptor.controls > 2 || dimension != (2U << descriptor.controls))
        return {};
    Eigen::Matrix2cd v;
    v << matrix(dimension - 2, dimension - 2), matrix(dimension - 2, dimension - 1),
        matrix(dimension - 1, dimension - 2), matrix(dimension - 1, dimension - 1);
    return CompileVerifiedSingle(v, descriptor.controls, descriptor);
}

// A compiled table as its cache hands it to a recorder, with the operation
// the recorder chose for it when it first saw it: a native kernel with the
// same action, or the table itself. A repeated gate is then recorded from the
// cache entry alone, without looking at the table again. The choice is
// forgotten whenever the entry takes another table.
struct CachedTransfer
{
    enum class Recording : unsigned char
    {
        Unknown,
        // The table itself, as the local operation of its number of qubits.
        Table,
        // A rotation about one axis, by the angle, cosine and sine below.
        Rotation,
        // A Clifford with a native kernel on one qubit, on the two qubits of
        // the table in their order, or on the two in exchanged order.
        Single,
        Pair,
        ExchangedPair
    };

    std::shared_ptr<const LocalTransfer> table;
    Recording recording = Recording::Unknown;
    OperationType operation = OperationType::X;
    // The table is that of a matrix that is not unitary (see LocalTransfer).
    bool nonunitary = false;
    // Of a rotation: whether it is an encoded quarter turn, and its numbers.
    bool quarter = false;
    double angle = 0., cosine = 0., sine = 0.;

    void Take(std::shared_ptr<const LocalTransfer> compiled)
    {
        table = std::move(compiled);
        recording = Recording::Unknown;
    }
};

// The value bits of the entries of a matrix are the key of its table. They
// are taken row by row, the real part of an entry before its imaginary part,
// and conjugation changes the sign of zero without changing the operator, so
// a negative zero counts as zero. The functions below hash the entries of a
// matrix on the given number of qubits, compare them with stored words and
// store them.
//
// The hash runs in four lanes, each taking every fourth entry, so that it is
// known a few multiplications after the last entry is read, however many
// there are. The comparison reads the matrix itself: a lookup writes nothing,
// and nothing it reads has just been written, so the table of a repeated gate
// is found without waiting on a chain of multiplications or on the store
// buffer.
template <class Matrix>
void MatrixEntryWords(const Matrix &matrix, unsigned qubits, unsigned entry, uint64_t &real, uint64_t &imaginary)
{
    const std::complex<double> value = matrix(entry >> qubits, entry & ((1U << qubits) - 1U));
    const double parts[2] = {value.real(), value.imag()};
    std::memcpy(&real, &parts[0], sizeof(real));
    std::memcpy(&imaginary, &parts[1], sizeof(imaginary));
}

template <class Matrix> uint64_t MatrixEntriesHash(const Matrix &matrix, unsigned qubits, uint64_t seed)
{
    constexpr uint64_t negativeZero = UINT64_C(0x8000000000000000), multiplier = UINT64_C(0x9E3779B97F4A7C15);
    const auto take = [&](unsigned entry, uint64_t &lane) {
        uint64_t real, imaginary;
        MatrixEntryWords(matrix, qubits, entry, real, imaginary);
        if (real == negativeZero)
            real = 0;
        if (imaginary == negativeZero)
            imaginary = 0;
        lane = (lane ^ real) * multiplier + imaginary;
    };
    uint64_t first = seed, second = UINT64_C(0xC2B2AE3D27D4EB4F), third = UINT64_C(0x165667B19E3779F9),
             fourth = UINT64_C(0x27D4EB2F165667C5);
    // Four, sixteen or sixty-four entries.
    for (unsigned entry = 0, count = 1U << (2 * qubits); entry < count; entry += 4)
    {
        take(entry, first);
        take(entry + 1, second);
        take(entry + 2, third);
        take(entry + 3, fourth);
    }
    const uint64_t mixed = (first ^ ((second << 17) | (second >> 47))) + (third ^ ((fourth << 41) | (fourth >> 23)));
    return (mixed ^ (mixed >> 29)) * multiplier;
}

template <class Matrix> bool SameMatrixEntries(const Matrix &matrix, unsigned qubits, const uint64_t *words)
{
    constexpr uint64_t negativeZero = UINT64_C(0x8000000000000000);
    for (unsigned entry = 0, count = 1U << (2 * qubits); entry < count; ++entry)
    {
        uint64_t real, imaginary;
        MatrixEntryWords(matrix, qubits, entry, real, imaginary);
        // A stored word is never a negative zero.
        if (words[2 * entry] != real && !(real == negativeZero && words[2 * entry] == 0))
            return false;
        if (words[2 * entry + 1] != imaginary && !(imaginary == negativeZero && words[2 * entry + 1] == 0))
            return false;
    }
    return true;
}

template <class Matrix> void StoreMatrixEntries(const Matrix &matrix, unsigned qubits, uint64_t *words)
{
    constexpr uint64_t negativeZero = UINT64_C(0x8000000000000000);
    for (unsigned entry = 0, count = 1U << (2 * qubits); entry < count; ++entry)
    {
        uint64_t real, imaginary;
        MatrixEntryWords(matrix, qubits, entry, real, imaginary);
        words[2 * entry] = real == negativeZero ? 0 : real;
        words[2 * entry + 1] = imaginary == negativeZero ? 0 : imaginary;
    }
}

// Where the hash of a matrix starts from: its number of qubits, or the number
// of controls put on it, and the kind of candidate it came with.
inline uint64_t CandidateSeed(unsigned number, const Gates::UnitaryGateDescriptor &descriptor)
{
    return number | (static_cast<uint64_t>(descriptor.kind) << 8) | (static_cast<uint64_t>(descriptor.adjoint) << 16);
}

// The stored key of the matrix cache. Only the words the number of qubits
// uses are written, compared or copied: a one-qubit entry does not touch the
// room of three qubits.
struct MatrixWords
{
    unsigned qubits = 0;
    std::array<uint64_t, 128> words;

    MatrixWords() = default;

    MatrixWords(const MatrixWords &other) : qubits(other.qubits)
    {
        std::copy_n(other.words.begin(), Used(), words.begin());
    }

    MatrixWords &operator=(const MatrixWords &other)
    {
        if (this != &other)
        {
            qubits = other.qubits;
            std::copy_n(other.words.begin(), Used(), words.begin());
        }
        return *this;
    }

    // The number of words in use: none in an entry that holds no matrix yet.
    size_t Used() const
    {
        return qubits ? size_t(2) << (2 * qubits) : 0;
    }
};

// One entry of the matrix cache: the table, the matrix by value bits and the
// candidate it came with.
struct MatrixTransferSlot : CachedTransfer
{
    MatrixWords key;
    Gates::UnitaryGateDescriptor descriptor;
};

// One cache per thread for every matrix type: the key is the matrix content.
inline TransferCache<MatrixTransferSlot> &MatrixTransferCache()
{
    thread_local TransferCache<MatrixTransferSlot> cache;
    return cache;
}

// A matrix that is exactly a*I + b*P for one Pauli string P that moves at least
// one bit: every diagonal entry is a, and row j has one other entry, b times
// the entry of P, in column j ^ x. This is what a rotation about a Pauli
// product looks like when it is given as a plain matrix (Rx, Ry, RXX, RYY,
// RZX, ...). Its transfer follows from a and b alone:
//   Q commuting with P:      Q -> (|a|^2 + |b|^2) Q + 2 Re(conj(a) b) QP
//   Q anticommuting with P:  Q -> (|a|^2 - |b|^2) Q + 2i Im(conj(a) b) QP
// which costs four sums for the whole table, where CompileMatrix forms every
// coefficient from the entries. The rule for a coefficient is the same: it is
// dropped when, and only when, it does not exceed the rounding its own sum can
// carry, and nothing is rescaled. Empty when the matrix has not this form.
template <class Matrix> std::shared_ptr<const LocalTransfer> CompilePauliPair(const Matrix &matrix, unsigned qubits)
{
    using Complex = std::complex<double>;
    const unsigned dimension = 1U << qubits;
    const Complex a = matrix(0, 0);
    // The other entry of row 0 gives the X mask.
    unsigned x = 0;
    for (unsigned col = 1; col < dimension; ++col)
        if (matrix(0, col) != Complex(0.))
        {
            if (x)
                return {};
            x = col;
        }
    if (!x)
        return {};
    // P|k> = i^(|x & z| + 2 |k & z|) |k ^ x>: the entry of row x (k = 0) is b i^|x & z|, and
    // the row that reaches a single bit k tells whether that bit is in the Z mask.
    const Complex reference = matrix(x, 0);
    unsigned z = 0;
    for (unsigned q = 0; q < qubits; ++q)
    {
        const Complex entry = matrix((1U << q) ^ x, 1U << q);
        if (entry == -reference && entry != reference)
            z |= 1U << q;
        else if (entry != reference)
            return {};
    }
    for (unsigned row = 0; row < dimension; ++row)
        for (unsigned col = 0; col < dimension; ++col)
        {
            const Complex entry = matrix(row, col);
            if (col == row)
            {
                if (entry != a)
                    return {};
            }
            else if (col == (row ^ x))
            {
                if (entry != ((Clifford::detail::Popcount(col & z) & 1U) ? -reference : reference))
                    return {};
            }
            else if (entry != Complex(0.))
                return {};
        }
    if (!std::isfinite(a.real()) || !std::isfinite(a.imag()) || !std::isfinite(reference.real()) ||
        !std::isfinite(reference.imag()))
        return {};
    // b = reference / i^|x & z|: a power of i only exchanges and negates the parts.
    Complex b = reference;
    for (unsigned k = Clifford::detail::Popcount(x & z) & 3U; k; --k)
        b = Complex(b.imag(), -b.real());
    constexpr double half = .5 * std::numeric_limits<double>::epsilon();
    const auto kept = [&](const RoundedSum &sum) {
        return std::abs(sum.value) > half * sum.Rounding() ? sum.value : 0.;
    };
    RoundedSum norm, difference, real, imaginary;
    for (double part : {a.real(), a.imag()})
    {
        norm.Add(part, part);
        difference.Add(part, part);
    }
    for (double part : {b.real(), b.imag()})
    {
        norm.Add(part, part);
        difference.Add(-part, part);
    }
    real.Add(a.real(), b.real());
    real.Add(a.imag(), b.imag());
    imaginary.Add(a.real(), b.imag());
    imaginary.Add(-a.imag(), b.real());
    const double same = kept(norm), opposite = kept(difference), symmetric = 2. * kept(real),
                 skew = 2. * kept(imaginary);
    auto result = std::make_shared<LocalTransfer>(qubits);
    // The identity input gives the Pauli coefficients of M^dagger M.
    result->nonunitary = !(std::abs(norm.value - 1.) <= 1e-12) || !(2. * std::abs(real.value) <= 1e-12);
    const detail::LocalPauli axis{uint8_t(x), uint8_t(z)};
    std::array<double, 64> values{};
    for (unsigned input = 0; input < dimension * dimension; ++input)
    {
        const auto pauli = DecodeLocal(input, qubits);
        const unsigned output = EncodeLocal({uint8_t(pauli.x ^ x), uint8_t(pauli.z ^ z)}, qubits);
        // Q P = i^k (the Pauli of the output label); k is even exactly when Q commutes with P.
        const unsigned k = detail::LocalPauliSum::ProductPhase(pauli, axis);
        if (k & 1U)
        {
            values[input] = opposite;
            values[output] = k == 1 ? -skew : skew;
        }
        else
        {
            values[input] = same;
            values[output] = k == 0 ? symmetric : -symmetric;
        }
        result->TakeColumn(input, values, (uint64_t(1) << input) | (uint64_t(1) << output));
    }
    return result;
}

// Compile the table of a matrix the cache does not hold, and keep it. This is
// where a matrix that is not finite is refused: such a matrix is never kept,
// so a lookup cannot find it.
template <class Matrix>
CachedTransfer &CompileMatrixTransfer(const Matrix &matrix, unsigned qubits,
                                      const Gates::UnitaryGateDescriptor &descriptor, uint64_t hash)
{
    const auto key = detail::MatrixKey(matrix);
    std::shared_ptr<const LocalTransfer> table;
    if (descriptor && Gates::VerifyUnitaryGateDescriptor(matrix, descriptor))
        table = CompileVerifiedUnitary(matrix, qubits, descriptor);
    if (!table)
    {
        if (const auto *canonical = detail::FindCanonicalClifford(matrix, key))
            table = CompileCanonicalClifford(*canonical, qubits);
        else if (auto pair = CompilePauliPair(matrix, qubits))
            table = std::move(pair);
        else
        {
            const auto structure = Gates::ClassifyGateMatrix(matrix);
            table = qubits > 1 && structure.controls == qubits - 1 ? CompileControlledBlock(matrix, qubits, structure)
                                                                   : CompileMatrix(matrix, qubits);
        }
    }
    auto &slot = MatrixTransferCache().Replace(hash);
    slot.key.qubits = qubits;
    StoreMatrixEntries(matrix, qubits, slot.key.words.data());
    slot.descriptor = descriptor;
    slot.Take(std::move(table));
    return slot;
}

// The table of a gate given by its matrix on one to three qubits, local bit 0
// first. A candidate descriptor is used only if the matrix is exactly the one
// its constructor builds; any other matrix, a mutated typed gate included, is
// compiled from its entries. Repeated matrices share one immutable table.
// This form returns the cache entry itself, for a recorder that may not need
// to own the table: the reference holds until this thread compiles again.
template <class Matrix>
CachedTransfer &CachedMatrixTransfer(const Matrix &matrix, unsigned qubits,
                                     const Gates::UnitaryGateDescriptor &descriptor = {})
{
    if (!qubits || qubits > 3 || matrix.rows() != (1U << qubits) || matrix.cols() != matrix.rows())
        throw std::invalid_argument("Gate requires a finite square matrix on one to three qubits");

    // The kind of the candidate is hashed with the entries: a typed gate and
    // the same matrix given without a candidate, which are two entries, then
    // have homes of their own.
    const uint64_t hash = MatrixEntriesHash(matrix, qubits, CandidateSeed(qubits, descriptor));
    if (auto *slot = MatrixTransferCache().Find(hash, [&](const MatrixTransferSlot &candidate) {
            return candidate.key.qubits == qubits && SameMatrixEntries(matrix, qubits, candidate.key.words.data()) &&
                   candidate.descriptor == descriptor;
        }))
        return *slot;
    return CompileMatrixTransfer(matrix, qubits, descriptor, hash);
}

template <class Matrix>
std::shared_ptr<const LocalTransfer> MatrixTransfer(const Matrix &matrix, unsigned qubits,
                                                    const Gates::UnitaryGateDescriptor &descriptor = {})
{
    return CachedMatrixTransfer(matrix, qubits, descriptor).table;
}

// One entry of the controlled cache: the table, the single-qubit matrix by
// value bits, its candidate and the number of controls.
struct ControlledTransferSlot : CachedTransfer
{
    std::array<uint64_t, 8> bits{};
    Gates::UnitaryGateDescriptor descriptor;
    unsigned controls = 0;
};

inline TransferCache<ControlledTransferSlot> &ControlledTransferCache()
{
    thread_local TransferCache<ControlledTransferSlot> cache;
    return cache;
}

// The table of a single-qubit matrix under one or two controls: the target is
// local bit 0 and the controls follow, as in the matrix of ControlledGate and
// DoublyControlledGate. It is the table MatrixTransfer gives for that matrix,
// keyed by the four entries of the single-qubit matrix, so that recording
// such a gate neither builds nor hashes the controlled matrix. This form
// returns the cache entry, which holds until this thread compiles again.
template <class Matrix>
CachedTransfer &CachedControlledTransfer(const Matrix &single, unsigned controls,
                                         const Gates::UnitaryGateDescriptor &descriptor = {})
{
    if (single.rows() != 2 || single.cols() != 2)
        throw std::invalid_argument("Controlled gate requires a finite 2 by 2 matrix");
    if (controls != 1 && controls != 2)
        throw std::invalid_argument("A controlled local gate has one or two controls");
    const uint64_t hash = MatrixEntriesHash(single, 1, CandidateSeed(controls, descriptor));
    auto &cache = ControlledTransferCache();
    if (auto *slot = cache.Find(hash, [&](const ControlledTransferSlot &candidate) {
            return candidate.controls == controls && SameMatrixEntries(single, 1, candidate.bits.data()) &&
                   candidate.descriptor == descriptor;
        }))
        return *slot;
    // A matrix that is not finite is never kept, so it always arrives here.
    if (!single.allFinite())
        throw std::invalid_argument("Controlled gate requires a finite 2 by 2 matrix");
    Eigen::Matrix2cd v;
    v << single(0, 0), single(0, 1), single(1, 0), single(1, 1);
    std::shared_ptr<const LocalTransfer> table;
    if (descriptor && !descriptor.controls && Gates::VerifyUnitaryGateDescriptor(v, descriptor))
        table = CompileVerifiedControlled(v, controls, descriptor);
    if (!table)
    {
        const unsigned dimension = 2U << controls;
        Gates::UnitaryLocalMatrix matrix = Gates::UnitaryLocalMatrix::Identity(dimension, dimension);
        matrix.block(dimension - 2, dimension - 2, 2, 2) = v;
        table = MatrixTransfer(matrix, controls + 1);
    }
    auto &slot = cache.Replace(hash);
    StoreMatrixEntries(single, 1, slot.bits.data());
    slot.descriptor = descriptor;
    slot.controls = controls;
    slot.Take(std::move(table));
    return slot;
}

template <class Matrix>
std::shared_ptr<const LocalTransfer> ControlledTransfer(const Matrix &single, unsigned controls,
                                                        const Gates::UnitaryGateDescriptor &descriptor = {})
{
    return CachedControlledTransfer(single, controls, descriptor).table;
}

inline std::shared_ptr<const LocalTransfer> ECRTransfer()
{
    static const auto table =
        CompileExact(detail::LocalPauliSum::Pauli(0, 'X') +
                         detail::LocalPauliSum::Pauli(0, 'Y', -1.) * detail::LocalPauliSum::Pauli(1, 'X'),
                     2, .5);
    return table;
}

inline std::shared_ptr<const LocalTransfer> DCXTransfer()
{
    using Sum = detail::LocalPauliSum;
    static const auto table =
        CompileExact((Sum::Projector(1, false) + Sum::Projector(1, true) * Sum::Pauli(0, 'X')) *
                         (Sum::Projector(0, false) + Sum::Projector(0, true) * Sum::Pauli(1, 'X')),
                     2);
    return table;
}

// One entry of the cache of two-qubit interactions: the two Pauli letters (or
// the sign of an XX+YY, XX-YY rotation) and the exact bits of the two angles.
struct InteractionSlot
{
    char first = 0, second = 0;
    uint64_t angle = 0, beta = 0;
    std::shared_ptr<const LocalTransfer> table;
};

// The key of an interaction, which also tells whether a slot holds it.
struct InteractionKey
{
    char first, second;
    uint64_t angle, beta, hash;

    InteractionKey(char one, char two, double theta, double phase) : first(one), second(two)
    {
        std::memcpy(&angle, &theta, sizeof(angle));
        std::memcpy(&beta, &phase, sizeof(beta));
        const uint64_t product =
            (angle ^ (uint64_t(first) << 8) ^ uint64_t(second)) * UINT64_C(0x9E3779B97F4A7C15) ^ beta;
        // The fold lets angles in arithmetic progression spread over the cache.
        hash = (product ^ (product >> 29)) * UINT64_C(0x9E3779B97F4A7C15);
    }

    bool operator()(const InteractionSlot &slot) const
    {
        return slot.first == first && slot.second == second && slot.angle == angle && slot.beta == beta;
    }
};

// An interaction that the entry at the home of its key does not hold: look
// among the other entries, then compile and keep the table.
inline std::shared_ptr<const LocalTransfer> InteractionTransferMiss(TransferCache<InteractionSlot> &cache,
                                                                    const InteractionKey &key, double angle,
                                                                    double beta)
{
    if (const InteractionSlot *found = cache.Moved(key.hash, key))
        return found->table;
    auto table = key.first == '+' || key.first == '-' ? CompileXXYY(key.first == '+', angle, beta)
                                                      : CompilePairRotation(key.first, key.second, angle);
    InteractionSlot &slot = cache.Replace(key.hash);
    slot.first = key.first;
    slot.second = key.second;
    slot.angle = key.angle;
    slot.beta = key.beta;
    slot.table = table;
    return table;
}

// Small per-thread cache for repeated variational angles. Exact bit keys keep
// nearby rotations distinct, and tables are immutable across recorded gates.
inline std::shared_ptr<const LocalTransfer> InteractionTransfer(char first, char second, double angle, double beta = 0.)
{
    thread_local TransferCache<InteractionSlot> cache;
    if (angle == 0.)
        angle = 0.;
    if (beta == 0.)
        beta = 0.;
    const InteractionKey key(first, second, angle, beta);
    if (const InteractionSlot *home = cache.Home(key.hash, key))
        return home->table;
    return InteractionTransferMiss(cache, key, angle, beta);
}

inline std::shared_ptr<const LocalTransfer> FixedTransfer(OperationType type)
{
    constexpr double halfPi = 1.57079632679489661923132169163975144;
    using Sum = detail::LocalPauliSum;
    switch (type)
    {
    case OperationType::CS: {
        static const auto map = CompileCP(halfPi);
        return map;
    }
    case OperationType::CSDAG: {
        static const auto map = CompileCP(-halfPi);
        return map;
    }
    case OperationType::CSX:
    case OperationType::CSXDAG: {
        const auto build = [halfPi](bool inverse) {
            const double sign = inverse ? -1. : 1.;
            return CompileControlled(
                {std::complex<double>(.5, .5 * sign), std::complex<double>(.5, -.5 * sign), 0., 0.},
                AxisRotation(0, Trig(sign * halfPi)), {0., 1., 1.});
        };
        static const auto forward = build(false), inverse = build(true);
        return type == OperationType::CSX ? forward : inverse;
    }
    case OperationType::CH: {
        static const auto map =
            CompileControlled({0., std::sqrt(.5), 0., std::sqrt(.5)},
                              Rotation3{{{{0., 0., 1.}}, {{0., -1., 0.}}, {{1., 0., 0.}}}}, {1., 2., 1.});
        return map;
    }
    case OperationType::CCX: {
        static const auto map = CompileExact(Sum::Identity() + Sum::Projector(1, true) * Sum::Projector(2, true) *
                                                                   (Sum::Pauli(0, 'X') + Sum::Identity(-1.)),
                                             3);
        return map;
    }
    case OperationType::CSWAP: {
        static const auto map = [] {
            auto swap = Sum::Identity(.5);
            for (char axis : {'X', 'Y', 'Z'})
                swap = swap + Sum::Pauli(0, axis, .5) * Sum::Pauli(1, axis);
            return CompileExact(Sum::Projector(2, false) + Sum::Projector(2, true) * swap, 3);
        }();
        return map;
    }
    default:
        throw std::invalid_argument("Not a fixed local Pauli gate");
    }
}

// The table of a parameterized gate from its angles.
inline std::shared_ptr<const LocalTransfer> CompileParameterized(OperationType type, double a, double b, double c,
                                                                 double d)
{
    switch (type)
    {
    case OperationType::U:
        return CompileU(a, b, c);
    case OperationType::CU:
        return CompileCU(a, b, c, d);
    case OperationType::CRX:
        return CompileCR(0, a);
    case OperationType::CRY:
        return CompileCR(1, a);
    case OperationType::CRZ:
        return CompileCR(2, a);
    case OperationType::CP:
        return CompileCP(a);
    default:
        throw std::invalid_argument("Not a parameterized local Pauli gate");
    }
}

// One entry of the cache of parameterized gates: the gate and the exact bits
// of its four angles.
struct ParameterizedSlot
{
    OperationType type = OperationType::X;
    std::array<uint64_t, 4> bits{};
    std::shared_ptr<const LocalTransfer> table;
};

// The key of a parameterized gate, which also tells whether a slot holds it.
// The words stay apart: they are hashed and compared one by one.
struct ParameterizedKey
{
    OperationType type;
    uint64_t a, b, c, d, hash;

    ParameterizedKey(OperationType gate, double first, double second, double third, double fourth) : type(gate)
    {
        std::memcpy(&a, &first, sizeof(a));
        std::memcpy(&b, &second, sizeof(b));
        std::memcpy(&c, &third, sizeof(c));
        std::memcpy(&d, &fourth, sizeof(d));
        // Two products side by side and one after a fold: angles in arithmetic
        // progression, whose bit patterns are one too, still spread over the
        // sets of the cache.
        const uint64_t mixed =
            ((a ^ ((b << 32) | (b >> 32))) * UINT64_C(0x9E3779B97F4A7C15)) ^
            ((c ^ ((d << 32) | (d >> 32)) ^ static_cast<uint64_t>(type)) * UINT64_C(0xC2B2AE3D27D4EB4F));
        hash = (mixed ^ (mixed >> 29)) * UINT64_C(0x9E3779B97F4A7C15);
    }

    bool operator()(const ParameterizedSlot &slot) const
    {
        return slot.type == type && slot.bits[0] == a && slot.bits[1] == b && slot.bits[2] == c && slot.bits[3] == d;
    }

    bool HasNegativeZero() const
    {
        constexpr uint64_t negativeZero = UINT64_C(0x8000000000000000);
        return (a == negativeZero) | (b == negativeZero) | (c == negativeZero) | (d == negativeZero);
    }
};

// A parameterized gate that the entry at the home of its key does not hold:
// look among the other entries, then compile and keep the table. A negative
// zero is the same angle as zero and is never stored: such a key is looked up
// as zero here, so keys without one, nearly all of them, pay one test for
// the rule.
inline std::shared_ptr<const LocalTransfer> ParameterizedTransferMiss(TransferCache<ParameterizedSlot> &cache,
                                                                      OperationType type, double a, double b, double c,
                                                                      double d)
{
    ParameterizedKey key(type, a, b, c, d);
    const bool negativeZero = key.HasNegativeZero();
    if (negativeZero)
        key = ParameterizedKey(type, a == 0. ? 0. : a, b == 0. ? 0. : b, c == 0. ? 0. : c, d == 0. ? 0. : d);
    if (const ParameterizedSlot *found = negativeZero ? cache.Find(key.hash, key) : cache.Moved(key.hash, key))
        return found->table;
    auto table = CompileParameterized(type, a, b, c, d);
    ParameterizedSlot &slot = cache.Replace(key.hash);
    slot.type = type;
    slot.bits = {{key.a, key.b, key.c, key.d}};
    slot.table = table;
    return table;
}

// The table of a parameterized gate; angles not used by the gate must be zero.
// Circuits often repeat angles (QAOA layers, Trotter steps), so each thread
// keeps the last tables it compiled in a small cache, keyed by the exact
// angle bits. A hit shares the immutable table.
inline std::shared_ptr<const LocalTransfer> ParameterizedTransfer(OperationType type, double a, double b = 0.,
                                                                  double c = 0., double d = 0.)
{
    thread_local TransferCache<ParameterizedSlot> cache;
    const ParameterizedKey key(type, a, b, c, d);
    if (const ParameterizedSlot *home = cache.Home(key.hash, key))
        return home->table;
    return ParameterizedTransferMiss(cache, type, a, b, c, d);
}

// ---- Noise channels ----
//
// A channel N(rho) = sum_k E_k rho E_k^dagger acts on an observable with its
// adjoint, P -> sum_k E_k^dagger P E_k. Trace preservation makes the adjoint
// unital: column 0 of every table below is exactly the identity, so a channel
// outside the support of a term leaves it alone and stays out of its cone.

// Eigenvalues of a channel that is diagonal in the Pauli basis, indexed by the
// local label. Entry 0 belongs to the identity and is exactly 1.
using Eigenvalues = std::array<double, 16>;

// Transfer columns of a one- or two-qubit channel: [input label][output label].
using ChannelTransfer = std::array<std::array<double, 16>, 16>;

// Public channel arguments index each qubit as I, X, Y, Z (0..3), with q0 in
// the low base-4 digit. Local labels interleave x,z bits: I=0, X=1, Z=2, Y=3.
inline unsigned LabelFromPauliIndex(unsigned index, unsigned qubits)
{
    static constexpr unsigned char code[4] = {0, 1, 3, 2};
    unsigned label = 0;
    for (unsigned q = 0; q < qubits; ++q)
        label |= unsigned(code[(index >> (2 * q)) & 3U]) << (2 * q);
    return label;
}

inline bool LabelsAnticommute(unsigned a, unsigned b)
{
    constexpr unsigned xBits = 0x55U;
    unsigned parity = (a & (b >> 1) & xBits) ^ ((a >> 1) & b & xBits);
    parity ^= parity >> 4;
    parity ^= parity >> 2;
    return parity & 1U;
}

// lambda_a = 1 - 2 * (probability of the errors anticommuting with a). This
// form stays accurate for weak noise and never uses the identity probability.
// probabilities holds 4^qubits values in the public I, X, Y, Z index order.
template <size_t N> Eigenvalues PauliChannelEigenvalues(const std::array<double, N> &probabilities, unsigned qubits)
{
    const unsigned count = 1U << (2 * qubits);
    double total = 0.;
    for (unsigned index = 0; index < count; ++index)
    {
        const double p = probabilities[index];
        if (!(p >= 0. && p <= 1.))
            throw std::invalid_argument("Pauli channel probabilities must lie in [0, 1]");
        total += p;
    }
    if (!(std::abs(total - 1.) <= 1E-9))
        throw std::invalid_argument("Pauli channel probabilities must sum to one");
    Eigenvalues lambda;
    lambda.fill(1.);
    for (unsigned label = 1; label < count; ++label)
    {
        double flipped = 0.;
        for (unsigned index = 1; index < count; ++index)
            if (LabelsAnticommute(label, LabelFromPauliIndex(index, qubits)))
                flipped += probabilities[index];
        lambda[label] = std::max(-1., 1. - 2. * flipped);
    }
    return lambda;
}

// Row-major matrix of at most two qubits.
using LocalMatrix = std::array<std::complex<double>, 16>;

// The matrix of a local label, with qubit 0 in the low bit of the basis index.
inline LocalMatrix LabelMatrix(unsigned label, unsigned qubits)
{
    const unsigned dim = 1U << qubits;
    LocalMatrix matrix{};
    for (unsigned col = 0; col < dim; ++col)
    {
        unsigned row = col;
        std::complex<double> phase = 1.;
        for (unsigned q = 0; q < qubits; ++q)
        {
            const bool x = (label >> (2 * q)) & 1U, z = (label >> (2 * q + 1)) & 1U, bit = (col >> q) & 1U;
            if (x)
                row ^= 1U << q;
            if (x && z)
                phase *= std::complex<double>(0., bit ? -1. : 1.);
            else if (z && bit)
                phase = -phase;
        }
        matrix[row * dim + col] = phase;
    }
    return matrix;
}

// R[a][b] = Re Tr(P_b sum_k E_k^dagger P_a E_k) / 2^qubits, with the first
// Kraus index bit on q0. Rejects operators that are not finite, not square of
// the right size, or not trace preserving. Entries within rounding of 0 or
// +-1 are snapped: a Pauli or Clifford channel then keeps its exact structure.
template <class Matrix> ChannelTransfer KrausTransfer(const std::vector<Matrix> &kraus, unsigned qubits)
{
    if (kraus.empty())
        throw std::invalid_argument("A channel must contain at least one Kraus operator");
    const unsigned dim = 1U << qubits, count = dim * dim;
    std::vector<LocalMatrix> operators;
    operators.reserve(kraus.size());
    LocalMatrix completeness{};
    for (const auto &E : kraus)
    {
        if (static_cast<size_t>(E.rows()) != dim || static_cast<size_t>(E.cols()) != dim)
            throw std::invalid_argument("All Kraus operators must have the same dimensions");
        LocalMatrix e{};
        for (unsigned r = 0; r < dim; ++r)
            for (unsigned c = 0; c < dim; ++c)
            {
                const std::complex<double> value = E(r, c);
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::invalid_argument("Kraus operators must contain only finite values");
                e[r * dim + c] = value;
            }
        for (unsigned r = 0; r < dim; ++r)
            for (unsigned c = 0; c < dim; ++c)
                for (unsigned k = 0; k < dim; ++k)
                    completeness[r * dim + c] += std::conj(e[k * dim + r]) * e[k * dim + c];
        operators.push_back(e);
    }
    double errorSquared = 0.;
    for (unsigned r = 0; r < dim; ++r)
        for (unsigned c = 0; c < dim; ++c)
            errorSquared += std::norm(completeness[r * dim + c] - std::complex<double>(r == c ? 1. : 0.));
    const double tolerance = 1E-10 * dim;
    if (!std::isfinite(errorSquared) || errorSquared > tolerance * tolerance)
        throw std::invalid_argument("Kraus operators do not define a trace-preserving channel");

    std::array<LocalMatrix, 16> basis;
    for (unsigned label = 0; label < count; ++label)
        basis[label] = LabelMatrix(label, qubits);
    ChannelTransfer transfer{};
    transfer[0][0] = 1.;
    constexpr double snap = 1E-14;
    for (unsigned input = 1; input < count; ++input)
    {
        LocalMatrix evolved{};
        for (const auto &e : operators)
        {
            LocalMatrix product{};
            for (unsigned r = 0; r < dim; ++r)
                for (unsigned c = 0; c < dim; ++c)
                    for (unsigned k = 0; k < dim; ++k)
                        product[r * dim + c] += basis[input][r * dim + k] * e[k * dim + c];
            for (unsigned r = 0; r < dim; ++r)
                for (unsigned c = 0; c < dim; ++c)
                    for (unsigned k = 0; k < dim; ++k)
                        evolved[r * dim + c] += std::conj(e[k * dim + r]) * product[k * dim + c];
        }
        for (unsigned output = 0; output < count; ++output)
        {
            std::complex<double> trace = 0.;
            for (unsigned r = 0; r < dim; ++r)
                for (unsigned k = 0; k < dim; ++k)
                    trace += basis[output][r * dim + k] * evolved[k * dim + r];
            double value = trace.real() / dim;
            if (std::abs(value) <= snap)
                value = 0.;
            else if (std::abs(value - 1.) <= snap)
                value = 1.;
            else if (std::abs(value + 1.) <= snap)
                value = -1.;
            transfer[input][output] = value;
        }
    }
    return transfer;
}

// A channel ready to record: its eigenvalues when the transfer matrix is
// diagonal (a Pauli channel, however it was specified), a local table otherwise.
struct CompiledChannel
{
    bool diagonal = true;
    Eigenvalues eigenvalues;
    std::shared_ptr<const LocalTransfer> table;
};

inline CompiledChannel CompileChannel(const ChannelTransfer &transfer, unsigned qubits)
{
    const unsigned count = 1U << (2 * qubits);
    CompiledChannel result;
    result.eigenvalues.fill(1.);
    for (unsigned input = 0; input < count && result.diagonal; ++input)
        for (unsigned output = 0; output < count; ++output)
            if (input != output && transfer[input][output] != 0.)
            {
                result.diagonal = false;
                break;
            }
    if (result.diagonal)
    {
        for (unsigned input = 1; input < count; ++input)
            result.eigenvalues[input] = transfer[input][input];
        return result;
    }
    auto table = std::make_shared<LocalTransfer>(static_cast<int>(qubits));
    for (unsigned input = 0; input < count; ++input)
        table->Column(input, transfer[input]);
    result.table = std::move(table);
    return result;
}

// X,Y -> coherence X,Y and Z -> shift I + scale Z. This covers amplitude
// damping, its finite-temperature form and thermal relaxation.
inline CompiledChannel CompileRelaxation(double coherence, double scale, double shift)
{
    ChannelTransfer transfer{};
    transfer[0][0] = 1.;
    transfer[1][1] = transfer[3][3] = coherence;
    transfer[2][2] = scale;
    transfer[2][0] = shift;
    return CompileChannel(transfer, 1);
}

} // namespace PauliDetail
} // namespace QC
