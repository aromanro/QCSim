#pragma once

#include "LocalGateExtensions.h"
#include "LocalPauliSum.h"
#include "QuantumGate.h"
#include <memory>
#include <vector>

namespace QC
{
namespace detail
{
// Exact content recognition, not a near-unitarity test. Known constructor
// encodings and canonical Clifford matrices with a common eighth-turn global
// phase are accepted. Noncanonical coefficients retain the generic matrix
// path. Raw values cannot distinguish a scaled matrix that lands bit-for-bit
// on an accepted constructor encoding from that constructor's matrix.
inline std::shared_ptr<const LocalGate> &CanonicalCliffordSlot(size_t index)
{
    // Stable owners retain only the canonical circuits actually used. A cold
    // H/ECR call need not allocate or compile hundreds of unrelated gates.
    static std::array<std::shared_ptr<const LocalGate>, 226> slots;
    return slots[index];
}

template <class Matrix, class Expected, class Compile>
const LocalGate *MatchCanonicalClifford(const Matrix &matrix, const Expected &expected, size_t index, unsigned arity,
                                        Compile compile)
{
    if (!(matrix.array() == expected.array()).all())
        return nullptr;
    auto &slot = CanonicalCliffordSlot(index);
    auto gate = std::atomic_load(&slot);
    if (!gate)
    {
        auto fresh = std::make_shared<LocalGate>(arity);
        compile(*fresh);
        gate = std::move(fresh);
        std::shared_ptr<const LocalGate> empty;
        if (!std::atomic_compare_exchange_strong(&slot, &empty, gate))
            gate = std::move(empty);
    }
    return gate.get();
}

template <class Matrix> const LocalGate *SearchCanonicalClifford(const Matrix &matrix)
{
    constexpr double halfPi = 1.57079632679489661923132169163975144;
    const unsigned arity = matrix.rows() == 2 ? 1 : 2;
    const auto match = [&](const auto &expected, size_t index, auto compile) {
        return MatchCanonicalClifford(matrix, expected, index, arity, compile);
    };
    const LocalGate *found;
    if (arity == 1)
    {
        if ((found = match(Gates::HadamardGate<>().getRawOperatorMatrix(), 3, [](LocalGate &g) { g.H(0); })))
            return found;
        if ((found = match(Gates::HyGate<>().getRawOperatorMatrix(), 4, [](LocalGate &g) { g.Sdg(0).H(0).S(0); })))
            return found;
        if ((found = match(Gates::SGate<>().getRawOperatorMatrix(), 5, [](LocalGate &g) { g.S(0); })))
            return found;
        if ((found = match(Gates::SDGGate<>().getRawOperatorMatrix(), 6, [](LocalGate &g) { g.Sdg(0); })))
            return found;
        if ((found = match(Gates::PauliXGate<>().getRawOperatorMatrix(), 7,
                           [&](LocalGate &g) { g.Rotate(0, 'X', 2 * halfPi); })))
            return found;
        if ((found = match(Gates::PauliYGate<>().getRawOperatorMatrix(), 8,
                           [&](LocalGate &g) { g.Rotate(0, 'Y', 2 * halfPi); })))
            return found;
        if ((found = match(Gates::PauliZGate<>().getRawOperatorMatrix(), 9,
                           [&](LocalGate &g) { g.Rotate(0, 'Z', 2 * halfPi); })))
            return found;
        if ((found = match(Gates::SquareRootNOTGate<>().getRawOperatorMatrix(), 10,
                           [&](LocalGate &g) { g.Rotate(0, 'X', halfPi); })))
            return found;
        if ((found = match(Gates::SquareRootNOTDagGate<>().getRawOperatorMatrix(), 11,
                           [&](LocalGate &g) { g.Rotate(0, 'X', -halfPi); })))
            return found;
        if ((found = match(Gates::SplitterGate<>().getRawOperatorMatrix(), 12,
                           [&](LocalGate &g) { g.Rotate(0, 'X', -halfPi); })))
            return found;
        if ((found = match(Gates::SplitterGate<>().getRawOperatorMatrix().adjoint(), 224,
                           [&](LocalGate &g) { g.Rotate(0, 'X', halfPi); })))
            return found;
        for (int turns = -8; turns <= 8; ++turns)
        {
            const double angle = turns * halfPi;
            const size_t slot = 21 + 4 * (turns + 8);
            if ((found = match(Gates::RxGate<>(angle).getRawOperatorMatrix(), slot,
                               [&](LocalGate &g) { g.Rotate(0, 'X', angle); })))
                return found;
            if ((found = match(Gates::RyGate<>(angle).getRawOperatorMatrix(), slot + 1,
                               [&](LocalGate &g) { g.Rotate(0, 'Y', angle); })))
                return found;
            if ((found = match(Gates::RzGate<>(angle).getRawOperatorMatrix(), slot + 2,
                               [&](LocalGate &g) { g.Rotate(0, 'Z', angle); })))
                return found;
            if ((found = match(Gates::PhaseShiftGate<>(angle).getRawOperatorMatrix(), slot + 3,
                               [&](LocalGate &g) { g.Rotate(0, 'Z', angle); })))
                return found;
        }
    }
    else
    {
        const auto swap = [](LocalGate &g) { g.CX(0, 1).CX(1, 0).CX(0, 1); };
        if ((found = match(Gates::CNOTGate<>().getRawOperatorMatrix(), 13, [](LocalGate &g) { g.CX(0, 1); })))
            return found;
        if ((found = match(Gates::ControlledYGate<>().getRawOperatorMatrix(), 14,
                           [](LocalGate &g) { g.Sdg(0).CX(0, 1).S(0); })))
            return found;
        if ((found = match(Gates::ControlledZGate<>().getRawOperatorMatrix(), 15, [](LocalGate &g) { g.CZ(0, 1); })))
            return found;
        if ((found = match(Gates::SwapGate<>().getRawOperatorMatrix(), 16, swap)))
            return found;
        if ((found = match(Gates::iSwapGate<>().getRawOperatorMatrix(), 17, [&](LocalGate &g) {
                 g.CZ(0, 1);
                 swap(g);
                 g.S(0).S(1);
             })))
            return found;
        if ((found = match(Gates::iSwapDagGate<>().getRawOperatorMatrix(), 18, [&](LocalGate &g) {
                 g.CZ(0, 1);
                 swap(g);
                 g.Sdg(0).Sdg(1);
             })))
            return found;
        if ((found = match(Gates::ECRGate<>().getRawOperatorMatrix(), 19,
                           [&](LocalGate &g) { g.S(0).Rotate(1, 'X', halfPi).CX(1, 0).Rotate(0, 'X', 2 * halfPi); })))
            return found;
        if ((found = match(Gates::DCXGate<>().getRawOperatorMatrix(), 20, [](LocalGate &g) { g.CX(1, 0).CX(0, 1); })))
            return found;
        if ((found = match(Gates::DCXGate<>().getRawOperatorMatrix().adjoint(), 225,
                           [](LocalGate &g) { g.CX(0, 1).CX(1, 0); })))
            return found;
        for (int turns = -8; turns <= 8; ++turns)
        {
            const double angle = turns * halfPi;
            const size_t slot = 89 + 5 * (turns + 8);
            if ((found = match(Gates::RXXGate<>(angle).getRawOperatorMatrix(), slot,
                               [&](LocalGate &g) { AppendPairRotation(g, 'X', 'X', angle); })))
                return found;
            if ((found = match(Gates::RYYGate<>(angle).getRawOperatorMatrix(), slot + 1,
                               [&](LocalGate &g) { AppendPairRotation(g, 'Y', 'Y', angle); })))
                return found;
            if ((found = match(Gates::RZZGate<>(angle).getRawOperatorMatrix(), slot + 2,
                               [&](LocalGate &g) { AppendPairRotation(g, 'Z', 'Z', angle); })))
                return found;
            if ((found = match(Gates::RZXGate<>(angle).getRawOperatorMatrix(), slot + 3,
                               [&](LocalGate &g) { AppendPairRotation(g, 'Z', 'X', angle); })))
                return found;
            if ((found = match(Gates::RXYGate<>(angle).getRawOperatorMatrix(), slot + 4,
                               [&](LocalGate &g) { AppendPairRotation(g, 'X', 'Y', angle); })))
                return found;
        }
        for (bool plus : {false, true})
            for (int t = -2; t <= 2; ++t)
                for (int b = -2; b <= 2; ++b)
                {
                    const double theta = t * (2 * halfPi), beta = b * halfPi;
                    const size_t slot = 174 + (plus ? 25 : 0) + 5 * (t + 2) + (b + 2);
                    const Eigen::MatrixXcd expected = plus ? Gates::XXPlusYYGate<>(theta, beta).getRawOperatorMatrix()
                                                           : Gates::XXMinusYYGate<>(theta, beta).getRawOperatorMatrix();
                    found = match(expected, slot, [&](LocalGate &g) {
                        const unsigned q = plus ? 0 : 1;
                        g.Rotate(q, 'Z', plus ? beta : -beta);
                        AppendPairRotation(g, 'X', 'X', .5 * theta);
                        AppendPairRotation(g, 'Y', 'Y', plus ? .5 * theta : -.5 * theta);
                        g.Rotate(q, 'Z', plus ? -beta : beta);
                    });
                    if (found)
                        return found;
                }
    }
    return nullptr;
}

struct SignedCliffordImage
{
    uint8_t x = 0, z = 0;
    bool negative = false;
};

// Multiplication by powers of i uses component swaps, so even subnormal
// entries participate in exact comparisons without products or squared norms.
inline std::complex<double> CliffordPhase(std::complex<double> value, unsigned phase)
{
    switch (phase & 3U)
    {
    case 0:
        return value;
    case 1:
        return {-value.imag(), value.real()};
    case 2:
        return {-value.real(), -value.imag()};
    default:
        return {value.imag(), -value.real()};
    }
}

// Recover the image of each X/Z generator from U P = Q U. Unlike evaluating
// U P U^dagger, this has no floating accumulation and cannot hide a tiny
// non-Clifford perturbation. The separate exact coefficient guard restricts
// the magnitude to the library's constructor encodings, without a tolerance.
template <class Matrix>
bool CanonicalCliffordImages(const Matrix &matrix, unsigned arity, std::array<SignedCliffordImage, 6> &images)
{
    const unsigned dimension = 1U << arity;
    // Controlled non-Clifford gates often have an identity first column and
    // a different final coefficient. Reject that mismatch before a full scan.
    unsigned first = 0;
    while (first < dimension && matrix(first, 0) == std::complex<double>(0.))
        ++first;
    if (first == dimension)
        return false;
    const std::complex<double> firstValue = matrix(first, 0), lastValue = matrix(dimension - 1, dimension - 1);
    // Pick one global coefficient family. SX, for example, has balanced
    // dyadic components (+/-1 +/-i)/2 rather than pure real/imaginary entries.
    // Mixing families entry by entry would not represent a common phase.
    const bool balanced = firstValue.real() != 0. && firstValue.imag() != 0.;
    const auto matchesAmplitude = [balanced](std::complex<double> value, double amplitude) {
        if (balanced)
            return std::abs(value.real()) == amplitude && std::abs(value.imag()) == amplitude;
        return (value.imag() == 0. && std::abs(value.real()) == amplitude) ||
               (value.real() == 0. && std::abs(value.imag()) == amplitude);
    };
    const double firstAmplitude = firstValue.real() == 0. ? std::abs(firstValue.imag()) : std::abs(firstValue.real());
    if (!matchesAmplitude(firstValue, firstAmplitude) ||
        (lastValue != std::complex<double>(0.) && !matchesAmplitude(lastValue, firstAmplitude)))
        return false;
    unsigned support = 0;
    for (unsigned row = 0; row < dimension; ++row)
        support += matrix(row, 0) != std::complex<double>(0.);
    if (!support || (support & (support - 1)))
        return false;
    // Constructors use both 1/sqrt(d) and sqrt(1/d), which differ by one ULP
    // for d=2 and d=8. Select one encoding for the entire matrix; accepting
    // either value independently per entry would hide genuine perturbations.
    // Test the legacy encoding first, preserving its one-square-root path.
    const double divisor = double(support) * (balanced ? 2. : 1.);
    if (firstAmplitude != 1. / std::sqrt(divisor) && firstAmplitude != std::sqrt(1. / divisor))
        return false;
    for (unsigned col = 0; col < dimension; ++col)
    {
        unsigned count = 0;
        for (unsigned row = 0; row < dimension; ++row)
        {
            const std::complex<double> value = matrix(row, col);
            if (value == std::complex<double>(0.))
                continue;
            if (!matchesAmplitude(value, firstAmplitude))
                return false;
            ++count;
        }
        if (count != support)
            return false;
    }
    for (unsigned generator = 0; generator < 2 * arity; ++generator)
    {
        const bool inputZ = generator >= arity;
        const unsigned bit = 1U << (generator % arity);
        bool found = false;
        for (unsigned x = 0; x < dimension && !found; ++x)
            for (unsigned z = 0; z < dimension && !found; ++z)
                for (unsigned negative = 0; negative < 2 && !found; ++negative)
                {
                    bool equal = true;
                    const unsigned phase = Clifford::detail::Popcount(x & z) + 2 * negative;
                    for (unsigned row = 0; row < dimension && equal; ++row)
                    {
                        const unsigned rowPhase = phase + 2 * (Clifford::detail::Popcount((row ^ x) & z) & 1U);
                        for (unsigned col = 0; col < dimension && equal; ++col)
                        {
                            const auto left =
                                inputZ ? ((col & bit) ? -matrix(row, col) : matrix(row, col)) : matrix(row, col ^ bit);
                            equal = left == CliffordPhase(matrix(row ^ x, col), rowPhase);
                        }
                    }
                    if (equal)
                    {
                        images[generator] = {uint8_t(x), uint8_t(z), bool(negative)};
                        found = true;
                    }
                }
        if (!found)
            return false;
    }
    return true;
}

inline std::shared_ptr<const LocalGate> CompileCliffordImages(unsigned arity, std::array<SignedCliffordImage, 6> images)
{
    uint64_t key = arity;
    for (unsigned i = 0; i < 2 * arity; ++i)
        key = (key << 7) | images[i].x | (uint64_t(images[i].z) << 3) | (uint64_t(images[i].negative) << 6);

    // Reuse equivalent signed maps without retaining the potentially very
    // large three-qubit Clifford group. Callers that retain a compilation hold
    // a shared owner, so evicting either bounded cache remains safe.
    struct Slot
    {
        uint64_t key = 0;
        std::shared_ptr<const LocalGate> gate;
    };

    thread_local std::array<Slot, 64> circuits;
    auto &slot = circuits[(key * UINT64_C(11400714819323198485)) >> 58];
    if (slot.gate && slot.key == key)
        return slot.gate;

    using Op = LocalGate::CliffordOp;
    using Kind = Op::Kind;
    std::vector<Op> reduction;
    const auto append = [&](Kind kind, unsigned q, unsigned other = 0) {
        reduction.push_back({kind, uint8_t(q), uint8_t(other)});
        const unsigned bit = 1U << q, otherBit = 1U << other;
        for (unsigned i = 0; i < 2 * arity; ++i)
        {
            auto &image = images[i];
            const bool x = (image.x & bit) != 0, z = (image.z & bit) != 0;
            switch (kind)
            {
            case Kind::H:
                image.negative ^= x && z;
                if (x != z)
                {
                    image.x ^= bit;
                    image.z ^= bit;
                }
                break;
            case Kind::Sdg:
                image.negative ^= x && !z;
                if (x)
                    image.z ^= bit;
                break;
            case Kind::X:
                image.negative ^= z;
                break;
            case Kind::Z:
                image.negative ^= x;
                break;
            case Kind::CX: {
                const bool controlX = (image.x & otherBit) != 0, controlZ = (image.z & otherBit) != 0;
                image.negative ^= controlX && z && (x == controlZ);
                if (controlX)
                    image.x ^= bit;
                if (z)
                    image.z ^= otherBit;
                break;
            }
            default:
                assert(false);
            }
        }
    };
    for (unsigned q = 0; q < arity; ++q)
    {
        const auto &zImage = images[arity + q];
        unsigned pivot = q;
        while (!(((zImage.x | zImage.z) >> pivot) & 1U))
            ++pivot;
        if (pivot != q)
        {
            append(Kind::CX, q, pivot);
            append(Kind::CX, pivot, q);
            append(Kind::CX, q, pivot);
        }
        for (unsigned j = q; j < arity; ++j)
            if ((zImage.x >> j) & 1U)
            {
                if ((zImage.z >> j) & 1U)
                    append(Kind::Sdg, j);
                append(Kind::H, j);
            }
        for (unsigned j = q + 1; j < arity; ++j)
            if ((zImage.z >> j) & 1U)
                append(Kind::CX, q, j);
        if (zImage.negative)
            append(Kind::X, q);

        const auto &xImage = images[q];
        if ((xImage.z >> q) & 1U)
            append(Kind::Sdg, q);
        for (unsigned j = q + 1; j < arity; ++j)
        {
            if ((xImage.z >> j) & 1U)
            {
                if ((xImage.x >> j) & 1U)
                    append(Kind::Sdg, j);
                else
                    append(Kind::H, j);
            }
            if ((xImage.x >> j) & 1U)
                append(Kind::CX, j, q);
        }
        if (xImage.negative)
            append(Kind::Z, q);
    }
    auto gate = std::make_shared<LocalGate>(arity);
    constexpr double pi = 3.14159265358979323846264338327950288;
    for (auto it = reduction.rbegin(); it != reduction.rend(); ++it)
        switch (it->kind)
        {
        case Kind::H:
            gate->H(it->first);
            break;
        case Kind::Sdg:
            gate->S(it->first);
            break;
        case Kind::X:
            gate->Rotate(it->first, 'X', pi);
            break;
        case Kind::Z:
            gate->Rotate(it->first, 'Z', pi);
            break;
        case Kind::CX:
            gate->CX(it->first, it->second);
            break;
        default:
            assert(false);
        }
    slot = {key, gate};
    return gate;
}

template <class Matrix> std::shared_ptr<const LocalGate> SynthesizeCanonicalClifford(const Matrix &matrix)
{
    const unsigned arity = matrix.rows() == 2 ? 1 : matrix.rows() == 4 ? 2 : 3;
    std::array<SignedCliffordImage, 6> images;
    return CanonicalCliffordImages(matrix, arity, images) ? CompileCliffordImages(arity, images) : nullptr;
}

// The fixed catalog already owns its entries for program lifetime. An alias
// with no control block gives the owned API the same allocation-free hot path.
inline std::shared_ptr<const LocalGate> StaticCliffordReference(const LocalGate *gate)
{
    return {std::shared_ptr<const LocalGate>{}, gate};
}

// One cache per thread for every matrix type: the key holds the contents.
inline LocalMatrixCache<std::shared_ptr<const LocalGate>> &CanonicalCliffordCache()
{
    thread_local LocalMatrixCache<std::shared_ptr<const LocalGate>> cache;
    return cache;
}

template <class Matrix>
const std::shared_ptr<const LocalGate> &CachedCanonicalClifford(const Matrix &matrix, const LocalMatrixKey &key)
{
    auto &cache = CanonicalCliffordCache();
    const unsigned set = cache.Set(key.hash, UINT64_C(0x9E3779B97F4A7C15));
    if (auto *slot = cache.Find(key, set, [](const std::shared_ptr<const LocalGate> &) { return true; }))
        return slot->value;
    std::shared_ptr<const LocalGate> found;
    // Eigen's isIdentity(0) compares squared magnitudes, which can underflow
    // for genuine tiny entries. Recognition requires exact scalar equality.
    bool identity = true;
    for (Eigen::Index row = 0; row < matrix.rows() && identity; ++row)
        for (Eigen::Index col = 0; col < matrix.cols() && identity; ++col)
            identity = matrix(row, col) == std::complex<double>(row == col ? 1. : 0.);
    if (identity)
    {
        if (matrix.rows() == 2)
        {
            static const LocalGate identity(1);
            found = StaticCliffordReference(&identity);
        }
        else if (matrix.rows() == 4)
        {
            static const LocalGate identity(2);
            found = StaticCliffordReference(&identity);
        }
        else
        {
            static const LocalGate identity(3);
            found = StaticCliffordReference(&identity);
        }
    }
    else if (matrix.rows() == 2 && matrix(0, 0) == std::complex<double>(1. / std::sqrt(2.)) &&
             matrix(0, 1) == matrix(0, 0) && matrix(1, 0) == matrix(0, 0) && matrix(1, 1) == -matrix(0, 0))
    {
        static const LocalGate hadamard = [] {
            LocalGate gate(1);
            gate.H(0);
            return gate;
        }();
        found = StaticCliffordReference(&hadamard);
    }
    else
    {
        // The exact synthesis decides. It rejects a matrix that is not a
        // canonical Clifford at its first coefficients, where the catalog
        // would compare it with every named gate. A Clifford that is one of
        // the named gates keeps the catalog's short circuit, which also
        // lives for the whole program; any other gets the synthesized one.
        const unsigned arity = matrix.rows() == 2 ? 1 : matrix.rows() == 4 ? 2 : 3;
        std::array<SignedCliffordImage, 6> images;
        if (CanonicalCliffordImages(matrix, arity, images))
        {
            if (arity < 3)
                found = StaticCliffordReference(SearchCanonicalClifford(matrix));
            if (!found)
                found = CompileCliffordImages(arity, images);
        }
    }
    auto &slot = cache.Insert(key, set);
    slot.value = std::move(found);
    return slot.value;
}

template <class Matrix>
std::shared_ptr<const LocalGate> FindCanonicalCliffordOwned(const Matrix &matrix, const LocalMatrixKey &key)
{
    return CachedCanonicalClifford(matrix, key);
}

template <class Matrix> std::shared_ptr<const LocalGate> FindCanonicalCliffordOwned(const Matrix &matrix)
{
    return FindCanonicalCliffordOwned(matrix, MatrixKey(matrix));
}

// Immediate-use bridge for consumers that finish compiling/replaying before
// another recognition call. A synthesized raw pointer can expire on cache
// eviction; consumers retaining it must use FindCanonicalCliffordOwned.
template <class Matrix> const LocalGate *FindCanonicalClifford(const Matrix &matrix, const LocalMatrixKey &key)
{
    return CachedCanonicalClifford(matrix, key).get();
}

template <class Matrix> const LocalGate *FindCanonicalClifford(const Matrix &matrix)
{
    return CachedCanonicalClifford(matrix, MatrixKey(matrix)).get();
}
} // namespace detail
} // namespace QC
