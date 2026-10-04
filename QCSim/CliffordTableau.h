#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "BitOps.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace QC
{
namespace Clifford
{
namespace detail
{

struct BitReference
{
    Word *word;
    Word mask;

    operator bool() const noexcept
    {
        return (*word & mask) != 0;
    }

    BitReference &operator=(bool value) noexcept
    {
        *word = (*word & ~mask) | (value ? mask : 0);
        return *this;
    }

    BitReference &operator=(const BitReference &other) noexcept
    {
        return *this = bool(other);
    }
};

template <bool Const> struct BitSpan
{
    using Pointer = std::conditional_t<Const, const Word *, Word *>;
    Pointer words;
    size_t bits;

    size_t size() const noexcept
    {
        return bits;
    }

    auto operator[](size_t bit) const noexcept
    {
        assert(bit < bits);
        if constexpr (Const)
            return (words[bit / 64] & (Word(1) << (bit % 64))) != 0;
        else
            return BitReference{words + bit / 64, Word(1) << (bit % 64)};
    }
};

template <bool Const> struct SignReference
{
    std::conditional_t<Const, const uint8_t *, uint8_t *> value;

    operator bool() const noexcept
    {
        return *value != 0;
    }

    SignReference &operator=(bool sign) noexcept
    {
        static_assert(!Const, "Read-only row");
        *value = sign;
        return *this;
    }

    SignReference &operator=(const SignReference &other) noexcept
    {
        return *this = bool(other);
    }

    SignReference &operator^=(bool flip) noexcept
    {
        return *this = (bool(*this) != flip);
    }
};

// Views never own storage and must not outlive their tableau. Rows and signs
// have separate storage so independent rows can be updated by separate threads.
template <bool Const> struct TableauRow
{
    BitSpan<Const> X, Z;
    SignReference<Const> PhaseSign;
    // Copy construction aliases a complete view; assigning its proxy members
    // would rebind the bits but write through the old sign pointer.
    TableauRow &operator=(const TableauRow &) = delete;

    size_t Words() const noexcept
    {
        return X.bits / 64 + (X.bits % 64 != 0);
    }

    size_t GetNrQubits() const noexcept
    {
        return X.bits;
    }

    template <bool C> bool operator==(const TableauRow<C> &other) const noexcept
    {
        if (X.bits != other.X.bits || bool(PhaseSign) != bool(other.PhaseSign))
            return false;
        for (size_t w = 0; w < Words(); ++w)
            if (X.words[w] != other.X.words[w] || Z.words[w] != other.Z.words[w])
                return false;
        return true;
    }

    void Clear() noexcept
    {
        static_assert(!Const, "Read-only row");
        std::fill_n(X.words, Words(), Word(0));
        std::fill_n(Z.words, Words(), Word(0));
        PhaseSign = false;
    }

    template <bool C> void CopyFrom(const TableauRow<C> &source) noexcept
    {
        assert(X.bits == source.X.bits);
        for (size_t w = 0; w < Words(); ++w)
        {
            X.words[w] = source.X.words[w];
            Z.words[w] = source.Z.words[w];
        }
        PhaseSign = bool(source.PhaseSign);
    }

    template <bool C> bool Anticommutes(const TableauRow<C> &other) const noexcept
    {
        Word parity = 0;
        for (size_t w = 0; w < Words(); ++w)
            parity ^= (X.words[w] & other.Z.words[w]) ^ (Z.words[w] & other.X.words[w]);
        return (Popcount(parity) & 1) != 0;
    }

    bool HasX() const noexcept
    {
        for (size_t w = 0; w < Words(); ++w)
            if (X.words[w])
                return true;
        return false;
    }

    // P(x,z) = (-1)^sign i^popcount(x&z) X^x Z^z. The extra phase
    // allows products of anticommuting rows when conjugating inverse images.
    template <bool C> void Multiply(const TableauRow<C> &right, unsigned extraPhase = 0) noexcept
    {
        unsigned phase = 2 * unsigned(bool(PhaseSign) != bool(right.PhaseSign)) + extraPhase;
        for (size_t w = 0; w < Words(); ++w)
        {
            const Word x = X.words[w], z = Z.words[w], rx = right.X.words[w], rz = right.Z.words[w];
            const Word nx = x ^ rx, nz = z ^ rz;
            phase += Popcount(x & z) + Popcount(rx & rz) + 2 * Popcount(z & rx) - Popcount(nx & nz);
            X.words[w] = nx;
            Z.words[w] = nz;
        }
        assert((phase & 1) == 0);
        PhaseSign = (phase & 2) != 0;
    }
};

inline void SwapRows(TableauRow<false> left, TableauRow<false> right) noexcept
{
    assert(left.GetNrQubits() == right.GetNrQubits());
    for (size_t w = 0; w < left.Words(); ++w)
    {
        std::swap(left.X.words[w], right.X.words[w]);
        std::swap(left.Z.words[w], right.Z.words[w]);
    }
    const bool sign = left.PhaseSign;
    left.PhaseSign = bool(right.PhaseSign);
    right.PhaseSign = sign;
}

// All X/Z rows live in one allocation; no per-generator heap buffers.
class PackedTableau
{
  public:
    PackedTableau() = default;

    explicit PackedTableau(size_t rows) : PackedTableau(rows, rows)
    {
    }

    PackedTableau(size_t rows, size_t qubits) : rows(rows), qubits(qubits), words(qubits / 64 + (qubits % 64 != 0))
    {
        if (words > bits.max_size() / 2 || (words && rows > bits.max_size() / (2 * words)))
            throw std::length_error("Clifford tableau exceeds container capacity");
        bits.resize(2 * rows * words);
        signs.resize(rows);
    }

    PackedTableau(const PackedTableau &) = default;

    PackedTableau(PackedTableau &&other) noexcept
    {
        swap(other);
    }

    PackedTableau &operator=(const PackedTableau &other)
    {
        if (this != &other)
        {
            PackedTableau copy(other);
            swap(copy);
        }
        return *this;
    }

    PackedTableau &operator=(PackedTableau &&other) noexcept
    {
        swap(other);
        return *this;
    }

    size_t size() const noexcept
    {
        return rows;
    }

    size_t GetNrQubits() const noexcept
    {
        return qubits;
    }

    bool empty() const noexcept
    {
        return rows == 0;
    }

    bool HasSameShape(const PackedTableau &other) const noexcept
    {
        return rows == other.rows && qubits == other.qubits;
    }

    void CopyFrom(const PackedTableau &other) noexcept
    {
        assert(HasSameShape(other));
        if (this == &other)
            return;
        std::copy(other.bits.begin(), other.bits.end(), bits.begin());
        std::copy(other.signs.begin(), other.signs.end(), signs.begin());
    }

    void clear() noexcept
    {
        PackedTableau empty;
        swap(empty);
    }

    void Clear() noexcept
    {
        std::fill(bits.begin(), bits.end(), Word(0));
        std::fill(signs.begin(), signs.end(), uint8_t(0));
    }

    void swap(PackedTableau &other) noexcept
    {
        std::swap(rows, other.rows);
        std::swap(qubits, other.qubits);
        std::swap(words, other.words);
        bits.swap(other.bits);
        signs.swap(other.signs);
    }

    TableauRow<false> operator[](size_t row) noexcept
    {
        assert(row < rows);
        Word *x = words ? bits.data() + 2 * row * words : nullptr;
        return {{x, qubits}, {words ? x + words : nullptr, qubits}, {signs.data() + row}};
    }

    TableauRow<true> operator[](size_t row) const noexcept
    {
        assert(row < rows);
        const Word *x = words ? bits.data() + 2 * row * words : nullptr;
        return {{x, qubits}, {words ? x + words : nullptr, qubits}, {signs.data() + row}};
    }

    void SwapRows(size_t a, size_t b) noexcept
    {
        if (a == b)
            return;
        for (size_t w = 0; w < 2 * words; ++w)
            std::swap(bits[2 * a * words + w], bits[2 * b * words + w]);
        std::swap(signs[a], signs[b]);
    }

    bool operator==(const PackedTableau &other) const noexcept
    {
        return rows == other.rows && qubits == other.qubits && bits == other.bits && signs == other.signs;
    }

  private:
    size_t rows = 0, qubits = 0, words = 0;
    std::vector<Word> bits;
    std::vector<uint8_t> signs;
};

// Measuring physical Z_q is random exactly when its image has an X part; the
// pivot is the image's lowest logical X bit.
template <bool C> bool FindPivot(const TableauRow<C> &row, size_t &pivot) noexcept
{
    for (size_t w = 0; w < row.Words(); ++w)
        if (row.X.words[w])
        {
            pivot = 64 * w + TrailingZero(row.X.words[w]);
            return true;
        }
    return false;
}

// Multiply the image of the physical one-qubit Pauli ('X', 'Y' or 'Z') on
// qubit q into result, using Y = i X Z.
inline void MultiplyImage(TableauRow<false> result, const PackedTableau &inverseX, const PackedTableau &inverseZ,
                          size_t q, char pauli) noexcept
{
    if (pauli != 'Z')
        result.Multiply(inverseX[q]);
    if (pauli != 'X')
        result.Multiply(inverseZ[q], pauli == 'Y' ? 1 : 0);
}

// The images form a signed symplectic basis exactly when they have the
// canonical Pauli commutation relations.
inline bool IsSymplectic(const PackedTableau &inverseX, const PackedTableau &inverseZ) noexcept
{
    const size_t n = inverseZ.size();
    for (size_t left = 0; left < n; ++left)
    {
        for (size_t right = left; right < n; ++right)
            if (inverseX[left].Anticommutes(inverseX[right]) || inverseZ[left].Anticommutes(inverseZ[right]))
                return false;
        for (size_t right = 0; right < n; ++right)
            if (inverseX[left].Anticommutes(inverseZ[right]) != (left == right))
                return false;
    }
    return true;
}

// Packed rows need far less work than a bit-by-bit kernel. Keep at least 128
// rows per worker below the caller's OpenMP limit, and stay serial below 512.
inline int CollapseWorkers(size_t qubits, bool enableMultithreading) noexcept
{
#ifdef _OPENMP
    if (enableMultithreading && qubits >= 512)
        return static_cast<int>(std::min(qubits / 128, size_t(omp_get_max_threads())));
#else
    (void)qubits;
    (void)enableMultithreading;
#endif
    return 1;
}

// The inverse images x[q] = U^dagger X_q U and z[q] = U^dagger Z_q U of a
// Clifford map U. Appending G to U maps every image P to U^dagger G^dagger P G U,
// so each gate touches a constant number of rows. The map does not own storage;
// callers validate qubits and invalidate their own caches.
class InverseMap
{
  public:
    InverseMap(PackedTableau &inverseX, PackedTableau &inverseZ) noexcept : x(inverseX), z(inverseZ)
    {
        assert(x.size() == z.size());
    }

    void SetIdentity() noexcept
    {
        x.Clear();
        z.Clear();
        for (size_t q = 0; q < z.size(); ++q)
        {
            x[q].X[q] = true;
            z[q].Z[q] = true;
        }
    }

    void ApplyH(size_t q) noexcept
    {
        SwapRows(x[q], z[q]);
    }

    void ApplyS(size_t q) noexcept
    {
        x[q].Multiply(z[q], 3);
    }

    void ApplySdg(size_t q) noexcept
    {
        x[q].Multiply(z[q], 1);
    }

    void ApplyX(size_t q) noexcept
    {
        z[q].PhaseSign ^= true;
    }

    void ApplyY(size_t q) noexcept
    {
        x[q].PhaseSign ^= true;
        z[q].PhaseSign ^= true;
    }

    void ApplyZ(size_t q) noexcept
    {
        x[q].PhaseSign ^= true;
    }

    void ApplySx(size_t q) noexcept
    {
        z[q].Multiply(x[q], 3);
    }

    void ApplySxDag(size_t q) noexcept
    {
        z[q].Multiply(x[q], 1);
    }

    void ApplyK(size_t q) noexcept
    {
        z[q].Multiply(x[q], 3);
        x[q].PhaseSign ^= true;
    }

    void ApplyCX(size_t target, size_t control) noexcept
    {
        x[control].Multiply(x[target]);
        z[target].Multiply(z[control]);
    }

    void ApplyCY(size_t target, size_t control) noexcept
    {
        ApplySdg(target);
        ApplyCX(target, control);
        ApplyS(target);
    }

    void ApplyCZ(size_t target, size_t control) noexcept
    {
        x[target].Multiply(z[control]);
        x[control].Multiply(z[target]);
    }

    void ApplySwap(size_t a, size_t b) noexcept
    {
        SwapRows(x[a], x[b]);
        SwapRows(z[a], z[b]);
    }

    void ApplyISwap(size_t a, size_t b) noexcept
    {
        ApplyS(a);
        ApplyS(b);
        ApplyCZ(a, b);
        ApplySwap(a, b);
    }

    void ApplyISwapDag(size_t a, size_t b) noexcept
    {
        ApplySdg(a);
        ApplySdg(b);
        ApplyCZ(a, b);
        ApplySwap(a, b);
    }

    // One signed pi/2 rotation exp(-+i pi/4 P) about the one-qubit Pauli P
    // ('X', 'Y' or 'Z'), up to global phase. Images of the operators that
    // anticommute with P become i^e image(P) * row = i^(e+2) row * image(P),
    // with e = 1 for the rotation and e = 3 for its inverse. scratch holds
    // image(Y); the X and Z images are rows already.
    void ApplyQuarterTurn(size_t q, char axis, bool inverse, TableauRow<false> scratch) noexcept
    {
        const unsigned extraPhase = inverse ? 1 : 3;
        if (axis == 'X')
            z[q].Multiply(x[q], extraPhase);
        else if (axis == 'Z')
            x[q].Multiply(z[q], extraPhase);
        else
        {
            scratch.Clear();
            MultiplyImage(scratch, x, z, q, 'Y');
            x[q].Multiply(scratch, extraPhase);
            z[q].Multiply(scratch, extraPhase);
        }
    }

    // Rebase after a random measurement of physical Z_qubit, whose image has
    // logical X at pivot. Afterwards that image is logical Z_pivot and the
    // outcome is folded into the signs, so a zero logical input stays zero.
    // With outcome false this is the outcome-independent change of basis
    // used by frames that record outcomes in their component labels.
    // measured is scratch storage of the same width.
    void CollapseZ(size_t qubit, size_t pivot, bool outcome, TableauRow<false> measured, int workers) noexcept
    {
        measured.CopyFrom(z[qubit]);
        unsigned measuredY = 0;
        for (size_t w = 0; w < measured.Words(); ++w)
            measuredY += Popcount(measured.X.words[w] & measured.Z.words[w]);
        const size_t pivotWord = pivot / 64;
        const Word mask = Word(1) << (pivot % 64);
        const auto update = [&](TableauRow<false> row, bool anticommutes) {
            // Only the image of physical X_qubit anticommutes with measured Z.
            // Rows without logical X_p need only their two pivot bits changed.
            if (!row.X[pivot])
            {
                row.X[pivot] = anticommutes;
                row.Z[pivot] = false;
                return;
            }
            const bool sign = bool(row.PhaseSign) ^ bool(measured.PhaseSign) ^ anticommutes ^ outcome;
            unsigned phase = 2 * unsigned(sign) + measuredY;
            for (size_t w = 0; w < row.Words(); ++w)
            {
                const Word rx = row.X.words[w], rz = row.Z.words[w];
                Word nx = rx ^ measured.X.words[w], nz = rz ^ measured.Z.words[w];
                if (w == pivotWord)
                {
                    nx = (nx & ~mask) | (anticommutes ? mask : 0);
                    nz |= mask;
                }
                // Convert between Hermitian-Pauli and ordered X/Z phases;
                // absorb the observed outcome into the new logical-Z signs.
                phase += Popcount(rx & rz) - Popcount(nx & nz) + 2 * Popcount(rx & measured.Z.words[w]);
                row.X.words[w] = nx;
                row.Z.words[w] = nz;
            }
            assert((phase & 1) == 0);
            row.PhaseSign = (phase & 2) != 0;
        };
        const size_t n = z.size();
        if (workers <= 1)
        {
            for (size_t q = 0; q < n; ++q)
            {
                update(x[q], q == qubit);
                update(z[q], false);
            }
            return;
        }
#pragma omp parallel for num_threads(workers)
        for (long long q = 0; q < static_cast<long long>(n); ++q)
        {
            update(x[q], size_t(q) == qubit);
            update(z[q], false);
        }
    }

  private:
    PackedTableau &x;
    PackedTableau &z;
};

} // namespace detail
} // namespace Clifford
} // namespace QC
