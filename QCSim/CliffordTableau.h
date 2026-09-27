#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace QC { namespace Clifford { namespace detail {

using Word = uint64_t;
inline unsigned Popcount(Word value) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && defined(__AVX2__)
    return static_cast<unsigned>(__popcnt64(value));
#elif defined(_MSC_VER) && defined(_M_IX86) && defined(__AVX2__)
    return __popcnt(static_cast<unsigned>(value)) + __popcnt(static_cast<unsigned>(value >> 32));
#elif defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_popcountll(value));
#else
    value -= (value >> 1) & UINT64_C(0x5555555555555555);
    value = (value & UINT64_C(0x3333333333333333)) + ((value >> 2) & UINT64_C(0x3333333333333333));
    value = (value + (value >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f);
    return static_cast<unsigned>((value * UINT64_C(0x0101010101010101)) >> 56);
#endif
}

inline unsigned TrailingZero(Word value) noexcept
{
    assert(value != 0);
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned long bit;
    _BitScanForward64(&bit, value);
    return bit;
#elif defined(_MSC_VER) && defined(_M_IX86)
    unsigned long bit;
    if (_BitScanForward(&bit, static_cast<unsigned long>(value))) return bit;
    _BitScanForward(&bit, static_cast<unsigned long>(value >> 32));
    return bit + 32;
#elif defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_ctzll(value));
#else
    unsigned bit = 0;
    while ((value & 1) == 0) { value >>= 1; ++bit; }
    return bit;
#endif
}

struct BitReference
{
    Word* word;
    Word mask;
    operator bool() const noexcept { return (*word & mask) != 0; }
    BitReference& operator=(bool value) noexcept
    { *word = (*word & ~mask) | (value ? mask : 0); return *this; }
    BitReference& operator=(const BitReference& other) noexcept { return *this = bool(other); }
};

template<bool Const> struct BitSpan
{
    using Pointer = std::conditional_t<Const, const Word*, Word*>;
    Pointer words;
    size_t bits;
    size_t size() const noexcept { return bits; }
    auto operator[](size_t bit) const noexcept
    {
        assert(bit < bits);
        if constexpr (Const) return (words[bit / 64] & (Word(1) << (bit % 64))) != 0;
        else return BitReference{words + bit / 64, Word(1) << (bit % 64)};
    }
};

template<bool Const> struct SignReference
{
    std::conditional_t<Const, const uint8_t*, uint8_t*> value;
    operator bool() const noexcept { return *value != 0; }
    SignReference& operator=(bool sign) noexcept { static_assert(!Const, "Read-only row"); *value = sign; return *this; }
    SignReference& operator=(const SignReference& other) noexcept { return *this = bool(other); }
    SignReference& operator^=(bool flip) noexcept { return *this = (bool(*this) != flip); }
};

// Views never own storage and must not outlive their tableau. Rows and signs
// have separate storage so independent rows can be updated by separate threads.
template<bool Const> struct TableauRow
{
    BitSpan<Const> X, Z;
    SignReference<Const> PhaseSign;
    // Copy construction aliases a complete view; assigning its proxy members
    // would rebind the bits but write through the old sign pointer.
    TableauRow& operator=(const TableauRow&) = delete;
    size_t Words() const noexcept { return X.bits / 64 + (X.bits % 64 != 0); }
    size_t GetNrQubits() const noexcept { return X.bits; }
    template<bool C> bool operator==(const TableauRow<C>& other) const noexcept
    {
        if (X.bits != other.X.bits || bool(PhaseSign) != bool(other.PhaseSign)) return false;
        for (size_t w = 0; w < Words(); ++w)
            if (X.words[w] != other.X.words[w] || Z.words[w] != other.Z.words[w]) return false;
        return true;
    }
    void Clear() noexcept
    {
        static_assert(!Const, "Read-only row");
        std::fill_n(X.words, Words(), Word(0));
        std::fill_n(Z.words, Words(), Word(0));
        PhaseSign = false;
    }
    template<bool C> void CopyFrom(const TableauRow<C>& source) noexcept
    {
        assert(X.bits == source.X.bits);
        for (size_t w = 0; w < Words(); ++w) { X.words[w] = source.X.words[w]; Z.words[w] = source.Z.words[w]; }
        PhaseSign = bool(source.PhaseSign);
    }
    template<bool C> bool Anticommutes(const TableauRow<C>& other) const noexcept
    {
        Word parity = 0;
        for (size_t w = 0; w < Words(); ++w)
            parity ^= (X.words[w] & other.Z.words[w]) ^ (Z.words[w] & other.X.words[w]);
        return (Popcount(parity) & 1) != 0;
    }
    bool HasX() const noexcept
    {
        for (size_t w = 0; w < Words(); ++w) if (X.words[w]) return true;
        return false;
    }
    // P(x,z) = (-1)^sign i^popcount(x&z) X^x Z^z. The extra phase
    // allows products of anticommuting rows when conjugating inverse images.
    template<bool C> void Multiply(const TableauRow<C>& right, unsigned extraPhase = 0) noexcept
    {
        unsigned phase = 2 * unsigned(bool(PhaseSign) != bool(right.PhaseSign)) + extraPhase;
        for (size_t w = 0; w < Words(); ++w)
        {
            const Word x = X.words[w], z = Z.words[w], rx = right.X.words[w], rz = right.Z.words[w];
            const Word nx = x ^ rx, nz = z ^ rz;
            phase += Popcount(x & z) + Popcount(rx & rz) + 2 * Popcount(z & rx) - Popcount(nx & nz);
            X.words[w] = nx; Z.words[w] = nz;
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
    explicit PackedTableau(size_t rows) : PackedTableau(rows, rows) {}
    PackedTableau(size_t rows, size_t qubits) : rows(rows), qubits(qubits), words(qubits / 64 + (qubits % 64 != 0))
    {
        if (words > bits.max_size() / 2 || (words && rows > bits.max_size() / (2 * words)))
            throw std::length_error("Clifford tableau exceeds container capacity");
        bits.resize(2 * rows * words);
        signs.resize(rows);
    }
    PackedTableau(const PackedTableau&) = default;
    PackedTableau(PackedTableau&& other) noexcept { swap(other); }
    PackedTableau& operator=(const PackedTableau& other)
    {
        if (this != &other) { PackedTableau copy(other); swap(copy); }
        return *this;
    }
    PackedTableau& operator=(PackedTableau&& other) noexcept { swap(other); return *this; }
    size_t size() const noexcept { return rows; }
    size_t GetNrQubits() const noexcept { return qubits; }
    bool empty() const noexcept { return rows == 0; }
    bool HasSameShape(const PackedTableau& other) const noexcept
    { return rows == other.rows && qubits == other.qubits; }
    void CopyFrom(const PackedTableau& other) noexcept
    {
        assert(HasSameShape(other));
        if (this == &other) return;
        std::copy(other.bits.begin(), other.bits.end(), bits.begin());
        std::copy(other.signs.begin(), other.signs.end(), signs.begin());
    }
    void clear() noexcept { PackedTableau empty; swap(empty); }
    void Clear() noexcept { std::fill(bits.begin(), bits.end(), Word(0)); std::fill(signs.begin(), signs.end(), uint8_t(0)); }
    void swap(PackedTableau& other) noexcept
    {
        std::swap(rows, other.rows); std::swap(qubits, other.qubits); std::swap(words, other.words);
        bits.swap(other.bits); signs.swap(other.signs);
    }
    TableauRow<false> operator[](size_t row) noexcept
    {
        assert(row < rows);
        Word* x = words ? bits.data() + 2 * row * words : nullptr;
        return {{x, qubits}, {words ? x + words : nullptr, qubits}, {signs.data() + row}};
    }
    TableauRow<true> operator[](size_t row) const noexcept
    {
        assert(row < rows);
        const Word* x = words ? bits.data() + 2 * row * words : nullptr;
        return {{x, qubits}, {words ? x + words : nullptr, qubits}, {signs.data() + row}};
    }
    void SwapRows(size_t a, size_t b) noexcept
    {
        if (a == b) return;
        for (size_t w = 0; w < 2 * words; ++w) std::swap(bits[2 * a * words + w], bits[2 * b * words + w]);
        std::swap(signs[a], signs[b]);
    }
    bool operator==(const PackedTableau& other) const noexcept
    { return rows == other.rows && qubits == other.qubits && bits == other.bits && signs == other.signs; }
private:
    size_t rows = 0, qubits = 0, words = 0;
    std::vector<Word> bits;
    std::vector<uint8_t> signs;
};

}}}
