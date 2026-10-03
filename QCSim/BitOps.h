#pragma once

#include <cassert>
#include <cstdint>
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

}}}
