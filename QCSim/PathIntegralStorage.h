#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#define QC_PATH_INTEGRAL_COMPACT_STORAGE 1

namespace QC
{
namespace PathIntegral
{

// Owning interchange value. Stored amplitude keys use only their active words;
// this fixed-capacity type remains useful for input, output and recursive paths.
struct FastVectorBool
{
    static constexpr size_t MaxWords = 16;
    FastVectorBool() = default;

    explicit FastVectorBool(size_t bits) : nBits(CheckedSize(bits))
    {
    }

    explicit FastVectorBool(const std::vector<bool> &bits) : FastVectorBool(bits.size())
    {
        for (size_t i = 0; i < nBits; ++i)
            if (bits[i])
                set(i, true);
    }

    FastVectorBool(const uint64_t *data, size_t bits) : FastVectorBool(bits)
    {
        std::copy_n(data, nWords(), words.begin());
    }

    bool get(size_t i) const
    {
        assert(i < nBits);
        return ((words[i / 64] >> (i % 64)) & 1) != 0;
    }

    void set(size_t i, bool value)
    {
        assert(i < nBits);
        const uint64_t mask = uint64_t{1} << (i % 64);
        if (value)
            words[i / 64] |= mask;
        else
            words[i / 64] &= ~mask;
    }

    size_t size() const
    {
        return nBits;
    }

    size_t nWords() const
    {
        return (nBits + 63) / 64;
    }

    const std::array<uint64_t, MaxWords> &getWords() const
    {
        return words;
    }

    std::vector<bool> toVector() const
    {
        std::vector<bool> result(nBits);
        for (size_t i = 0; i < nBits; ++i)
            result[i] = get(i);
        return result;
    }

    bool operator==(const FastVectorBool &other) const
    {
        return nBits == other.nBits && std::equal(words.begin(), words.begin() + nWords(), other.words.begin());
    }

    bool operator!=(const FastVectorBool &other) const
    {
        return !(*this == other);
    }

  private:
    static size_t CheckedSize(size_t bits)
    {
        if (bits > MaxWords * 64)
            throw std::length_error("Path integral states support at most 1024 qubits");
        return bits;
    }

    std::array<uint64_t, MaxWords> words{};
    size_t nBits = 0;
};

struct FastVectorBoolHash
{
    size_t operator()(const FastVectorBool &state) const
    {
        size_t seed = 0;
        for (size_t i = state.nWords(); i > 0; --i)
            seed ^= std::hash<uint64_t>{}(state.getWords()[i - 1]) + (seed << 6);
        return seed;
    }
};

// Views returned by iteration are valid until a structural container mutation.
// Convert to FastVectorBool to keep an owning key beyond that point.
class StateView
{
  public:
    StateView(const uint64_t *data, size_t bits) : words(data), nBits(bits)
    {
    }

    bool get(size_t i) const
    {
        assert(i < nBits);
        return ((words[i / 64] >> (i % 64)) & 1) != 0;
    }

    size_t size() const
    {
        return nBits;
    }

    size_t nWords() const
    {
        return (nBits + 63) / 64;
    }

    const uint64_t *getWords() const
    {
        return words;
    }

    operator FastVectorBool() const
    {
        return FastVectorBool(words, nBits);
    }

    std::vector<bool> toVector() const
    {
        std::vector<bool> result(nBits);
        for (size_t i = 0; i < nBits; ++i)
            result[i] = get(i);
        return result;
    }

  private:
    const uint64_t *words;
    size_t nBits;
};

namespace Detail
{
class MutableStateView : public StateView
{
  public:
    MutableStateView(uint64_t *data, size_t bits) : StateView(data, bits), words(data)
    {
    }

    void set(size_t i, bool value)
    {
        assert(i < size());
        const uint64_t mask = uint64_t{1} << (i % 64);
        if (value)
            words[i / 64] |= mask;
        else
            words[i / 64] &= ~mask;
    }

  private:
    uint64_t *words;
};
} // namespace Detail

class PathIntegralSimulator;

// A flat sparse amplitude map. Width is selected on the first insertion after
// clear(), shared by every key, and checked on subsequent insertions. Keys and
// amplitudes are contiguous; buckets hold indices, never separately allocated
// nodes. The common <=64-qubit case stores one uint64_t per key.
class AmplitudeMap
{
    friend class PathIntegralSimulator;
    using Complex = std::complex<double>;
    static constexpr size_t Deleted = std::numeric_limits<size_t>::max();

  public:
    AmplitudeMap() = default;

    AmplitudeMap(const AmplitudeMap &other)
        : values(other.values), keys(other.keys), nBits(other.nBits), wordCount(other.wordCount),
          occupied(other.occupied)
    {
        // A snapshot should not inherit a peak-sized index after pruning.
        const size_t needed = TableSize(size());
        if (Excess(other.buckets.size(), needed))
        {
            buckets.resize(needed);
            Reindex();
        }
        else
            buckets = other.buckets;
    }

    AmplitudeMap(AmplitudeMap &&) noexcept = default;
    AmplitudeMap &operator=(AmplitudeMap &&) noexcept = default;

    AmplitudeMap &operator=(const AmplitudeMap &other)
    {
        if (this != &other)
        {
            AmplitudeMap copy(other);
            swap(copy);
        }
        return *this;
    }

    template <bool Const> class Iterator
    {
        friend class AmplitudeMap;
        using Owner = std::conditional_t<Const, const AmplitudeMap, AmplitudeMap>;
        Owner *owner = nullptr;
        size_t index = 0;

        Iterator(Owner *map, size_t i) : owner(map), index(i)
        {
        }

      public:
        using iterator_category = std::input_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = std::pair<FastVectorBool, Complex>;
        using reference = std::pair<StateView, std::conditional_t<Const, const Complex &, Complex &>>;

        struct Arrow
        {
            reference value;

            const reference *operator->() const
            {
                return &value;
            }
        };

        using pointer = Arrow;
        Iterator() = default;

        reference operator*() const
        {
            return {owner->State(index), owner->values[index]};
        }

        Arrow operator->() const
        {
            return {**this};
        }

        Iterator &operator++()
        {
            ++index;
            return *this;
        }

        Iterator operator++(int)
        {
            auto old = *this;
            ++*this;
            return old;
        }

        bool operator==(const Iterator &other) const
        {
            return owner == other.owner && index == other.index;
        }

        bool operator!=(const Iterator &other) const
        {
            return !(*this == other);
        }
    };

    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;

    size_t size() const
    {
        return values.size();
    }

    bool empty() const
    {
        return values.empty();
    }

    size_t QubitCount() const
    {
        return nBits;
    }

    size_t max_size() const
    {
        return std::min({values.max_size(), keys.max_size() / (wordCount ? wordCount : 16), buckets.max_size() / 2});
    }

    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
        return iterator(this, size());
    }

    const_iterator begin() const
    {
        return const_iterator(this, 0);
    }

    const_iterator end() const
    {
        return const_iterator(this, size());
    }

    const_iterator cbegin() const
    {
        return begin();
    }

    const_iterator cend() const
    {
        return end();
    }

    void clear()
    {
        const size_t previousSize = size(), previousWords = wordCount;
        values.clear();
        keys.clear();
        // Forget the index without touching all its old buckets. A later
        // reserve/insertion initializes only the index it actually needs.
        buckets.clear();
        if (Excess(values.capacity(), previousSize))
            std::vector<Complex>().swap(values);
        if (Excess(keys.capacity() / (previousWords ? previousWords : 1), previousSize))
            std::vector<uint64_t>().swap(keys);
        if (Excess(buckets.capacity(), TableSize(previousSize)))
            std::vector<size_t>().swap(buckets);
        nBits = wordCount = occupied = 0;
    }

    void Release()
    {
        AmplitudeMap{}.swap(*this);
    }

    void swap(AmplitudeMap &other) noexcept
    {
        values.swap(other.values);
        keys.swap(other.keys);
        buckets.swap(other.buckets);
        std::swap(nBits, other.nBits);
        std::swap(wordCount, other.wordCount);
        std::swap(occupied, other.occupied);
    }

    void reserve(size_t count)
    {
        if (count > max_size())
            throw std::length_error("Path integral amplitude capacity exceeded");
        values.reserve(count);
        if (wordCount)
            keys.reserve(count * wordCount);
        const size_t tableSize = TableSize(count);
        if (tableSize > buckets.size())
            Rehash(tableSize);
    }

    template <class Key> Complex &operator[](const Key &key)
    {
        EnsureWidth(key.size());
        if (buckets.empty())
            Rehash(8);
        size_t slot = Slot(key);
        if (buckets[slot] && buckets[slot] != Deleted)
            return values[buckets[slot] - 1];
        if (occupied + 1 > buckets.size() / 2)
        {
            Rehash(size() + 1 > buckets.size() / 2 ? buckets.size() * 2 : buckets.size());
            slot = Slot(key);
        }
        // Reserve both arrays before changing their sizes, keeping insertion
        // consistent if allocation fails. Complex and word copies do not throw.
        if (values.size() == values.capacity() || keys.capacity() < (size() + 1) * wordCount)
        {
            const size_t capacity = std::max<size_t>(8, size() * 2 + 1);
            if (capacity > max_size())
                throw std::length_error("Path integral amplitude capacity exceeded");
            values.reserve(capacity);
            keys.reserve(capacity * wordCount);
        }
        const size_t index = size();
        keys.resize((index + 1) * wordCount);
        for (size_t w = 0; w < wordCount; ++w)
            keys[index * wordCount + w] = key.getWords()[w];
        values.emplace_back(0., 0.);
        if (buckets[slot] == 0)
            ++occupied;
        buckets[slot] = index + 1;
        return values.back();
    }

    template <class Key> iterator find(const Key &key)
    {
        return iterator(this, FindIndex(key));
    }

    template <class Key> const_iterator find(const Key &key) const
    {
        return const_iterator(this, FindIndex(key));
    }

    template <class Key> Complex &at(const Key &key)
    {
        const size_t i = FindIndex(key);
        if (i == size())
            throw std::out_of_range("Missing path integral amplitude");
        return values[i];
    }

    template <class Key> const Complex &at(const Key &key) const
    {
        const size_t i = FindIndex(key);
        if (i == size())
            throw std::out_of_range("Missing path integral amplitude");
        return values[i];
    }

    iterator erase(iterator position)
    {
        assert(position.owner == this && position.index < size());
        const size_t index = position.index, last = size() - 1;
        buckets[Slot(State(index))] = Deleted;
        if (index != last)
        {
            const size_t lastSlot = Slot(State(last));
            std::copy_n(keys.data() + last * wordCount, wordCount, keys.data() + index * wordCount);
            values[index] = values[last];
            buckets[lastSlot] = index + 1;
        }
        values.pop_back();
        keys.resize(size() * wordCount);
        return iterator(this, index);
    }

  private:
    static bool Excess(size_t capacity, size_t needed)
    {
        return capacity / 8 >= std::max<size_t>(8, needed);
    }

    void CompactIfSparse()
    {
        if (Excess(values.capacity(), size()) || Excess(keys.capacity() / (wordCount ? wordCount : 1), size()) ||
            Excess(buckets.capacity(), TableSize(size())))
        {
            // Reclaiming spare capacity is optional. Allocation pressure must
            // not turn a completed gate/collapse into a reported failure.
            try
            {
                AmplitudeMap compact(*this);
                swap(compact);
            }
            catch (const std::bad_alloc &)
            {
                // The original arrays and index remain intact.
            }
        }
    }

    void ReleaseIfOversized(size_t count, size_t width)
    {
        const size_t words = std::max<size_t>(1, (width + 63) / 64);
        if (Excess(values.capacity(), count) || Excess(keys.capacity() / words, count) ||
            Excess(buckets.capacity(), TableSize(count)))
            Release();
    }

    void PrepareOutput(size_t width, size_t count)
    {
        if (width > FastVectorBool::MaxWords * 64)
            throw std::length_error("Path integral states support at most 1024 qubits");
        const size_t words = std::max<size_t>(1, (width + 63) / 64);
        if (count > std::min({values.max_size(), keys.max_size() / words, buckets.max_size() / 2}))
            throw std::length_error("Path integral amplitude capacity exceeded");
        const size_t neededTable = TableSize(count);
        const bool releaseTable = Excess(buckets.capacity(), neededTable);
        const size_t tableSize = releaseTable ? neededTable : std::max(neededTable, buckets.size());
        // Size scratch storage against the next output, not its previous
        // population (or the zero population immediately after clear()).
        values.clear();
        keys.clear();
        buckets.clear();
        nBits = width;
        wordCount = words;
        occupied = 0;
        if (Excess(values.capacity(), count))
            std::vector<Complex>().swap(values);
        if (Excess(keys.capacity() / words, count))
            std::vector<uint64_t>().swap(keys);
        if (releaseTable)
            std::vector<size_t>().swap(buckets);
        values.reserve(count);
        buckets.resize(tableSize);
        keys.reserve(count * words);
    }

    template <class Key> void InsertUnique(const Key &key, const Complex &amplitude)
    {
        // Only used for distinct output groups, in a freshly prepared map.
        assert(key.size() == nBits && occupied == size());
        if (size() + 1 > buckets.size() / 2)
            Rehash(buckets.size() * 2);
        if (values.size() == values.capacity() || keys.capacity() < (size() + 1) * wordCount)
        {
            const size_t capacity = std::max<size_t>(8, size() * 2 + 1);
            if (capacity > max_size())
                throw std::length_error("Path integral amplitude capacity exceeded");
            values.reserve(capacity);
            keys.reserve(capacity * wordCount);
        }
        const size_t index = size();
        keys.resize((index + 1) * wordCount);
        std::copy_n(key.getWords(), wordCount, keys.data() + index * wordCount);
        values.push_back(amplitude);
        size_t slot = Hash(key) & (buckets.size() - 1);
        while (buckets[slot])
            slot = (slot + 1) & (buckets.size() - 1);
        buckets[slot] = index + 1;
        ++occupied;
    }

    StateView State(size_t index) const
    {
        return StateView(keys.data() + index * wordCount, nBits);
    }

    static uint64_t Mix(uint64_t x)
    {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    template <class Key> size_t Hash(const Key &key) const
    {
        uint64_t hash = 0;
        for (size_t w = 0; w < wordCount; ++w)
            hash = Mix(hash ^ key.getWords()[w]);
        return static_cast<size_t>(hash);
    }

    template <class Key> size_t Slot(const Key &key) const
    {
        size_t slot = Hash(key) & (buckets.size() - 1), firstDeleted = Deleted;
        for (;;)
        {
            const size_t entry = buckets[slot];
            if (!entry)
                return firstDeleted == Deleted ? slot : firstDeleted;
            if (entry == Deleted)
            {
                if (firstDeleted == Deleted)
                    firstDeleted = slot;
            }
            else
            {
                bool equal = true;
                for (size_t w = 0; w < wordCount; ++w)
                    if (keys[(entry - 1) * wordCount + w] != key.getWords()[w])
                    {
                        equal = false;
                        break;
                    }
                if (equal)
                    return slot;
            }
            slot = (slot + 1) & (buckets.size() - 1);
        }
    }

    template <class Key> size_t FindIndex(const Key &key) const
    {
        if (empty() || key.size() != nBits)
            return size();
        const size_t entry = buckets[Slot(key)];
        return entry && entry != Deleted ? entry - 1 : size();
    }

    void EnsureWidth(size_t width)
    {
        if (width > FastVectorBool::MaxWords * 64)
            throw std::length_error("Path integral states support at most 1024 qubits");
        if (wordCount && nBits != width)
        {
            if (!empty())
                throw std::invalid_argument("Amplitude keys must have the same qubit count");
            clear();
        }
        if (!wordCount)
        {
            nBits = width;
            wordCount = std::max<size_t>(1, (width + 63) / 64);
            if (values.capacity() > keys.max_size() / wordCount)
                throw std::length_error("Path integral key capacity exceeded");
            keys.reserve(values.capacity() * wordCount);
        }
    }

    static size_t TableSize(size_t count)
    {
        if (count > std::numeric_limits<size_t>::max() / 4)
            throw std::length_error("Path integral hash capacity exceeded");
        size_t result = 8;
        while (result / 2 < count)
            result *= 2;
        return result;
    }

    void Rehash(size_t tableSize)
    {
        buckets.resize(tableSize);
        Reindex();
    }

    void Reindex()
    {
        std::fill(buckets.begin(), buckets.end(), 0);
        for (size_t i = 0; i < size(); ++i)
        {
            size_t slot = Hash(State(i)) & (buckets.size() - 1);
            while (buckets[slot])
                slot = (slot + 1) & (buckets.size() - 1);
            buckets[slot] = i + 1;
        }
        occupied = size();
    }

    // Gate kernels only use nonthrowing callbacks. Injective transforms allow
    // key changes in place: rebuild the index after all entries have moved.
    template <class Function> void Transform(bool keysChange, Function &&function, bool shrinkIndex = false)
    {
        static_assert(noexcept(function(std::declval<Detail::MutableStateView &>(), std::declval<Complex &>())),
                      "An in-place amplitude transform must not throw");
        size_t output = 0;
        const size_t oldSize = size();
        for (size_t i = 0; i < oldSize; ++i)
        {
            Detail::MutableStateView state(keys.data() + i * wordCount, nBits);
            if (!function(state, values[i]))
                continue;
            if (i != output)
            {
                std::copy_n(keys.data() + i * wordCount, wordCount, keys.data() + output * wordCount);
                values[output] = values[i];
            }
            ++output;
        }
        values.resize(output);
        keys.resize(output * wordCount);
        if (keysChange || output != oldSize)
        {
            // Measurement knows the surviving population. Shrinking the
            // logical index before rebuilding avoids clearing the old peak;
            // CompactIfSparse separately releases excess physical capacity.
            if (shrinkIndex)
                buckets.resize(std::min(buckets.size(), TableSize(output)));
            Reindex();
            CompactIfSparse();
        }
    }

    std::vector<Complex> values;
    std::vector<uint64_t> keys;
    std::vector<size_t> buckets;
    size_t nBits = 0, wordCount = 0, occupied = 0;
};

} // namespace PathIntegral
} // namespace QC
