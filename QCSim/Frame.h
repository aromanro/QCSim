#pragma once

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "CliffordTableau.h"

namespace QC
{

// Logical computational-basis labels stored component-major in one packed
// allocation. This avoids vector<bool>'s proxy access and one heap allocation
// per component while retaining a small, value-like public type.
class PackedComponentLabels
{
  public:
    using Word = uint64_t;
    static constexpr size_t BitsPerWord = 64;

    PackedComponentLabels() = default;

    PackedComponentLabels(size_t bitsPerLabel, size_t nrLabels)
        : nrBits(bitsPerLabel), wordsPerLabel(WordsFor(bitsPerLabel)), nrLabels(nrLabels),
          words(nrLabels * wordsPerLabel, Word(0))
    {
    }

    size_t size() const noexcept
    {
        return nrLabels;
    }

    size_t GetNrBits() const noexcept
    {
        return nrBits;
    }

    size_t GetNrWords() const noexcept
    {
        return wordsPerLabel;
    }

    const Word *LabelWords(size_t label) const noexcept
    {
        return wordsPerLabel == 0 ? nullptr : words.data() + label * wordsPerLabel;
    }

    Word *LabelWords(size_t label) noexcept
    {
        return wordsPerLabel == 0 ? nullptr : words.data() + label * wordsPerLabel;
    }

    bool Get(size_t label, size_t bit) const noexcept
    {
        return (LabelWords(label)[bit / BitsPerWord] & (Word(1) << (bit % BitsPerWord))) != 0;
    }

    void Set(size_t label, size_t bit, bool value) noexcept
    {
        Word &word = LabelWords(label)[bit / BitsPerWord];
        const Word mask = Word(1) << (bit % BitsPerWord);
        if (value)
            word |= mask;
        else
            word &= ~mask;
    }

    void reserve(size_t labels)
    {
        words.reserve(labels * wordsPerLabel);
    }

    void resize(size_t labels)
    {
        words.resize(labels * wordsPerLabel, Word(0));
        nrLabels = labels;
    }

    void clear() noexcept
    {
        words.clear();
        nrLabels = 0;
    }

    void Reset(size_t bitsPerLabel)
    {
        nrBits = bitsPerLabel;
        wordsPerLabel = WordsFor(bitsPerLabel);
        nrLabels = 0;
        words.clear();
    }

    void swap(PackedComponentLabels &other) noexcept
    {
        using std::swap;
        swap(nrBits, other.nrBits);
        swap(wordsPerLabel, other.wordsPerLabel);
        swap(nrLabels, other.nrLabels);
        words.swap(other.words);
    }

    void CopyLabel(size_t destination, size_t source) noexcept
    {
        if (destination == source || wordsPerLabel == 0)
            return;
        std::memmove(LabelWords(destination), LabelWords(source), wordsPerLabel * sizeof(Word));
    }

    void Append(const Word *labelWords)
    {
        if (wordsPerLabel != 0)
            words.insert(words.end(), labelWords, labelWords + wordsPerLabel);
        ++nrLabels;
    }

    void AppendXor(const Word *labelWords, const Word *xorMask)
    {
        const size_t oldWords = words.size();
        words.resize(oldWords + wordsPerLabel);
        for (size_t word = 0; word < wordsPerLabel; ++word)
            words[oldWords + word] = labelWords[word] ^ xorMask[word];
        ++nrLabels;
    }

    bool operator==(const PackedComponentLabels &other) const noexcept
    {
        return nrBits == other.nrBits && nrLabels == other.nrLabels && words == other.words;
    }

    bool operator!=(const PackedComponentLabels &other) const noexcept
    {
        return !(*this == other);
    }

  private:
    static size_t WordsFor(size_t bits) noexcept
    {
        return (bits + BitsPerWord - 1) / BitsPerWord;
    }

    size_t nrBits = 0;
    size_t wordsPerLabel = 0;
    size_t nrLabels = 0;
    std::vector<Word> words;
};

// A reusable open-addressed index over PackedComponentLabels. Slots contain
// component indices, so label-buffer reallocations do not invalidate it.
class PackedComponentIndex
{
  public:
    using Word = PackedComponentLabels::Word;
    static constexpr size_t NotFound = std::numeric_limits<size_t>::max();

    void Build(const PackedComponentLabels &labels, size_t expectedLabels = 0)
    {
        const size_t required = std::max(labels.size(), expectedLabels);
        const size_t slotsNeeded = SlotsFor(required);
        if (slots.size() < slotsNeeded)
            slots.resize(slotsNeeded, NotFound);
        else if (slots.size() / slotsNeeded > 4)
            // Retain the allocation but drop a stale peak logical size. This keeps
            // later rebuild clears proportional to a frame that collapsed sharply.
            slots.resize(slotsNeeded);
        std::fill(slots.begin(), slots.end(), NotFound);
        for (size_t component = 0; component < labels.size(); ++component)
            Insert(labels, component);
    }

    size_t FindXor(const PackedComponentLabels &labels, const Word *source, const Word *xorMask) const noexcept
    {
        if (slots.empty())
            return NotFound;
        const size_t nrWords = labels.GetNrWords();
        size_t slot = static_cast<size_t>(HashXor(source, xorMask, nrWords)) & (slots.size() - 1);
        for (;;)
        {
            const size_t component = slots[slot];
            if (component == NotFound)
                return NotFound;
            if (EqualXor(labels.LabelWords(component), source, xorMask, nrWords))
                return component;
            slot = (slot + 1) & (slots.size() - 1);
        }
    }

    void Insert(const PackedComponentLabels &labels, size_t component) noexcept
    {
        const Word *candidate = labels.LabelWords(component);
        size_t slot = static_cast<size_t>(Hash(candidate, labels.GetNrWords())) & (slots.size() - 1);
        while (slots[slot] != NotFound)
            slot = (slot + 1) & (slots.size() - 1);
        slots[slot] = component;
    }

    size_t Capacity() const noexcept
    {
        return slots.size() / 2;
    }

  private:
    static uint64_t Mix(uint64_t value) noexcept
    {
        value ^= value >> 30;
        value *= UINT64_C(0xbf58476d1ce4e5b9);
        value ^= value >> 27;
        value *= UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31);
    }

    static uint64_t Hash(const Word *words, size_t nrWords) noexcept
    {
        uint64_t hash = Mix(UINT64_C(0x9e3779b97f4a7c15) ^ nrWords);
        for (size_t word = 0; word < nrWords; ++word)
            hash = Mix(hash ^ Mix(words[word] + word));
        return hash;
    }

    static uint64_t HashXor(const Word *left, const Word *right, size_t nrWords) noexcept
    {
        uint64_t hash = Mix(UINT64_C(0x9e3779b97f4a7c15) ^ nrWords);
        for (size_t word = 0; word < nrWords; ++word)
            hash = Mix(hash ^ Mix((left[word] ^ right[word]) + word));
        return hash;
    }

    static bool EqualXor(const Word *candidate, const Word *source, const Word *xorMask, size_t nrWords) noexcept
    {
        for (size_t word = 0; word < nrWords; ++word)
            if (candidate[word] != (source[word] ^ xorMask[word]))
                return false;
        return true;
    }

    static size_t SlotsFor(size_t labels)
    {
        const size_t minimumLabels = std::max<size_t>(labels, 1);
        if (minimumLabels > std::numeric_limits<size_t>::max() / 2)
            throw std::length_error("Too many packed component labels");
        const size_t minimumSlots = 2 * minimumLabels;
        size_t slotsNeeded = 8;
        while (slotsNeeded < minimumSlots)
        {
            if (slotsNeeded > std::numeric_limits<size_t>::max() / 2)
                throw std::length_error("Packed component index is too large");
            slotsNeeded *= 2;
        }
        return slotsNeeded;
    }

    std::vector<size_t> slots;
};

// The inverse Clifford action U^dagger P U of a frame's basis. Rows, phase
// arithmetic, gate updates and measurement rebasing are shared with the
// Clifford stabilizer simulator (see Clifford::detail::InverseMap); only the
// redundant forward tableau is deliberately not stored. Both images live in
// contiguous packed tableaux, so copies have constant allocation count and
// assignments between equal shapes reuse the existing storage.
class CliffordBasisMap
{
    using Tableau = Clifford::detail::PackedTableau;
    using Word = Clifford::detail::Word;

  public:
    using Row = Clifford::detail::TableauRow<true>;

    CliffordBasisMap() = delete;

    explicit CliffordBasisMap(size_t qubits) : inverseX(qubits), inverseZ(qubits), work(WorkRows, qubits)
    {
        Map().SetIdentity();
    }

    // Scratch rows carry no state and are not copied.
    CliffordBasisMap(const CliffordBasisMap &other)
        : inverseX(other.inverseX), inverseZ(other.inverseZ), work(WorkRows, other.GetNrQubits())
    {
    }

    CliffordBasisMap &operator=(const CliffordBasisMap &other)
    {
        if (this == &other)
            return *this;
        // SaveState/RestoreState copy equal shapes without allocating.
        if (inverseX.HasSameShape(other.inverseX))
        {
            inverseX.CopyFrom(other.inverseX);
            inverseZ.CopyFrom(other.inverseZ);
        }
        else
        {
            CliffordBasisMap copy(other);
            swap(copy);
        }
        return *this;
    }

    CliffordBasisMap(CliffordBasisMap &&) noexcept = default;
    CliffordBasisMap &operator=(CliffordBasisMap &&) noexcept = default;

    void swap(CliffordBasisMap &other) noexcept
    {
        inverseX.swap(other.inverseX);
        inverseZ.swap(other.inverseZ);
        work.swap(other.work);
    }

    size_t GetNrQubits() const noexcept
    {
        return inverseZ.size();
    }

    // Callers validate qubit indices.
    void ApplyH(size_t qubit) noexcept
    {
        Map().ApplyH(qubit);
    }

    void ApplyS(size_t qubit) noexcept
    {
        Map().ApplyS(qubit);
    }

    void ApplySdg(size_t qubit) noexcept
    {
        Map().ApplySdg(qubit);
    }

    void ApplyX(size_t qubit) noexcept
    {
        Map().ApplyX(qubit);
    }

    void ApplyY(size_t qubit) noexcept
    {
        Map().ApplyY(qubit);
    }

    void ApplyZ(size_t qubit) noexcept
    {
        Map().ApplyZ(qubit);
    }

    void ApplySx(size_t qubit) noexcept
    {
        Map().ApplySx(qubit);
    }

    void ApplySxDag(size_t qubit) noexcept
    {
        Map().ApplySxDag(qubit);
    }

    void ApplyK(size_t qubit) noexcept
    {
        Map().ApplyK(qubit);
    }

    void ApplyCX(size_t target, size_t control) noexcept
    {
        Map().ApplyCX(target, control);
    }

    void ApplyCY(size_t target, size_t control) noexcept
    {
        Map().ApplyCY(target, control);
    }

    void ApplyCZ(size_t target, size_t control) noexcept
    {
        Map().ApplyCZ(target, control);
    }

    void ApplySwap(size_t qubit1, size_t qubit2) noexcept
    {
        Map().ApplySwap(qubit1, qubit2);
    }

    void ApplyISwap(size_t qubit1, size_t qubit2) noexcept
    {
        Map().ApplyISwap(qubit1, qubit2);
    }

    void ApplyISwapDag(size_t qubit1, size_t qubit2) noexcept
    {
        Map().ApplyISwapDag(qubit1, qubit2);
    }

    // Apply one signed pi/2 rotation exp(-+i pi/4 P) about the physical
    // one-qubit Pauli P ('X', 'Y' or 'Z') to the moving basis. The map tracks
    // conjugation only; callers retain the canonical rotation's global phase
    // in component amplitudes.
    void ApplyQuarterTurn(size_t physicalQubit, char axis, bool inverse = false) noexcept
    {
        Map().ApplyQuarterTurn(physicalQubit, axis, inverse, work[ImageRow]);
    }

    // U^dagger X_q U and U^dagger Z_q U are stored rows. The views stay valid
    // until the map is assigned or destroyed, and see later updates.
    Row ImageX(size_t physicalQubit) const noexcept
    {
        return inverseX[physicalQubit];
    }

    Row ImageZ(size_t physicalQubit) const noexcept
    {
        return inverseZ[physicalQubit];
    }

    // Build the image of a product of one-qubit physical Paulis in scratch
    // storage. Image() is overwritten by the next build.
    void ClearImage() const noexcept
    {
        work[ImageRow].Clear();
    }

    void MultiplyImage(size_t physicalQubit, char pauli) const noexcept
    {
        MultiplyImageInto(work[ImageRow], physicalQubit, pauli);
    }

    // The same product into caller storage of the map's width.
    void MultiplyImageInto(Clifford::detail::TableauRow<false> image, size_t physicalQubit, char pauli) const noexcept
    {
        Clifford::detail::MultiplyImage(image, inverseX, inverseZ, physicalQubit, pauli);
    }

    // All inverse Z images, U^dagger Z_q U for every physical qubit q.
    const Tableau &ZImages() const noexcept
    {
        return inverseZ;
    }

    // Accumulate the factors' rows by XOR alone, without phase arithmetic.
    // The X part decides whether the product is off-diagonal in the
    // logical basis. If no factor row has an X part, the rows are commuting
    // diagonal Paulis and the accumulated Z part and sign are exact.
    // Returns whether a factor row had an X part.
    bool XorImage(size_t physicalQubit, char pauli) const noexcept
    {
        auto image = work[ImageRow];
        bool hasX = false;
        const auto accumulate = [&image, &hasX](const Row &row) noexcept {
            Word x = 0;
            for (size_t word = 0; word < image.Words(); ++word)
            {
                x |= row.X.words[word];
                image.X.words[word] ^= row.X.words[word];
                image.Z.words[word] ^= row.Z.words[word];
            }
            image.PhaseSign ^= bool(row.PhaseSign);
            hasX |= x != 0;
        };
        if (pauli != 'Z')
            accumulate(inverseX[physicalQubit]);
        if (pauli != 'X')
            accumulate(inverseZ[physicalQubit]);
        return hasX;
    }

    Row Image() const noexcept
    {
        return static_cast<const Tableau &>(work)[ImageRow];
    }

    // Collapse after a random measurement of physical Z with the given
    // outcome. With Q = U^dagger Z U = phase X^x Z^z, pivot is the lowest bit
    // of x (Clifford::detail::FindPivot). label is the measured component's
    // logical basis label; afterwards it is label xor x (if it contained the
    // pivot) with the pivot set to the outcome.
    void MeasureZ(size_t physicalQubit, size_t pivot, bool outcome, Word *label, bool multithreading) noexcept
    {
        const auto measured = inverseZ[physicalQubit];
        const size_t pivotWord = pivot / BitsPerWord;
        const Word mask = Word(1) << (pivot % BitsPerWord);
        if (label[pivotWord] & mask)
            for (size_t word = 0; word < measured.Words(); ++word)
                label[word] ^= measured.X.words[word];
        label[pivotWord] = outcome ? (label[pivotWord] | mask) : (label[pivotWord] & ~mask);
        RebaseZ(physicalQubit, pivot, multithreading);
    }

    // Rebase only the map, with the same pivot. Component labels carry
    // outcomes, so the map changes basis with outcome zero. This is used
    // directly when post-measurement labels were built in the new basis.
    void RebaseZ(size_t physicalQubit, size_t pivot, bool multithreading) noexcept
    {
        assert(inverseZ[physicalQubit].X[pivot]);
        Map().CollapseZ(physicalQubit, pivot, false, work[MeasuredRow],
                        Clifford::detail::CollapseWorkers(GetNrQubits(), multithreading));
    }

    // The inverse X/Z images form a signed symplectic basis exactly when
    // they have the canonical Pauli commutation relations. This validates
    // the complete inverse-only representation without keeping a redundant
    // forward tableau solely as a consistency oracle.
    bool IsConsistent() const noexcept
    {
        return Clifford::detail::IsSymplectic(inverseX, inverseZ);
    }

  private:
    static constexpr size_t BitsPerWord = 64;

    enum : size_t
    {
        MeasuredRow = 0,
        ImageRow,
        WorkRows
    };

    Clifford::detail::InverseMap Map() noexcept
    {
        return {inverseX, inverseZ};
    }

    Tableau inverseX;
    Tableau inverseZ;
    // Measurement and transform scratch; instances are not thread-safe.
    mutable Tableau work;
};

// Performance-oriented frame used by ExtendedStabilizer.  Component signs
// are logical computational-basis labels b in the moving Clifford basis U,
// so a component represents amplitude * U|b>.  The signed Clifford map owns
// the complete phase convention; no unpacked stabilizers, row maps, or phase
// correction vectors are duplicated here.
class ExtendedFrame
{
  public:
    using Word = PackedComponentLabels::Word;

    explicit ExtendedFrame(size_t nrQubits)
        : amplitudes(1, std::complex<double>(1.0, 0.0)), signs(nrQubits, 1), cliffordBasis(nrQubits),
          nextSigns(nrQubits, 0)
    {
    }

    ExtendedFrame(const ExtendedFrame &other)
        : amplitudes(other.amplitudes), signs(other.signs), cliffordBasis(other.cliffordBasis),
          nextSigns(other.GetNrQubits(), 0)
    {
    }

    ExtendedFrame &operator=(const ExtendedFrame &other)
    {
        if (this == &other)
            return *this;
        amplitudes = other.amplitudes;
        signs = other.signs;
        cliffordBasis = other.cliffordBasis;
        componentIndexValid = false;
        nextAmplitudes.clear();
        nextSigns.Reset(other.GetNrQubits());
        measurementPairs.clear();
        componentOrderWorkspace.clear();
        componentMagnitudeWorkspace.clear();
        return *this;
    }

    ExtendedFrame(ExtendedFrame &&) noexcept = default;
    ExtendedFrame &operator=(ExtendedFrame &&) noexcept = default;

    size_t GetNrQubits() const noexcept
    {
        return cliffordBasis.GetNrQubits();
    }

    size_t GetFrameSize() const noexcept
    {
        return amplitudes.size();
    }

    void EnsureComponentIndex(size_t expectedComponents = 0) const
    {
        if (!componentIndexValid || componentIndex.Capacity() < expectedComponents)
        {
            componentIndex.Build(signs, expectedComponents);
            componentIndexValid = true;
        }
    }

    size_t FindXorComponent(const Word *label, const Word *xorMask) const noexcept
    {
        return componentIndex.FindXor(signs, label, xorMask);
    }

    void InvalidateComponentIndex() const noexcept
    {
        componentIndexValid = false;
    }

    std::vector<std::complex<double>> amplitudes;
    PackedComponentLabels signs;
    CliffordBasisMap cliffordBasis;

  private:
    friend class ExtendedStabilizer;

    mutable PackedComponentIndex componentIndex;
    mutable bool componentIndexValid = false;
    std::vector<std::complex<double>> nextAmplitudes;
    PackedComponentLabels nextSigns;
    std::vector<std::pair<size_t, size_t>> measurementPairs;
    std::vector<size_t> componentOrderWorkspace;
    std::vector<double> componentMagnitudeWorkspace;
};

} // namespace QC
