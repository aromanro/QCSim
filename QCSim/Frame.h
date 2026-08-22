#pragma once

#include <algorithm>
#include <complex>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include "PauliStringXZ.h"

namespace QC {

	enum class NormalizationGateType : unsigned char {
		H = 0,
		S,
		CX,
		CZ
	};

	struct NormalizationGate {
		NormalizationGateType type = NormalizationGateType::H;
		size_t qubit1 = 0;
		size_t qubit2 = 0;
	};

	// A bounded, phase-aware representation of the inverse Clifford action
	// U^dagger P U.  Only the inverse X/Z images are needed by simulation; the
	// redundant forward tableau is deliberately not stored.  All Pauli rows
	// share contiguous packed buffers, so copying/restoring has constant
	// allocation count and Clifford updates avoid vector<bool> proxy traffic.
	class CliffordBasisMap {
		using Word = uint64_t;

		enum class Table : size_t {
			InverseX = 0,
			InverseZ,
			Count
		};

	public:
		// A non-owning packed Pauli view.  It stays valid until the owning map is
		// assigned or destroyed; Clifford mutations update the viewed contents,
		// including its sign bit.
		class PackedPauliView {
		public:
			PackedPauliView() = default;

			size_t GetNrQubits() const noexcept { return nrQubits; }
			size_t GetNrWords() const noexcept { return nrWords; }
			const Word* GetXWords() const noexcept { return xWords; }
			const Word* GetZWords() const noexcept { return zWords; }
			bool GetPhaseSign() const noexcept
			{
				return signWord != nullptr && (*signWord & signMask) != 0;
			}

			bool X(size_t qubit) const noexcept
			{
				return (xWords[qubit / BitsPerWord] & (Word(1) << (qubit % BitsPerWord))) != 0;
			}

			bool Z(size_t qubit) const noexcept
			{
				return (zWords[qubit / BitsPerWord] & (Word(1) << (qubit % BitsPerWord))) != 0;
			}

			bool HasX() const noexcept
			{
				for (size_t word = 0; word < nrWords; ++word)
					if (xWords[word] != 0)
						return true;
				return false;
			}

		private:
			friend class CliffordBasisMap;
			static constexpr size_t BitsPerWord = 64;

			PackedPauliView(const Word* x, const Word* z, size_t words,
				size_t qubits, const Word* packedSign, Word packedSignMask) noexcept
				: xWords(x), zWords(z), nrWords(words), nrQubits(qubits),
				signWord(packedSign), signMask(packedSignMask)
			{
			}

			const Word* xWords = nullptr;
			const Word* zWords = nullptr;
			size_t nrWords = 0;
			size_t nrQubits = 0;
			const Word* signWord = nullptr;
			Word signMask = 0;
		};

		CliffordBasisMap() = delete;

		explicit CliffordBasisMap(size_t qubits)
			: nrQubits(qubits),
			wordsPerRow(WordsFor(qubits)),
			signWords(WordsFor(qubits)),
			bits(static_cast<size_t>(Table::Count) * 2 * qubits * wordsPerRow, 0),
			signs(static_cast<size_t>(Table::Count) * signWords, 0),
			scratch(4 * wordsPerRow, 0)
		{
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
			{
				SetRowBit(Table::InverseX, qubit, false, qubit, true);
				SetRowBit(Table::InverseZ, qubit, true, qubit, true);
			}
		}

		size_t GetNrQubits() const noexcept
		{
			return nrQubits;
		}

		CliffordBasisMap(const CliffordBasisMap& other)
			: nrQubits(other.nrQubits),
			wordsPerRow(other.wordsPerRow),
			signWords(other.signWords),
			bits(other.bits),
			signs(other.signs),
			scratch(4 * other.wordsPerRow, 0)
		{
		}

		CliffordBasisMap& operator=(const CliffordBasisMap& other)
		{
			if (this == &other)
				return *this;

			nrQubits = other.nrQubits;
			wordsPerRow = other.wordsPerRow;
			signWords = other.signWords;
			// vector assignment reuses capacity for the common SaveState/RestoreState
			// case, while only these two buffers carry persistent map state.
			bits = other.bits;
			signs = other.signs;
			scratch.resize(4 * wordsPerRow);
			return *this;
		}

		CliffordBasisMap(CliffordBasisMap&& other) noexcept
			: nrQubits(std::exchange(other.nrQubits, 0)),
			wordsPerRow(std::exchange(other.wordsPerRow, 0)),
			signWords(std::exchange(other.signWords, 0)),
			bits(std::move(other.bits)),
			signs(std::move(other.signs)),
			scratch(std::move(other.scratch))
		{
			other.bits.clear();
			other.signs.clear();
			other.scratch.clear();
		}

		CliffordBasisMap& operator=(CliffordBasisMap&& other) noexcept
		{
			if (this == &other)
				return *this;

			nrQubits = std::exchange(other.nrQubits, 0);
			wordsPerRow = std::exchange(other.wordsPerRow, 0);
			signWords = std::exchange(other.signWords, 0);
			bits = std::move(other.bits);
			signs = std::move(other.signs);
			scratch = std::move(other.scratch);
			other.bits.clear();
			other.signs.clear();
			other.scratch.clear();
			return *this;
		}

		void ApplyH(size_t qubit)
		{
			SwapRows(Table::InverseX, qubit, Table::InverseZ, qubit);
		}

		void ApplyS(size_t qubit)
		{
			// S^dagger X S = -Y = -i XZ.
			HermitianProductInto(Table::InverseX, qubit,
				Table::InverseX, qubit, Table::InverseZ, qubit, 3);
		}

		void ApplyX(size_t qubit)
		{
			ToggleSign(Table::InverseZ, qubit);
		}

		void ApplyY(size_t qubit)
		{
			ToggleSign(Table::InverseX, qubit);
			ToggleSign(Table::InverseZ, qubit);
		}

		void ApplyZ(size_t qubit)
		{
			ToggleSign(Table::InverseX, qubit);
		}

		void ApplyCX(size_t target, size_t control)
		{
			// CX is self-inverse.  Its inverse conjugation maps
			// X_control -> X_control X_target and
			// Z_target  -> Z_control Z_target.
			HermitianProductInto(Table::InverseX, control,
				Table::InverseX, control, Table::InverseX, target);
			HermitianProductInto(Table::InverseZ, target,
				Table::InverseZ, control, Table::InverseZ, target);
		}

		PauliStringXZWithSign TransformToBasis(const PauliStringXZ& physicalPauli) const
		{
			if (physicalPauli.GetNrQubits() != nrQubits)
				throw std::invalid_argument("Pauli and Clifford basis dimensions do not match");
			if (nrQubits == 0)
				return PauliStringXZWithSign(0);

			auto& workspace = ThreadWorkspace(wordsPerRow);
			Word* physicalX = workspace.data();
			Word* physicalZ = physicalX + wordsPerRow;
			Word* logicalX = physicalZ + wordsPerRow;
			Word* logicalZ = logicalX + wordsPerRow;
			PackPauli(physicalPauli, physicalX, physicalZ);
			bool resultSign = false;
			TransformPacked(physicalX, physicalZ, false,
				logicalX, logicalZ, resultSign);
			return UnpackPauli(logicalX, logicalZ, resultSign);
		}

		PauliStringXZWithSign TransformToBasis(
			const PauliStringXZWithSign& physicalPauli) const
		{
			if (physicalPauli.GetNrQubits() != nrQubits)
				throw std::invalid_argument("Pauli and Clifford basis dimensions do not match");
			if (nrQubits == 0)
			{
				PauliStringXZWithSign result(0);
				result.PhaseSign = physicalPauli.PhaseSign;
				return result;
			}

			auto& workspace = ThreadWorkspace(wordsPerRow);
			Word* physicalX = workspace.data();
			Word* physicalZ = physicalX + wordsPerRow;
			Word* logicalX = physicalZ + wordsPerRow;
			Word* logicalZ = logicalX + wordsPerRow;
			PackPauli(physicalPauli, physicalX, physicalZ);
			bool resultSign = false;
			TransformPacked(physicalX, physicalZ, physicalPauli.PhaseSign,
				logicalX, logicalZ, resultSign);
			return UnpackPauli(logicalX, logicalZ, resultSign);
		}

		// Allocation-free views for the common single-qubit observables.  In
		// particular U^dagger Z_q U is already the q-th inverse-Z row.
		PackedPauliView TransformXToBasisPacked(size_t physicalQubit) const
		{
			ValidateQubit(physicalQubit);
			return View(Table::InverseX, physicalQubit);
		}

		PackedPauliView TransformZToBasisPacked(size_t physicalQubit) const
		{
			ValidateQubit(physicalQubit);
			return View(Table::InverseZ, physicalQubit);
		}

		PauliStringXZWithSign TransformXToBasis(size_t physicalQubit) const
		{
			ValidateQubit(physicalQubit);
			return UnpackPauli(View(Table::InverseX, physicalQubit));
		}

		PauliStringXZWithSign TransformZToBasis(size_t physicalQubit) const
		{
			ValidateQubit(physicalQubit);
			return UnpackPauli(View(Table::InverseZ, physicalQubit));
		}

		// Rebase the inverse Clifford map after a nondeterministic measurement
		// of a positive physical Pauli. logicalSigns are eigenvalue bits for the
		// conceptual signed forward-Z generators and are updated to label the
		// selected ray in the new basis; those forward rows need not be stored.
		void MeasurePauli(const PauliStringXZ& measuredPauli, bool outcome,
			std::vector<bool>& logicalSigns)
		{
			if (measuredPauli.GetNrQubits() != nrQubits
				|| logicalSigns.size() != nrQubits)
				throw std::invalid_argument(
					"Measurement and Clifford basis dimensions do not match");

			PackPauli(measuredPauli, Scratch(0), Scratch(1));
			bool logicalSign = false;
			TransformPacked(Scratch(0), Scratch(1), false,
				Scratch(2), Scratch(3), logicalSign);
			MeasurePacked(Scratch(2), Scratch(3), logicalSign,
				outcome, logicalSigns);
		}

		// Allocation-free specialization used by computational-basis sampling.
		void MeasureZ(size_t physicalQubit, bool outcome,
			std::vector<bool>& logicalSigns)
		{
			ValidateQubit(physicalQubit);
			if (logicalSigns.size() != nrQubits)
				throw std::invalid_argument(
					"Measurement and Clifford basis dimensions do not match");

			CopyWords(XRow(Table::InverseZ, physicalQubit), Scratch(0));
			CopyWords(ZRow(Table::InverseZ, physicalQubit), Scratch(1));
			const bool logicalSign = GetSign(Table::InverseZ, physicalQubit);
			MeasurePacked(Scratch(0), Scratch(1), logicalSign,
				outcome, logicalSigns);
		}

		// The inverse X/Z images form a signed symplectic basis exactly when
		// they have the canonical Pauli commutation relations.  This validates
		// the complete inverse-only representation without keeping a redundant
		// forward tableau solely as a consistency oracle.
		bool IsConsistent() const
		{
			for (size_t left = 0; left < nrQubits; ++left)
			{
				for (size_t right = left; right < nrQubits; ++right)
					if (RowsAnticommute(Table::InverseX, left,
						Table::InverseX, right)
						|| RowsAnticommute(Table::InverseZ, left,
							Table::InverseZ, right))
						return false;

				for (size_t right = 0; right < nrQubits; ++right)
					if (RowsAnticommute(Table::InverseX, left,
						Table::InverseZ, right) != (left == right))
						return false;
			}
			return true;
		}

	private:
		static constexpr size_t BitsPerWord = 64;

		static size_t WordsFor(size_t bits) noexcept
		{
			return (bits + BitsPerWord - 1) / BitsPerWord;
		}

		static std::vector<Word>& ThreadWorkspace(size_t words)
		{
			thread_local std::vector<Word> workspace;
			workspace.resize(4 * words);
			return workspace;
		}

		static int Mod4(int value) noexcept
		{
			const int result = value % 4;
			return result < 0 ? result + 4 : result;
		}

		static unsigned PopCount(Word value) noexcept
		{
#ifdef _MSC_VER
			return static_cast<unsigned>(__popcnt64(value));
#else
			return static_cast<unsigned>(__builtin_popcountll(value));
#endif
		}

		static int CountAndMod4(const Word* left, const Word* right,
			size_t words) noexcept
		{
			unsigned count = 0;
			for (size_t word = 0; word < words; ++word)
				count += PopCount(left[word] & right[word]);
			return static_cast<int>(count & 3U);
		}

		static int CountXorAndMod4(const Word* leftX, const Word* rightX,
			const Word* leftZ, const Word* rightZ, size_t words) noexcept
		{
			unsigned count = 0;
			for (size_t word = 0; word < words; ++word)
				count += PopCount((leftX[word] ^ rightX[word])
					& (leftZ[word] ^ rightZ[word]));
			return static_cast<int>(count & 3U);
		}

		static bool ParityAnd(const Word* left, const Word* right,
			size_t words) noexcept
		{
			unsigned parity = 0;
			for (size_t word = 0; word < words; ++word)
				parity ^= PopCount(left[word] & right[word]) & 1U;
			return parity != 0;
		}

		static bool GetBit(const Word* words, size_t bit) noexcept
		{
			return (words[bit / BitsPerWord]
				& (Word(1) << (bit % BitsPerWord))) != 0;
		}

		static void SetBit(Word* words, size_t bit, bool value) noexcept
		{
			const Word mask = Word(1) << (bit % BitsPerWord);
			Word& word = words[bit / BitsPerWord];
			if (value)
				word |= mask;
			else
				word &= ~mask;
		}

		void ValidateQubit(size_t qubit) const
		{
			if (qubit >= nrQubits)
				throw std::out_of_range("Qubit index is outside the Clifford basis");
		}

		size_t PlaneOffset(Table table, bool zPlane) const noexcept
		{
			return (static_cast<size_t>(table) * 2 + (zPlane ? 1 : 0))
				* nrQubits * wordsPerRow;
		}

		Word* XRow(Table table, size_t row) noexcept
		{
			return bits.data() + PlaneOffset(table, false) + row * wordsPerRow;
		}

		const Word* XRow(Table table, size_t row) const noexcept
		{
			return bits.data() + PlaneOffset(table, false) + row * wordsPerRow;
		}

		Word* ZRow(Table table, size_t row) noexcept
		{
			return bits.data() + PlaneOffset(table, true) + row * wordsPerRow;
		}

		const Word* ZRow(Table table, size_t row) const noexcept
		{
			return bits.data() + PlaneOffset(table, true) + row * wordsPerRow;
		}

		Word* Scratch(size_t plane) noexcept
		{
			return scratch.data() + plane * wordsPerRow;
		}

		void SetRowBit(Table table, size_t row, bool zPlane,
			size_t qubit, bool value) noexcept
		{
			SetBit(zPlane ? ZRow(table, row) : XRow(table, row), qubit, value);
		}

		bool GetSign(Table table, size_t row) const noexcept
		{
			const size_t offset = static_cast<size_t>(table) * signWords;
			return GetBit(signs.data() + offset, row);
		}

		void SetSign(Table table, size_t row, bool value) noexcept
		{
			const size_t offset = static_cast<size_t>(table) * signWords;
			SetBit(signs.data() + offset, row, value);
		}

		void ToggleSign(Table table, size_t row) noexcept
		{
			SetSign(table, row, !GetSign(table, row));
		}

		void ClearWords(Word* destination) const noexcept
		{
			std::fill_n(destination, wordsPerRow, Word(0));
		}

		void CopyWords(const Word* source, Word* destination) const noexcept
		{
			std::copy_n(source, wordsPerRow, destination);
		}

		void PackPauli(const PauliStringXZ& pauli, Word* x, Word* z) const
		{
			ClearWords(x);
			ClearWords(z);
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
			{
				if (pauli.X[qubit]) SetBit(x, qubit, true);
				if (pauli.Z[qubit]) SetBit(z, qubit, true);
			}
		}

		PauliStringXZWithSign UnpackPauli(const Word* x, const Word* z,
			bool sign) const
		{
			PauliStringXZWithSign result(nrQubits);
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
			{
				result.X[qubit] = GetBit(x, qubit);
				result.Z[qubit] = GetBit(z, qubit);
			}
			result.PhaseSign = sign;
			return result;
		}

		PauliStringXZWithSign UnpackPauli(const PackedPauliView& view) const
		{
			return UnpackPauli(view.GetXWords(), view.GetZWords(),
				view.GetPhaseSign());
		}

		PackedPauliView View(Table table, size_t row) const noexcept
		{
			const size_t signOffset = static_cast<size_t>(table) * signWords;
			return PackedPauliView(XRow(table, row), ZRow(table, row),
				wordsPerRow, nrQubits,
				signs.data() + signOffset + row / BitsPerWord,
				Word(1) << (row % BitsPerWord));
		}

		void MultiplyRightRaw(Word* leftX, Word* leftZ, int& phase,
			Table rightTable, size_t rightRow) const noexcept
		{
			const Word* rightX = XRow(rightTable, rightRow);
			const Word* rightZ = ZRow(rightTable, rightRow);
			phase = Mod4(phase
				+ CountAndMod4(rightX, rightZ, wordsPerRow)
				+ (GetSign(rightTable, rightRow) ? 2 : 0)
				+ (ParityAnd(leftZ, rightX, wordsPerRow) ? 2 : 0));
			for (size_t word = 0; word < wordsPerRow; ++word)
			{
				leftX[word] ^= rightX[word];
				leftZ[word] ^= rightZ[word];
			}
		}

		void TransformPacked(const Word* physicalX, const Word* physicalZ,
			bool physicalSign, Word* logicalX, Word* logicalZ,
			bool& logicalSign) const
		{
			ClearWords(logicalX);
			ClearWords(logicalZ);
			int phase = CountAndMod4(physicalX, physicalZ, wordsPerRow)
				+ (physicalSign ? 2 : 0);

			// P(x,z) = i^|x&z| X^x Z^z, so all X images precede all Z images.
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
				if (GetBit(physicalX, qubit))
					MultiplyRightRaw(logicalX, logicalZ, phase,
						Table::InverseX, qubit);
			for (size_t qubit = 0; qubit < nrQubits; ++qubit)
				if (GetBit(physicalZ, qubit))
					MultiplyRightRaw(logicalX, logicalZ, phase,
						Table::InverseZ, qubit);

			const int signPhase = Mod4(phase
				- CountAndMod4(logicalX, logicalZ, wordsPerRow));
			if (signPhase != 0 && signPhase != 2)
				throw std::logic_error(
					"A Clifford map produced a non-Hermitian Pauli");
			logicalSign = signPhase == 2;
		}

		void HermitianProductInto(Table destinationTable, size_t destinationRow,
			Table leftTable, size_t leftRow, Table rightTable, size_t rightRow,
			int extraIExponent = 0)
		{
			const Word* leftX = XRow(leftTable, leftRow);
			const Word* leftZ = ZRow(leftTable, leftRow);
			const Word* rightX = XRow(rightTable, rightRow);
			const Word* rightZ = ZRow(rightTable, rightRow);

			const int rawPhase = CountAndMod4(leftX, leftZ, wordsPerRow)
				+ (GetSign(leftTable, leftRow) ? 2 : 0)
				+ extraIExponent
				+ CountAndMod4(rightX, rightZ, wordsPerRow)
				+ (GetSign(rightTable, rightRow) ? 2 : 0)
				+ (ParityAnd(leftZ, rightX, wordsPerRow) ? 2 : 0);
			const int resultY = CountXorAndMod4(leftX, rightX,
				leftZ, rightZ, wordsPerRow);
			const int signPhase = Mod4(rawPhase - resultY);
			if (signPhase != 0 && signPhase != 2)
				throw std::logic_error(
					"A Clifford map produced a non-Hermitian Pauli product");

			Word* destinationX = XRow(destinationTable, destinationRow);
			Word* destinationZ = ZRow(destinationTable, destinationRow);
			for (size_t word = 0; word < wordsPerRow; ++word)
			{
				destinationX[word] = leftX[word] ^ rightX[word];
				destinationZ[word] = leftZ[word] ^ rightZ[word];
			}
			SetSign(destinationTable, destinationRow, signPhase == 2);
		}

		void SwapRows(Table leftTable, size_t leftRow,
			Table rightTable, size_t rightRow) noexcept
		{
			Word* leftX = XRow(leftTable, leftRow);
			Word* leftZ = ZRow(leftTable, leftRow);
			Word* rightX = XRow(rightTable, rightRow);
			Word* rightZ = ZRow(rightTable, rightRow);
			for (size_t word = 0; word < wordsPerRow; ++word)
			{
				std::swap(leftX[word], rightX[word]);
				std::swap(leftZ[word], rightZ[word]);
			}
			const bool sign = GetSign(leftTable, leftRow);
			SetSign(leftTable, leftRow, GetSign(rightTable, rightRow));
			SetSign(rightTable, rightRow, sign);
		}

		void ChangeInverseBasisInPlace(Table table, size_t row,
			const Word* measuredX, const Word* measuredZ, bool measuredSign,
			size_t pivot, int measuredY)
		{
			Word* x = XRow(table, row);
			Word* z = ZRow(table, row);
			const bool alphaPivot =
				ParityAnd(x, measuredZ, wordsPerRow)
				!= ParityAnd(z, measuredX, wordsPerRow);
			const bool oldAlphaPivot = GetBit(x, pivot);
			const int oldExponent = CountAndMod4(x, z, wordsPerRow)
				+ (GetSign(table, row) ? 2 : 0);

			if (oldAlphaPivot)
				for (size_t word = 0; word < wordsPerRow; ++word)
				{
					x[word] ^= measuredX[word];
					z[word] ^= measuredZ[word];
				}
			SetBit(x, pivot, alphaPivot);
			SetBit(z, pivot, oldAlphaPivot);

			const bool newXPivot = alphaPivot;
			const bool measuredZPivot = GetBit(measuredZ, pivot);
			const bool t = newXPivot
				!= ParityAnd(x, measuredZ, wordsPerRow)
				!= (newXPivot && measuredZPivot);
			int newExponent = CountAndMod4(x, z, wordsPerRow);
			if (oldAlphaPivot)
				newExponent += measuredY + (measuredSign ? 2 : 0)
					+ (t ? 2 : 0);

			const int signPhase = Mod4(oldExponent - newExponent);
			if (signPhase != 0 && signPhase != 2)
				throw std::logic_error(
					"Measurement produced an invalid inverse Clifford phase");
			SetSign(table, row, signPhase == 2);
		}

		bool RowsAnticommute(Table leftTable, size_t leftRow,
			Table rightTable, size_t rightRow) const noexcept
		{
			return ParityAnd(XRow(leftTable, leftRow),
				ZRow(rightTable, rightRow), wordsPerRow)
				!= ParityAnd(ZRow(leftTable, leftRow),
					XRow(rightTable, rightRow), wordsPerRow);
		}

		void MeasurePacked(const Word* logicalX, const Word* logicalZ,
			bool logicalSign,
			bool outcome, std::vector<bool>& logicalSigns)
		{
			size_t pivot = nrQubits;
			for (size_t logical = 0; logical < nrQubits; ++logical)
				if (GetBit(logicalX, logical))
				{
					pivot = logical;
					break;
				}

			if (pivot == nrQubits)
				throw std::logic_error(
					"Cannot rebase a deterministic Pauli measurement");

			const int logicalY = CountAndMod4(logicalX, logicalZ, wordsPerRow);
			for (size_t row = 0; row < nrQubits; ++row)
				ChangeInverseBasisInPlace(Table::InverseX, row,
					logicalX, logicalZ, logicalSign, pivot, logicalY);
			for (size_t row = 0; row < nrQubits; ++row)
				ChangeInverseBasisInPlace(Table::InverseZ, row,
					logicalX, logicalZ, logicalSign, pivot, logicalY);

			const bool oldPivotSign = logicalSigns[pivot];
			for (size_t logical = 0; logical < nrQubits; ++logical)
			{
				if (logical == pivot) continue;

				if (GetBit(logicalX, logical))
					logicalSigns[logical] =
						logicalSigns[logical] != oldPivotSign;
			}

			logicalSigns[pivot] = outcome;
		}

		size_t nrQubits = 0;
		size_t wordsPerRow = 0;
		size_t signWords = 0;
		std::vector<Word> bits;
		std::vector<Word> signs;
		std::vector<Word> scratch;
	};

	// Performance-oriented frame used by ExtendedStabilizer.  Component signs
	// are logical computational-basis labels b in the moving Clifford basis U,
	// so a component represents amplitude * U|b>.  The signed Clifford map owns
	// the complete phase convention; no unpacked stabilizers, row maps, or phase
	// correction vectors are duplicated here.
	class ExtendedFrame {
	public:
		using Signs = std::vector<bool>;

		explicit ExtendedFrame(size_t nrQubits)
			: amplitudes(1, std::complex<double>(1.0, 0.0)),
			signs(1, Signs(nrQubits, false)),
			cliffordBasis(nrQubits)
		{
		}

		size_t GetNrQubits() const noexcept
		{
			return cliffordBasis.GetNrQubits();
		}

		size_t GetFrameSize() const noexcept
		{
			return amplitudes.size();
		}

		std::vector<std::complex<double>> amplitudes;
		std::vector<Signs> signs;
		CliffordBasisMap cliffordBasis;
	};

}
