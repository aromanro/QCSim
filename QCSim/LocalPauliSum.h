#pragma once

#include <cassert>
#include <vector>

#include "CliffordTableau.h"
#include "LocalPauli.h"

namespace QC
{
namespace detail
{

// A gate on at most three local qubits, compiled as C S: the Pauli sum S
// acts first and the local Clifford C, recorded as gates, acts last.
// Building appends operators in time order. Each Clifford is moved past
// every later operator by conjugation, and rotations at multiples of pi/2
// are Cliffords, so S holds only the non-Clifford part. C is tracked by a
// local inverse tableau, always as wide as the largest gate, so one object
// can be reset and reused without allocating. Global phase is not tracked.
class LocalGate
{
  public:
    using Complex = LocalPauliSum::Complex;

    struct CliffordOp
    {
        enum class Kind : uint8_t
        {
            H,
            S,
            Sdg,
            X,
            Y,
            Z,
            CX,
            CZ,
            QuarterTurn,
            InverseQuarterTurn
        };
        Kind kind;
        uint8_t first = 0;  // the qubit, or the target of CX/CZ
        uint8_t second = 0; // the control of CX/CZ
        char axis = 'Z';    // the Pauli of a quarter turn
    };

    explicit LocalGate(size_t qubits)
        : inverseX(LocalPauli::MaxQubits), inverseZ(LocalPauli::MaxQubits), scratch(1, LocalPauli::MaxQubits)
    {
        Reset(qubits);
    }

    // Start an empty gate on the given number of local qubits.
    void Reset(size_t qubits)
    {
        assert(qubits >= 1 && qubits <= LocalPauli::MaxQubits);
        nrQubits = qubits;
        Clifford::detail::InverseMap(inverseX, inverseZ).SetIdentity();
        sum = LocalPauliSum::Identity();
        cliffords.clear();
    }

    LocalGate &H(size_t qubit)
    {
        return Append({CliffordOp::Kind::H, uint8_t(qubit)});
    }
    LocalGate &S(size_t qubit)
    {
        return Append({CliffordOp::Kind::S, uint8_t(qubit)});
    }
    LocalGate &Sdg(size_t qubit)
    {
        return Append({CliffordOp::Kind::Sdg, uint8_t(qubit)});
    }
    LocalGate &CX(size_t target, size_t control)
    {
        return Append({CliffordOp::Kind::CX, uint8_t(target), uint8_t(control)});
    }
    LocalGate &CZ(size_t target, size_t control)
    {
        return Append({CliffordOp::Kind::CZ, uint8_t(target), uint8_t(control)});
    }

    // exp(-i angle P / 2) for the Pauli P ('X', 'Y' or 'Z') on a local qubit.
    LocalGate &Rotate(size_t qubit, char axis, double angle)
    {
        long long quarterTurns = 0;
        if (TryGetQuarterTurns(angle, quarterTurns))
        {
            const uint8_t q = uint8_t(qubit);
            switch (((quarterTurns % 4) + 4) % 4)
            {
            case 0:
                return *this;
            case 1:
                return Append({CliffordOp::Kind::QuarterTurn, q, 0, axis});
            case 2:
                return Append({PauliKind(axis), q}); // -iP
            default:
                return Append({CliffordOp::Kind::InverseQuarterTurn, q, 0, axis});
            }
        }
        const double half = 0.5 * angle;
        LocalPauliSum rotation = LocalPauliSum::Identity(std::cos(half));
        AddConjugated(rotation, LocalPauli::Single(qubit, axis), Complex(0.0, -std::sin(half)));
        sum = rotation * sum;
        return *this;
    }

    // An exact operator on the local qubits, acting after everything so far.
    LocalGate &Multiply(const LocalPauliSum &op)
    {
        LocalPauliSum conjugated;
        op.ForEachTerm([&](LocalPauli pauli, Complex coefficient) { AddConjugated(conjugated, pauli, coefficient); });
        sum = conjugated * sum;
        return *this;
    }

    size_t GetNrQubits() const noexcept
    {
        return nrQubits;
    }
    const LocalPauliSum &Sum() const noexcept
    {
        return sum;
    }
    const std::vector<CliffordOp> &Cliffords() const noexcept
    {
        return cliffords;
    }

    // Apply the Clifford part, in order, to any map with the frame map's gate
    // names; qubits maps local qubits to the map's qubits.
    template <class Map> void ReplayCliffords(Map &map, const size_t *qubits) const
    {
        for (const auto &op : cliffords)
            Apply(map, op, qubits);
    }

  private:
    // The local tableau with the frame map's single quarter-turn signature.
    struct LocalMap
    {
        Clifford::detail::InverseMap map;
        Clifford::detail::TableauRow<false> scratch;
        void ApplyH(size_t q) noexcept
        {
            map.ApplyH(q);
        }
        void ApplyS(size_t q) noexcept
        {
            map.ApplyS(q);
        }
        void ApplySdg(size_t q) noexcept
        {
            map.ApplySdg(q);
        }
        void ApplyX(size_t q) noexcept
        {
            map.ApplyX(q);
        }
        void ApplyY(size_t q) noexcept
        {
            map.ApplyY(q);
        }
        void ApplyZ(size_t q) noexcept
        {
            map.ApplyZ(q);
        }
        void ApplyCX(size_t t, size_t c) noexcept
        {
            map.ApplyCX(t, c);
        }
        void ApplyCZ(size_t t, size_t c) noexcept
        {
            map.ApplyCZ(t, c);
        }
        void ApplyQuarterTurn(size_t q, char axis, bool inverse) noexcept
        {
            map.ApplyQuarterTurn(q, axis, inverse, scratch);
        }
    };

    static CliffordOp::Kind PauliKind(char axis) noexcept
    {
        return axis == 'X' ? CliffordOp::Kind::X : axis == 'Y' ? CliffordOp::Kind::Y : CliffordOp::Kind::Z;
    }

    template <class Map> static void Apply(Map &map, const CliffordOp &op, const size_t *qubits)
    {
        const size_t first = qubits[op.first], second = qubits[op.second];
        switch (op.kind)
        {
        case CliffordOp::Kind::H:
            map.ApplyH(first);
            break;
        case CliffordOp::Kind::S:
            map.ApplyS(first);
            break;
        case CliffordOp::Kind::Sdg:
            map.ApplySdg(first);
            break;
        case CliffordOp::Kind::X:
            map.ApplyX(first);
            break;
        case CliffordOp::Kind::Y:
            map.ApplyY(first);
            break;
        case CliffordOp::Kind::Z:
            map.ApplyZ(first);
            break;
        case CliffordOp::Kind::CX:
            map.ApplyCX(first, second);
            break;
        case CliffordOp::Kind::CZ:
            map.ApplyCZ(first, second);
            break;
        case CliffordOp::Kind::QuarterTurn:
            map.ApplyQuarterTurn(first, op.axis, false);
            break;
        case CliffordOp::Kind::InverseQuarterTurn:
            map.ApplyQuarterTurn(first, op.axis, true);
            break;
        }
    }

    LocalGate &Append(const CliffordOp &op)
    {
        cliffords.push_back(op);
        static constexpr size_t identity[LocalPauli::MaxQubits] = {0, 1, 2};
        LocalMap local{{inverseX, inverseZ}, scratch[0]};
        Apply(local, op, identity);
        return *this;
    }

    // Add coefficient C^dagger P C, for the Clifford C applied so far.
    void AddConjugated(LocalPauliSum &target, LocalPauli pauli, Complex coefficient)
    {
        auto image = scratch[0];
        image.Clear();
        for (size_t qubit = 0; qubit < nrQubits; ++qubit)
        {
            const char factor = pauli.Factor(qubit);
            if (factor != 'I')
                Clifford::detail::MultiplyImage(image, inverseX, inverseZ, qubit, factor);
        }
        const LocalPauli conjugated{uint8_t(image.X.words[0]), uint8_t(image.Z.words[0])};
        target.Add(conjugated, image.PhaseSign ? -coefficient : coefficient);
    }

    size_t nrQubits = 0;
    Clifford::detail::PackedTableau inverseX;
    Clifford::detail::PackedTableau inverseZ;
    Clifford::detail::PackedTableau scratch;
    LocalPauliSum sum;
    std::vector<CliffordOp> cliffords;
};

} // namespace detail
} // namespace QC
