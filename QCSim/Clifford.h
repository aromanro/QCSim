#pragma once

#include "LocalCliffordRecognition.h"
#include "LocalPauli.h"
#include "StabilizerState.h"
#include <memory>
#include <string>

namespace QC
{
namespace Clifford
{

class StabilizerSimulator : public StabilizerState
{
  public:
    StabilizerSimulator() = default;

    explicit StabilizerSimulator(size_t n) : StabilizerState(n)
    {
    }

    // Gate updates are shared with the extended stabilizer's frame map, see
    // detail::InverseMap. Gates touch a constant number of packed rows, with no
    // parallel launch. Diagonal gates leave the computational-basis cache valid.
    void ApplyH(size_t q)
    {
        BeginGate(q);
        Map().ApplyH(q);
    }

    void ApplyS(size_t q)
    {
        ValidateQubit(q);
        Map().ApplyS(q);
    }

    void ApplySdg(size_t q)
    {
        ValidateQubit(q);
        Map().ApplySdg(q);
    }

    void ApplyX(size_t q)
    {
        ValidateQubit(q);
        Map().ApplyX(q);
        FlipDistributionBit(q);
    }

    void ApplyY(size_t q)
    {
        ValidateQubit(q);
        Map().ApplyY(q);
        FlipDistributionBit(q);
    }

    void ApplyZ(size_t q)
    {
        ValidateQubit(q);
        Map().ApplyZ(q);
    }

    void ApplySx(size_t q)
    {
        BeginGate(q);
        Map().ApplySx(q);
    }

    void ApplySxDag(size_t q)
    {
        BeginGate(q);
        Map().ApplySxDag(q);
    }

    void ApplyK(size_t q)
    {
        BeginGate(q);
        Map().ApplyK(q);
    }

    // Only exactly encoded Clifford angles are accepted; nearby rotations are
    // not rounded onto the tableau. Global phase is immaterial to this backend.
    void ApplyRx(size_t q, double angle)
    {
        ApplySingleRotation(q, 'X', angle);
    }

    void ApplyRy(size_t q, double angle)
    {
        ApplySingleRotation(q, 'Y', angle);
    }

    void ApplyRz(size_t q, double angle)
    {
        ApplySingleRotation(q, 'Z', angle);
    }

    void ApplyPhase(size_t q, double angle)
    {
        ApplySingleRotation(q, 'Z', angle);
    }

    // The generic bridge recognizes exact matrix contents of the supported
    // canonical Clifford gates. It deliberately rejects unrecognized matrices
    // rather than rounding a near-Clifford operator onto this backend. A
    // compact gate on at most three qubits is recognized through its dense
    // matrix, whose local bits are the given qubits in order (for an applied
    // compact gate: its target, then its controls; a diagonal's table bits).
    // Nothing changes before the gate is known to be a Clifford.
    template <class MatrixClass>
    void ApplyGate(const Gates::QuantumGateWithOp<MatrixClass> &source, size_t q0, size_t q1 = 0, size_t q2 = 0)
    {
        const size_t qubits[] = {q0, q1, q2};
        if (const Gates::StructuredGate *compact = source.getStructuredGate())
        {
            const Gates::UnitaryLocalMatrix matrix = compact->DenseMatrix<Gates::UnitaryLocalMatrix>(3);
            ApplyCliffordMatrix(matrix, compact->Arity(), qubits);
            return;
        }
        const size_t count = source.getQubitsNumber();
        if (count < 1 || count > 3)
            throw std::invalid_argument("Clifford gates require a matrix on one to three qubits");
        ApplyCliffordMatrix(source.getRawOperatorMatrix(), count, qubits);
    }

    template <class MatrixClass> void ApplyGate(const Gates::AppliedGate<MatrixClass> &gate)
    {
        ApplyGate(gate, gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
    }

    void ApplyCX(size_t target, size_t control)
    {
        BeginPair(target, control);
        Map().ApplyCX(target, control);
    }

    void ApplyCY(size_t target, size_t control)
    {
        BeginPair(target, control);
        Map().ApplyCY(target, control);
    }

    void ApplyCZ(size_t target, size_t control)
    {
        ValidatePair(target, control);
        Map().ApplyCZ(target, control);
    }

    void ApplySwap(size_t a, size_t b)
    {
        ValidatePair(a, b, true);
        if (a == b)
            return;
        InvalidateDistribution();
        Map().ApplySwap(a, b);
    }

    void ApplyISwap(size_t a, size_t b)
    {
        BeginPair(a, b);
        Map().ApplyISwap(a, b);
    }

    void ApplyISwapDag(size_t a, size_t b)
    {
        BeginPair(a, b);
        Map().ApplyISwapDag(a, b);
    }

    // Arguments are the low and high local gate bits, respectively.
    void ApplyDCX(size_t q0, size_t q1)
    {
        BeginPair(q0, q1);
        Map().ApplyCX(q1, q0);
        Map().ApplyCX(q0, q1);
    }

    void ApplyECR(size_t q0, size_t q1)
    {
        BeginPair(q0, q1);
        Map().ApplyS(q0);
        Map().ApplySx(q1);
        Map().ApplyCX(q1, q0);
        Map().ApplyX(q0);
    }

    // Pauli-product rotations are Clifford exactly at multiples of pi/2.
    // Other angles are rejected before changing the tableau or its caches.
    void ApplyRXX(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'X', 'X', angle);
    }

    void ApplyRYY(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'Y', 'Y', angle);
    }

    void ApplyRZZ(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'Z', 'Z', angle);
    }

    void ApplyRZX(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'Z', 'X', angle);
    }

    void ApplyRXY(size_t a, size_t b, double angle)
    {
        ApplyPairRotation(a, b, 'X', 'Y', angle);
    }

    void ApplyXXPlusYY(size_t a, size_t b, double theta, double beta = 0.)
    {
        ApplyXXYY(a, b, theta, beta, true);
    }

    void ApplyXXMinusYY(size_t a, size_t b, double theta, double beta = 0.)
    {
        ApplyXXYY(a, b, theta, beta, false);
    }

    double ExpectationValue(const std::string &pauliString) const
    {
        if (pauliString.size() > getNrQubits())
            throw std::invalid_argument("Pauli string exceeds the number of qubits");
        std::vector<std::pair<size_t, char>> positions;
        std::pair<size_t, char> first{0, 'I'};
        for (size_t q = 0; q < pauliString.size(); ++q)
        {
            char c;
            switch (pauliString[q])
            {
            case 'I':
            case 'i':
                continue;
            case 'X':
            case 'x':
                c = 'X';
                break;
            case 'Y':
            case 'y':
                c = 'Y';
                break;
            case 'Z':
            case 'z':
                c = 'Z';
                break;
            default:
                throw std::runtime_error("Invalid operator in the Pauli string");
            }
            if (first.second == 'I')
                first = {q, c};
            else
            {
                if (positions.empty())
                {
                    positions.reserve(pauliString.size());
                    positions.push_back(first);
                }
                positions.emplace_back(q, c);
            }
        }
        if (first.second == 'I')
            return 1.0;
        // Sparse one-qubit observables need no allocation or scratch product.
        if (positions.empty())
        {
            if (first.second != 'Y')
            {
                const auto row = first.second == 'X' ? inverseX[first.first] : inverseZ[first.first];
                return row.HasX() ? 0.0 : (row.PhaseSign ? -1.0 : 1.0);
            }
            const auto x = inverseX[first.first], z = inverseZ[first.first];
            unsigned phase = 1 + 2 * unsigned(bool(x.PhaseSign) != bool(z.PhaseSign));
            for (size_t w = 0; w < x.Words(); ++w)
            {
                if (x.X.words[w] != z.X.words[w])
                    return 0.0;
                phase += detail::Popcount(x.X.words[w] & x.Z.words[w]) + detail::Popcount(z.X.words[w] & z.Z.words[w]) +
                         2 * detail::Popcount(x.Z.words[w] & z.X.words[w]);
            }
            assert((phase & 1) == 0);
            return (phase & 2) ? -1.0 : 1.0;
        }
        detail::PackedTableau work(1, getNrQubits());
        auto result = work[0];
        bool allDiagonal = true, diagonalSign = false;
        const auto accumulate = [&](auto row) {
            diagonalSign ^= bool(row.PhaseSign);
            for (size_t w = 0; w < result.Words(); ++w)
            {
                const auto x = row.X.words[w];
                result.X.words[w] ^= x;
                if (x)
                    allDiagonal = false;
            }
        };
        for (const auto &op : positions)
        {
            if (op.second != 'Z')
                accumulate(inverseX[op.first]);
            if (op.second != 'X')
                accumulate(inverseZ[op.first]);
        }
        if (result.HasX())
            return 0.0;
        if (allDiagonal)
            return diagonalSign ? -1.0 : 1.0;
        result.Clear();
        for (const auto &op : positions)
            detail::MultiplyImage(result, inverseX, inverseZ, op.first, op.second);
        return result.PhaseSign ? -1.0 : 1.0;
    }

    std::unique_ptr<StabilizerSimulator> Clone() const
    {
        return std::make_unique<StabilizerSimulator>(*this);
    }

  private:
    template <class Matrix> void ApplyCliffordMatrix(const Matrix &matrix, size_t count, const size_t *qubits)
    {
        for (size_t i = 0; i < count; ++i)
        {
            ValidateQubit(qubits[i]);
            for (size_t j = 0; j < i; ++j)
                if (qubits[i] == qubits[j])
                    throw std::invalid_argument("Clifford gate requires distinct qubits");
        }
        if (static_cast<size_t>(matrix.rows()) != (size_t(1) << count) || matrix.cols() != matrix.rows() ||
            !matrix.allFinite())
            throw std::invalid_argument("Clifford gate matrix has invalid dimensions or coefficients");
        const auto *compiled = QC::detail::FindCanonicalClifford(matrix);
        if (!compiled)
            throw std::invalid_argument("Gate matrix is not a recognized canonical Clifford");
        using Kind = QC::detail::LocalGate::CliffordOp::Kind;
        for (const auto &op : compiled->Cliffords())
        {
            const size_t first = qubits[op.first], second = qubits[op.second];
            switch (op.kind)
            {
            case Kind::H:
                ApplyH(first);
                break;
            case Kind::S:
                ApplyS(first);
                break;
            case Kind::Sdg:
                ApplySdg(first);
                break;
            case Kind::X:
                ApplyX(first);
                break;
            case Kind::Y:
                ApplyY(first);
                break;
            case Kind::Z:
                ApplyZ(first);
                break;
            case Kind::CX:
                ApplyCX(first, second);
                break;
            case Kind::CZ:
                ApplyCZ(first, second);
                break;
            case Kind::QuarterTurn:
            case Kind::InverseQuarterTurn:
                ApplySingleRotation(first, op.axis,
                                    op.kind == Kind::QuarterTurn ? 1.57079632679489661923 : -1.57079632679489661923);
                break;
            }
        }
    }

    void ApplySingleRotation(size_t q, char axis, double angle)
    {
        ValidateQubit(q);
        long long turns = 0;
        if (!QC::detail::TryGetQuarterTurns(angle, turns))
            throw std::invalid_argument("Clifford rotation angle must be an exact multiple of pi/2");
        const unsigned count = unsigned((turns % 4 + 4) % 4);
        if (!count)
            return;
        if (axis == 'Z')
        {
            if (count == 1)
                Map().ApplyS(q);
            else if (count == 2)
                Map().ApplyZ(q);
            else
                Map().ApplySdg(q);
            return;
        }
        if (count == 2)
        {
            if (axis == 'X')
                ApplyX(q);
            else
                ApplyY(q);
            return;
        }
        InvalidateDistribution();
        if (axis == 'Y')
            Map().ApplySdg(q);
        if (count == 1)
            Map().ApplySx(q);
        else
            Map().ApplySxDag(q);
        if (axis == 'Y')
            Map().ApplyS(q);
    }

    void ApplyXXYY(size_t a, size_t b, double theta, double beta, bool plus)
    {
        ValidatePair(a, b);
        long long turns = 0, phaseTurns = 0;
        if (!std::isfinite(theta) || !std::isfinite(beta) || !QC::detail::TryGetQuarterTurns(.5 * theta, turns))
            throw std::invalid_argument("Clifford XX+/-YY angle must be an exact multiple of pi");
        const unsigned count = unsigned((turns % 4 + 4) % 4);
        if (count & 1U)
        {
            if (!QC::detail::TryGetQuarterTurns(beta, phaseTurns))
                throw std::invalid_argument("Clifford XX+/-YY phase must be an exact multiple of pi/2");
        }
        else
        {
            if (count == 2)
            {
                Map().ApplyZ(a);
                Map().ApplyZ(b);
            }
            return;
        }
        InvalidateDistribution();
        const size_t phaseQubit = plus ? a : b;
        unsigned phase = unsigned((phaseTurns % 4 + 4) % 4);
        if (!plus)
            phase = (4 - phase) & 3U;
        for (unsigned k = 0; k < phase; ++k)
            Map().ApplyS(phaseQubit);
        ApplyPairRotation(a, b, 'X', 'X', .5 * theta);
        ApplyPairRotation(a, b, 'Y', 'Y', (plus ? .5 : -.5) * theta);
        for (unsigned k = 0; k < phase; ++k)
            Map().ApplySdg(phaseQubit);
    }

    void ApplyPairRotation(size_t a, size_t b, char first, char second, double angle)
    {
        ValidatePair(a, b);
        long long turns = 0;
        if (!QC::detail::TryGetQuarterTurns(angle, turns))
            throw std::invalid_argument("Clifford rotation angle must be an exact multiple of pi/2");
        const unsigned count = unsigned((turns % 4 + 4) % 4);
        if (!count)
            return;
        InvalidateDistribution();
        const size_t qubits[] = {a, b};
        const char axes[] = {first, second};
        for (unsigned i = 0; i < 2; ++i)
        {
            if (axes[i] == 'Y')
                Map().ApplySdg(qubits[i]);
            if (axes[i] != 'Z')
                Map().ApplyH(qubits[i]);
        }
        Map().ApplyCX(b, a);
        for (unsigned i = 0; i < count; ++i)
            Map().ApplyS(b);
        Map().ApplyCX(b, a);
        for (unsigned i = 0; i < 2; ++i)
        {
            if (axes[i] != 'Z')
                Map().ApplyH(qubits[i]);
            if (axes[i] == 'Y')
                Map().ApplyS(qubits[i]);
        }
    }

    void BeginGate(size_t q)
    {
        ValidateQubit(q);
        InvalidateDistribution();
    }

    void BeginPair(size_t a, size_t b)
    {
        ValidatePair(a, b);
        InvalidateDistribution();
    }
};

} // namespace Clifford
} // namespace QC
