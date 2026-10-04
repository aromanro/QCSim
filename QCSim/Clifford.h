#pragma once

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
