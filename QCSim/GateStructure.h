#pragma once

#include <cstddef>
#include <complex>
#include <cstdint>

namespace QC {
namespace Gates {

// Describes the active block after removing controls on any local bits. It contains no
// coefficients: callers may supply the conjugated matrix (density-matrix rows).
struct GateStructure
{
    enum class Kind : uint8_t { Dense, Diagonal, Antidiagonal, Swap, ISwap, ISwapDag, Permutation };
    Kind kind = Kind::Dense;
    uint8_t controls = 0;
    uint8_t controlMask = 0;
    uint8_t controlValue = 0;
    uint32_t permutation = 0; // Three bits per active output row: its source column.

    unsigned Source(unsigned row) const { return (permutation >> (3 * row)) & 7u; }
};

template<class MatrixClass>
GateStructure ClassifyGateMatrix(const MatrixClass& matrix)
{
    using Scalar = typename MatrixClass::Scalar;
    using Kind = GateStructure::Kind;
    const size_t dimension = static_cast<size_t>(matrix.rows());
    GateStructure result;
    if (matrix.rows() != matrix.cols() || (dimension != 2 && dimension != 4 && dimension != 8)) return result;
    if (dimension == 2)
    {
        if (matrix(0, 1) == Scalar(0.) && matrix(1, 0) == Scalar(0.)) result.kind = Kind::Diagonal;
        else if (matrix(0, 0) == Scalar(0.) && matrix(1, 1) == Scalar(0.)) result.kind = Kind::Antidiagonal;
        return result;
    }

    // Every deviation from identity must be in the active sector. A bit is a
    // control exactly when all such row AND column indices share its value.
    // Exact comparisons are essential: small cross terms must never be lost.
    unsigned common = dimension > 2 ? unsigned(dimension - 1) : 0;
    unsigned anchor = unsigned(dimension - 1);
    bool found = false;
    for (unsigned r = 0; r < dimension && common; ++r)
        for (unsigned c = 0; c < dimension && common; ++c)
            if (matrix(r, c) != Scalar(r == c ? 1. : 0.))
            {
                if (!found) { anchor = r; found = true; }
                common &= ~(r ^ anchor) & ~(c ^ anchor);
            }
    // Keep at least one target, including identity and single-entry phases.
    if (common == dimension - 1) common &= ~1u;
    result.controlMask = uint8_t(common);
    result.controlValue = uint8_t(anchor & common);
    for (unsigned mask = common; mask; mask &= mask - 1) ++result.controls;
    size_t indices[8]{}, active = 0;
    for (size_t i = 0; i < dimension; ++i)
        if ((i & common) == result.controlValue) indices[active++] = i;
    if (active == 2)
    {
        if (matrix(indices[0], indices[1]) == Scalar(0.) && matrix(indices[1], indices[0]) == Scalar(0.)) result.kind = Kind::Diagonal;
        else if (matrix(indices[0], indices[0]) == Scalar(0.) && matrix(indices[1], indices[1]) == Scalar(0.)) result.kind = Kind::Antidiagonal;
        return result;
    }

    bool diagonal = true, antidiagonal = true, permutation = true;
    unsigned usedColumns = 0;
    for (size_t r = 0; r < active; ++r)
    {
        size_t source = active;
        for (size_t c = 0; c < active; ++c)
            if (matrix(indices[r], indices[c]) != Scalar(0.))
            {
                if (r != c) diagonal = false;
                if (r + c != active - 1) antidiagonal = false;
                if (source != active) permutation = false;
                source = c;
            }
        if (source == active || (usedColumns & (1u << source))) permutation = false;
        else
        {
            usedColumns |= 1u << source;
            result.permutation |= uint32_t(source) << (3 * r);
        }
    }
    if (diagonal) result.kind = Kind::Diagonal;
    else if (antidiagonal) result.kind = Kind::Antidiagonal;
    else if (permutation)
    {
        result.kind = Kind::Permutation;
        if (active == 4 && result.Source(0) == 0 && result.Source(1) == 2 && result.Source(2) == 1 && result.Source(3) == 3
            && matrix(indices[0], indices[0]) == Scalar(1.) && matrix(indices[3], indices[3]) == Scalar(1.))
        {
            const Scalar phase = matrix(indices[1], indices[2]);
            if (phase == matrix(indices[2], indices[1]) && (phase == Scalar(1.) || phase == Scalar(0., 1.) || phase == Scalar(0., -1.)))
                result.kind = phase == Scalar(1.) ? Kind::Swap : (phase == Scalar(0., 1.) ? Kind::ISwap : Kind::ISwapDag);
        }
    }
    return result;
}

}
}
