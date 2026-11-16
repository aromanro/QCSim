#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

#include <Eigen/Dense>

namespace QC
{
namespace TensorNetworks
{
namespace ThreeSite
{
using Matrix = Eigen::MatrixXcd;
using Index = Eigen::Index;

// Multiply by 2^exponent without first forming an overflowing reciprocal for
// subnormal inputs. Ordinary exponents retain Eigen's vectorized multiplication.
template <class Derived> void ScalePowerOfTwo(Eigen::MatrixBase<Derived> &matrix, int exponent)
{
    if (exponent >= -1022 && exponent <= 1023)
    {
        matrix *= std::ldexp(1., exponent);
        return;
    }
    for (Index column = 0; column < matrix.cols(); ++column)
        for (Index row = 0; row < matrix.rows(); ++row)
            if constexpr (Eigen::NumTraits<typename Derived::Scalar>::IsComplex)
            {
                const auto value = matrix(row, column);
                matrix(row, column) = {std::ldexp(value.real(), exponent), std::ldexp(value.imag(), exponent)};
            }
            else
                matrix(row, column) = std::ldexp(matrix(row, column), exponent);
}

inline void Validate(const Matrix &gate, Index count, Index target, Index control1, Index control2,
                     bool adjacent = false)
{
    if (gate.rows() != 8 || gate.cols() != 8 || !gate.allFinite())
        throw std::invalid_argument("Three-qubit gate must be a finite 8x8 matrix");
    std::array<Index, 3> sites{target, control1, control2};
    std::sort(sites.begin(), sites.end());
    if (sites[0] < 0 || sites[2] >= count)
        throw std::invalid_argument("Qubit index out of bounds");
    if (sites[0] == sites[1] || sites[1] == sites[2])
        throw std::invalid_argument("Three-qubit gates require distinct qubits");
    if (adjacent && sites[2] - sites[0] != 2)
        throw std::invalid_argument("Three-qubit gates need to act on consecutive sites");
}

inline bool IsUnitary(const Matrix &gate)
{
    // A unitary coefficient cannot exceed one. Check before forming U^dagger U
    // so extreme, finite filters cannot overflow during classification.
    return gate.real().cwiseAbs().maxCoeff() <= 1. + 1E-12 && gate.imag().cwiseAbs().maxCoeff() <= 1. + 1E-12 &&
           (gate.adjoint() * gate).isIdentity(1E-12);
}

// Tensor memory order is (left, physical, right), including ket+2*bra for MPO.
// Return rows (left, physical1, physical2), columns (physical3, right).
// Choose the contraction tree from the actual four bond dimensions.
template <int D, class Tensor> Matrix Contract(const Tensor &a, const Tensor &b, const Tensor &c)
{
    constexpr int Last = D == 2 ? 2 : 3;
    const Index L = a.dimension(0), M = a.dimension(Last), N = b.dimension(Last), R = c.dimension(Last);
    const double leftCost = double(D * D) * L * N * (M + double(D) * R);
    const double rightCost = double(D * D) * M * R * (N + double(D) * L);
    const Eigen::Map<const Matrix> am(a.data(), D * L, M), cm(c.data(), N, D * R);
    Matrix result(D * D * L, D * R);
    if (leftCost < rightCost || (leftCost == rightCost && L * N <= M * R))
    {
        const Eigen::Map<const Matrix> bm(b.data(), M, D * N);
        const Matrix pair = am * bm;
        const Eigen::Map<const Matrix> pairView(pair.data(), D * D * L, N);
        result.noalias() = pairView * cm;
    }
    else
    {
        const Eigen::Map<const Matrix> bm(b.data(), D * M, N);
        const Matrix pair = bm * cm;
        const Eigen::Map<const Matrix> pairView(pair.data(), M, D * D * R);
        Eigen::Map<Matrix> resultView(result.data(), D * L, D * D * R);
        resultView.noalias() = am * pairView;
    }
    return result;
}

inline Eigen::Matrix<std::complex<double>, 8, 8> OrderedGate(const Matrix &gate, Index first, Index target,
                                                             Index control1, Index control2)
{
    Eigen::Matrix<std::complex<double>, 8, 8> ordered;
    std::array<int, 8> permutation{};
    for (int i = 0; i < 8; ++i)
        permutation[i] = ((i >> (target - first)) & 1) | (((i >> (control1 - first)) & 1) << 1) |
                         (((i >> (control2 - first)) & 1) << 2);
    for (int col = 0; col < 8; ++col)
        for (int row = 0; row < 8; ++row)
            ordered(row, col) = gate(permutation[row], permutation[col]);
    return ordered;
}

// Used only before rescaling an extreme operator. Coefficients in unoccupied
// input columns cannot contribute, but their magnitude must not determine the
// scale of the occupied columns (for example 1e-310 beside an unused 1e300).
// Support is exact: no small amplitude or matrix coefficient is thresholded.
template <int D>
void RestrictToInputSupport(Matrix &gate, const Matrix &theta, Index first, Index target, Index control1,
                            Index control2, Index L, Index R)
{
    std::array<bool, 8> occupied{};
    for (int physical = 0; physical < D * D * D; ++physical)
    {
        const int p1 = physical % D, p2 = (physical / D) % D, p3 = physical / (D * D);
        const int ket = (p1 & 1) + 2 * (p2 & 1) + 4 * (p3 & 1);
        const int bra = D == 2 ? ket : (p1 >> 1) + 2 * (p2 >> 1) + 4 * (p3 >> 1);
        if (occupied[ket] && occupied[bra])
            continue;
        bool found = false;
        for (Index r = 0; r < R && !found; ++r)
            for (Index l = 0; l < L; ++l)
                if (theta(l + L * (p1 + D * p2), p3 + D * r) != std::complex<double>(0.))
                {
                    found = true;
                    break;
                }
        if (found)
            occupied[ket] = occupied[bra] = true;
    }
    for (int physical = 0; physical < 8; ++physical)
        if (!occupied[physical])
        {
            const int column = ((physical >> (target - first)) & 1) | (((physical >> (control1 - first)) & 1) << 1) |
                               (((physical >> (control2 - first)) & 1) << 2);
            gate.col(column).setZero();
        }
}

inline int OperatorPhysicalIndex(int ket, int bra)
{
    return (ket & 1) + 2 * (bra & 1) + 4 * ((ket >> 1) & 1) + 8 * ((bra >> 1) & 1) + 16 * ((ket >> 2) & 1) +
           32 * ((bra >> 2) & 1);
}

// Each column here is a contiguous length-L physical block. Stage changed
// outputs before overwriting any input; untouched controlled sectors need no
// copy. Reuse scratch over R, keeping small bonds on the stack kernel below.
template <int D>
void ApplySlices(Matrix &theta, const Eigen::Matrix<std::complex<double>, 8, 8> &ordered,
                 const std::array<std::array<int, 8>, 8> &sources, const std::array<int, 8> &counts, bool onePerRow,
                 bool sparse, Index L, Index R)
{
    std::array<bool, 8> identity{};
    for (int row = 0; row < 8; ++row)
        identity[row] = counts[row] == 1 && sources[row][0] == row && ordered(row, row) == std::complex<double>(1.);
    if (std::all_of(identity.begin(), identity.end(), [](bool value) { return value; }))
        return;
    Matrix output(L, D * D * D), ket;
    if constexpr (D == 4)
        if (!onePerRow)
            ket.resize(L, 64);
    for (Index r = 0; r < R; ++r)
    {
        Eigen::Map<Matrix> input(theta.data() + D * D * D * L * r, L, D * D * D);
        if constexpr (D == 2)
        {
            if (sparse)
            {
                for (int row = 0; row < 8; ++row)
                    if (!identity[row])
                    {
                        output.col(row).setZero();
                        for (int entry = 0; entry < counts[row]; ++entry)
                        {
                            const int source = sources[row][entry];
                            output.col(row) += ordered(row, source) * input.col(source);
                        }
                    }
            }
            else
                output.noalias() = input * ordered.transpose();
            for (int row = 0; row < 8; ++row)
                if (!identity[row])
                    input.col(row) = output.col(row);
        }
        else
        {
            if (onePerRow)
            {
                for (int bra = 0; bra < 8; ++bra)
                    for (int row = 0; row < 8; ++row)
                        if (!(identity[row] && identity[bra]))
                            output.col(OperatorPhysicalIndex(row, bra)) =
                                ordered(row, sources[row][0]) *
                                input.col(OperatorPhysicalIndex(sources[row][0], sources[bra][0])) *
                                std::conj(ordered(bra, sources[bra][0]));
            }
            else
            {
                for (int bra = 0; bra < 8; ++bra)
                    for (int row = 0; row < 8; ++row)
                        if (!identity[row])
                        {
                            auto column = ket.col(OperatorPhysicalIndex(row, bra));
                            column.setZero();
                            for (int entry = 0; entry < counts[row]; ++entry)
                            {
                                const int source = sources[row][entry];
                                column += ordered(row, source) * input.col(OperatorPhysicalIndex(source, bra));
                            }
                        }
                for (int bra = 0; bra < 8; ++bra)
                    for (int row = 0; row < 8; ++row)
                        if (!(identity[row] && identity[bra]))
                        {
                            auto column = output.col(OperatorPhysicalIndex(row, bra));
                            column.setZero();
                            for (int entry = 0; entry < counts[bra]; ++entry)
                            {
                                const int source = sources[bra][entry];
                                if (identity[row])
                                    column +=
                                        std::conj(ordered(bra, source)) * input.col(OperatorPhysicalIndex(row, source));
                                else
                                    column +=
                                        std::conj(ordered(bra, source)) * ket.col(OperatorPhysicalIndex(row, source));
                            }
                        }
            }
            for (int bra = 0; bra < 8; ++bra)
                for (int row = 0; row < 8; ++row)
                    if (!(identity[row] && identity[bra]))
                        input.col(OperatorPhysicalIndex(row, bra)) = output.col(OperatorPhysicalIndex(row, bra));
        }
    }
}

template <int D>
void Apply(Matrix &theta, const Matrix &gate, Index first, Index target, Index control1, Index control2, Index L,
           Index R)
{
    const auto ordered = OrderedGate(gate, first, target, control1, control2);
    std::array<std::array<int, 8>, 8> sources{};
    std::array<int, 8> counts{};
    int nonzeros = 0;
    for (int row = 0; row < 8; ++row)
        for (int col = 0; col < 8; ++col)
            if (ordered(row, col) != std::complex<double>(0.))
            {
                sources[row][counts[row]++] = col;
                ++nonzeros;
            }
    const bool onePerRow = std::all_of(counts.begin(), counts.end(), [](int count) { return count == 1; });
    const bool sparse = nonzeros <= 24;
    // Measured crossover includes R=1: allocations lose on tiny blocks. Dense
    // MPOs keep the efficient fixed 8x8 ket/bra products at every bond size.
    if ((D == 2 && ((sparse && L >= 4) || (!sparse && L >= 8))) || (D == 4 && sparse && L >= 8))
    {
        ApplySlices<D>(theta, ordered, sources, counts, onePerRow, sparse, L, R);
        return;
    }
    for (Index r = 0; r < R; ++r)
        for (Index l = 0; l < L; ++l)
        {
            if constexpr (D == 2)
            {
                Eigen::Matrix<std::complex<double>, 8, 1> input;
                for (int p = 0; p < 8; ++p)
                    input[p] = theta(l + L * (p & 3), (p >> 2) + 2 * r);
                Eigen::Matrix<std::complex<double>, 8, 1> output;
                if (sparse)
                    for (int row = 0; row < 8; ++row)
                    {
                        output[row] = 0.;
                        for (int entry = 0; entry < counts[row]; ++entry)
                        {
                            const int col = sources[row][entry];
                            output[row] += ordered(row, col) * input[col];
                        }
                    }
                else
                    output.noalias() = ordered * input;
                for (int p = 0; p < 8; ++p)
                    theta(l + L * (p & 3), (p >> 2) + 2 * r) = output[p];
            }
            else
            {
                Eigen::Matrix<std::complex<double>, 8, 8> input;
                for (int bra = 0; bra < 8; ++bra)
                    for (int ket = 0; ket < 8; ++ket)
                    {
                        const int p1 = (ket & 1) + 2 * (bra & 1);
                        const int p2 = ((ket >> 1) & 1) + 2 * ((bra >> 1) & 1);
                        const int p3 = (ket >> 2) + 2 * (bra >> 2);
                        input(ket, bra) = theta(l + L * (p1 + 4 * p2), p3 + 4 * r);
                    }
                // Diagonal/permutation gates need one multiply per leg, and
                // controlled gates use only their exact nonzero entries.
                Eigen::Matrix<std::complex<double>, 8, 8> output;
                if (onePerRow)
                    for (int bra = 0; bra < 8; ++bra)
                        for (int ket = 0; ket < 8; ++ket)
                            output(ket, bra) = ordered(ket, sources[ket][0]) * input(sources[ket][0], sources[bra][0]) *
                                               std::conj(ordered(bra, sources[bra][0]));
                else if (sparse)
                {
                    Eigen::Matrix<std::complex<double>, 8, 8> ket;
                    ket.setZero();
                    output.setZero();
                    for (int row = 0; row < 8; ++row)
                        for (int entry = 0; entry < counts[row]; ++entry)
                        {
                            const int source = sources[row][entry];
                            ket.row(row) += ordered(row, source) * input.row(source);
                        }
                    for (int col = 0; col < 8; ++col)
                        for (int entry = 0; entry < counts[col]; ++entry)
                        {
                            const int source = sources[col][entry];
                            output.col(col) += std::conj(ordered(col, source)) * ket.col(source);
                        }
                }
                else
                {
                    // Two 8x8 contractions; never build a dense 64x64 superoperator.
                    const Eigen::Matrix<std::complex<double>, 8, 8> ket = ordered * input;
                    output.noalias() = ket * ordered.adjoint();
                }
                for (int bra = 0; bra < 8; ++bra)
                    for (int ketIndex = 0; ketIndex < 8; ++ketIndex)
                    {
                        const int p1 = (ketIndex & 1) + 2 * (bra & 1);
                        const int p2 = ((ketIndex >> 1) & 1) + 2 * ((bra >> 1) & 1);
                        const int p3 = (ketIndex >> 2) + 2 * (bra >> 2);
                        theta(l + L * (p1 + 4 * p2), p3 + 4 * r) = output(ketIndex, bra);
                    }
            }
        }
}

struct Split
{
    Matrix left, right;
    Eigen::VectorXd lambda;
    bool truncated = false;
};

inline Matrix Weighted(const Matrix &theta, const Eigen::VectorXd *lambda, Index L)
{
    Matrix result = theta;
    if (lambda)
        for (Index col = 0; col < result.cols(); ++col)
            for (Index offset = 0; offset < result.rows(); offset += L)
                result.col(col).segment(offset, L).array() *= lambda->array();
    return result;
}
} // namespace ThreeSite
} // namespace TensorNetworks
} // namespace QC
