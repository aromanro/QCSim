#pragma once

#include "StructuredGate.h"

#include <Eigen/Eigen>

#include <algorithm>
#include <array>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace QC
{
namespace TensorNetworks
{

// The operator of a diagonal or multiplexed gate on many qubits, as a chain of site matrices.
//
// A table with one two-valued index per position is a product of matrices,
//   t[i_0 ... i_(m-1)] = G_0[i_0] G_1[i_1] ... G_(m-1)[i_(m-1)],
// one pair of matrices per position (a tensor train). The sizes of the matrices, the ranks of
// the table, say how much the gate couples the two sides of each cut: a product of one-qubit
// phases has rank one everywhere, a sum of nearest-neighbour couplings rank two, a table without
// structure up to 2^(m/2). Read along the chain, the pairs are the site tensors of the operator:
// a diagonal applies G[s] where the site holds s, and a multiplexed rotation, whose matrix is
// cos(theta/2) 1 - i sin(theta/2) P, takes the cosines and the sines as the two values of one
// more index, sitting at its target. Sites between the qubits of the gate pass the bond on.
//
// The ranks are found with singular value decompositions. A singular value is taken as zero
// when it does not exceed the larger dimension of its matrix times the machine epsilon times
// the largest one, the usual numerical rank: what is left out is of the size of the rounding
// of the table itself.
struct ChainOperator
{
    using Matrix = Eigen::MatrixXcd;
    // The matrices of one site by (output bit, input bit): operators[2 * out + in], each of
    // size (rank of the left bond) x (rank of the right bond). An empty one is zero.
    using Site = std::array<Matrix, 4>;

    Eigen::Index first = 0;
    std::vector<Site> sites;
    Eigen::Index maxRank = 1;
};

namespace TableTrainDetail
{
// The tensor train of `table`, whose index holds position 0 in its highest bit. The ranks may
// not exceed rankLimit.
inline std::vector<std::array<Eigen::MatrixXcd, 2>> DecomposeTable(std::vector<std::complex<double>> table,
                                                                   size_t positions, Eigen::Index rankLimit)
{
    using Matrix = Eigen::MatrixXcd;
    std::vector<std::array<Matrix, 2>> cores(positions);
    Eigen::Index rank = 1;
    // rest(row, col): what remains to be split, rows = (left bond, next bit), row major in the table.
    Eigen::Index cols = static_cast<Eigen::Index>(table.size() / 2);
    Matrix rest = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        table.data(), 2, cols);
    std::vector<std::complex<double>>().swap(table);
    for (size_t position = 0; position + 1 < positions; ++position)
    {
        const Eigen::JacobiSVD<Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(rest);
        const auto &values = svd.singularValues();
        const double floor = values.size() ? values[0] * double(std::max(rest.rows(), rest.cols())) *
                                                 std::numeric_limits<double>::epsilon()
                                           : 0.;
        Eigen::Index kept = 0;
        while (kept < values.size() && values[kept] > floor)
            ++kept;
        if (kept == 0)
            kept = 1; // a table of zeros
        if (kept > rankLimit)
            throw std::invalid_argument("The table of this gate couples the chain too strongly for a tensor network: "
                                        "its rank along the chain exceeds the limit");
        for (int bit = 0; bit < 2; ++bit)
        {
            cores[position][bit].resize(rank, kept);
            for (Eigen::Index a = 0; a < rank; ++a)
                cores[position][bit].row(a) = svd.matrixU().row(2 * a + bit).head(kept);
        }
        // (singular values) (V^dagger): kept x cols, then the next bit moves from the columns to the rows.
        const Matrix carried = values.head(kept).asDiagonal() * svd.matrixV().leftCols(kept).adjoint();
        const Eigen::Index half = cols / 2;
        Matrix next(2 * kept, half);
        for (Eigen::Index b = 0; b < kept; ++b)
            for (int bit = 0; bit < 2; ++bit)
                next.row(2 * b + bit) = carried.row(b).segment(bit * half, half);
        rest = std::move(next);
        rank = kept;
        cols = half;
    }
    for (int bit = 0; bit < 2; ++bit)
    {
        cores[positions - 1][bit].resize(rank, 1);
        for (Eigen::Index a = 0; a < rank; ++a)
            cores[positions - 1][bit](a, 0) = rest(2 * a + bit, 0);
    }
    return cores;
}
} // namespace TableTrainDetail

// The chain operator of a diagonal or of a multiplexed rotation. siteOf(qubit) is the position
// of a qubit of the gate in the chain; the qubits must sit at distinct sites.
template <class SiteOf>
ChainOperator MakeChainOperator(const Gates::StructuredGate &gate, const SiteOf &siteOf, Eigen::Index rankLimit)
{
    using Complex = std::complex<double>;
    using Kind = Gates::StructuredGate::Kind;
    if (gate.kind == Kind::Controlled)
        throw std::invalid_argument("A controlled gate has its own chain operator");
    const bool multiplexed = gate.kind == Kind::MultiplexedRotation;
    // A Z multiplexer is the diagonal with exp(-+ i theta / 2) on the two values of its target.
    const bool channel = multiplexed && gate.axis != Gates::RotationAxis::Z;
    const size_t tableBits = gate.qubits.size() + (multiplexed ? 1 : 0);
    if (tableBits < 2 || tableBits >= std::numeric_limits<size_t>::digits)
        throw std::invalid_argument("A chain operator needs a table on two or more qubits");
    // (site, table bit); the last table bit of a multiplexer is its target.
    std::vector<std::pair<Eigen::Index, unsigned>> positions;
    for (size_t k = 0; k < gate.qubits.size(); ++k)
        positions.emplace_back(static_cast<Eigen::Index>(siteOf(gate.qubits[k])), static_cast<unsigned>(k));
    if (multiplexed)
        positions.emplace_back(static_cast<Eigen::Index>(siteOf(gate.target)),
                               static_cast<unsigned>(gate.qubits.size()));
    std::sort(positions.begin(), positions.end());
    for (size_t p = 1; p < positions.size(); ++p)
        if (positions[p].first == positions[p - 1].first)
            throw std::invalid_argument("Gate qubits must sit at distinct sites");

    const size_t entries = size_t{1} << tableBits, selectors = size_t{1} << gate.qubits.size();
    std::vector<std::array<double, 2>> rotation;
    if (multiplexed)
    {
        rotation.resize(selectors);
        for (size_t s = 0; s < selectors; ++s)
            QC::detail::RotationHalfAngleCosSin(gate.angles[s], rotation[s][0], rotation[s][1]);
    }
    // The table in the order of the chain: position 0 in the highest bit of the index.
    std::vector<Complex> ordered(entries);
    for (size_t index = 0; index < entries; ++index)
    {
        size_t source = 0;
        for (size_t p = 0; p < tableBits; ++p)
            source |= ((index >> (tableBits - 1 - p)) & 1) << positions[p].second;
        if (!multiplexed)
            ordered[index] = gate.diagonal[source];
        else
        {
            const size_t selector = source & (selectors - 1);
            const bool second = (source >> gate.qubits.size()) & 1;
            if (channel)
                ordered[index] = rotation[selector][second ? 1 : 0];
            else
                ordered[index] =
                    Complex(rotation[selector][0], second ? rotation[selector][1] : -rotation[selector][1]);
        }
    }
    const auto cores = TableTrainDetail::DecomposeTable(std::move(ordered), tableBits, rankLimit);

    ChainOperator result;
    result.first = positions.front().first;
    const Eigen::Index last = positions.back().first;
    result.sites.resize(static_cast<size_t>(last - result.first + 1));
    size_t next = 0;
    Eigen::Index bond = 1;
    for (Eigen::Index site = result.first; site <= last; ++site)
    {
        ChainOperator::Site &operators = result.sites[static_cast<size_t>(site - result.first)];
        if (positions[next].first != site)
        {
            // not a qubit of the gate: the bond passes through
            operators[0] = operators[3] = ChainOperator::Matrix::Identity(bond, bond);
            continue;
        }
        const auto &core = cores[next];
        if (channel && positions[next].second == gate.qubits.size())
        {
            // cos 1 - i sin P on the target: -iX = [[0, -i], [-i, 0]], -iY = [[0, -1], [1, 0]]
            operators[0] = operators[3] = core[0];
            if (gate.axis == Gates::RotationAxis::X)
                operators[1] = operators[2] = Complex(0., -1.) * core[1];
            else
            {
                operators[1] = -core[1];
                operators[2] = core[1];
            }
        }
        else
        {
            operators[0] = core[0];
            operators[3] = core[1];
        }
        bond = core[0].cols();
        result.maxRank = std::max(result.maxRank, bond);
        ++next;
    }
    return result;
}

} // namespace TensorNetworks
} // namespace QC
