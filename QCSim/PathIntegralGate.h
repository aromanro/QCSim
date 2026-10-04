#pragma once

#include "QuantumGate.h"
#include <algorithm>
#include <array>
#include <complex>
#include <limits>
#include <stdexcept>

namespace QC
{
namespace PathIntegral
{
namespace Detail
{

// This is built from the current matrix and simulator threshold, not the
// legacy gate predicates. Each input column lists only retained transitions.
struct SparseGate
{
    std::array<size_t, 3> qubits;
    std::array<unsigned, 8> counts{};
    std::array<std::array<unsigned, 8>, 8> rows{};
    std::array<std::array<std::complex<double>, 8>, 8> values{};
    unsigned arity;
    unsigned maxFanout = 0;
    bool diagonal = true;
    bool injective = true;
    bool distinctQubits = true;

    SparseGate(const QC::Gates::AppliedGate<> &gate, double epsilon)
        : qubits{gate.getQubit1(), gate.getQubit2(), gate.getQubit3()},
          arity(static_cast<unsigned>(gate.getQubitsNumber()))
    {
        if (arity == 0 || arity > 3)
            throw std::invalid_argument("Path integral gates must act on one to three qubits");
        const auto &matrix = gate.getRawOperatorMatrix();
        const unsigned dimension = 1u << arity;
        if (matrix.rows() != dimension || matrix.cols() != dimension)
            throw std::invalid_argument("Invalid path integral gate dimensions");
        unsigned usedRows = 0;
        for (unsigned column = 0; column < dimension; ++column)
        {
            for (unsigned row = 0; row < dimension; ++row)
            {
                const auto value = matrix(row, column);
                if (!(std::norm(value) > epsilon))
                    continue;
                const unsigned index = counts[column]++;
                rows[column][index] = row;
                values[column][index] = value;
                diagonal = diagonal && row == column;
                if (index || (usedRows & (1u << row)))
                    injective = false;
                usedRows |= 1u << row;
            }
            maxFanout = std::max(maxFanout, counts[column]);
        }
        // Repeated target indices can destroy injectivity in the global state.
        for (unsigned i = 0; i < arity; ++i)
            for (unsigned j = i + 1; j < arity; ++j)
                if (qubits[i] == qubits[j])
                    injective = distinctQubits = false;
    }

    template <unsigned Arity, class State> unsigned Column(const State &state) const
    {
        unsigned column = 0;
        for (unsigned bit = 0; bit < Arity; ++bit)
            column |= unsigned(state.get(qubits[bit])) << bit;
        return column;
    }

    template <unsigned Arity, class State> void SetRow(State &state, unsigned row) const
    {
        for (unsigned bit = 0; bit < Arity; ++bit)
            state.set(qubits[bit], ((row >> bit) & 1u) != 0);
    }

    size_t ReserveSize(size_t inputSize, size_t width, size_t maximum, unsigned estimate = 0) const
    {
        // Dense transitions often reconverge. Start with a bounded estimate;
        // geometric growth still handles genuinely expanding support.
        const size_t fanout = estimate ? std::min(estimate, maxFanout) : std::min<unsigned>(maxFanout, 2);
        size_t result = fanout && inputSize <= maximum / fanout ? inputSize * fanout : maximum;
        if (width < std::numeric_limits<size_t>::digits)
            result = std::min(result, size_t{1} << width);
        return result;
    }
};

} // namespace Detail
} // namespace PathIntegral
} // namespace QC
