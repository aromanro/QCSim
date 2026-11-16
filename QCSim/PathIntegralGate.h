#pragma once

#include "QuantumGate.h"
#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace QC
{
namespace PathIntegral
{
namespace Detail
{

inline bool RetainGateCoefficient(std::complex<double> value, double epsilon)
{
    // Exact mode retains representable amplitudes even when squaring would
    // underflow. Positive thresholds retain the existing probability rule.
    // A retained coefficient is decided by the first comparison alone, as
    // with a plain threshold; only a dropped one asks which mode applies.
    return std::norm(value) > epsilon || (epsilon == 0. && value != std::complex<double>(0.));
}

using CompactGate = Gates::AppliedGateDetail::StructuredPayload;

// How a compact operation moves basis states under the given threshold: the
// largest number of outputs of one input, whether no input changes, and
// whether no two inputs share an output. With `transposed` the answer is for
// the transposed operation, a backward step. A multiplexed rotation is decided
// by its first angles when they branch, and otherwise by the extremes derived
// once from its table, never by another scan of the table.
inline void ClassifyCompact(const CompactGate &gate, double epsilon, bool transposed, unsigned &maxFanout,
                            bool &diagonal, bool &injective)
{
    using Kind = Gates::StructuredGate::Kind;
    if (gate.OperationKind() == Kind::Diagonal)
    {
        maxFanout = 1;
        diagonal = injective = true;
        return;
    }
    if (gate.OperationKind() == Kind::Controlled)
    {
        const auto &m = gate.matrix;
        const bool a = RetainGateCoefficient(m[0], epsilon), b = RetainGateCoefficient(m[transposed ? 2 : 1], epsilon),
                   c = RetainGateCoefficient(m[transposed ? 1 : 2], epsilon), d = RetainGateCoefficient(m[3], epsilon);
        diagonal = !b && !c;
        injective = !(a && c) && !(b && d) && !(a && b) && !(c && d);
        maxFanout =
            std::max<unsigned>(!gate.qubits.empty(), std::max(unsigned(a) + unsigned(c), unsigned(b) + unsigned(d)));
        return;
    }
    // The transpose of a rotation only changes the sign of a sine.
    const bool z = gate.Axis() == Gates::RotationAxis::Z;
    bool anyBoth = false, anySine = false, anyKept = false;
    const Gates::StructuredGateInfo *info = gate.Info(false);
    if (!info)
    {
        const size_t size = gate.TableSize(), scanned = std::min<size_t>(size, 16);
        for (size_t selector = 0; selector < scanned; ++selector)
        {
            double c, s;
            gate.Rotation(selector, nullptr, c, s);
            if (z)
                anyKept = anyKept || RetainGateCoefficient({c, s}, epsilon);
            else
            {
                const bool keepC = RetainGateCoefficient(c, epsilon), keepS = RetainGateCoefficient(s, epsilon);
                anyBoth = anyBoth || (keepC && keepS);
                anySine = anySine || keepS;
                anyKept = anyKept || keepC || keepS;
            }
        }
        if (scanned < size && !(z ? anyKept : anyBoth))
            info = gate.Info(true);
    }
    if (info)
    {
        const bool exact = epsilon == 0.;
        if (z)
            anyKept = exact || info->maxNormSquare > epsilon;
        else
        {
            anyBoth = exact ? info->anyBothNonzero : info->maxMinSquare > epsilon;
            anySine = exact ? info->anySineNonzero : info->maxSineSquare > epsilon;
            anyKept = exact || info->maxMaxSquare > epsilon;
        }
    }
    diagonal = z || !anySine;
    injective = z || !anyBoth;
    maxFanout = !z && anyBoth ? 2 : unsigned(anyKept);
}

// The largest number of outputs of one basis state under the gate, or under
// its transpose, read from the matrix or the compact operation without
// building transition tables or a transposed gate.
inline unsigned GateFanout(const QC::Gates::AppliedGate<> &gate, double epsilon, bool transposed)
{
    if (const CompactGate *compact = gate.getCompactPayload())
    {
        unsigned fanout = 0;
        bool diagonal, injective;
        ClassifyCompact(*compact, epsilon, transposed, fanout, diagonal, injective);
        return fanout;
    }
    const auto &matrix = gate.getRawOperatorMatrix();
    const Eigen::Index outerSize = transposed ? matrix.rows() : matrix.cols();
    const Eigen::Index innerSize = transposed ? matrix.cols() : matrix.rows();
    unsigned fanout = 0;
    for (Eigen::Index outer = 0; outer < outerSize; ++outer)
    {
        unsigned count = 0;
        for (Eigen::Index inner = 0; inner < innerSize; ++inner)
            count += RetainGateCoefficient(transposed ? matrix(outer, inner) : matrix(inner, outer), epsilon);
        fanout = std::max(fanout, count);
    }
    return fanout;
}

// This is built from the current matrix and simulator threshold, not the
// legacy gate predicates. Each input column lists only retained transitions.
struct SparseGate
{
    // What every step reads comes first, on the cache lines of the counts.
    std::array<unsigned, 8> counts{};
    std::array<size_t, 3> qubits;
    // Set for a compact operation on no qubit or on more than three, which is
    // applied from its tables. The circuit owns this immutable payload. Circuit
    // copies share ownership, so cached transitions remain valid across
    // simulator copies as well.
    const CompactGate *structured = nullptr;
    unsigned arity;
    unsigned maxFanout = 0;
    bool diagonal = true;
    bool injective = true;
    bool distinctQubits = true;
    // Only the first counts[column] transitions of a column are ever written
    // or read. The tables are therefore left unset: clearing them costs more
    // than listing a small gate, and a table is built for every gate applied
    // outside a compiled circuit. A coefficient is its real and imaginary part.
    std::array<std::array<unsigned, 8>, 8> rows;
    std::array<std::array<std::array<double, 2>, 8>, 8> values;

    SparseGate(const QC::Gates::AppliedGate<> &gate, double epsilon)
        : qubits{gate.getQubit1(), gate.getQubit2(), gate.getQubit3()}
    {
        if (const CompactGate *compact = gate.getCompactPayload())
        {
            BuildCompact(*compact, epsilon);
            return;
        }
        arity = static_cast<unsigned>(gate.getQubitsNumber());
        if (arity == 0 || arity > 3)
            throw std::invalid_argument("Path integral gates must act on one to three qubits");
        const auto &matrix = gate.getRawOperatorMatrix();
        const unsigned dimension = 1u << arity;
        if (matrix.rows() != dimension || matrix.cols() != dimension)
            throw std::invalid_argument("Invalid path integral gate dimensions");
        // Most entries of a gate matrix are zero, so the test of one entry is
        // the cost that matters here. It is a single comparison in loops of
        // constant bounds: one instance per matrix size and per retention rule.
        const std::complex<double> *const data = matrix.data();
        const auto list = [&](const auto &keep) {
            switch (arity)
            {
            case 1:
                ListDense<2>(data, keep);
                break;
            case 2:
                ListDense<4>(data, keep);
                break;
            default:
                ListDense<8>(data, keep);
                break;
            }
        };
        if (epsilon == 0.)
            list([](const std::complex<double> &value) {
                // value != 0 as one integer test: some bit other than a sign is set.
                uint64_t parts[2];
                std::memcpy(parts, &value, sizeof(parts));
                return ((parts[0] | parts[1]) << 1) != 0;
            });
        else
            list([epsilon](const std::complex<double> &value) { return std::norm(value) > epsilon; });
        // Repeated target indices can destroy injectivity in the global state.
        for (unsigned i = 0; i < arity; ++i)
            for (unsigned j = i + 1; j < arity; ++j)
                if (qubits[i] == qubits[j])
                    injective = distinctQubits = false;
    }

    // The transitions of `source` whose row has the bits of `value` at the
    // positions of `mask`: the table for a step whose outputs are fixed there.
    // Only a gate that lists its transitions has any.
    SparseGate(const SparseGate &source, unsigned mask, unsigned value)
        : qubits(source.qubits), structured(source.structured), arity(source.arity), diagonal(source.diagonal),
          injective(source.injective), distinctQubits(source.distinctQubits)
    {
        const unsigned dimension = structured ? 0u : 1u << arity;
        for (unsigned column = 0; column < dimension; ++column)
        {
            unsigned retained = 0;
            for (unsigned transition = 0; transition < source.counts[column]; ++transition)
                if ((source.rows[column][transition] & mask) == value)
                {
                    rows[column][retained] = source.rows[column][transition];
                    values[column][retained++] = source.values[column][transition];
                }
            counts[column] = retained;
            maxFanout = std::max(maxFanout, retained);
        }
    }

    // Copies carry the unset transitions as raw bytes.
    SparseGate(const SparseGate &other) noexcept
        : counts(other.counts), qubits(other.qubits), structured(other.structured), arity(other.arity),
          maxFanout(other.maxFanout), diagonal(other.diagonal), injective(other.injective),
          distinctQubits(other.distinctQubits)
    {
        std::memcpy(&rows, &other.rows, sizeof(rows));
        std::memcpy(&values, &other.values, sizeof(values));
    }

    SparseGate &operator=(const SparseGate &other) noexcept
    {
        if (this != &other)
        {
            qubits = other.qubits;
            counts = other.counts;
            arity = other.arity;
            maxFanout = other.maxFanout;
            diagonal = other.diagonal;
            injective = other.injective;
            distinctQubits = other.distinctQubits;
            structured = other.structured;
            std::memcpy(&rows, &other.rows, sizeof(rows));
            std::memcpy(&values, &other.values, sizeof(values));
        }
        return *this;
    }

    std::complex<double> Value(unsigned column, unsigned transition) const
    {
        return {values[column][transition][0], values[column][transition][1]};
    }

    // The retained entries of a Dimension x Dimension matrix, column by column
    // and in increasing row order. The bounds are constants so that the rows
    // of a column are tested by straight-line code.
    template <unsigned Dimension, class Keep> void ListDense(const std::complex<double> *data, const Keep &keep)
    {
        constexpr bool rowMajor = Eigen::MatrixXcd::IsRowMajor != 0;
        unsigned usedRows = 0;
        for (unsigned column = 0; column < Dimension; ++column)
        {
            for (unsigned row = 0; row < Dimension; ++row)
            {
                const std::complex<double> value = data[rowMajor ? row * Dimension + column : column * Dimension + row];
                if (!keep(value))
                    continue;
                const unsigned index = counts[column]++;
                rows[column][index] = row;
                values[column][index] = {value.real(), value.imag()};
                diagonal = diagonal && row == column;
                if (index || (usedRows & (1u << row)))
                    injective = false;
                usedRows |= 1u << row;
            }
            maxFanout = std::max(maxFanout, counts[column]);
        }
    }

    // A compact operation on one to three qubits is listed like a dense gate
    // and runs through the same kernels. Local bit 0 is its target and the
    // controls follow; a diagonal keeps the order of its table bits, which is
    // the order of the gate's three qubit indices. The coefficients are those
    // of the compact operation itself. Any other compact operation is applied
    // from its tables and only classified here. Kept out of line: a dense
    // gate pays one test of its payload pointer for this path.
#if defined(_MSC_VER)
    __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    __attribute__((noinline))
#endif
    void
    BuildCompact(const CompactGate &gate, double epsilon)
    {
        arity = static_cast<unsigned>(gate.Arity());
        if (arity == 0 || arity > 3)
        {
            structured = &gate;
            ClassifyCompact(gate, epsilon, false, maxFanout, diagonal, injective);
            return;
        }
        unsigned usedRows = 0;
        const auto add = [&](unsigned column, unsigned row, const std::complex<double> &value) {
            if (!RetainGateCoefficient(value, epsilon))
                return;
            const unsigned index = counts[column]++;
            rows[column][index] = row;
            values[column][index] = {value.real(), value.imag()};
            diagonal = diagonal && row == column;
            if (index || (usedRows & (1u << row)))
                injective = false;
            usedRows |= 1u << row;
        };
        const unsigned dimension = 1u << arity;
        if (gate.OperationKind() == Gates::StructuredGate::Kind::Diagonal)
            for (unsigned column = 0; column < dimension; ++column)
                add(column, column, gate.Diagonal(column));
        else
        {
            // At most four angles: keep their cosines and sines with the tables.
            const auto *coefficients =
                gate.Coefficients(gate.OperationKind() == Gates::StructuredGate::Kind::MultiplexedRotation);
            for (unsigned selector = 0; selector < dimension / 2; ++selector)
            {
                const auto block = gate.Block(selector, coefficients);
                for (unsigned bit = 0; bit < 2; ++bit)
                {
                    add(2 * selector + bit, 2 * selector, block[bit]);
                    add(2 * selector + bit, 2 * selector + 1, block[2 + bit]);
                }
            }
        }
        for (unsigned column = 0; column < dimension; ++column)
            maxFanout = std::max(maxFanout, counts[column]);
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
