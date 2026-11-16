#pragma once

#include "GateAngles.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace QC
{
namespace Gates
{

enum class RotationAxis
{
    X,
    Y,
    Z
};

// The checks of a compact operation, shared by StructuredGate and its view. With
// tables = false the entries of a diagonal or of an angle table are not read:
// for tables that were validated when their gate was made, only the qubits and
// the sizes remain to be checked, which costs the number of qubits.
template <class Operation> void ValidateStructuredOperation(const Operation &operation, size_t nrQubits, bool tables);

// A compact applied operation. Local bit i is qubits[i]; for controlled operations
// qubits contains only the controls. No exponentially sized dense matrix is stored.
struct StructuredGate
{
    enum class Kind
    {
        Diagonal,
        Controlled,
        MultiplexedRotation
    };
    Kind kind = Kind::Diagonal;
    std::vector<size_t> qubits;
    size_t target = 0;
    size_t controlState = std::numeric_limits<size_t>::max();
    std::vector<std::complex<double>> diagonal;
    std::array<std::complex<double>, 4> matrix{}; // row major
    std::vector<double> angles;
    RotationAxis axis = RotationAxis::X;

    void Validate(size_t nrQubits) const
    {
        ValidateStructuredOperation(*this, nrQubits, true);
    }

    StructuredGate Adjoint() const
    {
        StructuredGate inverse(*this);
        if (kind == Kind::Diagonal)
            for (auto &value : inverse.diagonal)
                value = std::conj(value);
        else if (kind == Kind::Controlled)
        {
            inverse.matrix = {std::conj(matrix[0]), std::conj(matrix[2]), std::conj(matrix[1]), std::conj(matrix[3])};
        }
        else
            for (double &angle : inverse.angles)
                angle = -angle;
        return inverse;
    }

    // Number of qubits the operation acts on, the target included.
    size_t Arity() const
    {
        return qubits.size() + (kind == Kind::Diagonal ? 0 : 1);
    }

    // The qubits in the local bit order of DenseMatrix: the target first, then
    // the controls; for a diagonal, the order of its table bits.
    std::vector<size_t> OrderedQubits() const
    {
        std::vector<size_t> result;
        result.reserve(Arity());
        if (kind != Kind::Diagonal)
            result.push_back(target);
        result.insert(result.end(), qubits.begin(), qubits.end());
        return result;
    }

    // The dense matrix of a compact operation on few qubits, for the backends
    // that take local matrices. Local bit i is OrderedQubits()[i], the order in
    // which ApplyGate(gate, q0, q1, q2) takes its qubits. The coefficients are
    // those of the compact kernels: the same quarter-turn values, no rounding
    // beyond the trigonometry itself. The caller bounds the size.
    template <class MatrixClass> MatrixClass DenseMatrix(size_t maxQubits = 3) const
    {
        const size_t arity = Arity();
        if (arity == 0 || arity > maxQubits || arity >= std::numeric_limits<size_t>::digits)
            throw std::invalid_argument("This compact gate acts on too many qubits for a dense matrix here");
        const size_t dimension = size_t{1} << arity;
        MatrixClass result = MatrixClass::Zero(dimension, dimension);
        if (kind == Kind::Diagonal)
        {
            for (size_t i = 0; i < dimension; ++i)
                result(i, i) = diagonal[i];
            return result;
        }
        const size_t selected =
            controlState == std::numeric_limits<size_t>::max() ? (size_t{1} << qubits.size()) - 1 : controlState;
        for (size_t selector = 0; selector < dimension / 2; ++selector)
        {
            std::array<std::complex<double>, 4> block{{1., 0., 0., 1.}};
            if (kind == Kind::Controlled)
            {
                if (selector == selected)
                    block = matrix;
            }
            else
            {
                double c, s;
                detail::RotationHalfAngleCosSin(angles[selector], c, s);
                if (axis == RotationAxis::X)
                    block = {c, std::complex<double>(0., -s), std::complex<double>(0., -s), c};
                else if (axis == RotationAxis::Y)
                    block = {c, -s, s, c};
                else
                    block = {std::complex<double>(c, -s), 0., 0., std::complex<double>(c, s)};
            }
            for (size_t row = 0; row < 2; ++row)
                for (size_t col = 0; col < 2; ++col)
                    result(2 * selector + row, 2 * selector + col) = block[2 * row + col];
        }
        return result;
    }
};

// What a kernel reads of a compact operation: the fields of StructuredGate with
// the two tables as (pointer, size). An entry point that is handed tables by
// its caller applies them through a view and copies nothing; a StructuredGate
// converts to its view. A view borrows the qubit list and the tables: it must
// not outlive them.
struct StructuredGateView
{
    template <class T> struct Table
    {
        const T *first = nullptr;
        size_t count = 0;

        size_t size() const
        {
            return count;
        }

        const T *data() const
        {
            return first;
        }

        const T *begin() const
        {
            return first;
        }

        const T *end() const
        {
            return first + count;
        }

        const T &operator[](size_t index) const
        {
            return first[index];
        }
    };

    StructuredGate::Kind kind = StructuredGate::Kind::Diagonal;
    const std::vector<size_t> &qubits;
    size_t target = 0;
    size_t controlState = std::numeric_limits<size_t>::max();
    Table<std::complex<double>> diagonal;
    std::array<std::complex<double>, 4> matrix{}; // row major
    Table<double> angles;
    RotationAxis axis = RotationAxis::X;

    StructuredGateView(const StructuredGate &gate)
        : kind(gate.kind), qubits(gate.qubits), target(gate.target), controlState(gate.controlState),
          diagonal{gate.diagonal.data(), gate.diagonal.size()}, matrix(gate.matrix),
          angles{gate.angles.data(), gate.angles.size()}, axis(gate.axis)
    {
    }

    // A diagonal: values[s] multiplies the basis states whose bit qubits[i] is bit i of s.
    static StructuredGateView Diagonal(const std::vector<std::complex<double>> &values,
                                       const std::vector<size_t> &qubits)
    {
        StructuredGateView view(qubits);
        view.diagonal = {values.data(), values.size()};
        return view;
    }

    // A row-major 2x2 matrix on the target for one state of the controls.
    static StructuredGateView Controlled(const std::array<std::complex<double>, 4> &matrix, size_t target,
                                         const std::vector<size_t> &controls, size_t controlState)
    {
        StructuredGateView view(controls);
        view.kind = StructuredGate::Kind::Controlled;
        view.target = target;
        view.controlState = controlState;
        view.matrix = matrix;
        return view;
    }

    // A rotation of the target by angles[s] for the control state s.
    static StructuredGateView Multiplexed(RotationAxis axis, const std::vector<double> &angles, size_t target,
                                          const std::vector<size_t> &controls)
    {
        StructuredGateView view(controls);
        view.kind = StructuredGate::Kind::MultiplexedRotation;
        view.axis = axis;
        view.target = target;
        view.angles = {angles.data(), angles.size()};
        return view;
    }

    void Validate(size_t nrQubits) const
    {
        ValidateStructuredOperation(*this, nrQubits, true);
    }

    // The same without reading the table entries, for tables validated before.
    void ValidateQubits(size_t nrQubits) const
    {
        ValidateStructuredOperation(*this, nrQubits, false);
    }

    // An owning copy, for a caller that keeps the operation (a recording).
    StructuredGate Materialize() const
    {
        StructuredGate gate;
        gate.kind = kind;
        gate.qubits = qubits;
        gate.target = target;
        gate.controlState = controlState;
        gate.diagonal.assign(diagonal.begin(), diagonal.end());
        gate.matrix = matrix;
        gate.angles.assign(angles.begin(), angles.end());
        gate.axis = axis;
        return gate;
    }

  private:
    explicit StructuredGateView(const std::vector<size_t> &list) : qubits(list)
    {
    }
};

template <class Operation> void ValidateStructuredOperation(const Operation &operation, size_t nrQubits, bool tables)
{
    using Kind = StructuredGate::Kind;
    const auto &qubits = operation.qubits;
    if (qubits.size() >= std::numeric_limits<size_t>::digits)
        throw std::invalid_argument("Too many gate qubits");
    for (size_t i = 0; i < qubits.size(); ++i)
    {
        const size_t q = qubits[i];
        if (q >= nrQubits)
            throw std::invalid_argument("Qubit number is too high");
        if (std::find(qubits.begin(), qubits.begin() + i, q) != qubits.begin() + i)
            throw std::invalid_argument("Gate qubits must be distinct");
    }
    const size_t dimension = size_t{1} << qubits.size();
    if (operation.kind == Kind::Diagonal)
    {
        if (operation.diagonal.size() != dimension)
            throw std::invalid_argument("Diagonal size must be 2 to the number of qubits");
        if (tables)
            for (const auto &value : operation.diagonal)
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::invalid_argument("Diagonal entries must be finite");
        return;
    }
    if (operation.target >= nrQubits || std::find(qubits.begin(), qubits.end(), operation.target) != qubits.end())
        throw std::invalid_argument("Target must be in range and distinct from controls");
    if (operation.kind == Kind::Controlled)
    {
        if (operation.controlState != std::numeric_limits<size_t>::max() && operation.controlState >= dimension)
            throw std::invalid_argument("Control state does not fit its controls");
        for (const auto &value : operation.matrix)
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                throw std::invalid_argument("Controlled matrix entries must be finite");
    }
    else if (operation.kind == Kind::MultiplexedRotation)
    {
        if (operation.axis != RotationAxis::X && operation.axis != RotationAxis::Y && operation.axis != RotationAxis::Z)
            throw std::invalid_argument("Invalid rotation axis");
        if (operation.angles.size() != dimension)
            throw std::invalid_argument("Angle table size must be 2 to the number of controls");
        if (tables)
            for (double angle : operation.angles)
                if (!std::isfinite(angle))
                    throw std::invalid_argument("Rotation angles must be finite");
    }
    else
        throw std::invalid_argument("Invalid structured operation");
}

// Where the bits of a selector live in a basis state kept as 64-bit words:
// qubit q is bit q % 64 of word q / 64, and selector bit i belongs to
// qubits[i]. Consecutive selector bits that are consecutive bits of one word
// are read with a single shift and mask, and a required control pattern is
// compared one word at a time, so neither costs one test per qubit.
struct StructuredSelector
{
    struct Run
    {
        size_t word;
        uint64_t mask;
        unsigned shift;
        unsigned out;
    };

    struct Match
    {
        size_t word;
        uint64_t mask;
        uint64_t value;
    };

    StructuredSelector() = default;

    // wanted is the selector value that Matches accepts.
    StructuredSelector(const std::vector<size_t> &qubits, size_t wanted)
    {
        for (size_t i = 0; i < qubits.size();)
        {
            const size_t first = qubits[i], word = first / 64;
            size_t length = 1;
            while (i + length < qubits.size() && qubits[i + length] == first + length && (first + length) / 64 == word)
                ++length;
            runs.push_back(
                {word, length == 64 ? ~uint64_t{0} : (uint64_t{1} << length) - 1, unsigned(first % 64), unsigned(i)});
            i += length;
        }
        for (size_t i = 0; i < qubits.size(); ++i)
        {
            const size_t word = qubits[i] / 64;
            const uint64_t bit = uint64_t{1} << (qubits[i] % 64);
            auto found =
                std::find_if(matches.begin(), matches.end(), [word](const Match &m) { return m.word == word; });
            if (found == matches.end())
            {
                matches.push_back({word, 0, 0});
                found = matches.end() - 1;
            }
            found->mask |= bit;
            if ((wanted >> i) & 1)
                found->value |= bit;
        }
    }

    // The words must cover every qubit; callers check the register size once.
    size_t Gather(const uint64_t *words) const noexcept
    {
        size_t selector = 0;
        for (const Run &run : runs)
            selector |= size_t((words[run.word] >> run.shift) & run.mask) << run.out;
        return selector;
    }

    bool Matches(const uint64_t *words) const noexcept
    {
        for (const Match &match : matches)
            if ((words[match.word] & match.mask) != match.value)
                return false;
        return true;
    }

    std::vector<Run> runs;
    std::vector<Match> matches;
};

// The half-angle cosine and sine of every angle of a multiplexed rotation.
using StructuredRotationCoefficients = std::vector<std::array<double, 2>>;

inline StructuredRotationCoefficients RotationCoefficients(const StructuredGate &gate)
{
    StructuredRotationCoefficients coefficients(gate.angles.size());
    for (size_t i = 0; i < gate.angles.size(); ++i)
        detail::RotationHalfAngleCosSin(gate.angles[i], coefficients[i][0], coefficients[i][1]);
    return coefficients;
}

// What one pass over the table of a compact operation tells about it, without
// keeping anything of the table's size: the extremes that decide, for any
// retention threshold, how many outputs one basis state can reach. A backend
// that needs them computes them once and keeps them with the tables;
// applying a StructuredGate directly never requires them.
struct StructuredGateInfo
{
    explicit StructuredGateInfo(const StructuredGate &gate)
    {
        if (gate.kind == StructuredGate::Kind::Diagonal)
        {
            anyNonzero = std::any_of(gate.diagonal.begin(), gate.diagonal.end(), [](const std::complex<double> &value) {
                return value != std::complex<double>(0.);
            });
            return;
        }
        if (gate.kind != StructuredGate::Kind::MultiplexedRotation)
            return;
        for (size_t i = 0; i < gate.angles.size(); ++i)
        {
            double c, s;
            detail::RotationHalfAngleCosSin(gate.angles[i], c, s);
            // The squares as a retention test of one real coefficient forms them.
            const double c2 = std::norm(std::complex<double>(c)), s2 = std::norm(std::complex<double>(s));
            maxMinSquare = std::max(maxMinSquare, std::min(c2, s2));
            maxMaxSquare = std::max(maxMaxSquare, std::max(c2, s2));
            maxSineSquare = std::max(maxSineSquare, s2);
            maxNormSquare = std::max(maxNormSquare, std::norm(std::complex<double>(c, s)));
            anyBothNonzero = anyBothNonzero || (c != 0. && s != 0.);
            anySineNonzero = anySineNonzero || s != 0.;
        }
    }

    // Largest over the selectors of min(c^2, s^2), of max(c^2, s^2), of s^2 and of c^2 + s^2,
    // for the half-angle cosine c and sine s of a multiplexed rotation.
    double maxMinSquare = 0., maxMaxSquare = 0., maxSineSquare = 0., maxNormSquare = 0.;
    bool anyBothNonzero = false, anySineNonzero = false; // exact tests of the same coefficients
    bool anyNonzero = false;                             // a diagonal with a nonzero entry
};

} // namespace Gates
} // namespace QC
