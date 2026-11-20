#pragma once

#define _USE_MATH_DEFINES
#include <Eigen/Eigen>
#include <math.h>

#include "GateAngles.h"
#include "GateStructure.h"
#include "StructuredGate.h"
#include "UnitaryGate.h"
#include <array>
#include <atomic>
#include <cassert>
#include <cstring>
#include <memory>
#include <type_traits>
#include <typeinfo>

namespace QC
{
namespace Gates
{

// Whether a matrix is unitary, to the rounding a product of many gates carries (1e-12, the tolerance of the
// three-qubit gates of the matrix product state): its columns are orthonormal. Spelled out on the entries, so that
// the small matrices of the usual gates need no temporary.
template <class MatrixClass> bool IsUnitaryMatrix(const MatrixClass &matrix)
{
    const Eigen::Index size = matrix.rows();
    if (size != matrix.cols())
        return false;
    if (size == 2) // every one-qubit matrix: written out
    {
        const std::complex<double> a = matrix(0, 0), b = matrix(0, 1), c = matrix(1, 0), d = matrix(1, 1);
        return std::abs(std::norm(a) + std::norm(c) - 1.) <= 1E-12 &&
               std::abs(std::norm(b) + std::norm(d) - 1.) <= 1E-12 &&
               std::norm(std::conj(a) * b + std::conj(c) * d) <= 1E-24;
    }
    for (Eigen::Index j = 0; j < size; ++j)
        for (Eigen::Index k = j; k < size; ++k)
        {
            std::complex<double> product = 0.;
            for (Eigen::Index i = 0; i < size; ++i)
                product += std::conj(std::complex<double>(matrix(i, j))) * std::complex<double>(matrix(i, k));
            // A NaN fails either comparison.
            if (!(j == k ? std::abs(product.real() - 1.) <= 1E-12 : std::norm(product) <= 1E-24))
                return false;
        }
    return true;
}

template <class MatrixClass = Eigen::MatrixXcd> class QuantumGate
{
  public:
    // controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for
    // gates with a lower number of qubits
    virtual MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0,
                                          size_t controllingQubit2 = 0) const = 0;
    virtual ~QuantumGate() = default;

    virtual size_t getQubitsNumber() const
    {
        return 0; // 0 at this point means 'unknown', should be overriden if needed
    }

    // the following two functions could be used for some optimizatiions
    // only for reducing a constant... if it matters at one point I'll implement them
    // see QubitRegister::applyGate for details of the optimizations
    virtual bool isControlled() const
    {
        return false;
    }

    virtual bool isControlQubit(size_t /*qubit*/) const
    {
        return false;
    }

    virtual bool isDiagonal() const
    {
        return false;
    }

    virtual bool isAntidiagonal() const
    {
        return false;
    }

    virtual bool isSwapGate() const
    {
        return false;
    }

    virtual bool IsISwapGate() const
    {
        return false;
    }

    virtual bool IsISwapDagGate() const
    {
        return false;
    }

    // Whether the matrix of the gate is unitary. The gates defined by name (X, H, CX, the rotations and so on) are
    // by construction, and nothing is tested for them. A gate built on a matrix that was handed over (AppliedGate)
    // tests it, once; an operator that is no gate (a projection) says so. A simulator for which it matters asks: the
    // matrix product state applies what is not unitary as a filter, see MPSSimulatorImpl::ApplyGate.
    virtual bool isUnitary() const
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class QuantumGateWithOp : public QuantumGate<MatrixClass>
{
  public:
    using BaseClass = QuantumGate<MatrixClass>;

    QuantumGateWithOp(const QuantumGateWithOp &other) : operatorMat(other.operatorMat)
    {
    }

    QuantumGateWithOp(QuantumGateWithOp &&other) : operatorMat(std::move(other.operatorMat))
    {
    }

    QuantumGateWithOp(const MatrixClass &U) : operatorMat(U)
    {
    }

    QuantumGateWithOp(MatrixClass &&U) : operatorMat(std::move(U))
    {
    }

    QuantumGateWithOp &operator=(const QuantumGateWithOp &other)
    {
        if (this == &other)
            return *this;
        if (other.getStructuredGate())
            throw std::invalid_argument("Assign compact gates through AppliedGate to preserve their representation");
        operatorMat = other.operatorMat;
        OperatorChanged();
        return *this;
    }

    QuantumGateWithOp &operator=(QuantumGateWithOp &&other)
    {
        if (this == &other)
            return *this;
        if (other.getStructuredGate())
            throw std::invalid_argument("Assign compact gates through AppliedGate to preserve their representation");
        operatorMat = std::move(other.operatorMat);
        OperatorChanged();
        return *this;
    }

    QuantumGateWithOp &operator=(const MatrixClass &U)
    {
        operatorMat = U;
        OperatorChanged();
        return *this;
    }

    QuantumGateWithOp &operator=(MatrixClass &&U)
    {
        operatorMat = std::move(U);
        OperatorChanged();
        return *this;
    }

    virtual GateStructure getStructure() const
    {
        return ClassifyGateMatrix(operatorMat);
    }

    const MatrixClass &getRawOperatorMatrix() const
    {
        // Compact AppliedGate payloads use an empty matrix sentinel. Keep the
        // original nonvirtual dense accessor: only that exceptional sentinel
        // needs a representation query, including calls through a base reference.
        // The sentinel is a matrix without rows, so a dense gate pays one
        // comparison of a value its caller reads anyway.
        if (operatorMat.rows() == 0 && getStructuredGate())
            throw std::invalid_argument(
                "A compact applied gate has no dense operator matrix; use a structured gate backend");
        return operatorMat;
    }

    virtual UnitaryGateDescriptor getUnitaryDescriptor() const
    {
        return {};
    }

    virtual const StructuredGate *getStructuredGate() const
    {
        return nullptr;
    }

    size_t getQubitsNumber() const override
    {
        // Common arities avoid a transcendental call on every AppliedGate.
        const auto rows = operatorMat.rows();
        if (rows == 2)
            return 1;
        if (rows == 4)
            return 2;
        if (rows == 8)
            return 3;
        return rows > 0 ? static_cast<size_t>(log2(rows)) : 0;
    }

    void setOperator(const MatrixClass &U)
    {
        assert(U.rows() == U.cols());
        // assert(U.rows() == 1ULL << getQubitsNumber());

        operatorMat = U;
        OperatorChanged();
    }

  protected:
    // Derived caches refresh eagerly, including mutations through a base reference.
    virtual void OperatorChanged()
    {
    }

    MatrixClass operatorMat;
};

template <class MatrixClass = Eigen::MatrixXcd> class SingleQubitGate : public QuantumGateWithOp<MatrixClass>
{
  public:
    using BaseClass = QuantumGateWithOp<MatrixClass>;

    SingleQubitGate() noexcept : BaseClass(MatrixClass::Zero(2, 2))
    {
    }

    SingleQubitGate(const SingleQubitGate &other) noexcept : BaseClass(other.operatorMat)
    {
    }

    SingleQubitGate(SingleQubitGate &&other) noexcept : BaseClass(std::move(other.operatorMat))
    {
    }

    SingleQubitGate(const MatrixClass &U) noexcept : BaseClass(U)
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 2);
    }

    SingleQubitGate(MatrixClass &&U) noexcept : BaseClass(std::move(U))
    {
        assert(BaseClass::operatorMat.rows() == BaseClass::operatorMat.cols());
        assert(BaseClass::operatorMat.rows() == 2);
    }

    SingleQubitGate &operator=(const SingleQubitGate &other)
    {
        assert(other.operatorMat.rows() == 2);
        assert(other.operatorMat.cols() == 2);
        BaseClass::operator=(other.operatorMat);
        return *this;
    }

    SingleQubitGate &operator=(SingleQubitGate &&other) noexcept
    {
        assert(other.operatorMat.rows() == 2);
        assert(other.operatorMat.cols() == 2);
        BaseClass::operator=(std::move(other.operatorMat));
        return *this;
    }

    SingleQubitGate &operator=(const MatrixClass &U)
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 2);
        BaseClass::operator=(U);
        return *this;
    }

    SingleQubitGate &operator=(MatrixClass &&U) noexcept
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 2);
        BaseClass::operator=(std::move(U));
        return *this;
    }

    size_t getQubitsNumber() const override
    {
        return 1;
    }

    // controllingQubit1 and controllingQubit2 are ignored, they will be used for two (only controllingQubit1) and three
    // qubit gates this is not used anymore, instead there is a more efficient implementation in
    // QubitRegister::applyGate even this matrix could be constructed more efficiently in a similar manner (see
    // QubitRegister::applyGate for details) but I won't bother, I'll keep this different implementation because it
    // might be more clear - construction by tensor product and also could be useful for debugging in the case the
    // optimized version has a problem
    MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0,
                                  size_t controllingQubit2 = 0) const override
    {
        const size_t nrBasisStates = 1ULL << nrQubits;
        MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

        const size_t qubitBit = 1ULL << qubit;

        // since this can be quite condensed, here is a description:
        // this just computes the tensor product between the 2x2 operator matrix for the qubit and the identity
        // operators for the other qubits

        // for the two and three qubits operators it's more complex but analogous

        for (size_t i = 0; i < nrBasisStates; ++i)
        {
            const size_t ind1 = i | qubitBit;
            for (size_t j = 0; j < nrBasisStates; ++j)
                if (ind1 == (j | qubitBit)) // this is just a lot of delta 'functions' for the 'other' qubits, the bit
                                            // was forced to 1 for the qubit with this operator
                    extOperatorMat(i, j) = BaseClass::operatorMat(
                        i & qubitBit ? 1 : 0,
                        j & qubitBit ? 1 : 0); // pick the correct matrix element for the qubit operator
        }

        return extOperatorMat;
    };
};

template <class MatrixClass = Eigen::MatrixXcd> class TwoQubitsGate : public QuantumGateWithOp<MatrixClass>
{
  public:
    using BaseClass = QuantumGateWithOp<MatrixClass>;

    TwoQubitsGate() noexcept : BaseClass(MatrixClass::Identity(4, 4))
    {
    }

    TwoQubitsGate(const TwoQubitsGate &other) noexcept : BaseClass(other.operatorMat)
    {
    }

    TwoQubitsGate(TwoQubitsGate &&other) noexcept : BaseClass(std::move(other.operatorMat))
    {
    }

    TwoQubitsGate(const MatrixClass &U) noexcept : BaseClass(U)
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 4);
    }

    TwoQubitsGate(MatrixClass &&U) noexcept : BaseClass(std::move(U))
    {
        assert(BaseClass::operatorMat.rows() == BaseClass::operatorMat.cols());
        assert(BaseClass::operatorMat.rows() == 4);
    }

    TwoQubitsGate &operator=(const TwoQubitsGate &other)
    {
        assert(other.operatorMat.rows() == 4);
        assert(other.operatorMat.cols() == 4);
        BaseClass::operator=(other.operatorMat);
        return *this;
    }

    TwoQubitsGate &operator=(TwoQubitsGate &&other) noexcept
    {
        assert(other.operatorMat.rows() == 4);
        assert(other.operatorMat.cols() == 4);
        BaseClass::operator=(std::move(other.operatorMat));
        return *this;
    }

    TwoQubitsGate &operator=(const MatrixClass &U)
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 4);
        BaseClass::operator=(U);
        return *this;
    }

    TwoQubitsGate &operator=(MatrixClass &&U) noexcept
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 4);
        BaseClass::operator=(std::move(U));
        return *this;
    }

    size_t getQubitsNumber() const override
    {
        return 2;
    }

    // controllingQubit2 is ignored, it will be used for three qubit gates only
    // this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
    // even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for
    // details) but I won't bother, I'll keep this different implementation because it might be more clear -
    // construction by tensor product and also could be useful for debugging in the case the optimized version has a
    // problem
    MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0,
                                  size_t controllingQubit2 = 0) const override
    {
        assert(qubit != controllingQubit1);

        const size_t nrBasisStates = 1ULL << nrQubits;
        MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

        const size_t qubitBit = 1ULL << qubit;
        const size_t ctrlQubitBit = 1ULL << controllingQubit1;
        const size_t mask = qubitBit | ctrlQubitBit;

        // computing the tensor product between the gate matrix and identity operators for the other qubits

        for (size_t i = 0; i < nrBasisStates; ++i)
        {
            const size_t ind1 = i | mask;
            for (size_t j = 0; j < nrBasisStates; ++j)
                if (ind1 == (j | mask)) // the delta 'function'
                    extOperatorMat(i, j) = BaseClass::operatorMat((i & ctrlQubitBit ? 2 : 0) | (i & qubitBit ? 1 : 0),
                                                                  (j & ctrlQubitBit ? 2 : 0) | (j & qubitBit ? 1 : 0));
        }

        return extOperatorMat;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class ThreeQubitsGate : public QuantumGateWithOp<MatrixClass>
{
  public:
    using BaseClass = QuantumGateWithOp<MatrixClass>;

    ThreeQubitsGate() noexcept : BaseClass(MatrixClass::Identity(8, 8))
    {
    }

    ThreeQubitsGate(const ThreeQubitsGate &other) noexcept : BaseClass(other.operatorMat)
    {
    }

    ThreeQubitsGate(ThreeQubitsGate &&other) noexcept : BaseClass(std::move(other.operatorMat))
    {
    }

    ThreeQubitsGate(const MatrixClass &U) noexcept : BaseClass(U)
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 8);
    }

    ThreeQubitsGate(MatrixClass &&U) noexcept : BaseClass(std::move(U))
    {
        assert(BaseClass::operatorMat.rows() == BaseClass::operatorMat.cols());
        assert(BaseClass::operatorMat.rows() == 8);
    }

    ThreeQubitsGate &operator=(const ThreeQubitsGate &other)
    {
        assert(other.operatorMat.rows() == 8);
        assert(other.operatorMat.cols() == 8);
        BaseClass::operator=(other.operatorMat);
        return *this;
    }

    ThreeQubitsGate &operator=(ThreeQubitsGate &&other) noexcept
    {
        assert(other.operatorMat.rows() == 8);
        assert(other.operatorMat.cols() == 8);
        BaseClass::operator=(std::move(other.operatorMat));
        return *this;
    }

    ThreeQubitsGate &operator=(const MatrixClass &U)
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 8);
        BaseClass::operator=(U);
        return *this;
    }

    ThreeQubitsGate &operator=(MatrixClass &&U) noexcept
    {
        assert(U.rows() == U.cols());
        assert(U.rows() == 8);
        BaseClass::operator=(std::move(U));
        return *this;
    }

    size_t getQubitsNumber() const override
    {
        return 3;
    }

    // this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
    // even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for
    // details) but I won't bother, I'll keep this different implementation because it might be more clear -
    // construction by tensor product and also could be useful for debugging in the case the optimized version has a
    // problem
    MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0,
                                  size_t controllingQubit2 = 0) const override
    {
        assert(qubit != controllingQubit1 && controllingQubit1 != controllingQubit2);

        const size_t nrBasisStates = 1ULL << nrQubits;
        MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

        const size_t qubitBit = 1ULL << qubit;
        const size_t qubitBit2 = 1ULL << controllingQubit1;
        const size_t ctrlQubitBit = 1ULL << controllingQubit2;
        const size_t mask = qubitBit | qubitBit2 | ctrlQubitBit;

        // computing the tensor product between the gate matrix and identity operators for the other qubits

        for (size_t i = 0; i < nrBasisStates; ++i)
        {
            const size_t ind1 = i | mask;
            for (size_t j = 0; j < nrBasisStates; ++j)
                if (ind1 == (j | mask)) // the delta 'function'
                    extOperatorMat(i, j) = BaseClass::operatorMat(
                        (i & ctrlQubitBit ? 4 : 0) | (i & qubitBit2 ? 2 : 0) | (i & qubitBit ? 1 : 0),
                        (j & ctrlQubitBit ? 4 : 0) | (j & qubitBit2 ? 2 : 0) | (j & qubitBit ? 1 : 0));
        }

        return extOperatorMat;
    }
};

namespace AppliedGateDetail
{
// What an applied gate knows of its matrix: whether it is unitary, once it was tested or told. Copies carry it, and
// it is atomic because a gate that threads share may be asked by several of them (they would all find the same).
struct KnownUnitarity
{
    enum : uint8_t
    {
        NotTested,
        Unitary,
        NotUnitary
    };

    std::atomic<uint8_t> value{NotTested};
    KnownUnitarity() = default;

    KnownUnitarity(const KnownUnitarity &other) noexcept : value(other.value.load(std::memory_order_relaxed))
    {
    }

    KnownUnitarity &operator=(const KnownUnitarity &other) noexcept
    {
        value.store(other.value.load(std::memory_order_relaxed), std::memory_order_relaxed);
        return *this;
    }
};

// Dense raw gates keep a null pointer. Compact operations and explicitly
// typed provenance share the existing pointer-sized representation field.
struct Payload
{
    enum class Kind
    {
        Structured,
        Unitary
    };

    explicit Payload(Kind value) : kind(value)
    {
    }

    const Kind kind;
};

// The tables of a compact operation, immutable once validated. Every gate made
// from it by copying, by remapping its qubits or by taking its adjoint or its
// transpose shares this object, and with it the values derived from the tables.
struct StructuredCore
{
    explicit StructuredCore(StructuredGate &&value) : operation(std::move(value))
    {
    }

    StructuredCore(const StructuredCore &) = delete;
    StructuredCore &operator=(const StructuredCore &) = delete;

    ~StructuredCore()
    {
        delete info.load(std::memory_order_relaxed);
        delete coefficients.load(std::memory_order_relaxed);
    }

    // What a pass over the table tells about the operation, computed by the
    // first caller that asks to build it. Concurrent first calls may both
    // compute; one result is kept.
    const StructuredGateInfo *Info(bool build) const
    {
        return Derived(info, build, [this] { return new StructuredGateInfo(operation); });
    }

    // The cosines and sines of a multiplexed rotation, as large as its angle
    // table: built only for a caller that will read a good part of them.
    const StructuredRotationCoefficients *Coefficients(bool build) const
    {
        return Derived(coefficients, build,
                       [this] { return new StructuredRotationCoefficients(RotationCoefficients(operation)); });
    }

    const StructuredGate operation;

  private:
    template <class Value, class Make>
    static const Value *Derived(std::atomic<const Value *> &slot, bool build, const Make &make)
    {
        const Value *current = slot.load(std::memory_order_acquire);
        if (!current && build)
        {
            const Value *fresh = make();
            if (slot.compare_exchange_strong(current, fresh, std::memory_order_acq_rel, std::memory_order_acquire))
                current = fresh;
            else
                delete fresh;
        }
        return current;
    }

    mutable std::atomic<const StructuredGateInfo *> info{nullptr};
    mutable std::atomic<const StructuredRotationCoefficients *> coefficients{nullptr};
};

// A compact operation as one applied gate sees it: the shared tables, the
// qubits this gate acts on, the target matrix of a controlled operation, and
// whether the tables are read conjugated (the adjoint of a diagonal, or of a
// multiplexed rotation: its angles negated). Remapping qubits and taking the
// adjoint or the transpose therefore never copy or rescan a table.
struct StructuredPayload : Payload
{
    // A new operation; the caller has validated it.
    explicit StructuredPayload(StructuredGate &&value)
        : Payload(Kind::Structured), core(std::make_shared<const StructuredCore>(std::move(value))),
          qubits(core->operation.qubits), target(core->operation.target), matrix(core->operation.matrix),
          conjugated(false)
    {
        Bind();
    }

    // The tables of another payload under a different binding.
    StructuredPayload(const StructuredPayload &source, std::vector<size_t> boundQubits, size_t boundTarget,
                      const std::array<std::complex<double>, 4> &boundMatrix, bool conjugate)
        : Payload(Kind::Structured), core(source.core), qubits(std::move(boundQubits)), target(boundTarget),
          matrix(boundMatrix), conjugated(conjugate)
    {
        Bind();
    }

    StructuredPayload(const StructuredPayload &) = delete;
    StructuredPayload &operator=(const StructuredPayload &) = delete;

    ~StructuredPayload()
    {
        delete materialized.load(std::memory_order_relaxed);
    }

    StructuredGate::Kind OperationKind() const
    {
        return core->operation.kind;
    }

    RotationAxis Axis() const
    {
        return core->operation.axis;
    }

    size_t Arity() const
    {
        return qubits.size() + (OperationKind() == StructuredGate::Kind::Diagonal ? 0 : 1);
    }

    // Local bit order: the target first, then the controls; a diagonal keeps
    // the order of its table bits.
    size_t Qubit(size_t local) const
    {
        if (OperationKind() == StructuredGate::Kind::Diagonal)
            return qubits[local];
        return local == 0 ? target : qubits[local - 1];
    }

    // The selector for which the matrix of a controlled operation applies.
    size_t ControlValue() const
    {
        const size_t state = core->operation.controlState;
        return state == std::numeric_limits<size_t>::max() ? (size_t{1} << qubits.size()) - 1 : state;
    }

    std::complex<double> Diagonal(size_t selector) const
    {
        const std::complex<double> value = core->operation.diagonal[selector];
        return conjugated ? std::conj(value) : value;
    }

    size_t TableSize() const
    {
        return OperationKind() == StructuredGate::Kind::Diagonal ? core->operation.diagonal.size()
                                                                 : core->operation.angles.size();
    }

    // Half-angle cosine and sine of the rotation for one selector, from the
    // coefficient table when the caller has it.
    void Rotation(size_t selector, const StructuredRotationCoefficients *coefficients, double &c, double &s) const
    {
        if (coefficients)
        {
            c = (*coefficients)[selector][0];
            s = (*coefficients)[selector][1];
        }
        else
            detail::RotationHalfAngleCosSin(core->operation.angles[selector], c, s);
        if (conjugated)
            s = -s;
    }

    // The row-major 2x2 matrix applied to the target for one selector.
    std::array<std::complex<double>, 4> Block(size_t selector,
                                              const StructuredRotationCoefficients *coefficients = nullptr) const
    {
        using Complex = std::complex<double>;
        if (OperationKind() == StructuredGate::Kind::Controlled)
            return selector == ControlValue() ? matrix : std::array<Complex, 4>{{1., 0., 0., 1.}};
        double c, s;
        Rotation(selector, coefficients, c, s);
        if (Axis() == RotationAxis::X)
            return {{c, Complex(0., -s), Complex(0., -s), c}};
        if (Axis() == RotationAxis::Y)
            return {{c, -s, s, c}};
        return {{Complex(c, -s), 0., 0., Complex(c, s)}};
    }

    // The operation as the kernels read it: the shared tables under this gate's
    // binding, with nothing copied or built. For a conjugated payload the
    // tables are those of the adjoint (a diagonal read conjugated, the angles
    // negated), which the caller hands to the kernel together with the view;
    // the matrix of a controlled operation is already the one to apply.
    StructuredGateView View() const
    {
        const StructuredGate &tables = core->operation;
        if (tables.kind == StructuredGate::Kind::Diagonal)
            return StructuredGateView::Diagonal(tables.diagonal, qubits);
        if (tables.kind == StructuredGate::Kind::Controlled)
            return StructuredGateView::Controlled(matrix, target, qubits, tables.controlState);
        return StructuredGateView::Multiplexed(tables.axis, tables.angles, target, qubits);
    }

    const StructuredGateInfo *Info(bool build) const
    {
        return core->Info(build);
    }

    const StructuredRotationCoefficients *Coefficients(bool build) const
    {
        return core->Coefficients(build);
    }

    // The operation as a StructuredGate. A gate that still has the binding it
    // was created with returns the shared tables themselves; a remapped,
    // adjoint or transposed gate builds its own copy the first time a backend
    // that reads a StructuredGate asks for it.
    const StructuredGate &Gate() const
    {
        return plain ? core->operation : Materialized();
    }

    const std::shared_ptr<const StructuredCore> core;
    const std::vector<size_t> qubits; // the controls, or the table bits of a diagonal
    const size_t target;
    const std::array<std::complex<double>, 4> matrix; // row major
    const bool conjugated;
    bool plain = false;    // exactly the operation the tables were created with
    bool distinct = false; // no qubit appears twice
    // A multiplexer with so few angles that their cosines and sines are worth
    // keeping from the first use on, whatever reads them.
    bool fewAngles = false;
    size_t maxQubit = 0;
    StructuredSelector selector;

  private:
    void Bind()
    {
        const StructuredGate &operation = core->operation;
        plain = !conjugated && target == operation.target && matrix == operation.matrix && qubits == operation.qubits;
        distinct = true;
        for (size_t i = 0; i < qubits.size(); ++i)
        {
            maxQubit = std::max(maxQubit, qubits[i]);
            if (std::find(qubits.begin(), qubits.begin() + i, qubits[i]) != qubits.begin() + i)
                distinct = false;
        }
        if (operation.kind != StructuredGate::Kind::Diagonal)
        {
            maxQubit = std::max(maxQubit, target);
            if (std::find(qubits.begin(), qubits.end(), target) != qubits.end())
                distinct = false;
        }
        selector = StructuredSelector(qubits, ControlValue());
        fewAngles = operation.kind == StructuredGate::Kind::MultiplexedRotation && operation.angles.size() <= 64;
    }

    // Kept out of line: the accessors of a dense gate reach this only through
    // a test of its payload, and must stay as small as they were.
#if defined(_MSC_VER)
    __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    __attribute__((noinline))
#endif
    const StructuredGate &
    Materialized() const
    {
        const StructuredGate *current = materialized.load(std::memory_order_acquire);
        if (!current)
        {
            auto fresh = std::make_unique<StructuredGate>(core->operation);
            fresh->qubits = qubits;
            fresh->target = target;
            fresh->matrix = matrix;
            if (conjugated)
            {
                for (auto &value : fresh->diagonal)
                    value = std::conj(value);
                for (double &angle : fresh->angles)
                    angle = -angle;
            }
            if (materialized.compare_exchange_strong(current, fresh.get(), std::memory_order_acq_rel,
                                                     std::memory_order_acquire))
                current = fresh.release();
        }
        return *current;
    }

    mutable std::atomic<const StructuredGate *> materialized{nullptr};
};

struct UnitaryPayload : Payload
{
    explicit UnitaryPayload(const UnitaryGateDescriptor &value) : Payload(Kind::Unitary), descriptor(value)
    {
    }

    UnitaryGateDescriptor descriptor;
};

inline std::shared_ptr<const Payload> InternUnitary(const UnitaryGateDescriptor &descriptor)
{
    if (!descriptor)
        return {};

    struct Slot
    {
        UnitaryGateDescriptor descriptor;
        std::shared_ptr<const Payload> payload;
    };

    thread_local std::array<Slot, 16> cache;
    uint64_t hash = static_cast<uint64_t>(descriptor.kind);
    for (double angle : descriptor.angles)
    {
        if (angle == 0.)
            angle = 0.;
        uint64_t bits;
        std::memcpy(&bits, &angle, sizeof(bits));
        hash = (hash ^ bits) * UINT64_C(0x9E3779B97F4A7C15);
    }
    hash = (hash ^ descriptor.controls ^ (uint64_t(descriptor.first) << 8) ^ (uint64_t(descriptor.second) << 16) ^
            (uint64_t(descriptor.adjoint) << 24)) *
           UINT64_C(0x9E3779B97F4A7C15);
    auto &slot = cache[hash >> 60];
    if (slot.payload && slot.descriptor == descriptor)
        return slot.payload;
    auto payload = std::make_shared<const UnitaryPayload>(descriptor);
    slot = {descriptor, payload};
    return payload;
}
} // namespace AppliedGateDetail

// used for recording the applied gates, to be applied again later
// also for 'uncompute'
template <class MatrixClass = Eigen::MatrixXcd> class AppliedGate : public Gates::QuantumGateWithOp<MatrixClass>
{
  public:
    using BaseClass = Gates::QuantumGateWithOp<MatrixClass>;

    AppliedGate() : Gates::QuantumGateWithOp<MatrixClass>(MatrixClass::Zero(1, 1)), q1(0), q2(0), q3(0)
    {
        OperatorChanged();
    }

    AppliedGate(const AppliedGate &other)
        : Gates::QuantumGateWithOp<MatrixClass>(other.operatorMat), q1(other.q1), q2(other.q2), q3(other.q3),
          payload(other.payload)
    {
        structure = other.CurrentStructure();
    }

    AppliedGate(AppliedGate &&other) noexcept
        : Gates::QuantumGateWithOp<MatrixClass>(std::move(other.operatorMat)), q1(other.q1), q2(other.q2), q3(other.q3),
          payload(std::move(other.payload))
    {
        // A derived source may have changed its protected matrix directly.
        // Its matrix has already moved, so classify the destination instead.
        structure = other.IsExactAppliedGate() ? other.structure : ClassifyGateMatrix(BaseClass::operatorMat);
    }

    AppliedGate(const MatrixClass &op, size_t q1 = 0, size_t q2 = 0, size_t q3 = 0)
        : Gates::QuantumGateWithOp<MatrixClass>(op), q1(q1), q2(q2), q3(q3)
    {
        OperatorChanged();
    }

    AppliedGate(MatrixClass &&op, size_t q1 = 0, size_t q2 = 0, size_t q3 = 0) noexcept
        : Gates::QuantumGateWithOp<MatrixClass>(std::move(op)), q1(q1), q2(q2), q3(q3)
    {
        OperatorChanged();
    }

    // Explicit typed construction preserves the constructor candidate. Raw
    // matrix overloads remain untrusted; generic bridges verify every entry
    // before using this candidate, including mutations in derived classes.
    template <class Gate, class = decltype(std::declval<const Gate &>().getRawOperatorMatrix()),
              std::enable_if_t<!std::is_base_of_v<AppliedGate, std::decay_t<Gate>>, int> = 0>
    explicit AppliedGate(const Gate &gate, size_t q1 = 0, size_t q2 = 0, size_t q3 = 0)
        : Gates::QuantumGateWithOp<MatrixClass>(MatrixClass(gate.getRawOperatorMatrix())), q1(q1), q2(q2), q3(q3)
    {
        OperatorChanged();
        const auto candidate = GetUnitaryGateDescriptor(gate);
        payload = AppliedGateDetail::InternUnitary(candidate);
    }

    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return payload && payload->kind == AppliedGateDetail::Payload::Kind::Unitary
                   ? static_cast<const AppliedGateDetail::UnitaryPayload &>(*payload).descriptor
                   : UnitaryGateDescriptor{};
    }

    explicit AppliedGate(StructuredGate operation) : BaseClass(MatrixClass()), q1(0), q2(0), q3(0)
    {
        if (BaseClass::operatorMat.size() != 0)
            throw std::invalid_argument("Compact AppliedGate storage requires a dynamically sized matrix type");
        operation.Validate(std::numeric_limits<size_t>::max());
        payload = std::make_shared<const AppliedGateDetail::StructuredPayload>(std::move(operation));
        RefreshQubits();
        structure = ClassifyGateMatrix(BaseClass::operatorMat);
    }

    static AppliedGate Diagonal(std::vector<std::complex<double>> values, std::vector<size_t> qubits)
    {
        StructuredGate operation;
        operation.diagonal = std::move(values);
        operation.qubits = std::move(qubits);
        return AppliedGate(std::move(operation));
    }

    static AppliedGate MultiControlled(const MatrixClass &matrix, size_t target, std::vector<size_t> controls,
                                       size_t controlState = std::numeric_limits<size_t>::max())
    {
        if (matrix.rows() != 2 || matrix.cols() != 2)
            throw std::invalid_argument("A controlled target matrix must be 2x2");
        StructuredGate operation;
        operation.kind = StructuredGate::Kind::Controlled;
        operation.target = target;
        operation.qubits = std::move(controls);
        operation.controlState = controlState;
        operation.matrix = {matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1)};
        return AppliedGate(std::move(operation));
    }

    static AppliedGate MultiControlled(const BaseClass &gate, size_t target, std::vector<size_t> controls,
                                       size_t controlState = std::numeric_limits<size_t>::max())
    {
        return MultiControlled(gate.getRawOperatorMatrix(), target, std::move(controls), controlState);
    }

    static AppliedGate MultiplexedRotation(RotationAxis axis, std::vector<double> angles, size_t target,
                                           std::vector<size_t> controls)
    {
        StructuredGate operation;
        operation.kind = StructuredGate::Kind::MultiplexedRotation;
        operation.axis = axis;
        operation.angles = std::move(angles);
        operation.target = target;
        operation.qubits = std::move(controls);
        return AppliedGate(std::move(operation));
    }

    // The compact operation as a StructuredGate, null for a dense gate. A gate
    // whose qubits were remapped, or that is the adjoint or the transpose of
    // another compact gate, shares that gate's tables; the first call builds
    // its own StructuredGate (one copy of the table), later calls return it.
    const StructuredGate *getStructuredGate() const override
    {
        return HasCompact() ? &Compact().Gate() : nullptr;
    }

    // The compact operation without building a StructuredGate, null for a
    // dense gate: the shared tables, this gate's qubits and the values
    // derived from the tables.
    const AppliedGateDetail::StructuredPayload *getCompactPayload() const
    {
        return HasCompact() ? &Compact() : nullptr;
    }

    size_t getQubitsNumber() const override
    {
        // A dense gate answers from its matrix alone, as its base class does;
        // only a matrix without rows can stand for a compact operation.
        const auto rows = BaseClass::operatorMat.rows();
        if (rows == 2)
            return 1;
        if (rows == 4)
            return 2;
        if (rows == 8)
            return 3;
        return HasCompact() ? Compact().Arity() : BaseClass::getQubitsNumber();
    }

    // Compact controls follow the target; diagonal qubits retain their table-bit order.
    std::vector<size_t> getQubits() const
    {
        if (HasCompact())
        {
            const auto &compact = Compact();
            std::vector<size_t> result;
            result.reserve(compact.Arity());
            if (compact.OperationKind() != StructuredGate::Kind::Diagonal)
                result.push_back(compact.target);
            result.insert(result.end(), compact.qubits.begin(), compact.qubits.end());
            return result;
        }
        const size_t count = getQubitsNumber();
        if (count > 3)
            throw std::invalid_argument("A legacy full-register matrix does not store individual qubit indices");
        const size_t indices[3] = {q1, q2, q3};
        return std::vector<size_t>(indices, indices + count);
    }

    // The complete qubit list, in the order of getQubits. This is the only way
    // to change a qubit beyond the third one of a compact gate. A compact list
    // must name distinct qubits; its tables are shared, not copied.
    void setQubits(const std::vector<size_t> &qubits)
    {
        if (qubits.size() != getQubitsNumber())
            throw std::invalid_argument("Qubit list must match the applied gate arity");
        if (HasCompact())
        {
            const bool diagonal = Compact().OperationKind() == StructuredGate::Kind::Diagonal;
            const size_t none = std::numeric_limits<size_t>::max();
            for (size_t i = diagonal ? 0 : 1; i < qubits.size(); ++i)
            {
                if (qubits[i] == none)
                    throw std::invalid_argument("Qubit number is too high");
                if (std::find(qubits.begin() + (diagonal ? 0 : 1), qubits.begin() + i, qubits[i]) != qubits.begin() + i)
                    throw std::invalid_argument("Gate qubits must be distinct");
            }
            if (!diagonal &&
                (qubits[0] == none || std::find(qubits.begin() + 1, qubits.end(), qubits[0]) != qubits.end()))
                throw std::invalid_argument("Target must be in range and distinct from controls");
            Rebind(qubits);
        }
        else
        {
            if (qubits.size() > 3)
                throw std::invalid_argument("Legacy dense applied gates store at most three indices");
            q1 = qubits.size() > 0 ? qubits[0] : 0;
            q2 = qubits.size() > 1 ? qubits[1] : 0;
            q3 = qubits.size() > 2 ? qubits[2] : 0;
        }
    }

    AppliedGate Adjoint() const
    {
        if (!payload)
            return AppliedGate(MatrixClass(BaseClass::operatorMat.adjoint()), q1, q2, q3);
        if (HasCompact())
        {
            // The tables stay shared: a controlled matrix is replaced by its
            // adjoint, a diagonal or an angle table is read conjugated.
            const auto &compact = Compact();
            const auto &m = compact.matrix;
            const bool controlled = compact.OperationKind() == StructuredGate::Kind::Controlled;
            return AppliedGate(std::make_shared<const AppliedGateDetail::StructuredPayload>(
                                   compact, compact.qubits, compact.target,
                                   std::array<std::complex<double>, 4>{
                                       {std::conj(m[0]), std::conj(m[2]), std::conj(m[1]), std::conj(m[3])}},
                                   controlled ? false : !compact.conjugated),
                               *this);
        }
        AppliedGate inverse(MatrixClass(BaseClass::operatorMat.adjoint()), q1, q2, q3);
        auto descriptor = getUnitaryDescriptor();
        if (descriptor)
        {
            descriptor.adjoint = !descriptor.adjoint;
            inverse.payload = AppliedGateDetail::InternUnitary(descriptor);
        }
        return inverse;
    }

    // The gate with the transposed matrix, on the same qubits. A compact gate
    // keeps sharing its tables: only a controlled matrix or the sign of the
    // angles of a Y multiplexer changes, every other compact operation is
    // symmetric and is returned as it is.
    AppliedGate Transposed() const
    {
        if (!HasCompact())
            return AppliedGate(MatrixClass(BaseClass::operatorMat.transpose()), q1, q2, q3);
        const auto &compact = Compact();
        const auto &m = compact.matrix;
        if (compact.OperationKind() == StructuredGate::Kind::Controlled)
        {
            if (m[1] == m[2])
                return *this;
            return AppliedGate(std::make_shared<const AppliedGateDetail::StructuredPayload>(
                                   compact, compact.qubits, compact.target,
                                   std::array<std::complex<double>, 4>{{m[0], m[2], m[1], m[3]}}, false),
                               *this);
        }
        if (compact.OperationKind() != StructuredGate::Kind::MultiplexedRotation || compact.Axis() != RotationAxis::Y)
            return *this;
        return AppliedGate(std::make_shared<const AppliedGateDetail::StructuredPayload>(
                               compact, compact.qubits, compact.target, m, !compact.conjugated),
                           *this);
    }

    // The assignments take the source's descriptor instead of classifying the matrix again. A derived
    // source is classified from its current matrix (read before a move); a derived destination still
    // gets its OperatorChanged hook.
    AppliedGate &operator=(const AppliedGate &other)
    {
        if (this == &other)
            return *this;
        BaseClass::operatorMat = other.operatorMat;
        AssignStructure(other.CurrentStructure());
        q1 = other.q1;
        q2 = other.q2;
        q3 = other.q3;
        payload = other.payload;
        return *this;
    }

    AppliedGate &operator=(AppliedGate &&other) noexcept
    {
        if (this == &other)
            return *this;
        const GateStructure otherStructure = other.CurrentStructure();
        BaseClass::operatorMat = std::move(other.operatorMat);
        AssignStructure(otherStructure);
        q1 = other.q1;
        q2 = other.q2;
        q3 = other.q3;
        payload = std::move(other.payload);
        return *this;
    }

    // Kernel metadata is separate from the legacy predicates. Other
    // callers (including MPS routing) rely on AppliedGate's generic flags.
    GateStructure getStructure() const override
    {
        return CurrentStructure();
    }

    size_t getQubit1() const
    {
        return q1;
    }

    size_t getQubit2() const
    {
        return q2;
    }

    size_t getQubit3() const
    {
        return q3;
    }

    void setQubit1(size_t q)
    {
        SetQubit(0, q);
    }

    void setQubit2(size_t q)
    {
        SetQubit(1, q);
    }

    void setQubit3(size_t q)
    {
        SetQubit(2, q);
    }

    bool isBranching() const
    {
        if (HasCompact())
            return getMaxBranching() > 1;
        const auto &mat = BaseClass::getRawOperatorMatrix();
        for (size_t i = 0; i < static_cast<size_t>(mat.rows()); ++i)
        {
            size_t nonZeroElements = 0;
            for (size_t j = 0; j < static_cast<size_t>(mat.cols()); ++j)
            {
                if (std::norm(mat(i, j)) > 1e-15)
                {
                    ++nonZeroElements;
                    if (nonZeroElements > 1)
                        return true;
                }
            }
        }

        return false;
    }

    size_t getMaxBranching() const
    {
        if (HasCompact())
        {
            const auto &compact = Compact();
            const StructuredGate &tables = compact.core->operation;
            // The first entries of a table usually decide. When they do not,
            // the answer comes from the values derived once from the whole
            // table and kept with it, not from another scan.
            constexpr size_t prefix = 16;
            if (tables.kind == StructuredGate::Kind::Diagonal)
            {
                const size_t scanned = std::min(prefix, tables.diagonal.size());
                for (size_t i = 0; i < scanned; ++i)
                    if (tables.diagonal[i] != std::complex<double>(0.))
                        return 1;
                return scanned < tables.diagonal.size() && compact.Info(true)->anyNonzero ? 1 : 0;
            }
            if (tables.kind == StructuredGate::Kind::Controlled)
            {
                const auto &m = compact.matrix;
                // This metadata describes backward branching, as for dense gates:
                // rows of U are columns of the transposed backward operation.
                const size_t first = (m[0] != std::complex<double>(0.)) + (m[1] != std::complex<double>(0.));
                const size_t second = (m[2] != std::complex<double>(0.)) + (m[3] != std::complex<double>(0.));
                return std::max(std::max(first, second), size_t(!compact.qubits.empty()));
            }
            if (tables.axis == RotationAxis::Z)
                return 1;
            const StructuredGateInfo *info = compact.Info(false);
            if (!info)
            {
                const size_t scanned = std::min(prefix, tables.angles.size());
                for (size_t i = 0; i < scanned; ++i)
                {
                    if (tables.angles[i] == 0.)
                        continue;
                    double c, s;
                    detail::RotationHalfAngleCosSin(tables.angles[i], c, s);
                    if (s != 0. && c != 0.)
                        return 2;
                }
                if (scanned == tables.angles.size())
                    return 1;
                info = compact.Info(true);
            }
            return info->anyBothNonzero ? 2 : 1;
        }
        const auto &mat = BaseClass::getRawOperatorMatrix();
        size_t maxBranching = 0;
        for (size_t i = 0; i < static_cast<size_t>(mat.rows()); ++i)
        {
            size_t nonZeroElements = 0;
            for (size_t j = 0; j < static_cast<size_t>(mat.cols()); ++j)
            {
                if (std::norm(mat(i, j)) > 1e-15)
                {
                    ++nonZeroElements;
                }
            }
            if (nonZeroElements > maxBranching)
                maxBranching = nonZeroElements;
        }
        return maxBranching;
    }

    bool isBranchingInState(Eigen::Index state, Eigen::Index &firstNextState) const
    {
        if (HasCompact())
        {
            const auto &compact = Compact();
            const size_t arity = compact.Arity();
            if (state < 0 || (arity < std::numeric_limits<size_t>::digits && size_t(state) >= (size_t{1} << arity)))
                throw std::invalid_argument("Local state is outside the applied gate");
            firstNextState = state;
            if (compact.OperationKind() == StructuredGate::Kind::Diagonal)
                return false;
            const size_t selector = size_t(state) >> 1;
            const size_t column = size_t(state) & 1;
            if (compact.OperationKind() == StructuredGate::Kind::Controlled ? selector != compact.ControlValue()
                                                                            : compact.Axis() == RotationAxis::Z)
                return false;
            // Column `column` of the target block for this selector.
            const auto block = compact.Block(selector, compact.Coefficients(false));
            const std::complex<double> first = block[column], second = block[2 + column];
            if (first != std::complex<double>(0.))
                firstNextState = state & ~Eigen::Index{1};
            else if (second != std::complex<double>(0.))
                firstNextState = state | Eigen::Index{1};
            return first != std::complex<double>(0.) && second != std::complex<double>(0.);
        }
        const auto &mat = BaseClass::getRawOperatorMatrix();
        size_t nonZeroElements = 0;

        assert(state >= 0 && state < mat.rows());

        for (Eigen::Index i = 0; i < mat.rows(); ++i)
        {
            if (std::norm(mat(i, state)) > 1e-15)
            {
                ++nonZeroElements;
                if (nonZeroElements > 1)
                    return true;
                else
                    firstNextState = i;
            }
        }

        return false;
    }

    // The matrix was handed over, so nothing is known of it until it is tested: once, on the first call, and again
    // only after the matrix changes. A caller that knows can say so and spare the test (setUnitary). A compact gate
    // has no matrix here; whoever applies its tables validates them.
    bool isUnitary() const override
    {
        uint8_t known = unitarity.value.load(std::memory_order_relaxed);
        if (known == AppliedGateDetail::KnownUnitarity::NotTested)
        {
            if (HasCompact())
                return true;
            known = IsUnitaryMatrix(BaseClass::operatorMat) ? AppliedGateDetail::KnownUnitarity::Unitary
                                                            : AppliedGateDetail::KnownUnitarity::NotUnitary;
            // Kept by the exact type only, which alone can guarantee that every matrix mutation uses
            // OperatorChanged (see CurrentStructure): a derived gate is tested on every call.
            if (IsExactAppliedGate())
                unitarity.value.store(known, std::memory_order_relaxed);
        }
        return known == AppliedGateDetail::KnownUnitarity::Unitary;
    }

    // For a caller that knows whether its matrix is unitary: the matrix is then not tested, until it changes.
    void setUnitary(bool unitary) noexcept
    {
        unitarity.value.store(unitary ? AppliedGateDetail::KnownUnitarity::Unitary
                                      : AppliedGateDetail::KnownUnitarity::NotUnitary,
                              std::memory_order_relaxed);
    }

  protected:
    void OperatorChanged() override
    {
        payload.reset();
        structure = ClassifyGateMatrix(BaseClass::operatorMat);
        unitarity.value.store(AppliedGateDetail::KnownUnitarity::NotTested, std::memory_order_relaxed);
    }

  private:
    // A gate derived from a compact one: the same tables under the binding of
    // the given payload. The legacy indices beyond the arity are kept.
    AppliedGate(std::shared_ptr<const AppliedGateDetail::StructuredPayload> compact, const AppliedGate &source)
        : BaseClass(MatrixClass()), q1(source.q1), q2(source.q2), q3(source.q3), payload(std::move(compact))
    {
        RefreshQubits();
        structure = ClassifyGateMatrix(BaseClass::operatorMat);
    }

    const AppliedGateDetail::StructuredPayload &Compact() const
    {
        return static_cast<const AppliedGateDetail::StructuredPayload &>(*payload);
    }

    bool HasCompact() const
    {
        // Derived gates may replace the protected matrix directly, or override
        // OperatorChanged. A populated matrix remains authoritative in that case.
        return payload && BaseClass::operatorMat.rows() == 0 &&
               payload->kind == AppliedGateDetail::Payload::Kind::Structured;
    }

    void RefreshQubits()
    {
        const auto &compact = Compact();
        const size_t arity = compact.Arity();
        if (arity > 0)
            q1 = compact.Qubit(0);
        if (arity > 1)
            q2 = compact.Qubit(1);
        if (arity > 2)
            q3 = compact.Qubit(2);
    }

    // The complete qubit list, in the order of getQubits, over the same tables.
    void Rebind(const std::vector<size_t> &qubits)
    {
        const auto &compact = Compact();
        const bool diagonal = compact.OperationKind() == StructuredGate::Kind::Diagonal;
        payload = std::make_shared<const AppliedGateDetail::StructuredPayload>(
            compact, std::vector<size_t>(qubits.begin() + (diagonal ? 0 : 1), qubits.end()),
            diagonal ? compact.target : qubits[0], compact.matrix, compact.conjugated);
        RefreshQubits();
    }

    // The single-index setters. A dense gate stores any of its three indices,
    // used or not. A compact gate does the same: one of its first three
    // qubits is replaced over the shared tables, and an index beyond its arity
    // is kept and ignored. Nothing is checked here, so that a remap may pass
    // through a list with a repeated qubit; a backend rejects the gate only if
    // it is applied in that state. Qubits beyond the third go through setQubits.
    void SetQubit(size_t index, size_t value)
    {
        if (!HasCompact() || index >= Compact().Arity())
        {
            (index == 0 ? q1 : (index == 1 ? q2 : q3)) = value;
            return;
        }
        auto qubits = getQubits();
        qubits[index] = value;
        Rebind(qubits);
    }

    bool IsExactAppliedGate() const
    {
#if defined(__cpp_rtti) || defined(_CPPRTTI)
        return typeid(*this) == typeid(AppliedGate);
#else
        return false;
#endif
    }

    GateStructure CurrentStructure() const
    {
        // Only the exact type can guarantee that every matrix mutation uses
        // OperatorChanged. Subclasses retain access to operatorMat; without
        // RTTI use the same safe, uncached path for all instances.
        return IsExactAppliedGate() ? structure : ClassifyGateMatrix(BaseClass::operatorMat);
    }

    // After the matrix was assigned: the exact type takes the known descriptor, a derived
    // destination (which may override OperatorChanged) is refreshed through its hook.
    void AssignStructure(const GateStructure &known)
    {
        if (IsExactAppliedGate())
        {
            structure = known;
            unitarity.value.store(AppliedGateDetail::KnownUnitarity::NotTested, std::memory_order_relaxed);
        }
        else
            OperatorChanged();
    }

    // don't use it!
    MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0,
                                  size_t controllingQubit2 = 0) const override
    {
        if (HasCompact())
            throw std::invalid_argument("A compact applied gate cannot be expanded through the legacy three-qubit API");
        // this is a hack used only for trying out the old way... without OPTIMIZED_TENSOR_PRODUCT, to not be used
        // otherwise!
        // const size_t N = log2(BaseClass::operatorMat.rows());
        const size_t N =
            BaseClass::getQubitsNumber(); // no need to use log2, there is a function in the base class that computes it
        if (N == 1)
        {
            SingleQubitGate<MatrixClass> gate(BaseClass::operatorMat);
            return gate.getOperatorMatrix(nrQubits, qubit);
        }
        else if (N == 2)
        {
            TwoQubitsGate<MatrixClass> gate(BaseClass::operatorMat);
            return gate.getOperatorMatrix(nrQubits, qubit, controllingQubit1);
        }
        else if (N == 3)
        {
            ThreeQubitsGate<MatrixClass> gate(BaseClass::operatorMat);
            return gate.getOperatorMatrix(nrQubits, qubit, controllingQubit1, controllingQubit2);
        }

        return BaseClass::getRawOperatorMatrix();
    }

    GateStructure structure;
    size_t q1;
    size_t q2;
    size_t q3;
    std::shared_ptr<const AppliedGateDetail::Payload> payload;
    mutable AppliedGateDetail::KnownUnitarity unitarity;
};

template <class MatrixClass = Eigen::MatrixXcd> class HadamardGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::H, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    HadamardGate()
    {
        // 1./sqrt(2) * (X + Z)
        static const double norm = 1. / sqrt(2.);
        OpClass::operatorMat(0, 0) = norm;
        OpClass::operatorMat(0, 1) = norm;
        OpClass::operatorMat(1, 0) = norm;
        OpClass::operatorMat(1, 1) = -norm;
    }
};

// hadamard can be used to switch back and forth to X basis
// this gate is for switching to Y basis and back
// also known as K gate
template <class MatrixClass = Eigen::MatrixXcd> class HyGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Hy, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    HyGate()
    {
        // 1./sqrt(2) * (Y + Z)
        static const double norm = 1. / sqrt(2.);
        OpClass::operatorMat(0, 0) = norm;
        OpClass::operatorMat(0, 1) = std::complex<double>(0, -norm);
        OpClass::operatorMat(1, 0) = std::complex<double>(0, norm);
        OpClass::operatorMat(1, 1) = -norm;
    }
};

// also known as the V gate or Phase gate
template <class MatrixClass = Eigen::MatrixXcd> class SGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::S, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    SGate()
    {
        OpClass::operatorMat(0, 0) = 1;
        OpClass::operatorMat(1, 1) = std::complex<double>(0, 1);
    }

    bool isDiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class SDGGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Sdg, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    SDGGate()
    {
        OpClass::operatorMat(0, 0) = 1;
        OpClass::operatorMat(1, 1) = std::complex<double>(0, -1);
    }

    bool isDiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class TGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::T, {}, 0, {{recordedPhase}});
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    TGate()
    {
        OpClass::operatorMat(0, 0) = 1;
        recordedPhase = std::polar<double>(1, M_PI / 4.);
        OpClass::operatorMat(1, 1) = recordedPhase;
    }

    bool isDiagonal() const override
    {
        return true;
    }

  private:
    std::complex<double> recordedPhase;
};

template <class MatrixClass = Eigen::MatrixXcd> class TDGGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Tdg, {}, 0, {{recordedPhase}});
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    TDGGate()
    {
        OpClass::operatorMat(0, 0) = 1;
        recordedPhase = std::polar<double>(1, -M_PI / 4.);
        OpClass::operatorMat(1, 1) = recordedPhase;
    }

    bool isDiagonal() const override
    {
        return true;
    }

  private:
    std::complex<double> recordedPhase;
};

template <class MatrixClass = Eigen::MatrixXcd> class PhaseShiftGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Phase, {{recordedTheta}}, 0,
                                            {{recordedPhase}});
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    PhaseShiftGate(double theta = 0)
    {
        OpClass::operatorMat(0, 0) = 1;
        SetPhaseShift(theta);
    }

    void SetPhaseShift(double theta)
    {
        recordedTheta = theta;
        recordedPhase = UnitaryPhase(theta);
        OpClass::operatorMat(1, 1) = recordedPhase;
    }

    bool isDiagonal() const override
    {
        return true;
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedPhase = 1.;
};

// also known as the flip gate or not gate
template <class MatrixClass = Eigen::MatrixXcd> class PauliXGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::X, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    PauliXGate()
    {
        OpClass::operatorMat(0, 1) = 1;
        OpClass::operatorMat(1, 0) = 1;
    }

    bool isAntidiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class PauliYGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Y, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    PauliYGate()
    {
        OpClass::operatorMat(0, 1) = std::complex<double>(0, -1);
        OpClass::operatorMat(1, 0) = std::complex<double>(0, 1);
    }

    bool isAntidiagonal() const override
    {
        return true;
    }
};

// also known as the phase flip gate
template <class MatrixClass = Eigen::MatrixXcd> class PauliZGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Z, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    PauliZGate()
    {
        OpClass::operatorMat(0, 0) = 1;
        OpClass::operatorMat(1, 1) = -1;
    }

    bool isDiagonal() const override
    {
        return true;
    }
};

// sx gate, also known as the square root of x gate
template <class MatrixClass = Eigen::MatrixXcd> class SquareRootNOTGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Sx, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    SquareRootNOTGate()
    {
        OpClass::operatorMat(0, 0) = std::complex<double>(0.5, 0.5);
        OpClass::operatorMat(0, 1) = std::complex<double>(0.5, -0.5);
        OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
        OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class SquareRootNOTDagGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Sxdg, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    SquareRootNOTDagGate()
    {
        OpClass::operatorMat(0, 0) = std::complex<double>(0.5, -0.5);
        OpClass::operatorMat(0, 1) = std::complex<double>(0.5, 0.5);
        OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
        OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class SplitterGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Splitter, {}, 0};
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    SplitterGate()
    {
        static const double norm = 1. / sqrt(2.);
        OpClass::operatorMat(0, 0) = norm;
        OpClass::operatorMat(0, 1) = std::complex<double>(0, norm);
        OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
        OpClass::operatorMat(1, 1) = norm;
    }
};

// rotation gates

template <class MatrixClass = Eigen::MatrixXcd> class RotationGate : public SingleQubitGate<MatrixClass>
{
  public:
    using BaseClass = SingleQubitGate<MatrixClass>;

    virtual void SetTheta(double theta) = 0;
};

template <class MatrixClass = Eigen::MatrixXcd> class RxGate : public RotationGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Rx, {{recordedTheta}}, 0, {{recordedCosSin}});
    }

    using BaseClass = RotationGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    RxGate(double theta = 0)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta) override
    {
        recordedTheta = theta;
        double cosine, sine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        recordedCosSin = {cosine, sine};

        OpClass::operatorMat(0, 0) = std::complex<double>(cosine, 0);
        OpClass::operatorMat(0, 1) = std::complex<double>(0, -sine);
        OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
        OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> class RyGate : public RotationGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Ry, {{recordedTheta}}, 0, {{recordedCosSin}});
    }

    using BaseClass = RotationGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    RyGate(double theta = 0)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta) override
    {
        recordedTheta = theta;
        double cosine, sine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        recordedCosSin = {cosine, sine};

        OpClass::operatorMat(0, 0) = std::complex<double>(cosine, 0);
        OpClass::operatorMat(0, 1) = std::complex<double>(-sine, 0);
        OpClass::operatorMat(1, 0) = std::complex<double>(sine, 0);
        OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> class RzGate : public RotationGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Rz, {{recordedTheta}}, 0, {{recordedCosSin}});
    }

    using BaseClass = RotationGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    RzGate(double theta = 0)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta) override
    {
        recordedTheta = theta;
        double cosine, sine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        recordedCosSin = {cosine, sine};

        OpClass::operatorMat(0, 0) = std::complex<double>(cosine, -sine);
        OpClass::operatorMat(1, 1) = std::complex<double>(cosine, sine);
    }

    bool isDiagonal() const override
    {
        return true;
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

// generic rotation gate
template <class MatrixClass = Eigen::MatrixXcd> class UGate : public SingleQubitGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::U, recordedAngles, 0, recordedCoefficients);
    }

    using BaseClass = SingleQubitGate<MatrixClass>;
    using OpClass = typename SingleQubitGate<MatrixClass>::BaseClass;

    UGate(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
    {
        SetParams(theta, phi, lambda, gamma);
    }

    virtual void SetParams(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
    {
        recordedAngles = {theta, phi, lambda, gamma};
        const auto coefficients = UnitaryEulerCoefficients(theta, phi, lambda, gamma);
        recordedCoefficients = coefficients;
        OpClass::operatorMat(0, 0) = coefficients[0];
        OpClass::operatorMat(0, 1) = coefficients[1];
        OpClass::operatorMat(1, 0) = coefficients[2];
        OpClass::operatorMat(1, 1) = coefficients[3];
    }

  private:
    std::array<double, 4> recordedAngles{};
    std::array<std::complex<double>, 4> recordedCoefficients{{1., 0., 0., 1.}};
};

} // namespace Gates
} // namespace QC
