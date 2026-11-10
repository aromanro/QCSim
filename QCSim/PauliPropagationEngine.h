#pragma once

#include "PauliPropOper.h"
#include "ThreadPool.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <typeinfo>

namespace QC
{
namespace PauliDetail
{

// The public vector<bool> Pauli API is unchanged. Propagation uses word masks.
template <size_t N> struct Term
{
    std::array<uint64_t, N> X{}, Z{};
    double Coefficient = 1.;

    explicit Term(size_t = 0)
    {
    }
};

template <> struct Term<0>
{
    std::vector<uint64_t> X, Z;
    double Coefficient = 1.;

    explicit Term(size_t qubits = 0) : X((qubits + 63) / 64), Z(X.size())
    {
    }
};

inline size_t Popcount(uint64_t x)
{
#if defined(_MSC_VER) && defined(_M_X64)
    return static_cast<size_t>(__popcnt64(x));
#elif defined(__GNUC__) || defined(__clang__)
    return static_cast<size_t>(__builtin_popcountll(x));
#else
    size_t n = 0;
    for (; x; x &= x - 1)
        ++n;
    return n;
#endif
}

template <class T> size_t Weight(const T &t)
{
    size_t result = 0;
    for (size_t w = 0; w < t.X.size(); ++w)
        result += Popcount(t.X[w] | t.Z[w]);
    return result;
}

template <class T> double Expectation(const T &t)
{
    for (auto x : t.X)
        if (x)
            return 0.;
    return t.Coefficient;
}

template <class T> T Pack(const PauliStringXZWithCoefficient &p, size_t qubits)
{
    T t(qubits);
    t.Coefficient = p.Coefficient;
    for (size_t q = 0; q < std::min(qubits, p.X.size()); ++q)
    {
        if (p.X[q])
            t.X[q / 64] |= uint64_t(1) << (q % 64);
        if (p.Z[q])
            t.Z[q / 64] |= uint64_t(1) << (q % 64);
    }
    return t;
}

template <class T> PauliStringXZWithCoefficient Unpack(const T &t, size_t qubits)
{
    PauliStringXZWithCoefficient p(qubits);
    p.Coefficient = t.Coefficient;
    for (size_t q = 0; q < qubits; ++q)
    {
        p.X[q] = (t.X[q / 64] >> (q % 64)) & 1;
        p.Z[q] = (t.Z[q / 64] >> (q % 64)) & 1;
    }
    return p;
}

// Primitive gates have no per-gate heap allocation. The payload owns an
// immutable local table, the eigenvalues of a two-qubit diagonal channel, or
// a custom operator (userDefined); the operation type tells the first two
// apart. Sharing this slot keeps the primitive descriptor's original size.
struct Operation
{
    OperationType type;
    int q0, q1, q2 = 0;
    // A one-qubit diagonal channel keeps its eigenvalues here instead, by
    // local label: angle for X, sine for Z, cosine for Y.
    double angle = 0., sine = 0., cosine = 1., coefficient = 1.;
    bool projectOne = false;
    bool userDefined = false, clifford = false;
    std::shared_ptr<const void> payload;

    Operation(OperationType t = OperationType::X, int a = 0, int b = 0)
        : type(t), q0(a), q1(b), clifford(IsPauliClifford(t))
    {
    }

    const Operator *Custom() const
    {
        return userDefined ? static_cast<const Operator *>(payload.get()) : nullptr;
    }

    const LocalTransfer &Transfer() const
    {
        return *static_cast<const LocalTransfer *>(payload.get());
    }

    const Eigenvalues &Diagonal() const
    {
        return *static_cast<const Eigenvalues *>(payload.get());
    }

    // Eigenvalues of either diagonal channel, by local label.
    Eigenvalues ChannelEigenvalues() const
    {
        if (type == OperationType::PAULI_CHANNEL_2Q)
            return Diagonal();
        Eigenvalues lambda;
        lambda.fill(1.);
        lambda[1] = angle;
        lambda[2] = sine;
        lambda[3] = cosine;
        return lambda;
    }

    static Operation PauliChannel(int q, const Eigenvalues &lambda)
    {
        Operation op(OperationType::PAULI_CHANNEL_1Q, q);
        op.angle = lambda[1];
        op.sine = lambda[2];
        op.cosine = lambda[3];
        return op;
    }

    static Operation PauliChannel(int a, int b, const Eigenvalues &lambda)
    {
        Operation op(OperationType::PAULI_CHANNEL_2Q, a, b);
        op.payload = std::make_shared<const Eigenvalues>(lambda);
        return op;
    }

    static Operation Reset(int q, bool toOne)
    {
        Operation op(OperationType::RESET, q);
        op.projectOne = toOne;
        return op;
    }

    // The noise properties below derive from the type and its stored values,
    // so the descriptor and the cost of recording a gate are what they were
    // before noise existed. They are meaningful for built-in operations only.

    // A diagonal channel or reset: applied in place, it never adds a term.
    bool InPlace() const
    {
        return IsPauliDiagonalChannel(type) || type == OperationType::RESET;
    }

    // May set coefficients to zero; the engine then compacts the terms.
    bool Annihilates() const
    {
        switch (type)
        {
        case OperationType::PROJ:
        case OperationType::RESET:
            return true;
        case OperationType::PAULI_CHANNEL_1Q:
            return angle == 0. || sine == 0. || cosine == 0.;
        case OperationType::PAULI_CHANNEL_2Q: {
            const auto &lambda = Diagonal();
            return std::any_of(lambda.begin(), lambda.end(), [](double value) { return value == 0.; });
        }
        default:
            return false;
        }
    }

    // Applied in place by Apply, without ever appending to the branch buffer.
    // Consecutive operations of this kind share one pass over the terms. A
    // Clifford gate answers from the two flags it always read.
    bool NonBranching() const
    {
        return !userDefined && (clifford || InPlace());
    }

    int Qubit(int index) const
    {
        return index == 0 ? q0 : index == 1 ? q1 : q2;
    }

    static Operation Local(OperationType type, int a, int b, int c, std::shared_ptr<const LocalTransfer> table)
    {
        Operation op(type, a, b);
        op.q2 = c;
        op.clifford = table->clifford;
        op.payload = std::move(table);
        return op;
    }

    bool Clifford() const
    {
        return !userDefined && clifford;
    }

    static Operation Rotation(OperationType t, int q, double angle)
    {
        Operation op(t, q);
        op.angle = angle;
        // Exact quarter turns never branch and can share Clifford dispatch.
        if (detail::QuarterTurnCosSin(angle, op.cosine, op.sine))
            op.clifford = true;
        return op;
    }

    std::unique_ptr<Operator> Legacy() const
    {
        if (Custom())
            return Custom()->Clone();
        if (IsPauliLocalGate(type))
            return std::make_unique<OperatorLocal>(type, q0, q1, q2,
                                                   std::static_pointer_cast<const LocalTransfer>(payload));
#define QC_PP_OP(name)                                                                                                 \
    case OperationType::name:                                                                                          \
        return std::make_unique<Operator##name>(q0)
#define QC_PP_OP2(name)                                                                                                \
    case OperationType::name:                                                                                          \
        return std::make_unique<Operator##name>(q0, q1)
        switch (type)
        {
            QC_PP_OP(X);
            QC_PP_OP(Y);
            QC_PP_OP(Z);
            QC_PP_OP(H);
            QC_PP_OP(K);
            QC_PP_OP(S);
            QC_PP_OP(SDG);
            QC_PP_OP(SX);
            QC_PP_OP(SXDG);
            QC_PP_OP2(CX);
            QC_PP_OP2(CY);
            QC_PP_OP2(CZ);
            QC_PP_OP2(SWAP);
            QC_PP_OP2(ISWAP);
            QC_PP_OP2(ISWAPDG);
        case OperationType::PROJ:
            return std::make_unique<Projector>(q0, projectOne, coefficient);
        case OperationType::PAULI_CHANNEL_1Q:
            return std::make_unique<OperatorPauliChannel>(q0, ChannelEigenvalues());
        case OperationType::PAULI_CHANNEL_2Q:
            return std::make_unique<OperatorPauliChannel>(q0, q1, ChannelEigenvalues());
        case OperationType::RESET:
            return std::make_unique<OperatorReset>(q0, projectOne);
        case OperationType::RX:
            return std::make_unique<OperatorRX>(q0, angle);
        case OperationType::RY:
            return std::make_unique<OperatorRY>(q0, angle);
        case OperationType::RZ:
            return std::make_unique<OperatorRZ>(q0, angle);
        default:
            break;
        }
#undef QC_PP_OP
#undef QC_PP_OP2
        throw std::invalid_argument("Unknown Pauli operation");
    }

    static Operation Import(std::unique_ptr<Operator> p)
    {
        if (!p)
            throw std::invalid_argument("Null Pauli operation");
        Operation op(p->GetType(), p->GetQubit(0), p->GetNrQubits() > 1 ? p->GetQubit(1) : 0);
        op.q2 = p->GetNrQubits() > 2 ? p->GetQubit(2) : 0;
        if (IsPauliLocalGate(op.type) && typeid(*p) == typeid(OperatorLocal))
            return Local(op.type, op.q0, op.q1, op.q2, static_cast<const OperatorLocal &>(*p).GetTransfer());
        if (IsPauliDiagonalChannel(op.type) && typeid(*p) == typeid(OperatorPauliChannel))
        {
            const auto &lambda = static_cast<const OperatorPauliChannel &>(*p).GetEigenvalues();
            return op.type == OperationType::PAULI_CHANNEL_1Q ? PauliChannel(op.q0, lambda)
                                                              : PauliChannel(op.q0, op.q1, lambda);
        }
        if (op.type == OperationType::RESET && typeid(*p) == typeid(OperatorReset))
            return Reset(op.q0, static_cast<const OperatorReset &>(*p).IsToOne());
        bool builtin = false;
#define QC_PP_MATCH(name)                                                                                              \
    case OperationType::name:                                                                                          \
        builtin = typeid(*p) == typeid(Operator##name);                                                                \
        break
        switch (op.type)
        {
            QC_PP_MATCH(X);
            QC_PP_MATCH(Y);
            QC_PP_MATCH(Z);
            QC_PP_MATCH(H);
            QC_PP_MATCH(K);
            QC_PP_MATCH(S);
            QC_PP_MATCH(SDG);
            QC_PP_MATCH(SX);
            QC_PP_MATCH(SXDG);
            QC_PP_MATCH(CX);
            QC_PP_MATCH(CY);
            QC_PP_MATCH(CZ);
            QC_PP_MATCH(SWAP);
            QC_PP_MATCH(ISWAP);
            QC_PP_MATCH(ISWAPDG);
            QC_PP_MATCH(RX);
            QC_PP_MATCH(RY);
            QC_PP_MATCH(RZ);
        case OperationType::PROJ:
            builtin = typeid(*p) == typeid(Projector);
            break;
        default:
            break;
        }
#undef QC_PP_MATCH
        if (!builtin)
        {
            op.userDefined = true;
            op.payload = std::shared_ptr<const Operator>(std::move(p));
        }
        else if (IsPauliRotation(op.type))
            op = Rotation(op.type, op.q0, static_cast<const OperatorRotation &>(*p).GetAngle());
        else if (op.type == OperationType::PROJ)
        {
            const auto &proj = static_cast<const Projector &>(*p);
            op.projectOne = proj.IsProjectOne();
            op.coefficient = proj.GetCoefficient();
        }
        return op;
    }
};

// Derive signed local Clifford permutations from the established public kernels
// once. This also pins the adjoint and target/control conventions in one place.
inline const std::array<std::array<unsigned char, 16>, 15> &CliffordMaps()
{
    static const auto maps = [] {
        std::array<std::array<unsigned char, 16>, 15> result{};
        for (size_t gate = 0; gate < result.size(); ++gate)
        {
            const auto op = Operation(static_cast<OperationType>(gate), 0, 1).Legacy();
            for (unsigned p = 0; p < 16; ++p)
            {
                PauliStringXZWithCoefficient t(2);
                PauliStringStorage unused;
                t.X[0] = p & 1;
                t.Z[0] = p & 2;
                t.X[1] = p & 4;
                t.Z[1] = p & 8;
                op->Apply(t, unused);
                result[gate][p] = static_cast<unsigned char>((t.X[0] ? 1 : 0) | (t.Z[0] ? 2 : 0) | (t.X[1] ? 4 : 0) |
                                                             (t.Z[1] ? 8 : 0) | (t.Coefficient < 0 ? 16 : 0));
            }
        }
        return result;
    }();
    return maps;
}

template <OperationType Axis, class T>
void Rotate(const Operation &op, std::vector<T> &terms, size_t begin, size_t end, std::vector<T> &extra)
{
    if (op.sine == 0. && op.cosine == 1.)
        return;
    const size_t w = static_cast<size_t>(op.q0) / 64;
    const uint64_t bit = uint64_t(1) << (op.q0 % 64);
    for (size_t i = begin; i < end; ++i)
    {
        auto &t = terms[i];
        if (t.Coefficient == 0.)
            continue;
        const bool x = (t.X[w] & bit) != 0, z = (t.Z[w] & bit) != 0;
        const bool anti = Axis == OperationType::RX ? z : Axis == OperationType::RZ ? x : x != z;
        if (!anti)
            continue;
        if (op.sine == 0.)
        {
            t.Coefficient *= op.cosine;
            continue;
        }
        const double sign = Axis == OperationType::RX   ? (x ? -1. : 1.)
                            : Axis == OperationType::RZ ? (z ? 1. : -1.)
                                                        : (x ? 1. : -1.);
        if (op.cosine == 0.)
        {
            t.Coefficient *= sign * op.sine;
            if constexpr (Axis != OperationType::RZ)
                t.X[w] ^= bit;
            if constexpr (Axis != OperationType::RX)
                t.Z[w] ^= bit;
        }
        else
        {
            const double coefficient = t.Coefficient * sign * op.sine;
            if (coefficient != 0.)
            {
                extra.push_back(t);
                auto &branch = extra.back();
                branch.Coefficient = coefficient;
                if constexpr (Axis != OperationType::RZ)
                    branch.X[w] ^= bit;
                if constexpr (Axis != OperationType::RX)
                    branch.Z[w] ^= bit;
            }
            t.Coefficient *= op.cosine;
        }
    }
}

// Word and bit positions of the K gate qubits, and the local table label of
// a term there, with x,z bits interleaved.
template <unsigned K>
void LocalPositions(const Operation &op, std::array<size_t, K> &words, std::array<unsigned, K> &shifts)
{
    for (unsigned q = 0; q < K; ++q)
    {
        words[q] = static_cast<size_t>(op.Qubit(q)) / 64;
        shifts[q] = unsigned(op.Qubit(q)) % 64;
    }
}

template <unsigned K, class T>
unsigned LocalLabel(const T &t, const std::array<size_t, K> &words, const std::array<unsigned, K> &shifts)
{
    unsigned input = 0;
    for (unsigned q = 0; q < K; ++q)
        input |= unsigned((t.X[words[q]] >> shifts[q]) & 1) << (2 * q) | unsigned((t.Z[words[q]] >> shifts[q]) & 1)
                                                                             << (2 * q + 1);
    return input;
}

template <unsigned K, bool SameWord, class T>
void ApplyLocalKernel(const Operation &op, std::vector<T> &terms, size_t begin, size_t end, std::vector<T> &extra,
                      bool preallocated, size_t write, bool append)
{
    const auto &table = op.Transfer();
    std::array<size_t, K> words;
    std::array<unsigned, K> shifts;
    LocalPositions<K>(op, words, shifts);
    std::array<uint64_t, K> bits;
    for (unsigned q = 0; q < K; ++q)
        bits[q] = uint64_t(1) << shifts[q];
    const auto label = [&](const T &t) { return LocalLabel<K>(t, words, shifts); };
    if (!preallocated)
        write = extra.size();
    if (!preallocated && table.maxOutputs > 1)
    {
        size_t count = extra.size();
        for (size_t i = begin; i < end; ++i)
            if (terms[i].Coefficient != 0.)
            {
                const unsigned input = label(terms[i]);
                const auto outputs = table.offsets[input + 1] - table.offsets[input];
                const size_t branches = outputs > 1 ? outputs - 1 : 0;
                if (branches > extra.max_size() - count)
                    throw std::length_error("Pauli expansion exceeds vector capacity");
                count += branches;
            }
        extra.resize(count);
    }
    std::array<uint64_t, 1U << (2 * K)> xChanges{}, zChanges{};
    const bool useMasks = SameWord && end - begin >= 32;
    if (useMasks)
        for (unsigned delta = 0; delta < xChanges.size(); ++delta)
            for (unsigned q = 0; q < K; ++q)
            {
                xChanges[delta] |= (uint64_t(0) - ((delta >> (2 * q)) & 1U)) & bits[q];
                zChanges[delta] |= (uint64_t(0) - ((delta >> (2 * q + 1)) & 1U)) & bits[q];
            }
    for (size_t i = begin; i < end; ++i)
    {
        auto &t = terms[i];
        if (t.Coefficient == 0.)
            continue;
        const unsigned input = label(t);
        if ((table.unchanged >> input) & 1U)
            continue;
        const auto first = table.offsets[input], last = table.offsets[input + 1];
        const double coefficient = t.Coefficient;
        const uint64_t xWord = t.X[words[0]], zWord = t.Z[words[0]];
        const auto set = [&](T &out, const LocalTransfer::Entry &e) {
            out.Coefficient = coefficient * e.coefficient;
            const unsigned delta = input ^ e.pauli;
            if constexpr (SameWord)
            {
                if (useMasks)
                {
                    out.X[words[0]] = xWord ^ xChanges[delta];
                    out.Z[words[0]] = zWord ^ zChanges[delta];
                    return;
                }
                uint64_t x = 0, z = 0;
                for (unsigned q = 0; q < K; ++q)
                {
                    x |= (uint64_t(0) - ((delta >> (2 * q)) & 1U)) & bits[q];
                    z |= (uint64_t(0) - ((delta >> (2 * q + 1)) & 1U)) & bits[q];
                }
                out.X[words[0]] = xWord ^ x;
                out.Z[words[0]] = zWord ^ z;
            }
            else
            {
                for (unsigned q = 0; q < K; ++q)
                {
                    out.X[words[q]] ^= (uint64_t(0) - ((delta >> (2 * q)) & 1U)) & bits[q];
                    out.Z[words[q]] ^= (uint64_t(0) - ((delta >> (2 * q + 1)) & 1U)) & bits[q];
                }
            }
        };
        for (auto j = first + 1; j < last; ++j)
        {
            if constexpr (std::is_same_v<T, Term<0>> || std::is_same_v<T, Term<4>> || std::is_same_v<T, Term<8>>)
            {
                if (append)
                {
                    extra.push_back(t);
                    set(extra.back(), table.entries[j]);
                    continue;
                }
            }
            auto &out = extra[write++];
            out = t;
            set(out, table.entries[j]);
        }
        if (first == last)
            t.Coefficient = 0.;
        else
            set(t, table.entries[first]);
    }
}

template <class T>
void ApplyLocal(const Operation &op, std::vector<T> &terms, size_t begin, size_t end, std::vector<T> &extra,
                bool preallocated = false, size_t write = 0, bool append = false)
{
    if (op.Transfer().qubits == 1)
        ApplyLocalKernel<1, true>(op, terms, begin, end, extra, preallocated, write, append);
    else if (op.Transfer().qubits == 2)
    {
        if (op.q0 / 64 == op.q1 / 64)
            ApplyLocalKernel<2, true>(op, terms, begin, end, extra, preallocated, write, append);
        else
            ApplyLocalKernel<2, false>(op, terms, begin, end, extra, preallocated, write, append);
    }
    else
    {
        if (op.q0 / 64 == op.q1 / 64 && op.q0 / 64 == op.q2 / 64)
            ApplyLocalKernel<3, true>(op, terms, begin, end, extra, preallocated, write, append);
        else
            ApplyLocalKernel<3, false>(op, terms, begin, end, extra, preallocated, write, append);
    }
}

template <class T>
void Apply(const Operation &op, std::vector<T> &terms, size_t begin, size_t end, std::vector<T> &extra, size_t qubits)
{
    if (op.Custom())
    {
        for (size_t i = begin; i < end; ++i)
        {
            auto p = Unpack(terms[i], qubits);
            PauliStringStorage output;
            op.Custom()->Apply(p, output);
            terms[i] = Pack<T>(p, qubits);
            for (const auto &child : output)
                extra.push_back(Pack<T>(child, qubits));
        }
        return;
    }
    const size_t w = static_cast<size_t>(op.q0) / 64;
    const uint64_t bit = uint64_t(1) << (op.q0 % 64);
    switch (op.type)
    {
    case OperationType::U:
    case OperationType::CU:
    case OperationType::CRX:
    case OperationType::CRY:
    case OperationType::CRZ:
    case OperationType::CP:
    case OperationType::CS:
    case OperationType::CSDAG:
    case OperationType::CSX:
    case OperationType::CSXDAG:
    case OperationType::CH:
    case OperationType::CCX:
    case OperationType::CSWAP:
    case OperationType::CHANNEL_1Q:
    case OperationType::CHANNEL_2Q:
        ApplyLocal(op, terms, begin, end, extra);
        return;
    case OperationType::PAULI_CHANNEL_1Q: {
        // Branch-free: the identity label multiplies by exactly 1.
        const unsigned shift = unsigned(op.q0) % 64;
        const double lambda[4] = {1., op.angle, op.sine, op.cosine};
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            t.Coefficient *= lambda[unsigned((t.X[w] >> shift) & 1) | unsigned((t.Z[w] >> shift) & 1) << 1];
        }
        return;
    }
    case OperationType::PAULI_CHANNEL_2Q: {
        const auto &lambda = op.Diagonal();
        const size_t w1 = static_cast<size_t>(op.q1) / 64;
        const unsigned shift = unsigned(op.q0) % 64, shift1 = unsigned(op.q1) % 64;
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            const unsigned label = unsigned((t.X[w] >> shift) & 1) | unsigned((t.Z[w] >> shift) & 1) << 1 |
                                   unsigned((t.X[w1] >> shift1) & 1) << 2 | unsigned((t.Z[w1] >> shift1) & 1) << 3;
            t.Coefficient *= lambda[label];
        }
        return;
    }
    case OperationType::RESET:
        // The adjoint of an unrecorded reset: X and Y vanish, Z becomes +-I.
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            if (t.X[w] & bit)
                t.Coefficient = 0.;
            else if (t.Z[w] & bit)
            {
                t.Z[w] &= ~bit;
                if (op.projectOne)
                    t.Coefficient = -t.Coefficient;
            }
        }
        return;
    case OperationType::X:
        for (size_t i = begin; i < end; ++i)
            if (terms[i].Z[w] & bit)
                terms[i].Coefficient = -terms[i].Coefficient;
        return;
    case OperationType::Z:
        for (size_t i = begin; i < end; ++i)
            if (terms[i].X[w] & bit)
                terms[i].Coefficient = -terms[i].Coefficient;
        return;
    case OperationType::Y:
        for (size_t i = begin; i < end; ++i)
            if ((terms[i].X[w] ^ terms[i].Z[w]) & bit)
                terms[i].Coefficient = -terms[i].Coefficient;
        return;
    case OperationType::H:
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            if (t.X[w] & t.Z[w] & bit)
                t.Coefficient = -t.Coefficient;
            const uint64_t delta = (t.X[w] ^ t.Z[w]) & bit;
            t.X[w] ^= delta;
            t.Z[w] ^= delta;
        }
        return;
    case OperationType::S:
    case OperationType::SDG:
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            const uint64_t phaseZ = op.type == OperationType::S ? ~t.Z[w] : t.Z[w];
            if (t.X[w] & phaseZ & bit)
                t.Coefficient = -t.Coefficient;
            t.Z[w] ^= t.X[w] & bit;
        }
        return;
    case OperationType::CX: {
        const size_t controlWord = static_cast<size_t>(op.q1) / 64;
        const uint64_t controlBit = uint64_t(1) << (op.q1 % 64);
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            const bool xc = (t.X[controlWord] & controlBit) != 0, zt = (t.Z[w] & bit) != 0;
            const bool xt = (t.X[w] & bit) != 0, zc = (t.Z[controlWord] & controlBit) != 0;
            if (xc && zt && xt == zc)
                t.Coefficient = -t.Coefficient;
            t.X[w] ^= xc ? bit : 0;
            t.Z[controlWord] ^= zt ? controlBit : 0;
        }
        return;
    }
    case OperationType::RX:
        Rotate<OperationType::RX>(op, terms, begin, end, extra);
        return;
    case OperationType::RY:
        Rotate<OperationType::RY>(op, terms, begin, end, extra);
        return;
    case OperationType::RZ:
        Rotate<OperationType::RZ>(op, terms, begin, end, extra);
        return;
    case OperationType::PROJ:
        for (size_t i = begin; i < end; ++i)
        {
            auto &t = terms[i];
            if (t.X[w] & bit)
            {
                t.Coefficient = 0.;
                continue;
            }
            t.Coefficient *= op.coefficient;
            if (t.Coefficient == 0.)
                continue;
            extra.push_back(t);
            auto &branch = extra.back();
            branch.Z[w] ^= bit;
            if (op.projectOne)
                branch.Coefficient = -branch.Coefficient;
        }
        return;
    default:
        break;
    }
    const bool two = PauliOperationArity(op.type) == 2;
    const size_t w1 = two ? static_cast<size_t>(op.q1) / 64 : w;
    const uint64_t bit1 = two ? uint64_t(1) << (op.q1 % 64) : 0;
    const auto &map = CliffordMaps()[static_cast<size_t>(op.type)];
    for (size_t i = begin; i < end; ++i)
    {
        auto &t = terms[i];
        const unsigned input = ((t.X[w] & bit) ? 1 : 0) | ((t.Z[w] & bit) ? 2 : 0) | ((t.X[w1] & bit1) ? 4 : 0) |
                               ((t.Z[w1] & bit1) ? 8 : 0);
        const unsigned out = map[input];
        if (out & 16)
            t.Coefficient = -t.Coefficient;
        t.X[w] = (t.X[w] & ~bit) | ((out & 1) ? bit : 0);
        t.Z[w] = (t.Z[w] & ~bit) | ((out & 2) ? bit : 0);
        if (two)
        {
            t.X[w1] = (t.X[w1] & ~bit1) | ((out & 4) ? bit1 : 0);
            t.Z[w1] = (t.Z[w1] & ~bit1) | ((out & 8) ? bit1 : 0);
        }
    }
}

inline uint64_t Mix(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

template <class T> size_t Hash(const T &t)
{
    uint64_t h = UINT64_C(0x9e3779b97f4a7c15);
    for (size_t w = 0; w < t.X.size(); ++w)
        h = Mix(h ^ Mix(t.X[w]) ^ (Mix(t.Z[w]) + UINT64_C(0x9e3779b97f4a7c15)));
    return static_cast<size_t>(h);
}

template <class T> struct Workspace
{
    std::vector<T> terms;
    std::vector<std::vector<T>> extra;
    std::vector<size_t> hashSlots;
    // Cumulative scratch diagnostic: lets regressions detect repeated growth
    // after logical shrinking without timing-dependent assertions.
    size_t dedupRehashCount = 0;
};

struct Settings
{
    size_t qubits = 0, weight = std::numeric_limits<size_t>::max();
    double coefficient = 0.;
    int trims = std::numeric_limits<int>::max(), dedup = std::numeric_limits<int>::max();
    size_t parallelThreshold = 16384, batch = 4096, sumThreshold = 65536, sumBatch = 16384;
    ThreadPool<> *pool = nullptr;

    bool Dedup(size_t index) const
    {
        return dedup != std::numeric_limits<int>::max() && index % dedup == 0;
    }

    bool Trim(size_t index) const
    {
        return trims != std::numeric_limits<int>::max() && index % trims == 0;
    }
};

template <class T> void Trim(std::vector<T> &terms, const Settings &s)
{
    size_t write = 0;
    for (size_t read = 0; read < terms.size(); ++read)
    {
        auto &t = terms[read];
        if (std::abs(t.Coefficient) > s.coefficient && (s.weight >= s.qubits || Weight(t) <= s.weight))
        {
            if (read != write)
                terms[write] = std::move(t);
            ++write;
        }
    }
    terms.resize(write);
}

template <class T> size_t InitialDedupSlots(const std::vector<T> &terms, bool retainedCapacity = false)
{
    if (terms.size() < 512)
        return 64;
    // Sampling repeated hashes estimates the number of distinct labels without
    // reserving for every input (which can be almost entirely duplicates).
    // This is only a capacity hint; collisions and biased samples cannot affect
    // correctness because the table still grows at its usual load factor.
    constexpr size_t sampleCount = 64;
    std::array<size_t, sampleCount> hashes;
    size_t collisions = 0;
    for (size_t i = 0; i < sampleCount; ++i)
    {
        const size_t index = static_cast<size_t>(Mix(i + UINT64_C(0x9e3779b97f4a7c15))) % terms.size();
        hashes[i] = Hash(terms[index]);
        for (size_t previous = 0; previous < i; ++previous)
            collisions += hashes[previous] == hashes[i];
    }
    const size_t estimate = !collisions || (retainedCapacity && collisions <= 1)
                                ? terms.size()
                                : (sampleCount * (sampleCount - 1) / 2) / collisions;
    const size_t bounded = std::min<size_t>(32768, std::min(terms.size(), std::max<size_t>(32, estimate)));
    size_t slots = 64;
    while (slots / 2 < bounded)
        slots *= 2;
    return slots;
}

template <class T> void Deduplicate(Workspace<T> &ws, const Settings &s)
{
    auto &terms = ws.terms;
    auto &slots = ws.hashSlots;
    constexpr size_t empty = std::numeric_limits<size_t>::max();
    if (slots.empty())
        slots.assign(InitialDedupSlots(terms), empty);
    else
    {
        // Keep physical capacity for a later large query, but avoid clearing
        // its old peak on every small query. The hot insert loop below stays
        // unchanged: it needs no per-term hash cache or occupied-slot list.
        if (terms.size() < slots.size() / 8)
        {
            size_t logical = 64;
            while (logical / 2 < terms.size())
                logical *= 2;
            slots.resize(std::min(slots.size(), logical));
        }
        else if (slots.size() < 1024 && terms.size() / slots.size() > 8)
            slots.resize(std::max(slots.size(), InitialDedupSlots(terms, slots.capacity() > slots.size())));
        else if (slots.capacity() > slots.size() && terms.size() > slots.size() / 2)
        {
            // A preceding small frontier may have shrunk only the logical
            // table. Jump back toward the sampled requirement using retained
            // storage, rather than rehashing through every intermediate size.
            // Steady duplicate tables (size == capacity) do not pay this scan.
            size_t logical = InitialDedupSlots(terms, true);
            while (logical > slots.capacity())
                logical /= 2;
            slots.resize(std::max(slots.size(), logical));
        }
        std::fill(slots.begin(), slots.end(), empty);
    }
    size_t write = 0;
    for (size_t read = 0; read < terms.size(); ++read)
    {
        auto &t = terms[read];
        if (s.weight < s.qubits && Weight(t) > s.weight)
            continue;
        if ((write + 1) * 2 > slots.size())
        {
            slots.assign(slots.size() * 2, empty);
            for (size_t j = 0; j < write; ++j)
            {
                size_t h = Hash(terms[j]) & (slots.size() - 1);
                while (slots[h] != empty)
                    h = (h + 1) & (slots.size() - 1);
                slots[h] = j;
            }
            ++ws.dedupRehashCount;
        }
        size_t h = Hash(t) & (slots.size() - 1);
        while (slots[h] != empty && !(terms[slots[h]].X == t.X && terms[slots[h]].Z == t.Z))
            h = (h + 1) & (slots.size() - 1);
        if (slots[h] != empty)
            terms[slots[h]].Coefficient += t.Coefficient;
        else
        {
            slots[h] = write;
            if (read != write)
                terms[write] = std::move(t);
            ++write;
        }
    }
    terms.resize(write);
    Trim(terms, s);
}

inline size_t Parts(size_t size, size_t threshold, size_t grain, ThreadPool<> *pool)
{
    if (!pool || pool->IsWorkerThread() || size < threshold)
        return 1;
    return std::max<size_t>(1, std::min(pool->GetThreadCount() + 1, size / grain));
}

inline size_t DefaultWorkerCount(size_t hardwareThreads)
{
    return hardwareThreads > 1 ? hardwareThreads - 1 : 0;
}

// Up to four work items per participant balance variable core speeds without
// allocating a branch buffer/future for every term when batch is very small.
inline size_t Grain(size_t size, size_t parts, size_t batch)
{
    return parts == 1 ? std::max<size_t>(1, size) : std::max(batch, 1 + (size - 1) / (parts * 4));
}

inline size_t ChunkCount(size_t size, size_t grain)
{
    return size ? 1 + (size - 1) / grain : 1;
}

// Every submitted job is drained before references to the query's scratch die,
// including allocation/enqueue errors and exceptions thrown by a custom gate.
// Only participants are enqueued; they pick up indexed work items dynamically.
// Output buffers belong to work-item indices, never to completion order.
template <class F> void Chunks(size_t size, size_t parts, ThreadPool<> *pool, F &&fn, size_t grain = 0)
{
    if (parts == 1)
    {
        fn(0, 0, size);
        return;
    }
    std::vector<std::future<double>> jobs;
    jobs.reserve(parts - 1);
    std::exception_ptr error;
    if (!grain)
        grain = 1 + (size - 1) / parts;
    const size_t count = ChunkCount(size, grain);
    std::atomic<size_t> next{0};
    const auto run = [&] {
        std::exception_ptr failure;
        for (;;)
        {
            const size_t job = next.fetch_add(1, std::memory_order_relaxed);
            if (job >= count)
                break;
            const size_t begin = job * grain;
            const size_t end = begin + std::min(grain, size - begin);
            try
            {
                fn(job, begin, end);
            }
            catch (...)
            {
                if (!failure)
                    failure = std::current_exception();
            }
        }
        if (failure)
            std::rethrow_exception(failure);
    };
    try
    {
        for (size_t j = 1; j < parts; ++j)
        {
            jobs.push_back(pool->Enqueue([&] {
                run();
                return 0.;
            }));
        }
        run();
    }
    catch (...)
    {
        error = std::current_exception();
    }
    for (auto &job : jobs)
        try
        {
            job.get();
        }
        catch (...)
        {
            if (!error)
                error = std::current_exception();
        }
    if (error)
        std::rethrow_exception(error);
}

// Native columns tell us the exact output size. Reserve disjoint ranges before
// launching workers, then write branches straight into the final vector. This
// avoids a second branch buffer and a second pass copying the full register.
template <unsigned K, class T>
size_t LocalBranchCount(const std::vector<T> &terms, const Operation &op, size_t begin, size_t end)
{
    const auto &table = op.Transfer();
    std::array<size_t, K> words;
    std::array<unsigned, K> shifts;
    LocalPositions<K>(op, words, shifts);
    size_t count = 0;
    for (size_t i = begin; i < end; ++i)
    {
        const auto &t = terms[i];
        if (t.Coefficient == 0.)
            continue;
        const unsigned input = LocalLabel<K>(t, words, shifts);
        const auto outputs = table.offsets[input + 1] - table.offsets[input];
        count += outputs > 1 ? outputs - 1 : 0;
    }
    return count;
}

template <class T> void ExecuteLocal(std::vector<T> &terms, const Operation &op, const Settings &s)
{
    const auto &table = op.Transfer();
    const size_t size = terms.size();
    const size_t parts = Parts(size, s.parallelThreshold, s.batch, s.pool);
    const size_t grain = parts == 1 ? size : Grain(size, parts, s.batch);
    const size_t chunks = parts == 1 ? 1 : ChunkCount(size, grain);
    std::vector<size_t> offsets(parts == 1 ? 0 : chunks);
    size_t total = size;
    if (table.maxOutputs > 1 && size > (terms.max_size() - size) / (table.maxOutputs - 1))
        throw std::length_error("Pauli expansion exceeds vector capacity");
    const auto countBranches = [&](size_t begin, size_t end) {
        return table.qubits == 1   ? LocalBranchCount<1>(terms, op, begin, end)
               : table.qubits == 2 ? LocalBranchCount<2>(terms, op, begin, end)
                                   : LocalBranchCount<3>(terms, op, begin, end);
    };
    if (parts > 1 && size >= 65536)
    {
        Chunks(
            size, parts, s.pool,
            [&](size_t job, size_t begin, size_t end) { offsets[job] = countBranches(begin, end); }, grain);
        // Counts use worker-independent chunk indices. The prefix sum fixes
        // exactly the same final branch ranges as the sequential count pass.
        for (auto &offset : offsets)
        {
            const size_t count = offset;
            offset = total;
            total += count;
        }
    }
    else
    {
        for (size_t job = 0; job < chunks; ++job)
        {
            if (parts > 1)
                offsets[job] = total;
            const size_t begin = job * grain;
            const size_t end = std::min(size, (job + 1) * grain);
            total += countBranches(begin, end);
        }
    }
    if constexpr (std::is_same_v<T, Term<0>> || std::is_same_v<T, Term<4>> || std::is_same_v<T, Term<8>>)
    {
        if (parts == 1)
        {
            // Copy-construct wide terms directly instead of initializing words
            // that will be overwritten. Reserving first keeps sources stable.
            terms.reserve(total);
            ApplyLocal(op, terms, 0, size, terms, true, size, true);
            return;
        }
    }
    terms.resize(total);
    if (parts == 1)
        ApplyLocal(op, terms, 0, size, terms, true, size);
    else
        Chunks(
            size, parts, s.pool,
            [&](size_t job, size_t begin, size_t end) { ApplyLocal(op, terms, begin, end, terms, true, offsets[job]); },
            grain);
}

template <class T> bool RotationBranches(const Operation &op, const T &term)
{
    if (term.Coefficient == 0. || op.sine == 0. || op.cosine == 0. || term.Coefficient * op.sine == 0.)
        return false;
    const size_t word = static_cast<size_t>(op.q0) / 64;
    const uint64_t bit = uint64_t(1) << (op.q0 % 64);
    const bool x = (term.X[word] & bit) != 0, z = (term.Z[word] & bit) != 0;
    return op.type == OperationType::RX ? z : op.type == OperationType::RZ ? x : x != z;
}

// Compose the adjoint action on the four local Pauli labels. Intermediate
// coefficients are merged on this one qubit before the large frontier is
// traversed. Only exact scalar queries use this: public expansions and cutoff
// paths retain their original per-gate representation and checkpoint behavior.
inline std::shared_ptr<const LocalTransfer> CompileScalarBlock(const std::vector<Operation> &operations, size_t low,
                                                               size_t high)
{
    auto table = std::make_shared<LocalTransfer>(1);
    std::vector<Term<1>> terms, branches;
    terms.reserve(4);
    branches.reserve(4);
    for (unsigned input = 0; input < 4; ++input)
    {
        terms.clear();
        terms.emplace_back();
        terms.back().X[0] = input & 1;
        terms.back().Z[0] = (input >> 1) & 1;
        std::array<double, 4> coefficients{};
        for (size_t next = high + 1; next > low;)
        {
            auto operation = operations[--next];
            operation.q0 = 0;
            branches.clear();
            Apply(operation, terms, 0, terms.size(), branches, 1);
            coefficients.fill(0.);
            for (const auto &term : terms)
                coefficients[unsigned(term.X[0] | (term.Z[0] << 1))] += term.Coefficient;
            for (const auto &term : branches)
                coefficients[unsigned(term.X[0] | (term.Z[0] << 1))] += term.Coefficient;
            terms.clear();
            for (unsigned output = 0; output < 4; ++output)
                if (coefficients[output] != 0.)
                {
                    terms.emplace_back();
                    terms.back().X[0] = output & 1;
                    terms.back().Z[0] = (output >> 1) & 1;
                    terms.back().Coefficient = coefficients[output];
                }
        }
        // A block is nonempty, so the last gate populated all four values.
        table->Column(input, coefficients);
    }
    return table;
}

// A conservative cone is computed once per query, rather than rescanning the
// frontier after each gate. Include zero-coefficient labels: custom operators
// may inspect or revive them. An empty vector means every gate is needed.
template <class T>
std::vector<unsigned char> BackwardCone(const std::vector<T> &terms, const std::vector<Operation> &operations,
                                        size_t qubits)
{
    if (terms.empty() || operations.empty())
        return {};
    auto support = terms.front().X;
    std::fill(support.begin(), support.end(), 0);
    size_t occupied = 0;
    for (const auto &term : terms)
    {
        for (size_t w = 0; w < support.size(); ++w)
        {
            const uint64_t added = (term.X[w] | term.Z[w]) & ~support[w];
            occupied += Popcount(added);
            support[w] |= added;
        }
        if (occupied == qubits)
            return {};
    }
    std::vector<unsigned char> needed;
    for (size_t next = operations.size(); next;)
    {
        const size_t index = --next;
        const auto &op = operations[index];
        // A custom action has no declared support contract. Nothing before it
        // can be excluded from the cone, even if its nominal qubits are local.
        if (op.Custom())
            break;
        bool intersects = op.type == OperationType::PROJ;
        // Public local tables can describe nonunital actions. Only skip their
        // identity input when the table promises the exact identity action.
        if (IsPauliLocalGate(op.type) && !(op.Transfer().unchanged & 1))
            intersects = true;
        const int arity = PauliOperationArity(op.type);
        for (int q = 0; q < arity; ++q)
        {
            const size_t bit = static_cast<size_t>(op.Qubit(q));
            intersects = intersects || ((support[bit / 64] >> (bit % 64)) & 1);
        }
        if (!intersects)
        {
            if (needed.empty())
                needed.assign(operations.size(), 1);
            needed[index] = 0;
            continue;
        }
        for (int q = 0; q < arity; ++q)
        {
            const size_t bit = static_cast<size_t>(op.Qubit(q));
            const uint64_t mask = uint64_t(1) << (bit % 64);
            if (!(support[bit / 64] & mask))
            {
                support[bit / 64] |= mask;
                ++occupied;
            }
        }
        if (occupied == qubits)
            break;
    }
    return needed;
}

enum class CombinationSafety : unsigned char
{
    Unchecked,
    Allowed,
    Forbidden
};

inline bool CanCombine(const std::vector<Operation> &operations, CombinationSafety &safety)
{
    if (safety == CombinationSafety::Unchecked)
    {
        // Contraction and fusion can cancel finite branches before a later
        // nonfinite rotation. Inspect once, only when combining can pay.
        const bool finite = std::all_of(operations.begin(), operations.end(), [](const Operation &op) {
            return !IsPauliRotation(op.type) || (std::isfinite(op.cosine) && std::isfinite(op.sine));
        });
        safety = finite ? CombinationSafety::Allowed : CombinationSafety::Forbidden;
    }
    return safety == CombinationSafety::Allowed;
}

// Once no earlier gate touches these qubits, evaluate their factors against
// |0>. This is for exact scalar queries only: a returned operator expansion or
// coefficient/weight approximation must retain its original representation.
template <class T, class Mask>
void ContractFinished(Workspace<T> &ws, const Mask &mask, const Settings &s, const std::vector<Operation> &operations,
                      CombinationSafety &safety)
{
    auto &terms = ws.terms;
    size_t write = 0;
    // Bound any merged sum, with headroom for rounding. Large finite
    // coefficients must stay separate until earlier rotations scale them.
    const double coefficientLimit =
        (std::numeric_limits<double>::max() / 2.) / static_cast<double>(std::max<size_t>(1, terms.size()));
    bool safeCoefficients = true;
    for (size_t read = 0; read < terms.size(); ++read)
    {
        auto &term = terms[read];
        bool zero = false;
        for (size_t w = 0; w < mask.size(); ++w)
            zero = zero || bool(term.X[w] & mask[w]);
        if (zero)
            continue;
        safeCoefficients = safeCoefficients && std::abs(term.Coefficient) <= coefficientLimit;
        if (read != write)
            terms[write] = std::move(term);
        ++write;
    }
    terms.resize(write);
    // Keep finished Z labels too when combining is unsafe: removing them can
    // make a later scheduled deduplication merge otherwise distinct branches.
    if (!safeCoefficients || (terms.size() > 1 && !CanCombine(operations, safety)))
        return;
    bool changedLabels = false;
    for (auto &term : terms)
        for (size_t w = 0; w < mask.size(); ++w)
        {
            changedLabels = changedLabels || bool(term.Z[w] & mask[w]);
            term.Z[w] &= ~mask[w];
        }
    if (changedLabels && terms.size() > 1)
        Deduplicate(ws, s);
}

template <class T> double SumExpectations(const std::vector<T> &terms, const Settings &s)
{
    const size_t parts = Parts(terms.size(), s.sumThreshold, s.sumBatch, s.pool);
    if (parts == 1)
    {
        double result = 0.;
        for (const auto &t : terms)
            result += Expectation(t);
        return result;
    }
    std::vector<double> sums(parts);
    Chunks(terms.size(), parts, s.pool, [&](size_t job, size_t begin, size_t end) {
        double sum = 0.;
        for (size_t i = begin; i < end; ++i)
            sum += Expectation(terms[i]);
        sums[job] = sum;
    });
    double result = 0.;
    for (double sum : sums)
        result += sum;
    return result;
}

template <class T>
double Execute(Workspace<T> &ws, const std::vector<Operation> &operations, const Settings &s,
               bool preserveExpansion = true)
{
    auto &terms = ws.terms;
    if (terms.empty())
        return 0.;
    if (operations.empty())
        return SumExpectations(terms, s);
    const bool mayContract = !preserveExpansion && s.coefficient == 0. && s.weight >= s.qubits;
    constexpr size_t contractionFrontierThreshold = 64;
    bool attemptedContract = false, contract = false;
    CombinationSafety combinationSafety = CombinationSafety::Unchecked;
    std::vector<size_t> firstTouch;
    std::decay_t<decltype(terms.front().X)> finished{};
    const auto needed = BackwardCone(terms, operations, s.qubits);
    const auto finish = [&](size_t low, size_t high) {
        if (!mayContract)
            return;
        if (!attemptedContract)
        {
            // Planning cannot pay for a tiny nonbranching query. Wait until
            // the frontier grows, while still bounding independent branching.
            if (terms.size() < contractionFrontierThreshold)
                return;
            attemptedContract = true;
            if (std::any_of(operations.begin(), operations.end(),
                            [](const Operation &op) { return bool(op.Custom()); }))
                return;
            contract = true;
            firstTouch.assign(s.qubits, operations.size());
            for (size_t i = 0; i < operations.size(); ++i)
                for (int q = 0; q < PauliOperationArity(operations[i].type); ++q)
                {
                    const size_t bit = static_cast<size_t>(operations[i].Qubit(q));
                    if (firstTouch[bit] == operations.size())
                        firstTouch[bit] = i;
                }
            finished = terms.front().X;
            std::fill(finished.begin(), finished.end(), 0);
            bool any = false;
            // Initialization may happen after many groups. Include every bit
            // already finished, as well as bits the circuit never touches.
            for (size_t q = 0; q < s.qubits; ++q)
                if (firstTouch[q] >= low)
                {
                    finished[q / 64] |= uint64_t(1) << (q % 64);
                    any = true;
                }
            if (any)
                ContractFinished(ws, finished, s, operations, combinationSafety);
            return;
        }
        if (!contract)
            return;
        std::fill(finished.begin(), finished.end(), 0);
        bool any = false;
        for (size_t index = low; index <= high; ++index)
        {
            if (!needed.empty() && !needed[index])
                continue;
            for (int q = 0; q < PauliOperationArity(operations[index].type); ++q)
            {
                const size_t bit = static_cast<size_t>(operations[index].Qubit(q));
                if (firstTouch[bit] == index)
                {
                    finished[bit / 64] |= uint64_t(1) << (bit % 64);
                    any = true;
                }
            }
        }
        if (any)
            ContractFinished(ws, finished, s, operations, combinationSafety);
    };
    if (terms.size() >= contractionFrontierThreshold)
        finish(operations.size(), operations.size());
    size_t firstCustom = operations.size();
    bool checkedCustom = false;
    for (size_t next = operations.size(); next && !terms.empty();)
    {
        const size_t high = next - 1;
        // Compile same-qubit scalar blocks between the original checkpoints.
        // Two or more branching rotations justify one table construction and
        // one frontier traversal; tiny, expansion, approximation and custom
        // queries keep the original per-gate path.
        if (!preserveExpansion && contract && terms.size() >= 16384 && !operations[high].Custom() &&
            !operations[high].Clifford() && IsPauliRotation(operations[high].type) &&
            CanCombine(operations, combinationSafety))
        {
            const int qubit = operations[high].q0;
            bool branches = false;
            for (size_t sample = 0; sample < 32 && !branches; ++sample)
                branches = RotationBranches(operations[high], terms[sample * (terms.size() / 32)]);
            if (branches)
            {
                size_t cursor = next, blockLow = high;
                unsigned rotations = 0;
                while (cursor && high - (cursor - 1) < 32)
                {
                    const size_t index = cursor - 1;
                    const auto &candidate = operations[index];
                    // Diagonal channels and resets between the rotations of one
                    // qubit join the block; only branching rotations justify it.
                    if (candidate.Custom() || candidate.q0 != qubit || PauliOperationArity(candidate.type) != 1 ||
                        IsPauliLocalGate(candidate.type) ||
                        (!candidate.NonBranching() && !IsPauliRotation(candidate.type)))
                        break;
                    rotations += !candidate.NonBranching();
                    blockLow = index;
                    cursor = index;
                    if (s.Dedup(index) || s.Trim(index) || firstTouch[static_cast<size_t>(qubit)] == index)
                        break;
                }
                // Nonfinite coefficients must follow the original arithmetic
                // and zero/NaN removal order: composing first can hide an
                // intermediate Inf-Inf that the unfused path would trim.
                if (rotations >= 2 && std::all_of(terms.begin(), terms.end(),
                                                  [](const T &term) { return std::isfinite(term.Coefficient); }))
                {
                    const auto table = CompileScalarBlock(operations, blockLow, high);
                    const auto fused = Operation::Local(OperationType::U, qubit, 0, 0, table);
                    ExecuteLocal(terms, fused, s);
                    if (s.Dedup(blockLow))
                        Deduplicate(ws, s);
                    else if (s.Trim(blockLow))
                        Trim(terms, s);
                    finish(blockLow, high);
                    next = blockLow;
                    continue;
                }
            }
        }
        size_t low = high;
        const auto &op = operations[high];
        // Clifford gates, diagonal channels and resets never branch: a run of
        // them is one pass over the terms. Terms a run sets to zero are
        // removed after it, as after a projector. The scan reads, for a
        // Clifford gate, only the two flags it read before noise existed.
        const bool run = op.NonBranching();
        bool compact = !op.Custom() && op.Annihilates();
        if (run)
            while (low && !s.Dedup(low) && !s.Trim(low))
            {
                const auto &previous = operations[low - 1];
                if (previous.userDefined)
                    break;
                if (!previous.clifford)
                {
                    if (!previous.InPlace())
                        break;
                    compact = compact || previous.Annihilates();
                }
                --low;
            }
        bool anyNeeded = needed.empty();
        if (!anyNeeded)
            for (size_t index = low; index <= high; ++index)
                anyNeeded = anyNeeded || bool(needed[index]);
        if (!anyNeeded)
        {
            // Keep the original checkpoint indices even when their gates do
            // nothing to the observable. Approximate runs depend on them.
            if (s.Dedup(low))
                Deduplicate(ws, s);
            else if (s.Trim(low))
                Trim(terms, s);
            next = low;
            continue;
        }
        if (!op.Clifford() && !op.Custom() && IsPauliLocalGate(op.type))
        {
            ExecuteLocal(terms, op, s);
            if (s.Dedup(high))
                Deduplicate(ws, s);
            else if (s.Trim(high))
                Trim(terms, s);
            finish(high, high);
            next = high;
            continue;
        }
        const size_t size = terms.size();
        const size_t threshold =
            run ? std::max<size_t>(1, s.parallelThreshold / (high - low + 1)) : s.parallelThreshold;
        const size_t parts = op.Custom() ? 1 : Parts(size, threshold, s.batch, s.pool);
        // Avoid chunk-count division on the common sequential path.
        const size_t grain = parts == 1 ? size : Grain(size, parts, s.batch);
        const size_t chunks = parts == 1 ? 1 : ChunkCount(size, grain);
        if (ws.extra.size() < chunks)
            ws.extra.resize(chunks);
        for (size_t j = 0; j < chunks; ++j)
            ws.extra[j].clear();
        const auto applyChunk = [&](size_t job, size_t begin, size_t end) {
            auto &extra = ws.extra[job];
            if (!run && (op.Custom() || !IsPauliLocalGate(op.type)))
                extra.reserve(end - begin);
            for (size_t pos = high + 1; pos > low;)
            {
                --pos;
                if (needed.empty() || needed[pos])
                    Apply(operations[pos], terms, begin, end, extra, s.qubits);
            }
        };
        // Keep the serial gate loop independent of parallel dispatch and
        // exception-draining machinery, including in compiler inlining decisions.
        if (parts == 1)
            applyChunk(0, 0, size);
        else
            Chunks(size, parts, s.pool, applyChunk, grain);
        if (compact)
        {
            // Compact only after every worker has released its vector range.
            // Custom Apply implementations may observe or revive zero terms;
            // retain those terms while such an operation remains ahead.
            if (!checkedCustom)
            {
                firstCustom =
                    static_cast<size_t>(std::find_if(operations.begin(), operations.end(),
                                                     [](const Operation &gate) { return bool(gate.Custom()); }) -
                                        operations.begin());
                checkedCustom = true;
            }
            if (high < firstCustom)
                terms.erase(
                    std::remove_if(terms.begin(), terms.end(), [](const T &term) { return term.Coefficient == 0.; }),
                    terms.end());
        }
        size_t total = terms.size();
        for (size_t j = 0; j < chunks; ++j)
            total += ws.extra[j].size();
        terms.reserve(total);
        // Original terms followed by branches in input order, regardless of
        // worker count. Global coefficient aggregation precedes truncation.
        for (size_t j = 0; j < chunks; ++j)
            for (auto &child : ws.extra[j])
                terms.push_back(std::move(child));
        if (s.Dedup(low))
            Deduplicate(ws, s);
        else if (s.Trim(low))
            Trim(terms, s);
        finish(low, high);
        next = low;
    }
    return SumExpectations(terms, s);
}

} // namespace PauliDetail
} // namespace QC
