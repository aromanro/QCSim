#pragma once

#include "PauliPropagationEngine.h"
#include <numeric>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

namespace QC
{
// Kept for callers using the public Pauli-string representation.
struct PauliStringHash
{
    size_t operator()(const PauliStringXZWithCoefficient &p) const
    {
        return std::hash<std::vector<bool>>{}(p.X) ^ (std::hash<std::vector<bool>>{}(p.Z) << 1);
    }
};
struct PauliStringEqual
{
    bool operator()(const PauliStringXZWithCoefficient &a, const PauliStringXZWithCoefficient &b) const
    {
        return a.X == b.X && a.Z == b.Z;
    }
};

// Propagate observables backwards through U: <0|U^dagger P U|0>.
// Gates are recorded in statevector order, then visited in reverse order.
// Each kernel applies the adjoint action U^dagger P U, so a recorded S gate
// maps X to -Y (S^dagger X S), whereas SDG maps X to Y. RX/RY/RZ use the
// same convention. Controlled gates store the target in q0 and control in q1.
// Queries own their scratch buffers. Recorded circuits are shared until mutated.
class PauliPropagator
{
    template <class F> auto WithTerm(F &&f) const
    {
        if (nrQubits <= 64)
            return f(PauliDetail::Term<1>(nrQubits));
        if (nrQubits <= 128)
            return f(PauliDetail::Term<2>(nrQubits));
        if (nrQubits <= 256)
            return f(PauliDetail::Term<4>(nrQubits));
        return f(PauliDetail::Term<0>(nrQubits));
    }

  public:
    PauliPropagator() : operations(std::make_shared<Program>())
    {
        std::random_device rd;
        rng.seed(rd());
    }
    void SetSeed(uint64_t seed)
    {
        std::seed_seq seq{uint32_t(seed), uint32_t(seed >> 32)};
        rng.seed(seq);
    }
    int GetNrQubits() const
    {
        return nrQubits;
    }
    void SetNrQubits(int n)
    {
        if (n < 0)
            throw std::invalid_argument("Negative Pauli register width");
        for (const auto &op : *operations)
        {
            for (int q = 0; q < PauliOperationArity(op.type); ++q)
                if (op.Qubit(q) >= n)
                    throw std::invalid_argument("Pauli circuit exceeds register width");
        }
        nrQubits = n;
    }
    double GetCoefficientThreshold() const
    {
        return coefThreshold;
    }
    void SetCoefficientThreshold(double value)
    {
        if (!std::isfinite(value) || value < 0.)
            throw std::invalid_argument("Invalid coefficient threshold");
        coefThreshold = value;
    }
    size_t GetPauliWeightThreshold() const
    {
        return pauliWeightThreshold;
    }
    void SetPauliWeightThreshold(size_t value)
    {
        pauliWeightThreshold = value;
    }
    int StepsBetweenTrims() const
    {
        return stepsBetweenTrims;
    }
    void SetStepsBetweenTrims(int value)
    {
        Positive(value);
        stepsBetweenTrims = value;
    }
    int StepsBetweenDeduplication() const
    {
        return stepsBetweenDeduplication;
    }
    void SetStepsBetweenDeduplication(int value)
    {
        Positive(value);
        stepsBetweenDeduplication = value;
    }
    size_t GetParallelThreshold() const
    {
        return parallelThreshold;
    }
    void SetParallelThreshold(size_t value)
    {
        Positive(value);
        parallelThreshold = value;
    }
    template <class Integer, std::enable_if_t<std::is_integral_v<Integer> && std::is_signed_v<Integer>, int> = 0>
    void SetParallelThreshold(Integer value)
    {
        SetParallelThreshold(PositiveSize(value));
    }
    size_t GetBatchSize() const
    {
        return batchSize;
    }
    void SetBatchSize(size_t value)
    {
        Positive(value);
        batchSize = value;
    }
    template <class Integer, std::enable_if_t<std::is_integral_v<Integer> && std::is_signed_v<Integer>, int> = 0>
    void SetBatchSize(Integer value)
    {
        SetBatchSize(PositiveSize(value));
    }
    size_t GetParallelThresholdForSum() const
    {
        return parallelThresholdSum;
    }
    void SetParallelThresholdForSum(size_t value)
    {
        Positive(value);
        parallelThresholdSum = value;
    }
    template <class Integer, std::enable_if_t<std::is_integral_v<Integer> && std::is_signed_v<Integer>, int> = 0>
    void SetParallelThresholdForSum(Integer value)
    {
        SetParallelThresholdForSum(PositiveSize(value));
    }
    size_t GetBatchSizeForSum() const
    {
        return batchSizeSum;
    }
    void SetBatchSizeForSum(size_t value)
    {
        Positive(value);
        batchSizeSum = value;
    }
    template <class Integer, std::enable_if_t<std::is_integral_v<Integer> && std::is_signed_v<Integer>, int> = 0>
    void SetBatchSizeForSum(Integer value)
    {
        SetBatchSizeForSum(PositiveSize(value));
    }
    size_t GetSamplingCacheMaxNodes() const
    {
        return samplingCacheMaxNodes;
    }
    void SetSamplingCacheMaxNodes(size_t value)
    {
        samplingCacheMaxNodes = value;
    }
    size_t GetSavePosition() const
    {
        return pos;
    }
    void SetSavePosition(size_t value)
    {
        if (value > operations->size())
            throw std::invalid_argument("Invalid Pauli save position");
        pos = value;
    }
    void EnableParallel(size_t workers = 0)
    {
        // The caller participates. Unknown/single-thread hardware stays serial.
        if (!workers)
            workers = PauliDetail::DefaultWorkerCount(std::thread::hardware_concurrency());
        if (!workers)
        {
            DisableParallel();
            return;
        }
        if (!threadPool || threadPool->GetThreadCount() != workers)
            threadPool = std::make_unique<ThreadPool<>>(workers);
    }
    void DisableParallel()
    {
        threadPool.reset();
    }
    bool IsParallelEnabled() const
    {
        return bool(threadPool);
    }
    size_t GetThreadCount() const
    {
        return threadPool ? threadPool->GetThreadCount() : 0;
    }
    void SaveState()
    {
        pos = operations->size();
    }
    void RestoreState()
    {
        pos = std::min(pos, operations->size());
        if (pos < operations->size())
            MutableOperations().resize(pos);
    }
    void ClearOperations()
    {
        operations = std::make_shared<Program>();
        pos = 0;
    }

    void ApplyX(int q)
    {
        Add(OperationType::X, q);
    }
    void ApplyY(int q)
    {
        Add(OperationType::Y, q);
    }
    void ApplyZ(int q)
    {
        Add(OperationType::Z, q);
    }
    void ApplyH(int q)
    {
        Add(OperationType::H, q);
    }
    void ApplyK(int q)
    {
        Add(OperationType::K, q);
    }
    void ApplyS(int q)
    {
        Add(OperationType::S, q);
    }
    void ApplySDG(int q)
    {
        Add(OperationType::SDG, q);
    }
    void ApplySX(int q)
    {
        Add(OperationType::SX, q);
    }
    void ApplySXDG(int q)
    {
        Add(OperationType::SXDG, q);
    }
    void ApplyCX(int control, int target)
    {
        AddTwo(OperationType::CX, target, control);
    }
    void ApplyCY(int control, int target)
    {
        AddTwo(OperationType::CY, target, control);
    }
    void ApplyCZ(int control, int target)
    {
        AddTwo(OperationType::CZ, target, control);
    }
    void ApplySWAP(int a, int b)
    {
        AddTwo(OperationType::SWAP, a, b);
    }
    void ApplyISWAP(int a, int b)
    {
        AddTwo(OperationType::ISWAP, a, b);
    }
    void ApplyISWAPDG(int a, int b)
    {
        AddTwo(OperationType::ISWAPDG, a, b);
    }
    void ApplyRX(int q, double angle)
    {
        AddRotation(OperationType::RX, q, angle);
    }
    void ApplyRY(int q, double angle)
    {
        AddRotation(OperationType::RY, q, angle);
    }
    void ApplyRZ(int q, double angle)
    {
        AddRotation(OperationType::RZ, q, angle);
    }

    // Native gates are single recorded operations. Trimming/deduplication
    // intervals therefore count complete gates, including at special angles.
    // Controlled-gate arguments follow this class's control/target convention.
    void ApplyU(int q, double theta, double phi, double lambda, double gamma = 0.)
    {
        if (!std::isfinite(gamma))
            throw std::invalid_argument("Gate angle must be finite");
        AddLocal(OperationType::U, q, 0, 0, PauliDetail::ParameterizedTransfer(OperationType::U, theta, phi, lambda));
    }
    void ApplyCU(int control, int target, double theta, double phi, double lambda, double gamma = 0.)
    {
        AddLocal(OperationType::CU, target, control, 0,
                 PauliDetail::ParameterizedTransfer(OperationType::CU, theta, phi, lambda, gamma));
    }
    void ApplyCRX(int control, int target, double angle)
    {
        AddLocal(OperationType::CRX, target, control, 0, PauliDetail::ParameterizedTransfer(OperationType::CRX, angle));
    }
    void ApplyCRY(int control, int target, double angle)
    {
        AddLocal(OperationType::CRY, target, control, 0, PauliDetail::ParameterizedTransfer(OperationType::CRY, angle));
    }
    void ApplyCRZ(int control, int target, double angle)
    {
        AddLocal(OperationType::CRZ, target, control, 0, PauliDetail::ParameterizedTransfer(OperationType::CRZ, angle));
    }
    void ApplyCP(int control, int target, double angle)
    {
        AddLocal(OperationType::CP, target, control, 0, PauliDetail::ParameterizedTransfer(OperationType::CP, angle));
    }
    void ApplyCS(int control, int target)
    {
        AddFixed(OperationType::CS, target, control);
    }
    void ApplyCSDAG(int control, int target)
    {
        AddFixed(OperationType::CSDAG, target, control);
    }
    void ApplyCSX(int control, int target)
    {
        AddFixed(OperationType::CSX, target, control);
    }
    void ApplyCSXDAG(int control, int target)
    {
        AddFixed(OperationType::CSXDAG, target, control);
    }
    void ApplyCH(int control, int target)
    {
        AddFixed(OperationType::CH, target, control);
    }
    void ApplyCCX(int control1, int control2, int target)
    {
        AddFixed(OperationType::CCX, target, control1, control2);
    }
    void ApplyCSwap(int control, int target1, int target2)
    {
        AddFixed(OperationType::CSWAP, target1, target2, control);
    }

    std::vector<std::unique_ptr<Operator>> GetOperations() const
    {
        std::vector<std::unique_ptr<Operator>> result;
        result.reserve(operations->size());
        for (const auto &op : *operations)
            result.push_back(op.Legacy());
        return result;
    }
    void SetOperations(std::vector<std::unique_ptr<Operator>> &&input)
    {
        auto program = std::make_shared<Program>();
        program->reserve(input.size());
        for (auto &op : input)
        {
            auto compiled = PauliDetail::Operation::Import(std::move(op));
            for (int q = 0; q < PauliOperationArity(compiled.type); ++q)
            {
                CheckQubit(compiled.Qubit(q));
                for (int other = 0; other < q; ++other)
                    if (compiled.Qubit(q) == compiled.Qubit(other))
                        throw std::invalid_argument("Repeated gate qubit");
            }
            program->push_back(std::move(compiled));
        }
        operations = std::move(program);
        pos = std::min(pos, operations->size());
    }
    // Share the circuit and its checkpoint, including when the source is *this.
    // Subsequent appends/projectors/restore are isolated by copy-on-write.
    void ShareOperationsFrom(const PauliPropagator &source)
    {
        if (nrQubits != source.nrQubits)
            throw std::invalid_argument("Pauli clone width mismatch");
        // Custom operators may contain mutable state; retain their Clone contract.
        const bool custom = std::any_of(source.operations->begin(), source.operations->end(),
                                        [](const PauliDetail::Operation &op) { return bool(op.Custom()); });
        if (custom)
            SetOperations(source.GetOperations());
        else
            operations = source.operations;
        pos = source.pos;
    }

    double ExpectationValue(const std::string &pauli) const
    {
        if (pauli.size() > static_cast<size_t>(nrQubits))
            throw std::invalid_argument("Pauli observable exceeds register width");
        return WithTerm([&](auto term) {
            using T = decltype(term);
            for (size_t q = 0; q < pauli.size(); ++q)
            {
                const char c = pauli[q];
                const uint64_t bit = uint64_t(1) << (q % 64);
                if (c == 'X' || c == 'x' || c == 'Y' || c == 'y')
                    term.X[q / 64] |= bit;
                if (c == 'Z' || c == 'z' || c == 'Y' || c == 'y')
                    term.Z[q / 64] |= bit;
            }
            PauliDetail::Workspace<T> ws;
            ws.terms.push_back(std::move(term));
            return PauliDetail::Execute(ws, *operations, Settings());
        });
    }
    double ExpectationValue(const PauliStringXZWithCoefficient &pauli) const
    {
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            PauliDetail::Workspace<T> ws;
            ws.terms.push_back(PauliDetail::Pack<T>(pauli, nrQubits));
            return PauliDetail::Execute(ws, *operations, Settings());
        });
    }
    double ExpectationValue(PauliStringXZWithCoefficient &&pauli) const
    {
        return ExpectationValue(static_cast<const PauliStringXZWithCoefficient &>(pauli));
    }
    double ExpectationValue(const PauliStringStorage &input) const
    {
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            PauliDetail::Workspace<T> ws;
            ws.terms.reserve(input.size());
            for (const auto &t : input)
                ws.terms.push_back(PauliDetail::Pack<T>(t, nrQubits));
            return PauliDetail::Execute(ws, *operations, Settings());
        });
    }
    // The storage overloads continue to return the propagated expansion in the
    // supplied public buffer. Ordinary scalar queries avoid materializing it.
    double ExpectationValue(const std::string &pauli, PauliStringStorage &output) const
    {
        if (pauli.size() > static_cast<size_t>(nrQubits))
            throw std::invalid_argument("Pauli observable exceeds register width");
        PauliStringXZWithCoefficient p(nrQubits);
        for (size_t q = 0; q < pauli.size(); ++q)
        {
            const char c = pauli[q];
            p.X[q] = c == 'X' || c == 'x' || c == 'Y' || c == 'y';
            p.Z[q] = c == 'Z' || c == 'z' || c == 'Y' || c == 'y';
        }
        return ExpectationValue(std::move(p), output);
    }
    double ExpectationValue(const PauliStringXZWithCoefficient &p, PauliStringStorage &output) const
    {
        auto copy = p;
        return ExpectationValue(std::move(copy), output);
    }
    double ExpectationValue(PauliStringXZWithCoefficient &&p, PauliStringStorage &output) const
    {
        p.Resize(nrQubits);
        output.push_back(std::move(p));
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            PauliDetail::Workspace<T> ws;
            ws.terms.reserve(output.size());
            for (const auto &t : output)
                ws.terms.push_back(PauliDetail::Pack<T>(t, nrQubits));
            const double result = PauliDetail::Execute(ws, *operations, Settings());
            output.clear();
            output.reserve(ws.terms.size());
            for (const auto &t : ws.terms)
                output.push_back(PauliDetail::Unpack(t, nrQubits));
            return result;
        });
    }
    double Probability0(int q) const
    {
        return 0.5 * (1. + BoundedZ(ZExpectation(q)));
    }
    double Probability1(int q) const
    {
        return 0.5 * (1. - BoundedZ(ZExpectation(q)));
    }
    double Probability0(int q, PauliStringStorage &output) const
    {
        CheckQubit(q);
        PauliStringXZWithCoefficient p(nrQubits);
        p.Z[q] = true;
        return 0.5 * (1. + BoundedZ(ExpectationValue(std::move(p), output)));
    }
    double Probability1(int q, PauliStringStorage &output) const
    {
        CheckQubit(q);
        PauliStringXZWithCoefficient p(nrQubits);
        p.Z[q] = true;
        return 0.5 * (1. - BoundedZ(ExpectationValue(std::move(p), output)));
    }
    std::vector<bool> Measure(const std::vector<int> &qubits)
    {
        for (int q : qubits)
            CheckQubit(q);
        std::vector<bool> result;
        result.reserve(qubits.size());
        for (int q : qubits)
        {
            const double p1 = Probability1(q);
            const bool one = uniformZeroOne(rng) < p1;
            result.push_back(one);
            AddProjector(q, one, 0.5 / (one ? p1 : 1. - p1));
        }
        return result;
    }
    double Probability(size_t outcome)
    {
        if (!nrQubits || (nrQubits < int(8 * sizeof(size_t)) && outcome >= (size_t(1) << nrQubits)))
            return 0.;
        const size_t savedSize = operations->size();
        double result = 1.;
        try
        {
            for (int q = 0; q < nrQubits; ++q)
            {
                const bool one = (outcome & 1) != 0;
                outcome >>= 1;
                const double prob = one ? Probability1(q) : Probability0(q);
                if (prob <= 0.)
                {
                    result = 0.;
                    break;
                }
                result *= prob;
                if (q + 1 < nrQubits)
                    AddProjector(q, one, 0.5 / prob);
            }
        }
        catch (...)
        {
            MutableOperations().resize(savedSize);
            throw;
        }
        MutableOperations().resize(savedSize);
        return result;
    }
    std::vector<bool> Sample(const std::vector<int> &qubits)
    {
        CheckSampleQubits(qubits);
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            Sampler<T> sampler;
            return SampleOne(qubits, sampler, 0);
        });
    }
    std::unordered_map<std::vector<bool>, size_t> SampleCounts(const std::vector<int> &qubits, size_t shots)
    {
        CheckSampleQubits(qubits);
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            Sampler<T> sampler;
            const size_t limit = shots > 1 ? samplingCacheMaxNodes : 0;
            if (limit)
                sampler.nodes.emplace_back();
            std::unordered_map<std::vector<bool>, size_t> counts;
            for (size_t shot = 0; shot < shots; ++shot)
                ++counts[SampleOne(qubits, sampler, limit)];
            return counts;
        });
    }

  private:
    using Program = std::vector<PauliDetail::Operation>;
    template <class T> static void Positive(T value)
    {
        if (value <= 0)
            throw std::invalid_argument("Pauli interval/batch must be positive");
    }
    template <class T> static size_t PositiveSize(T value)
    {
        Positive(value);
        if (static_cast<uintmax_t>(value) > std::numeric_limits<size_t>::max())
            throw std::invalid_argument("Pauli interval/batch exceeds size_t");
        return static_cast<size_t>(value);
    }
    // Truncation need not leave a physical distribution. Clamp finite estimates
    // to valid conditional probabilities; reject undefined normalization rather
    // than silently sampling from NaN. This does not bound approximation error.
    static double BoundedZ(double numerator, double denominator = 1.)
    {
        if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator <= 0.)
            throw std::domain_error("Invalid Pauli probability normalization");
        return std::clamp(numerator, -denominator, denominator);
    }
    void CheckQubit(int q) const
    {
        if (q < 0 || q >= nrQubits)
            throw std::out_of_range("Pauli qubit outside register");
    }
    void CheckSampleQubits(const std::vector<int> &qubits) const
    {
        std::unordered_set<int> seen;
        for (int q : qubits)
        {
            CheckQubit(q);
            if (!seen.insert(q).second)
                throw std::invalid_argument("Repeated sampled qubit");
        }
    }
    Program &MutableOperations()
    {
        if (operations.use_count() != 1)
            operations = std::make_shared<Program>(*operations);
        return *operations;
    }
    void Add(OperationType type, int q)
    {
        CheckQubit(q);
        MutableOperations().emplace_back(type, q);
    }
    void AddTwo(OperationType type, int a, int b)
    {
        CheckQubit(a);
        CheckQubit(b);
        if (a == b)
            throw std::invalid_argument("Repeated gate qubit");
        MutableOperations().emplace_back(type, a, b);
    }
    void AddRotation(OperationType type, int q, double angle)
    {
        CheckQubit(q);
        MutableOperations().push_back(PauliDetail::Operation::Rotation(type, q, angle));
    }
    void AddProjector(int q, bool one, double coefficient)
    {
        PauliDetail::Operation op(OperationType::PROJ, q);
        op.projectOne = one;
        op.coefficient = coefficient;
        MutableOperations().push_back(std::move(op));
    }
    void AddLocal(OperationType type, int a, int b, int c, std::shared_ptr<const PauliDetail::LocalTransfer> table)
    {
        const int qubits[] = {a, b, c};
        for (int q = 0; q < PauliOperationArity(type); ++q)
        {
            CheckQubit(qubits[q]);
            for (int other = 0; other < q; ++other)
                if (qubits[q] == qubits[other])
                    throw std::invalid_argument("Repeated gate qubit");
        }
        MutableOperations().push_back(PauliDetail::Operation::Local(type, a, b, c, std::move(table)));
    }
    void AddFixed(OperationType type, int a, int b, int c = 0)
    {
        AddLocal(type, a, b, c, PauliDetail::FixedTransfer(type));
    }
    PauliDetail::Settings Settings() const
    {
        PauliDetail::Settings s;
        s.qubits = nrQubits;
        s.weight = pauliWeightThreshold;
        s.coefficient = coefThreshold;
        s.trims = stepsBetweenTrims;
        s.dedup = stepsBetweenDeduplication;
        s.parallelThreshold = parallelThreshold;
        s.batch = batchSize;
        s.sumThreshold = parallelThresholdSum;
        s.sumBatch = batchSizeSum;
        s.pool = threadPool.get();
        return s;
    }
    double ZExpectation(int q) const
    {
        CheckQubit(q);
        return WithTerm([&](auto term) {
            using T = decltype(term);
            term.Z[q / 64] |= uint64_t(1) << (q % 64);
            PauliDetail::Workspace<T> ws;
            ws.terms.push_back(std::move(term));
            return PauliDetail::Execute(ws, *operations, Settings());
        });
    }
    static constexpr size_t noNode = std::numeric_limits<size_t>::max();
    struct SampleNode
    {
        double numerator = 0.;
        size_t child[2] = {noNode, noNode};
        bool ready = false;
    };
    template <class T> struct Sampler
    {
        PauliDetail::Workspace<T> workspace;
        std::vector<T> prefix, branch;
        std::vector<SampleNode> nodes;
    };
    template <class T> std::vector<bool> SampleOne(const std::vector<int> &qubits, Sampler<T> &sampler, size_t limit)
    {
        auto &start = sampler.prefix;
        auto &branch = sampler.branch;
        auto &ws = sampler.workspace;
        start.clear();
        start.emplace_back(nrQubits);
        std::vector<bool> result;
        result.reserve(qubits.size());
        double denominator = 1.;
        size_t built = 0, node = limit ? 0 : noNode;
        const auto settings = Settings();
        const auto makeBranch = [&](int q) {
            branch.clear();
            branch.reserve(start.size());
            const uint64_t bit = uint64_t(1) << (q % 64);
            for (const auto &t : start)
            {
                branch.push_back(t);
                branch.back().Z[q / 64] |= bit;
            }
        };
        for (size_t depth = 0; depth < qubits.size(); ++depth)
        {
            double numerator;
            if (node != noNode && sampler.nodes[node].ready)
                numerator = sampler.nodes[node].numerator;
            else
            {
                // Materialize the projector expansion only on a cache miss.
                // Replaying the prefix retains the original sampling trim rules.
                while (built < depth)
                {
                    makeBranch(qubits[built]);
                    if (pauliWeightThreshold < static_cast<size_t>(nrQubits) && built % stepsBetweenTrims == 0)
                        PauliDetail::Trim(branch, settings);
                    start.reserve(start.size() + branch.size());
                    for (auto &t : branch)
                    {
                        if (result[built])
                            t.Coefficient = -t.Coefficient;
                        start.push_back(std::move(t));
                    }
                    ++built;
                }
                makeBranch(qubits[depth]);
                ws.terms = branch;
                numerator = PauliDetail::Execute(ws, *operations, settings);
                if (node != noNode)
                {
                    sampler.nodes[node].numerator = numerator;
                    sampler.nodes[node].ready = true;
                }
            }
            // Update the prefix normalization with the same bounded numerator
            // used for the draw, so it cannot drift onto an impossible branch.
            numerator = BoundedZ(numerator, denominator);
            const double p1 = 0.5 * (1. - numerator / denominator);
            const bool one = uniformZeroOne(rng) < p1;
            result.push_back(one);
            if (depth + 1 == qubits.size())
                break;
            denominator += one ? -numerator : numerator;
            if (node != noNode)
            {
                size_t child = sampler.nodes[node].child[one ? 1 : 0];
                if (child == noNode && sampler.nodes.size() < limit)
                {
                    child = sampler.nodes.size();
                    sampler.nodes[node].child[one ? 1 : 0] = child;
                    sampler.nodes.emplace_back();
                }
                node = child;
            }
        }
        return result;
    }

    int nrQubits = 0;
    int stepsBetweenTrims = std::numeric_limits<int>::max();
    int stepsBetweenDeduplication = std::numeric_limits<int>::max();
    std::shared_ptr<Program> operations;
    size_t pos = 0;
    double coefThreshold = 0.;
    size_t pauliWeightThreshold = std::numeric_limits<size_t>::max();
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniformZeroOne{0., 1.};
    std::unique_ptr<ThreadPool<>> threadPool;
    size_t parallelThresholdSum = 65536, batchSizeSum = 16384;
    size_t parallelThreshold = 16384, batchSize = 4096;
    size_t samplingCacheMaxNodes = 8192;
};
} // namespace QC
