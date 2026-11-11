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
        else if (nrQubits <= 128)
            return f(PauliDetail::Term<2>(nrQubits));
        else if (nrQubits <= 256)
            return f(PauliDetail::Term<4>(nrQubits));
        else if (nrQubits <= 512)
            return f(PauliDetail::Term<8>(nrQubits));
        return f(PauliDetail::Term<0>(nrQubits));
    }

  public:
    PauliPropagator() : operations(std::make_shared<Program>())
    {
        std::random_device rdv;
        rng.seed(rdv());
    }

    // Seeds the engine with the value itself. A std::seed_seq costs about ten
    // microseconds per call, too much for callers that reseed for every shot.
    void SetSeed(uint64_t seed)
    {
        rng.seed(seed);
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
        if (nrQubits != n)
            InvalidateSampling();
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
        if (coefThreshold != value)
            InvalidateSampling();
        coefThreshold = value;
    }

    size_t GetPauliWeightThreshold() const
    {
        return pauliWeightThreshold;
    }

    void SetPauliWeightThreshold(size_t value)
    {
        if (pauliWeightThreshold != value)
            InvalidateSampling();
        pauliWeightThreshold = value;
    }

    int StepsBetweenTrims() const
    {
        return stepsBetweenTrims;
    }

    void SetStepsBetweenTrims(int value)
    {
        Positive(value);
        if (stepsBetweenTrims != value)
            InvalidateSampling();
        stepsBetweenTrims = value;
    }

    int StepsBetweenDeduplication() const
    {
        return stepsBetweenDeduplication;
    }

    void SetStepsBetweenDeduplication(int value)
    {
        Positive(value);
        if (stepsBetweenDeduplication != value)
            InvalidateSampling();
        stepsBetweenDeduplication = value;
    }

    size_t GetParallelThreshold() const
    {
        return parallelThreshold;
    }

    void SetParallelThreshold(size_t value)
    {
        Positive(value);
        if (parallelThreshold != value)
            InvalidateSampling();
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
        if (batchSize != value)
            InvalidateSampling();
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
        if (parallelThresholdSum != value)
            InvalidateSampling();
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
        if (batchSizeSum != value)
            InvalidateSampling();
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
        if (samplingCacheMaxNodes != value)
            InvalidateSampling();
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
        {
            InvalidateSampling();
            threadPool = std::make_unique<ThreadPool<>>(workers);
        }
    }

    void DisableParallel()
    {
        if (threadPool)
            InvalidateSampling();
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
        importedSamplingTables = false;
        InvalidateSampling();
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

    // ---- Noise channels ----
    //
    // Each call records one operation. The observable is propagated through
    // the adjoint (Heisenberg) channel, so the recorded circuit gives the exact
    // ensemble average of the noisy state: no trajectories are sampled, and
    // expectation values, probabilities and samples all describe the mixed
    // state. Pauli channels and phase damping are diagonal in the Pauli basis
    // and only scale terms. Damping channels send Z to a mix of I and Z and can
    // add a term. A reset never adds one.

    // (1 - px - py - pz) rho + px X rho X + py Y rho Y + pz Z rho Z
    void ApplyPauliChannel(int q, double px, double py, double pz)
    {
        const double identity = 1. - (px + py + pz);
        if (!(identity >= -1E-9))
            throw std::invalid_argument("Pauli error probabilities must sum to at most one");
        const std::array<double, 4> probabilities{{std::max(0., identity), px, py, pz}};
        AddChannel(Diagonal(PauliDetail::PauliChannelEigenvalues(probabilities, 1)), q, 0, 1);
    }

    // probabilities[p0 + 4 * p1] is the probability of the error p0 on q0 and p1
    // on q1, each indexed I, X, Y, Z = 0..3. The sixteen values sum to one.
    void ApplyTwoQubitPauliChannel(int q0, int q1, const std::array<double, 16> &probabilities)
    {
        AddChannel(Diagonal(PauliDetail::PauliChannelEigenvalues(probabilities, 2)), q0, q1, 2);
    }

    void ApplyBitFlipNoise(int q, double p)
    {
        ApplyPauliChannel(q, p, 0., 0.);
    }

    void ApplyPhaseFlipNoise(int q, double p)
    {
        ApplyPauliChannel(q, 0., 0., p);
    }

    // (1 - p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
    void ApplyDepolarizingNoise(int q, double p)
    {
        ApplyPauliChannel(q, p / 3., p / 3., p / 3.);
    }

    // (1 - p) rho + p/15 (sum over the fifteen non-identity two-qubit Paulis)
    void ApplyTwoQubitDepolarizingNoise(int q0, int q1, double p)
    {
        CheckProbability(p, "Depolarizing probability");
        std::array<double, 16> probabilities;
        probabilities.fill(p / 15.);
        probabilities[0] = 1. - p;
        ApplyTwoQubitPauliChannel(q0, q1, probabilities);
    }

    // Coherences are multiplied by sqrt(1 - gamma); populations are unchanged.
    void ApplyPhaseDamping(int q, double gamma)
    {
        CheckProbability(gamma, "Phase-damping probability");
        const double coherence = std::sqrt(1. - gamma);
        AddChannel(PauliDetail::CompileRelaxation(coherence, 1., 0.), q, 0, 1);
    }

    // Relaxation with probability gamma towards the population
    // excitedPopulation of |1> (zero: plain |1> -> |0> amplitude damping).
    void ApplyAmplitudeDamping(int q, double gamma, double excitedPopulation = 0.)
    {
        CheckProbability(gamma, "Amplitude-damping probability");
        CheckProbability(excitedPopulation, "Excited-state population");
        AddChannel(
            PauliDetail::CompileRelaxation(std::sqrt(1. - gamma), 1. - gamma, gamma * (1. - 2. * excitedPopulation)), q,
            0, 1);
    }

    // Populations relax towards excitedPopulation as exp(-duration / t1) and
    // coherences decay as exp(-duration / t2). Requires t2 <= 2 * t1; either
    // time may be infinite. No clamping: t1 < t2 <= 2 * t1 is represented exactly.
    void ApplyThermalRelaxation(int q, double duration, double t1, double t2, double excitedPopulation = 0.)
    {
        if (!std::isfinite(duration) || duration < 0.)
            throw std::invalid_argument("Thermal-relaxation duration must be finite and nonnegative");
        if (!(t1 > 0.) || !(t2 > 0.))
            throw std::invalid_argument("T1 and T2 must be positive");
        CheckProbability(excitedPopulation, "Excited-state population");
        if (std::isinf(t2) ? !std::isinf(t1) : (std::isfinite(t1) && t2 > 2. * t1 * (1. + 1E-12)))
            throw std::invalid_argument("Thermal relaxation requires T2 <= 2*T1");
        const double coherence = std::isinf(t2) ? 1. : std::exp(-duration / t2);
        const double gamma = std::isinf(t1) ? 0. : -std::expm1(-duration / t1);
        AddChannel(PauliDetail::CompileRelaxation(coherence, 1. - gamma, gamma * (1. - 2. * excitedPopulation)), q, 0,
                   1);
    }

    // An arbitrary one- or two-qubit channel from its Kraus operators (2x2 or
    // 4x4 matrices with operator()(row, col), rows() and cols(); for two qubits
    // q0 is the low bit of the matrix index). The operators must be trace
    // preserving. A channel whose transfer matrix turns out diagonal is
    // recorded as a Pauli channel and does not branch.
    template <class Matrix> void ApplyKrausChannel(const std::vector<Matrix> &kraus, int q0, int q1 = 0)
    {
        if (kraus.empty())
            throw std::invalid_argument("A channel must contain at least one Kraus operator");
        const auto dimension = static_cast<size_t>(kraus.front().rows());
        if (dimension != 2 && dimension != 4)
            throw std::invalid_argument("Kraus operators must act on one or two qubits");
        const unsigned qubits = dimension == 2 ? 1 : 2;
        AddChannel(PauliDetail::CompileChannel(PauliDetail::KrausTransfer(kraus, qubits), qubits), q0, q1, qubits);
    }

    // Reset to |0> (or |1>) without recording an outcome. Unlike Measure, this
    // is deterministic: it draws no random number and inserts no projector.
    void ApplyReset(int q, bool toOne = false)
    {
        CheckQubit(q);
        MutableOperations().push_back(PauliDetail::Operation::Reset(q, toOne));
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
        bool externalTables = false;
        program->reserve(input.size());
        for (auto &op : input)
        {
            auto compiled = PauliDetail::Operation::Import(std::move(op));
            externalTables = externalTables || IsPauliLocalGate(compiled.type);
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
        importedSamplingTables = externalTables;
        InvalidateSampling();
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
        {
            operations = source.operations;
            importedSamplingTables = source.importedSamplingTables;
            InvalidateSampling();
        }
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
            WorkspaceLease<T> lease;
            auto &ws = lease.Get();
            ws.terms.push_back(std::move(term));
            return PauliDetail::Execute(ws, *operations, Settings(), false);
        });
    }

    double ExpectationValue(const PauliStringXZWithCoefficient &pauli) const
    {
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            WorkspaceLease<T> lease;
            auto &ws = lease.Get();
            ws.terms.push_back(PauliDetail::Pack<T>(pauli, nrQubits));
            return PauliDetail::Execute(ws, *operations, Settings(), false);
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
            WorkspaceLease<T> lease;
            auto &ws = lease.Get();
            ws.terms.reserve(input.size());
            for (const auto &t : input)
                ws.terms.push_back(PauliDetail::Pack<T>(t, nrQubits));
            return PauliDetail::Execute(ws, *operations, Settings(), false);
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
            WorkspaceLease<T> lease;
            auto &ws = lease.Get();
            ws.terms.reserve(output.size());
            for (const auto &t : output)
                ws.terms.push_back(PauliDetail::Pack<T>(t, nrQubits));
            const double result = PauliDetail::Execute(ws, *operations, Settings(), true);
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
        if (nrQubits <= maximumMomentQubits && CanUseMoments(static_cast<size_t>(nrQubits)))
        {
            std::vector<int> qubits(static_cast<size_t>(nrQubits));
            std::iota(qubits.begin(), qubits.end(), 0);
            return WithTerm([&](auto tag) {
                using T = decltype(tag);
                Sampler<T> local;
                auto &sampler = SamplingContext(qubits, local, 0);
                SamplingScratchGuard<T> guard(sampler);
                if (sampler.moments.empty())
                {
                    for (size_t index = 0; index < sampler.distinctOutcomeCount; ++index)
                        if (sampler.requestedOutcomes[index] == outcome)
                            return sampler.requestedProbabilities[index];
                    // Preserve the efficient normalized-projector algorithm
                    // for isolated outcomes. Cache its scalar answer, then
                    // promote a genuine multi-outcome workload to full moments.
                    if (sampler.distinctOutcomeCount + 1 < sampler.requestedOutcomes.size())
                    {
                        const double probability = UncachedProbabilityPreservingContext(outcome);
                        const size_t index = sampler.distinctOutcomeCount++;
                        sampler.requestedOutcomes[index] = outcome;
                        sampler.requestedProbabilities[index] = probability;
                        return probability;
                    }
                    PrepareMoments(sampler, qubits.size());
                }
                if (sampler.probabilities.empty())
                {
                    for (size_t subset = 0; subset < sampler.moments.size(); ++subset)
                        Moment(qubits, sampler, subset);
                    auto probabilities = sampler.moments;
                    // A Walsh-Hadamard transform produces every outcome from
                    // the same Z-subset moments in O(n * 2^n), once per circuit.
                    for (size_t stride = 1; stride < probabilities.size(); stride *= 2)
                        for (size_t block = 0; block < probabilities.size(); block += 2 * stride)
                            for (size_t index = block; index < block + stride; ++index)
                            {
                                const double a = probabilities[index], b = probabilities[index + stride];
                                probabilities[index] = a + b;
                                probabilities[index + stride] = a - b;
                            }
                    for (auto &probability : probabilities)
                    {
                        probability /= probabilities.size();
                        if (!std::isfinite(probability))
                            throw std::domain_error("Invalid Pauli probability normalization");
                        probability = std::clamp(probability, 0., 1.);
                    }
                    sampler.probabilities = std::move(probabilities);
                }
                return sampler.probabilities[outcome];
            });
        }
        return UncachedProbability(outcome);
    }

    std::vector<bool> Sample(const std::vector<int> &qubits)
    {
        CheckSampleQubits(qubits);
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            Sampler<T> local;
            auto &sampler = SamplingContext(qubits, local, 0);
            SamplingScratchGuard<T> guard(sampler);
            return SampleOne(qubits, sampler, sampler.limit);
        });
    }

    std::unordered_map<std::vector<bool>, size_t> SampleCounts(const std::vector<int> &qubits, size_t shots)
    {
        CheckSampleQubits(qubits);
        if (!shots)
            return {};
        return WithTerm([&](auto tag) {
            using T = decltype(tag);
            Sampler<T> local;
            auto &sampler = SamplingContext(qubits, local, shots > 1 ? EffectiveSamplingLimit() : 0, shots >= 16);
            SamplingScratchGuard<T> guard(sampler);
            std::unordered_map<std::vector<bool>, size_t> counts;
            for (size_t shot = 0; shot < shots; ++shot)
                ++counts[SampleOne(qubits, sampler, sampler.limit)];
            return counts;
        });
    }

  private:
    using Program = std::vector<PauliDetail::Operation>;

    // A thread owns one reusable workspace per packed width. A custom operator
    // can reenter a query on this thread; that nested query gets local scratch.
    // Completed jobs are drained by Execute before the lease is returned.
    template <class T> class WorkspaceLease
    {
        struct Cached
        {
            PauliDetail::Workspace<T> workspace;
            bool busy = false;
        };

        static Cached &Cache()
        {
            static thread_local Cached value;
            return value;
        }

      public:
        WorkspaceLease()
        {
            auto &candidate = Cache();
            if (!candidate.busy)
            {
                candidate.busy = true;
                cached = &candidate;
                value = &candidate.workspace;
            }
            else
                value = &local;
        }

        WorkspaceLease(const WorkspaceLease &) = delete;
        WorkspaceLease &operator=(const WorkspaceLease &) = delete;

        ~WorkspaceLease()
        {
            if (!cached)
                return;
            auto &ws = cached->workspace;
            // Bound per-thread retention after a large one-off query. Clearing
            // terms destroys heap masks of dynamically packed wide labels.
            constexpr size_t maximumBytes = 16 * 1024 * 1024;
            size_t bytes = ws.terms.capacity() * sizeof(T) + ws.extra.capacity() * sizeof(std::vector<T>) +
                           ws.hashSlots.capacity() * sizeof(size_t);
            for (const auto &extra : ws.extra)
                bytes += extra.capacity() * sizeof(T);
            if (bytes > maximumBytes)
                ws = PauliDetail::Workspace<T>{};
            else
            {
                ws.terms.clear();
                for (auto &extra : ws.extra)
                    extra.clear();
            }
            cached->busy = false;
        }

        PauliDetail::Workspace<T> &Get()
        {
            return *value;
        }

      private:
        PauliDetail::Workspace<T> local;
        Cached *cached = nullptr;
        PauliDetail::Workspace<T> *value = nullptr;
    };

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
        // Up to eight readout sites avoid hash allocations. Larger subsets
        // keep the original path: measured fresh-instance costs favor it.
        // Preserve left-to-right range/duplicate exception order in both paths.
        if (qubits.size() <= 8)
        {
            for (size_t index = 0; index < qubits.size(); ++index)
            {
                CheckQubit(qubits[index]);
                for (size_t previous = 0; previous < index; ++previous)
                    if (qubits[previous] == qubits[index])
                        throw std::invalid_argument("Repeated sampled qubit");
            }
            return;
        }
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
        InvalidateSampling();
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

    static void CheckProbability(double value, const char *name)
    {
        if (!(value >= 0. && value <= 1.))
            throw std::invalid_argument(std::string(name) + " must lie in [0, 1]");
    }

    static PauliDetail::CompiledChannel Diagonal(const PauliDetail::Eigenvalues &eigenvalues)
    {
        PauliDetail::CompiledChannel channel;
        channel.eigenvalues = eigenvalues;
        return channel;
    }

    void AddChannel(const PauliDetail::CompiledChannel &channel, int a, int b, unsigned qubits)
    {
        CheckQubit(a);
        if (qubits == 2)
        {
            CheckQubit(b);
            if (a == b)
                throw std::invalid_argument("Repeated gate qubit");
        }
        if (!channel.diagonal)
            MutableOperations().push_back(
                PauliDetail::Operation::Local(qubits == 1 ? OperationType::CHANNEL_1Q : OperationType::CHANNEL_2Q, a,
                                              qubits == 2 ? b : 0, 0, channel.table));
        else if (qubits == 1)
            MutableOperations().push_back(PauliDetail::Operation::PauliChannel(a, channel.eigenvalues));
        else
            MutableOperations().push_back(PauliDetail::Operation::PauliChannel(a, b, channel.eigenvalues));
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
            WorkspaceLease<T> lease;
            auto &ws = lease.Get();
            ws.terms.push_back(std::move(term));
            return PauliDetail::Execute(ws, *operations, Settings(), false);
        });
    }

    static constexpr size_t noNode = std::numeric_limits<size_t>::max();

    struct SampleNode
    {
        double numerator = 0.;
        size_t child[2] = {noNode, noNode};
        bool ready = false;
    };

    struct SamplingCacheBase
    {
        virtual ~SamplingCacheBase() = default;
        std::weak_ptr<const Program> program;
        std::vector<int> qubits;
        uint64_t revision = 0;
        int width = 0;
        size_t limit = 0;
    };

    // Temporary projectors must not destroy the persistent context being
    // filled by Probability(). Holding the original program also isolates all
    // temporary mutations through copy-on-write, including exceptional exits.
    struct ScopedProbabilityState
    {
        PauliPropagator &owner;
        std::shared_ptr<Program> program;
        std::unique_ptr<SamplingCacheBase> cache;
        uint64_t revision;
        int eligibility;
        bool externalTables;
        size_t position;

        explicit ScopedProbabilityState(PauliPropagator &p)
            : owner(p), program(p.operations), cache(std::move(p.samplingCache)), revision(p.samplingRevision),
              eligibility(p.samplingEligibility), externalTables(p.importedSamplingTables), position(p.pos)
        {
        }

        ~ScopedProbabilityState()
        {
            owner.operations = std::move(program);
            owner.samplingCache = std::move(cache);
            owner.samplingRevision = revision;
            owner.samplingEligibility = eligibility;
            owner.importedSamplingTables = externalTables;
            owner.pos = position;
        }
    };

    double UncachedProbability(size_t outcome)
    {
        const size_t savedSize = operations->size();
        double result = 1.;
        try
        {
            for (int q = 0; q < nrQubits; ++q)
            {
                const bool one = (outcome & 1) != 0;
                outcome >>= 1;
                const double probability = one ? Probability1(q) : Probability0(q);
                if (probability <= 0.)
                {
                    result = 0.;
                    break;
                }
                result *= probability;
                if (q + 1 < nrQubits)
                    AddProjector(q, one, 0.5 / probability);
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

    double UncachedProbabilityPreservingContext(size_t outcome)
    {
        ScopedProbabilityState state(*this);
        return UncachedProbability(outcome);
    }

    template <class T> struct Sampler : SamplingCacheBase
    {
        PauliDetail::Workspace<T> workspace;
        std::vector<T> prefix, branch;
        std::vector<SampleNode> nodes;
        std::vector<double> moments;
        std::vector<unsigned char> momentReady;
        std::vector<double> probabilities;
        std::array<size_t, 16> requestedOutcomes{};
        std::array<double, 16> requestedProbabilities{};
        size_t distinctOutcomeCount = 0;
    };

    static constexpr int maximumMomentQubits = 13;

    void InvalidateSampling()
    {
        ++samplingRevision;
        samplingEligibility = -1;
        samplingCache.reset();
    }

    size_t EffectiveSamplingLimit() const
    {
        // A cache is an acceleration only. Bound retention even when a caller
        // provides an excessive node count; the uncached path remains valid.
        return std::min(samplingCacheMaxNodes, size_t{1} << 18);
    }

    int SamplingEligibility()
    {
        if (samplingEligibility >= 0)
            return samplingEligibility;
        samplingEligibility = importedSamplingTables ? 0 : 2;
        for (const auto &op : *operations)
        {
            if (op.Custom())
                return samplingEligibility = 0;
            if (op.type == OperationType::PROJ && samplingEligibility)
                samplingEligibility = 1;
        }
        return samplingEligibility;
    }

    bool CanUseMoments(size_t count)
    {
        return count <= static_cast<size_t>(maximumMomentQubits) && coefThreshold == 0. &&
               pauliWeightThreshold >= static_cast<size_t>(nrQubits) &&
               (size_t{1} << count) <= EffectiveSamplingLimit() && SamplingEligibility() == 2;
    }

    template <class T>
    Sampler<T> &SamplingContext(const std::vector<int> &qubits, Sampler<T> &local, size_t fallbackLimit,
                                bool requestMoments = false)
    {
        const size_t limit = EffectiveSamplingLimit();
        if (!limit || qubits.size() > 4096 || !SamplingEligibility())
        {
            local.limit = fallbackLimit;
            if (fallbackLimit)
                local.nodes.emplace_back();
            return local;
        }
        if (!samplingCache || samplingCache->revision != samplingRevision || samplingCache->width != nrQubits ||
            samplingCache->qubits != qubits || samplingCache->limit != limit ||
            samplingCache->program.lock().get() != operations.get())
        {
            auto replacement = std::make_unique<Sampler<T>>();
            replacement->program = operations;
            replacement->qubits = qubits;
            replacement->revision = samplingRevision;
            replacement->width = nrQubits;
            replacement->limit = limit;
            replacement->nodes.emplace_back();
            samplingCache = std::move(replacement);
        }
        // WithTerm selects a unique packed type for each width interval. Width
        // and circuit/configuration revisions are part of the context key.
        auto &sampler = static_cast<Sampler<T> &>(*samplingCache);
        if (requestMoments && CanUseMoments(qubits.size()))
            PrepareMoments(sampler, qubits.size());
        return sampler;
    }

    template <class T> static void PrepareMoments(Sampler<T> &sampler, size_t count)
    {
        if (!sampler.moments.empty())
            return;
        // Commit both arrays together, keeping a valid prefix cache if either
        // allocation fails during a late promotion of a persistent context.
        std::vector<double> values(size_t{1} << count);
        std::vector<unsigned char> ready(values.size());
        values[0] = 1.;
        ready[0] = 1;
        sampler.moments = std::move(values);
        sampler.momentReady = std::move(ready);
    }

    template <class T> static void TrimSamplingScratch(Sampler<T> &sampler)
    {
        auto &ws = sampler.workspace;
        constexpr size_t maximumBytes = 16 * 1024 * 1024;
        size_t bytes = (ws.terms.capacity() + sampler.prefix.capacity() + sampler.branch.capacity()) * sizeof(T) +
                       ws.extra.capacity() * sizeof(std::vector<T>) + ws.hashSlots.capacity() * sizeof(size_t);
        for (const auto &extra : ws.extra)
            bytes += extra.capacity() * sizeof(T);
        if (bytes > maximumBytes)
        {
            ws = PauliDetail::Workspace<T>{};
            sampler.prefix = std::vector<T>{};
            sampler.branch = std::vector<T>{};
        }
        else
        {
            ws.terms.clear();
            for (auto &extra : ws.extra)
                extra.clear();
            sampler.prefix.clear();
            sampler.branch.clear();
        }
    }

    template <class T> struct SamplingScratchGuard
    {
        Sampler<T> &sampler;

        explicit SamplingScratchGuard(Sampler<T> &value) : sampler(value)
        {
        }

        ~SamplingScratchGuard()
        {
            TrimSamplingScratch(sampler);
        }
    };

    template <class T> double Moment(const std::vector<int> &qubits, Sampler<T> &sampler, size_t subset)
    {
        if (sampler.momentReady[subset])
            return sampler.moments[subset];
        auto &ws = sampler.workspace;
        ws.terms.clear();
        ws.terms.emplace_back(nrQubits);
        for (size_t bit = 0; bit < qubits.size(); ++bit)
            if (subset & (size_t{1} << bit))
                ws.terms[0].Z[qubits[bit] / 64] |= uint64_t{1} << (qubits[bit] % 64);
        const double value = PauliDetail::Execute(ws, *operations, Settings(), false);
        sampler.moments[subset] = value;
        sampler.momentReady[subset] = 1;
        return value;
    }

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
        size_t chosenBits = 0;
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
            else if (!sampler.moments.empty())
            {
                numerator = 0.;
                for (size_t subset = 0; subset < (size_t{1} << depth); ++subset)
                {
                    const double value = Moment(qubits, sampler, subset | (size_t{1} << depth));
                    numerator +=
                        (PauliDetail::Popcount(static_cast<uint64_t>(subset & chosenBits)) & 1) ? -value : value;
                }
                if (node != noNode)
                {
                    sampler.nodes[node].numerator = numerator;
                    sampler.nodes[node].ready = true;
                }
            }
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
                numerator = PauliDetail::Execute(ws, *operations, settings, false);
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
            if (one && !sampler.moments.empty())
                chosenBits |= size_t{1} << depth;
            if (depth + 1 == qubits.size())
                break;
            denominator += one ? -numerator : numerator;
            if (node != noNode)
            {
                size_t child = sampler.nodes[node].child[one ? 1 : 0];
                if (child == noNode && sampler.nodes.size() < limit)
                {
                    child = sampler.nodes.size();
                    sampler.nodes.emplace_back();
                    // Persist the edge only after allocation succeeds: this
                    // sampler may be reused after a failed sampling request.
                    sampler.nodes[node].child[one ? 1 : 0] = child;
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
    // Nonconst sampling already mutates rng. Its probability cache is local to
    // this instance; const expectation queries continue to use thread leases.
    uint64_t samplingRevision = 0;
    int samplingEligibility = -1;
    bool importedSamplingTables = false;
    std::unique_ptr<SamplingCacheBase> samplingCache;
};
} // namespace QC
