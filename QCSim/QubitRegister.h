#pragma once

#include "QubitRegisterCalculator.h"

#include <cstdint>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them
// usually from left to right

namespace QC
{

template <class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd>
class QubitRegister : public QubitRegisterCalculator<VectorClass, MatrixClass>,
                      public StructuredGateOperations<QubitRegister<VectorClass, MatrixClass>, MatrixClass>
{
  public:
    using GateClass = Gates::QuantumGateWithOp<MatrixClass>;
    using BaseClass = QubitRegisterCalculator<VectorClass, MatrixClass>;
    using StructuredGateOperations<QubitRegister<VectorClass, MatrixClass>, MatrixClass>::ApplyDiagonalGate;

    QubitRegister(size_t N = 3, unsigned int addseed = 0)
        : NrQubits(N), NrBasisStates(1ULL << NrQubits), registerStorage(VectorClass::Zero(NrBasisStates)),
          uniformZeroOne(0, 1), recordGates(false)
    {
        assert(N > 0);

        if (addseed == 0)
        {
            std::random_device rdl;
            addseed = rdl();
        }

        const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
        std::seed_seq seed{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(seed);

        registerStorage(0) = 1;
    }

    // this is a special constructor, I need it for something in a derived work (closed source)
    QubitRegister(size_t N, VectorClass &v, unsigned int addseed = 0)
        : NrQubits(N), NrBasisStates(1ULL << NrQubits), uniformZeroOne(0, 1), recordGates(false)
    {
        assert(N > 0);
        registerStorage.swap(v);

        if (addseed == 0)
        {
            std::random_device rd;
            addseed = rd();
        }

        const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
        std::seed_seq seed{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(seed);
    }

    // Seeds the engine with the value itself. A std::seed_seq costs about ten
    // microseconds per call, too much for callers that reseed for every shot.
    void SetSeed(uint64_t theSeed)
    {
        rng.seed(theSeed);
    }

    size_t getNrQubits() const
    {
        return NrQubits;
    };

    size_t getNrBasisStates() const
    {
        return NrBasisStates;
    };

    std::complex<double> getBasisStateAmplitude(size_t State) const
    {
        if (State >= NrBasisStates)
            return 0;

        return registerStorage(State);
    }

    double getBasisStateProbability(size_t State) const
    {
        if (State >= NrBasisStates)
            return 0;

        return std::norm(registerStorage(State));
    }

    void setToBasisState(size_t State)
    {
        if (State >= NrBasisStates)
            return;

        Clear();
        registerStorage(State) = 1;
    }

    void setToQubitState(size_t q)
    {
        if (q >= NrQubits)
            return;

        Clear();
        registerStorage(1ULL << q) = 1;
    }

    // measurement should give either all 0 or all 1
    void setToCatState()
    {
        Clear();
        static const double OneOverSqrt2 = 1. / sqrt(2.);

        registerStorage(0) = OneOverSqrt2;
        registerStorage(NrBasisStates - 1) = OneOverSqrt2;
    }

    void Reset()
    {
        setToBasisState(0);
    }

    // all states have equal amplitude, so measurement should give any state with equal probability
    void setToEqualSuperposition()
    {
        registerStorage.setConstant(1. / sqrt(NrBasisStates));
    }

    // to be able to set them all, after setting them, call Normalize
    void setRawAmplitude(size_t State, std::complex<double> val)
    {
        if (State >= NrBasisStates)
            return;

        registerStorage(State) = val;
    }

    void Clear()
    {
        registerStorage.setZero();
    }

    void Normalize()
    {
        const double norm = registerStorage.norm();
        if (norm < 1E-20)
            return;

        registerStorage *= 1. / norm;
    }

    // to be able to compare different results
    void AdjustPhaseAndNormalize()
    {
        std::complex<double> v0 = registerStorage[0];
        double av0 = abs(v0);

        if (av0 >= 1E-5)
        {
            for (size_t i = 0; i < getNrBasisStates(); ++i)
                registerStorage[i] /= v0;
        }
        else
        {
            v0 = registerStorage[NrBasisStates >> 1];
            av0 = abs(v0);

            if (av0 >= 1E-5)
            {
                for (size_t i = 0; i < getNrBasisStates(); ++i)
                    registerStorage[i] /= v0;
            }
            else
            {
                v0 = registerStorage[NrBasisStates - 1];
                av0 = abs(v0);

                if (av0 >= 1E-5)
                {
                    for (size_t i = 0; i < getNrBasisStates(); ++i)
                        registerStorage[i] /= v0;
                }
            }
        }

        Normalize();
    }

    size_t MeasureAll()
    {
        const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy
        const size_t state =
            BaseClass::SampleBasisState(NrBasisStates, registerStorage, prob, NrBasisStates - 1, UseMultithreading());

        setToBasisState(state); // collapse

        return state;
    }

    // shortcut for measuring a single qubit
    size_t MeasureQubit(size_t qubit)
    {
        return Measure(qubit, qubit);
    }

    // measure a 'subregister' as a separate register
    // can measure a single qubit, if firstQubit == secondQubit
    // will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on
    // position 0 and so on) so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other
    // measured ones are zero, 2 means that the next one 1 one and all others are zero and so on

    size_t Measure(size_t firstQubit, size_t secondQubit)
    {
        const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy

        if (firstQubit == secondQubit)
        {
            if (!UseMultithreading())
                return BaseClass::MeasureQubit(NrBasisStates, registerStorage, firstQubit, prob);

            return BaseClass::MeasureQubitOmp(NrBasisStates, registerStorage, firstQubit, prob);
        }

        if (!UseMultithreading())
            return BaseClass::Measure(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);

        return BaseClass::MeasureOmp(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
    }

    std::map<size_t, size_t> RepeatedMeasure(size_t nrTimes = 1000)
    {
        if (nrTimes == 0)
            return {}; // nothing to measure

        std::map<size_t, size_t> measurements;

        if (nrTimes == 1) // shortcut for a single measurement
        {
            const size_t meas = MeasureNoCollapse();
            ++measurements[meas];
            return measurements;
        }

        // a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N),
        // where N is the number of states), there is an even O(1) method (see
        // https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

        std::vector<double> probabilities(registerStorage.size());

        double accum = 0;
        for (size_t i = 0; i < static_cast<size_t>(registerStorage.size()); ++i)
        {
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
            ++i;
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
        }

        for (size_t shot = 0; shot < nrTimes; ++shot)
        {
            const double prob = 1. - uniformZeroOne(rng);
            const size_t meas =
                std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

            ++measurements[meas];
        }

        return measurements;
    }

    std::unordered_map<size_t, size_t> RepeatedMeasureUnordered(size_t nrTimes = 1000)
    {
        if (nrTimes == 0)
            return {}; // nothing to measure

        std::unordered_map<size_t, size_t> measurements;

        if (nrTimes == 1) // shortcut for a single measurement
        {
            const size_t meas = MeasureNoCollapse();
            ++measurements[meas];
            return measurements;
        }

        // a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N),
        // where N is the number of states), there is an even O(1) method (see
        // https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

        std::vector<double> probabilities(registerStorage.size());

        double accum = 0;
        for (size_t i = 0; i < (size_t)registerStorage.size(); ++i)
        {
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
            ++i;
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
        }

        for (size_t shot = 0; shot < nrTimes; ++shot)
        {
            const double prob = 1. - uniformZeroOne(rng);
            const size_t meas =
                std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

            ++measurements[meas];
        }

        return measurements;
    }

    std::map<size_t, size_t> RepeatedMeasure(size_t firstQubit, size_t secondQubit, size_t nrTimes = 1000)
    {
        if (nrTimes == 0)
            return {}; // nothing to measure

        std::map<size_t, size_t> measurements;

        const size_t secondQubitp1 = secondQubit + 1;
        const size_t firstPartMask = (1ULL << firstQubit) - 1;
        const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;

        if (nrTimes == 1) // shortcut for a single measurement
        {
            const size_t meas = MeasureNoCollapse(firstQubit, secondQubit);
            ++measurements[(meas & measuredPartMask) >> firstQubit];
            return measurements;
        }

        // a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N),
        // where N is the number of states), there is an even O(1) method (see
        // https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

        std::vector<double> probabilities(registerStorage.size());

        double accum = 0;
        for (size_t i = 0; i < static_cast<size_t>(registerStorage.size()); ++i)
        {
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
            ++i;
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
        }

        for (size_t shot = 0; shot < nrTimes; ++shot)
        {
            const double prob = 1. - uniformZeroOne(rng);
            const size_t meas =
                std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

            ++measurements[(meas & measuredPartMask) >> firstQubit];
        }

        return measurements;
    }

    std::unordered_map<size_t, size_t> RepeatedMeasureUnordered(size_t firstQubit, size_t secondQubit,
                                                                size_t nrTimes = 1000)
    {
        if (nrTimes == 0)
            return {}; // nothing to measure

        std::unordered_map<size_t, size_t> measurements;

        const size_t secondQubitp1 = secondQubit + 1;
        const size_t firstPartMask = (1ULL << firstQubit) - 1;
        const size_t measuredPartMask = (1ULL << secondQubitp1) - 1 - firstPartMask;

        if (nrTimes == 1) // shortcut for a single measurement
        {
            const size_t meas = MeasureNoCollapse(firstQubit, secondQubit);
            ++measurements[(meas & measuredPartMask) >> firstQubit];
            return measurements;
        }

        // a faster sampling way (O(n) where n is the number of qubits, except the preprocessing phase which is O(N),
        // where N is the number of states), there is an even O(1) method (see
        // https://en.wikipedia.org/wiki/Alias_method), but I won't bother here, it's not a so spectacular improvement

        std::vector<double> probabilities(registerStorage.size());

        double accum = 0;
        for (size_t i = 0; i < registerStorage.size(); ++i)
        {
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
            ++i;
            accum += std::norm(registerStorage[i]);
            probabilities[i] = accum;
            if (accum > 1.0 - std::numeric_limits<double>::epsilon())
            {
                probabilities.resize(i + 1);
                break;
            }
        }

        for (size_t shot = 0; shot < nrTimes; ++shot)
        {
            const double prob = 1. - uniformZeroOne(rng);
            const size_t meas =
                std::lower_bound(probabilities.begin(), probabilities.end(), prob) - probabilities.begin();

            ++measurements[(meas & measuredPartMask) >> firstQubit];
        }

        return measurements;
    }

    // controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for
    // gates with a lower number of qubits
    void ApplyGate(const GateClass &gate, size_t qubit, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0)
    {
        const auto &matrix = gate.getRawOperatorMatrix();
        const auto dimension = matrix.rows();
        if (dimension != matrix.cols() || (dimension != 2 && dimension != 4 && dimension != 8))
            throw std::invalid_argument("Statevector gates must have a 2x2, 4x4 or 8x8 operator matrix");
        const size_t gateQubits = gate.getQubitsNumber();
        const size_t matrixQubits = dimension == 2 ? 1 : (dimension == 4 ? 2 : 3);
        if (gateQubits != matrixQubits)
            throw std::invalid_argument("Gate arity does not match its operator matrix");

        CheckQubits(gate, qubit, controllingQubit1, controllingQubit2, gateQubits);

#define OPTIMIZED_TENSOR_PRODUCT 1
#ifdef OPTIMIZED_TENSOR_PRODUCT
        assert(gateQubits > 0 && gateQubits <= 3);

        const std::array<size_t, 3> bits{size_t{1} << qubit, gateQubits > 1 ? size_t{1} << controllingQubit1 : 0,
                                         gateQubits > 2 ? size_t{1} << controllingQubit2 : 0};
        BaseClass::ApplyGateInPlace(registerStorage, matrix, gate.getStructure(), bits,
                                    static_cast<unsigned>(gateQubits), NrBasisStates, UseMultithreading());
#else
        registerStorage =
            gate.getOperatorMatrix(NrQubits, qubit, controllingQubit1, controllingQubit2) * registerStorage;
#endif

        if (recordGates)
            computeGates.emplace_back(
                Gates::AppliedGate(gate.getRawOperatorMatrix(), qubit, controllingQubit1, controllingQubit2));
    }

    void ApplyGate(const Gates::AppliedGate<MatrixClass> &gate)
    {
        if (const auto *payload = gate.getCompactPayload())
        {
            // The kernels read the shared tables of the gate under its own qubits:
            // a remapped or adjoint gate needs no copy of them. The tables were
            // validated when the gate was made and cannot change afterwards, so
            // only the qubits are checked; the cosines and sines of a
            // multiplexer are computed once and stay with the tables.
            const bool multiplexed = payload->OperationKind() == Gates::StructuredGate::Kind::MultiplexedRotation;
            ApplyStructuredGateImpl(payload->View(), &gate, true, payload->conjugated,
                                    multiplexed ? payload->Coefficients(true) : nullptr);
        }
        else if (const auto *compact = gate.getStructuredGate())
            ApplyStructuredGateImpl(*compact, &gate, false);
        else
            ApplyGate(gate, gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
    }

    void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>> &gates)
    {
        for (const auto &gate : gates)
            ApplyGate(gate);
    }

    void ApplyStructuredGate(const Gates::StructuredGate &gate)
    {
        ApplyStructuredGateImpl(gate, nullptr, false);
    }

    // The same for an operation whose tables stay with the caller: nothing is
    // copied unless the register is recording its gates.
    void ApplyStructuredGateView(const Gates::StructuredGateView &gate)
    {
        ApplyStructuredGateImpl(gate, nullptr, false);
    }

    void ApplyOperatorMatrix(const MatrixClass &m)
    {
        registerStorage = m * registerStorage;

        if (recordGates)
            computeGates.emplace_back(Gates::AppliedGate(m, 0, 1, 2));
    }

    const VectorClass &getRegisterStorage() const
    {
        return registerStorage;
    }

    void setRegisterStorage(const VectorClass &vals)
    {
        if (registerStorage.size() != vals.size())
            return;

        registerStorage = vals;
        Normalize();
    }

    // warning, you should be sure that the vector is normalized and has the proper size
    void setRegisterStorageFastNoNormalize(VectorClass &vals)
    {
        registerStorage.swap(vals);
    }

    // to check how well the computed state matches some 'exact' known one
    double stateFidelity(const VectorClass &state) const
    {
        if (registerStorage.size() != state.size())
            return 0;

        const std::complex<double> p = (registerStorage.adjoint() * state)(0);

        return norm(p);
    }

    void ComputeStart()
    {
        recordGates = true;
        computeGates.clear();
    }

    void ComputeEnd()
    {
        recordGates = false;
    }

    void ComputeClear()
    {
        computeGates.clear();
    }

    // applies again the recorded gates
    // with this the same operations can be repeated several times
    void Compute()
    {
        // avoid recording the gates again if somehow the user forgot to call ComputeEnd
        const bool recordSave = recordGates;
        recordGates = false;

        try
        {
            for (const auto &gate : computeGates)
                ReplayRecordedGate(gate, false);
        }
        catch (...)
        {
            recordGates = recordSave;
            throw;
        }

        recordGates = recordSave;
    }

    // undoes the recorded gates
    // the operations are unitary, so U^-1 = U^t and (U1 * U2)^t = U2^t * U1^t
    void Uncompute()
    {
        const bool recordSave = recordGates;
        recordGates = false;

        try
        {
            for (auto gate = computeGates.rbegin(); gate != computeGates.rend(); ++gate)
                ReplayRecordedGate(*gate, true);
        }
        catch (...)
        {
            recordGates = recordSave;
            throw;
        }

        recordGates = recordSave;
    }

    double GetQubitProbability(size_t qubit) const
    {
        if (!UseMultithreading())
            return BaseClass::GetQubitProbability(NrBasisStates, registerStorage, qubit);

        return BaseClass::GetQubitProbabilityOmp(NrBasisStates, registerStorage, qubit);
    }

    // Reduced density matrix of one, two or three qubits, in one read-only pass. The first qubit is the least
    // significant bit of the row and of the column index.
    MatrixClass LocalDensityMatrix(size_t qubit) const
    {
        return LocalDensity<1>({qubit});
    }

    MatrixClass LocalDensityMatrix(size_t qubit, size_t qubit2) const
    {
        return LocalDensity<2>({qubit, qubit2});
    }

    MatrixClass LocalDensityMatrix(size_t qubit, size_t qubit2, size_t qubit3) const
    {
        return LocalDensity<3>({qubit, qubit2, qubit3});
    }

    // Raw squared-amplitude sums in the same basis order as LocalDensityMatrix. Their sum is the state's
    // squared norm, including for an unnormalized state.
    std::array<double, 2> LocalPopulations(size_t qubit) const
    {
        return Populations<1>({qubit});
    }

    std::array<double, 4> LocalPopulations(size_t qubit, size_t qubit2) const
    {
        return Populations<2>({qubit, qubit2});
    }

    std::array<double, 8> LocalPopulations(size_t qubit, size_t qubit2, size_t qubit3) const
    {
        return Populations<3>({qubit, qubit2, qubit3});
    }

    // Applies a 2x2, 4x4 or 8x8 matrix, which need not be unitary, to one, two or three qubits and returns the
    // squared norm of the resulting, unnormalized state. The first qubit is the least significant bit of the
    // matrix index. A trajectory can apply the no-jump operators of a block and obtain their probability in
    // the same pass. A very small squared norm is recomputed from the output to recover underflowed terms.
    double ApplyLocalMapNorm(const MatrixClass &matrix, size_t qubit)
    {
        return LocalNorm<1>(matrix, {qubit});
    }

    double ApplyLocalMapNorm(const MatrixClass &matrix, size_t qubit, size_t qubit2)
    {
        return LocalNorm<2>(matrix, {qubit, qubit2});
    }

    double ApplyLocalMapNorm(const MatrixClass &matrix, size_t qubit, size_t qubit2, size_t qubit3)
    {
        return LocalNorm<3>(matrix, {qubit, qubit2, qubit3});
    }

    // One step of a quantum trajectory through the channel rho -> sum_k E_k rho E_k^dagger: applies one Kraus
    // operator, chosen with the probability the state gives it, and renormalizes. Averaging the trajectories
    // reproduces the channel. The operators are 2x2, 4x4 or 8x8, on one, two or three qubits respectively.
    // The first target is the low bit of the matrix index. `random` is uniform in [0, 1).
    //
    // The probabilities tr(E_k rho E_k^dagger) come from the reduced density matrix of the targets. When every
    // operator has at most one nonzero per row (as in amplitude damping), only the target populations are
    // needed. Real operators need only the real part of the reduced density matrix. In the ordinary range,
    // one pass supplies the weights for all operators and one gate application produces the selected branch.
    // Trace preservation is not checked: the weights are taken relative to their sum, and the state leaves with
    // unit norm. Returns the index of the operator applied.
    //
    // Precision. A weight that is a small part of what its terms add up to is read again from the amplitudes (see
    // UnresolvedWeight). What the amplitudes cannot tell from their own rounding counts as no weight: one below
    // (64 epsilon)^2 times the bound on those terms, that is, an amplitude left below about 1.4e-14 of the
    // amplitudes that cancelled to give it. An operator without weight is never drawn. When none is left the call
    // throws std::domain_error and preserves the register. The diagnostic includes the largest directly measured
    // weight and its cutoff when available. The cutoff also applies to exact binary inputs: roundoff in a state
    // produced by gates can leave a residual of the same size for a branch whose true weight is zero.
    // When amplitudes or coefficients reach the limits of the scalar range,
    // powers of two keep the products and normalization representable. In that case the cumulative draw is
    // exact for the computed weights; the dot products giving those weights still have ordinary roundoff.
    size_t ApplyKrausOperator(const std::vector<MatrixClass> &kraus, double random, size_t qubit,
                              size_t controllingQubit1 = 0)
    {
        if (!kraus.empty() && kraus.front().rows() == 2)
            return ApplyKraus<1>(kraus, random, {qubit});
        return ApplyKraus<2>(kraus, random, {qubit, controllingQubit1});
    }

    size_t ApplyKrausOperator(const std::vector<MatrixClass> &kraus, double random, size_t qubit, size_t qubit2,
                              size_t qubit3)
    {
        return ApplyKraus<3>(kraus, random, {qubit, qubit2, qubit3});
    }

    // The same, drawing from this register's own random stream.
    size_t SampleKrausOperator(const std::vector<MatrixClass> &kraus, size_t qubit, size_t controllingQubit1 = 0)
    {
        return ApplyKrausOperator(kraus, uniformZeroOne(rng), qubit, controllingQubit1);
    }

    size_t SampleKrausOperator(const std::vector<MatrixClass> &kraus, size_t qubit, size_t qubit2, size_t qubit3)
    {
        return ApplyKrausOperator(kraus, uniformZeroOne(rng), qubit, qubit2, qubit3);
    }

  private:
    template <unsigned Qubits>
    size_t ApplyKraus(const std::vector<MatrixClass> &kraus, double random, const std::array<size_t, Qubits> &qubits)
    {
        constexpr Eigen::Index dimension = Eigen::Index{1} << Qubits;
        if (kraus.empty())
            throw std::invalid_argument("A channel must contain at least one Kraus operator");
        for (const auto &e : kraus)
            if (e.rows() != dimension || e.cols() != dimension || !e.allFinite())
                throw std::invalid_argument(
                    "Kraus operators must be finite square matrices matching their target qubits");

        // Stack storage for the usual few operators: the weight of each, and its reach when the weight is in doubt.
        double fewWeights[16], fewReaches[16];
        std::vector<double> manyValues(kraus.size() > 16 ? 2 * kraus.size() : 0);
        double *const weights = kraus.size() > 16 ? manyValues.data() : fewWeights;
        double *const reaches = kraus.size() > 16 ? manyValues.data() + kraus.size() : fewReaches;
        const bool usePopulations = std::all_of(kraus.begin(), kraus.end(), [](const MatrixClass &op) {
            // Distinct input columns have disjoint row support, so E^dagger E is diagonal without cancellation.
            for (Eigen::Index row = 0; row < op.rows(); ++row)
            {
                bool occupied = false;
                for (Eigen::Index col = 0; col < op.cols(); ++col)
                    if (op(row, col) != typename MatrixClass::Scalar(0))
                    {
                        if (occupied)
                            return false;
                        occupied = true;
                    }
            }
            return true;
        });
        // A population rounded to zero can hide a finite contribution from a large coefficient.
        // With bounded Kraus coefficients its contribution is itself below the weight range, and
        // only the tiny-weight or boundary checks below need to recover it.
        const auto unreliablePopulation = [&](double p, size_t c) {
            if (p != 0.)
                return !std::isfinite(p) || p < std::numeric_limits<double>::min();
            return std::any_of(kraus.begin(), kraus.end(), [&](const MatrixClass &op) {
                return op.col(c).cwiseAbs().maxCoeff() > 1.;
            });
        };
        MatrixClass rho;
        if (usePopulations)
        {
            const auto p = Populations<Qubits>(qubits);
            for (size_t c = 0; c < dimension; ++c)
                if (unreliablePopulation(p[c], c))
                    return ApplyScaledKraus<Qubits>(kraus, random, qubits);
            KrausPopulationWeights(kraus, p, weights);
        }
        else
        {
            const bool real = std::all_of(kraus.begin(), kraus.end(),
                                          [](const MatrixClass &op) { return (op.imag().array() == 0.).all(); });
            rho = real ? LocalDensity<Qubits, true>(qubits) : LocalDensity<Qubits>(qubits);
            for (size_t c = 0; c < dimension; ++c)
                if (unreliablePopulation(double(rho(c, c).real()), c))
                    return ApplyScaledKraus<Qubits>(kraus, random, qubits);
        }
        double total = 0., doubt = 0.;
        for (size_t k = 0; k < kraus.size(); ++k)
        {
            if (usePopulations)
            {
                reaches[k] = -1.; // A sum of nonnegative terms cannot cancel.
                total += weights[k];
                continue;
            }
            // The weight, and what the terms of its sum can add up to: a bound on them, d * sum_a |column a of E|^2
            // rho(a, a), by Cauchy-Schwarz. The entries of rho are known to their own rounding, so the weight is
            // known to that rounding of its reach. A weight that is a small fraction of its reach is what the
            // cancellation of the terms left (see UnresolvedWeight), whatever its sign. The reach is d times the
            // terms with a = b, so it comes out of the same pass.
            double weight = 0., diagonal = 0.;
            for (Eigen::Index row = 0; row < dimension; ++row)
                for (Eigen::Index a = 0; a < dimension; ++a)
                    for (Eigen::Index b = 0; b < dimension; ++b)
                    {
                        const double term =
                            std::real(std::complex<double>(kraus[k](row, a)) * std::complex<double>(rho(a, b)) *
                                      std::conj(std::complex<double>(kraus[k](row, b))));
                        weight += term;
                        if (a == b)
                            diagonal += term;
                    }
            const double reach = static_cast<double>(dimension) * diagonal;
            if (!std::isfinite(weight) || !std::isfinite(reach))
                return ApplyScaledKraus<Qubits>(kraus, random, qubits);
            weights[k] = weight > 0. ? weight : 0.;
            if (weight >= UnresolvedWeight * reach)
                reaches[k] = -1.; // resolved
            else
            {
                reaches[k] = reach;
                doubt += DensityRounding * reach;
            }
            total += weights[k];
        }

        size_t chosen = kraus.size();
        bool boundaryAtZero = false;
        // Draws the operator for `random`; says whether the threshold is clear of the ends of its interval by more
        // than the weights in doubt can move them.
        const auto draw = [&]() {
            const double threshold = random * total;
            double accumulated = 0., lower = 0.;
            boundaryAtZero = false;
            for (size_t k = 0; k < kraus.size(); ++k)
            {
                if (weights[k] <= 0.)
                    continue;
                chosen = k; // also the fallback if rounding leaves the threshold above the last sum
                lower = accumulated;
                accumulated += weights[k];
                if (threshold < accumulated)
                {
                    if (threshold == lower && random > 0.)
                        boundaryAtZero = std::any_of(weights, weights + kraus.size(), [&](double w) {
                            return w <= std::numeric_limits<double>::epsilon() * total;
                        });
                    break;
                }
            }
            return threshold - lower > 2. * doubt && accumulated - threshold > 2. * doubt;
        };
        bool drawn = false;
        // The largest weight read from the amplitudes, and the cutoff of its operator: for the message, if none stays.
        double largestFound = 0., cutoffOfLargest = 0.;
        if (doubt > 0.)
        {
            // Some weights are not resolved. Usually that changes nothing: they are too small to move the draw, and
            // it is settled as it stands. Otherwise each of them is read from the amplitudes, without changing the
            // state, before anything is drawn: a weight that is there is then found whatever the rounding of rho made
            // of it, and one that is not there is not drawn.
            drawn = std::isfinite(total) && total > 0. && draw();
            if (!drawn)
            {
                total = 0.;
                for (size_t k = 0; k < kraus.size(); ++k)
                {
                    if (reaches[k] >= 0.)
                    {
                        const double exact = LocalMapWeight<Qubits>(kraus[k], qubits);
                        const double cutoff = AmplitudeRounding * AmplitudeRounding * reaches[k];
                        if (exact > largestFound)
                        {
                            largestFound = exact;
                            cutoffOfLargest = cutoff;
                        }
                        weights[k] = exact >= cutoff ? exact : 0.;
                    }
                    total += weights[k];
                }
                doubt = 0.;
            }
        }
        if (!drawn)
        {
            // Tiny squared norms can round to zero or retain only a few significant bits. This can
            // affect a draw even in a trace-preserving family if the first interval is that small.
            constexpr double minimum = std::numeric_limits<double>::min();
            if (!std::isfinite(total) || (total <= minimum && largestFound == 0.) ||
                (random * total <= minimum && weights[0] <= minimum && total > 0.))
                return ApplyScaledKraus<Qubits>(kraus, random, qubits);
            if (total <= 0.)
            {
                // A positive weight that was read from the amplitudes and fell below its cutoff is reported: it is
                // not told from the rounding of amplitudes that cancel, and counts as none.
                if (!(largestFound > 0.))
                    throw std::domain_error("No Kraus operator has any probability on this state");
                std::ostringstream message;
                message << "No Kraus operator has any probability on this state: the largest weight found, "
                        << std::setprecision(3) << largestFound << ", is below " << cutoffOfLargest
                        << ", the rounding of amplitudes that cancel, and counts as none";
                throw std::domain_error(message.str());
            }
            draw();
        }
        if (boundaryAtZero || weights[chosen] < std::numeric_limits<double>::min())
            return ApplyScaledKraus<Qubits>(kraus, random, qubits);

        // An operator whose weight was in doubt takes the state almost to nothing: what it leaves is the difference
        // of nearly equal amplitudes. It is applied as it is, as a gate, and the state is then scaled by the norm it
        // actually has, as a second gate. A recording holds the two, and replays exactly what was done here.
        const size_t qubit2 = Qubits > 1 ? qubits[1] : 0, qubit3 = Qubits > 2 ? qubits[2] : 0;
        if (reaches[chosen] >= 0.)
        {
            ApplyGate(Gates::AppliedGate<MatrixClass>(kraus[chosen], qubits[0], qubit2, qubit3));
            const double norm = registerStorage.squaredNorm();
            if (!(norm > 0.) || !std::isfinite(norm))
                throw std::domain_error("The Kraus operator drawn leaves nothing of this state");
            ApplyGate(
                Gates::AppliedGate<MatrixClass>(MatrixClass(MatrixClass::Identity(2, 2) / std::sqrt(norm)), qubits[0]));
            return chosen;
        }

        // The scaled operator goes straight to the local-map kernel: no gate object, and a diagonal or
        // single-entry operator (the no-jump and jump operators of a damping channel) costs what it holds.
        const double scale = 1. / std::sqrt(weights[chosen]);
        std::array<std::complex<double>, dimension * dimension> coefficients;
        for (Eigen::Index row = 0; row < dimension; ++row)
            for (Eigen::Index col = 0; col < dimension; ++col)
                coefficients[static_cast<size_t>(row * dimension + col)] =
                    std::complex<double>(kraus[chosen](row, col)) * scale;
        // An unoccupied input column can hold a finite coefficient too large to premultiply by scale.
        if (scale > 1. && !std::all_of(coefficients.begin(), coefficients.end(), [](std::complex<double> z) {
                return std::isfinite(z.real()) && std::isfinite(z.imag());
            }))
            return ApplyScaledKraus<Qubits>(kraus, random, qubits);
        BaseClass::template ApplyLocalMap<Qubits>(registerStorage, coefficients.data(), LocalBits<Qubits>(qubits),
                                                  NrBasisStates, UseMultithreading());
        if (recordGates)
            computeGates.emplace_back(
                Gates::AppliedGate<MatrixClass>(MatrixClass(kraus[chosen] * scale), qubits[0], qubit2, qubit3));
        return chosen;
    }

    // Exceptional range only. Factor a power of two from each populated input column and a second
    // power from each operator's image. The matrix and state seen by the ordinary dot products then
    // have bounded entries, even when the original probability lies outside the double range.
    template <unsigned Qubits>
    size_t ApplyScaledKraus(const std::vector<MatrixClass> &kraus, double random,
                            const std::array<size_t, Qubits> &qubits)
    {
        constexpr size_t dimension = size_t{1} << Qubits;
        const auto bits = LocalBits<Qubits>(qubits);
        const auto column = [&](size_t index) {
            size_t c = 0;
            for (unsigned b = 0; b < Qubits; ++b)
                if (index & bits[b])
                    c |= size_t{1} << b;
            return c;
        };
        const auto magnitude = [](std::complex<double> z) { return std::max(std::abs(z.real()), std::abs(z.imag())); };
        const auto shift = [](std::complex<double> z, int exponent) {
            return std::complex<double>(std::ldexp(z.real(), exponent), std::ldexp(z.imag(), exponent));
        };
        if (!registerStorage.allFinite())
            throw std::domain_error("Kraus weights require a finite state");
        std::array<double, dimension> maxima{}, populations{};
        std::array<int, dimension> columns{};
        for (size_t i = 0; i < NrBasisStates; ++i)
        {
            const size_t c = column(i);
            maxima[c] = std::max(maxima[c], magnitude(registerStorage[i]));
        }
        for (size_t c = 0; c < dimension; ++c)
            if (maxima[c] > 0.)
                columns[c] = std::ilogb(maxima[c]);
        VectorClass scaled = registerStorage;
        for (size_t i = 0; i < NrBasisStates; ++i)
        {
            const size_t c = column(i);
            const auto z = shift(registerStorage[i], -columns[c]);
            scaled[i] = typename VectorClass::Scalar(z);
            populations[c] += std::norm(z);
        }

        constexpr int precision = std::numeric_limits<double>::digits;
        struct Weight
        {
            std::uint64_t significand = 0;
            int exponent = 0, imageExponent = 0;
        };
        std::vector<Weight> weights(kraus.size());
        int smallestExponent = std::numeric_limits<int>::max();
        int largestExponent = std::numeric_limits<int>::min();
        VectorClass image;
        std::array<std::complex<double>, dimension * dimension> coefficients;
        for (size_t k = 0; k < kraus.size(); ++k)
        {
            auto &w = weights[k];
            int imageExponent = std::numeric_limits<int>::min();
            for (size_t c = 0; c < dimension; ++c)
                if (maxima[c] > 0.)
                    for (size_t r = 0; r < dimension; ++r)
                    {
                        const double m = magnitude(kraus[k](r, c));
                        if (m > 0.)
                            imageExponent = std::max(imageExponent, columns[c] + std::ilogb(m));
                    }
            if (imageExponent == std::numeric_limits<int>::min())
                continue;
            w.imageExponent = imageExponent;
            double diagonal = 0.;
            for (size_t r = 0; r < dimension; ++r)
                for (size_t c = 0; c < dimension; ++c)
                {
                    const auto z = maxima[c] > 0. ? shift(kraus[k](r, c), columns[c] - imageExponent)
                                                  : std::complex<double>{};
                    coefficients[r * dimension + c] = z;
                    diagonal += std::norm(z) * populations[c];
                }
            image = scaled;
            const double weight = BaseClass::template ApplyLocalMapNorm<Qubits>(
                image, coefficients.data(), bits, NrBasisStates, UseMultithreading());
            if (weight <= 0. || weight < (AmplitudeRounding * AmplitudeRounding * dimension) * diagonal)
                continue;
            w.significand = static_cast<std::uint64_t>(std::ldexp(std::frexp(weight, &w.exponent), precision));
            w.exponent += 2 * imageExponent - precision;
            smallestExponent = std::min(smallestExponent, w.exponent);
            largestExponent = std::max(largestExponent, w.exponent);
        }
        if (largestExponent == std::numeric_limits<int>::min())
            throw std::domain_error("No Kraus operator has any probability on this state");

        // Each computed weight and the draw are dyadic. Compare the cumulative weight with
        // random * total in base 2^16, exactly for those values. Floating addition would discard
        // a tiny interval beside a large one and could move a draw across their shared boundary.
        int drawExponent = 0;
        const double fraction = std::frexp(std::min(1., std::max(0., random)), &drawExponent);
        const auto draw = static_cast<std::uint64_t>(std::ldexp(fraction, precision));
        drawExponent -= precision;
        const int low = smallestExponent + drawExponent;
        const int high = largestExponent + precision + std::numeric_limits<size_t>::digits;
        using Digits = std::vector<std::uint16_t>;
        Digits threshold((high - low) / 16 + 2), cumulative(threshold.size());
        const auto addProduct = [&](Digits &sum, const Weight &w, std::uint64_t multiplier, int shift) {
            for (std::uint64_t x = w.significand; x; x >>= 16, shift += 16)
            {
                int bit = w.exponent + shift - low;
                for (std::uint64_t y = multiplier; y; y >>= 16, bit += 16)
                {
                    std::uint64_t carry = ((x & 0xffffu) * (y & 0xffffu)) << (bit % 16);
                    for (size_t pos = static_cast<size_t>(bit / 16); carry; ++pos)
                    {
                        carry += sum[pos];
                        sum[pos] = static_cast<std::uint16_t>(carry);
                        carry >>= 16;
                    }
                }
            }
        };
        for (const auto &w : weights)
            addProduct(threshold, w, draw, drawExponent);
        size_t chosen = 0;
        for (size_t k = 0; k < weights.size(); ++k)
        {
            const auto &w = weights[k];
            if (w.significand == 0)
                continue;
            chosen = k;
            addProduct(cumulative, w, 1, 0);
            if (std::lexicographical_compare(threshold.rbegin(), threshold.rend(), cumulative.rbegin(), cumulative.rend()))
                break;
        }
        const auto &w = weights[chosen];
        MatrixClass bounded(dimension, dimension);
        for (size_t c = 0; c < dimension; ++c)
        {
            // Factor empty columns too, so a recording still represents the original linear operator
            // when replayed on a different state.
            if (maxima[c] == 0.)
            {
                double m = 0.;
                for (size_t r = 0; r < dimension; ++r)
                    m = std::max(m, magnitude(kraus[chosen](r, c)));
                if (m > 0.)
                    columns[c] = w.imageExponent - std::ilogb(m);
            }
            for (size_t r = 0; r < dimension; ++r)
                bounded(r, c) = typename MatrixClass::Scalar(shift(kraus[chosen](r, c), columns[c] - w.imageExponent));
        }
        const size_t q2 = Qubits > 1 ? qubits[1] : 0, q3 = Qubits > 2 ? qubits[2] : 0;
        // Apply to scratch before changing the register or its recording. In particular, keep the
        // normalization after the dot products: scaling their coefficients first can magnify rounding
        // when the selected branch almost cancels.
        image = registerStorage;
        const std::array<size_t, 3> gateBits{bits[0], Qubits > 1 ? bits[1] : 0, Qubits > 2 ? bits[2] : 0};
        std::vector<Gates::AppliedGate<MatrixClass>> recording;
        const auto apply = [&](const MatrixClass &matrix) {
            Gates::AppliedGate<MatrixClass> gate(matrix, qubits[0], q2, q3);
            BaseClass::ApplyGateInPlace(image, matrix, gate.getStructure(), gateBits, Qubits, NrBasisStates,
                                        UseMultithreading());
            if (recordGates)
                recording.emplace_back(std::move(gate));
        };
        // Split large exponents so every recorded diagonal is representable in the matrix scalar.
        constexpr int step = std::numeric_limits<typename MatrixClass::RealScalar>::max_exponent / 2;
        while (std::any_of(columns.begin(), columns.end(), [](int e) { return e != 0; }))
        {
            MatrixClass diagonal = MatrixClass::Zero(dimension, dimension);
            for (size_t c = 0; c < dimension; ++c)
            {
                const int e = std::clamp(columns[c], -step, step);
                diagonal(c, c) = std::ldexp(1., -e);
                columns[c] -= e;
            }
            apply(diagonal);
        }
        apply(bounded);
        const double norm = image.stableNorm();
        if (!(norm > 0.) || !std::isfinite(norm))
            throw std::domain_error("The Kraus operator drawn leaves nothing of this state");
        apply(MatrixClass::Identity(dimension, dimension) / norm);
        computeGates.insert(computeGates.end(), std::make_move_iterator(recording.begin()),
                            std::make_move_iterator(recording.end()));
        registerStorage.swap(image);
        return chosen;
    }

  public:
    void SaveState()
    {
        savedStateStorage = registerStorage;
    }

    void RestoreState()
    {
        if (savedStateStorage.size() == 0)
            return;
        registerStorage = savedStateStorage;
    }

    void RestoreStateDestructive()
    {
        if (savedStateStorage.size() == 0)
            return;
        registerStorage.swap(savedStateStorage);
        savedStateStorage.resize(0);
    }

    // the following ones should be used for 'repeated measurements' that avoid reexecuting the circuit each time
    size_t MeasureNoCollapse()
    {
        const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy

        return BaseClass::SampleBasisState(NrBasisStates, registerStorage, prob, 0, UseMultithreading());
    }

    // does not check the gates, that's why it returns a complex number
    // the caller should ensure the hermicity and extract the real part
    std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>> &gates)
    {
        if (gates.empty())
            return 1.;

        // Pauli strings (the common case) are computed in a single read-only pass, without copying the state
        size_t xMask = 0;
        size_t zMask = 0;
        size_t nrY = 0;
        if (GetPauliStringMasks(gates, xMask, zMask, nrY))
            return BaseClass::PauliExpectationValue(NrBasisStates, registerStorage, xMask, zMask, nrY,
                                                    UseMultithreading());

        VectorClass savedState = registerStorage;
        const bool recording = recordGates;
        recordGates = false;
        try
        {
            ApplyGates(gates);
            const std::complex<double> result = (savedState.adjoint() * registerStorage)(0);
            registerStorage.swap(savedState);
            recordGates = recording;
            return result;
        }
        catch (...)
        {
            registerStorage.swap(savedState);
            recordGates = recording;
            throw;
        }
    }

    std::unique_ptr<QubitRegister<VectorClass, MatrixClass>> Clone() const
    {
        auto qr = std::make_unique<QubitRegister<VectorClass, MatrixClass>>(1);
        qr->SetMultithreading(BaseClass::GetMultithreading());
        qr->NrQubits = NrQubits;
        qr->NrBasisStates = NrBasisStates;
        qr->registerStorage = registerStorage;
        qr->savedStateStorage = savedStateStorage;
        qr->computeGates = computeGates;
        qr->recordGates = recordGates;

        return qr;
    }

    MatrixClass getDensityMatrix() const
    {
        return registerStorage * registerStorage.adjoint();
    }

  protected:
    void ApplyStructuredGateImpl(const Gates::StructuredGateView &gate, const Gates::AppliedGate<MatrixClass> *recorded,
                                 bool tablesValidated, bool adjointTables = false,
                                 const Gates::StructuredRotationCoefficients *coefficients = nullptr)
    {
        if (tablesValidated)
            gate.ValidateQubits(NrQubits);
        else
            gate.Validate(NrQubits);
        // Allocate a recording before touching the state; remove it if a kernel's
        // temporary allocation fails. Compact operations remain compact in replay.
        if (recordGates)
        {
            if (recorded)
                computeGates.push_back(*recorded);
            else
                computeGates.emplace_back(gate.Materialize());
        }
        try
        {
            BaseClass::ApplyStructuredGateInPlace(registerStorage, gate, NrBasisStates, UseMultithreading(), 0, false,
                                                  coefficients, adjointTables);
        }
        catch (...)
        {
            if (recordGates)
                computeGates.pop_back();
            throw;
        }
    }

    void ReplayRecordedGate(const Gates::AppliedGate<MatrixClass> &gate, bool adjoint)
    {
        if (gate.getStructuredGate())
        {
            // The adjoint of a compact gate shares its tables.
            if (adjoint)
                ApplyGate(gate.Adjoint());
            else
                ApplyGate(gate);
        }
        else if (gate.getQubitsNumber() > 3)
        {
            if (adjoint)
                ApplyOperatorMatrix(gate.getRawOperatorMatrix().adjoint());
            else
                ApplyOperatorMatrix(gate.getRawOperatorMatrix());
        }
        else if (adjoint)
            ApplyGate(Gates::AppliedGate<MatrixClass>(gate.getRawOperatorMatrix().adjoint(), gate.getQubit1(),
                                                      gate.getQubit2(), gate.getQubit3()));
        else
            ApplyGate(gate);
    }

    bool UseMultithreading() const
    {
        return BaseClass::GetMultithreading() && NrBasisStates >= BaseClass::GetParallelMinBasisStates();
    }

    template <unsigned Qubits> double LocalNorm(const MatrixClass &matrix, const std::array<size_t, Qubits> &qubits)
    {
        constexpr Eigen::Index dimension = Eigen::Index{1} << Qubits;
        if (matrix.rows() != dimension || matrix.cols() != dimension || !matrix.allFinite())
            throw std::invalid_argument("The matrix must be finite and match its qubits");
        const auto bits = LocalBits<Qubits>(qubits);
        std::array<std::complex<double>, (size_t{1} << (2 * Qubits))> coefficients;
        for (Eigen::Index row = 0; row < dimension; ++row)
            for (Eigen::Index col = 0; col < dimension; ++col)
                coefficients[static_cast<size_t>(row * dimension + col)] = std::complex<double>(matrix(row, col));
        const double norm = BaseClass::template ApplyLocalMapNorm<Qubits>(registerStorage, coefficients.data(), bits,
                                                                          NrBasisStates, UseMultithreading());
        if (recordGates)
            computeGates.emplace_back(Gates::AppliedGate<MatrixClass>(matrix, qubits[0], Qubits > 1 ? qubits[1] : 0,
                                                                      Qubits > 2 ? qubits[2] : 0));
        return norm;
    }

    // The squared norm the state would have after the matrix acts on the qubits, with nothing written: the sum over
    // the register of |sum_c matrix(r, c) amplitude(c)|^2. A weight taken from the reduced density matrix is the
    // difference of sums over the whole register and cannot be smaller than their rounding; this one is a sum of
    // squares of what the matrix leaves of each group of amplitudes, and holds down to the rounding of the amplitudes
    // themselves. The first qubit is the least significant bit of the matrix index; the qubits were checked by
    // LocalDensity. Partial sums are taken over fixed blocks of indices and added in order, so the value does not
    // depend on the number of threads.
    template <unsigned Qubits>
    double LocalMapWeight(const MatrixClass &matrix, const std::array<size_t, Qubits> &qubits) const
    {
        constexpr size_t dimension = size_t{1} << Qubits;
        std::array<size_t, Qubits> sorted = qubits;
        std::sort(sorted.begin(), sorted.end());
        std::array<size_t, dimension> offsets{};
        for (size_t c = 0; c < dimension; ++c)
            for (unsigned b = 0; b < Qubits; ++b)
                if ((c >> b) & 1)
                    offsets[c] |= size_t{1} << qubits[b];
        std::array<std::complex<double>, dimension * dimension> coefficients;
        for (size_t row = 0; row < dimension; ++row)
            for (size_t col = 0; col < dimension; ++col)
                coefficients[row * dimension + col] =
                    std::complex<double>(matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)));

        const long long groups = static_cast<long long>(NrBasisStates >> Qubits);
        constexpr long long blockSize = 2048;
        const long long blocks = (groups + blockSize - 1) / blockSize;
        std::vector<double> partial(static_cast<size_t>(blocks), 0.);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(ParallelTeamSize(NrBasisStates)) if (UseMultithreading())
#endif
        for (long long block = 0; block < blocks; ++block)
        {
            double sum = 0.;
            const long long end = std::min(groups, (block + 1) * blockSize);
            for (long long group = block * blockSize; group < end; ++group)
            {
                // Deposit the group index into the bits that are not targets.
                size_t base = static_cast<size_t>(group);
                for (unsigned b = 0; b < Qubits; ++b)
                {
                    const size_t low = (size_t{1} << sorted[b]) - 1;
                    base = (base & low) | ((base & ~low) << 1);
                }
                std::array<std::complex<double>, dimension> amplitudes;
                for (size_t c = 0; c < dimension; ++c)
                    amplitudes[c] = std::complex<double>(registerStorage[static_cast<Eigen::Index>(base | offsets[c])]);
                for (size_t row = 0; row < dimension; ++row)
                {
                    std::complex<double> left = 0.;
                    for (size_t col = 0; col < dimension; ++col)
                        left += coefficients[row * dimension + col] * amplitudes[col];
                    sum += std::norm(left);
                }
            }
            partial[static_cast<size_t>(block)] = sum;
        }
        double weight = 0.;
        for (const double value : partial)
            weight += value;
        return weight;
    }

    // Keep the multiplication order of the density path. Squaring the coefficient first could overflow for an
    // incomplete Kraus family on an unnormalized state even when its branch weight is representable.
    template <size_t Dimension>
    static void KrausPopulationWeights(const std::vector<MatrixClass> &kraus, const std::array<double, Dimension> &p,
                                       double *weights)
    {
        for (size_t k = 0; k < kraus.size(); ++k)
        {
            double weight = 0.;
            for (size_t row = 0; row < Dimension; ++row)
                for (size_t col = 0; col < Dimension; ++col)
                {
                    const std::complex<double> value(kraus[k](row, col));
                    if (value != std::complex<double>{})
                    {
                        weight += std::real((value * p[col]) * std::conj(value));
                        break;
                    }
                }
            weights[k] = weight;
        }
    }

    template <unsigned Qubits>
    std::array<size_t, Qubits> LocalBits(const std::array<size_t, Qubits> &qubits) const
    {
        std::array<size_t, Qubits> bits{};
        for (unsigned k = 0; k < Qubits; ++k)
        {
            if (qubits[k] >= NrQubits)
                throw std::invalid_argument("Qubit number is too high");
            for (unsigned j = 0; j < k; ++j)
                if (qubits[j] == qubits[k])
                    throw std::invalid_argument("Qubits must be different");
            bits[k] = size_t{1} << qubits[k];
        }
        return bits;
    }

    template <unsigned Qubits>
    std::array<double, (size_t{1} << Qubits)> Populations(const std::array<size_t, Qubits> &qubits) const
    {
        return BaseClass::template LocalPopulations<Qubits>(registerStorage, LocalBits<Qubits>(qubits),
                                                             NrBasisStates, UseMultithreading());
    }

    template <unsigned Qubits, bool RealPartOnly = false>
    MatrixClass LocalDensity(const std::array<size_t, Qubits> &qubits) const
    {
        const auto bits = LocalBits<Qubits>(qubits);
        const auto rho = BaseClass::template LocalDensityMatrix<Qubits, RealPartOnly>(
            registerStorage, bits, NrBasisStates, UseMultithreading());
        constexpr Eigen::Index dimension = Eigen::Index{1} << Qubits;
        MatrixClass result(dimension, dimension);
        for (Eigen::Index row = 0; row < dimension; ++row)
            for (Eigen::Index col = 0; col < dimension; ++col)
                result(row, col) = typename MatrixClass::Scalar(rho[static_cast<size_t>(row * dimension + col)]);
        return result;
    }

    inline void CheckQubits(const GateClass & /*gate*/, size_t qubit, size_t controllingQubit1,
                            size_t controllingQubit2, size_t gateQubits) const
    {
        if (NrQubits == 0)
            throw std::invalid_argument("Qubit number is zero");
        else if (NrQubits <= qubit)
            throw std::invalid_argument("Qubit number is too high");
        else if (gateQubits == 2)
        {
            if (NrQubits <= controllingQubit1)
                throw std::invalid_argument("Controlling qubit number is too high");
            else if (qubit == controllingQubit1)
                throw std::invalid_argument("Qubit and controlling qubit are the same");
        }
        else if (gateQubits == 3)
        {
            if (NrQubits <= controllingQubit1 || NrQubits <= controllingQubit2)
                throw std::invalid_argument("Controlling qubit number is too high");
            else if (qubit == controllingQubit1 || qubit == controllingQubit2 || controllingQubit1 == controllingQubit2)
                throw std::invalid_argument("Qubits must be different");
        }
    }

    // Recognizes a product of single qubit X, Y, Z (or identity) gates, each on a different qubit.
    // Exact comparisons only, anything else (including a repeated qubit) is left to the generic path.
    bool GetPauliStringMasks(const std::vector<Gates::AppliedGate<MatrixClass>> &gates, size_t &xMask, size_t &zMask,
                             size_t &nrY) const
    {
        const std::complex<double> zero(0., 0.);
        const std::complex<double> one(1., 0.);
        const std::complex<double> i(0., 1.);

        size_t usedQubits = 0;
        for (const auto &gate : gates)
        {
            if (gate.getStructuredGate())
                return false;
            const MatrixClass &m = gate.getRawOperatorMatrix();
            if (m.rows() != 2 || m.cols() != 2)
                return false;

            const size_t qubit = gate.getQubit1();
            if (qubit >= NrQubits)
                return false;

            const size_t qubitBit = 1ULL << qubit;
            if (usedQubits & qubitBit)
                return false;
            usedQubits |= qubitBit;

            if (m(0, 1) == zero && m(1, 0) == zero && m(0, 0) == one)
            {
                if (m(1, 1) == one)
                    continue;        // identity
                if (m(1, 1) == -one) // Z
                {
                    zMask |= qubitBit;
                    continue;
                }
            }
            else if (m(0, 0) == zero && m(1, 1) == zero)
            {
                if (m(0, 1) == one && m(1, 0) == one) // X
                {
                    xMask |= qubitBit;
                    continue;
                }
                if (m(0, 1) == -i && m(1, 0) == i) // Y
                {
                    xMask |= qubitBit;
                    zMask |= qubitBit;
                    ++nrY;
                    continue;
                }
            }

            return false;
        }

        return true;
    }

    // shortcut for measuring a single qubit
    size_t MeasureNoCollapse(size_t qubit)
    {
        return MeasureNoCollapse(qubit, qubit);
    }

    // measure a 'subregister' as a separate register
    // can measure a single qubit, if firstQubit == secondQubit
    // will return a 'state' as if the measured sequence is in a separate register (that is, the 'firstQubit' is on
    // position 0 and so on) so 0 means that all measured qubits are zero, 1 means that firstQubit is 1 and all other
    // measured ones are zero, 2 means that the next one 1 one and all others are zero and so on

    size_t MeasureNoCollapse(size_t firstQubit, size_t secondQubit)
    {
        const double prob = 1. - uniformZeroOne(rng); // this excludes 0 as probabiliy

        if (firstQubit == secondQubit)
            return BaseClass::MeasureQubitNoCollapse(NrBasisStates, registerStorage, firstQubit, prob);

        return BaseClass::MeasureNoCollapse(NrBasisStates, registerStorage, firstQubit, secondQubit, prob);
    }

    size_t NrQubits;
    size_t NrBasisStates;

    VectorClass registerStorage;

    VectorClass savedStateStorage;

    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniformZeroOne;

    std::vector<Gates::AppliedGate<MatrixClass>> computeGates;

    bool recordGates;

    // The weight tr(E rho E^dagger) of a Kraus operator is in doubt when it is below this fraction of its reach, the
    // bound on the terms it sums (see ApplyKrausOperator): it is then what their cancellation left, and its relative
    // error is the rounding of rho over this fraction. The usual operators of a noise channel (diagonal, a single
    // entry per column, a multiple of a unitary) have a weight of half or a quarter of their reach on any state.
    static constexpr double UnresolvedWeight = 1e-3;
    // The relative rounding assumed in the entries of the reduced density matrix, with room for a large register:
    // times the reach of an operator, how far the weight computed from it can be from the true one.
    static constexpr double DensityRounding = 1e-12;
    // The relative rounding of an amplitude after a local map. A weight read from the amplitudes that is below its
    // square times the reach is the rounding of amplitudes that cancel exactly, and counts as none.
    static constexpr double AmplitudeRounding = 64. * 2.220446049250313e-16;
};

} // namespace QC
