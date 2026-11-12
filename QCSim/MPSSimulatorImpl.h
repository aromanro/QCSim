#pragma once

#include "MPSSimulatorBase.h"
#include "TensorPauli.h"

#include <iomanip>
#include <sstream>

#define USE_FAST_SVD 1

namespace QC
{

namespace TensorNetworks
{

class MPSSimulatorImpl : public MPSSimulatorBase
{
  public:
    friend class MPSSimulator;

    MPSSimulatorImpl(size_t N, unsigned int addseed = 0) : MPSSimulatorBase(N, addseed)
    {
#ifdef USE_FAST_SVD
        // 16 is Eigen's default. With the Vidal form (dividing by the lambdas) it caused precision issues in tests
        // against statevector, so it used to be 64, but with Hastings' method it's as precise as Jacobi and much faster
        // for bond dimensions between 8 and 32
        SVD.setSwitchSize(blockSizeLimit); // lower sizes will use Jacobi
#endif
    }

    void ApplyGate(const Gates::AppliedGate<MatrixClass> &gate) override
    {
        ApplyGate(gate, gate.getQubit1(), gate.getQubit2());
    }

    // three qubit gates not supported, convert them to two qubit gates
    // also two qubit gates need to act on adjacent qubits
    // don't try to apply a gate that doesn't satisfy these conditions
    // use swap gates to move qubits around
    // this is wrapped up into a higher level simulator that swaps the qubits for you and maps them to minimize swaps
    // so this should not be used directly
    void ApplyGate(const GateClass &gate, IndexType qubit, IndexType controllingQubit1 = 0) override
    {
        if (gate.getQubitsNumber() > 2)
            throw std::invalid_argument("Three qubit gates not supported");
        else if ((qubit < 0 || qubit >= static_cast<IndexType>(gammas.size())) ||
                 (gate.getQubitsNumber() == 2 &&
                  (controllingQubit1 < 0 || controllingQubit1 >= static_cast<IndexType>(gammas.size()))))
            throw std::invalid_argument("Qubit index out of bounds");
        else if (gate.getQubitsNumber() == 2 &&
                 std::abs(static_cast<int>(qubit) - static_cast<int>(controllingQubit1)) != 1)
            throw std::invalid_argument("Two qubit gates need to act on adjacent qubits");
        else if (!gate.getRawOperatorMatrix().allFinite())
            throw std::invalid_argument("Gate matrix contains a non-finite value");

        // Environments kept only for reads (see EnvironmentsForRead) end with them: a gate goes the way it always
        // went on a chain that holds no Kraus operator.
        if (gaugeForReads)
            LeaveGaugeMode();

        if (gate.getQubitsNumber() == 1)
            ApplySingleQubitGate(gate, qubit);
        else
            ApplyTwoQubitGate(gate, qubit, controllingQubit1);
    }

    void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>> &gates) override
    {
        for (const auto &gate : gates)
            ApplyGate(gate);
    }

    // Walks over the chain and, wherever a bond dimension is bigger than the currently set
    // bond dimension limit, contracts the two neighbour sites and re-splits them with a
    // truncated SVD to bring the bond dimension down to the limit. No gate is applied.
    // Useful after lowering the bond dimension limit with setLimitBondDimension once the
    // bonds have already grown beyond the new limit.
    void Trim() override
    {
        if (!limitSize)
            return; // nothing to trim against
        EnsureCanonical();

        for (IndexType qubit1 = 0; qubit1 < static_cast<IndexType>(lambdas.size()); ++qubit1)
        {
            if (lambdas[qubit1].size() <= chi)
                continue;

            // contract the two neighbour sites (no gate is applied), then re-split with a truncated SVD
            const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);
            const MatrixClass thetaMatrix = ReshapeTheta(theta);

            DecomposeAndSetGammas(thetaMatrix, qubit1, qubit1 + 1);
        }
    }

    void ReCanonicalize() override
    {
        RestoreCanonicalForm(true);
    }

    // The body of ReCanonicalize. Without `compress` the state is left as it is: the sweeps only drop the directions
    // the decompositions find without weight, and none of the cuts the user asked for (bond limit, singular value
    // threshold) is applied. That is what a read needs, see EnsureCanonical.
    void RestoreCanonicalForm(bool compress)
    {
        // The sweeps below rebuild everything from the tensors: they need no environment and leave none stale.
        LeaveGaugeMode();
        if (lambdas.empty())
        {
            canonicalFormExact = true;
            return;
        }

        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            bool applyUserCompression = compress;
            bool reduced;
            do
            {
                const auto previousDimensions = getBondDimensions();
                // Repair the basis without trusting old Schmidt weights or dropping directions.
                for (IndexType site = static_cast<IndexType>(gammas.size()) - 1; site > 0; --site)
                    RightCanonicalizeSite(site);

                // The right block is now orthonormal. Starting at bond zero makes each
                // newly computed lambda the valid left environment of the following bond.
                for (IndexType bond = 0; bond < static_cast<IndexType>(lambdas.size()); ++bond)
                {
                    const MatrixClass theta = ReshapeTheta(ContractTwoQubits(bond));
                    DecomposeAndSetGammasImpl(theta, bond, bond + 1, applyUserCompression, lambdas[bond].size());
                }

                // Truncation can stale earlier Schmidt spectra. Refresh until a pass
                // changes no dimensions; no bond can grow, so this terminates.
                reduced = getBondDimensions() != previousDimensions;
                // A relative cutoff must also hold for the refreshed spectra. Discarded
                // weight is a per-compression budget: do not spend it again on refresh.
                applyUserCompression = compress && truncationMode == TruncationMode::RelativeToMax;
            } while (reduced);
        });
        // The last pass changed no dimension, so it discarded nothing: every lambda is the Schmidt spectrum of its
        // bond and every tensor to the right of a bond is orthonormal.
        canonicalFormExact = true;
    }

    // false if measured 0, true if measured 1
    bool MeasureQubit(IndexType qubit) override
    {
        if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");

        const double rndVal = 1. - uniformZeroOne(rng);
        // On a chain that holds Kraus operators a measurement is one more of them: the channel of the two
        // projectors, with the same rule for the draw (0 when the number drawn is below its probability). The site
        // is split again where the projector left it (see RegaugeSite); nothing sweeps the chain.
        if (gaugeMode && !gaugeForReads)
            return ApplyKrausOperator({zeroProjection.getRawOperatorMatrix(), oneProjection.getRawOperatorMatrix()},
                                      rndVal, qubit) == 1;

        // The probability is taken from the environments of the site when the canonical form is approximate (see
        // GetProbability); the collapse and what follows it stay on the canonical path.
        const double prob0 = ValidMeasurementProbability(GetProbability(qubit));
        if (gaugeForReads)
            LeaveGaugeMode();
        const bool zeroMeasured = PickZeroMeasurement(rndVal, prob0);

        MatrixClass projMat = zeroMeasured ? zeroProjection.getRawOperatorMatrix() * 1. / sqrt(prob0)
                                           : oneProjection.getRawOperatorMatrix() * 1. / sqrt(1. - prob0);
        const QC::Gates::SingleQubitGate projOp(std::move(projMat));

        ApplySingleQubitGate(projOp, qubit);

        // propagate to the other qubits to the right and left until the end or there is no entanglement with the next
        // one

        for (IndexType q = qubit; q < static_cast<IndexType>(lambdas.size()); ++q)
        {
            if (lambdas[q].size() == 1)
                break;
            ApplyTwoQubitGate(projOp, q, q + 1, true);
        }

        for (IndexType q = qubit; q > 0; --q)
        {
            const IndexType q1 = q - 1;
            if (lambdas[q1].size() == 1)
                break;
            ApplyTwoQubitGate(projOp, q1, q, true);
        }

        return !zeroMeasured;
    }

    std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType> &qubits) override
    {
        if (qubits.empty())
            return {};

        std::unordered_map<IndexType, bool> res;
        if (qubits.size() == 1)
        {
            const IndexType qubit = *qubits.begin();
            res[qubit] = MeasureQubit(qubit);
            return res;
        }
        if (*qubits.cbegin() < 0 || *qubits.crbegin() >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        // On a chain that holds Kraus operators, one by one, each as the channel of its projectors (see MeasureQubit).
        if (gaugeMode && !gaugeForReads)
        {
            for (const auto qubit : qubits)
                res[qubit] = MeasureQubit(qubit);
            return res;
        }

        const auto lastIt = --qubits.end();

        for (auto it = qubits.begin(); it != lastIt;)
        {
            const auto qubit = *it;
            ++it;
            const auto nextQubit = *it;

            const double rndVal = 1. - uniformZeroOne(rng);
            const double prob0 = ValidMeasurementProbability(GetProbability(qubit));
            if (gaugeForReads)
                LeaveGaugeMode();
            const bool zeroMeasured = PickZeroMeasurement(rndVal, prob0);

            res[qubit] = !zeroMeasured;

            MatrixClass projMat = zeroMeasured ? zeroProjection.getRawOperatorMatrix() * 1. / sqrt(prob0)
                                               : oneProjection.getRawOperatorMatrix() * 1. / sqrt(1. - prob0);
            const QC::Gates::SingleQubitGate projOp(std::move(projMat));

            ApplySingleQubitGate(projOp, qubit);

            for (IndexType q = qubit; q < nextQubit; ++q)
            {
                if (lambdas[q].size() == 1)
                    break;
                ApplyTwoQubitGate(projOp, q, q + 1, true);
            }

            for (IndexType q = qubit; q > 0; --q)
            {
                const IndexType q1 = q - 1;
                if (lambdas[q1].size() == 1)
                    break;
                ApplyTwoQubitGate(projOp, q1, q, true);
            }
        }

        const auto lastQubit = *lastIt;
        res[lastQubit] = MeasureQubit(lastQubit);

        return res;
    }

    // One-site operator products may be non-Hermitian, so the result is complex.
    // Only the affected physical interval is contracted.
    std::complex<double> ExpectationValue(const std::vector<Gates::AppliedGate<MatrixClass>> &gates) override
    {
        if (gates.empty())
            return 1.;
        // The environments it needs, when the canonical form does not give them, are set up where the interval of
        // the operators is known (see ExpectationValueImpl).
        std::complex<double> result;
        RunMaybeSingleThreaded(enableMultithreading, [&]() { result = ExpectationValueImpl(gates); });
        return result;
    }

    std::vector<std::complex<double>> ExpectationValues(const std::vector<std::string> &paulis) const
    {
        if (paulis.empty())
            return {};
        std::vector<std::complex<double>> result;
        RunMaybeSingleThreaded(enableMultithreading, [&]() { result = ExpectationValuesImpl(paulis); });
        return result;
    }

    // Sampling needs the block to the right of every site. On an exact canonical form it is orthonormal; otherwise
    // its environments are computed once and serve every sample until the chain changes (see EnvironmentsForRead).
    std::unordered_map<IndexType, bool> MeasureNoCollapse() override
    {
        if (gammas.empty())
            return {};
        PrepareSampling();
        const size_t nrQubits = gammas.size();
        return MeasureNoCollapse(static_cast<IndexType>(nrQubits - 1));
    }

    // Brings every right environment up to date, where the chain keeps or needs them, inside the threading guard:
    // the samples then only read them.
    void PrepareSampling()
    {
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            EnvironmentsForRead(0);
            if (gaugeMode && !lambdas.empty())
                RightEnvironment(0);
        });
    }

    std::unordered_map<IndexType, bool> MeasureNoCollapse(const std::set<IndexType> &qubits) override
    {
        if (qubits.empty())
            return {};
        if (*qubits.cbegin() < 0 || *qubits.crbegin() >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        PrepareSampling();

        return MeasureNoCollapse(*qubits.crbegin());
    }

    // The base class answers from the canonical form: the squared lambdas to the left of the site, orthonormal
    // tensors to its right. When the chain keeps explicit environments, or when gates have left that form
    // approximate by discarding Schmidt values, the probability is taken from the environments of the site instead.
    double GetProbability(IndexType qubit, bool zeroVal = true) const override
    {
        if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        if (!gaugeMode && canonicalFormExact)
            return MPSSimulatorBase::GetProbability(qubit, zeroVal);

        MPSSimulatorImpl *const self = const_cast<MPSSimulatorImpl *>(this);
        Eigen::Matrix2cd rho;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            self->EnvironmentsForRead(qubit);
            rho = self->SiteDensity(qubit);
        });
        return ClampProbability(rho(zeroVal ? 0 : 1, zeroVal ? 0 : 1).real());
    }

    // A saved state holds no environments: it is taken on the canonical path.
    std::shared_ptr<MPSSimulatorStateInterface> getState() const override
    {
        EnsureCanonical();
        return MPSSimulatorBase::getState();
    }

    void setState(const std::shared_ptr<MPSSimulatorStateInterface> &state) override
    {
        if (state)
            LeaveGaugeMode();
        MPSSimulatorBase::setState(state);
    }

    void setStateDestructive(std::shared_ptr<MPSSimulatorStateInterface> &state) override
    {
        if (state)
            LeaveGaugeMode();
        MPSSimulatorBase::setStateDestructive(state);
    }

    void Clear() override
    {
        LeaveGaugeMode();
        MPSSimulatorBase::Clear();
    }

    void InitOnesState() override
    {
        LeaveGaugeMode();
        MPSSimulatorBase::InitOnesState();
    }

    // setToQubitState and setToBasisState are not overridden: the base class resets the chain through Clear(), which
    // ends the gauge mode above, and only once it has accepted its argument. A basis state that does not fit the
    // register is refused (or ignored) before anything changes, and the chain then keeps its environments too.

    // -----------------------------------------------------------------------------------------------------------------
    // Quantum trajectories: one Kraus operator of a channel, drawn with the probability the state gives it.
    //
    // A non-unitary one-site operator takes the chain out of the canonical form every other operation here relies on:
    // the tensor it acts on is no longer right-orthonormal, and the stored lambdas no longer describe the blocks to
    // the left of the bonds on its right. Restoring the form costs a sweep of decompositions over the chain, as a
    // measurement does, and a noisy circuit applies such an operator after every gate.
    //
    // Instead, from the first such operator on, the chain keeps its environments explicitly ('gauge mode'): for each
    // bond the Gram matrix of the block to its left and of the block to its right, cached and recomputed only across
    // the sites that changed since they were last used. Then
    //  - the reduced density matrix of a site, which gives the probability of every Kraus operator, costs a few
    //    matrix products, and applying the operator drawn leaves all the other tensors as they are;
    //  - a two-qubit gate weighs its decomposition by the two environments instead of the left lambda, so that its
    //    singular values are the Schmidt values of the bond and cutting them is the optimal truncation, exactly as
    //    in the canonical form (see DecomposeWithEnvironments);
    //  - a unitary one-qubit gate changes no environment.
    // No operation sweeps the chain. Whatever assumes the canonical form (measurement, expectation values, saved
    // states) restores it first with ReCanonicalize, which ends the gauge mode.

    // Reduced density matrix of one site: rho(s, s') = <s| rho |s'>. The state is not changed. On a chain whose gates
    // have discarded Schmidt values the canonical form no longer gives the environments of the site: they are
    // computed from the tensors and kept for the reads that follow, as for the other reads (see
    // EnvironmentsForRead), so the gates after this one are those of a chain that was not read.
    MatrixClass LocalDensityMatrix(IndexType qubit)
    {
        if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        Eigen::Matrix2cd rho;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            EnvironmentsForRead(qubit);
            rho = SiteDensity(qubit);
        });
        return rho;
    }

    // One step of a trajectory through the channel rho -> sum_k E_k rho E_k^dagger on one site: applies one of the
    // 2x2 operators, chosen with the probability the state gives it, and renormalizes. `random` is uniform in
    // [0, 1). Trace preservation is not checked: the weights are taken relative to their sum. Returns the index of
    // the operator applied.
    //
    // Precision, as for the statevector (QubitRegister::ApplyKrausOperator). A weight that is a small part of what
    // its terms add up to is read again from the tensors, on the canonical form, and one below (64 epsilon)^2
    // times the bound on those terms counts as none: it is what the rounding of amplitudes that cancel leaves.
    // Such an operator is never drawn, and when none is left the call throws std::domain_error with the largest
    // weight found and the bound it fell below. This refuses, too, a branch that small which is exact.
    size_t ApplyKrausOperator(const std::vector<MatrixClass> &kraus, double random, IndexType qubit)
    {
        if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        if (kraus.empty())
            throw std::invalid_argument("A channel must contain at least one Kraus operator");
        for (const auto &e : kraus)
            if (e.rows() != 2 || e.cols() != 2 || !e.allFinite())
                throw std::invalid_argument("Matrix product state Kraus operators must be finite 2x2 matrices");

        size_t chosen = 0;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            // The canonical form gives the environments of the site only while it is exact (see LocalDensityMatrix).
            if (!gaugeMode && !canonicalFormExact)
                EnterGaugeModeNormalized(qubit);
            const Eigen::Matrix2cd rho = SiteDensity(qubit);
            // Stack storage for the usual few operators: the weight of each, and its reach when it is in doubt.
            double fewWeights[16], fewReaches[16];
            std::vector<double> manyValues(kraus.size() > 16 ? 2 * kraus.size() : 0);
            double *const weights = kraus.size() > 16 ? manyValues.data() : fewWeights;
            double *const reaches = kraus.size() > 16 ? manyValues.data() + kraus.size() : fewReaches;
            double total = 0., doubt = 0.;
            for (size_t k = 0; k < kraus.size(); ++k)
            {
                const Eigen::Matrix2cd e = kraus[k];
                const double weight = (e * rho * e.adjoint()).trace().real();
                // The bound on the terms of that trace, 2 sum_a |column a of E|^2 rho(a, a): see UnresolvedWeight.
                const double reach =
                    2. * (e.col(0).squaredNorm() * rho(0, 0).real() + e.col(1).squaredNorm() * rho(1, 1).real());
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

            // Whether the weights were read on the canonical form, which the draw below then leaves to be applied.
            bool canonical = false;
            if (doubt == 0.)
            {
                if (!std::isfinite(total) || total <= 0.)
                    throw std::domain_error("No Kraus operator has any probability on this state");
                DrawOperator(weights, kraus.size(), total, 0., random, chosen);
            }
            else if (!(std::isfinite(total) && total > 0. &&
                       DrawOperator(weights, kraus.size(), total, doubt, random, chosen) && reaches[chosen] < 0.))
            {
                // A weight in doubt decides the draw, or is the one drawn. The environments are Gram matrices: they
                // hold a weight to the rounding of their largest entry, and the norm of what a nearly orthogonal
                // operator leaves is below that. The canonical form is restored instead, where the weights are
                // sums of squares of the tensor of the site (see CanonicalSiteWeight); they are all read there
                // before anything is drawn, so that a weight that is there is found whatever the rounding made of
                // it, and one that is not there is not drawn.
                RestoreCanonicalForm(false);
                total = 0.;
                double largest = 0., cutoffOfLargest = 0.;
                for (size_t k = 0; k < kraus.size(); ++k)
                {
                    double reach = 0.;
                    const double weight = CanonicalSiteWeight(qubit, kraus[k], reach);
                    const double cutoff = AmplitudeRounding * AmplitudeRounding * reach;
                    if (weight > largest)
                    {
                        largest = weight;
                        cutoffOfLargest = cutoff;
                    }
                    weights[k] = weight >= cutoff ? weight : 0.;
                    total += weights[k];
                }
                if (!std::isfinite(total) || total <= 0.)
                    ThrowNoProbability(largest, cutoffOfLargest);
                DrawOperator(weights, kraus.size(), total, 0., random, chosen);
                canonical = true;
            }

            const MatrixClass scaled = kraus[chosen] / std::sqrt(weights[chosen]);
            const Eigen::Matrix2cd drawn = kraus[chosen];
            // Read from the density matrix the weights came from; after the weights were read again on the canonical
            // form the branch below is taken whatever the operator.
            const bool rankOne = !canonical && LeavesOneState(drawn, rho, weights[chosen]);
            if (canonical || (!rankOne && ProjectionLike(drawn)))
            {
                // A near-null operator, of which what was a small part of the state is now most, or one that nearly
                // projects without leaving the site in one state: the lambdas know a small weight only to the
                // rounding of the largest, and explicit environments would go on from Gram matrices that have all
                // but lost their rank. The canonical form is rebuilt from the tensors, which also normalizes the
                // state from them; the scale applied here only keeps the numbers in range. Both are rare: the first
                // needs a state nearly orthogonal to the operator, the second a damping of all but a hundredth of
                // an amplitude.
                ApplySingleQubitGate(gammas[qubit], scaled);
                NormalizeLoneSite();
                RestoreCanonicalForm(false);
                return;
            }
            // The environments are taken from the canonical form the first time, before the site changes. Those
            // kept for reads serve as well: from here on the chain holds a Kraus operator.
            EnterGaugeMode();
            gaugeForReads = false;
            ApplySingleQubitGate(gammas[qubit], scaled);
            // Every environment that contains the site is stale: the left blocks of the bonds from this site on,
            // and the right blocks of the bonds before it.
            leftValid = std::min(leftValid, qubit);
            rightValid = std::max(rightValid, qubit);
            // A projection (a jump, a reset, one outcome of a measurement) leaves the site in one state, whatever
            // the rest holds: the site is split again where the projection left it.
            if (rankOne)
                RegaugeSite(qubit, drawn);
        });
        return chosen;
    }

    // The same step for 4x4 operators on two neighbouring sites, `qubit` the low bit of the matrix index and
    // `controllingQubit1` the high one. The two-site tensor is built once. With the environments around it, it
    // gives the reduced density matrix of the pair, and so the probability of every operator; the operator drawn
    // is applied to that same tensor, which is then split in one decomposition weighed by the environments, as
    // for a two-qubit gate. The step costs what a gate costs, whatever the number of operators. The precision of
    // the weights, and what is refused for want of one, are those of the one-site step.
    size_t ApplyKrausOperator(const std::vector<MatrixClass> &kraus, double random, IndexType qubit,
                              IndexType controllingQubit1)
    {
        const IndexType sites = static_cast<IndexType>(gammas.size());
        if (qubit < 0 || qubit >= sites || controllingQubit1 < 0 || controllingQubit1 >= sites)
            throw std::invalid_argument("Qubit index out of bounds");
        if (std::abs(qubit - controllingQubit1) != 1)
            throw std::invalid_argument("Two-qubit Kraus operators need to act on adjacent qubits");
        if (kraus.empty())
            throw std::invalid_argument("A channel must contain at least one Kraus operator");
        for (const auto &e : kraus)
            if (e.rows() != 4 || e.cols() != 4 || !e.allFinite())
                throw std::invalid_argument("Two-qubit Kraus operators must be finite 4x4 matrices");

        const IndexType qubit1 = std::min(qubit, controllingQubit1), qubit2 = qubit1 + 1;
        // The left site is the high bit of the operators' index, unless it is `qubit`.
        const bool reversed = qubit < controllingQubit1;
        size_t chosen = 0;
        double weight = 0.;
        bool decomposed = false, projects = false;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            // The environments are taken from the canonical form the first time, before the sites change. Those
            // kept for reads serve as well: from here on the chain holds a Kraus operator.
            EnterGaugeMode();
            gaugeForReads = false;
            const IndexType bonds = static_cast<IndexType>(lambdas.size());
            // theta(j * szl + l, k * szr + r), with j and k the states of the left and of the right site
            const MatrixClass theta = ReshapeTheta(ContractTwoQubits(qubit1));
            const IndexType szl = theta.rows() / 2, szr = theta.cols() / 2;
            MatrixClass weighted = theta;
            if (qubit1 > 0)
            {
                const MatrixClass &left = LeftEnvironment(qubit1 - 1);
                for (IndexType j = 0; j < 2; ++j)
                    weighted.middleRows(j * szl, szl).noalias() = left * theta.middleRows(j * szl, szl);
            }
            if (qubit2 < bonds)
            {
                const MatrixClass &right = RightEnvironment(qubit2);
                for (IndexType k = 0; k < 2; ++k)
                {
                    const MatrixClass block = weighted.middleCols(k * szr, szr) * right.transpose();
                    weighted.middleCols(k * szr, szr) = block;
                }
            }
            // Where the block of one state of the pair starts, the state given in the index of the operators.
            const auto row = [&](int index) { return (reversed ? index % 2 : index / 2) * szl; };
            const auto col = [&](int index) { return (reversed ? index / 2 : index % 2) * szr; };

            // rho(a, b) = <a| rho |b> = sum conj(theta_b) (G theta_a H^T), with G and H the environments
            Eigen::Matrix4cd rho;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    rho(a, b) = theta.block(row(b), col(b), szl, szr)
                                    .conjugate()
                                    .cwiseProduct(weighted.block(row(a), col(a), szl, szr))
                                    .sum();

            // The weight of each operator, and its reach when the weight is in doubt (see UnresolvedWeight).
            std::vector<double> weights(kraus.size()), reaches(kraus.size());
            double total = 0., doubt = 0.;
            for (size_t k = 0; k < kraus.size(); ++k)
            {
                const Eigen::Matrix4cd e = kraus[k];
                const double w = (e * rho * e.adjoint()).trace().real();
                double reach = 0.;
                for (int a = 0; a < 4; ++a)
                    reach += e.col(a).squaredNorm() * rho(a, a).real();
                reach *= 4.;
                weights[k] = w > 0. ? w : 0.;
                if (w >= UnresolvedWeight * reach)
                    reaches[k] = -1.; // resolved
                else
                {
                    reaches[k] = reach;
                    doubt += DensityRounding * reach;
                }
                total += weights[k];
            }

            if (doubt == 0.)
            {
                if (!std::isfinite(total) || total <= 0.)
                    throw std::domain_error("No Kraus operator has any probability on this state");
                DrawOperator(weights.data(), kraus.size(), total, 0., random, chosen);
            }
            else if (!(std::isfinite(total) && total > 0. &&
                       DrawOperator(weights.data(), kraus.size(), total, doubt, random, chosen) &&
                       reaches[chosen] < 0.))
            {
                // A weight in doubt decides the draw, or is the one drawn: as for one site, the canonical form is
                // restored and every weight read there, as a sum of squares of what the operator leaves of the
                // two-site tensor, before anything is drawn. The operator drawn is applied as it is, and the
                // canonical decomposition normalizes by the singular values of what it left.
                RestoreCanonicalForm(false);
                const MatrixClass exact = ReshapeTheta(ContractTwoQubits(qubit1));
                const IndexType exactL = exact.rows() / 2, exactR = exact.cols() / 2;
                const auto exactRow = [&](int index) { return (reversed ? index % 2 : index / 2) * exactL; };
                const auto exactCol = [&](int index) { return (reversed ? index / 2 : index % 2) * exactR; };
                // The squared norm of a matrix of the shape of theta, its rows weighed by the left lambdas.
                const auto squaredNorm = [&](const MatrixClass &matrix) {
                    double sum = 0.;
                    for (IndexType j = 0; j < 2; ++j)
                        for (IndexType l = 0; l < exactL; ++l)
                        {
                            const double lambda = qubit1 == 0 ? 1. : lambdas[qubit1 - 1][l];
                            sum += lambda * lambda * matrix.row(j * exactL + l).squaredNorm();
                        }
                    return sum;
                };
                const auto applyTo = [&](const Eigen::Matrix4cd &op) {
                    MatrixClass applied = MatrixClass::Zero(exact.rows(), exact.cols());
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b)
                            if (op(a, b) != std::complex<double>{})
                                applied.block(exactRow(a), exactCol(a), exactL, exactR) +=
                                    op(a, b) * exact.block(exactRow(b), exactCol(b), exactL, exactR);
                    return applied;
                };
                double populations[4];
                for (int a = 0; a < 4; ++a)
                {
                    MatrixClass part = MatrixClass::Zero(exact.rows(), exact.cols());
                    part.block(exactRow(a), exactCol(a), exactL, exactR) =
                        exact.block(exactRow(a), exactCol(a), exactL, exactR);
                    populations[a] = squaredNorm(part);
                }
                total = 0.;
                double largest = 0., cutoffOfLargest = 0.;
                for (size_t k = 0; k < kraus.size(); ++k)
                {
                    const Eigen::Matrix4cd e = kraus[k];
                    double reach = 0.;
                    for (int a = 0; a < 4; ++a)
                        reach += e.col(a).squaredNorm() * populations[a];
                    reach *= 4.;
                    const double w = squaredNorm(applyTo(e));
                    const double cutoff = AmplitudeRounding * AmplitudeRounding * reach;
                    if (w > largest)
                    {
                        largest = w;
                        cutoffOfLargest = cutoff;
                    }
                    weights[k] = w >= cutoff ? w : 0.;
                    total += weights[k];
                }
                if (!std::isfinite(total) || total <= 0.)
                    ThrowNoProbability(largest, cutoffOfLargest);
                DrawOperator(weights.data(), kraus.size(), total, 0., random, chosen);
                const Eigen::Matrix4cd drawn = kraus[chosen] / std::sqrt(weights[chosen]);
                DecomposeAndSetGammasImpl(applyTo(drawn), qubit1, qubit2, true);
                // What was a small part of the state may now be most of it, and the lambdas knew it only to the
                // rounding of the largest: the form is rebuilt once more, which normalizes from the tensors.
                RestoreCanonicalForm(false);
                decomposed = true;
                return;
            }
            weight = weights[chosen];

            const Eigen::Matrix4cd op = kraus[chosen] / std::sqrt(weight);
            MatrixClass applied = MatrixClass::Zero(theta.rows(), theta.cols());
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    if (op(a, b) != std::complex<double>{})
                        applied.block(row(a), col(a), szl, szr) += op(a, b) * theta.block(row(b), col(b), szl, szr);
            decomposed = DecomposeWithEnvironments(applied, qubit1, qubit2, (op.adjoint() * op).isIdentity(1e-10));
            projects = ProjectionLike(Eigen::Matrix4cd(kraus[chosen]));
        });
        if (!decomposed)
        {
            // An environment has lost rank: the gate path restores the canonical form for this update.
            const Gates::AppliedGate<MatrixClass> gate(MatrixClass(kraus[chosen] / std::sqrt(weight)), qubit,
                                                       controllingQubit1);
            ApplyTwoQubitGate(gate, qubit, controllingQubit1);
        }
        // As for one site: after an operator that nearly projects, the canonical form is rebuilt from the tensors.
        if (projects)
            RestoreCanonicalForm(false);
        return chosen;
    }

    // Whether the chain keeps explicit environments instead of the canonical form (see above).
    bool InGaugeMode() const
    {
        return gaugeMode;
    }

    // Computes <0|Psi> by contracting the current MPS chain with the zero state
    // useful for maestro/composer for checking some stuff
    std::complex<double> ProjectOnZero() const override
    {
        const size_t nrQubits = getNrQubits();
        if (nrQubits == 0)
            return 0.;

        static const Indexes product_dims{IntIndexPair(1, 0)};
        MatrixTensorType res = gammas[0].chip(0, 1);

        // the lambdas are already included in the B tensors
        for (size_t q = 1; q < nrQubits; ++q)
        {
            // why? Needs this intermediary variable here, not even calling eval() works if assigning directly to res
            MatrixTensorType tmp = res.contract(gammas[q].chip(0, 1), product_dims);
            res = std::move(tmp);
        }

        return res(0, 0);
    }

  private:
    std::complex<double> ExpectationValueImpl(const std::vector<Gates::AppliedGate<MatrixClass>> &gates)
    {
        if (gates.empty())
            return 1.;
        IndexType first = static_cast<IndexType>(gammas.size()), last = 0;
        for (const auto &gate : gates)
        {
            const IndexType q = gate.getQubit1();
            if (q < 0 || q >= static_cast<IndexType>(gammas.size()))
                throw std::invalid_argument("Qubit index out of bounds");
            first = std::min(first, q);
            last = std::max(last, q);
        }
        EnvironmentsForRead(first);
        std::vector<Eigen::Matrix2cd> operators(last - first + 1, Eigen::Matrix2cd::Identity());
        // Repeated operators at a site retain the caller's application order.
        for (const auto &gate : gates)
        {
            auto &op = operators[gate.getQubit1() - first];
            op = gate.getRawOperatorMatrix() * op;
        }
        ExpectationWorkspace workspace;
        for (IndexType q = first; q < last; ++q)
            ContractExpectationSite(q, operators[q - first], q == first, false, workspace);
        return ContractExpectationSite(last, operators[last - first], last == first, true, workspace);
    }

    std::vector<std::complex<double>> ExpectationValuesImpl(const std::vector<std::string> &paulis) const
    {
        struct Query
        {
            size_t index;
            IndexType first, last;
            std::string pauli;
        };

        std::vector<Query> queries;
        std::vector<std::complex<double>> result(paulis.size(), 1.);
        for (size_t i = 0; i < paulis.size(); ++i)
        {
            ValidatePauliString(paulis[i], gammas.size());
            std::string pauli = paulis[i];
            for (char &p : pauli)
                p = CanonicalPauli(p);
            const auto first = pauli.find_first_not_of('I');
            if (first == std::string::npos)
                continue;
            const auto last = pauli.find_last_not_of('I');
            queries.push_back({i, static_cast<IndexType>(first), static_cast<IndexType>(last), std::move(pauli)});
        }
        std::sort(queries.begin(), queries.end(), [](const Query &a, const Query &b) {
            return a.first != b.first ? a.first < b.first : a.pauli < b.pauli;
        });
        if (!queries.empty())
            const_cast<MPSSimulatorImpl *>(this)->EnvironmentsForRead(queries.front().first);
        ExpectationWorkspace workspace;
        std::vector<MatrixClass> prefixes;
        IndexType groupFirst = -1, validUntil = 0;
        for (size_t i = 0; i < queries.size(); ++i)
        {
            const auto &query = queries[i];
            if (i != 0 && query.pauli == queries[i - 1].pauli)
            {
                result[query.index] = result[queries[i - 1].index];
                continue;
            }
            if (query.first != groupFirst)
            {
                prefixes.clear();
                groupFirst = query.first;
                validUntil = groupFirst;
            }
            IndexType start = query.first;
            if (i != 0 && queries[i - 1].first == query.first)
                while (start < std::min(validUntil, query.last) && query.pauli[start] == queries[i - 1].pauli[start])
                    ++start;
            if (start != query.first)
                workspace.left = prefixes[start - query.first - 1];

            // Cache only prefixes that the next query will actually reuse.
            IndexType keepUntil = query.first;
            if (i + 1 < queries.size() && queries[i + 1].first == query.first)
                while (keepUntil < std::min(query.last, queries[i + 1].last) &&
                       query.pauli[keepUntil] == queries[i + 1].pauli[keepUntil])
                    ++keepUntil;
            size_t bytes = 0;
            IndexType boundedUntil = query.first;
            while (boundedUntil < keepUntil)
            {
                const size_t R = static_cast<size_t>(gammas[boundedUntil].dimension(2));
                const size_t available = ObservableCacheBytes - bytes;
                if (available < sizeof(MatrixClass) ||
                    R > (available - sizeof(MatrixClass)) / sizeof(std::complex<double>) / R)
                    break;
                bytes += sizeof(MatrixClass) + R * R * sizeof(std::complex<double>);
                ++boundedUntil;
            }
            // Drop deeper old environments before growing a different prefix.
            prefixes.resize(static_cast<size_t>(boundedUntil - query.first));
            for (IndexType q = start; q < query.last; ++q)
            {
                ContractExpectationSite(q, PauliOperator(query.pauli[q]), q == query.first, false, workspace);
                if (q < boundedUntil)
                    prefixes[q - query.first] = workspace.left;
            }
            result[query.index] = ContractExpectationSite(query.last, PauliOperator(query.pauli[query.last]),
                                                          query.last == query.first, true, workspace);
            validUntil = boundedUntil;
        }
        return result;
    }

    struct ExpectationWorkspace
    {
        MatrixClass left, next, ket0, ket1, work0, work1;
    };

    std::complex<double> ContractExpectationSite(IndexType q, const Eigen::Matrix2cd &op, bool first, bool last,
                                                 ExpectationWorkspace &w) const
    {
        const auto &g = gammas[q];
        const IndexType L = g.dimension(0), R = g.dimension(2);
        using Slice = Eigen::Map<const MatrixClass, 0, Eigen::OuterStride<>>;
        const Slice zero(g.data(), L, R, Eigen::OuterStride<>(2 * L));
        const Slice one(g.data() + L, L, R, Eigen::OuterStride<>(2 * L));
        const auto a = op(0, 0), b = op(0, 1), c = op(1, 0), d = op(1, 1);
        // With explicit environments the interval is closed by the Gram matrices of the blocks around it, where the
        // canonical form has the squared lambdas on its left and nothing to add on its right.
        const IndexType bonds = static_cast<IndexType>(lambdas.size());
        MPSSimulatorImpl *const self = const_cast<MPSSimulatorImpl *>(this);
        if (first && last && !gaugeMode)
        {
            std::complex<double> value = 0.;
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < L; ++l)
                {
                    const double weight = q == 0 ? 1. : lambdas[q - 1][l] * lambdas[q - 1][l];
                    const auto z = zero(l, r), o = one(l, r);
                    value += weight * (std::conj(z) * (a * z + b * o) + std::conj(o) * (c * z + d * o));
                }
            return value;
        }
        if (b == std::complex<double>{} && c == std::complex<double>{})
        {
            w.ket0 = a * zero;
            w.ket1 = d * one;
        }
        else if (a == std::complex<double>{} && d == std::complex<double>{})
        {
            w.ket0 = b * one;
            w.ket1 = c * zero;
        }
        else
        {
            w.ket0 = a * zero + b * one;
            w.ket1 = c * zero + d * one;
        }
        if (first)
        {
            if (q != 0 && gaugeMode)
            {
                const MatrixClass &left = self->LeftEnvironment(q - 1);
                w.work0.noalias() = left * w.ket0;
                w.work1.noalias() = left * w.ket1;
                w.ket0.swap(w.work0);
                w.ket1.swap(w.work1);
            }
            else if (q != 0)
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                    {
                        const double weight = lambdas[q - 1][l] * lambdas[q - 1][l];
                        w.ket0(l, r) *= weight;
                        w.ket1(l, r) *= weight;
                    }
        }
        else
        {
            w.work0.noalias() = w.left * w.ket0;
            w.work1.noalias() = w.left * w.ket1;
            w.ket0.swap(w.work0);
            w.ket1.swap(w.work1);
        }
        if (last)
        {
            if (gaugeMode && q < bonds)
            {
                const MatrixClass &right = self->RightEnvironment(q);
                w.work0.noalias() = w.ket0 * right.transpose();
                w.work1.noalias() = w.ket1 * right.transpose();
                w.ket0.swap(w.work0);
                w.ket1.swap(w.work1);
            }
            return (zero.conjugate().cwiseProduct(w.ket0)).sum() + (one.conjugate().cwiseProduct(w.ket1)).sum();
        }
        w.next.noalias() = zero.adjoint() * w.ket0;
        w.next.noalias() += one.adjoint() * w.ket1;
        w.left.swap(w.next);
        return 0.;
    }

    // Draws the operator for `random` among the weights. Returns whether the threshold is clear of the two ends of
    // its interval by more than twice `doubt`, which is what the weights in doubt can move them by.
    static bool DrawOperator(const double *weights, size_t count, double total, double doubt, double random,
                             size_t &chosen)
    {
        const double threshold = random * total;
        double accumulated = 0., lower = 0.;
        for (size_t k = 0; k < count; ++k)
        {
            if (weights[k] <= 0.)
                continue;
            chosen = k; // also the fallback if rounding leaves the threshold above the last sum
            lower = accumulated;
            accumulated += weights[k];
            if (threshold < accumulated)
                break;
        }
        return threshold - lower > 2. * doubt && accumulated - threshold > 2. * doubt;
    }

    // Whether an operator nearly projects: its smallest singular value is below about a hundredth of the others.
    // From the determinant, the product of the singular values, against the squared Frobenius norm, their sum of
    // squares: for two of them s and r s, |det| / norm^2 is r / (1 + r^2).
    static bool ProjectionLike(const Eigen::Matrix2cd &op)
    {
        return !(std::abs(op.determinant()) > 1e-2 * op.squaredNorm());
    }

    // Whether a one-site operator has rank one, to rounding: a projection in the wide sense, |a><b|.
    static bool RankOne(const Eigen::Matrix2cd &op)
    {
        return !(std::abs(op.determinant()) > 1e-14 * op.squaredNorm());
    }

    // The state a rank-one operator |a><b| leaves the site in: its column of larger norm, normalized (the columns
    // are multiples of each other). Zero for the zero operator.
    static Eigen::Vector2cd KeptDirection(const Eigen::Matrix2cd &op)
    {
        Eigen::Vector2cd a = op.col(0).squaredNorm() >= op.col(1).squaredNorm() ? op.col(0) : op.col(1);
        const double length = a.norm();
        if (length > 0.)
            a /= length;
        else
            a.setZero();
        return a;
    }

    // Whether an operator leaves the site in one state whatever the rest holds, so that the site can be split again
    // around that state (see RegaugeSite). An operator of rank one does, on any state. One that has rank one only
    // to rounding, as diag(1, 1e-15), does not on a state that holds little of what it keeps: what it leaves outside
    // its larger direction is then part of the result, or all of it. So the test is on the state. `rho` is the
    // density matrix of the site before the operator and `weight` = tr(op rho op^dagger); the part of that weight
    // outside the direction kept has to be below the rounding of the amplitudes.
    static bool LeavesOneState(const Eigen::Matrix2cd &op, const Eigen::Matrix2cd &rho, double weight)
    {
        if (!RankOne(op))
            return false;
        const Eigen::Vector2cd a = KeptDirection(op);
        // <a_perp| op, with a_perp = (-conj(a1), conj(a0)) the state orthogonal to a
        const Eigen::RowVector2cd outside = -a[1] * op.row(0) + a[0] * op.row(1);
        const double left = (outside * rho * outside.adjoint())(0, 0).real();
        return left <= AmplitudeRounding * AmplitudeRounding * weight;
    }

    // A chain of one qubit has no bond, so the sweeps that normalize from the tensors have nothing to do there: its
    // tensor is divided by its own norm.
    void NormalizeLoneSite()
    {
        if (!lambdas.empty() || gammas.empty())
            return;
        Eigen::Map<Eigen::ArrayXcd> values(gammas[0].data(), gammas[0].size());
        const double norm = std::sqrt(values.abs2().sum());
        if (norm > 0. && std::isfinite(norm))
            values /= norm;
    }

    // No operator is left to draw. When the weights were read from the tensors and the largest is positive, it is
    // reported with the cutoff it fell below: such a weight is not told from the rounding of amplitudes that cancel.
    [[noreturn]] static void ThrowNoProbability(double largest, double cutoff)
    {
        if (!(largest > 0.))
            throw std::domain_error("No Kraus operator has any probability on this state");
        std::ostringstream message;
        message << "No Kraus operator has any probability on this state: the largest weight found, "
                << std::setprecision(3) << largest << ", is below " << cutoff
                << ", the rounding of amplitudes that cancel, and counts as none";
        throw std::domain_error(message.str());
    }

    // After an operator |a><b| on `site`, in gauge mode: every slice of the tensor is a multiple of one matrix,
    // T(l, s, r) = a(s) m(l, r), and the site is in the state a whatever the rest holds. With m = U S V^dagger, U goes
    // into the tensor on the left, V^dagger into the one on the right, and the site keeps a(s) S: the same state,
    // with the two bonds of the site at the rank the projection left them and in bases in which their
    // environments are graded by S instead of having lost their rank. This is the local form of what the sweeps of
    // a measurement do on the canonical path, and costs one decomposition of a bond-sized matrix. The lambdas of
    // the two bonds only keep their dimensions here: nothing reads them in gauge mode.
    void RegaugeSite(IndexType site, const Eigen::Matrix2cd &op)
    {
        const IndexType sites = static_cast<IndexType>(gammas.size());
        const GammaType &gamma = gammas[site];
        const IndexType L = gamma.dimension(0), R = gamma.dimension(2);
        // a: the state the operator leaves the site in. The caller has checked that it does (see LeavesOneState).
        const Eigen::Vector2cd a = KeptDirection(op);
        if (!(a.squaredNorm() > 0.))
            return;
        const MatrixClass m = std::conj(a[0]) * Slice(gamma, 0) + std::conj(a[1]) * Slice(gamma, 1);
        Eigen::BDCSVD<MatrixClass> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if (svd.info() != Eigen::Success)
            return; // left as it is: a gate that meets the environments restores the canonical form
        const IndexType rank = std::max<IndexType>(1, svd.rank());
        const MatrixClass U = svd.matrixU().leftCols(rank), V = svd.matrixV().leftCols(rank);
        const auto S = svd.singularValues().head(rank);

        GammaType middle(rank, 2, rank);
        middle.setZero();
        for (IndexType s = 0; s < 2; ++s)
            for (IndexType k = 0; k < rank; ++k)
                middle(k, s, k) = a[s] * S[k];
        if (site > 0)
        {
            const GammaType &left = gammas[site - 1];
            const IndexType leftL = left.dimension(0);
            GammaType updated(leftL, 2, rank);
            for (IndexType s = 0; s < 2; ++s)
                Eigen::Map<MatrixClass, 0, Eigen::OuterStride<>>(updated.data() + s * leftL, leftL, rank,
                                                                 Eigen::OuterStride<>(2 * leftL))
                    .noalias() = Slice(left, s) * U;
            gammas[site - 1] = std::move(updated);
            lambdas[site - 1] = LambdaType::Ones(rank);
            bondDimensions.Update(site - 1, rank);
        }
        if (site + 1 < sites)
        {
            const GammaType &right = gammas[site + 1];
            const IndexType rightR = right.dimension(2);
            GammaType updated(rank, 2, rightR);
            for (IndexType s = 0; s < 2; ++s)
                Eigen::Map<MatrixClass, 0, Eigen::OuterStride<>>(updated.data() + s * rank, rank, rightR,
                                                                 Eigen::OuterStride<>(2 * rank))
                    .noalias() = V.adjoint() * Slice(right, s);
            gammas[site + 1] = std::move(updated);
            lambdas[site] = LambdaType::Ones(rank);
            bondDimensions.Update(site, rank);
        }
        // At an end of the chain the missing side has dimension one and its factor is a phase: it stays in the site.
        if (site == 0)
            for (IndexType s = 0; s < 2; ++s)
                for (IndexType k = 0; k < rank; ++k)
                    middle(k, s, k) *= U(0, k);
        if (site + 1 == sites)
            for (IndexType s = 0; s < 2; ++s)
                for (IndexType k = 0; k < rank; ++k)
                    middle(k, s, k) *= std::conj(V(0, k));
        gammas[site] = std::move(middle);
        // The two bonds have new bases: their environments, and those that contain the site, are computed again.
        leftValid = std::min(leftValid, std::max<IndexType>(site - 1, 0));
        rightValid = std::max(rightValid, std::min<IndexType>(site + 1, static_cast<IndexType>(lambdas.size())));
    }

    static bool ProjectionLike(const Eigen::Matrix4cd &op)
    {
        const double mean = op.squaredNorm() / 4.;
        return !(std::abs(op.determinant()) > 1e-2 * mean * mean);
    }

    // On an exact canonical form the block to the left of a site has the squared lambdas of its left bond as weights
    // and the block to its right is orthonormal, so the squared norm of the state after a one-site operator is
    // sum lambda(l)^2 |sum_a op(s, a) B(l, a, r)|^2: a sum of squares of what the operator leaves of the tensor, which
    // holds down to the rounding of the tensor itself. `reach` receives the bound on the terms of the same weight
    // taken from the density matrix of the site, 2 sum_a |column a of op|^2 rho(a, a).
    double CanonicalSiteWeight(IndexType site, const Eigen::Matrix2cd &op, double &reach) const
    {
        const GammaType &gamma = gammas[site];
        const IndexType L = gamma.dimension(0), R = gamma.dimension(2);
        double weight = 0., zeroPopulation = 0., onePopulation = 0.;
        for (IndexType r = 0; r < R; ++r)
            for (IndexType l = 0; l < L; ++l)
            {
                const double lambda2 = site == 0 ? 1. : lambdas[site - 1][l] * lambdas[site - 1][l];
                const std::complex<double> zero = gamma(l, 0, r), one = gamma(l, 1, r);
                weight += lambda2 *
                          (std::norm(op(0, 0) * zero + op(0, 1) * one) + std::norm(op(1, 0) * zero + op(1, 1) * one));
                zeroPopulation += lambda2 * std::norm(zero);
                onePopulation += lambda2 * std::norm(one);
            }
        reach = 2. * (op.col(0).squaredNorm() * zeroPopulation + op.col(1).squaredNorm() * onePopulation);
        return weight;
    }

    // M = B_i reshaped as (L, 2R). QR of M^dagger gives M = R^dagger Q^dagger:
    // keep Q^dagger at this site and move R^dagger left. No lambda or rank cutoff.
    void RightCanonicalizeSite(IndexType site)
    {
        const GammaType &gamma = gammas[site];
        const IndexType L = gamma.dimension(0);
        const IndexType R = gamma.dimension(2);
        const IndexType rank = std::min(L, 2 * R);
        MatrixClass adjoint(2 * R, L);
        for (IndexType physical = 0; physical < 2; ++physical)
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < L; ++l)
                    adjoint(physical * R + r, l) = std::conj(gamma(l, physical, r));

        Eigen::HouseholderQR<MatrixClass> qr(adjoint);
        const MatrixClass Q = qr.householderQ() * MatrixClass::Identity(2 * R, rank);
        const MatrixClass triangular = qr.matrixQR().topRows(rank).template triangularView<Eigen::Upper>();
        const GammaType &previous = gammas[site - 1];
        const IndexType previousL = previous.dimension(0);
        const Eigen::Map<const MatrixClass> previousMatrix(previous.data(), 2 * previousL, L);
        const MatrixClass left = previousMatrix * triangular.adjoint();
        GammaType newLeft(previousL, 2, rank);
        std::copy(left.data(), left.data() + left.size(), newLeft.data());
        GammaType newRight(rank, 2, R);
        for (IndexType physical = 0; physical < 2; ++physical)
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < rank; ++l)
                    newRight(l, physical, r) = std::conj(Q(physical * R + r, l));
        gammas[site - 1] = std::move(newLeft);
        gammas[site] = std::move(newRight);
        lambdas[site - 1] = LambdaType::Ones(rank); // dimensions only until the forward sweep
        bondDimensions.Update(site - 1, lambdas[site - 1].size());
    }

    // Brings a measurement probability into [0, 1]. One that is outside by more than rounding is reported on cerr,
    // with its value and what the settings say of the cause, and the simulation goes on with the clamped value.
    double ValidMeasurementProbability(double probability) const
    {
        constexpr double toleranceLow = 1E-9;
        constexpr double toleranceHigh = 1E-5;
        if (probability < -toleranceLow || probability > 1. + toleranceHigh)
            ReportInvalidProbability(probability);

        return ClampProbability(std::clamp(probability, 0., 1.));
    }

    // Only what is already known goes into the message: the probability, the truncation settings, the bond
    // dimensions and whether anything was truncated. The canonical form is exact when a probability is taken (see
    // EnsureCanonical), so a value out of range means that the state itself is off its norm.
    void ReportInvalidProbability(double probability) const
    {
        std::ostringstream message;
        message << "Invalid measurement probability produced by the MPS state: " << std::setprecision(12) << probability
                << std::setprecision(3);
        if (probability > 1.)
            message << " (above 1 by " << probability - 1. << "). ";
        else
            message << " (below 0 by " << -probability << "). ";
        // A truncation alone no longer gives this: the form is rebuilt, and the state normalized, before a read that
        // follows one. What was truncated is still said, as it tells how approximate the state is.
        message << "The state is off its norm, as after a gate that is not unitary. Truncation is not the cause: ";
        if (schmidtValuesDiscarded)
        {
            message << "Schmidt values have been discarded (";
            if (limitSize)
                message << "bond limit " << chi << ", ";
            if (limitEntanglement)
                message << "singular value threshold " << singularValueThreshold
                        << (truncationMode == TruncationMode::RelativeToMax ? " relative to the largest, "
                                                                            : " of discarded weight, ");
            message << "largest bond " << getMaxBondDimension() << "), but the canonical form was rebuilt after that.";
        }
        else
            message << "no Schmidt value has been discarded (largest bond " << getMaxBondDimension() << ").";
        std::cerr << message.str() << std::endl;
    }

    static bool PickZeroMeasurement(double rndVal, double prob0)
    {
        if (prob0 <= 0.)
            return false;
        if (prob0 >= 1.)
            return true;

        return rndVal < prob0;
    }

    static void SwapTheta(Eigen::Tensor<std::complex<double>, 4> &theta)
    {
        for (IndexType j = 0; j < theta.dimension(3); ++j)
            for (IndexType i = 0; i < theta.dimension(0); ++i)
                std::swap(theta(i, 0, 1, j), theta(i, 1, 0, j));
    }

    MatrixClass ConstructTheta(const GateClass &gate, IndexType qubit1, bool dontApplyGate, bool isSwapGate,
                               bool reversed)
    {
        MatrixClass thetaMatrix;

        // TODO: optimize for some gates? The most important would be the swap gate, as it's applied a lot to move
        // qubits near each other before applying other two qubit gates something like:
        if (dontApplyGate || isSwapGate)
        {
            Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);
            // eiter apply the swap gate, optimized
            // or this case is for dealing with entanglement with the qubits entangled with a measured one
            // see 'Measure' function for details

            // NOTE: with proper data structures this could be way more optimized
            // instead of those generic eigen tensors for example, a more optimized tensor for this would hold Eigen
            // matrices inside and the physical indices would select among them that way swapping will reduce to
            // swapping pointers, it would be very fast but I'm lazy so I use the generic implementation

            // qiskit aer does it, though

            // it's less important than it appears, because usually SVD will take the most time (unless severely limited
            // or having some particularly easy circuits)

            if (isSwapGate)
                SwapTheta(theta);

            thetaMatrix = ReshapeTheta(theta);
        }
        else
        {
            const Eigen::TensorFixedSize<std::complex<double>, Eigen::Sizes<2, 2, 2, 2>> U =
                GetTwoQubitsGateTensor(gate, reversed);
            const Eigen::Tensor<std::complex<double>, 4> thetabar = ConstructThetaBar(qubit1, U);

            thetaMatrix = ReshapeThetaBar(thetabar);
        }

        return thetaMatrix;
    }

    void ApplyTwoQubitGate(const GateClass &gate, IndexType qubit, IndexType controllingQubit1,
                           bool dontApplyGate = false)
    {
        // it's more complex than the single qubit gate
        // very shortly:
        // contract tensors for the two qubits (the B tensors already contain the lambdas, except the left one)
        // shape the gate into a tensor
        // contract the gate tensor with the two qubit tensor
        // apply SVD (with the left lambda multiplied in) to separate out the resulting tensor into the two qubit
        // tensors and the lambda between them

        IndexType qubit1 = controllingQubit1;
        IndexType qubit2 = qubit;
        bool reversed = false;

        if (qubit1 > qubit2)
        {
            std::swap(qubit1, qubit2);
            reversed = true;
        }

        const bool isSwapGate = gate.isSwapGate();
        if (!gaugeMode && !dontApplyGate && isSwapGate && TrySwapProductSites(qubit1))
            return;

        // the size is (2 * left bond dim) x (2 * right bond dim)
        const MatrixClass thetaMatrix = ConstructTheta(gate, qubit1, dontApplyGate, isSwapGate, reversed);

        if (gaugeMode)
        {
            const auto &op = gate.getRawOperatorMatrix();
            const bool unitary = dontApplyGate || isSwapGate || (op.adjoint() * op).isIdentity(1e-10);
            bool decomposed = false;
            RunMaybeSingleThreaded(enableMultithreading, [&]() {
                decomposed = DecomposeWithEnvironments(thetaMatrix, qubit1, qubit2, unitary);
            });
            if (decomposed)
                return;
            // An environment has lost rank (a projection was applied). The canonical form is restored for this
            // update, whose own environments are then exact, and the explicit ones go on from it. The update may
            // truncate, which leaves the left site off the canonical form; but the form was exact until then, so
            // it still describes the left blocks up to the bond of the update and the right blocks from it on.
            // The others, which contain the two sites, are recomputed.
            RestoreCanonicalForm(false);
            DecomposeAndSetGammas(ConstructTheta(gate, qubit1, dontApplyGate, isSwapGate, reversed), qubit1, qubit2);
            SetCanonicalEnvironments();
            leftValid = qubit1 + 1;
            rightValid = qubit1;
            return;
        }

        DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2);
    }

    // ---- gauge mode (see ApplyKrausOperator) ----

    using SliceMap = Eigen::Map<const MatrixClass, 0, Eigen::OuterStride<>>;

    // The matrix of a site tensor for one value of its physical index.
    static SliceMap Slice(const GammaType &gamma, IndexType physical)
    {
        const IndexType L = gamma.dimension(0);
        return SliceMap(gamma.data() + physical * L, L, gamma.dimension(2), Eigen::OuterStride<>(2 * L));
    }

    // The canonical form describes the same state: restoring it changes no observable. None of the cuts the user
    // asked for is applied here, so that reading a state (a probability, an expectation value, a saved state) leaves
    // it as it is; the gates and the measurements that follow apply them as they do on a chain that never left the
    // form.
    //
    // This is for what needs the chain on its canonical path: saving it, trimming it. It ends the gauge mode, by
    // restoring the form, or by dropping the environments when they were only kept for reads. It does not repair a
    // form that gates have left approximate by discarding Schmidt values (see canonicalFormExact): what reads such a
    // chain takes the environments from the tensors instead (see EnvironmentsForRead), which costs matrix products
    // where a repair costs a decomposition per site.
    void EnsureCanonical()
    {
        if (!gaugeMode)
            return;
        if (gaugeForReads)
            LeaveGaugeMode();
        else
            RestoreCanonicalForm(false);
    }

    void EnsureCanonical() const
    {
        const_cast<MPSSimulatorImpl *>(this)->EnsureCanonical();
    }

    // Makes the environments available to a read that does not change the chain (a probability, an expectation
    // value, a sample), with `site` one of the sites read. On an exact canonical form the lambdas and the
    // orthonormal tensors are the environments, and nothing is done. On a form that truncations have left
    // approximate they are computed from the tensors, by matrix products, and kept for the reads that follow: the
    // chain is then in gauge mode for reading only, and the next gate drops them and goes on as before.
    void EnvironmentsForRead(IndexType site)
    {
        if (gaugeMode || canonicalFormExact)
            return;
        EnterGaugeModeNormalized(site);
        gaugeForReads = true;
    }

    void LeaveGaugeMode()
    {
        gaugeMode = false;
        gaugeForReads = false;
        leftEnv.clear();
        rightEnv.clear();
        leftCanonical.clear();
        rightCanonical.clear();
        leftValid = rightValid = 0;
    }

    // The environments of an exact canonical form: the left block of a bond has the squared Schmidt values as its
    // Gram matrix and the right block, built of right-orthonormal tensors, the identity.
    void SetCanonicalEnvironments()
    {
        const IndexType bonds = static_cast<IndexType>(lambdas.size());
        leftEnv.resize(bonds);
        rightEnv.resize(bonds);
        for (IndexType bond = 0; bond < bonds; ++bond)
        {
            leftEnv[bond] = lambdas[bond].cwiseAbs2().template cast<std::complex<double>>().asDiagonal();
            rightEnv[bond] = MatrixClass::Identity(lambdas[bond].size(), lambdas[bond].size());
        }
        leftCanonical.assign(bonds, 1);
        rightCanonical.assign(bonds, 1);
        leftValid = bonds;
        rightValid = 0;
        gaugeMode = true;
    }

    // The environments come from the canonical form while it is exact. Once a gate has discarded Schmidt values it
    // no longer gives them (see canonicalFormExact), and they are all computed from the tensors, when first needed.
    void EnterGaugeMode()
    {
        if (gaugeMode)
            return;
        if (canonicalFormExact)
        {
            SetCanonicalEnvironments();
            return;
        }
        const IndexType bonds = static_cast<IndexType>(lambdas.size());
        leftEnv.resize(bonds);
        rightEnv.resize(bonds);
        leftCanonical.assign(bonds, 0);
        rightCanonical.assign(bonds, 0);
        leftValid = 0;
        rightValid = bonds;
        gaugeMode = true;
    }

    // Enters the gauge mode before a read of `site` on a chain whose canonical form is not exact. Such a chain is
    // not exactly normalized either: each truncation normalized it by lambdas that the next one left approximate.
    // The environments of the site, computed from the tensors, give the norm, and it is put into the tensor of the
    // site, which none of the environments computed so far contains. The gauge mode then starts, as it does from an
    // exact form, on a state of unit norm.
    void EnterGaugeModeNormalized(IndexType site)
    {
        EnterGaugeMode();
        const double norm = SiteDensity(site).trace().real();
        if (norm > 0. && std::isfinite(norm))
            Eigen::Map<Eigen::ArrayXcd>(gammas[site].data(), gammas[site].size()) /= std::sqrt(norm);
    }

    // G(r', r) = <L_r'|L_r> for the block of the sites up to `bond`: G_j = sum_s M_s^dagger G_(j-1) M_s.
    const MatrixClass &LeftEnvironment(IndexType bond)
    {
        for (IndexType b = leftValid; b <= bond; ++b)
        {
            const GammaType &gamma = gammas[b];
            const IndexType R = gamma.dimension(2);
            MatrixClass env = MatrixClass::Zero(R, R);
            for (IndexType s = 0; s < 2; ++s)
            {
                const SliceMap m = Slice(gamma, s);
                if (b == 0)
                    env.noalias() += m.adjoint() * m;
                else
                {
                    const MatrixClass weighted = leftEnv[b - 1] * m;
                    env.noalias() += m.adjoint() * weighted;
                }
            }
            leftEnv[b] = std::move(env);
            leftCanonical[b] = 0;
        }
        leftValid = std::max(leftValid, bond + 1);
        return leftEnv[bond];
    }

    // H(l', l) = <R_l'|R_l> for the block of the sites after `bond`: H_(j-1) = sum_s conj(M_s) H_j M_s^T.
    const MatrixClass &RightEnvironment(IndexType bond)
    {
        const IndexType bonds = static_cast<IndexType>(lambdas.size());
        for (IndexType b = rightValid - 1; b >= bond; --b)
        {
            const GammaType &gamma = gammas[b + 1];
            const IndexType L = gamma.dimension(0);
            MatrixClass env = MatrixClass::Zero(L, L);
            for (IndexType s = 0; s < 2; ++s)
            {
                const SliceMap m = Slice(gamma, s);
                if (b + 1 == bonds)
                    env.noalias() += m.conjugate() * m.transpose();
                else
                {
                    const MatrixClass weighted = m.conjugate() * rightEnv[b + 1];
                    env.noalias() += weighted * m.transpose();
                }
            }
            rightEnv[b] = std::move(env);
            rightCanonical[b] = 0;
        }
        rightValid = std::min(rightValid, bond);
        return rightEnv[bond];
    }

    // rho(s, s') = sum conj(M_s'(l', r')) (G M_s H^T)(l', r'), with G and H the environments of the site.
    Eigen::Matrix2cd SiteDensity(IndexType site)
    {
        const GammaType &gamma = gammas[site];
        const IndexType bonds = static_cast<IndexType>(lambdas.size());
        MatrixClass weighted[2];
        for (IndexType s = 0; s < 2; ++s)
        {
            const SliceMap m = Slice(gamma, s);
            if (site == 0)
                weighted[s] = m;
            else if (gaugeMode)
                weighted[s].noalias() = LeftEnvironment(site - 1) * m;
            else
                weighted[s] = lambdas[site - 1].cwiseAbs2().template cast<std::complex<double>>().asDiagonal() * m;
            if (gaugeMode && site < bonds)
            {
                const MatrixClass right = weighted[s] * RightEnvironment(site).transpose();
                weighted[s] = right;
            }
        }
        Eigen::Matrix2cd rho;
        for (IndexType s = 0; s < 2; ++s)
            for (IndexType s2 = 0; s2 < 2; ++s2)
                rho(s, s2) = Slice(gamma, s2).conjugate().cwiseProduct(weighted[s]).sum();
        return rho;
    }

    // DecomposeAndSetGammas for a chain in gauge mode: theta is weighed by factors of both environments instead of
    // the left lambda. With G = A^dagger A to the left and H = D^dagger D to the right, the norm of the state is the
    // Frobenius norm of theta' = (1 (x) A) theta (1 (x) D)^T, so the singular values of theta' are the Schmidt
    // values of the middle bond and cutting them is the optimal truncation. With theta' = U S W^dagger, the new
    // tensors are X = theta (1 (x) D)^T W, which is (1 (x) A)^-1 U S computed without inverting A, and
    // Y = W^dagger (1 (x) D)^-T. In the canonical form A is the diagonal of the left lambdas and D the identity,
    // and this is Hastings' update.
    //
    // The factors come from the environments, which are Gram matrices, and a Gram matrix holds a weight only to the
    // rounding of its entries. This update is therefore used only where its factors are as good as the canonical
    // ones, and returns false otherwise, with nothing changed, for the caller to restore the canonical form:
    //  - An environment still in its canonical form (the diagonal of the squared lambdas, the identity) is used as
    //    it is: its square root is the lambdas themselves.
    //  - A left environment accumulated through products is factored after scaling to a unit diagonal. A Gram
    //    matrix is graded like the weights of its bond, and its entries are known relative to the diagonal ones on
    //    their row and column, so the scaling takes the grading out exactly and what is left says how well the
    //    directions are told apart: a direction of small weight is as accurate as a large one, and the pivots of the
    //    scaled matrix fall only when the vectors of the left block have become nearly dependent, which is what a
    //    projection (a jump, a reset) does. Below ScaledPivotLimit the update is refused.
    //  - The right factor is also inverted. An accumulated right environment is scaled in the same way, the
    //    grading going into a diagonal that is inverted exactly; when the pivots of the scaled matrix spread beyond
    //    RightPivotLimit, or a direction has no weight at all, the update is refused.
    // No singular value is dropped beyond the numerical rank of the decomposition, as in the canonical form: when
    // the factors are accepted they resolve the small ones too.
    //
    // `unitary` says that theta holds a unitary gate (or none): without truncation the other environments then keep
    // their values.
    bool DecomposeWithEnvironments(const MatrixClass &thetaMatrix, IndexType qubit1, IndexType qubit2, bool unitary)
    {
        const IndexType bonds = static_cast<IndexType>(lambdas.size());
        const IndexType szl = thetaMatrix.rows() / 2, szr = thetaMatrix.cols() / 2;

        MatrixClass weighted = thetaMatrix;
        // The environments are brought up to date before their flags are read: recomputing one clears its flag.
        const MatrixClass *const leftGram = qubit1 > 0 ? &LeftEnvironment(qubit1 - 1) : nullptr;
        const MatrixClass *const rightGram = qubit2 < bonds ? &RightEnvironment(qubit2) : nullptr;
        if (leftGram && leftCanonical[qubit1 - 1])
        {
            const Eigen::VectorXd scale = leftGram->diagonal().real().cwiseMax(0.).cwiseSqrt();
            for (IndexType j = 0; j < 2; ++j)
                weighted.middleRows(j * szl, szl) =
                    scale.template cast<std::complex<double>>().asDiagonal() * thetaMatrix.middleRows(j * szl, szl);
        }
        else if (leftGram)
        {
            // G = S C S with S the square roots of its diagonal and C = P^T L D L^dagger P, so A = sqrt(D) L^dagger P S
            const MatrixClass &gram = *leftGram;
            Eigen::VectorXd diagonal = gram.diagonal().real();
            if (!diagonal.allFinite())
                return false;
            Eigen::VectorXd inverse(szl);
            for (IndexType i = 0; i < szl; ++i)
            {
                if (diagonal[i] > 0.)
                {
                    diagonal[i] = std::sqrt(diagonal[i]);
                    inverse[i] = 1. / diagonal[i];
                }
                else if (diagonal[i] == 0. && gram.row(i).cwiseAbs().maxCoeff() == 0.)
                    inverse[i] = 0.; // a direction without weight: it takes no part
                else
                    return false; // the rounding of a cancellation, not a weight
            }
            MatrixClass scaled = inverse.template cast<std::complex<double>>().asDiagonal() * gram *
                                 inverse.template cast<std::complex<double>>().asDiagonal();
            for (IndexType i = 0; i < szl; ++i)
                if (inverse[i] == 0.)
                    scaled(i, i) = 1.;
            const Eigen::LDLT<MatrixClass> factor(scaled);
            if (factor.info() != Eigen::Success)
                return false;
            const Eigen::VectorXd pivots = factor.vectorD().real();
            if (!pivots.allFinite() || !(pivots.minCoeff() > ScaledPivotLimit * pivots.maxCoeff()))
                return false;
            const Eigen::VectorXd scale = pivots.cwiseSqrt();
            for (IndexType j = 0; j < 2; ++j)
            {
                MatrixClass block =
                    diagonal.template cast<std::complex<double>>().asDiagonal() * thetaMatrix.middleRows(j * szl, szl);
                block = factor.transpositionsP() * block;
                block = factor.matrixU() * block;
                weighted.middleRows(j * szl, szl) = scale.template cast<std::complex<double>>().asDiagonal() * block;
            }
        }
        MatrixClass rightFactor; // D, upper triangular; left empty for a canonical right environment, the identity
        if (rightGram && !rightCanonical[qubit2])
        {
            // H = S C S with S the square roots of its diagonal and C = U^dagger U, so D = U S, still upper triangular.
            // The grading goes into S, which is inverted exactly; a direction without weight cannot be.
            Eigen::VectorXd diagonal = rightGram->diagonal().real();
            if (!diagonal.allFinite() || !(diagonal.minCoeff() > 0.))
                return false;
            diagonal = diagonal.cwiseSqrt();
            const Eigen::VectorXd inverse = diagonal.cwiseInverse();
            const MatrixClass scaled = inverse.template cast<std::complex<double>>().asDiagonal() * (*rightGram) *
                                       inverse.template cast<std::complex<double>>().asDiagonal();
            const Eigen::LLT<MatrixClass, Eigen::Upper> factor(scaled);
            if (factor.info() != Eigen::Success)
                return false;
            rightFactor = factor.matrixU();
            const auto pivots = rightFactor.diagonal().real();
            if (!rightFactor.allFinite() || !(pivots.minCoeff() > RightPivotLimit * pivots.maxCoeff()))
                return false;
            rightFactor = rightFactor * diagonal.template cast<std::complex<double>>().asDiagonal();
            for (IndexType k = 0; k < 2; ++k)
            {
                const MatrixClass block = weighted.middleCols(k * szr, szr) * rightFactor.transpose();
                weighted.middleCols(k * szr, szr) = block;
            }
        }

#ifdef USE_FAST_SVD
        const bool computeWithJacobi = weighted.rows() < blockSizeLimit && weighted.cols() < blockSizeLimit;
        if (computeWithJacobi)
            jacobiSVD.compute(weighted);
        else
            SVD.compute(weighted);
        if ((computeWithJacobi ? jacobiSVD.info() : SVD.info()) != Eigen::Success)
            throw std::runtime_error("MPS two-qubit SVD failed (non-finite or numerically invalid input)");
        const MatrixClass &VmatrixFull = computeWithJacobi ? jacobiSVD.matrixV() : SVD.matrixV();
        const LambdaType &SvaluesFull = computeWithJacobi ? jacobiSVD.singularValues() : SVD.singularValues();
        const IndexType numericalRank = computeWithJacobi ? jacobiSVD.rank() : SVD.rank();
#else
        jacobiSVD.compute(weighted);
        if (jacobiSVD.info() != Eigen::Success)
            throw std::runtime_error("MPS two-qubit SVD failed (non-finite or numerically invalid input)");
        const MatrixClass &VmatrixFull = jacobiSVD.matrixV();
        const LambdaType &SvaluesFull = jacobiSVD.singularValues();
        const IndexType numericalRank = jacobiSVD.rank();
#endif
        // The same cuts as in the canonical form: these are the Schmidt values of the bond.
        IndexType szm = limitEntanglement
                            ? ComputeCompressedRank(SvaluesFull, numericalRank, truncationMode, singularValueThreshold)
                            : numericalRank;
        if (szm == 0)
            szm = 1;
        const IndexType sz = limitSize ? std::min<IndexType>(chi, szm) : szm;

        LambdaType newLambda = SvaluesFull.head(sz);
        const double norm = newLambda.norm();
        if (norm > 0)
            newLambda /= norm;

        // mixed = (1 (x) D)^T W and rows = W^dagger (1 (x) D)^-T
        MatrixClass mixed = VmatrixFull.leftCols(sz);
        MatrixClass rows = mixed.adjoint();
        if (rightFactor.size() != 0)
        {
            const MatrixClass lower = rightFactor.transpose();
            for (IndexType k = 0; k < 2; ++k)
            {
                const MatrixClass block = lower * mixed.middleRows(k * szr, szr);
                mixed.middleRows(k * szr, szr) = block;
                auto columns = rows.middleCols(k * szr, szr);
                lower.template triangularView<Eigen::Lower>().template solveInPlace<Eigen::OnTheRight>(columns);
            }
            if (!rows.allFinite())
                return false;
        }

        // Build everything first, then commit together, as the canonical decomposition does.
        GammaType newLeft(szl, 2, sz), newRight(sz, 2, szr);
        Eigen::Map<MatrixClass> leftB(newLeft.data(), 2 * szl, sz);
        leftB.noalias() = thetaMatrix * mixed;
        if (norm > 0 && norm != 1.)
            leftB /= norm;
        for (IndexType k = 0; k < szr; ++k)
            for (IndexType j = 0; j < 2; ++j)
            {
                const IndexType jind = j * szr + k;
                for (IndexType i = 0; i < sz; ++i)
                    newRight(i, j, k) = rows(i, jind);
            }
        MatrixClass newLeftEnv = newLambda.cwiseAbs2().template cast<std::complex<double>>().asDiagonal();

        lambdas[qubit1] = std::move(newLambda);
        bondDimensions.Update(qubit1, lambdas[qubit1].size());
        gammas[qubit1] = std::move(newLeft);
        gammas[qubit2] = std::move(newRight);

        // In the new basis of the middle bond the two environments are those of the canonical form.
        leftEnv[qubit1] = std::move(newLeftEnv);
        rightEnv[qubit1] = MatrixClass::Identity(sz, sz);
        leftCanonical[qubit1] = rightCanonical[qubit1] = 1;
        if (sz < numericalRank)
            schmidtValuesDiscarded = true;
        if (unitary && sz >= numericalRank)
        {
            // Nothing else changed: a unitary on the two sites, of which every direction the decomposition found
            // was kept, leaves the left blocks that contain both, and the right blocks that contain both, with the
            // Gram matrices they had.
            leftValid = std::max(leftValid, qubit1 + 1);
            rightValid = std::min(rightValid, qubit1);
        }
        else
        {
            leftValid = qubit1 + 1;
            rightValid = qubit1;
        }
        return true;
    }

    // Explicit environments: leftEnv[bond] is current for bond < leftValid, rightEnv[bond] for bond >= rightValid.
    bool gaugeMode = false;
    // The gauge mode was entered by a read on a chain whose canonical form is approximate, and no Kraus operator
    // has been applied in it: the tensors are still those of the canonical path (see EnvironmentsForRead).
    bool gaugeForReads = false;
    std::vector<MatrixClass> leftEnv, rightEnv;
    IndexType leftValid = 0, rightValid = 0;
    // Whether a current environment still has its canonical form (the diagonal of the squared lambdas on the left,
    // the identity on the right), as set from that form or by a decomposition, and was not accumulated through
    // products: see the cut in DecomposeWithEnvironments.
    std::vector<char> leftCanonical, rightCanonical;
    // The smallest pivot, relative to the largest, of an accumulated left environment scaled to a unit diagonal for
    // which DecomposeWithEnvironments still uses it: the relative error of the weight of a direction is the rounding
    // of the environment over its pivot.
    static constexpr double ScaledPivotLimit = 1e-4;
    // The same for the pivots of the factor of an accumulated right environment, scaled likewise, which is inverted:
    // the rounding of the new right tensor grows as the inverse of this.
    static constexpr double RightPivotLimit = 1e-2;

    // The weight tr(E rho E^dagger) of a Kraus operator is in doubt when it is below this fraction of its reach, the
    // bound on the terms it sums, d sum_a |column a of E|^2 rho(a, a): it is then what their cancellation left, and
    // its relative error is the rounding of rho over this fraction. The usual operators of a noise channel
    // (diagonal, a single entry per column, a multiple of a unitary) have a weight of half or a quarter of their
    // reach on any state, so only an operator nearly orthogonal to the state is concerned.
    static constexpr double UnresolvedWeight = 1e-3;
    // The relative rounding assumed in the entries of a density matrix taken through the environments: times the
    // reach of an operator, how far the weight computed from it can be from the true one.
    static constexpr double DensityRounding = 1e-12;
    // The relative rounding of a tensor entry after a local operator. A weight read on the canonical form that is
    // below its square times the reach is the rounding of entries that cancel exactly, and counts as none.
    static constexpr double AmplitudeRounding = 64. * 2.220446049250313e-16;

    // For L = middle = R = 1, SWAP turns a * b^T into b * a^T. Split that
    // rank-one matrix directly, including the normalization done by the SVD path.
    // The middle lambda is already absorbed in a; it need not equal one.
    bool TrySwapProductSites(IndexType first)
    {
        auto &left = gammas[first];
        auto &right = gammas[first + 1];
        if (left.dimension(0) != 1 || left.dimension(2) != 1 || right.dimension(0) != 1 || right.dimension(2) != 1 ||
            lambdas[first].size() != 1)
            return false;

        const Eigen::Map<const Eigen::Vector2cd> a(left.data()), b(right.data());
        const double aNorm = a.norm(), bNorm = b.norm();
        const double leftWeight = first == 0 ? 1. : std::abs(lambdas[first - 1][0]);
        const double singularValue = leftWeight * aNorm * bNorm;
        // Retain the usual decomposition/error handling for zero or extreme data.
        const double squaredValue = singularValue * singularValue;
        if (!(aNorm > 0.) || !std::isnormal(squaredValue))
            return false;

        const double leftScale = aNorm / singularValue, rightScale = 1. / aNorm;
        const Eigen::Vector2cd newLeft = b * leftScale, newRight = a * rightScale;
        if (!newLeft.allFinite() || !newRight.allFinite())
            return false;

        // All potentially failing work precedes the commit. Tensor moves exchange
        // storage, and the existing bond vector and dimension summary stay size one.
        std::swap(left, right);
        if (leftScale != 1.)
            Eigen::Map<Eigen::Vector2cd>(left.data()) = newLeft;
        if (rightScale != 1.)
            Eigen::Map<Eigen::Vector2cd>(right.data()) = newRight;
        lambdas[first][0] = 1.;
        return true;
    }

    // SVD the (already built) theta matrix and write back the two new B tensors and the lambda
    // in between. Shared by two-qubit gates, Trim, and ReCanonicalize. User-requested
    // chi / entanglement cuts are applied only when applyUserCompression is true (the
    // default); ReCanonicalize disables them during its final gauge refresh (except relative-cutoff cleanup).
    //
    // This is Hastings' method: theta = B1 B2 (with the gate applied) does not contain the left lambda.
    // The SVD is done on lambda_left * theta = U S V^dagger, which is the two site tensor in the Schmidt basis of the
    // left bond, so the singular values are the Schmidt values of the middle bond. Then the new right B is V^dagger and
    // the new left B is theta * V (which is lambda_left^-1 * U * S, but computed without dividing). U is not needed at
    // all.
    void DecomposeAndSetGammas(const MatrixClass &thetaMatrix, IndexType qubit1, IndexType qubit2,
                               bool applyUserCompression = true)
    {
        // the SVD and the matrix product for the new left B are the parts Eigen parallelizes
        RunMaybeSingleThreaded(enableMultithreading,
                               [&]() { DecomposeAndSetGammasImpl(thetaMatrix, qubit1, qubit2, applyUserCompression); });
    }

    void DecomposeAndSetGammasImpl(const MatrixClass &thetaMatrix, IndexType qubit1, IndexType qubit2,
                                   bool applyUserCompression,
                                   IndexType maxBondDimension = std::numeric_limits<IndexType>::max())
    {
        const IndexType szl = qubit1 == 0 ? 1 : lambdas[qubit1 - 1].size();
        const IndexType szr = qubit2 == static_cast<IndexType>(lambdas.size()) ? 1 : lambdas[qubit2].size();

        assert(2 * szl == thetaMatrix.rows());
        assert(2 * szr == thetaMatrix.cols());

        // the rows are (physical index, left bond index) = j * szl + i
        MatrixClass weightedTheta;
        if (qubit1 != 0)
        {
            const LambdaType &leftLambda = lambdas[qubit1 - 1];
            weightedTheta.resize(thetaMatrix.rows(), thetaMatrix.cols());
            for (IndexType c = 0; c < thetaMatrix.cols(); ++c)
                for (IndexType j = 0; j < 2; ++j)
                {
                    const IndexType jszl = j * szl;
                    for (IndexType i = 0; i < szl; ++i)
                        weightedTheta(jszl + i, c) = thetaMatrix(jszl + i, c) * leftLambda[i];
                }
        }
        const MatrixClass &svdMatrix = qubit1 == 0 ? thetaMatrix : weightedTheta;

#ifdef USE_FAST_SVD
        const bool computeWithJacobi = svdMatrix.rows() < blockSizeLimit && svdMatrix.cols() < blockSizeLimit;
#endif

        // n x p is decomposed into U = n x n, singular vals diagonal matrix = n x p, V^t = p x p
        // the thin version U = n x min(n, p), singular vals diagonal matrix = min(n,p) x min(n,p), V^t = min(n,p) x p
        // only V is computed, U is not needed
#ifdef USE_FAST_SVD
        if (computeWithJacobi)
#endif
            jacobiSVD.compute(svdMatrix);
#ifdef USE_FAST_SVD
        else
            SVD.compute(svdMatrix);

        const MatrixClass &VmatrixFull = computeWithJacobi ? jacobiSVD.matrixV() : SVD.matrixV();
        const LambdaType &SvaluesFull = computeWithJacobi ? jacobiSVD.singularValues() : SVD.singularValues();

        // Eigen's own numerical rank (no threshold is set on the SVD objects): only the singular values
        // below diagSize * epsilon * sigma_max are dropped, they are not distinguishable from zero in double precision.
        // Any real cut is only the user's choice (setLimitEntanglement / setLimitBondDimension).
        // A failed decomposition (non-finite input) leaves the previous results in the reused SVD
        // objects; never read them.
        if ((computeWithJacobi ? jacobiSVD.info() : SVD.info()) != Eigen::Success)
            throw std::runtime_error("MPS two-qubit SVD failed (non-finite or numerically invalid input)");
        const IndexType numericalRank = computeWithJacobi ? jacobiSVD.rank() : SVD.rank();
#else
        if (jacobiSVD.info() != Eigen::Success)
            throw std::runtime_error("MPS two-qubit SVD failed (non-finite or numerically invalid input)");
        const MatrixClass &VmatrixFull = jacobiSVD.matrixV();
        const LambdaType &SvaluesFull = jacobiSVD.singularValues();

        const IndexType numericalRank = jacobiSVD.rank();
#endif
        // If user-requested compression is enabled, further reduce the rank according to
        // the configured truncation mode, applied on the (still descending-sorted) singular values.
        IndexType szm = (applyUserCompression && limitEntanglement)
                            ? ComputeCompressedRank(SvaluesFull, numericalRank, truncationMode, singularValueThreshold)
                            : numericalRank;

        if (szm == 0)
            szm = 1; // Shouldn't happen (unless some big limit was put on 'zero')!

        // A gauge-only update must not grow the incoming bond through roundoff.
        const IndexType sz =
            std::min(maxBondDimension, (applyUserCompression && limitSize) ? std::min<IndexType>(chi, szm) : szm);

        assert(sz <= VmatrixFull.cols());
        assert(2 * szr == VmatrixFull.rows());

        // Discarding Schmidt values here leaves the lambdas of the other bonds, and the tensors to the left of this
        // one, only approximately canonical.
        if (sz < numericalRank)
        {
            canonicalFormExact = false;
            if (applyUserCompression)
                schmidtValuesDiscarded = true;
        }

        // Build the new lambda and both B tensors first, then commit them together (noexcept moves),
        // so an allocation failure cannot leave the lambda paired with the old tensors.
        LambdaType newLambda = SvaluesFull.head(sz);
        assert(newLambda[0] != 0.);

        // the norm of the kept singular values is the norm of the state after truncation
        // the lambdas are normalized, and the left B has to be scaled the same way to keep the state normalized
        const double norm = newLambda.norm();
        if (norm > 0)
            newLambda /= norm;

        GammaType newLeft, newRight;
        SetNewGammas(thetaMatrix, VmatrixFull, norm > 0 ? 1. / norm : 1., szl, sz, szr, newLeft, newRight);
        lambdas[qubit1] = std::move(newLambda);
        bondDimensions.Update(qubit1, lambdas[qubit1].size());
        gammas[qubit1] = std::move(newLeft);
        gammas[qubit2] = std::move(newRight);
    }

    // Given the (descending-sorted, already limited to the numerical rank) singular values and a
    // user-requested compression threshold, returns how many of them to keep according to
    // the selected truncation mode. Always keeps at least the largest singular value.
    static IndexType ComputeCompressedRank(const LambdaType &sortedDescendingSVs, IndexType rank, TruncationMode mode,
                                           double threshold)
    {
        if (rank <= 1)
            return rank;

        if (mode == TruncationMode::RelativeToMax)
        {
            // Reproduces Eigen::SVDBase::rank() exactly: a singular value is kept iff it is
            // strictly greater than threshold * sigma_max (with the same guard against a
            // zero sigma_max that Eigen's own premultiplied threshold uses). This is this
            // simulator's original (pre-mode-switch) truncation behavior.
            const double premultipliedThreshold =
                std::max(threshold * sortedDescendingSVs[0], std::numeric_limits<double>::min());
            IndexType i = rank - 1;
            while (i >= 0 && sortedDescendingSVs[i] < premultipliedThreshold)
                --i;
            return i + 1;
        }

        // DiscardedWeight: discard the smallest singular values while the cumulative sum of
        // their squares, normalized by the sum of squares of all of them (at this bond),
        // stays below the threshold. Matches Qiskit Aer's reduce_zeros (see
        // build/qiskit-aer/src/simulators/matrix_product_state/svd.cpp in the maestro repo)
        // and ITensor's default 'cutoff'. Never discards the largest singular value.
        const double total = sortedDescendingSVs.head(rank).squaredNorm();
        if (total <= 0.)
            return rank;

        double discarded = 0.;
        IndexType keep = rank;
        for (IndexType i = rank - 1; i > 0; --i)
        {
            const double sq = sortedDescendingSVs[i] * sortedDescendingSVs[i];
            if ((discarded + sq) / total >= threshold)
                break;

            discarded += sq;
            keep = i;
        }

        return keep;
    }

    static TwoQubitsGateTensor GetTwoQubitsGateTensor(const GateClass &gate, bool reversed)
    {
        TwoQubitsGateTensor result;

        const auto &gateMat = gate.getRawOperatorMatrix();

        if (reversed)
            for (int q0l = 0; q0l < 2; ++q0l)
            {
                const int l0 = q0l << 1;
                for (int q0c = 0; q0c < 2; ++q0c)
                {
                    const int c0 = q0c << 1;
                    for (int q1l = 0; q1l < 2; ++q1l)
                        for (int q1c = 0; q1c < 2; ++q1c)
                            result(q1l, q0l, q1c, q0c) = gateMat(l0 | q1l, c0 | q1c);
                }
            }
        else
            for (int q0l = 0; q0l < 2; ++q0l) // ctrl qubit
            {
                const int l0 = q0l << 1;
                for (int q0c = 0; q0c < 2; ++q0c) // ctrl qubit
                {
                    const int c0 = q0c << 1;
                    for (int q1l = 0; q1l < 2; ++q1l)
                        for (int q1c = 0; q1c < 2; ++q1c)
                            result(q0l, q1l, q0c, q1c) = gateMat(l0 | q1l, c0 | q1c);
                }
            }

        return result;
    }

    // the left B is (theta * V) * scale, the right B is V^dagger
    // V is the full V from the SVD, only the first sz columns are used
    static void SetNewGammas(const MatrixClass &thetaMatrix, const MatrixClass &Vmatrix, double scale, IndexType szl,
                             IndexType sz, IndexType szr, GammaType &newLeft, GammaType &rightB)
    {
        // left site: the (szl, 2, sz) tensor has the same memory layout as the (2 * szl) x sz matrix with the rows j *
        // szl + i so the product can be written directly into it
        newLeft.resize(szl, 2, sz);
        Eigen::Map<MatrixClass> leftB(newLeft.data(), 2 * szl, sz);
        leftB.noalias() = thetaMatrix * Vmatrix.leftCols(sz);
        if (scale != 1.)
            leftB *= scale;

        // right site: the columns of V^dagger are (physical index, right bond index) = j * szr + k
        rightB.resize(sz, 2, szr);
        for (IndexType k = 0; k < szr; ++k)
            for (IndexType j = 0; j < 2; ++j)
            {
                const IndexType jind = j * szr + k;
                for (IndexType i = 0; i < sz; ++i)
                    rightB(i, j, k) = std::conj(Vmatrix(jind, i));
            }
    }

    // the B tensors already contain the lambdas: the middle one is in the left B, the right one in the right B
    // the left lambda is not included, it's needed only for the SVD (see DecomposeAndSetGammas)
    Eigen::Tensor<std::complex<double>, 4> ContractTwoQubits(IndexType qubit1) const
    {
        const IndexType qubit2 = qubit1 + 1;

        static const Indexes product_dims_int{IntIndexPair(2, 0)};

        assert(gammas[qubit1].dimension(2) == gammas[qubit2].dimension(0));

        // contract the first B with the next one
        // the resulting tensor has four legs, 1 and 2 are the physical ones

        return gammas[qubit1].contract(gammas[qubit2], product_dims_int);
    }

    // for the two qubit gate - U is the gate tensor
    Eigen::Tensor<std::complex<double>, 4> ConstructThetaBar(IndexType qubit1, const TwoQubitsGateTensor &U)
    {
        const Eigen::Tensor<std::complex<double>, 4> theta = ContractTwoQubits(qubit1);

        // apply evolution/gate operator
        using DimPair = Eigen::Tensor<std::complex<double>, 4>::DimensionPair;

        // from theta the physical indexes are contracted out
        // the last two become the physical indexes

        static const Eigen::array<DimPair, 2> product_dim{DimPair(1, 2), DimPair(2, 3)};

        // this applies the time evolution/gate operator U
        return theta.contract(U, product_dim);
    }

    MatrixClass ReshapeThetaBar(const Eigen::Tensor<std::complex<double>, 4> &theta)
    {
        // get it into a matrix for SVD - use JacobiSVD
        IndexType sz0 = theta.dimension(0);
        IndexType sz1 = theta.dimension(1);

        assert(sz0 > 0);
        assert(sz1 > 0);

        assert(2 == theta.dimension(2));
        assert(2 == theta.dimension(3));

        MatrixClass thetaMatrix(2 * sz0, 2 * sz1);

        for (IndexType l = 0; l < 2; ++l)
        {
            const IndexType lsz1 = l * sz1;
            for (IndexType k = 0; k < 2; ++k)
            {
                const IndexType ksz0 = k * sz0;
                for (IndexType j = 0; j < sz1; ++j)
                {
                    const IndexType lsz1j = lsz1 + j;
                    for (IndexType i = 0; i < sz0; ++i)
                        // 2, 0, 3, 1 - k & l from theta are the physical indexes
                        thetaMatrix(ksz0 + i, lsz1j) = theta(i, j, k, l);
                }
            }
        }

        return thetaMatrix;
    }

    static MatrixClass ReshapeTheta(const Eigen::Tensor<std::complex<double>, 4> &theta)
    {
        // get it into a matrix for SVD - use JacobiSVD
        IndexType sz0 = theta.dimension(0);
        assert(sz0 > 0);

        assert(2 == theta.dimension(1));
        assert(2 == theta.dimension(2));

        IndexType sz3 = theta.dimension(3);
        assert(sz3 > 0);

        MatrixClass thetaMatrix(2 * sz0, 2 * sz3);

        for (IndexType l = 0; l < sz3; ++l)
            for (IndexType k = 0; k < 2; ++k)
            {
                const IndexType ksz3l = k * sz3 + l;
                for (IndexType j = 0; j < 2; ++j)
                {
                    const IndexType jsz0 = j * sz0;
                    for (IndexType i = 0; i < sz0; ++i)
                        // j & k from theta are the physical indexes
                        thetaMatrix(jsz0 + i, ksz3l) = theta(i, j, k, l);
                }
            }

        return thetaMatrix;
    }

    std::unordered_map<IndexType, bool> MeasureNoCollapse(IndexType limit)
    {
        if (gammas.empty())
            return {};

        const IndexType limit1 = limit + 1;
        std::unordered_map<IndexType, bool> res;

        // the product of the B matrices (the lambdas are included in them) for the values measured so far,
        // divided by the square root of their probability: the sites to the right are in the right canonical
        // form, so the conditional probability of 0 is the squared norm of the zero branch of this row vector.
        // Keeping the row normalized avoids the joint probability of the prefix, which underflows for
        // long chains (2^-n for n unbiased qubits) and made every later outcome 1.
        Eigen::RowVectorXcd vec = Eigen::RowVectorXcd::Ones(1);

        // the matrix for a physical index p starts at offset p * dim1 in the tensor data, the columns are 2 * dim1
        // apart
        using SliceMap = Eigen::Map<const MatrixClass, 0, Eigen::OuterStride<>>;

        for (IndexType qubit = 0; qubit < limit1; ++qubit)
        {
            const IndexType dim1 = gammas[qubit].dimension(0);
            const IndexType dim2 = gammas[qubit].dimension(2);

            // 1. First, compute probability for measuring 0
            // zero matrix
            const SliceMap zeroMat(gammas[qubit].data(), dim1, dim2, Eigen::OuterStride<>(2 * dim1));
            Eigen::RowVectorXcd zeroVec = vec * zeroMat;

            // With explicit environments the block to the right of the site is not orthonormal: the squared norm
            // of a branch is taken in the Gram matrix of that block, v H^T v^dagger.
            const MatrixClass *const right =
                gaugeMode && qubit < static_cast<IndexType>(lambdas.size()) ? &RightEnvironment(qubit) : nullptr;
            const auto squaredNorm = [right](const Eigen::RowVectorXcd &branch) {
                return right ? (branch * right->transpose()).dot(branch).real() : branch.squaredNorm();
            };

            // the probability for the current qubit to be 0, conditioned on the values picked so far
            const double prob0 = ValidMeasurementProbability(squaredNorm(zeroVec));

            // 2. use that probability to measure the qubit
            const double rndVal = 1. - uniformZeroOne(rng);
            const bool zeroMeasured = PickZeroMeasurement(rndVal, prob0);
            res[qubit] = !zeroMeasured;

            // now update the vector, renormalized by the probability of the picked value
            if (zeroMeasured) // no need to compute it again if 0 was measured, it was already computed above
            {
                vec.swap(zeroVec);
                vec /= std::sqrt(prob0);
            }
            else
            {
                const SliceMap oneMat(gammas[qubit].data() + dim1, dim1, dim2, Eigen::OuterStride<>(2 * dim1));
                vec = vec * oneMat;
                // In the Gram matrix the branch is normalized by its own norm, which is 1 - prob0 only up to the
                // rounding of the environments.
                const double prob1 = right ? squaredNorm(vec) : 1. - prob0;
                vec /= std::sqrt(prob1 > 0. ? prob1 : 1. - prob0);
            }
        }

        return res;
    }

    // with Hastings' method the U matrix from the SVD is not needed, only V
#ifdef USE_FAST_SVD
    constexpr static IndexType blockSizeLimit = 16;
    Eigen::BDCSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinV> SVD;
#endif
    Eigen::JacobiSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinV> jacobiSVD;
};

} // namespace TensorNetworks
} // namespace QC
