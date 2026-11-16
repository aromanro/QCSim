#pragma once

#include "MPOSimulatorBase.h"
#include "StructuredGateOperations.h"
#include "TableTrain.h"
#include "ThreeSiteTensor.h"

#include <iostream>

#define USE_FAST_SVD 1

namespace QC
{

namespace TensorNetworks
{

// The actual MPO simulator working on adjacent qubits.
// Two qubit gates need to act on adjacent qubits; the MPOSimulator decorator
// (see MPOSimulator.h) takes care of swapping non adjacent qubits together.
// If compression is enabled, the SVD split below truncates the operator-space
// MPO, not a positivity-preserving purification. Therefore the result remains
// an approximation to the density matrix, but it is not guaranteed to be a valid
// density matrix after truncation.
class MPOSimulatorImpl : public MPOSimulatorBase,
                         public StructuredGateOperations<MPOSimulatorImpl, MPOSimulatorInterface::MatrixClass>
{
  public:
    friend class MPOSimulator;

    MPOSimulatorImpl(size_t N, unsigned int addseed = 0) : MPOSimulatorBase(N, addseed)
    {
#ifdef USE_FAST_SVD
        SVD.setSwitchSize(blockSizeLimit); // lower sizes will use Jacobi
#endif
    }

    void ApplyGate(const Gates::AppliedGate<MatrixClass> &gate) override
    {
        if (const auto *compact = gate.getStructuredGate())
            ApplyStructuredGate(*compact);
        else if (gate.getQubitsNumber() <= 2)
            ApplyOperator(gate);
        else
            ApplyGate(gate, gate.getQubit1(), gate.getQubit2(), gate.getQubit3());
    }

    // A compact gate. A multi-controlled one-qubit operator is applied in place for any number
    // of controls at any sites (see ApplyMultiControlledSites). A diagonal or multiplexed gate
    // on at most three qubits is expanded to its matrix, on sites that must then be consecutive
    // as for every local operator here; on more qubits it is applied in place as the chain
    // operator of its table, whose cost is set by how strongly the table couples the chain.
    void ApplyStructuredGate(const Gates::StructuredGate &gate)
    {
        gate.Validate(gammas.size());
        if (gate.kind == Gates::StructuredGate::Kind::Controlled && !gate.qubits.empty())
        {
            MatrixClass u(2, 2);
            u << gate.matrix[0], gate.matrix[1], gate.matrix[2], gate.matrix[3];
            std::vector<IndexType> controls;
            std::vector<bool> values;
            for (size_t k = 0; k < gate.qubits.size(); ++k)
            {
                controls.push_back(static_cast<IndexType>(gate.qubits[k]));
                values.push_back(gate.controlState == std::numeric_limits<size_t>::max() ||
                                 ((gate.controlState >> k) & 1));
            }
            ApplyMultiControlledSites(u, static_cast<IndexType>(gate.target), controls, values);
            return;
        }
        if (gate.Arity() > 3)
        {
            // A wide diagonal or multiplexed table is applied in place as the chain operator
            // of the table (see ApplyChainOperator).
            ApplyChainOperator(WideTableOperator(gate, [](size_t qubit) { return static_cast<IndexType>(qubit); }),
                               WideTableIsUnitary(gate));
            return;
        }
        const auto qubits = gate.OrderedQubits();
        if (qubits.empty())
        {
            // A table of one entry is a scalar: rho is multiplied by its squared modulus.
            ApplyOperator(Gates::AppliedGate<MatrixClass>(MatrixClass(gate.diagonal[0] * MatrixClass::Identity(2, 2)), 0));
            return;
        }
        ApplyGate(Gates::AppliedGate<MatrixClass>(gate.DenseMatrix<MatrixClass>(3), qubits[0],
                                                  qubits.size() > 1 ? qubits[1] : 0,
                                                  qubits.size() > 2 ? qubits[2] : 0));
    }

    // The largest rank along the chain that a wide diagonal or multiplexed table may have
    // here. The operator acts on the ket and on the bra side, so between two sites the tensors
    // grow by the square of that rank before they are cut again.
    static constexpr IndexType maxChainOperatorRank = 16;

    // The chain operator of a diagonal or multiplexed gate on more than three qubits (see
    // TableTrain.h); siteOf gives the site of a qubit.
    template <class SiteOf>
    static ChainOperator WideTableOperator(const Gates::StructuredGate &gate, const SiteOf &siteOf)
    {
        return MakeChainOperator(gate, siteOf, maxChainOperatorRank);
    }

    // A multiplexed rotation is unitary; a diagonal is when every entry has modulus one.
    static bool WideTableIsUnitary(const Gates::StructuredGate &gate)
    {
        if (gate.kind == Gates::StructuredGate::Kind::Diagonal)
            for (const auto &value : gate.diagonal)
                if (!(std::abs(std::norm(value) - 1.) <= 1E-12))
                    return false;
        return true;
    }

    // rho -> O rho O^dagger for an operator given as a chain of site matrices (see
    // TableTrain.h), from its first site to its last. Each tensor takes the matrices of its
    // site on the ket leg and their conjugates on the bra leg, which multiplies its two bonds
    // by the squared ranks of the operator; the bonds are then cut again as after a
    // multi-controlled gate.
    void ApplyChainOperator(const ChainOperator &chain, bool unitary)
    {
        const IndexType first = chain.first, last = first + static_cast<IndexType>(chain.sites.size()) - 1;
        if (chain.sites.empty() || first < 0 || last >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");
        InvalidateSamplingCache();
        for (IndexType site = first; site <= last; ++site)
        {
            const ChainOperator::Site &operators = chain.sites[static_cast<size_t>(site - first)];
            IndexType leftRank = 0, rightRank = 0;
            for (const auto &part : operators)
                if (part.size())
                {
                    leftRank = part.rows();
                    rightRank = part.cols();
                }
            const TensorType &gamma = gammas[site];
            const IndexType L = gamma.dimension(0), R = gamma.dimension(3);
            TensorType product(L * leftRank * leftRank, 2, 2, R * rightRank * rightRank);
            product.setZero();
            for (IndexType ketOut = 0; ketOut < 2; ++ketOut)
                for (IndexType ketIn = 0; ketIn < 2; ++ketIn)
                {
                    const auto &ketPart = operators[static_cast<size_t>(2 * ketOut + ketIn)];
                    if (!ketPart.size())
                        continue;
                    for (IndexType braOut = 0; braOut < 2; ++braOut)
                        for (IndexType braIn = 0; braIn < 2; ++braIn)
                        {
                            const auto &braPart = operators[static_cast<size_t>(2 * braOut + braIn)];
                            if (!braPart.size())
                                continue;
                            for (IndexType c = 0; c < rightRank; ++c)
                                for (IndexType c2 = 0; c2 < rightRank; ++c2)
                                    for (IndexType a = 0; a < leftRank; ++a)
                                        for (IndexType a2 = 0; a2 < leftRank; ++a2)
                                        {
                                            const std::complex<double> factor = ketPart(a, c) * std::conj(braPart(a2, c2));
                                            if (factor == std::complex<double>(0.))
                                                continue;
                                            const IndexType rowOffset = (a * leftRank + a2) * L,
                                                            colOffset = (c * rightRank + c2) * R;
                                            for (IndexType r = 0; r < R; ++r)
                                                for (IndexType l = 0; l < L; ++l)
                                                    product(rowOffset + l, ketOut, braOut, colOffset + r) +=
                                                        factor * gamma(l, ketIn, braIn, r);
                                        }
                        }
                }
            gammas[site] = std::move(product);
            if (site < last)
            {
                const IndexType cols = gammas[site].dimension(3);
                lambdas[site] = LambdaType::Ones(cols);
                bondDimensions.Update(site, cols);
            }
        }
        // As after a multi-controlled gate: a unitary keeps a canonical chain canonical.
        if (canonicalFormValid && unitary)
            RunMaybeSingleThreaded(enableMultithreading, [&]() {
                for (IndexType site = last; site > first; --site)
                    RightCanonicalizeSite(site);
            });
        else
            InvalidateCanonicalForm(first, last);
        bool truncated = false;
        for (IndexType bond = first; bond < last; ++bond)
        {
            const bool canonicalWeights = canonicalFormValid;
            if (!canonicalWeights)
                PrepareBond(bond);
            const MatrixClass theta = ReshapeTheta(ContractTwoQubits(bond));
            truncated = DecomposeAndSetGammas(theta, bond, bond + 1, true, canonicalWeights) || truncated;
        }
        ApplyPostTruncationPatches(truncated);
    }

    // Two- and three-qubit gates need to act on consecutive sites.
    // (the higher level MPOSimulator swaps qubits around to satisfy this)
    void ApplyGate(const GateClass &gate, IndexType qubit, IndexType controllingQubit1 = 0) override
    {
        ApplyOperator(gate, qubit, controllingQubit1);
    }

    void ApplyGate(const GateClass &gate, IndexType qubit, IndexType controllingQubit1,
                   IndexType controllingQubit2) override
    {
        if (gate.getQubitsNumber() == 3)
        {
            ThreeSite::Validate(gate.getRawOperatorMatrix(), static_cast<IndexType>(gammas.size()), qubit,
                                controllingQubit1, controllingQubit2, true);
            ApplyThreeQubitGate(gate.getRawOperatorMatrix(), qubit, controllingQubit1, controllingQubit2);
            return;
        }
        ApplyGate(gate, qubit, controllingQubit1);
    }

    void ApplyGates(const std::vector<Gates::AppliedGate<MatrixClass>> &gates) override
    {
        for (const auto &gate : gates)
            ApplyGate(gate);
    }

    // A one-qubit operator on `target`, applied where every control holds its value:
    // O = 1 + (product of the control projectors) (x) (U - 1), and rho -> O rho O^dagger.
    // The two terms make O an operator of bond dimension two along the chain, whatever the
    // number of controls and wherever they sit, so it is applied in place: no swaps, no
    // ancillas and no 2^k matrix. Between the first and the last of its sites each tensor
    // becomes the direct sum of what the four ket/bra combinations of the two terms make of
    // it, and one sweep brings the bonds down again, with the limits the user has set.
    // The sites are positions in the chain; U need not be unitary.
    void ApplyMultiControlledSites(const MatrixClass &u, IndexType target, const std::vector<IndexType> &controls,
                                   const std::vector<bool> &controlValues)
    {
        const IndexType sites = static_cast<IndexType>(gammas.size());
        if (u.rows() != 2 || u.cols() != 2 || !u.allFinite())
            throw std::invalid_argument("A multi-controlled target matrix must be a finite 2x2 matrix");
        if (controls.size() != controlValues.size())
            throw std::invalid_argument("Each control needs its value");
        if (target < 0 || target >= sites)
            throw std::invalid_argument("Qubit index out of bounds");
        // -1: the operator is the identity there; 0 or 1: a control with that value; 2: the target
        std::vector<int> role(static_cast<size_t>(sites), -1);
        role[target] = 2;
        IndexType first = target, last = target;
        for (size_t k = 0; k < controls.size(); ++k)
        {
            const IndexType control = controls[k];
            if (control < 0 || control >= sites)
                throw std::invalid_argument("Qubit index out of bounds");
            if (role[control] != -1)
                throw std::invalid_argument("Controls must be distinct and differ from the target");
            role[control] = controlValues[k] ? 1 : 0;
            first = std::min(first, control);
            last = std::max(last, control);
        }
        if (controls.empty())
        {
            ApplyValidatedOperator(u, 1, target, 0);
            return;
        }
        if (u.isIdentity(0.))
            return;

        InvalidateSamplingCache();
        const MatrixClass difference = u - MatrixClass::Identity(2, 2);
        const bool unitary = (u.adjoint() * u).isIdentity(1E-12);
        for (IndexType site = first; site <= last; ++site)
        {
            const TensorType &gamma = gammas[site];
            const IndexType L = gamma.dimension(0), R = gamma.dimension(3);
            // What the second term makes of the ket leg, of the bra leg, and of both;
            // where the operator is the identity, all of them are the tensor itself.
            TensorType ket, bra, both;
            const TensorType *blocks[4] = {&gamma, &gamma, &gamma, &gamma};
            if (role[site] != -1)
            {
                ket = TensorType(L, 2, 2, R);
                bra = TensorType(L, 2, 2, R);
                both = TensorType(L, 2, 2, R);
                blocks[1] = &ket;
                blocks[2] = &bra;
                blocks[3] = &both;
            }
            if (role[site] == 2)
            {
                ket.setZero();
                bra.setZero();
                both.setZero();
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType q = 0; q < 2; ++q)
                        for (IndexType k = 0; k < 2; ++k)
                            for (IndexType l = 0; l < L; ++l)
                            {
                                const std::complex<double> value = gamma(l, k, q, r);
                                for (IndexType out = 0; out < 2; ++out)
                                {
                                    ket(l, out, q, r) += difference(out, k) * value;
                                    bra(l, k, out, r) += std::conj(difference(out, q)) * value;
                                    for (IndexType outBra = 0; outBra < 2; ++outBra)
                                        both(l, out, outBra, r) +=
                                            difference(out, k) * std::conj(difference(outBra, q)) * value;
                                }
                            }
            }
            else if (role[site] != -1)
            {
                const IndexType selected = role[site];
                ket.setZero();
                bra.setZero();
                both.setZero();
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                    {
                        for (IndexType q = 0; q < 2; ++q)
                            ket(l, selected, q, r) = gamma(l, selected, q, r);
                        for (IndexType k = 0; k < 2; ++k)
                            bra(l, k, selected, r) = gamma(l, k, selected, r);
                        both(l, selected, selected, r) = gamma(l, selected, selected, r);
                    }
            }
            const IndexType rows = site == first ? L : 4 * L, cols = site == last ? R : 4 * R;
            TensorType sum(rows, 2, 2, cols);
            sum.setZero();
            for (IndexType block = 0; block < 4; ++block)
                CopyGammaBlock(*blocks[block], sum, site == first ? 0 : block * L, site == last ? 0 : block * R);
            gammas[site] = std::move(sum);
            if (site < last)
            {
                lambdas[site] = LambdaType::Ones(cols);
                bondDimensions.Update(site, cols);
            }
        }

        // A unitary leaves the weights of the bonds outside its sites as they were, so a chain
        // on its canonical form stays on it: the sites after the first are made orthonormal
        // from the right, then every bond is cut with its left lambda as usual. Otherwise the
        // sites join the window that the sweep below moves along.
        if (canonicalFormValid && unitary)
            RunMaybeSingleThreaded(enableMultithreading, [&]() {
                for (IndexType site = last; site > first; --site)
                    RightCanonicalizeSite(site);
            });
        else
            InvalidateCanonicalForm(first, last);
        bool truncated = false;
        for (IndexType bond = first; bond < last; ++bond)
        {
            const bool canonicalWeights = canonicalFormValid;
            if (!canonicalWeights)
                PrepareBond(bond);
            const MatrixClass theta = ReshapeTheta(ContractTwoQubits(bond));
            truncated = DecomposeAndSetGammas(theta, bond, bond + 1, true, canonicalWeights) || truncated;
        }
        ApplyPostTruncationPatches(truncated);
    }

    // Compresses the current state using either or both configured limits, without a gate.
    void Trim() override
    {
        ApplyPostTruncationPatches(Compress());
    }

    void ReCanonicalize() override
    {
        InvalidateCanonicalForm();
        const IndexType nrBonds = static_cast<IndexType>(lambdas.size());

        // QR transfers every component to the left neighbor, without consulting stale
        // lambdas or dropping small directions that a large left environment may amplify.
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            for (IndexType site = nrBonds; site > 0; --site)
                RightCanonicalizeSite(site);
        });

        // With orthonormal right environments, the left-to-right SVDs now see the
        // actual operator Schmidt weights. Only numerical rank filtering is allowed.
        for (IndexType bond = 0; bond < nrBonds; ++bond)
        {
            const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(bond);
            const MatrixClass thetaMatrix = ReshapeTheta(theta);
            DecomposeAndSetGammas(thetaMatrix, bond, bond + 1, false);
        }
        canonicalFormValid = true;
        centerFirst = centerLast = 0;
    }

    void Hermitize() override
    {
        const bool wasApplying = applyingPostTruncationPatches;
        ScopedFlag patchGuard(applyingPostTruncationPatches);
        InvalidateCanonicalForm();

        // the operator is the product of the site tensors, so the adjoint is the product of the adjoint sites
        std::vector<TensorType> adjointGammas;
        adjointGammas.reserve(gammas.size());
        for (const auto &gamma : gammas)
            adjointGammas.emplace_back(AdjointSite(gamma));

        AddState(lambdas, gammas, adjointGammas);
        bondDimensions.Refresh(lambdas);
        ScaleSite(0, 0.5);

        ReCanonicalize();
        // The patches are suspended here, so the limits are applied without them. The trace is the one of
        // (rho + rho^dagger) / 2 unless this cut something: only then is there a trace to restore, as after
        // a gate. An operator that is deliberately not normalized keeps its trace otherwise.
        const bool truncated = Compress();

        if (truncated && restoreTraceAfterTruncation && !wasApplying)
            RestoreTraceIfSafe();
    }

    void ApplyOperator(const Gates::AppliedGate<MatrixClass> &op) override
    {
        if (const auto *compact = op.getStructuredGate())
        {
            ApplyStructuredGate(*compact);
            return;
        }
        const size_t operatorQubits = op.getQubitsNumber();
        if (operatorQubits == 3)
        {
            ApplyOperator(op, CheckedThreeSiteQubit(op.getQubit1()), CheckedThreeSiteQubit(op.getQubit2()),
                          CheckedThreeSiteQubit(op.getQubit3()));
            return;
        }
        ValidateOperator(op, operatorQubits);
        ValidateQubits(operatorQubits, op.getQubit1(), op.getQubit2());

        ApplyValidatedOperator(op.getRawOperatorMatrix(), operatorQubits, static_cast<IndexType>(op.getQubit1()),
                               operatorQubits > 1 ? static_cast<IndexType>(op.getQubit2()) : 0);
    }

    // Multi-site operators need to act on adjacent qubits
    // (the higher level MPOSimulator swaps qubits around to satisfy this)
    void ApplyOperator(const GateClass &op, IndexType qubit, IndexType controllingQubit1 = 0) override
    {
        const size_t operatorQubits = ValidateOperator(op);
        ValidateQubits(operatorQubits, qubit, controllingQubit1);

        ApplyValidatedOperator(op.getRawOperatorMatrix(), operatorQubits, qubit, controllingQubit1);
    }

    void ApplyOperators(const std::vector<Gates::AppliedGate<MatrixClass>> &ops) override
    {
        for (const auto &op : ops)
            ApplyOperator(op);
    }

    void ApplyOperator(const GateClass &op, IndexType qubit, IndexType controllingQubit1,
                       IndexType controllingQubit2) override
    {
        // ApplyGate's three-site path already implements A rho A^dagger and
        // validates the finite 8x8 matrix and all three consecutive sites.
        ApplyGate(op, qubit, controllingQubit1, controllingQubit2);
    }

    void ApplyOperatorAndNormalize(const Gates::AppliedGate<MatrixClass> &op) override
    {
        const size_t operatorQubits = op.getQubitsNumber();
        if (operatorQubits == 3)
        {
            ApplyOperatorAndNormalize(op, CheckedThreeSiteQubit(op.getQubit1()), CheckedThreeSiteQubit(op.getQubit2()),
                                      CheckedThreeSiteQubit(op.getQubit3()));
            return;
        }
        ValidateOperator(op, operatorQubits);
        ValidateQubits(operatorQubits, op.getQubit1(), op.getQubit2());

        ApplyValidatedOperatorAndNormalize(op.getRawOperatorMatrix(), operatorQubits,
                                           static_cast<IndexType>(op.getQubit1()),
                                           operatorQubits > 1 ? static_cast<IndexType>(op.getQubit2()) : 0);
    }

    void ApplyOperatorAndNormalize(const GateClass &op, IndexType qubit, IndexType controllingQubit1 = 0) override
    {
        const size_t operatorQubits = ValidateOperator(op);
        ValidateQubits(operatorQubits, qubit, controllingQubit1);

        ApplyValidatedOperatorAndNormalize(op.getRawOperatorMatrix(), operatorQubits, qubit, controllingQubit1);
    }

    void ApplyOperatorAndNormalize(const GateClass &op, IndexType qubit, IndexType controllingQubit1,
                                   IndexType controllingQubit2) override
    {
        if (op.getQubitsNumber() != 3)
        {
            ApplyOperatorAndNormalize(op, qubit, controllingQubit1);
            return;
        }
        const auto &matrix = op.getRawOperatorMatrix();
        ThreeSite::Validate(matrix, static_cast<IndexType>(gammas.size()), qubit, controllingQubit1, controllingQubit2,
                            true);
        ApplyValidatedOperatorAndNormalize(matrix, 3, qubit, controllingQubit1, controllingQubit2);
    }

    void ApplyKrausOperators(const std::vector<Gates::AppliedGate<MatrixClass>> &ops) override
    {
        if (ops.empty())
            return;

        const size_t qubit = ops.front().getQubit1();
        const size_t controllingQubit1 = ops.front().getQubit2();
        const size_t operatorQubits = ValidateOperator(ops.front());

        // Validate the complete channel before changing any tensor.
        for (const auto &op : ops)
        {
            if (ValidateOperator(op) != operatorQubits)
                throw std::invalid_argument("All Kraus operators need to have the same arity");
        }

        ValidateQubits(operatorQubits, qubit, controllingQubit1);
        for (const auto &op : ops)
            if (op.getQubit1() != qubit || (operatorQubits > 1 && op.getQubit2() != controllingQubit1))
                throw std::invalid_argument("All Kraus operators need to act on the same qubits");

        ApplyKrausOperatorsImpl(ops, operatorQubits, static_cast<IndexType>(qubit),
                                operatorQubits > 1 ? static_cast<IndexType>(controllingQubit1) : 0);
    }

    void ApplyKrausOperators(const std::vector<MatrixClass> &ops, IndexType qubit,
                             IndexType controllingQubit1 = 0) override
    {
        if (ops.empty())
            return;

        const size_t operatorQubits = ValidateOperatorMatrix(ops.front());
        for (const auto &op : ops)
            if (ValidateOperatorMatrix(op) != operatorQubits)
                throw std::invalid_argument("All Kraus operators need to have the same arity");

        ValidateQubits(operatorQubits, qubit, controllingQubit1);
        ApplyKrausOperatorsImpl(ops, operatorQubits, qubit, controllingQubit1);
    }

    // A linear map on the operator space of one site, or of two neighbouring sites, given as a superoperator:
    // the 4x4 or 16x16 matrix that multiplies the local density blocks indexed ket + d * bra, with `qubit` the
    // low bit of both and `controllingQubit1` the high one. For Kraus operators E_k it is
    // sum_k conj(E_k) (x) E_k, and for a unitary U it is conj(U) (x) U. Maps applied one after another multiply
    // (later ones on the left), so a gate, its noise and the gates around them on the same sites cost one
    // update and one decomposition, where a gate followed by a two-qubit channel costs two. The matrix is
    // trusted: only its shape and its entries being finite are checked.
    void ApplySuperoperator(const MatrixClass &superoperator, IndexType qubit, IndexType controllingQubit1 = 0)
    {
        const size_t operatorQubits = superoperator.rows() == 4 ? 1 : superoperator.rows() == 16 ? 2 : 0;
        if (!operatorQubits || superoperator.cols() != superoperator.rows() || !superoperator.allFinite())
            throw std::invalid_argument("A superoperator must be a finite 4x4 or 16x16 matrix");
        ValidateQubits(operatorQubits, qubit, controllingQubit1);
        if (operatorQubits == 1)
            ApplySingleSiteSuperoperator(superoperator, qubit);
        else
            ApplyTwoSiteSuperoperator(superoperator, qubit, controllingQubit1);
    }

    // false if measured 0, true if measured 1
    bool MeasureQubit(IndexType qubit) override
    {
        if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
            throw std::invalid_argument("Qubit index out of bounds");

        const double prob0 = ResolvedProbability(ValidMeasurementProbability(GetProbability(qubit, true)));
        // the number is drawn whatever the probability, so that the sequence of a seed does not depend on it
        const double drawn = uniformZeroOne(rng);
        const bool zeroMeasured = prob0 >= 1. || drawn < prob0;

        // for a density matrix, the post measurement state is P rho P (suitably normalized)
        // P acts only on the measured site, so this is a purely local operation - no need to
        // propagate the collapse along the chain like the MPS simulator does
        const MatrixClass &projMat =
            zeroMeasured ? zeroProjection.getRawOperatorMatrix() : oneProjection.getRawOperatorMatrix();
        TensorType collapsedGamma = gammas[qubit];
        ApplySingleQubitGate(collapsedGamma, projMat);

        // Normalize by the actual post-projection trace, rather than by the probability
        // normalized against the input trace. This also handles deliberately unnormalized MPOs:
        // the outcome taken has a probability above the rounding, so its trace is not compared
        // with a fixed bound, which would refuse ordinary outcomes of an operator of small trace.
        // Only its size is checked: a truncated operator that is no longer Hermitian has a
        // complex trace here, and is divided by it as before.
        const std::complex<double> postTrace = ContractChain([this, qubit, &collapsedGamma](IndexType q) {
            return q == qubit ? TraceSiteTensor(collapsedGamma) : SiteTraceMatrix(q);
        });
        if (!IsFinite(postTrace) || !(std::abs(postTrace) >= std::numeric_limits<double>::min()))
            throw std::runtime_error("Cannot collapse an MPO state onto a zero-probability outcome");

        ScaleTensor(collapsedGamma, 1. / postTrace);
        gammas[qubit] = std::move(collapsedGamma);
        InvalidateCanonicalForm(qubit, qubit);

        return !zeroMeasured;
    }

    std::unordered_map<IndexType, bool> MeasureQubits(const std::set<IndexType> &qubits) override
    {
        for (const IndexType qubit : qubits)
            if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
                throw std::invalid_argument("Qubit index out of bounds");

        std::unordered_map<IndexType, bool> res;
        for (const auto qubit : qubits)
            res[qubit] = MeasureQubit(qubit);

        return res;
    }

  private:
    class ScopedFlag
    {
        bool &flag;
        const bool previous;

      public:
        explicit ScopedFlag(bool &value) : flag(value), previous(value)
        {
            flag = true;
        }

        ~ScopedFlag()
        {
            flag = previous;
        }

        ScopedFlag(const ScopedFlag &) = delete;
        ScopedFlag &operator=(const ScopedFlag &) = delete;
    };

    void PrepareBond(IndexType bond, IndexType last = -1)
    {
        if (last < 0)
            last = bond + 1;
        InvalidateSamplingCache();
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            while (centerFirst < bond)
            {
                LeftCanonicalizeSite(centerFirst);
                ++centerFirst;
                centerLast = std::max(centerLast, centerFirst);
            }
            while (centerLast > last)
            {
                RightCanonicalizeSite(centerLast);
                --centerLast;
                centerFirst = std::min(centerFirst, centerLast);
            }
        });
    }

    void LeftCanonicalizeSite(IndexType site)
    {
        const auto &gamma = gammas[site];
        const IndexType L = gamma.dimension(0), R = gamma.dimension(3);
        const IndexType rank = std::min(4 * L, R);
        const Eigen::Map<const MatrixClass> matrix(gamma.data(), 4 * L, R);
        Eigen::HouseholderQR<MatrixClass> qr(matrix);
        const MatrixClass Q = qr.householderQ() * MatrixClass::Identity(4 * L, rank);
        MatrixClass transfer = qr.matrixQR().topRows(rank).template triangularView<Eigen::Upper>();
        // R carries the norm of the site into the next one: keep it in range (see scaleExponent)
        int64_t removed = 0;
        RescaleIfOutOfRange(transfer, removed);
        const auto &next = gammas[site + 1];
        const IndexType nextR = next.dimension(3);
        const Eigen::Map<const MatrixClass> nextMatrix(next.data(), R, 4 * nextR);
        const MatrixClass right = transfer * nextMatrix;
        TensorType newLeft(L, 2, 2, rank), newRight(rank, 2, 2, nextR);
        std::copy(Q.data(), Q.data() + Q.size(), newLeft.data());
        std::copy(right.data(), right.data() + right.size(), newRight.data());
        LambdaType placeholder = LambdaType::Ones(rank);
        gammas[site] = std::move(newLeft);
        gammas[site + 1] = std::move(newRight);
        lambdas[site] = std::move(placeholder);
        bondDimensions.Update(site, lambdas[site].size());
        scaleExponent += removed;
    }

    // Exact physical reordering for an overlap copy, without gate contractions,
    // user compression, trace patches, callbacks or a global SVD sweep.
    void SwapSitesForOverlap(IndexType first)
    {
        if (canonicalFormValid)
            InvalidateCanonicalForm(0, 0);
        PrepareBond(first);
        if (TrySwapProductSites(first, false))
            return;
        const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(first);
        DecomposeAndSetGammas(ReshapeSwappedTheta(theta), first, first + 1, false, false);
    }

    void RightCanonicalizeSite(IndexType site)
    {
        const TensorType &gamma = gammas[site];
        const IndexType L = gamma.dimension(0);
        const IndexType R = gamma.dimension(3);
        const IndexType rank = std::min(L, 4 * R);
        MatrixClass adjoint(4 * R, L);
        for (IndexType ket = 0; ket < 2; ++ket)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                        adjoint((ket * 2 + bra) * R + r, l) = std::conj(gamma(l, ket, bra, r));

        Eigen::HouseholderQR<MatrixClass> qr(adjoint);
        const MatrixClass Q = qr.householderQ() * MatrixClass::Identity(4 * R, rank);
        const MatrixClass triangular = qr.matrixQR().topRows(rank).template triangularView<Eigen::Upper>();
        MatrixClass transfer = triangular.adjoint();
        // R^dagger carries the norm of the site into the previous one: keep it in range
        int64_t removed = 0;
        RescaleIfOutOfRange(transfer, removed);
        const TensorType &previous = gammas[site - 1];
        const IndexType previousL = previous.dimension(0);
        // Tensor storage groups (left, ket, bra) contiguously before the right bond.
        const Eigen::Map<const MatrixClass> previousMatrix(previous.data(), 4 * previousL, L);
        const MatrixClass left = previousMatrix * transfer;
        TensorType newLeft(previousL, 2, 2, rank);
        std::copy(left.data(), left.data() + left.size(), newLeft.data());
        TensorType newRight(rank, 2, 2, R);
        for (IndexType ket = 0; ket < 2; ++ket)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < rank; ++l)
                        newRight(l, ket, bra, r) = std::conj(Q((ket * 2 + bra) * R + r, l));
        LambdaType placeholder = LambdaType::Ones(rank);
        gammas[site - 1] = std::move(newLeft);
        gammas[site] = std::move(newRight);
        lambdas[site - 1] = std::move(placeholder);
        bondDimensions.Update(site - 1, lambdas[site - 1].size());
        scaleExponent += removed;
    }

    static bool IsFinite(const std::complex<double> &value)
    {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    }

    static IndexType CheckedThreeSiteQubit(size_t qubit)
    {
        if (qubit > static_cast<size_t>(std::numeric_limits<IndexType>::max()))
            throw std::invalid_argument("Qubit index out of bounds");
        return static_cast<IndexType>(qubit);
    }

    static size_t ValidateOperatorMatrix(const MatrixClass &matrix)
    {
        if (matrix.rows() <= 0 || matrix.cols() <= 0 || matrix.rows() != matrix.cols())
            throw std::invalid_argument("Operator matrix must be a non-empty square matrix");
        if (matrix.rows() != 2 && matrix.rows() != 4)
            throw std::invalid_argument("MPO operators must have dimensions 2x2 or 4x4");

        for (IndexType col = 0; col < matrix.cols(); ++col)
            for (IndexType row = 0; row < matrix.rows(); ++row)
                if (!IsFinite(matrix(row, col)))
                    throw std::invalid_argument("Operator matrix elements must be finite");

        return matrix.rows() == 2 ? 1 : 2;
    }

    static size_t ValidateOperator(const GateClass &op, size_t knownArity = 0)
    {
        const size_t matrixQubits = ValidateOperatorMatrix(op.getRawOperatorMatrix());
        if ((knownArity ? knownArity : op.getQubitsNumber()) != matrixQubits)
            throw std::invalid_argument("Operator matrix dimensions do not match its declared arity");
        return matrixQubits;
    }

    void ValidateQubits(size_t operatorQubits, IndexType qubit, IndexType controllingQubit1) const
    {
        const IndexType nrQubits = static_cast<IndexType>(gammas.size());
        if (qubit < 0 || qubit >= nrQubits ||
            (operatorQubits == 2 && (controllingQubit1 < 0 || controllingQubit1 >= nrQubits)))
            throw std::invalid_argument("Qubit index out of bounds");
        if (operatorQubits == 2 && qubit == controllingQubit1)
            throw std::invalid_argument("Two qubit operators require distinct qubits");
        if (operatorQubits == 2 && std::abs(qubit - controllingQubit1) != 1)
            throw std::invalid_argument("Two qubit operators need to act on adjacent qubits");
    }

    void ValidateQubits(size_t operatorQubits, size_t qubit, size_t controllingQubit1) const
    {
        const size_t nrQubits = gammas.size();
        if (qubit >= nrQubits || (operatorQubits == 2 && controllingQubit1 >= nrQubits))
            throw std::invalid_argument("Qubit index out of bounds");
        if (operatorQubits == 2 && qubit == controllingQubit1)
            throw std::invalid_argument("Two qubit operators require distinct qubits");
        if (operatorQubits == 2 && qubit + 1 != controllingQubit1 && controllingQubit1 + 1 != qubit)
            throw std::invalid_argument("Two qubit operators need to act on adjacent qubits");
    }

    static const MatrixClass &GetOperatorMatrix(const MatrixClass &op)
    {
        return op;
    }

    static const MatrixClass &GetOperatorMatrix(const GateClass &op)
    {
        return op.getRawOperatorMatrix();
    }

    static MatrixClass TraceSiteTensor(const TensorType &gamma)
    {
        const IndexType L = gamma.dimension(0);
        const IndexType R = gamma.dimension(3);
        MatrixClass result = MatrixClass::Zero(L, R);

        for (IndexType state = 0; state < 2; ++state)
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < L; ++l)
                    result(l, r) += gamma(l, state, state, r);

        return result;
    }

    static void ScaleTensor(TensorType &gamma, std::complex<double> factor)
    {
        for (IndexType r = 0; r < gamma.dimension(3); ++r)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType ket = 0; ket < 2; ++ket)
                    for (IndexType l = 0; l < gamma.dimension(0); ++l)
                        gamma(l, ket, bra, r) *= factor;
    }

    void ApplyValidatedOperator(const MatrixClass &op, size_t operatorQubits, IndexType qubit,
                                IndexType controllingQubit1)
    {
        InvalidateSamplingCache();
        if (operatorQubits == 1)
        {
            if (limitSize || limitEntanglement)
                ApplySingleQubitGateConservative(gammas[qubit], op);
            else
                ApplySingleQubitGate(gammas[qubit], op);
        }
        else
            ApplyTwoQubitGate(op, qubit, controllingQubit1);
        if (!(op.adjoint() * op).isIdentity(1E-12))
            InvalidateCanonicalForm(std::min(qubit, operatorQubits == 1 ? qubit : controllingQubit1),
                                    std::max(qubit, operatorQubits == 1 ? qubit : controllingQubit1));
    }

    void ApplyValidatedOperatorAndNormalize(const MatrixClass &op, size_t operatorQubits, IndexType qubit,
                                            IndexType controllingQubit1, IndexType controllingQubit2 = 0)
    {
        if (operatorQubits == 1)
        {
            ApplySingleQubitOperatorAndNormalize(op, qubit);
            return;
        }

        auto originalLambdas = lambdas;
        auto originalGammas = gammas;
        const auto originalMetadata = GetCanonicalMetadata();
        const int64_t originalExponent = scaleExponent;

        try
        {
            ApplyAndNormalizeWithoutBackup(op, operatorQubits, qubit, controllingQubit1, controllingQubit2);
        }
        catch (...)
        {
            lambdas.swap(originalLambdas);
            gammas.swap(originalGammas);
            scaleExponent = originalExponent;
            RestoreCanonicalMetadata(originalMetadata);
            throw;
        }
    }

    void ApplySingleQubitOperatorAndNormalize(const MatrixClass &op, IndexType qubit)
    {
        // Only the target tensor needs staging. Keep the chain and its metadata
        // untouched until the prospective trace has passed normalization checks.
        TensorType updatedGamma = gammas[qubit];
        if (limitSize || limitEntanglement)
            ApplySingleQubitGateConservative(updatedGamma, op);
        else
            ApplySingleQubitGate(updatedGamma, op);

        // The trace the operator leaves is judged against the trace it was applied to.
        std::complex<double> trace, traceBefore;
        if (limitSize || limitEntanglement)
        {
            traceBefore = Trace();
            trace = ContractChain([this, qubit, &updatedGamma](IndexType q) {
                return q == qubit ? TraceSiteTensor(updatedGamma) : SiteTraceMatrix(q);
            });
        }
        else
        {
            // Match Trace's fast contraction, including its summation order and
            // power-of-two scaling, without installing the tentative tensor. The two traces
            // share the chain up to the site, so the one before is carried along from there.
            MatrixClass left = MatrixClass::Ones(1, 1), next, leftBefore, nextBefore;
            int64_t exponent = scaleExponent, exponentBefore = scaleExponent;
            for (IndexType q = 0; q < static_cast<IndexType>(gammas.size()); ++q)
            {
                if (q == qubit)
                {
                    leftBefore = left;
                    exponentBefore = exponent;
                }
                if (q >= qubit)
                {
                    ContractPauliSite(gammas[q], 'I', leftBefore, nextBefore);
                    RescaleIfOutOfRange(nextBefore, exponentBefore);
                    leftBefore.swap(nextBefore);
                }
                ContractPauliSite(q == qubit ? updatedGamma : gammas[q], 'I', left, next);
                RescaleIfOutOfRange(next, exponent);
                left.swap(next);
            }
            trace = ScaleByPowerOfTwo(left(0, 0), exponent);
            traceBefore = ScaleByPowerOfTwo(leftBefore(0, 0), exponentBefore);
        }
        RequireResolvedWeight(trace, traceBefore);
        const bool unitary = (op.adjoint() * op).isIdentity(1E-12);

        // All allocating work is complete. Preserve the existing normalization
        // site and Schmidt-weight updates so canonical environments stay valid.
        gammas[qubit] = std::move(updatedGamma);
        if (!unitary)
            InvalidateCanonicalForm(qubit, qubit);
        ScaleSite(canonicalFormValid ? 0 : centerFirst, 1. / trace);
    }

    void ApplyAndNormalizeWithoutBackup(const MatrixClass &op, size_t operatorQubits, IndexType qubit,
                                        IndexType controllingQubit1, IndexType controllingQubit2 = 0)
    {
        const auto traceBefore = Trace();
        bool truncated = false;
        {
            // The patches wait until the weight is judged: restoring the trace after a truncation would put a
            // trace of one in the place of what the operator left, and that would be compared with the trace
            // before it, refusing an ordinary operator on an operator of large trace and passing a weight that
            // is only rounding.
            ScopedFlag deferred(deferPostTruncationPatches);
            truncatedWhileDeferred = false;
            if (operatorQubits == 3)
                ApplyThreeQubitGate(op, qubit, controllingQubit1, controllingQubit2);
            else
                ApplyValidatedOperator(op, operatorQubits, qubit, controllingQubit1);
            truncated = truncatedWhileDeferred;
        }
        auto trace = Trace();
        RequireResolvedWeight(trace, traceBefore);
        if (truncated)
        {
            // Then as after a gate, and the state is normalized by the trace the patches leave, as it was before.
            ApplyPostTruncationPatches(true);
            trace = Trace();
            if (!IsFinite(trace) || !(trace.real() >= std::numeric_limits<double>::min()))
                throw std::runtime_error("Cannot normalize an MPO state: its truncation left no positive trace");
        }
        ScaleSite(canonicalFormValid ? 0 : centerFirst, 1. / trace);
    }

    template <class OperatorsContainer>
    void ApplyKrausOperatorsImpl(const OperatorsContainer &ops, size_t operatorQubits, IndexType qubit,
                                 IndexType controllingQubit1)
    {
        CheckKrausCompleteness(ops);

        if (operatorQubits == 1)
            ApplySingleQubitChannel(ops, qubit);
        else
            ApplyTwoQubitChannel(ops, qubit, controllingQubit1);
    }

    template <class OperatorsContainer> void CheckKrausCompleteness(const OperatorsContainer &ops) const
    {
        if (krausCompletenessCheck == KrausCompletenessCheck::Ignore)
            return;

        const MatrixClass &first = GetOperatorMatrix(ops.front());
        MatrixClass sum = MatrixClass::Zero(first.rows(), first.cols());
        for (const auto &op : ops)
        {
            const MatrixClass &K = GetOperatorMatrix(op);
            sum.noalias() += K.adjoint() * K;
        }

        const double residual = (sum - MatrixClass::Identity(sum.rows(), sum.cols())).norm();
        if (residual <= krausCompletenessTolerance)
            return;

        if (krausCompletenessCheck == KrausCompletenessCheck::Warn)
        {
            std::cerr << "Kraus operators do not satisfy the completeness relation (residual " << residual << ")"
                      << std::endl;
            return;
        }

        throw std::invalid_argument("Kraus operators do not satisfy the completeness relation");
    }

    template <class OperatorsContainer>
    bool ShouldAggregateChannel(const OperatorsContainer &ops, IndexType L, IndexType R) const
    {
        if (limitSize || limitEntanglement || L * R < 16 || ops.size() < 2)
            return false;
        if (ops.size() > 2)
            return true;
        // Two dense operators also save a full tensor pass. Sparse pairs
        // (for example amplitude damping) keep the zero-skipping loops.
        for (const auto &op : ops)
        {
            const auto &matrix = GetOperatorMatrix(op);
            for (IndexType i = 0; i < matrix.size(); ++i)
                if (matrix.data()[i] == std::complex<double>{})
                    return false;
        }
        return true;
    }

    // Keep superoperator setup and its stack storage out of the fallback
    // channel loops. Inlining this work can penalize even two-Kraus channels.
    template <class OperatorsContainer>
    static EIGEN_DONT_INLINE void AccumulateSingleQubitChannel(const OperatorsContainer &ops,
                                                               const TensorType &original, TensorType &result)
    {
        const IndexType L = original.dimension(0), R = original.dimension(3);
        Eigen::Matrix<std::complex<double>, 4, 4> channel;
        channel.setZero();
        for (const auto &op : ops)
        {
            const auto &E = GetOperatorMatrix(op);
            for (int out = 0; out < 4; ++out)
                for (int in = 0; in < 4; ++in)
                    channel(out, in) += E(out % 2, in % 2) * std::conj(E(out / 2, in / 2));
        }
        for (int out = 0; out < 4; ++out)
            for (int in = 0; in < 4; ++in)
            {
                const auto factor = channel(out, in);
                if (factor == std::complex<double>{})
                    continue;
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                        result(l, out % 2, out / 2, r) += factor * original(l, in % 2, in / 2, r);
            }
    }

    template <class OperatorsContainer>
    static EIGEN_DONT_INLINE void AccumulateTwoQubitChannel(const OperatorsContainer &ops, bool reversed,
                                                            const Eigen::Tensor<std::complex<double>, 6> &theta,
                                                            Eigen::Tensor<std::complex<double>, 6> &result)
    {
        const IndexType L = theta.dimension(0), R = theta.dimension(5);
        Eigen::Matrix<std::complex<double>, 16, 16> channel;
        channel.setZero();
        // out/in = (two-site ket)*4 + (two-site bra).
        for (const auto &op : ops)
        {
            const auto E = GetTwoQubitsGateTensor(GetOperatorMatrix(op), reversed);
            for (int out = 0; out < 16; ++out)
                for (int in = 0; in < 16; ++in)
                    channel(out, in) += E(out / 8, (out / 4) % 2, in / 8, (in / 4) % 2) *
                                        std::conj(E((out / 2) % 2, out % 2, (in / 2) % 2, in % 2));
        }
        for (int out = 0; out < 16; ++out)
            for (int in = 0; in < 16; ++in)
            {
                const auto factor = channel(out, in);
                if (factor == std::complex<double>{})
                    continue;
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                        result(l, r, out / 8, (out / 4) % 2, (out / 2) % 2, out % 2) +=
                            factor * theta(l, in / 8, (in / 2) % 2, (in / 4) % 2, in % 2, r);
            }
    }

    // The site index ket + 2 * bra is the one the channel accumulation above uses.
    void ApplySingleSiteSuperoperator(const MatrixClass &channel, IndexType qubit)
    {
        InvalidateSamplingCache();
        const TensorType &original = gammas[qubit];
        const IndexType L = original.dimension(0), R = original.dimension(3);
        TensorType result(L, 2, 2, R);
        result.setZero();
        for (int out = 0; out < 4; ++out)
            for (int in = 0; in < 4; ++in)
            {
                const std::complex<double> factor = channel(out, in);
                if (factor == std::complex<double>{})
                    continue;
                for (IndexType r = 0; r < R; ++r)
                    for (IndexType l = 0; l < L; ++l)
                        result(l, out % 2, out / 2, r) += factor * original(l, in % 2, in / 2, r);
            }
        gammas[qubit] = std::move(result);
        InvalidateCanonicalForm(qubit, qubit);
    }

    // As a two-qubit channel, with the map given instead of built from Kraus operators.
    void ApplyTwoSiteSuperoperator(const MatrixClass &superoperator, IndexType qubit, IndexType controllingQubit1)
    {
        const bool canonicalWeights = canonicalFormValid;
        if (!canonicalWeights)
            PrepareBond(std::min(qubit, controllingQubit1));
        IndexType qubit1 = controllingQubit1;
        IndexType qubit2 = qubit;
        bool reversed = false;
        if (qubit1 > qubit2)
        {
            std::swap(qubit1, qubit2);
            reversed = true;
        }

        // The accumulation below packs (two-site ket) * 4 + (two-site bra), the left site being the high bit of
        // each. Unless reversed the left site is controllingQubit1, the high bit of the caller's index too.
        const auto callerIndex = [reversed](int packed) {
            const int ket1 = packed / 8, ket2 = (packed / 4) % 2, bra1 = (packed / 2) % 2, bra2 = packed % 2;
            const int ketLow = reversed ? ket1 : ket2, ketHigh = reversed ? ket2 : ket1;
            const int braLow = reversed ? bra1 : bra2, braHigh = reversed ? bra2 : bra1;
            return ketLow + 2 * ketHigh + 4 * (braLow + 2 * braHigh);
        };
        int caller[16];
        for (int packed = 0; packed < 16; ++packed)
            caller[packed] = callerIndex(packed);

        const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
        const IndexType L = theta.dimension(0);
        const IndexType R = theta.dimension(5);
        // One matrix product: the bond indices are the rows, the sixteen states of the pair the columns.
        MatrixClass gathered(L * R, 16), packed(16, 16);
        for (int in = 0; in < 16; ++in)
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < L; ++l)
                    gathered(l + L * r, in) = theta(l, in / 8, (in / 2) % 2, (in / 4) % 2, in % 2, r);
        for (int out = 0; out < 16; ++out)
            for (int in = 0; in < 16; ++in)
                packed(in, out) = superoperator(caller[out], caller[in]);
        // Eigen runs a product of this size on several threads unless told otherwise.
        MatrixClass mapped;
        RunMaybeSingleThreaded(enableMultithreading, [&]() { mapped.noalias() = gathered * packed; });
        Eigen::Tensor<std::complex<double>, 6> result(L, R, 2, 2, 2, 2);
        for (int out = 0; out < 16; ++out)
            for (IndexType r = 0; r < R; ++r)
                for (IndexType l = 0; l < L; ++l)
                    result(l, r, out / 8, (out / 4) % 2, (out / 2) % 2, out % 2) = mapped(l + L * r, out);

        const MatrixClass thetaMatrix = ReshapeThetaBar(result);
        ApplyPostTruncationPatches(DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2, true, canonicalWeights));
        InvalidateCanonicalForm(qubit1, qubit2);
    }

    template <class OperatorsContainer> void ApplySingleQubitChannel(const OperatorsContainer &ops, IndexType qubit)
    {
        InvalidateSamplingCache();
        const TensorType &original = gammas[qubit];
        const IndexType L = original.dimension(0);
        const IndexType R = original.dimension(3);
        TensorType result(L, 2, 2, R);
        result.setZero();

        // Small tensors, sparse pairs and user-compressed evolution retain
        // the original arithmetic. Dense pairs can also amortize the setup.
        if (ShouldAggregateChannel(ops, L, R))
        {
            AccumulateSingleQubitChannel(ops, original, result);
            gammas[qubit] = std::move(result);
            InvalidateCanonicalForm(qubit, qubit);
            return;
        }

        for (const auto &op : ops)
        {
            const MatrixClass &E = GetOperatorMatrix(op);
            for (IndexType ketOut = 0; ketOut < 2; ++ketOut)
                for (IndexType braOut = 0; braOut < 2; ++braOut)
                    for (IndexType ketIn = 0; ketIn < 2; ++ketIn)
                        for (IndexType braIn = 0; braIn < 2; ++braIn)
                        {
                            const std::complex<double> factor = E(ketOut, ketIn) * std::conj(E(braOut, braIn));
                            if (factor == std::complex<double>(0., 0.))
                                continue;

                            for (IndexType r = 0; r < R; ++r)
                                for (IndexType l = 0; l < L; ++l)
                                    result(l, ketOut, braOut, r) += factor * original(l, ketIn, braIn, r);
                        }
        }

        gammas[qubit] = std::move(result);
        InvalidateCanonicalForm(qubit, qubit);
    }

    template <class OperatorsContainer>
    void ApplyTwoQubitChannel(const OperatorsContainer &ops, IndexType qubit, IndexType controllingQubit1)
    {
        const bool canonicalWeights = canonicalFormValid;
        if (!canonicalWeights)
            PrepareBond(std::min(qubit, controllingQubit1));
        IndexType qubit1 = controllingQubit1;
        IndexType qubit2 = qubit;
        bool reversed = false;
        if (qubit1 > qubit2)
        {
            std::swap(qubit1, qubit2);
            reversed = true;
        }

        const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
        const IndexType L = theta.dimension(0);
        const IndexType R = theta.dimension(5);
        Eigen::Tensor<std::complex<double>, 6> result(L, R, 2, 2, 2, 2);
        result.setZero();

        if (ShouldAggregateChannel(ops, L, R))
        {
            AccumulateTwoQubitChannel(ops, reversed, theta, result);
            const MatrixClass thetaMatrix = ReshapeThetaBar(result);
            ApplyPostTruncationPatches(DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2, true, canonicalWeights));
            InvalidateCanonicalForm(qubit1, qubit2);
            return;
        }

        for (const auto &op : ops)
        {
            const TwoQubitsGateTensor E = GetTwoQubitsGateTensor(GetOperatorMatrix(op), reversed);
            for (IndexType ket1Out = 0; ket1Out < 2; ++ket1Out)
                for (IndexType ket2Out = 0; ket2Out < 2; ++ket2Out)
                    for (IndexType bra1Out = 0; bra1Out < 2; ++bra1Out)
                        for (IndexType bra2Out = 0; bra2Out < 2; ++bra2Out)
                            for (IndexType ket1In = 0; ket1In < 2; ++ket1In)
                                for (IndexType ket2In = 0; ket2In < 2; ++ket2In)
                                    for (IndexType bra1In = 0; bra1In < 2; ++bra1In)
                                        for (IndexType bra2In = 0; bra2In < 2; ++bra2In)
                                        {
                                            const std::complex<double> factor =
                                                E(ket1Out, ket2Out, ket1In, ket2In) *
                                                std::conj(E(bra1Out, bra2Out, bra1In, bra2In));
                                            if (factor == std::complex<double>(0., 0.))
                                                continue;

                                            for (IndexType r = 0; r < R; ++r)
                                                for (IndexType l = 0; l < L; ++l)
                                                    result(l, r, ket1Out, ket2Out, bra1Out, bra2Out) +=
                                                        factor * theta(l, ket1In, bra1In, ket2In, bra2In, r);
                                        }
        }

        const MatrixClass thetaMatrix = ReshapeThetaBar(result);
        ApplyPostTruncationPatches(DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2, true, canonicalWeights));
        InvalidateCanonicalForm(qubit1, qubit2);
    }

    ThreeSite::Split SplitThreeSite(const MatrixClass &theta, IndexType leftDimension, const LambdaType *leftLambda,
                                    int64_t &removed)
    {
        MatrixClass weighted;
        if (leftLambda)
        {
            // Only relative Schmidt weights matter for V and rank selection.
            // Remove an extreme common scale before multiplying, so a perfectly
            // representable weighted tensor cannot vanish through underflow.
            const int weightExponent = OutOfRangeExponent(leftLambda->cwiseAbs().maxCoeff());
            if (weightExponent)
            {
                LambdaType weights = *leftLambda;
                ThreeSite::ScalePowerOfTwo(weights, -weightExponent);
                weighted = ThreeSite::Weighted(theta, &weights, leftDimension);
            }
            else
                weighted = ThreeSite::Weighted(theta, leftLambda, leftDimension);
        }
        const MatrixClass &input = leftLambda ? weighted : theta;
#ifdef USE_FAST_SVD
        const bool computeWithJacobi = input.rows() < blockSizeLimit && input.cols() < blockSizeLimit;
        if (computeWithJacobi)
            jacobiSVD.compute(input);
        else
            SVD.compute(input);
        if ((computeWithJacobi ? jacobiSVD.info() : SVD.info()) != Eigen::Success)
            throw std::runtime_error("MPO three-site SVD failed");
        const auto &raw = computeWithJacobi ? jacobiSVD.singularValues() : SVD.singularValues();
        const auto &V = computeWithJacobi ? jacobiSVD.matrixV() : SVD.matrixV();
        const IndexType nonzero = computeWithJacobi ? jacobiSVD.nonzeroSingularValues() : SVD.nonzeroSingularValues();
        const double threshold = computeWithJacobi ? jacobiSVD.threshold() : SVD.threshold();
#else
        jacobiSVD.compute(input);
        if (jacobiSVD.info() != Eigen::Success)
            throw std::runtime_error("MPO three-site SVD failed");
        const auto &raw = jacobiSVD.singularValues();
        const auto &V = jacobiSVD.matrixV();
        const IndexType nonzero = jacobiSVD.nonzeroSingularValues();
        const double threshold = jacobiSVD.threshold();
#endif
        const int exponent = OutOfRangeExponent(raw.size() > 0 ? raw[0] : 0.);
        LambdaType values = raw;
        if (exponent)
            ThreeSite::ScalePowerOfTwo(values, -exponent);
        const IndexType rank = NumericalRank(values, nonzero, threshold);
        IndexType kept =
            limitEntanglement ? ComputeCompressedRank(values, rank, truncationMode, singularValueThreshold) : rank;
        kept = std::max<IndexType>(1, limitSize ? std::min(chi, kept) : kept);
        ThreeSite::Split result;
        result.lambda = values.head(kept);
        result.left.noalias() = theta * V.leftCols(kept);
        RescaleIfOutOfRange(result.left, removed);
        result.right = V.leftCols(kept).adjoint();
        result.truncated = kept < rank;
        return result;
    }

    void ApplyThreeQubitGate(const MatrixClass &gate, IndexType target, IndexType control1, IndexType control2)
    {
        const IndexType first = std::min({target, control1, control2});
        InvalidateSamplingCache();
        const bool canonicalWeights = canonicalFormValid;
        if (!canonicalWeights)
            PrepareBond(first, first + 2);
        bool truncated = false;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            const IndexType L = gammas[first].dimension(0), R = gammas[first + 2].dimension(3);
            // rho transforms with two factors of the operator. Keep their common
            // scale out of the contractions, even when its square is outside the
            // floating-point range; the existing MPO exponent stores it exactly.
            int gateExponent = OutOfRangeExponent(MaxAbs(gate));
            const bool extremeGate = gateExponent != 0;
            MatrixClass theta = ThreeSite::Contract<4>(gammas[first], gammas[first + 1], gammas[first + 2]);
            MatrixClass scaledGate;
            if (extremeGate)
            {
                scaledGate = gate;
                ThreeSite::RestrictToInputSupport<4>(scaledGate, theta, first, target, control1, control2, L, R);
                gateExponent = OutOfRangeExponent(MaxAbs(scaledGate));
                ThreeSite::ScalePowerOfTwo(scaledGate, -gateExponent);
            }
            const MatrixClass &operation = extremeGate ? scaledGate : gate;
            int64_t removed = 2 * int64_t(gateExponent);
            RescaleIfOutOfRange(theta, removed);
            ThreeSite::Apply<4>(theta, operation, first, target, control1, control2, L, R);
            RescaleIfOutOfRange(theta, removed);
            const LambdaType *weights = canonicalWeights && first > 0 ? &lambdas[first - 1] : nullptr;
            auto last = SplitThreeSite(theta, L, weights, removed);
            const IndexType lastBond = last.lambda.size();
            const MatrixClass pair = Eigen::Map<const MatrixClass>(last.left.data(), 4 * L, 4 * lastBond);
            auto left = SplitThreeSite(pair, L, weights, removed);
            const IndexType firstBond = left.lambda.size();
            TensorType a(L, 2, 2, firstBond), b(firstBond, 2, 2, lastBond), c(lastBond, 2, 2, R);
            std::copy(left.left.data(), left.left.data() + left.left.size(), a.data());
            std::copy(left.right.data(), left.right.data() + left.right.size(), b.data());
            std::copy(last.right.data(), last.right.data() + last.right.size(), c.data());
            gammas[first] = std::move(a);
            gammas[first + 1] = std::move(b);
            gammas[first + 2] = std::move(c);
            lambdas[first] = std::move(left.lambda);
            lambdas[first + 1] = std::move(last.lambda);
            bondDimensions.Update(first, firstBond);
            bondDimensions.Update(first + 1, lastBond);
            scaleExponent += removed;
            truncated = left.truncated || last.truncated;
            if (!canonicalWeights)
            {
                canonicalFormValid = false;
                centerFirst = centerLast = first;
            }
            else if (truncated)
                InvalidateCanonicalForm(first, first);
            if (extremeGate || !(gate.adjoint() * gate).isIdentity(1E-12))
                InvalidateCanonicalForm(first, first);
        });
        // Both splits must be committed before hermitization/trace restoration.
        ApplyPostTruncationPatches(truncated);
    }

    void ApplyTwoQubitGate(const MatrixClass &gate, IndexType qubit, IndexType controllingQubit1)
    {
        const bool canonicalWeights = canonicalFormValid;
        if (!canonicalWeights)
            PrepareBond(std::min(qubit, controllingQubit1));
        // contract the tensors for the two qubits (the B tensors already contain the lambdas, except the left one)
        // apply the gate on the kets and its conjugate on the bras (rho -> U rho U^dagger)
        // SVD the result (with the left lambda multiplied in) to separate the two qubit tensors and the new lambda in
        // between. If the SVD is truncated, this is an operator-space approximation only; it does not preserve
        // positivity of the represented density matrix.

        IndexType qubit1 = controllingQubit1;
        IndexType qubit2 = qubit;
        bool reversed = false;

        if (qubit1 > qubit2)
        {
            std::swap(qubit1, qubit2);
            reversed = true;
        }

        // (4 * leftBond) x (4 * rightBond) matrix, the physical dimension per site is 4 = 2 (ket) x 2 (bra)
        MatrixClass thetaMatrix;
        if (IsExactSwapMatrix(gate))
        {
            if (TrySwapProductSites(qubit1, canonicalWeights))
                return;
            // SWAP exchanges both ket and bra legs, independently of argument order.
            // Keep the usual SVD/compression path, but skip the gate contractions.
            thetaMatrix = ReshapeSwappedTheta(ContractTwoQubits(qubit1));
        }
        else
        {
            const TwoQubitsGateTensor U = GetTwoQubitsGateTensor(gate, reversed);
            thetaMatrix = limitSize || limitEntanglement ? ReshapeThetaBar(ConstructThetaBarConservative(qubit1, U))
                                                         : ConstructThetaMatrix(qubit1, U);
        }

        ApplyPostTruncationPatches(DecomposeAndSetGammas(thetaMatrix, qubit1, qubit2, true, canonicalWeights));
    }

    // Called with the same prepared environments as DecomposeAndSetGammas.
    // For L = middle = R = 1, the swapped matrix is b * a^T. Its right factor
    // is a / ||a|| and its left factor is b * ||a||: preserve the operator's
    // scale, including mixed/unnormalized states whose only singular value is not one.
    bool TrySwapProductSites(IndexType first, bool canonicalWeights)
    {
        auto &left = gammas[first];
        auto &right = gammas[first + 1];
        if (left.dimension(0) != 1 || left.dimension(3) != 1 || right.dimension(0) != 1 || right.dimension(3) != 1 ||
            lambdas[first].size() != 1)
            return false;

        const Eigen::Map<const Eigen::Vector4cd> a(left.data()), b(right.data());
        const double aNorm = a.norm(), bNorm = b.norm();
        const double leftWeight = canonicalWeights && first != 0 ? std::abs(lambdas[first - 1][0]) : 1.;
        const double singularValue = leftWeight * aNorm * bNorm;
        if (!(aNorm > 0.) || !(singularValue > 0.) || !std::isfinite(singularValue))
            return false;

        Eigen::Vector4cd newLeft = b * aNorm;
        const Eigen::Vector4cd newRight = a / aNorm;
        if (!newLeft.allFinite() || !newRight.allFinite())
            return false;
        int64_t removed = 0;
        RescaleIfOutOfRange(newLeft, removed);
        const int valuesExponent = OutOfRangeExponent(singularValue);
        const double newLambda = singularValue * std::ldexp(1., -valuesExponent);
        if (!newLeft.allFinite() || !std::isfinite(newLambda))
            return false;

        InvalidateSamplingCache();
        std::swap(left, right);
        if (aNorm != 1. || removed != 0)
            Eigen::Map<Eigen::Vector4cd>(left.data()) = newLeft;
        if (aNorm != 1.)
            Eigen::Map<Eigen::Vector4cd>(right.data()) = newRight;
        lambdas[first][0] = newLambda;
        scaleExponent += removed;
        // No dimension changes or user truncation: no summary update or patches.
        if (!canonicalWeights)
        {
            canonicalFormValid = false;
            centerFirst = centerLast = first;
        }
        return true;
    }

    // Applies the configured limits to every bond that needs it. Returns whether anything was cut.
    bool Compress()
    {
        if (!limitSize && !limitEntanglement)
            return false;
        if (!limitEntanglement && std::all_of(lambdas.begin(), lambdas.end(),
                                              [this](const LambdaType &lambda) { return lambda.size() <= chi; }))
            return false;

        bool truncated = false;
        for (IndexType qubit1 = 0; qubit1 < static_cast<IndexType>(lambdas.size()); ++qubit1)
        {
            // After a cut, update each intervening bond's environment even when it
            // already fits the cap; a later cut must not reuse its old Schmidt weights.
            if (!limitEntanglement && lambdas[qubit1].size() <= chi && !truncated)
                continue;

            const bool canonicalWeights = canonicalFormValid;
            if (!canonicalWeights)
                PrepareBond(qubit1);
            // Contract with orthonormal environments, then apply the requested cut.
            const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
            const MatrixClass thetaMatrix = ReshapeTheta(theta);

            truncated = DecomposeAndSetGammas(thetaMatrix, qubit1, qubit1 + 1, true, canonicalWeights) || truncated;
        }
        return truncated;
    }

    void ApplyPostTruncationPatches(bool truncated)
    {
        if (!truncated || applyingPostTruncationPatches)
            return;
        if (deferPostTruncationPatches)
        {
            // a normalized operator applies them itself, once it has judged the weight (ApplyAndNormalizeWithoutBackup)
            truncatedWhileDeferred = true;
            return;
        }
        if (!hermitizeAfterTruncation && !restoreTraceAfterTruncation)
            return;

        ScopedFlag patchGuard(applyingPostTruncationPatches);
        if (hermitizeAfterTruncation)
            Hermitize();
        if (restoreTraceAfterTruncation)
            RestoreTraceIfSafe();
    }

    // SVD the (already built) theta matrix and write back the two new site tensors and the
    // lambda in between. Shared by two-qubit gates, Trim, and ReCanonicalize. User-requested
    // chi / entanglement cuts are applied only when applyUserCompression is true (the
    // default); ReCanonicalize passes false so it is a gauge restore.
    // Returns true when a user-requested cut actually dropped singular values.
    //
    // This is Hastings' method: theta = B1 B2 (with the gate applied) does not contain the left lambda.
    // The SVD is done on lambda_left * theta = U S V^dagger (the two site tensor as it would be in the Vidal form),
    // then the new right B is V^dagger and the new left B is theta * V (which is lambda_left^-1 * U * S, but computed
    // without dividing by the lambdas). U is not needed at all.
    // Nothing is rescaled, so B1 B2 = theta exactly if nothing is cut, and Tr(rho) is preserved.
    bool DecomposeAndSetGammas(const MatrixClass &thetaMatrix, IndexType qubit1, IndexType qubit2,
                               bool applyUserCompression = true, bool canonicalWeights = true)
    {
        // the SVD and the matrix product for the new left B are the parts Eigen parallelizes
        bool truncated = false;
        RunMaybeSingleThreaded(enableMultithreading, [&]() {
            truncated = DecomposeAndSetGammasImpl(thetaMatrix, qubit1, qubit2, applyUserCompression, canonicalWeights);
        });

        return truncated;
    }

    bool DecomposeAndSetGammasImpl(const MatrixClass &thetaMatrix, IndexType qubit1, IndexType qubit2,
                                   bool applyUserCompression, bool canonicalWeights)
    {
        InvalidateSamplingCache();
        const IndexType L = qubit1 == 0 ? 1 : lambdas[qubit1 - 1].size();
        const IndexType R = qubit2 == static_cast<IndexType>(lambdas.size()) ? 1 : lambdas[qubit2].size();

        assert(4 * L == thetaMatrix.rows());
        assert(4 * R == thetaMatrix.cols());

        // the rows are (physical index, left bond index) = p * L + l
        MatrixClass weightedTheta;
        if (canonicalWeights && qubit1 != 0)
        {
            const LambdaType &leftLambda = lambdas[qubit1 - 1];
            weightedTheta.resize(thetaMatrix.rows(), thetaMatrix.cols());
            for (IndexType c = 0; c < thetaMatrix.cols(); ++c)
                for (IndexType p = 0; p < 4; ++p)
                {
                    const IndexType pL = p * L;
                    for (IndexType l = 0; l < L; ++l)
                        weightedTheta(pL + l, c) = thetaMatrix(pL + l, c) * leftLambda[l];
                }
        }
        const MatrixClass &svdMatrix = !canonicalWeights || qubit1 == 0 ? thetaMatrix : weightedTheta;

#ifdef USE_FAST_SVD
        const bool computeWithJacobi = svdMatrix.rows() < blockSizeLimit && svdMatrix.cols() < blockSizeLimit;
#endif

        // only V is computed, U is not needed
#ifdef USE_FAST_SVD
        if (computeWithJacobi)
#endif
            jacobiSVD.compute(svdMatrix);
#ifdef USE_FAST_SVD
        else
            SVD.compute(svdMatrix);

        // A failed decomposition (non-finite input) leaves the previous results in the reused SVD
        // objects; never read them.
        if ((computeWithJacobi ? jacobiSVD.info() : SVD.info()) != Eigen::Success)
            throw std::runtime_error("MPO two-site SVD failed (non-finite or numerically invalid input)");
        const MatrixClass &VmatrixFull = computeWithJacobi ? jacobiSVD.matrixV() : SVD.matrixV();
        const LambdaType &SvaluesFull = computeWithJacobi ? jacobiSVD.singularValues() : SVD.singularValues();
        const double rankThreshold = computeWithJacobi ? jacobiSVD.threshold() : SVD.threshold();
        const IndexType nonzeroValues =
            computeWithJacobi ? jacobiSVD.nonzeroSingularValues() : SVD.nonzeroSingularValues();
#else
        if (jacobiSVD.info() != Eigen::Success)
            throw std::runtime_error("MPO two-site SVD failed (non-finite or numerically invalid input)");
        const MatrixClass &VmatrixFull = jacobiSVD.matrixV();
        const LambdaType &SvaluesFull = jacobiSVD.singularValues();
        const double rankThreshold = jacobiSVD.threshold();
        const IndexType nonzeroValues = jacobiSVD.nonzeroSingularValues();
#endif
        // The singular values carry the norm of the operator (2^-n/2 for a maximally mixed
        // chain): squared weights and Eigen's DBL_MIN rank floor would fail far from 1. Rank
        // decisions use them scaled by a power of two, which changes nothing in ordinary ranges.
        const int valuesExponent = OutOfRangeExponent(SvaluesFull.size() > 0 ? SvaluesFull[0] : 0.);
        const LambdaType Svalues =
            valuesExponent != 0 ? LambdaType(SvaluesFull * std::ldexp(1., -valuesExponent)) : SvaluesFull;

        // Eigen's own numerical rank (no threshold is set on the SVD objects): only the singular values
        // below diagSize * epsilon * sigma_max are dropped, they are not distinguishable from zero in double precision.
        // Any real cut is only the user's choice (setLimitEntanglement / setLimitBondDimension).
        const IndexType numericalRank = NumericalRank(Svalues, nonzeroValues, rankThreshold);

        // If user-requested compression is enabled, further reduce the rank according to
        // the configured truncation mode, applied on the (still descending-sorted) singular values.
        IndexType szm = (applyUserCompression && limitEntanglement)
                            ? ComputeCompressedRank(Svalues, numericalRank, truncationMode, singularValueThreshold)
                            : numericalRank;

        if (szm == 0)
            szm = 1;

        const IndexType sz = (applyUserCompression && limitSize) ? std::min<IndexType>(chi, szm) : szm;
        const bool truncated = applyUserCompression && sz < numericalRank;

        assert(sz <= VmatrixFull.cols());
        assert(4 * R == VmatrixFull.rows());

        // the new lambda is the raw singular values: NOT normalized, like the B tensors they are not rescaled
        // (Tr(rho) is preserved when untruncated). Out of the window they are scaled by a power of two:
        // lambdas only weight later SVD inputs, so their common scale does not matter.
        // Everything is built first and committed with noexcept moves.
        LambdaType newLambda = Svalues.head(sz);
        TensorType newLeft, newRight;
        const int64_t removed =
            SetNewGammas(thetaMatrix, VmatrixFull, L, sz, R, newLeft, newRight, limitSize || limitEntanglement);
        lambdas[qubit1] = std::move(newLambda);
        bondDimensions.Update(qubit1, lambdas[qubit1].size());
        gammas[qubit1] = std::move(newLeft);
        gammas[qubit2] = std::move(newRight);
        scaleExponent += removed;
        if (!canonicalWeights)
        {
            canonicalFormValid = false;
            centerFirst = centerLast = qubit1;
        }
        else if (truncated)
            InvalidateCanonicalForm(qubit1, qubit1);
        return truncated;
    }

    // See MPSSimulatorImpl::ComputeCompressedRank for the full explanation of both modes -
    // identical logic, duplicated here per this codebase's per-file convention (no shared
    // header between the MPS and MPO impls today). Always keeps at least the largest
    // singular value.
    // Eigen's SVDBase::rank() on the given (possibly rescaled) singular values.
    static IndexType NumericalRank(const LambdaType &sortedDescendingSVs, IndexType nonzero, double threshold)
    {
        if (sortedDescendingSVs.size() == 0)
            return 0;
        const double premultipliedThreshold =
            std::max(sortedDescendingSVs[0] * threshold, std::numeric_limits<double>::min());
        IndexType i = std::min<IndexType>(nonzero, sortedDescendingSVs.size()) - 1;
        while (i >= 0 && sortedDescendingSVs[i] < premultipliedThreshold)
            --i;
        return i + 1;
    }

    static IndexType ComputeCompressedRank(const LambdaType &sortedDescendingSVs, IndexType rank, TruncationMode mode,
                                           double threshold)
    {
        if (rank <= 1)
            return rank;

        if (mode == TruncationMode::RelativeToMax)
        {
            const double premultipliedThreshold =
                std::max(threshold * sortedDescendingSVs[0], std::numeric_limits<double>::min());
            IndexType i = rank - 1;
            while (i >= 0 && sortedDescendingSVs[i] < premultipliedThreshold)
                --i;
            return i + 1;
        }

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

    static bool IsExactSwapMatrix(const MatrixClass &gate)
    {
        // Check coefficients, not gate metadata: operators may be scaled, phased
        // or nearly SWAP. Only the exact permutation permits a plain leg exchange.
        constexpr IndexType swappedIndex[4] = {0, 2, 1, 3};
        for (IndexType column = 0; column < 4; ++column)
            for (IndexType row = 0; row < 4; ++row)
                if (gate(row, column) != (row == swappedIndex[column] ? 1. : 0.))
                    return false;
        return true;
    }

    static TwoQubitsGateTensor GetTwoQubitsGateTensor(const MatrixClass &gateMat, bool reversed)
    {
        TwoQubitsGateTensor result;

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

    // contracts the two adjacent sites into a rank-6 tensor
    // with leg order (leftBond, ket1, bra1, ket2, bra2, rightBond)
    // the B tensors already contain the lambdas: the middle one is in the left B, the right one in the right B
    // the left lambda is not included, it's needed only for the SVD (see DecomposeAndSetGammas)
    Eigen::Tensor<std::complex<double>, 6> ContractTwoQubits(IndexType qubit1) const
    {
        const IndexType qubit2 = qubit1 + 1;

        static const Indexes contractMid{IntIndexPair(3, 0)};

        assert(gammas[qubit1].dimension(3) == gammas[qubit2].dimension(0));

        return gammas[qubit1].contract(gammas[qubit2], contractMid);
    }

    // Preserve the original summation order when user compression is enabled.
    Eigen::Tensor<std::complex<double>, 6> ConstructThetaBarConservative(IndexType qubit1, const TwoQubitsGateTensor &U)
    {
        const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);

        const TwoQubitsGateTensor Uconj = U.conjugate();

        // contract the kets: theta(ket1=1, ket2=3) with U(in1=2, in2=3)
        // (leftBond, ket1, bra1, ket2, bra2, rightBond) -> (leftBond, bra1, bra2, rightBond, ket1', ket2')
        static const Indexes2 ketDims{IntIndexPair(1, 2), IntIndexPair(3, 3)};
        const Eigen::Tensor<std::complex<double>, 6> theta2 = theta.contract(U, ketDims);

        // contract the bras: theta2(bra1=1, bra2=2) with U*(in1=2, in2=3)
        // (leftBond, bra1, bra2, rightBond, ket1', ket2') -> (leftBond, rightBond, ket1', ket2', bra1', bra2')
        static const Indexes2 braDims{IntIndexPair(1, 2), IntIndexPair(2, 3)};
        return theta2.contract(Uconj, braDims);
    }

    // Apply the gate to each 4x4 physical block, writing the SVD layout
    // directly. This avoids two rank-six contraction results and a reshape.
    MatrixClass ConstructThetaMatrix(IndexType qubit1, const TwoQubitsGateTensor &U)
    {
        const Eigen::Tensor<std::complex<double>, 6> theta = ContractTwoQubits(qubit1);
        const IndexType L = theta.dimension(0), R = theta.dimension(5);
        MatrixClass result(4 * L, 4 * R);
        const TwoQubitsGateTensor Uconj = U.conjugate();
        for (IndexType r = 0; r < R; ++r)
            for (IndexType l = 0; l < L; ++l)
            {
                std::complex<double> ketApplied[4][4];
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                    {
                        ketApplied[a][b] = 0.;
                        for (int c = 0; c < 4; ++c)
                            ketApplied[a][b] += U(a / 2, a % 2, c / 2, c % 2) * theta(l, c / 2, b / 2, c % 2, b % 2, r);
                    }
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                    {
                        std::complex<double> value = 0.;
                        for (int c = 0; c < 4; ++c)
                            value += ketApplied[a][c] * Uconj(b / 2, b % 2, c / 2, c % 2);
                        result((2 * (a / 2) + b / 2) * L + l, (2 * (a % 2) + b % 2) * R + r) = value;
                    }
            }
        return result;
    }

    // reshapes the rank-6 theta into a (4 * leftBond) x (4 * rightBond) matrix for the SVD
    // the combined physical index per site is p = ket * 2 + bra
    static MatrixClass ReshapeThetaBar(const Eigen::Tensor<std::complex<double>, 6> &theta)
    {
        // theta dims: 0=leftBond, 1=rightBond, 2=ket1, 3=ket2, 4=bra1, 5=bra2
        const IndexType L = theta.dimension(0);
        const IndexType R = theta.dimension(1);

        assert(L > 0);
        assert(R > 0);

        MatrixClass thetaMatrix(4 * L, 4 * R);

        for (IndexType ket1 = 0; ket1 < 2; ++ket1)
            for (IndexType bra1 = 0; bra1 < 2; ++bra1)
            {
                const IndexType p1L = (ket1 * 2 + bra1) * L;
                for (IndexType ket2 = 0; ket2 < 2; ++ket2)
                    for (IndexType bra2 = 0; bra2 < 2; ++bra2)
                    {
                        const IndexType p2R = (ket2 * 2 + bra2) * R;
                        for (IndexType r = 0; r < R; ++r)
                            for (IndexType l = 0; l < L; ++l)
                                thetaMatrix(p1L + l, p2R + r) = theta(l, r, ket1, ket2, bra1, bra2);
                    }
            }

        return thetaMatrix;
    }

    // reshapes the raw two-site contraction (leg order leftBond, ket1, bra1, ket2, bra2, rightBond,
    // as produced by ContractTwoQubits) into a (4 * leftBond) x (4 * rightBond) matrix for the SVD.
    // Same target layout as ReshapeThetaBar (physical index per site p = ket * 2 + bra), but for the
    // no-gate leg order used by Trim.
    static MatrixClass ReshapeTheta(const Eigen::Tensor<std::complex<double>, 6> &theta)
    {
        // theta dims: 0=leftBond, 1=ket1, 2=bra1, 3=ket2, 4=bra2, 5=rightBond
        const IndexType L = theta.dimension(0);
        const IndexType R = theta.dimension(5);

        assert(L > 0);
        assert(R > 0);

        MatrixClass thetaMatrix(4 * L, 4 * R);

        for (IndexType ket1 = 0; ket1 < 2; ++ket1)
            for (IndexType bra1 = 0; bra1 < 2; ++bra1)
            {
                const IndexType p1L = (ket1 * 2 + bra1) * L;
                for (IndexType ket2 = 0; ket2 < 2; ++ket2)
                    for (IndexType bra2 = 0; bra2 < 2; ++bra2)
                    {
                        const IndexType p2R = (ket2 * 2 + bra2) * R;
                        for (IndexType r = 0; r < R; ++r)
                            for (IndexType l = 0; l < L; ++l)
                                thetaMatrix(p1L + l, p2R + r) = theta(l, ket1, bra1, ket2, bra2, r);
                    }
            }

        return thetaMatrix;
    }

    // Fuse a physical SWAP into the SVD reshape: the old second site's ket/bra
    // pair becomes the row index, and the old first site's pair becomes the column.
    // This needs no gate multiplications or intermediate shuffled rank-6 tensor.
    static MatrixClass ReshapeSwappedTheta(const Eigen::Tensor<std::complex<double>, 6> &theta)
    {
        const IndexType L = theta.dimension(0);
        const IndexType R = theta.dimension(5);

        assert(L > 0);
        assert(R > 0);

        MatrixClass thetaMatrix(4 * L, 4 * R);
        for (IndexType ket1 = 0; ket1 < 2; ++ket1)
            for (IndexType bra1 = 0; bra1 < 2; ++bra1)
            {
                const IndexType p1R = (ket1 * 2 + bra1) * R;
                for (IndexType ket2 = 0; ket2 < 2; ++ket2)
                    for (IndexType bra2 = 0; bra2 < 2; ++bra2)
                    {
                        const IndexType p2L = (ket2 * 2 + bra2) * L;
                        for (IndexType r = 0; r < R; ++r)
                            for (IndexType l = 0; l < L; ++l)
                                thetaMatrix(p2L + l, p1R + r) = theta(l, ket1, bra1, ket2, bra2, r);
                    }
            }
        return thetaMatrix;
    }

    // the left B is theta * V, the right B is V^dagger
    // V is the full V from the SVD (4R x rank), only the first sz columns are used
    // Returns the power of two removed from the left B when it is out of range (see scaleExponent).
    static int64_t SetNewGammas(const MatrixClass &thetaMatrix, const MatrixClass &Vmatrix, IndexType L, IndexType sz,
                                IndexType R, TensorType &Btensor1, TensorType &Btensor2, bool conservative)
    {
        // (4L x sz), the rows are (physical index, left bond index) = p * L + l, with p = ket * 2 + bra
        int64_t removed = 0;
        if (conservative)
        {
            MatrixClass leftB = thetaMatrix * Vmatrix.leftCols(sz);
            RescaleIfOutOfRange(leftB, removed);

            Btensor1.resize(L, 2, 2, sz);
            for (IndexType m = 0; m < sz; ++m)
                for (IndexType bra = 0; bra < 2; ++bra)
                    for (IndexType ket = 0; ket < 2; ++ket)
                    {
                        const IndexType pL = (ket * 2 + bra) * L;
                        for (IndexType l = 0; l < L; ++l)
                            Btensor1(l, ket, bra, m) = leftB(pL + l, m);
                    }
        }
        else
        {
            Btensor1.resize(L, 2, 2, sz);
            Eigen::Map<MatrixClass> leftB(Btensor1.data(), 4 * L, sz);
            leftB.noalias() = thetaMatrix * Vmatrix.leftCols(sz);
            RescaleIfOutOfRange(leftB, removed);

            // SVD rows use 2*ket+bra; tensor memory uses ket+2*bra. Only the
            // two off-diagonal physical blocks need exchanging after the product.
            for (IndexType m = 0; m < sz; ++m)
                for (IndexType l = 0; l < L; ++l)
                    std::swap(leftB(L + l, m), leftB(2 * L + l, m));
        }

        // the columns of V^dagger are (physical index, right bond index) = p * R + r
        Btensor2.resize(sz, 2, 2, R);
        for (IndexType r = 0; r < R; ++r)
            for (IndexType bra = 0; bra < 2; ++bra)
                for (IndexType ket = 0; ket < 2; ++ket)
                {
                    const IndexType idx = (ket * 2 + bra) * R + r;
                    for (IndexType m = 0; m < sz; ++m)
                        Btensor2(m, ket, bra, r) = std::conj(Vmatrix(idx, m));
                }
        return removed;
    }

    // with Hastings' method the U matrix from the SVD is not needed, only V
#ifdef USE_FAST_SVD
    constexpr static IndexType blockSizeLimit = 16;
    Eigen::BDCSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinV> SVD;
#endif
    constexpr static double krausCompletenessTolerance = 1E-10;
    Eigen::JacobiSVD<MatrixClass, Eigen::DecompositionOptions::ComputeThinV> jacobiSVD;
    bool applyingPostTruncationPatches = false;
    // Set while a normalized operator is applied: a truncation is then only recorded, see ApplyPostTruncationPatches.
    bool deferPostTruncationPatches = false;
    bool truncatedWhileDeferred = false;
};

} // namespace TensorNetworks
} // namespace QC
