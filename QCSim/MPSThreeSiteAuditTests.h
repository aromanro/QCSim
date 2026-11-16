#pragma once

#include "MPSSimulator.h"
#include "QubitRegister.h"
#include "TensorNetworkRegressionTests.h"

namespace QC
{
namespace MPSThreeSiteAudit
{
using namespace TensorNetworkRegression;
using Vector = Eigen::VectorXcd;
using Simulator = TensorNetworks::MPSSimulatorImpl;

inline Vector ApplyDense(const Vector &state, const Matrix &gate, const std::vector<int> &sites)
{
    Vector result = Vector::Zero(state.size());
    for (Eigen::Index input = 0; input < state.size(); ++input)
    {
        Eigen::Index outside = input;
        int local = 0;
        for (size_t bit = 0; bit < sites.size(); ++bit)
        {
            local |= ((input >> sites[bit]) & 1) << bit;
            outside &= ~(Eigen::Index{1} << sites[bit]);
        }
        for (int output = 0; output < gate.rows(); ++output)
        {
            Eigen::Index global = outside;
            for (size_t bit = 0; bit < sites.size(); ++bit)
                global |= Eigen::Index((output >> bit) & 1) << sites[bit];
            result(global) += gate(output, local) * state(input);
        }
    }
    return result;
}

// Independent full-state Schmidt truncation. This uses U*S*V^dagger over all
// physical qubits, not the local Hastings factors or the simulator's rank helper.
inline Vector Truncate(const Vector &state, int bond, int cap, double threshold, bool discardedWeight)
{
    const Eigen::Index rows = Eigen::Index{1} << (bond + 1);
    const Matrix matrix = Eigen::Map<const Matrix>(state.data(), rows, state.size() / rows);
    Eigen::JacobiSVD<Matrix> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto values = svd.singularValues();
    int rank = static_cast<int>(values.size());
    while (rank > 1 && values[rank - 1] < 1E-13 * values[0])
        --rank;
    int keep = rank;
    if (discardedWeight)
    {
        const double budget = threshold * values.head(rank).squaredNorm();
        double removed = 0.;
        while (keep > 1 && removed + values[keep - 1] * values[keep - 1] <= budget)
        {
            removed += values[keep - 1] * values[keep - 1];
            --keep;
        }
    }
    else
        while (keep > 1 && values[keep - 1] < threshold * values[0])
            --keep;
    keep = std::min(keep, cap);
    const Matrix output =
        svd.matrixU().leftCols(keep) * values.head(keep).asDiagonal() * svd.matrixV().leftCols(keep).adjoint();
    Vector result = Eigen::Map<const Vector>(output.data(), output.size());
    result.normalize();
    return result;
}

inline void Prepare(Simulator &simulator)
{
    simulator.SetMultithreading(false);
    for (int layer = 0; layer < 3; ++layer)
        for (int site = 0; site + 1 < static_cast<int>(simulator.getNrQubits()); ++site)
            simulator.ApplyGate(Gate(Unitary(4, 73 + layer * 11 + site)), site, site + 1);
}

inline void CheckTruncation()
{
    const Matrix gate = Unitary(8, 197);
    for (bool discarded : {false, true})
        for (int cap : {1, 2, 3, 8})
            for (double threshold : {0.001, 0.18})
            {
                Simulator simulator(6);
                Prepare(simulator);
                // The reference below is the optimal cut of the state as it is, which needs
                // current environments once an earlier cut has left the lambdas approximate.
                simulator.setThreeQubitEnvironmentRefresh(true);
                Vector reference = simulator.getRegisterStorage();
                simulator.setLimitBondDimension(cap);
                simulator.setLimitEntanglement(threshold);
                simulator.setTruncationMode(discarded ? Simulator::TruncationMode::DiscardedWeight
                                                      : Simulator::TruncationMode::RelativeToMax);
                // Move the interval, revisit it, and cross both chain boundaries.
                // This checks exterior-cache reuse as well as both SVD cuts.
                for (int first : {1, 1, 2, 0, 3, 2})
                {
                    reference = ApplyDense(reference, gate, {first + 2, first, first + 1});
                    reference = Truncate(reference, first + 1, cap, threshold, discarded);
                    reference = Truncate(reference, first, cap, threshold, discarded);
                    simulator.ApplyGate(Gate(gate), first + 2, first, first + 1);
                    Close(simulator.getRegisterStorage(), reference,
                          "MPS double-cut result differs from full-state Schmidt truncation");
                    const auto dimensions = simulator.getBondDimensions();
                    Require(dimensions[first] <= cap && dimensions[first + 1] <= cap,
                            "MPS double cut exceeded bond cap");
                    for (int site : {0, 2, 5})
                    {
                        double probability = 0.;
                        for (Eigen::Index basis = 0; basis < reference.size(); ++basis)
                            if (!(basis & (Eigen::Index{1} << site)))
                                probability += std::norm(reference[basis]);
                        Close(simulator.GetProbability(site), probability,
                              "MPS double cut left stale probability environments");
                    }
                }
                simulator.dontLimitBondDimension();
                simulator.dontLimitEntanglement();
                const auto snapshot = simulator.getState();
                simulator.ApplyGate(Gates::CNOTGate<>(), 3, 4);
                reference = ApplyDense(reference, Gates::CNOTGate<>().getRawOperatorMatrix(), {3, 4});
                Close(simulator.getRegisterStorage(), reference, "MPS two-site gate after audited truncation mismatch");
                simulator.setState(snapshot);
                reference = simulator.getRegisterStorage();
                simulator.ApplyGate(Gate(gate), 1, 2, 3);
                reference = ApplyDense(reference, gate, {1, 2, 3});
                Close(simulator.getRegisterStorage(), reference, "MPS snapshot retained stale three-site environments");
                const bool outcome = simulator.MeasureQubit(2);
                for (Eigen::Index basis = 0; basis < reference.size(); ++basis)
                    if (bool(basis & 4) != outcome)
                        reference[basis] = 0.;
                reference.normalize();
                Close(simulator.getRegisterStorage(), reference,
                      "MPS measurement after double cut used stale environments");
                simulator.ApplyGate(Gate(gate), 0, 1, 2);
                reference = ApplyDense(reference, gate, {0, 1, 2});
                Close(simulator.getRegisterStorage(), reference, "MPS three-site update after measurement mismatch");
            }
}

// Without the environment refresh a unitary three-qubit gate treats a truncated chain as a
// two-qubit gate does: it weighs its cuts by the left lambda and stays on the canonical path.
inline void CheckCanonicalWeights()
{
    const Matrix gate = Unitary(8, 197);
    // No cut anywhere: both settings give the exact evolution.
    for (bool refresh : {false, true})
    {
        Simulator simulator(6);
        Prepare(simulator);
        simulator.setThreeQubitEnvironmentRefresh(refresh);
        Require(simulator.getThreeQubitEnvironmentRefresh() == refresh, "MPS three-site setting was not stored");
        Vector reference = simulator.getRegisterStorage();
        for (int first : {1, 3, 0, 2})
        {
            reference = ApplyDense(reference, gate, {first + 2, first, first + 1});
            simulator.ApplyGate(Gate(gate), first + 2, first, first + 1);
        }
        Close(simulator.getRegisterStorage(), reference, "MPS untruncated three-site evolution depends on the setting");
    }

    Vector exact;
    {
        Simulator simulator(8);
        Prepare(simulator);
        exact = simulator.getRegisterStorage();
        for (int first : {2, 4, 1, 5, 3, 0})
            exact = ApplyDense(exact, gate, {first + 2, first, first + 1});
    }
    double fidelity[2];
    for (bool refresh : {false, true})
    {
        Simulator simulator(8);
        Prepare(simulator);
        simulator.setThreeQubitEnvironmentRefresh(refresh);
        simulator.setLimitBondDimension(4);
        for (int first : {2, 4, 1, 5, 3, 0})
        {
            simulator.ApplyGate(Gate(gate), first + 2, first, first + 1);
            const auto dimensions = simulator.getBondDimensions();
            Require(dimensions[first] <= 4 && dimensions[first + 1] <= 4, "MPS three-site cut exceeded the bond cap");
        }
        // A read after the cuts sees a normalized state whose probabilities are its own.
        Vector state = simulator.getRegisterStorage();
        Require(std::abs(state.norm() - 1.) < 1E-9, "MPS three-site cuts left the state off its norm");
        for (int site : {0, 3, 7})
        {
            double probability = 0.;
            for (Eigen::Index basis = 0; basis < state.size(); ++basis)
                if (!(basis & (Eigen::Index{1} << site)))
                    probability += std::norm(state[basis]);
            Close(simulator.GetProbability(site), probability, "MPS three-site cut left stale probability environments");
        }
        fidelity[refresh] = std::norm(exact.dot(state));
        // With the limits lifted, what follows is exact for the state the chain holds:
        // approximate weights only ever chose what a cut discards.
        simulator.dontLimitBondDimension();
        state = simulator.getRegisterStorage();
        simulator.ApplyGate(Gates::CNOTGate<>(), 3, 4);
        state = ApplyDense(state, Gates::CNOTGate<>().getRawOperatorMatrix(), {3, 4});
        Close(simulator.getRegisterStorage(), state, "MPS two-site gate after three-site cuts mismatch");
        simulator.ApplyGate(Gate(gate), 5, 3, 4);
        state = ApplyDense(state, gate, {5, 3, 4});
        Close(simulator.getRegisterStorage(), state, "MPS three-site gate after three-site cuts mismatch");
        const bool outcome = simulator.MeasureQubit(4);
        for (Eigen::Index basis = 0; basis < state.size(); ++basis)
            if (bool(basis & 16) != outcome)
                state[basis] = 0.;
        state.normalize();
        Close(simulator.getRegisterStorage(), state, "MPS measurement after three-site cuts mismatch");
    }
    // The lambda weights are those of a two-qubit update: a little below the optimal cut.
    Require(fidelity[0] > 0. && fidelity[0] > fidelity[1] - 0.05,
            "MPS three-site cuts by canonical weights lost far more than the refreshed ones");

    // The setting follows a copy of the simulator.
    TensorNetworks::MPSSimulator routed(5);
    routed.setThreeQubitEnvironmentRefresh(true);
    Require(routed.Clone()->getThreeQubitEnvironmentRefresh(), "MPS three-site setting was lost by Clone");
}

inline void CheckExtremeFilters()
{
    const Matrix unitary = Unitary(8, 311);
    for (double scale : {1E-310, 1E-200, 1E200, 1E300})
    {
        Simulator simulator(5);
        Prepare(simulator);
        const Matrix filter = scale * unitary;
        // Normalize the represented coefficients without forming 1/scale.
        Matrix normalized = filter;
        for (Eigen::Index i = 0; i < normalized.size(); ++i)
            normalized.data()[i] = {filter.data()[i].real() / scale, filter.data()[i].imag() / scale};
        Vector expected = ApplyDense(simulator.getRegisterStorage(), normalized, {1, 2, 3});
        expected.normalize();
        simulator.ApplyGate(Gate(filter), 1, 2, 3);
        Close(simulator.getRegisterStorage(), expected,
              "MPS finite extreme-scale filter lost normalization or precision");
        simulator.ApplyGate(Gate(unitary), 2, 3, 4);
        expected = ApplyDense(expected, unitary, {2, 3, 4});
        Close(simulator.getRegisterStorage(), expected, "MPS extreme filter left invalid subsequent environments");
        TensorNetworks::MPSSimulator routed(5);
        routed.SetMultithreading(false);
        routed.ApplyGate(Gates::HadamardGate<>(), 2);
        expected = ApplyDense(routed.getRegisterStorage(), normalized, {0, 2, 4});
        expected.normalize();
        routed.ApplyGate(Gate(filter), 0, 2, 4);
        Close(routed.getRegisterStorage(), expected, "MPS wrapper rejected or damaged a nonzero extreme filter");
    }
    for (double amplitude : {1E-150, 1E-200, 1E-300})
    {
        Simulator simulator(5);
        simulator.SetMultithreading(false);
        simulator.ApplyGate(Gates::RyGate<>(2 * amplitude), 0);
        // The rare amplitude sits inside the selected interval, so the exact
        // rank cutoff on an earlier entangling gate cannot discard it first.
        Matrix projection = Matrix::Zero(8, 8);
        for (int i = 1; i < 8; i += 2)
            projection(i, i) = 1.;
        simulator.ApplyGate(Gate(projection), 0, 1, 2);
        Vector expected = Vector::Zero(32);
        expected[1] = 1.;
        Close(simulator.getRegisterStorage(), expected, "MPS rare nonzero branch underflowed during normalization");
    }
    Simulator simulator(5);
    Prepare(simulator);
    const Vector before = simulator.getRegisterStorage();
    Throws<std::domain_error>([&] { simulator.ApplyGate(Gate(Matrix::Zero(8, 8)), 1, 2, 3); },
                              "MPS zero three-site operator was accepted");
    Close(simulator.getRegisterStorage(), before, "MPS rejected zero operator changed state");
    TensorNetworks::MPSSimulator routed(5);
    routed.ApplyGate(Gates::HadamardGate<>(), 2);
    const auto map = routed.getQubitsMap();
    const Vector routedBefore = routed.getRegisterStorage();
    Throws<std::domain_error>([&] { routed.ApplyGate(Gate(Matrix::Zero(8, 8)), 0, 2, 4); },
                              "MPS routed zero operator was accepted");
    Require(routed.getQubitsMap() == map, "Rejected MPS zero operator changed routing");
    Close(routed.getRegisterStorage(), routedBefore, "Rejected routed MPS zero operator changed state");
    Simulator annihilated(5);
    Matrix impossible = Matrix::Identity(8, 8);
    impossible(0, 0) = 0.;
    const Vector unselected = annihilated.getRegisterStorage();
    Throws<std::domain_error>([&] { annihilated.ApplyGate(Gate(impossible), 1, 2, 3); },
                              "MPS zero-probability branch was accepted");
    Close(annihilated.getRegisterStorage(), unselected, "Rejected MPS zero-probability branch changed state");
    annihilated.ApplyGate(Gate(unitary), 1, 2, 3);
    Close(annihilated.getRegisterStorage(), ApplyDense(unselected, unitary, {1, 2, 3}),
          "Rejected MPS filter left stale environments");
}

inline void CheckRoutingAfterReads()
{
    for (int routing = 0; routing < 3; ++routing)
    {
        TensorNetworks::MPSSimulator read(6);
        read.SetMultithreading(false);
        read.setLimitBondDimension(2);
        for (int layer = 0; layer < 3; ++layer)
            for (int site = 0; site < 5; ++site)
                read.ApplyGate(Gate(Unitary(4, 401 + 7 * layer + site)), site, site + 1);
        read.SetUseOptimalMeetingPosition(routing != 0);
        if (routing == 2)
            read.SetMeetingPositionCallback([](const auto &) { return 2; });
        auto unread = read.Clone();
        int readNotifications = 0, unreadNotifications = 0;
        read.SetBondDimensionCallback([&](const auto &) { ++readNotifications; });
        unread->SetBondDimensionCallback([&](const auto &) { ++unreadNotifications; });
        // Truncation leaves approximate canonical weights. A probability read
        // creates temporary environments that must be dropped by the first
        // internal routing SWAP, just as by the old public gate entry point.
        (void)read.GetProbability(3);
        const Gate gate(Unitary(4, 439));
        read.ApplyGate(gate, 5, 0);
        unread->ApplyGate(gate, 5, 0);
        Close(read.getRegisterStorage(), unread->getRegisterStorage(),
              "MPS routing after a read changed truncated gate evolution");
        Require(read.getQubitsMap() == unread->getQubitsMap(), "MPS routing after a read changed logical mapping");
        Require(readNotifications > 1 && readNotifications == unreadNotifications,
                "MPS routing after a read changed bond callback boundaries");
        (void)read.GetProbability(1);
        read.MoveAtBeginningOfChain({3, 5});
        unread->MoveAtBeginningOfChain({3, 5});
        Close(read.getRegisterStorage(), unread->getRegisterStorage(),
              "MPS chain reordering after a read changed truncated state");
        Require(read.getQubitsMap() == unread->getQubitsMap() && readNotifications == unreadNotifications,
                "MPS chain reordering changed maps or callback boundaries");
    }
}

inline void CheckRejectedFilterRollback()
{
    struct Inspectable : Simulator
    {
        using Simulator::Simulator;

        auto RawState() const
        {
            return MPSSimulatorBase::getState();
        }
    };

    Matrix projector = Matrix::Zero(8, 8);
    projector(7, 7) = 1.;
    for (bool cached : {false, true})
    {
        Inspectable adjacent(5);
        adjacent.SetMultithreading(false);
        adjacent.ApplyGate(Gates::HadamardGate<>(), 0);
        adjacent.ApplyGate(Gates::CNOTGate<>(), 1, 0);
        if (cached)
        {
            Matrix k0 = Matrix::Identity(2, 2), k1 = Matrix::Zero(2, 2);
            k0(1, 1) = 0.9;
            k1(0, 1) = std::sqrt(0.19);
            adjacent.ApplyKrausOperator({k0, k1}, 0.2, 0);
            (void)adjacent.GetProbability(0);
        }
        const auto before = std::static_pointer_cast<TensorNetworks::MPSSimulatorBaseState>(adjacent.RawState());
        const Vector expected = adjacent.getRegisterStorage();
        Throws<std::domain_error>([&] { adjacent.ApplyGate(Gate(projector, 1, 2, 3)); },
                                  "Annihilating adjacent MPS filter was accepted");
        const auto after = std::static_pointer_cast<TensorNetworks::MPSSimulatorBaseState>(adjacent.RawState());
        Require(before->canonicalFormExact == after->canonicalFormExact &&
                    before->schmidtValuesDiscarded == after->schmidtValuesDiscarded,
                "Rejected MPS filter changed canonical metadata");
        for (size_t site = 0; site < before->gammas.size(); ++site)
        {
            const auto &a = before->gammas[site], &b = after->gammas[site];
            Require(a.dimensions() == b.dimensions() && std::equal(a.data(), a.data() + a.size(), b.data()),
                    "Rejected MPS filter changed a site tensor");
        }
        for (size_t bond = 0; bond < before->lambdas.size(); ++bond)
            Require(before->lambdas[bond] == after->lambdas[bond], "Rejected MPS filter changed Schmidt values");
        adjacent.ApplyGate(Gates::CNOTGate<>(), 2, 1);
        Close(adjacent.getRegisterStorage(), ApplyDense(expected, Gates::CNOTGate<>().getRawOperatorMatrix(), {2, 1}),
              "Rejected MPS filter left stale environments");
    }

    TensorNetworks::MPSSimulator routed(5);
    routed.SetMultithreading(false);
    routed.ApplyGate(Gates::HadamardGate<>(), 0);
    routed.ApplyGate(Gates::CNOTGate<>(), 1, 0);
    routed.SaveState();
    const auto map = routed.getQubitsMap(), dimensions = routed.getBondDimensions();
    const Vector routedBefore = routed.getRegisterStorage();
    int notifications = 0;
    routed.SetBondDimensionCallback([&](const auto &) { ++notifications; });
    routed.SetBondDimensionSummaryCallback([&](auto) { ++notifications; });
    routed.setLimitBondDimension(1); // routing would discard the Bell pair
    Throws<std::domain_error>([&] { routed.ApplyGate(Gate(projector, 0, 2, 4)); },
                              "Annihilating routed MPS filter was accepted");
    Require(routed.getQubitsMap() == map && routed.getBondDimensions() == dimensions && notifications == 0,
            "Rejected MPS filter exposed tentative routing or bond changes");
    Close(routed.getRegisterStorage(), routedBefore, "Rejected MPS filter retained routing truncation");
    routed.dontLimitBondDimension();
    routed.ApplyGate(Gates::HadamardGate<>(), 4);
    routed.RestoreState();
    Close(routed.getRegisterStorage(), routedBefore, "Rejected MPS filter overwrote a saved state");
    auto reference = routed.Clone();
    std::vector<std::vector<Simulator::IndexType>> observedDimensions, referenceDimensions;
    std::vector<Simulator::IndexType> observedMaxima, referenceMaxima;
    routed.SetBondDimensionCallback([&](const auto &value) {
        ++notifications;
        observedDimensions.push_back(value);
    });
    routed.SetBondDimensionSummaryCallback([&](auto value) {
        ++notifications;
        observedMaxima.push_back(value);
    });
    reference->SetBondDimensionCallback([&](const auto &value) { referenceDimensions.push_back(value); });
    reference->SetBondDimensionSummaryCallback([&](auto value) { referenceMaxima.push_back(value); });
    reference->ApplyGate(Gate(Matrix::Identity(8, 8), 0, 2, 4));
    routed.ApplyGate(Gate(Matrix(Matrix::Identity(8, 8) * 0.5), 0, 2, 4));
    Close(routed.getRegisterStorage(), routedBefore, "Successful MPS filter changed normalized state");
    Require(notifications == 6, "Successful MPS filter changed per-swap bond notification counts");
    Require(observedDimensions == referenceDimensions && observedMaxima == referenceMaxima,
            "Successful MPS filter changed per-swap bond notification payloads");
}

inline bool Run()
{
    try
    {
        CheckTruncation();
        CheckCanonicalWeights();
        CheckExtremeFilters();
        CheckRoutingAfterReads();
        CheckRejectedFilterRollback();
        std::cout << "MPS independent three-site truncation and scale audit: passed\n";
        return true;
    }
    catch (const std::exception &error)
    {
        std::cerr << "MPS independent three-site audit: " << error.what() << '\n';
        return false;
    }
}
} // namespace MPSThreeSiteAudit
} // namespace QC
