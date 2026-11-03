#pragma once

#include "PauliPropagator.h"
#include <cmath>
#include <iostream>
#include <map>
#include <random>

namespace QC
{
namespace PauliLightconeRegression
{
inline void Require(bool ok, const char *message)
{
    if (!ok)
        throw std::runtime_error(message);
}

inline void Close(double a, double b, const char *message)
{
    Require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) < 3e-10, message);
}

inline std::map<std::string, double> Canonical(const PauliStringStorage &terms)
{
    std::map<std::string, double> result;
    for (const auto &term : terms)
        result[term.PauliStringXZ::ToString()] += term.Coefficient;
    return result;
}

// Independent public virtual-operator evolution, with the original checkpoint
// indices and no cone or finished-qubit reduction.
inline PauliStringStorage Reference(const PauliPropagator &sim, PauliStringStorage terms)
{
    const auto ops = sim.GetOperations();
    for (size_t next = ops.size(); next;)
    {
        const size_t index = --next;
        PauliStringStorage extra;
        for (auto &term : terms)
            ops[index]->Apply(term, extra);
        for (auto &term : extra)
            terms.push_back(std::move(term));
        const bool dedup = sim.StepsBetweenDeduplication() != std::numeric_limits<int>::max() &&
                           index % sim.StepsBetweenDeduplication() == 0;
        const bool trim =
            sim.StepsBetweenTrims() != std::numeric_limits<int>::max() && index % sim.StepsBetweenTrims() == 0;
        if (dedup)
        {
            std::map<std::string, PauliStringXZWithCoefficient> unique;
            for (const auto &term : terms)
            {
                if (term.PauliWeight() > sim.GetPauliWeightThreshold())
                    continue;
                const auto key = term.PauliStringXZ::ToString();
                auto found = unique.find(key);
                if (found == unique.end())
                    unique.emplace(key, term);
                else
                    found->second.Coefficient += term.Coefficient;
            }
            terms.clear();
            for (auto &entry : unique)
                terms.push_back(std::move(entry.second));
        }
        if (trim || dedup)
            terms.erase(std::remove_if(terms.begin(), terms.end(),
                                       [&](const auto &term) {
                                           return std::abs(term.Coefficient) <= sim.GetCoefficientThreshold() ||
                                                  term.PauliWeight() > sim.GetPauliWeightThreshold();
                                       }),
                        terms.end());
    }
    return terms;
}

inline double Sum(const PauliStringStorage &terms)
{
    double result = 0.;
    for (const auto &term : terms)
        if (std::none_of(term.X.begin(), term.X.end(), [](bool bit) { return bit; }))
            result += term.Coefficient;
    return result;
}

inline void Check(const PauliPropagator &sim, const PauliStringXZWithCoefficient &input)
{
    const auto expected = Reference(sim, {input});
    Close(sim.ExpectationValue(input), Sum(expected), "lightcone scalar/reference mismatch");
    PauliStringStorage output;
    Close(sim.ExpectationValue(input, output), Sum(expected), "lightcone output/reference mismatch");
    const auto want = Canonical(expected), got = Canonical(output);
    for (const auto &entry : want)
    {
        const auto found = got.find(entry.first);
        Close(entry.second, found == got.end() ? 0. : found->second, "lightcone changed propagated expansion");
    }
    for (const auto &entry : got)
    {
        const auto found = want.find(entry.first);
        Close(entry.second, found == want.end() ? 0. : found->second, "lightcone added propagated expansion");
    }
}

class RemoteX : public OperatorX
{
  public:
    explicit RemoteX(int remote) : OperatorX(0), remote(remote)
    {
    }

    void Apply(PauliStringXZWithCoefficient &term, PauliStringStorage &) const override
    {
        term.X[remote] = !term.X[remote];
    }

    std::unique_ptr<Operator> Clone() const override
    {
        return std::make_unique<RemoteX>(*this);
    }

  private:
    int remote;
};

inline void Run()
{
    std::mt19937 random(731);
    for (int width : {4, 65, 129, 257})
        for (int trial = 0; trial < 18; ++trial)
        {
            PauliPropagator sim;
            sim.SetNrQubits(width);
            sim.SetStepsBetweenDeduplication(3);
            sim.SetStepsBetweenTrims(2);
            sim.SetCoefficientThreshold(trial % 3 == 0 ? .017 : 0.);
            sim.SetPauliWeightThreshold(trial % 4 == 0 ? 2 : std::numeric_limits<size_t>::max());
            const int active[] = {0, width / 2, width - 1};
            for (int gate = 0; gate < 12; ++gate)
            {
                const int ia = random() % 3, ib = (ia + 1 + random() % 2) % 3;
                const int a = active[ia], b = active[ib], c = active[3 - ia - ib];
                const double angle = .1 + (random() % 100) * .013;
                switch (random() % 8)
                {
                case 0:
                    sim.ApplyH(a);
                    break;
                case 1:
                    sim.ApplyRY(a, angle);
                    break;
                case 2:
                    sim.ApplyRZ(a, angle);
                    break;
                case 3:
                    sim.ApplyCX(a, b);
                    break;
                case 4:
                    sim.ApplyU(a, angle, -.2, .4);
                    break;
                case 5:
                    sim.ApplyCRY(a, b, angle);
                    break;
                case 6:
                    sim.ApplyCP(a, b, angle);
                    break;
                default:
                    sim.ApplyCSwap(a, b, c);
                    break;
                }
            }
            if (trial % 2 == 0)
            {
                auto ops = sim.GetOperations();
                // Outside the observable's initial support: projectors must
                // still contribute their probability and grow the cone.
                ops.push_back(std::make_unique<Projector>(1, trial % 4 == 0, .5));
                sim.SetOperations(std::move(ops));
            }
            PauliStringXZWithCoefficient input(width);
            for (int q : active)
            {
                input.X[q] = random() & 1;
                input.Z[q] = random() & 1;
            }
            Check(sim, input);
            if (trial == 1)
            {
                sim.EnableParallel(2);
                sim.SetParallelThreshold(1);
                sim.SetBatchSize(1);
                Check(sim, input);
            }
        }

    PauliPropagator independent;
    independent.SetNrQubits(129);
    PauliStringXZWithCoefficient product(129);
    double expected = 1.;
    for (int q : {0, 1, 63, 64, 65, 127, 128})
    {
        independent.ApplyRY(q, .37);
        product.X[q] = true;
        expected *= std::sin(.37);
    }
    Close(independent.ExpectationValue(product), expected, "finished-qubit product mismatch");
    product.X[7] = true;
    Close(independent.ExpectationValue(product), 0., "untouched X was not rejected");
    product.X[7] = false;
    product.Z[7] = true;
    Check(independent, product);
    independent.SaveState();
    PauliPropagator clone;
    clone.SetNrQubits(129);
    clone.ShareOperationsFrom(independent);
    independent.ApplyH(7);
    Check(independent, product);
    Check(clone, product);
    independent.RestoreState();
    Check(independent, product);
    PauliStringStorage initial(64, product);
    for (size_t i = 0; i < initial.size(); ++i)
    {
        initial[i].Coefficient = 1. / (i + 1);
        initial[i].X[7] = i % 2;
    }
    Close(independent.ExpectationValue(initial), Sum(Reference(independent, initial)),
          "initial large frontier contraction mismatch");

    PauliPropagator custom;
    custom.SetNrQubits(129);
    std::vector<std::unique_ptr<Operator>> ops;
    ops.push_back(std::make_unique<OperatorH>(128));
    ops.push_back(std::make_unique<RemoteX>(128));
    ops.push_back(std::make_unique<OperatorRY>(1, .71));
    custom.SetOperations(std::move(ops));
    PauliStringXZWithCoefficient identity(129);
    Close(custom.ExpectationValue(identity), 1., "custom remote support was pruned");
    Check(custom, identity);
    identity.X[128] = true;
    Close(custom.ExpectationValue(identity), 1., "custom never-touched X was rejected early");
    Check(custom, identity);

    // The contraction threshold is crossed after several branches. A custom
    // operator anywhere in the circuit must still prevent early contraction.
    PauliPropagator delayedCustom;
    delayedCustom.SetNrQubits(129);
    std::vector<std::unique_ptr<Operator>> delayedOps;
    delayedOps.push_back(std::make_unique<RemoteX>(128));
    delayedCustom.SetOperations(std::move(delayedOps));
    PauliStringXZWithCoefficient delayedInput(129);
    delayedInput.X[128] = true;
    for (int q = 1; q <= 6; ++q)
    {
        delayedCustom.ApplyRY(q, .37);
        delayedInput.X[q] = true;
    }
    Close(delayedCustom.ExpectationValue(delayedInput), std::pow(std::sin(.37), 6),
          "delayed contraction ignored a custom operator");
    Check(delayedCustom, delayedInput);

    // Imported local tables are public and need not represent unitary gates.
    // Their identity column must be checked before treating them as skippable.
    for (bool producesX : {false, true})
    {
        auto table = std::make_shared<PauliDetail::LocalTransfer>(1);
        for (unsigned input = 0; input < 4; ++input)
        {
            std::array<double, 4> values{};
            values[input == 0 && producesX ? 1 : input] = input == 0 && !producesX ? 2. : 1.;
            table->Column(input, values);
        }
        PauliPropagator nonunital;
        nonunital.SetNrQubits(129);
        std::vector<std::unique_ptr<Operator>> imported;
        imported.push_back(std::make_unique<OperatorLocal>(OperationType::U, 128, 0, 0, table));
        nonunital.SetOperations(std::move(imported));
        PauliStringXZWithCoefficient unit(129);
        Close(nonunital.ExpectationValue(unit), producesX ? 0. : 2., "nonunital identity action was skipped");
        Check(nonunital, unit);
    }

    // Contraction merges labels between checkpoints. It must not trim a NaN,
    // or a NaN made from Inf - Inf, that the exact unmerged sum would keep.
    for (bool nanAngle : {false, true})
    {
        PauliPropagator nonfinite;
        nonfinite.SetNrQubits(8);
        for (int q = 1; q <= 6; ++q)
            nonfinite.ApplyZ(q);
        nonfinite.ApplyRY(0, nanAngle ? std::numeric_limits<double>::quiet_NaN() : .3);
        PauliStringStorage frontier;
        for (int k = 0; k < 64; ++k)
        {
            PauliStringXZWithCoefficient term(8);
            term.Z[0] = true;
            for (int b = 0; b < 6; ++b)
                term.Z[1 + b] = (k >> b) & 1;
            term.Coefficient = 1. / 64;
            frontier.push_back(term);
        }
        if (!nanAngle)
        {
            // Never-touched qubit 7 is contracted first, giving these one label.
            frontier[0].Coefficient = std::numeric_limits<double>::infinity();
            frontier[1] = frontier[0];
            frontier[1].Z[7] = true;
            frontier[1].Coefficient = -std::numeric_limits<double>::infinity();
        }
        Require(std::isnan(Sum(Reference(nonfinite, frontier))), "nonfinite regression lost its reference NaN");
        Require(std::isnan(nonfinite.ExpectationValue(frontier)), "contraction trimmed a nonfinite expectation");
    }
    std::cout << "PASS lightcone scalar/expansion references, sparse/native/projector/custom circuits, "
                 "approximation checkpoints, wide masks and shared-circuit lifecycle\n";
}
} // namespace PauliLightconeRegression
} // namespace QC
