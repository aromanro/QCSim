#pragma once

#include "QuantumGate.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace QC
{
// The same public compact-gate API for statevector and density-matrix simulators.
// As with ApplyGate, controlled matrices need not be unitary. Uncompute assumes
// unitarity and applies the adjoint. Angles use the exp(-i theta P / 2) convention.
template <class Simulator, class MatrixClass> class StructuredGateOperations
{
  public:
    static constexpr size_t AllControls = std::numeric_limits<size_t>::max();

    // The table is applied where it is: a simulator with a view entry point
    // copies it only if it records the gate.
    void ApplyDiagonalGate(const std::vector<std::complex<double>> &values, const std::vector<size_t> &qubits)
    {
        ApplyView(Gates::StructuredGateView::Diagonal(values, qubits));
    }

    void ApplyMultiControlledGate(const MatrixClass &matrix, size_t target, const std::vector<size_t> &controls,
                                  size_t controlState = AllControls)
    {
        if (matrix.rows() != 2 || matrix.cols() != 2)
            throw std::invalid_argument("A multi-controlled target matrix must be 2x2");
        ApplyView(Gates::StructuredGateView::Controlled({matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1)}, target,
                                                        controls, controlState));
    }

    void ApplyMultiControlledGate(const Gates::QuantumGateWithOp<MatrixClass> &gate, size_t target,
                                  const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        if (gate.getQubitsNumber() != 1)
            throw std::invalid_argument("A multi-controlled target must be a one-qubit gate");
        ApplyMultiControlledGate(gate.getRawOperatorMatrix(), target, controls, controlState);
    }

    void ApplyMCX(size_t target, const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        ApplyMultiControlledGate(Gates::PauliXGate<MatrixClass>(), target, controls, controlState);
    }

    void ApplyMCZ(size_t target, const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        ApplyMultiControlledGate(Gates::PauliZGate<MatrixClass>(), target, controls, controlState);
    }

    void ApplyMCP(double angle, size_t target, const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        if (!std::isfinite(angle))
            throw std::invalid_argument("Gate angle must be finite");
        double c, s;
        detail::QuarterTurnCosSin(angle, c, s);
        ApplyControlledValues({1., 0., 0., std::complex<double>(c, s)}, target, controls, controlState);
    }

    void ApplyMCRX(double angle, size_t target, const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        ApplyControlledRotation(Gates::RotationAxis::X, angle, target, controls, controlState);
    }

    void ApplyMCRY(double angle, size_t target, const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        ApplyControlledRotation(Gates::RotationAxis::Y, angle, target, controls, controlState);
    }

    void ApplyMCRZ(double angle, size_t target, const std::vector<size_t> &controls, size_t controlState = AllControls)
    {
        ApplyControlledRotation(Gates::RotationAxis::Z, angle, target, controls, controlState);
    }

    void ApplyMCU(double theta, double phi, double lambda, size_t target, const std::vector<size_t> &controls,
                  double gamma = 0., size_t controlState = AllControls)
    {
        ApplyMultiControlledGate(Gates::UGate<MatrixClass>(theta, phi, lambda, gamma), target, controls, controlState);
    }

    void ApplyUniformlyControlledRotation(Gates::RotationAxis axis, const std::vector<double> &angles, size_t target,
                                          const std::vector<size_t> &controls)
    {
        ApplyView(Gates::StructuredGateView::Multiplexed(axis, angles, target, controls));
    }

    void ApplyMultiplexedRX(const std::vector<double> &angles, size_t target, const std::vector<size_t> &controls)
    {
        ApplyUniformlyControlledRotation(Gates::RotationAxis::X, angles, target, controls);
    }

    void ApplyMultiplexedRY(const std::vector<double> &angles, size_t target, const std::vector<size_t> &controls)
    {
        ApplyUniformlyControlledRotation(Gates::RotationAxis::Y, angles, target, controls);
    }

    void ApplyMultiplexedRZ(const std::vector<double> &angles, size_t target, const std::vector<size_t> &controls)
    {
        ApplyUniformlyControlledRotation(Gates::RotationAxis::Z, angles, target, controls);
    }

  private:
    void ApplyControlledValues(const std::array<std::complex<double>, 4> &matrix, size_t target,
                               const std::vector<size_t> &controls, size_t controlState)
    {
        ApplyView(Gates::StructuredGateView::Controlled(matrix, target, controls, controlState));
    }

    void ApplyControlledRotation(Gates::RotationAxis axis, double angle, size_t target,
                                 const std::vector<size_t> &controls, size_t controlState)
    {
        if (!std::isfinite(angle))
            throw std::invalid_argument("Gate angle must be finite");
        double c, s;
        detail::RotationHalfAngleCosSin(angle, c, s);
        using Complex = std::complex<double>;
        if (axis == Gates::RotationAxis::X)
            ApplyControlledValues({c, Complex(0., -s), Complex(0., -s), c}, target, controls, controlState);
        else if (axis == Gates::RotationAxis::Y)
            ApplyControlledValues({c, -s, s, c}, target, controls, controlState);
        else
            ApplyControlledValues({Complex(c, -s), 0., 0., Complex(c, s)}, target, controls, controlState);
    }

    // A simulator that can apply borrowed tables takes the view; any other
    // one is handed an owning copy, as before.
    template <class S>
    static auto ApplyViewTo(S &simulator, const Gates::StructuredGateView &gate, int)
        -> decltype(simulator.ApplyStructuredGateView(gate))
    {
        return simulator.ApplyStructuredGateView(gate);
    }

    template <class S> static void ApplyViewTo(S &simulator, const Gates::StructuredGateView &gate, long)
    {
        simulator.ApplyStructuredGate(gate.Materialize());
    }

    void ApplyView(const Gates::StructuredGateView &gate)
    {
        ApplyViewTo(*static_cast<Simulator *>(this), gate, 0);
    }
};

} // namespace QC
