#pragma once

#include "GateAngles.h"
#include "SimpleGates.h"
#include <stdexcept>
#include <type_traits>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them
// usually from left to right

namespace QC
{
namespace Gates
{

template <class MatrixClass = Eigen::MatrixXcd> class SwapGate : public TwoQubitsGate<MatrixClass>
{
  public:
    using BaseClass = TwoQubitsGate<MatrixClass>;
    using OpClass = typename TwoQubitsGate<MatrixClass>::BaseClass;

    SwapGate()
    {
        OpClass::operatorMat(1, 1) = 0;
        OpClass::operatorMat(2, 2) = 0;
        OpClass::operatorMat(1, 2) = 1;
        OpClass::operatorMat(2, 1) = 1;
    }

    bool isSwapGate() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class iSwapGate : public TwoQubitsGate<MatrixClass>
{
  public:
    using BaseClass = TwoQubitsGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    iSwapGate()
    {
        OpClass::operatorMat(1, 1) = 0;
        OpClass::operatorMat(2, 2) = 0;
        OpClass::operatorMat(1, 2) = std::complex<double>(0, 1);
        OpClass::operatorMat(2, 1) = std::complex<double>(0, 1);
    }

    bool IsISwapGate() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class iSwapDagGate : public TwoQubitsGate<MatrixClass>
{
  public:
    using BaseClass = TwoQubitsGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    iSwapDagGate()
    {
        OpClass::operatorMat(1, 1) = 0;
        OpClass::operatorMat(2, 2) = 0;
        OpClass::operatorMat(1, 2) = std::complex<double>(0, -1);
        OpClass::operatorMat(2, 1) = std::complex<double>(0, -1);
    }

    bool IsISwapDagGate() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class DecrementGate : public TwoQubitsGate<MatrixClass>
{
  public:
    using BaseClass = TwoQubitsGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    DecrementGate()
    {
        OpClass::operatorMat(0, 0) = 0;
        OpClass::operatorMat(1, 1) = 0;
        OpClass::operatorMat(2, 2) = 0;
        OpClass::operatorMat(3, 3) = 0;

        OpClass::operatorMat(3, 0) = 1;

        OpClass::operatorMat(0, 1) = 1;
        OpClass::operatorMat(1, 2) = 1;
        OpClass::operatorMat(2, 3) = 1;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class TwoQubitsControlledGate : public TwoQubitsGate<MatrixClass>
{
  public:
    using BaseClass = TwoQubitsGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    TwoQubitsControlledGate() {};

    bool isControlled() const override
    {
        return true;
    }

    bool isControlQubit(size_t qubit) const override
    {
        return qubit == 0;
    }

    void SetOperation(const MatrixClass &U)
    {
        assert(U.rows() == 2 && U.cols() == 2);
        OpClass::operatorMat.block(2, 2, 2, 2) = U;
    }
};

// also named controlled X gate
template <class MatrixClass = Eigen::MatrixXcd> class CNOTGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::X, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    CNOTGate()
    {
        /*
        MatrixClass U = MatrixClass::Zero(2, 2);

        U(0, 1) = 1;
        U(1, 0) = 1;

        TwoQubitsControlledGate<MatrixClass>::SetOperation(U);
        */

        // the commented code above is correct, but this one should be faster
        OpClass::operatorMat(2, 2) = 0;
        OpClass::operatorMat(3, 3) = 0;
        OpClass::operatorMat(2, 3) = 1;
        OpClass::operatorMat(3, 2) = 1;
    }

    // this reffers only to the controlled block, not the whole matrix
    bool isAntidiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class ControlledYGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Y, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledYGate()
    {
        OpClass::operatorMat(2, 2) = 0;
        OpClass::operatorMat(3, 3) = 0;
        OpClass::operatorMat(2, 3) = std::complex(0., -1.);
        OpClass::operatorMat(3, 2) = std::complex(0., 1.);
    }

    // this reffers only to the controlled block, not the whole matrix
    bool isAntidiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class ControlledZGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Z, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledZGate()
    {
        OpClass::operatorMat(3, 3) = -1;
    }

    bool isDiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd>
class ControlledHadamardGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::H, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledHadamardGate()
    {
        static const double norm = 1. / sqrt(2.);
        OpClass::operatorMat(2, 2) = norm;
        OpClass::operatorMat(2, 3) = norm;
        OpClass::operatorMat(3, 2) = norm;
        OpClass::operatorMat(3, 3) = -norm;
    }
};

template <class MatrixClass = Eigen::MatrixXcd>
class ControlledSquareRootNOTGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Sx, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledSquareRootNOTGate()
    {
        OpClass::operatorMat(2, 2) = std::complex<double>(0.5, 0.5);
        OpClass::operatorMat(2, 3) = std::complex<double>(0.5, -0.5);
        OpClass::operatorMat(3, 2) = OpClass::operatorMat(2, 3);
        OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
    }
};

template <class MatrixClass = Eigen::MatrixXcd>
class ControlledSquareRootNOTDagGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Sxdg, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledSquareRootNOTDagGate()
    {
        OpClass::operatorMat(2, 2) = std::complex<double>(0.5, -0.5);
        OpClass::operatorMat(2, 3) = std::complex<double>(0.5, 0.5);
        OpClass::operatorMat(3, 2) = OpClass::operatorMat(2, 3);
        OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
    }
};

// controlled-V gate
template <class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::S, {}, 1};
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledPhaseGate()
    {
        OpClass::operatorMat(3, 3) = std::complex<double>(0, 1);
    }

    bool isDiagonal() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd>
class ControlledPhaseShiftGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Phase, {{recordedTheta}}, 1,
                                            {{recordedPhase}});
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledPhaseShiftGate(double theta = 0)
    {
        SetPhaseShift(theta);
    }

    void SetPhaseShift(double theta)
    {
        recordedTheta = theta;
        recordedPhase = UnitaryPhase(theta);
        OpClass::operatorMat(3, 3) = recordedPhase;
    }

    bool isDiagonal() const override
    {
        return true;
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedPhase = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> class ControlledUGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::U, recordedAngles, 1, recordedCoefficients);
    }

    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledUGate(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
    {
        SetParams(theta, phi, lambda, gamma);
    }

    virtual void SetParams(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
    {
        recordedAngles = {theta, phi, lambda, gamma};
        const auto coefficients = UnitaryEulerCoefficients(theta, phi, lambda, gamma);
        recordedCoefficients = coefficients;
        OpClass::operatorMat(2, 2) = coefficients[0];
        OpClass::operatorMat(2, 3) = coefficients[1];
        OpClass::operatorMat(3, 2) = coefficients[2];
        OpClass::operatorMat(3, 3) = coefficients[3];
    }

  private:
    std::array<double, 4> recordedAngles{};
    std::array<std::complex<double>, 4> recordedCoefficients{{1., 0., 0., 1.}};
};

template <class MatrixClass = Eigen::MatrixXcd>
class ControlledRotationGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    using BaseClass = TwoQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    virtual void SetTheta(double theta) = 0;
};

template <class MatrixClass = Eigen::MatrixXcd> class ControlledRxGate : public ControlledRotationGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Rx, {{recordedTheta}}, 1,
                                            {{recordedCosSin}});
    }

    using BaseClass = ControlledRotationGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledRxGate(double theta = 0)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta) override
    {
        recordedTheta = theta;
        double cosine, sine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        recordedCosSin = {cosine, sine};

        OpClass::operatorMat(2, 2) = std::complex<double>(cosine, 0);
        OpClass::operatorMat(2, 3) = std::complex<double>(0, -sine);
        OpClass::operatorMat(3, 2) = OpClass::operatorMat(2, 3);
        OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> class ControlledRyGate : public ControlledRotationGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Ry, {{recordedTheta}}, 1,
                                            {{recordedCosSin}});
    }

    using BaseClass = ControlledRotationGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledRyGate(double theta = 0)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta) override
    {
        recordedTheta = theta;
        double cosine, sine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        recordedCosSin = {cosine, sine};

        OpClass::operatorMat(2, 2) = std::complex<double>(cosine, 0);
        OpClass::operatorMat(2, 3) = std::complex<double>(-sine, 0);
        OpClass::operatorMat(3, 2) = std::complex<double>(sine, 0);
        OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> class ControlledRzGate : public ControlledRotationGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::Rz, {{recordedTheta}}, 1,
                                            {{recordedCosSin}});
    }

    using BaseClass = ControlledRotationGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ControlledRzGate(double theta = 0)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta) override
    {
        recordedTheta = theta;
        double cosine, sine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        recordedCosSin = {cosine, sine};

        OpClass::operatorMat(2, 2) = std::complex<double>(cosine, -sine);
        OpClass::operatorMat(3, 3) = std::complex<double>(cosine, sine);
    }

    bool isDiagonal() const override
    {
        return true;
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> class ThreeQubitsControlledGate : public ThreeQubitsGate<MatrixClass>
{
  public:
    using BaseClass = ThreeQubitsGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass;

    ThreeQubitsControlledGate() {};

    bool isControlled() const override
    {
        return true;
    }

    bool isControlQubit(size_t qubit) const override
    {
        return qubit == 0;
    }

    void SetOperation(const MatrixClass &U)
    {
        assert(U.rows() == 4 && U.cols() == 4);
        OpClass::operatorMat.block(4, 4, 4, 4) = U;
    }
};

// also known as CCNOT or CCX
template <class MatrixClass = Eigen::MatrixXcd> class ToffoliGate : public ThreeQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::X, {}, 2};
    }

    using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    ToffoliGate()
    {
        OpClass::operatorMat(6, 6) = 0;
        OpClass::operatorMat(7, 7) = 0;
        OpClass::operatorMat(6, 7) = 1;
        OpClass::operatorMat(7, 6) = 1;
    }

    bool isControlQubit(size_t qubit) const override
    {
        return qubit <= 1;
    }

    // this reffers only to the controlled block, not the whole matrix
    bool isAntidiagonal() const override
    {
        return true;
    }
};

// also known as CSWAP
template <class MatrixClass = Eigen::MatrixXcd> class FredkinGate : public ThreeQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Fredkin};
    }

    using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    FredkinGate()
    {
        OpClass::operatorMat(5, 5) = 0;
        OpClass::operatorMat(6, 6) = 0;
        OpClass::operatorMat(5, 6) = 1;
        OpClass::operatorMat(6, 5) = 1;
    }

    bool isSwapGate() const override
    {
        return true;
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class CCZGate : public ThreeQubitsControlledGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return {UnitaryGateDescriptor::Kind::Z, {}, 2};
    }

    using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
    using OpClass = typename BaseClass::BaseClass::BaseClass;

    CCZGate()
    {
        OpClass::operatorMat(7, 7) = -1;
    }

    bool isControlQubit(size_t qubit) const override
    {
        return qubit <= 1;
    }

    bool isDiagonal() const override
    {
        return true;
    }
};

// Local bit 0 is the first ApplyGate argument. For asymmetric interactions
// the first Pauli acts on that bit: RZX(theta, q0, q1) = exp(-i theta X1 Z0 / 2).
template <char First, char Second, class MatrixClass = Eigen::MatrixXcd>
class PauliProductRotationGate : public TwoQubitsGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(UnitaryGateDescriptor::Kind::PairRotation, {{recordedTheta}}, 0,
                                            {{recordedCosSin}}, First, Second);
    }

    explicit PauliProductRotationGate(double theta = 0.)
    {
        SetTheta(theta);
    }

    void SetTheta(double theta)
    {
        if (!std::isfinite(theta))
            throw std::invalid_argument("Gate angle must be finite");
        recordedTheta = theta;
        this->operatorMat.setZero();
        double c, s;
        detail::RotationHalfAngleCosSin(theta, c, s);
        recordedCosSin = {c, s};
        for (unsigned col = 0; col < 4; ++col)
        {
            unsigned row = col;
            std::complex<double> phase(0., -s);
            const char axes[] = {First, Second};
            for (unsigned q = 0; q < 2; ++q)
            {
                const bool bit = (col >> q) & 1U;
                if (axes[q] != 'Z')
                    row ^= 1U << q;
                if (axes[q] == 'Y')
                    phase *= std::complex<double>(0., bit ? -1. : 1.);
                else if (axes[q] == 'Z' && bit)
                    phase = -phase;
            }
            this->operatorMat(col, col) = c;
            this->operatorMat(row, col) += phase;
        }
    }

    bool isDiagonal() const override
    {
        return First == 'Z' && Second == 'Z';
    }

  private:
    double recordedTheta = 0.;
    std::complex<double> recordedCosSin = 1.;
};

template <class MatrixClass = Eigen::MatrixXcd> using RXXGate = PauliProductRotationGate<'X', 'X', MatrixClass>;
template <class MatrixClass = Eigen::MatrixXcd> using RYYGate = PauliProductRotationGate<'Y', 'Y', MatrixClass>;
template <class MatrixClass = Eigen::MatrixXcd> using RZZGate = PauliProductRotationGate<'Z', 'Z', MatrixClass>;
template <class MatrixClass = Eigen::MatrixXcd> using RZXGate = PauliProductRotationGate<'Z', 'X', MatrixClass>;
template <class MatrixClass = Eigen::MatrixXcd> using RXYGate = PauliProductRotationGate<'X', 'Y', MatrixClass>;

// exp[-i theta (XX +/- YY)/4], conjugated by the beta phase. Only the
// odd-parity (+) or even-parity (-) two-dimensional sector is changed.
template <bool Plus, class MatrixClass = Eigen::MatrixXcd> class XXYYGate : public TwoQubitsGate<MatrixClass>
{
  public:
    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return UnitaryGateDescriptor::Built(Plus ? UnitaryGateDescriptor::Kind::XXPlusYY
                                                 : UnitaryGateDescriptor::Kind::XXMinusYY,
                                            {{recordedTheta, recordedBeta}}, 0, recordedEntries);
    }

    explicit XXYYGate(double theta = 0., double beta = 0.)
    {
        SetParams(theta, beta);
    }

    void SetParams(double theta, double beta = 0.)
    {
        if (!std::isfinite(theta) || !std::isfinite(beta))
            throw std::invalid_argument("Gate angles must be finite");
        recordedTheta = theta;
        recordedBeta = beta;
        this->operatorMat.setIdentity();
        constexpr unsigned a = Plus ? 1 : 0, b = Plus ? 2 : 3;
        double cosine, sine, betaCosine, betaSine;
        detail::RotationHalfAngleCosSin(theta, cosine, sine);
        detail::QuarterTurnCosSin(beta, betaCosine, betaSine);
        const std::complex<double> s(0., -sine);
        recordedEntries = {cosine, s * std::complex<double>(betaCosine, -betaSine),
                           s * std::complex<double>(betaCosine, betaSine), 0.};
        this->operatorMat(a, a) = this->operatorMat(b, b) = cosine;
        this->operatorMat(a, b) = recordedEntries[1];
        this->operatorMat(b, a) = recordedEntries[2];
    }

  private:
    double recordedTheta = 0., recordedBeta = 0.;
    // cosine and the two off-diagonal entries of the active sector, as stored
    std::array<std::complex<double>, 4> recordedEntries{{1., 0., 0., 0.}};
};
template <class MatrixClass = Eigen::MatrixXcd> using XXPlusYYGate = XXYYGate<true, MatrixClass>;
template <class MatrixClass = Eigen::MatrixXcd> using XXMinusYYGate = XXYYGate<false, MatrixClass>;

template <class MatrixClass = Eigen::MatrixXcd> class ECRGate : public TwoQubitsGate<MatrixClass>
{
  public:
    ECRGate()
    {
        this->operatorMat.setZero();
        const double r = std::sqrt(.5);
        for (unsigned col = 0; col < 4; ++col)
        {
            this->operatorMat(col ^ 1U, col) = r;
            this->operatorMat(col ^ 3U, col) = std::complex<double>(0., (col & 1U) ? r : -r);
        }
    }
};

template <class MatrixClass = Eigen::MatrixXcd> class DCXGate : public TwoQubitsGate<MatrixClass>
{
  public:
    DCXGate()
    {
        this->operatorMat.setZero();
        for (unsigned col = 0; col < 4; ++col)
        {
            unsigned row = col ^ ((col & 1U) << 1);
            row ^= (row >> 1) & 1U;
            this->operatorMat(row, col) = 1.;
        }
    }
};

// Generic controlled forms preserve the single-qubit matrix's global phase:
// a phase that is global before adding controls becomes a relative phase.
template <class MatrixClass = Eigen::MatrixXcd> class ControlledGate : public TwoQubitsControlledGate<MatrixClass>
{
  public:
    ControlledGate() = default;

    explicit ControlledGate(const QuantumGateWithOp<MatrixClass> &gate)
    {
        SetOperation(gate);
    }

    explicit ControlledGate(const MatrixClass &matrix)
    {
        SetOperation(matrix);
    }

    void SetOperation(const MatrixClass &matrix)
    {
        if (matrix.rows() != 2 || matrix.cols() != 2 || !matrix.allFinite())
            throw std::invalid_argument("Controlled gate requires a finite 2 by 2 matrix");
        descriptor = {};
        this->operatorMat.setIdentity();
        this->operatorMat.block(2, 2, 2, 2) = matrix;
    }

    bool isDiagonal() const override
    {
        return this->operatorMat(2, 3) == std::complex<double>(0.) &&
               this->operatorMat(3, 2) == std::complex<double>(0.);
    }

    // As for the legacy controlled gates, this describes the active block.
    bool isAntidiagonal() const override
    {
        return this->operatorMat(2, 2) == std::complex<double>(0.) &&
               this->operatorMat(3, 3) == std::complex<double>(0.);
    }

    template <class Gate, class = decltype(std::declval<const Gate &>().getRawOperatorMatrix())>
    explicit ControlledGate(const Gate &gate)
    {
        SetOperation(gate);
    }

    template <class Gate> auto SetOperation(const Gate &gate) -> decltype(gate.getRawOperatorMatrix(), void())
    {
        SetOperation(gate.getRawOperatorMatrix());
        auto candidate = GetUnitaryGateDescriptor(gate);
        if (candidate && candidate.controls == 0)
        {
            candidate.controls = 1;
            descriptor = candidate;
        }
    }

    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return descriptor;
    }

  private:
    UnitaryGateDescriptor descriptor;
};

template <class MatrixClass = Eigen::MatrixXcd>
class DoublyControlledGate : public ThreeQubitsControlledGate<MatrixClass>
{
  public:
    DoublyControlledGate() = default;

    explicit DoublyControlledGate(const QuantumGateWithOp<MatrixClass> &gate)
    {
        SetOperation(gate);
    }

    explicit DoublyControlledGate(const MatrixClass &matrix)
    {
        SetOperation(matrix);
    }

    void SetOperation(const MatrixClass &matrix)
    {
        if (matrix.rows() != 2 || matrix.cols() != 2 || !matrix.allFinite())
            throw std::invalid_argument("Doubly controlled gate requires a finite 2 by 2 matrix");
        descriptor = {};
        this->operatorMat.setIdentity();
        this->operatorMat.block(6, 6, 2, 2) = matrix;
    }

    bool isControlQubit(size_t qubit) const override
    {
        return qubit < 2;
    }

    bool isDiagonal() const override
    {
        return this->operatorMat(6, 7) == std::complex<double>(0.) &&
               this->operatorMat(7, 6) == std::complex<double>(0.);
    }

    bool isAntidiagonal() const override
    {
        return this->operatorMat(6, 6) == std::complex<double>(0.) &&
               this->operatorMat(7, 7) == std::complex<double>(0.);
    }

    template <class Gate, class = decltype(std::declval<const Gate &>().getRawOperatorMatrix())>
    explicit DoublyControlledGate(const Gate &gate)
    {
        SetOperation(gate);
    }

    template <class Gate> auto SetOperation(const Gate &gate) -> decltype(gate.getRawOperatorMatrix(), void())
    {
        SetOperation(gate.getRawOperatorMatrix());
        auto candidate = GetUnitaryGateDescriptor(gate);
        if (candidate && candidate.controls == 0)
        {
            candidate.controls = 2;
            descriptor = candidate;
        }
    }

    UnitaryGateDescriptor getUnitaryDescriptor() const override
    {
        return descriptor;
    }

  private:
    UnitaryGateDescriptor descriptor;
};

template <template <class> class SingleGate, bool Double, class MatrixClass = Eigen::MatrixXcd>
class ControlledSingleQubitGate
    : public std::conditional_t<Double, DoublyControlledGate<MatrixClass>, ControlledGate<MatrixClass>>
{
  public:
    using BaseClass = std::conditional_t<Double, DoublyControlledGate<MatrixClass>, ControlledGate<MatrixClass>>;

    template <class... Args>
    explicit ControlledSingleQubitGate(Args... args) : BaseClass(SingleGate<MatrixClass>(args...))
    {
    }

    void SetTheta(double theta)
    {
        this->SetOperation(SingleGate<MatrixClass>(theta));
    }

    void SetPhaseShift(double theta)
    {
        this->SetOperation(SingleGate<MatrixClass>(theta));
    }
};

template <class M = Eigen::MatrixXcd> using ControlledHyGate = ControlledSingleQubitGate<HyGate, false, M>;
template <class M = Eigen::MatrixXcd> using ControlledSGate = ControlledSingleQubitGate<SGate, false, M>;
template <class M = Eigen::MatrixXcd> using ControlledSDGGate = ControlledSingleQubitGate<SDGGate, false, M>;
template <class M = Eigen::MatrixXcd> using ControlledTGate = ControlledSingleQubitGate<TGate, false, M>;
template <class M = Eigen::MatrixXcd> using ControlledTDGGate = ControlledSingleQubitGate<TDGGate, false, M>;
template <class M = Eigen::MatrixXcd> using ControlledSplitterGate = ControlledSingleQubitGate<SplitterGate, false, M>;
template <class M = Eigen::MatrixXcd> using CSGate = ControlledSGate<M>;
template <class M = Eigen::MatrixXcd> using CSDGGate = ControlledSDGGate<M>;
template <class M = Eigen::MatrixXcd> using CCXGate = ToffoliGate<M>;
template <class M = Eigen::MatrixXcd> using CCYGate = ControlledSingleQubitGate<PauliYGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCHGate = ControlledSingleQubitGate<HadamardGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCHyGate = ControlledSingleQubitGate<HyGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCSGate = ControlledSingleQubitGate<SGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCSDGGate = ControlledSingleQubitGate<SDGGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCTGate = ControlledSingleQubitGate<TGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCTDGGate = ControlledSingleQubitGate<TDGGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCSXGate = ControlledSingleQubitGate<SquareRootNOTGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCSXDGGate = ControlledSingleQubitGate<SquareRootNOTDagGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCSplitterGate = ControlledSingleQubitGate<SplitterGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCRxGate = ControlledSingleQubitGate<RxGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCRyGate = ControlledSingleQubitGate<RyGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCRzGate = ControlledSingleQubitGate<RzGate, true, M>;
template <class M = Eigen::MatrixXcd> using CCPhaseShiftGate = ControlledSingleQubitGate<PhaseShiftGate, true, M>;

template <class M = Eigen::MatrixXcd> class CCUGate : public DoublyControlledGate<M>
{
  public:
    CCUGate(double theta = 0., double phi = 0., double lambda = 0., double gamma = 0.)
    {
        SetParams(theta, phi, lambda, gamma);
    }

    void SetParams(double theta, double phi = 0., double lambda = 0., double gamma = 0.)
    {
        if (!std::isfinite(theta) || !std::isfinite(phi) || !std::isfinite(lambda) || !std::isfinite(gamma))
            throw std::invalid_argument("Gate angles must be finite");
        const UGate<M> gate(theta, phi, lambda, gamma);
        this->SetOperation(gate);
    }
};

} // namespace Gates
} // namespace QC
