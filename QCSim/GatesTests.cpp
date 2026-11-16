#define _USE_MATH_DEFINES
#include <math.h>

#include <array>
#include <iostream>
#include <iterator>
#include <map>
#include <vector>

#include "AppliedUnitaryRegressionTests.h"
#include "DensityMatrix.h"
#include "ExtendedUnitaryRegressionTests.h"
#include "GateCoverageRegressionTests.h"
#include "PathIntegral.h"
#include "PauliGatePrecisionRegressionTests.h"
#include "QuantumGate.h"
#include "Tests.h"

bool checkSingleQubitGates1()
{
    if (QC::Gates::HadamardGate hadamard; !checkUnitary(hadamard.getRawOperatorMatrix()))
    {
        std::cout << "The Hadamard gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::HyGate hy; !checkUnitary(hy.getRawOperatorMatrix()))
    {
        std::cout << "The Hy gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::SGate phase; !checkUnitary(phase.getRawOperatorMatrix()))
    {
        std::cout << "The Phase/S gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::SDGGate sdg; !checkUnitary(sdg.getRawOperatorMatrix()))
    {
        std::cout << "The SDG gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::TGate t; !checkUnitary(t.getRawOperatorMatrix()))
    {
        std::cout << "The T gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::TDGGate tdg; !checkUnitary(tdg.getRawOperatorMatrix()))
    {
        std::cout << "The TDG gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::PhaseShiftGate phaseShift(M_PI / 4); !checkUnitary(phaseShift.getRawOperatorMatrix()))
    {
        std::cout << "The Phase Shift gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::PauliXGate x; !checkUnitary(x.getRawOperatorMatrix()))
    {
        std::cout << "The Pauli X gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::PauliYGate y; !checkUnitary(y.getRawOperatorMatrix()))
    {
        std::cout << "The Pauli Y gate is not unitary!" << std::endl;
        return false;
    }

    return true;
}

bool checkSingleQubitGates2()
{
    if (QC::Gates::PauliZGate z; !checkUnitary(z.getRawOperatorMatrix()))
    {
        std::cout << "The Pauli Z gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::SquareRootNOTGate notSquareRoot; !checkUnitary(notSquareRoot.getRawOperatorMatrix()))
    {
        std::cout << "The Squared Root NOT gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::SquareRootNOTDagGate notSquareRootDag; !checkUnitary(notSquareRootDag.getRawOperatorMatrix()))
    {
        std::cout << "The Squared Root NOT Dagger gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::SplitterGate splitter; !checkUnitary(splitter.getRawOperatorMatrix()))
    {
        std::cout << "The Splitter gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::RxGate rx(M_PI / 4); !checkUnitary(rx.getRawOperatorMatrix()))
    {
        std::cout << "The Rx gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::RyGate ry(M_PI / 4); !checkUnitary(ry.getRawOperatorMatrix()))
    {
        std::cout << "The Ry gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::RzGate rz(M_PI / 4); !checkUnitary(rz.getRawOperatorMatrix()))
    {
        std::cout << "The Rz gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::UGate u(M_PI / 4, M_PI / 5, M_PI / 7); !checkUnitary(u.getRawOperatorMatrix()))
    {
        std::cout << "The U gate is not unitary!" << std::endl;
        return false;
    }

    return true;
}

bool checkDoubleQubitGates1()
{
    if (QC::Gates::SwapGate swap; !checkUnitary(swap.getRawOperatorMatrix()))
    {
        std::cout << "The Swap gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::iSwapGate iswap; !checkUnitary(iswap.getRawOperatorMatrix()))
    {
        std::cout << "The iSwap gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::DecrementGate decrement; !checkUnitary(decrement.getRawOperatorMatrix()))
    {
        std::cout << "The Decrement gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::CNOTGate cnot; !checkUnitary(cnot.getRawOperatorMatrix()))
    {
        std::cout << "The CNOT gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledYGate cy; !checkUnitary(cy.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Y gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledPhaseGate cphase; !checkUnitary(cphase.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Phase gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledPhaseShiftGate cphaseShift(M_PI / 3); !checkUnitary(cphaseShift.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Phase Shift gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledZGate cz; !checkUnitary(cz.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Z gate is not unitary!" << std::endl;
        return false;
    }

    return true;
}

bool checkDoubleQubitGates2()
{
    if (QC::Gates::ControlledHadamardGate ch; !checkUnitary(ch.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Hadamard gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledSquareRootNOTGate cnotSquareRoot; !checkUnitary(cnotSquareRoot.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Squared Root NOT gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledSquareRootNOTDagGate cnotSquareRootDag;
        !checkUnitary(cnotSquareRootDag.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Squared Root NOT Dagger gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledUGate<> cu(M_PI / 4, M_PI / 3); !checkUnitary(cu.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled U gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledRxGate crx(M_PI / 5); !checkUnitary(crx.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Rx gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledRyGate cry(M_PI / 3); !checkUnitary(cry.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Ry gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::ControlledRzGate crz(M_PI / 4); !checkUnitary(crz.getRawOperatorMatrix()))
    {
        std::cout << "The Controlled Rz gate is not unitary!" << std::endl;
        return false;
    }

    return true;
}

bool checkTripleQubitGates()
{
    if (QC::Gates::ToffoliGate toffoli; !checkUnitary(toffoli.getRawOperatorMatrix()))
    {
        std::cout << "The Toffoli gate is not unitary!" << std::endl;
        return false;
    }

    if (QC::Gates::FredkinGate fredkin; !checkUnitary(fredkin.getRawOperatorMatrix()))
    {
        std::cout << "The Fredkin gate is not unitary!" << std::endl;
        return false;
    }

    return true;
}

bool checkDoubleQubitGates()
{
    return checkDoubleQubitGates1() && checkDoubleQubitGates2();
}

bool checkSingleQubitGates()
{
    return checkSingleQubitGates1() && checkSingleQubitGates2();
}

bool checkExpandedLocalGateBackends()
{
    using namespace QC;
    using namespace QC::Gates;
    std::vector<AppliedGate<>> preparation;
    for (size_t q = 0; q < 4; ++q)
        preparation.emplace_back(UGate<>(0.3 + q * 0.17, -0.7 + q * 0.19, 0.2).getRawOperatorMatrix(), q);
    preparation.emplace_back(CNOTGate<>().getRawOperatorMatrix(), 3, 0);
    Eigen::VectorXcd initial = Eigen::VectorXcd::Zero(16);
    initial[0] = 1.;
    for (const auto &gate : preparation)
    {
        std::vector<unsigned> qubits{static_cast<unsigned>(gate.getQubit1())};
        if (gate.getQubitsNumber() > 1)
            qubits.push_back(static_cast<unsigned>(gate.getQubit2()));
        initial = GateCoverageRegression::Apply(initial, gate.getRawOperatorMatrix(), qubits);
    }
    const std::vector<Eigen::MatrixXcd> matrices{
        RZZGate<>(0.37).getRawOperatorMatrix(),
        RXXGate<>(-0.21).getRawOperatorMatrix(),
        RYYGate<>(0.53).getRawOperatorMatrix(),
        RZXGate<>(-0.79).getRawOperatorMatrix(),
        RXYGate<>(0.43).getRawOperatorMatrix(),
        XXPlusYYGate<>(0.51, -0.29).getRawOperatorMatrix(),
        XXMinusYYGate<>(-0.63, 0.31).getRawOperatorMatrix(),
        ECRGate<>().getRawOperatorMatrix(),
        DCXGate<>().getRawOperatorMatrix(),
        CSGate<>().getRawOperatorMatrix(),
        CSDGGate<>().getRawOperatorMatrix(),
        iSwapGate<>().getRawOperatorMatrix(),
        CCUGate<>(0.43, -0.71, 0.19, 0.23).getRawOperatorMatrix(),
        DoublyControlledGate<>(HyGate<>()).getRawOperatorMatrix(),
        DoublyControlledGate<>(SquareRootNOTDagGate<>()).getRawOperatorMatrix()};
    for (const auto &matrix : matrices)
    {
        std::array<unsigned, 3> positions{0, 2, 3};
        do
        {
            const size_t arity = matrix.rows() == 8 ? 3 : 2;
            const std::vector<unsigned> qubits(positions.begin(), positions.begin() + arity);
            const Eigen::VectorXcd expected = GateCoverageRegression::Apply(initial, matrix, qubits);
            const AppliedGate<> applied(matrix, positions[0], positions[1], positions[2]);
            QubitRegister<> state(4, 31);
            state.setRegisterStorage(initial);
            state.ApplyGate(applied);
            DensityMatrix<> density(4, 31);
            density.setDensityMatrix((initial * initial.adjoint()).eval());
            density.ApplyGate(applied);
            if ((state.getRegisterStorage() - expected).norm() > 1E-12 ||
                (density.getDensityMatrix() - expected * expected.adjoint()).norm() > 1E-12)
                return false;
            auto circuit = preparation;
            circuit.push_back(applied);
            PathIntegral::PathIntegralSimulator path;
            path.SetTrimValue(0.);
            path.SetCircuit(circuit);
            for (size_t output = 0; output < 16; ++output)
            {
                std::vector<bool> bits(4);
                for (size_t q = 0; q < 4; ++q)
                    bits[q] = ((output >> q) & 1) != 0;
                if (std::abs(path.Propagate(bits) - expected[output]) > 1E-12)
                    return false;
            }
        } while (std::next_permutation(positions.begin(), positions.end()));
    }
    std::cout << "Expanded gate statevector, density-matrix and path-integral comparisons passed.\n";
    return true;
}

bool checkGates()
{
    std::cout << "\nChecking unitarity of the quantum gates..." << std::endl;

    if (!checkSingleQubitGates() || !checkDoubleQubitGates() || !checkTripleQubitGates())
        return false;

    if (!QC::GateCoverageRegression::Run())
        return false;

    try
    {
        QC::PauliGatePrecisionRegression::Run();
        QC::AppliedUnitaryRegression::Run();
        QC::ExtendedUnitaryRegression::Run();
    }
    catch (const std::exception &error)
    {
        std::cerr << "Pauli/extended gate precision regression: " << error.what() << '\n';
        return false;
    }

    if (!checkExpandedLocalGateBackends())
        return false;

    std::cout << "OK!" << std::endl;

    return true;
}
