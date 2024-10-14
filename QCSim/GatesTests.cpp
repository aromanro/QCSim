#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>
#include <iterator>
#include <map>

#include "Tests.h"
#include "QuantumGate.h"

bool checkSingleQubitGates()
{
	if (QC::Gates::HadamardGate hadamard; !checkUnitary(hadamard.getRawOperatorMatrix())) {
		std::cout << "The Hadamard gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::HyGate hy; !checkUnitary(hy.getRawOperatorMatrix())) {
		std::cout << "The Hy gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::SGate phase; !checkUnitary(phase.getRawOperatorMatrix())) {
		std::cout << "The Phase/S gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::SDGGate sdg; !checkUnitary(sdg.getRawOperatorMatrix())) {
		std::cout << "The SDG gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::TGate t; !checkUnitary(t.getRawOperatorMatrix())) {
		std::cout << "The T gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::TDGGate tdg; !checkUnitary(tdg.getRawOperatorMatrix())) {
		std::cout << "The TDG gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::PhaseShiftGate phaseShift(M_PI / 4); !checkUnitary(phaseShift.getRawOperatorMatrix())) {
		std::cout << "The Phase Shift gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::PauliXGate x; !checkUnitary(x.getRawOperatorMatrix())) {
		std::cout << "The Pauli X gate is not unitary!" << std::endl;
		return false;
	}
	
	if (QC::Gates::PauliYGate y; !checkUnitary(y.getRawOperatorMatrix())) {
		std::cout << "The Pauli Y gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::PauliZGate z; !checkUnitary(z.getRawOperatorMatrix())) {
		std::cout << "The Pauli Z gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::SquareRootNOTGate notSquareRoot; !checkUnitary(notSquareRoot.getRawOperatorMatrix())) {
		std::cout << "The Squared Root NOT gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::SquareRootNOTDagGate notSquareRootDag; !checkUnitary(notSquareRootDag.getRawOperatorMatrix())) {
		std::cout << "The Squared Root NOT Dagger gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::SplitterGate splitter; !checkUnitary(splitter.getRawOperatorMatrix())) {
		std::cout << "The Splitter gate is not unitary!" << std::endl;
		return false;
	}
	
	if (QC::Gates::RxGate rx(M_PI / 4); !checkUnitary(rx.getRawOperatorMatrix())) {
		std::cout << "The Rx gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::RyGate ry(M_PI / 4); !checkUnitary(ry.getRawOperatorMatrix())) {
		std::cout << "The Ry gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::RzGate rz(M_PI / 4); !checkUnitary(rz.getRawOperatorMatrix())) {
		std::cout << "The Rz gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::UGate u(M_PI / 4, M_PI / 5, M_PI / 7); !checkUnitary(u.getRawOperatorMatrix())) {
		std::cout << "The U gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkDoubleQubitGates1()
{
	if (QC::Gates::SwapGate swap; !checkUnitary(swap.getRawOperatorMatrix())) {
		std::cout << "The Swap gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::iSwapGate iswap; !checkUnitary(iswap.getRawOperatorMatrix())) {
		std::cout << "The iSwap gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::DecrementGate decrement; !checkUnitary(decrement.getRawOperatorMatrix())) {
		std::cout << "The Decrement gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::CNOTGate cnot; !checkUnitary(cnot.getRawOperatorMatrix())) {
		std::cout << "The CNOT gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledYGate cy; !checkUnitary(cy.getRawOperatorMatrix())) {
		std::cout << "The Controlled Y gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledPhaseGate cphase; !checkUnitary(cphase.getRawOperatorMatrix())) {
		std::cout << "The Controlled Phase gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledPhaseShiftGate cphaseShift(M_PI / 3); !checkUnitary(cphaseShift.getRawOperatorMatrix())) {
		std::cout << "The Controlled Phase Shift gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledZGate cz; !checkUnitary(cz.getRawOperatorMatrix())) {
		std::cout << "The Controlled Z gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkDoubleQubitGates2()
{
	if (QC::Gates::ControlledHadamardGate ch; !checkUnitary(ch.getRawOperatorMatrix())) {
		std::cout << "The Controlled Hadamard gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledSquareRootNOTGate cnotSquareRoot; !checkUnitary(cnotSquareRoot.getRawOperatorMatrix())) {
		std::cout << "The Controlled Squared Root NOT gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledSquareRootNOTDagGate cnotSquareRootDag; !checkUnitary(cnotSquareRootDag.getRawOperatorMatrix())) {
		std::cout << "The Controlled Squared Root NOT Dagger gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledUGate<> cu(M_PI / 4, M_PI / 3); !checkUnitary(cu.getRawOperatorMatrix())) {
		std::cout << "The Controlled U gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledRxGate crx(M_PI / 5); !checkUnitary(crx.getRawOperatorMatrix())) {
		std::cout << "The Controlled Rx gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledRyGate cry(M_PI / 3); !checkUnitary(cry.getRawOperatorMatrix())) {
		std::cout << "The Controlled Ry gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::ControlledRzGate crz(M_PI / 4); !checkUnitary(crz.getRawOperatorMatrix())) {
		std::cout << "The Controlled Rz gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkTripleQubitGates()
{
	if (QC::Gates::ToffoliGate toffoli; !checkUnitary(toffoli.getRawOperatorMatrix())) {
		std::cout << "The Toffoli gate is not unitary!" << std::endl;
		return false;
	}

	if (QC::Gates::FredkinGate fredkin; !checkUnitary(fredkin.getRawOperatorMatrix())) {
		std::cout << "The Fredkin gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkDoubleQubitGates()
{
	return checkDoubleQubitGates1() && checkDoubleQubitGates2();
}

bool checkGates()
{
	std::cout << "\nChecking unitarity of the quantum gates..." << std::endl;

	if (!checkSingleQubitGates() || !checkDoubleQubitGates() || !checkTripleQubitGates()) return false;

	std::cout << "OK!" << std::endl;

	return true;
}
