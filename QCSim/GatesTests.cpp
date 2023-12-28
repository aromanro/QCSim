#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>
#include <iterator>
#include <map>

#include "Tests.h"
#include "QuantumGate.h"

bool checkSingleQubitGates()
{
	QC::Gates::HadamardGate hadamard;
	if (!checkUnitary(hadamard.getRawOperatorMatrix())) {
		std::cout << "The Hadamard gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::HyGate hy;
	if (!checkUnitary(hy.getRawOperatorMatrix())) {
		std::cout << "The Hy gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::PhaseGate phase;
	if (!checkUnitary(phase.getRawOperatorMatrix())) {
		std::cout << "The Phase gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::PhaseShiftGate phaseShift;
	if (!checkUnitary(phaseShift.getRawOperatorMatrix())) {
		std::cout << "The Phase Shift gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::PauliXGate x;
	if (!checkUnitary(x.getRawOperatorMatrix())) {
		std::cout << "The Pauli X gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::PauliYGate y;
	if (!checkUnitary(y.getRawOperatorMatrix())) {
		std::cout << "The Pauli Y gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::PauliZGate z;
	if (!checkUnitary(z.getRawOperatorMatrix())) {
		std::cout << "The Pauli Z gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::SquareRootNOTGate notSquareRoot;
	if (!checkUnitary(notSquareRoot.getRawOperatorMatrix())) {
		std::cout << "The Squared Root NOT gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::SplitterGate splitter;
	if (!checkUnitary(splitter.getRawOperatorMatrix())) {
		std::cout << "The Splitter gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::RxGate rx;
	if (!checkUnitary(rx.getRawOperatorMatrix())) {
		std::cout << "The Rx gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::RyGate ry;
	if (!checkUnitary(ry.getRawOperatorMatrix())) {
		std::cout << "The Ry gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::RzGate rz;
	if (!checkUnitary(rz.getRawOperatorMatrix())) {
		std::cout << "The Rz gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::UGate u;
	u.SetParams(M_PI / 4, M_PI / 5, M_PI / 7);
	if (!checkUnitary(u.getRawOperatorMatrix())) {
		std::cout << "The U gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkDoubleQubitGates()
{
	QC::Gates::SwapGate swap;
	if (!checkUnitary(swap.getRawOperatorMatrix())) {
		std::cout << "The Swap gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::iSwapGate iswap;
	if (!checkUnitary(iswap.getRawOperatorMatrix())) {
		std::cout << "The iSwap gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::DecrementGate decrement;
	if (!checkUnitary(decrement.getRawOperatorMatrix())) {
		std::cout << "The Decrement gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::CNOTGate cnot;
	if (!checkUnitary(cnot.getRawOperatorMatrix())) {
		std::cout << "The CNOT gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::ControlledPhaseGate cphase;
	if (!checkUnitary(cphase.getRawOperatorMatrix())) {
		std::cout << "The Controlled Phase gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::ControlledPhaseShiftGate cphaseShift;
	if (!checkUnitary(cphaseShift.getRawOperatorMatrix())) {
		std::cout << "The Controlled Phase Shift gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::ControlledZGate cz;
	if (!checkUnitary(cz.getRawOperatorMatrix())) {
		std::cout << "The Controlled Z gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkTripleQubitGates()
{
	QC::Gates::ToffoliGate toffoli;
	if (!checkUnitary(toffoli.getRawOperatorMatrix())) {
		std::cout << "The Toffoli gate is not unitary!" << std::endl;
		return false;
	}

	QC::Gates::FredkinGate fredkin;
	if (!checkUnitary(fredkin.getRawOperatorMatrix())) {
		std::cout << "The Fredkin gate is not unitary!" << std::endl;
		return false;
	}

	return true;
}

bool checkGates()
{
	std::cout << "\nChecking unitarity of the quantum gates..." << std::endl;

	if (!checkSingleQubitGates() || !checkDoubleQubitGates() || !checkTripleQubitGates()) return false;

	std::cout << "OK!" << std::endl;

	return true;
}
