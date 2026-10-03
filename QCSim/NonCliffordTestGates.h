#pragma once

// The statevector side of the randomized non-Clifford gate tests, shared by the
// extended stabilizer and Pauli propagator suites.

#include <cmath>
#include <random>

#include "Tests.h"
#include "QubitRegister.h"

// Codes 0-14 are Clifford gates and 15-17 rotations, as in GetGate; 18-30 are
// the other supported non-Clifford gates. qubit1 is the target (the first
// target of CSwap), qubit2 the control (the second target of CSwap) and qubit3
// the second control (the control of CSwap). U and CU derive their further
// angles from angle.
enum NonCliffordGateCode : int {
	CodeRx = 15, CodeU = 18, CodeCU, CodeCRx, CodeCRy, CodeCRz, CodeCP, CodeCS,
	CodeCSdg, CodeCSx, CodeCSxDag, CodeCH, CodeCCX, CodeCSwap
};

struct NonCliffordGateAngles {
	double theta, phi, lambda, gamma;
};

inline NonCliffordGateAngles DeriveGateAngles(double angle)
{
	return { angle, 0.7 * angle + 0.3, 0.5 - 0.4 * angle, 0.25 * angle };
}

inline void ApplyStatevectorGate(QC::QubitRegister<>& qubitRegister, int code,
	size_t qubit1, size_t qubit2, double angle, size_t qubit3)
{
	if (code < CodeU)
	{
		auto gate = GetGate(code, angle);
		qubitRegister.ApplyGate(*gate, qubit1, qubit2);
		return;
	}
	const auto angles = DeriveGateAngles(angle);
	const double pi = std::acos(-1.0);
	switch (code)
	{
	case CodeU:
		qubitRegister.ApplyGate(QC::Gates::UGate<>(angles.theta, angles.phi, angles.lambda), qubit1);
		break;
	case CodeCU:
		qubitRegister.ApplyGate(QC::Gates::ControlledUGate<>(angles.theta, angles.phi,
			angles.lambda, angles.gamma), qubit1, qubit2);
		break;
	case CodeCRx: qubitRegister.ApplyGate(QC::Gates::ControlledRxGate<>(angle), qubit1, qubit2); break;
	case CodeCRy: qubitRegister.ApplyGate(QC::Gates::ControlledRyGate<>(angle), qubit1, qubit2); break;
	case CodeCRz: qubitRegister.ApplyGate(QC::Gates::ControlledRzGate<>(angle), qubit1, qubit2); break;
	case CodeCP: qubitRegister.ApplyGate(QC::Gates::ControlledPhaseShiftGate<>(angle), qubit1, qubit2); break;
	case CodeCS: qubitRegister.ApplyGate(QC::Gates::ControlledPhaseShiftGate<>(0.5 * pi), qubit1, qubit2); break;
	case CodeCSdg: qubitRegister.ApplyGate(QC::Gates::ControlledPhaseShiftGate<>(-0.5 * pi), qubit1, qubit2); break;
	case CodeCSx: qubitRegister.ApplyGate(QC::Gates::ControlledSquareRootNOTGate<>(), qubit1, qubit2); break;
	case CodeCSxDag: qubitRegister.ApplyGate(QC::Gates::ControlledSquareRootNOTDagGate<>(), qubit1, qubit2); break;
	case CodeCH: qubitRegister.ApplyGate(QC::Gates::ControlledHadamardGate<>(), qubit1, qubit2); break;
	case CodeCCX: qubitRegister.ApplyGate(QC::Gates::ToffoliGate<>(), qubit1, qubit2, qubit3); break;
	default: qubitRegister.ApplyGate(QC::Gates::FredkinGate<>(), qubit1, qubit2, qubit3); break;
	}
}

// Any supported non-Clifford gate that fits the register: a rotation, U, a
// controlled two-qubit gate or a three-qubit gate.
inline int RandomNonCliffordCode(std::mt19937& generator, size_t nrQubits)
{
	const int last = nrQubits >= 3 ? CodeCSwap : (nrQubits == 2 ? CodeCH : CodeU);
	return std::uniform_int_distribution<int>(CodeRx, last)(generator);
}

// A third qubit distinct from two others in a register of at least three.
inline size_t ThirdQubit(size_t nrQubits, size_t qubit1, size_t qubit2)
{
	size_t qubit3 = (qubit2 + 1) % nrQubits;
	while (qubit3 == qubit1 || qubit3 == qubit2) qubit3 = (qubit3 + 1) % nrQubits;
	return qubit3;
}
