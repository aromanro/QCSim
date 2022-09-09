#pragma once

#include <Eigen/eigen>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {

	class QuantumGate
	{
	public:
		virtual Eigen::MatrixXcd getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const = 0;
	};

	class QuantumGateWithOp : public QuantumGate
	{
	protected:
		Eigen::MatrixXcd operatorMat;
	};

	class SingleQubitGate : public QuantumGateWithOp
	{
	public:
		SingleQubitGate();

		// controllingQubit is ignored, it will be used for two qubit gates
		Eigen::MatrixXcd getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override;
	};

	class HadamardGate : public SingleQubitGate
	{
	public:
		HadamardGate();
	};

	class PhaseShiftGate : public SingleQubitGate
	{
	public:
		PhaseShiftGate(double theta = 0);

		void SetPhaseShift(double theta);
	};

	class PauliXGate : public SingleQubitGate
	{
	public:
		PauliXGate();
	};

	class PauliYGate : public SingleQubitGate
	{
	public:
		PauliYGate();
	};

	class PauliZGate : public SingleQubitGate
	{
	public:
		PauliZGate();
	};

	class TwoQubitsGate : public QuantumGateWithOp
	{
	public:
		TwoQubitsGate();
	};

	class SwapGate : public TwoQubitsGate
	{
	public:
		SwapGate();
	};

	class TwoQubitsControlledGate : public TwoQubitsGate
	{
	public:
		TwoQubitsControlledGate();

		void SetOperation(const Eigen::MatrixXcd& U);

		Eigen::MatrixXcd getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override;
	};

	class CNOTGate : public TwoQubitsControlledGate
	{
	public:
		CNOTGate();
	};

	class ControlledPhaseShiftGate : public TwoQubitsControlledGate
	{
	public:
		ControlledPhaseShiftGate(double theta = 0);

		void SetPhaseShift(double theta);
	};

	class ControlledZGate : public TwoQubitsControlledGate
	{
	public:
		ControlledZGate();
	};
}


