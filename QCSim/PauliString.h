#pragma once


namespace PauliString {

	// used by PauliDecomposedHamiltonianSimulation and Variational Quantum Eigensolver
	class PauliString {
	public:
		enum class PauliOp : unsigned char
		{
			opZ = 0,
			opX = 1,
			opY = 2
		};

		PauliString(int nrQubits) : ops(nrQubits, PauliOp::opZ), coeff(1.0) {}

		void setOperatorForQubit(unsigned int qubit, PauliOp op)
		{
			if (qubit >= ops.size()) return;

			ops[qubit] = op;
		}

		PauliOp getOperatorForQubit(unsigned int qubit) const
		{
			if (qubit >= ops.size()) return PauliOp::opZ;

			return ops[qubit];
		}

		double getCoefficient() const
		{
			return coeff;
		}

		void setCoefficient(double c)
		{
			coeff = c;
		}

	private:
		std::vector<PauliOp> ops;
		double coeff;
	};
}