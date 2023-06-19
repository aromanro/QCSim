#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NQubitsControlledQuantumGate
	{
	public:
		using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		NQubitsControlledQuantumGate(unsigned int N, const MatrixClass& op, unsigned int startQubit = 1, unsigned int controllingQubit = 0)
			: BaseClass(N, startQubit, startQubit + getOperatorQubitsNumber(op)), U(op), nrQubits(N), ctrlQubit(controllingQubit)
		{
		}

		unsigned int Execute(RegisterClass& reg) override
		{
			reg.ApplyOperatorMatrix(getOperatorMatrix());

			return 0; // no measurement, so the return should be ignored
		}

		MatrixClass getOperatorMatrix() const
		{
			const unsigned int startQubit = BaseClass::getStartQubit();
			const unsigned int endQubit = BaseClass::getEndQubit();
			const unsigned int nrBasisStates = 1u << nrQubits;
			const unsigned int ctrlQubitBit = 1u << ctrlQubit;

			MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);
			MatrixClass Unity = MatrixClass::Identity(U.rows(), U.cols());

			unsigned int mask = 0;
			for (unsigned int i = startQubit; i <= endQubit; ++i)
				mask |= 1u << i;

			for (unsigned int i = 0; i < nrBasisStates; ++i)
			{
				const unsigned int ind1 = i | mask;
				for (unsigned int j = 0; j < nrBasisStates; ++j)
					if (ind1 == (j | mask)) // the delta 'function'
					{
						const unsigned int ind1op = (i & mask) >> startQubit;
						const unsigned int ind2op = (j & mask) >> startQubit;

						if ((i & ctrlQubitBit) == ctrlQubitBit && (j & ctrlQubitBit) == ctrlQubitBit)
							extOperatorMat(i, j) = U(ind1op, ind2op);
						else if ((i & ctrlQubitBit) == 0 && (j & ctrlQubitBit) == 0)
							extOperatorMat(i, j) = Unity(ind1op, ind2op);
						// else leave it zero, the off diagonal blocks are zero and the target matrix is already initialized to zero
					}
			}

			return extOperatorMat;
		}

		const MatrixClass& getRawOperatorMatrix() const
		{
			return U;
		}

		unsigned int getRawOperatorQubitsNumber() const
		{
			return getOperatorQubitsNumber(U);
		};

		unsigned int getNrQubits() const
		{
			return nrQubits;
		}

	protected:
		static unsigned int getOperatorQubitsNumber(const MatrixClass& op)
		{
			unsigned int sz = static_cast<unsigned int>(op.rows());

			unsigned int res = 0;
			while (sz) {
				++res;
				sz >>= 1;
			}

			return res;
		};

		MatrixClass U;
		unsigned int nrQubits;
		unsigned int ctrlQubit;
	};

}


