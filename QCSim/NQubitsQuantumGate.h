#pragma once

#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NQubitsQuantumGate : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			NQubitsQuantumGate(unsigned int N, const MatrixClass& op, unsigned int startQubit = 0)
				: BaseClass(N, startQubit, startQubit + getOperatorQubitsNumber(op)), U(op), nrQubits(N)
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
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				unsigned int mask = 0;
				for (unsigned int i = startQubit; i <= endQubit; ++i)
					mask |= 1u << i;

				for (unsigned int i = 0; i < nrBasisStates; ++i)
				{
					const unsigned int ind1 = i | mask;
					for (unsigned int j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
							extOperatorMat(i, j) = U((i & mask) >> startQubit, (j & mask) >> startQubit);
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

		protected:
			static unsigned int getOperatorQubitsNumber(const MatrixClass& op)
			{
				unsigned int sz = static_cast<unsigned int>(op.rows() - 1);

				unsigned int res = 0;
				while (sz) {
					++res;
					sz >>= 1;
				}

				return res;
			};

			MatrixClass U;
			unsigned int nrQubits;
		};

	}

}


