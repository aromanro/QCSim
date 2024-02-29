#pragma once

#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NQubitsQuantumGate : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			NQubitsQuantumGate(size_t N, const MatrixClass& op, size_t startQubit = 0)
				: BaseClass(N, startQubit, startQubit + getOperatorQubitsNumber(op)), U(op), nrQubits(N)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				reg.ApplyOperatorMatrix(getOperatorMatrix());

				return 0; // no measurement, so the return should be ignored
			}

			MatrixClass getOperatorMatrix() const
			{
				const size_t startQubit = BaseClass::getStartQubit();
				const size_t endQubit = BaseClass::getEndQubit();
				const size_t nrBasisStates = 1u << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				size_t mask = 0;
				for (size_t i = startQubit; i <= endQubit; ++i)
					mask |= 1u << i;

				for (size_t i = 0; i < nrBasisStates; ++i)
				{
					const size_t ind1 = i | mask;
					for (size_t j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
							extOperatorMat(i, j) = U((i & mask) >> startQubit, (j & mask) >> startQubit);
				}

				return extOperatorMat;
			}

			const MatrixClass& getRawOperatorMatrix() const
			{
				return U;
			}

			size_t getRawOperatorQubitsNumber() const
			{
				return getOperatorQubitsNumber(U);
			};

		protected:
			static size_t getOperatorQubitsNumber(const MatrixClass& op)
			{
				size_t sz = static_cast<size_t>(op.rows() - 1);

				size_t res = 0;
				while (sz) {
					++res;
					sz >>= 1;
				}

				return res;
			};

			MatrixClass U;
			size_t nrQubits;
		};

	}

}


