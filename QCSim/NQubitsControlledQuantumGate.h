#pragma once


#include "QuantumAlgorithm.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class NQubitsControlledQuantumGate : public QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>
		{
		public:
			using BaseClass = QuantumSubAlgorithmOnSubregister<VectorClass, MatrixClass>;
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			NQubitsControlledQuantumGate(size_t N, const MatrixClass& op, size_t startQubit = 1, size_t controllingQubit = 0)
				: BaseClass(N, startQubit, startQubit + getOperatorQubitsNumber(op)), U(op), nrQubits(N), ctrlQubit(controllingQubit)
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
				const size_t nrBasisStates = 1ULL << nrQubits;
				const size_t ctrlQubitBit = 1ULL << ctrlQubit;

				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);
				MatrixClass Unity = MatrixClass::Identity(U.rows(), U.cols());

				size_t mask = 0;
				for (size_t i = startQubit; i <= endQubit; ++i)
					mask |= 1ULL << i;

				for (size_t i = 0; i < nrBasisStates; ++i)
				{
					const size_t ind1 = i | mask;
					for (size_t j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
						{
							const size_t ind1op = (i & mask) >> startQubit;
							const size_t ind2op = (j & mask) >> startQubit;

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

			size_t getRawOperatorQubitsNumber() const
			{
				return getOperatorQubitsNumber(U);
			};

			size_t getNrQubits() const
			{
				return nrQubits;
			}

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
			size_t ctrlQubit;
		};

	}

}


