#pragma once

#include "QubitRegister.h"

namespace QC {
	
	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumAlgorithm
	{
	public:
		QuantumAlgorithm(unsigned int N = 3, int addseed = 0)
			: reg(N, addseed)
		{
		}

		virtual unsigned int Execute() = 0;

		static unsigned int getQubitState(unsigned int q)
		{
			return 1u << q;
		}

		const QubitRegister<VectorClass, MatrixClass>& getRegister() const
		{
			return reg;
		}

	protected:
		QubitRegister<VectorClass, MatrixClass> reg;
	};

}
