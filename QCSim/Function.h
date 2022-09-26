#pragma once

#include <Eigen/eigen>

#include "QubitRegister.h"

namespace QC {

	// a function may be implemented out of quantum gates or have it a 'black box' - just construct its matrix and apply it on the register

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class Function
	{
	public:
		virtual void Apply(QubitRegister<VectorClass, MatrixClass>& reg) = 0;
	};

}


