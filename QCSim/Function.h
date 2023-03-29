#pragma once

#include <Eigen/eigen>

#include "QubitRegister.h"

namespace QC {

	// a function may be implemented out of quantum gates or have it a 'black box' - just construct its matrix and apply it on the register

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class Function
	{
	public:
		using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

		virtual ~Function() {}

		virtual void Apply(RegisterClass& reg) = 0;
	};

}


