#pragma once

#include "QubitRegister.h"

namespace QC {
	
	class QuantumAlgorithm
	{
	public:
		QuantumAlgorithm(int N = 3, int addseed = 0);

		virtual unsigned int Execute() = 0;

		QubitRegister reg;
	};

}
