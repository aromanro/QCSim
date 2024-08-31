#pragma once

#include "MPSSimulatorImpl.h"

#include <unordered_map>

namespace QC
{

	namespace TensorNetworks
	{

		// this is going to allow two qubit gates to be applied on qubits that are not adjacent
		class MPSSimulator : public MPSSimulatorInterface
		{
		public:
			MPSSimulator() = delete;

			MPSSimulator(size_t N, int addseed = 0)
				: impl(N, addseed)
			{
			}

			// TODO: implement it!

		private:
			MPSSimulatorImpl impl;
			std::unordered_map<IndexType, IndexType> qubitsMap;
			QC::Gates::SwapGate<MatrixClass> swapGate;
		};

	}

}