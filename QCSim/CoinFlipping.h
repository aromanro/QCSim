#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

#include "Tests.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Games {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class CoinFlipping : 
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		CoinFlipping(int addseed = 0)
			: BaseClass(1, addseed)
		{
			BaseClass::setToBasisState(0);

			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);
		}

		size_t Execute() override
		{
			BaseClass::setToBasisState(0);

			BaseClass::ApplyGate(hadamardGate, 0);

			if (dist_bool(rng))
				BaseClass::ApplyGate(flipGate, 0);

			BaseClass::ApplyGate(hadamardGate, 0);

			return BaseClass::Measure();
		}

	protected:
		QC::Gates::PauliXGate<MatrixClass> flipGate;
		QC::Gates::HadamardGate<MatrixClass> hadamardGate;

		std::mt19937_64 rng;
		std::bernoulli_distribution dist_bool;
	};

}

