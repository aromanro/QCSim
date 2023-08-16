#pragma once

#include "QuantumAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace MachineLearning {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QMeansClustering2D :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		struct DataPoint {
			double x = 0.;
			double y = 0.;
			int cluster = -1;
		};

		QMeansClustering2D(int addseed = 0)
		: BaseClass(3, addseed)
		{
		}

		// TODO: implement it
	};

}


