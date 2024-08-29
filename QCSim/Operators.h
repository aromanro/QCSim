#pragma once

#include "QuantumGate.h"


namespace QC {
	namespace Operators {

		template<class MatrixClass = Eigen::MatrixXcd> class SigmaPlus : public Gates::SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = Gates::SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SigmaPlus()
			{
				OpClass::operatorMat(1, 0) = 1;
			}

			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SigmaMinus : public Gates::SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = Gates::SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SigmaMinus()
			{
				OpClass::operatorMat(0, 1) = 1;
			}

			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ZeroProjection : public Gates::SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = Gates::SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			ZeroProjection()
			{
				OpClass::operatorMat(0, 0) = 1;
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class OneProjection : public Gates::SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = Gates::SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			OneProjection()
			{
				OpClass::operatorMat(1, 1) = 1;
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

	}
}

