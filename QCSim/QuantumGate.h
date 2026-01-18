#pragma once

#include "SimpleGates.h"

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {
	namespace Gates {

		template<class MatrixClass = Eigen::MatrixXcd> class SwapGate : public TwoQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsGate<MatrixClass>;
			using OpClass = typename TwoQubitsGate<MatrixClass>::BaseClass;

			SwapGate()
			{
				OpClass::operatorMat(1, 1) = 0;
				OpClass::operatorMat(2, 2) = 0;
				OpClass::operatorMat(1, 2) = 1;
				OpClass::operatorMat(2, 1) = 1;
			}

			bool isSwapGate() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class iSwapGate : public TwoQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			iSwapGate()
			{
				OpClass::operatorMat(1, 1) = 0;
				OpClass::operatorMat(2, 2) = 0;
				OpClass::operatorMat(1, 2) = std::complex<double>(0, 1);
				OpClass::operatorMat(2, 1) = std::complex<double>(0, 1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class iSwapDagGate : public TwoQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			iSwapDagGate()
			{
				OpClass::operatorMat(1, 1) = 0;
				OpClass::operatorMat(2, 2) = 0;
				OpClass::operatorMat(1, 2) = std::complex<double>(0, -1);
				OpClass::operatorMat(2, 1) = std::complex<double>(0, -1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class DecrementGate : public TwoQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			DecrementGate()
			{
				OpClass::operatorMat(0, 0) = 0;
				OpClass::operatorMat(1, 1) = 0;
				OpClass::operatorMat(2, 2) = 0;
				OpClass::operatorMat(3, 3) = 0;

				OpClass::operatorMat(3, 0) = 1;

				OpClass::operatorMat(0, 1) = 1;
				OpClass::operatorMat(1, 2) = 1;
				OpClass::operatorMat(2, 3) = 1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class TwoQubitsControlledGate : public TwoQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			TwoQubitsControlledGate() {};

			bool isControlled() const override
			{
				return true;
			}

			bool isControlQubit(size_t qubit) const override
			{
				return qubit == 0;
			}

			void SetOperation(const MatrixClass& U)
			{
				assert(U.rows() == 2 && U.cols() == 2);
				OpClass::operatorMat.block(2, 2, 2, 2) = U;
			}
		};


		// also named controlled X gate
		template<class MatrixClass = Eigen::MatrixXcd> class CNOTGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			CNOTGate()
			{
				/*
				MatrixClass U = MatrixClass::Zero(2, 2);

				U(0, 1) = 1;
				U(1, 0) = 1;

				TwoQubitsControlledGate<MatrixClass>::SetOperation(U);
				*/

				// the commented code above is correct, but this one should be faster
				OpClass::operatorMat(2, 2) = 0;
				OpClass::operatorMat(3, 3) = 0;
				OpClass::operatorMat(2, 3) = 1;
				OpClass::operatorMat(3, 2) = 1;
			}

			// this reffers only to the controlled block, not the whole matrix
			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledYGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledYGate()
			{
				OpClass::operatorMat(2, 2) = 0;
				OpClass::operatorMat(3, 3) = 0;
				OpClass::operatorMat(2, 3) = std::complex(0., -1.);
				OpClass::operatorMat(3, 2) = std::complex(0., 1.);
			}

			// this reffers only to the controlled block, not the whole matrix
			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledZGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledZGate()
			{
				OpClass::operatorMat(3, 3) = -1;
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledHadamardGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledHadamardGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(2, 2) = norm;
				OpClass::operatorMat(2, 3) = norm;
				OpClass::operatorMat(3, 2) = norm;
				OpClass::operatorMat(3, 3) = -norm;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledSquareRootNOTGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledSquareRootNOTGate()
			{
				OpClass::operatorMat(2, 2) = std::complex<double>(0.5, 0.5);
				OpClass::operatorMat(2, 3) = std::complex<double>(0.5, -0.5);
				OpClass::operatorMat(3, 2) = OpClass::operatorMat(2, 3);
				OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledSquareRootNOTDagGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledSquareRootNOTDagGate()
			{
				OpClass::operatorMat(2, 2) = std::complex<double>(0.5, -0.5);
				OpClass::operatorMat(2, 3) = std::complex<double>(0.5, 0.5);
				OpClass::operatorMat(3, 2) = OpClass::operatorMat(2, 3);
				OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
			}
		};

		// controlled-V gate
		template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledPhaseGate()
			{
				OpClass::operatorMat(3, 3) = std::complex<double>(0, 1);
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseShiftGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledPhaseShiftGate(double theta = 0)
			{
				SetPhaseShift(theta);
			}

			void SetPhaseShift(double theta)
			{
				OpClass::operatorMat(3, 3) = std::polar(1., theta);
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledUGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledUGate(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
			{
				SetParams(theta, phi, lambda, gamma);
			}

			virtual void SetParams(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(2, 2) = std::polar(1., gamma) * cos(t2);
				OpClass::operatorMat(2, 3) = -std::polar(1., gamma + lambda) * sin(t2);
				OpClass::operatorMat(3, 2) = std::polar(1., gamma + phi) * sin(t2);
				OpClass::operatorMat(3, 3) = std::polar(1., gamma + phi + lambda) * cos(t2);
			}
		};


		template<class MatrixClass = Eigen::MatrixXcd> class ControlledRotationGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			virtual void SetTheta(double theta) = 0;
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledRxGate : public ControlledRotationGate<MatrixClass>
		{
		public:
			using BaseClass = ControlledRotationGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledRxGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(2, 2) = std::complex<double>(cos(t2), 0);
				OpClass::operatorMat(2, 3) = std::complex<double>(0, -sin(t2));
				OpClass::operatorMat(3, 2) = OpClass::operatorMat(2, 3);
				OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledRyGate : public ControlledRotationGate<MatrixClass>
		{
		public:
			using BaseClass = ControlledRotationGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledRyGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(2, 2) = std::complex<double>(cos(t2), 0);
				OpClass::operatorMat(2, 3) = std::complex<double>(-sin(t2), 0);
				OpClass::operatorMat(3, 2) = std::complex<double>(sin(t2), 0);
				OpClass::operatorMat(3, 3) = OpClass::operatorMat(2, 2);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledRzGate : public ControlledRotationGate<MatrixClass>
		{
		public:
			using BaseClass = ControlledRotationGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ControlledRzGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(2, 2) = std::polar(1., -t2);
				OpClass::operatorMat(3, 3) = std::polar(1., t2);
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};



		template<class MatrixClass = Eigen::MatrixXcd> class ThreeQubitsControlledGate : public ThreeQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			ThreeQubitsControlledGate() {};

			bool isControlled() const override
			{
				return true;
			}

			bool isControlQubit(size_t qubit) const override
			{
				return qubit == 0;
			}

			void SetOperation(const MatrixClass& U)
			{
				assert(U.rows() == 4 && U.cols() == 4);
				OpClass::operatorMat.block(4, 4, 4, 4) = U;
			}
		};

		// also known as CCNOT or CCX
		template<class MatrixClass = Eigen::MatrixXcd> class ToffoliGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			ToffoliGate()
			{
				OpClass::operatorMat(6, 6) = 0;
				OpClass::operatorMat(7, 7) = 0;
				OpClass::operatorMat(6, 7) = 1;
				OpClass::operatorMat(7, 6) = 1;
			}

			bool isControlQubit(size_t qubit) const override
			{
				return qubit <= 1;
			}

			// this reffers only to the controlled block, not the whole matrix
			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		// also known as CSWAP
		template<class MatrixClass = Eigen::MatrixXcd> class FredkinGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			FredkinGate()
			{
				OpClass::operatorMat(5, 5) = 0;
				OpClass::operatorMat(6, 6) = 0;
				OpClass::operatorMat(5, 6) = 1;
				OpClass::operatorMat(6, 5) = 1;
			}

			bool isSwapGate() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class CCZGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			CCZGate()
			{
				OpClass::operatorMat(7, 7) = -1;
			}

			bool isControlQubit(size_t qubit) const override
			{
				return qubit <= 1;
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

	}
}


