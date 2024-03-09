#pragma once

#include <Eigen/Eigen>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {
	namespace Gates {

		template<class MatrixClass = Eigen::MatrixXcd> class QuantumGate
		{
		public:
			// controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for gates with a lower number of qubits
			virtual MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const = 0;
			virtual ~QuantumGate() = default;

			virtual size_t getQubitsNumber() const
			{
				return 0; // 0 at this point means 'unknown', should be overriden if needed
			}

			// the following two functions could be used for some optimizatiions
			// only for reducing a constant... if it matters at one point I'll implement them
			// see QubitRegister::applyGate for details of the optimizations
			virtual bool isControlled() const
			{
				return false;
			}

			virtual bool isControlQubit(size_t qubit) const
			{
				return false;
			}

			virtual bool isDiagonal() const
			{
				return false;
			}

			virtual bool isAntidiagonal() const
			{
				return false;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class QuantumGateWithOp : public QuantumGate<MatrixClass>
		{
		public:
			using BaseClass = QuantumGate<MatrixClass>;

			QuantumGateWithOp(const MatrixClass& U)
				: operatorMat(U)
			{
			}

			const MatrixClass& getRawOperatorMatrix() const
			{
				return operatorMat;
			}

			size_t getQubitsNumber() const override
			{
				size_t sz = static_cast<size_t>(operatorMat.rows() - 1);

				size_t res = 0;
				while (sz) {
					++res;
					sz >>= 1;
				}

				return res;
			}

			void setOperator(const MatrixClass& U)
			{
				assert(U.rows() == U.cols());
				assert(U.rows() == 1ULL << getQubitsNumber());

				operatorMat = U;
			}

		protected:
			MatrixClass operatorMat;
		};


		template<class MatrixClass = Eigen::MatrixXcd> class SingleQubitGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = QuantumGateWithOp<MatrixClass>;

			SingleQubitGate()
				: BaseClass(MatrixClass::Zero(2, 2))
			{
			}

			SingleQubitGate(const MatrixClass& U)
				: BaseClass(U)
			{
				assert(U.rows() == U.cols());
				assert(U.rows() == 2);
			}

			size_t getQubitsNumber() const override
			{
				return 1;
			}

			// controllingQubit1 and controllingQubit2 are ignored, they will be used for two (only controllingQubit1) and three qubit gates
			// this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
			// even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for details)
			// but I won't bother, I'll keep this different implementation because it might be more clear - construction by tensor product
			// and also could be useful for debugging in the case the optimized version has a problem
			MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
			{
				const size_t nrBasisStates = 1ULL << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				const size_t qubitBit = 1ULL << qubit;

				// since this can be quite condensed, here is a description:
				// this just computes the tensor product between the 2x2 operator matrix for the qubit and the identity operators for the other qubits

				// for the two and three qubits operators it's more complex but analogous

				for (size_t i = 0; i < nrBasisStates; ++i)
				{
					const size_t ind1 = i | qubitBit;
					for (size_t j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | qubitBit)) // this is just a lot of delta 'functions' for the 'other' qubits, the bit was forced to 1 for the qubit with this operator
							extOperatorMat(i, j) = BaseClass::operatorMat(i & qubitBit ? 1 : 0, j & qubitBit ? 1 : 0); // pick the correct matrix element for the qubit operator 
				}

				return extOperatorMat;
			};
		};

		template<class MatrixClass = Eigen::MatrixXcd> class TwoQubitsGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = QuantumGateWithOp<MatrixClass>;

			TwoQubitsGate()
				: BaseClass(MatrixClass::Identity(4, 4))
			{
			}

			TwoQubitsGate(const MatrixClass& U)
				: BaseClass(U)
			{
				assert(U.rows() == U.cols());
				assert(U.rows() == 4);
			}

			size_t getQubitsNumber() const override
			{
				return 2;
			}

			// controllingQubit2 is ignored, it will be used for three qubit gates only
			// this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
			// even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for details)
			// but I won't bother, I'll keep this different implementation because it might be more clear - construction by tensor product
			// and also could be useful for debugging in the case the optimized version has a problem
			MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
			{
				assert(qubit != controllingQubit1);

				const size_t nrBasisStates = 1ULL << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				const size_t qubitBit = 1ULL << qubit;
				const size_t ctrlQubitBit = 1ULL << controllingQubit1;
				const size_t mask = qubitBit | ctrlQubitBit;

				// computing the tensor product between the gate matrix and identity operators for the other qubits

				for (size_t i = 0; i < nrBasisStates; ++i)
				{
					const size_t ind1 = i | mask;
					for (size_t j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
							extOperatorMat(i, j) = BaseClass::operatorMat((i & ctrlQubitBit ? 2 : 0) | (i & qubitBit ? 1 : 0), (j & ctrlQubitBit ? 2 : 0) | (j & qubitBit ? 1 : 0));
				}

				return extOperatorMat;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ThreeQubitsGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = QuantumGateWithOp<MatrixClass>;

			ThreeQubitsGate()
				: BaseClass(MatrixClass::Identity(8, 8))
			{
			}

			ThreeQubitsGate(const MatrixClass& U)
				: BaseClass(U)
			{
				assert(U.rows() == U.cols());
				assert(U.rows() == 8);
			}

			size_t getQubitsNumber() const override
			{
				return 3;
			}

			// this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
			// even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for details)
			// but I won't bother, I'll keep this different implementation because it might be more clear - construction by tensor product
			// and also could be useful for debugging in the case the optimized version has a problem
			MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
			{
				assert(qubit != controllingQubit1 && controllingQubit1 != controllingQubit2);

				const size_t nrBasisStates = 1ULL << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				const size_t qubitBit = 1ULL << qubit;
				const size_t qubitBit2 = 1ULL << controllingQubit1;
				const size_t ctrlQubitBit = 1ULL << controllingQubit2;
				const size_t mask = qubitBit | qubitBit2 | ctrlQubitBit;

				// computing the tensor product between the gate matrix and identity operators for the other qubits

				for (size_t i = 0; i < nrBasisStates; ++i)
				{
					const size_t ind1 = i | mask;
					for (size_t j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
							extOperatorMat(i, j) = BaseClass::operatorMat((i & ctrlQubitBit ? 4 : 0) | (i & qubitBit2 ? 2 : 0) | (i & qubitBit ? 1 : 0), (j & ctrlQubitBit ? 4 : 0) | (j & qubitBit2 ? 2 : 0) | (j & qubitBit ? 1 : 0));
				}

				return extOperatorMat;
			}
		};

		// used for recording the applied gates, to be applied again later
		// also for 'uncompute'
		template<class MatrixClass = Eigen::MatrixXcd> class AppliedGate : public Gates::QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = Gates::QuantumGateWithOp<MatrixClass>;

			AppliedGate()
				: Gates::QuantumGateWithOp<MatrixClass>(MatrixClass::Zero(1, 1)), q1(0), q2(0), q3(0)
			{
			}

			AppliedGate(const MatrixClass& op, size_t q1 = 0, size_t q2 = 0, size_t q3 = 0)
				: Gates::QuantumGateWithOp<MatrixClass>(op), q1(q1), q2(q2), q3(q3)
			{
			}

			size_t getQubit1() const { return q1; }
			size_t getQubit2() const { return q2; }
			size_t getQubit3() const { return q3; }

			void setQubit1(size_t q) { q1 = q; }
			void setQubit2(size_t q) { q2 = q; }
			void setQubit3(size_t q) { q3 = q; }

		private:
			// don't use it!
			MatrixClass getOperatorMatrix(size_t nrQubits, size_t qubit = 0, size_t controllingQubit1 = 0, size_t controllingQubit2 = 0) const override
			{
				// this is a hack used only for trying out the old way... without OPTIMIZED_TENSOR_PRODUCT, to not be used otherwise!
				//const size_t N = log2(BaseClass::operatorMat.rows());
				const size_t N = BaseClass::getQubitsNumber(); // no need to use log2, there is a function in the base class that computes it
				if (N == 1)
				{
					SingleQubitGate<MatrixClass> gate(BaseClass::operatorMat);
					return gate.getOperatorMatrix(nrQubits, qubit);
				}
				else if (N == 2)
				{
					TwoQubitsGate<MatrixClass> gate(BaseClass::operatorMat);
					return gate.getOperatorMatrix(nrQubits, qubit, controllingQubit1);
				}
				else if (N == 3)
				{
					ThreeQubitsGate<MatrixClass> gate(BaseClass::operatorMat);
					return gate.getOperatorMatrix(nrQubits, qubit, controllingQubit1, controllingQubit2);
				}

				return BaseClass::getRawOperatorMatrix();
			}

			size_t q1;
			size_t q2;
			size_t q3;
		};



		template<class MatrixClass = Eigen::MatrixXcd> class HadamardGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			HadamardGate()
			{
				// 1./sqrt(2) * (X + Z)
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = norm;
				OpClass::operatorMat(1, 0) = norm;
				OpClass::operatorMat(1, 1) = -norm;
			}
		};

		// hadamard can be used to switch back and forth to X basis
		// this gate is for switching to Y basis and back
		// also known as K gate
		template<class MatrixClass = Eigen::MatrixXcd> class HyGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			HyGate()
			{
				// 1./sqrt(2) * (Y + Z)
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -norm);
				OpClass::operatorMat(1, 0) = std::complex<double>(0, norm);
				OpClass::operatorMat(1, 1) = -norm;
			}
		};

		// also known as the V gate
		template<class MatrixClass = Eigen::MatrixXcd> class PhaseGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			PhaseGate()
			{
				OpClass::operatorMat(0, 0) = 1;
				OpClass::operatorMat(1, 1) = std::complex<double>(0, 1);
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PhaseShiftGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			PhaseShiftGate(double theta = 0)
			{
				OpClass::operatorMat(0, 0) = 1;
				SetPhaseShift(theta);
			}

			void SetPhaseShift(double theta)
			{
				OpClass::operatorMat(1, 1) = std::polar(1., theta);
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		// also known as the flip gate or not gate
		template<class MatrixClass = Eigen::MatrixXcd> class PauliXGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			PauliXGate()
			{
				OpClass::operatorMat(0, 1) = 1;
				OpClass::operatorMat(1, 0) = 1;
			}

			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PauliYGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			PauliYGate()
			{
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -1);
				OpClass::operatorMat(1, 0) = std::complex<double>(0, 1);
			}

			bool isAntidiagonal() const override
			{
				return true;
			}
		};

		// also known as the phase flip gate
		template<class MatrixClass = Eigen::MatrixXcd> class PauliZGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			PauliZGate()
			{
				OpClass::operatorMat(0, 0) = 1;
				OpClass::operatorMat(1, 1) = -1;
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SigmaPlusGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SigmaPlusGate()
			{
				OpClass::operatorMat(0, 1) = 1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SigmaMinusGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SigmaMinusGate()
			{
				OpClass::operatorMat(1, 0) = 1;
			}
		};

		// sx gate, also known as the square root of x gate
		template<class MatrixClass = Eigen::MatrixXcd> class SquareRootNOTGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SquareRootNOTGate()
			{	
				OpClass::operatorMat(0, 0) = std::complex<double>(0.5, 0.5);
				OpClass::operatorMat(0, 1) = std::complex<double>(0.5, -0.5);
				OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
				OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SquareRootNOTDagGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SquareRootNOTDagGate()
			{
				OpClass::operatorMat(0, 0) = std::complex<double>(0.5, -0.5);
				OpClass::operatorMat(0, 1) = std::complex<double>(0.5, 0.5);
				OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
				OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SplitterGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass;

			SplitterGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = std::complex<double>(0, norm);
				OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
				OpClass::operatorMat(1, 1) = norm;
			}
		};

		// rotation gates

		template<class MatrixClass = Eigen::MatrixXcd> class RotationGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			
			virtual void SetTheta(double theta) = 0;
		};

		template<class MatrixClass = Eigen::MatrixXcd> class RxGate : public RotationGate<MatrixClass>
		{
		public:
			using BaseClass = RotationGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			RxGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(0, 0) = std::complex<double>(cos(t2), 0);
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -sin(t2));
				OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
				OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class RyGate : public RotationGate<MatrixClass>
		{
		public:
			using BaseClass = RotationGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			RyGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(0, 0) = std::complex<double>(cos(t2), 0);
				OpClass::operatorMat(0, 1) = std::complex<double>(-sin(t2), 0);
				OpClass::operatorMat(1, 0) = std::complex<double>(sin(t2), 0);
				OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class RzGate : public RotationGate<MatrixClass>
		{
		public:
			using BaseClass = RotationGate<MatrixClass>;
			using OpClass = typename BaseClass::BaseClass::BaseClass;

			RzGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(0, 0) = std::polar(1., -t2);
				OpClass::operatorMat(1, 1) = std::polar(1., t2);
			}

			bool isDiagonal() const override
			{
				return true;
			}
		};

		// generic rotation gate
		template<class MatrixClass = Eigen::MatrixXcd> class UGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = typename SingleQubitGate<MatrixClass>::BaseClass;

			UGate(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
			{
				SetParams(theta, phi, lambda, gamma);
			}

			virtual void SetParams(double theta = 0, double phi = 0, double lambda = 0, double gamma = 0)
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(0, 0) = std::polar(1., gamma) * std::complex<double>(cos(t2), 0.);
				OpClass::operatorMat(0, 1) = -std::polar(1., gamma + lambda) * sin(t2);
				OpClass::operatorMat(1, 0) = std::polar(1., gamma + phi) * sin(t2);
				OpClass::operatorMat(1, 1) = std::polar(1., gamma + phi + lambda) * cos(t2);
			}
		};




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

				OpClass::operatorMat(2, 2) = std::polar(1., gamma) * std::complex<double>(cos(t2), 0.);
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
		};

	}
}


