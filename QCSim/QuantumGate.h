#pragma once

#include <Eigen/eigen>

// Qubits are numbered from right to left, starting with zero, this might be confusing, since notation numbers them usually from left to right

namespace QC {
	namespace Gates {

		template<class MatrixClass = Eigen::MatrixXcd> class QuantumGate
		{
		public:
			// controllingQubit1 is for two qubit gates and controllingQubit2 is for three qubit gates, they are ignored for gates with a lower number of qubits
			virtual MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const = 0;
			virtual ~QuantumGate() {};

			virtual unsigned int getQubitsNumber() const
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

			virtual bool isControlQubit(unsigned int qubit) const
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
			};

			const MatrixClass& getRawOperatorMatrix() const
			{
				return operatorMat;
			};

			unsigned int getQubitsNumber() const override
			{
				unsigned int sz = static_cast<unsigned int>(operatorMat.rows());

				unsigned int res = 0;
				while (sz) {
					++res;
					sz >>= 1;
				}

				return res;
			};

			void setOperator(const MatrixClass& U)
			{
				assert(U.rows() == U.cols());

				operatorMat = U;
			}

		protected:
			MatrixClass operatorMat;
		};

		// used for recording the applied gates, to be applied again later
		// also for 'uncompute'
		template<class MatrixClass = Eigen::MatrixXcd> class AppliedGate : public Gates::QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = Gates::QuantumGateWithOp<MatrixClass>;

			AppliedGate(const MatrixClass& op, unsigned int q1 = 0, unsigned int q2 = 0, unsigned int q3 = 0)
				: Gates::QuantumGateWithOp<MatrixClass>(op), q1(q1), q2(q2), q3(q3)
			{
			}

			unsigned int getQubit1() const { return q1; }
			unsigned int getQubit2() const { return q2; }
			unsigned int getQubit3() const { return q3; }

		private:
			// don't use it!
			MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const override
			{
				return BaseClass::getRawOperatorMatrix();
			}

			unsigned int q1;
			unsigned int q2;
			unsigned int q3;
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SingleQubitGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = QuantumGateWithOp<MatrixClass>;

			SingleQubitGate()
				: BaseClass(MatrixClass::Zero(2, 2))
			{
			};

			SingleQubitGate(const MatrixClass& U)
				: BaseClass(U)
			{
			};

			unsigned int getQubitsNumber() const override
			{
				return 1;
			}

			// controllingQubit1 and controllingQubit2 are ignored, they will be used for two (only controllingQubit1) and three qubit gates
			// this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
			// even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for details)
			// but I won't bother, I'll keep this different implementation because it might be more clear - construction by tensor product
			// and also could be useful for debugging in the case the optimized version has a problem
			MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const override
			{
				const unsigned int nrBasisStates = 1u << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				const unsigned int qubitBit = 1u << qubit;

				// since this can be quite condensed, here is a description:
				// this just computes the tensor product between the 2x2 operator matrix for the qubit and the identity operators for the other qubits

				// for the two and three qubits operators it's more complex but analogous

				for (unsigned int i = 0; i < nrBasisStates; ++i)
				{
					const unsigned int ind1 = i | qubitBit;
					for (unsigned int j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | qubitBit)) // this is just a lot of delta 'functions' for the 'other' qubits, the bit was forced to 1 for the qubit with this operator
							extOperatorMat(i, j) = BaseClass::operatorMat(i & qubitBit ? 1 : 0, j & qubitBit ? 1 : 0); // pick the correct matrix element for the qubit operator 
				}

				return extOperatorMat;
			};
		};

		template<class MatrixClass = Eigen::MatrixXcd> class HadamardGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			HadamardGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = norm;
				OpClass::operatorMat(1, 0) = norm;
				OpClass::operatorMat(1, 1) = -norm;
			}
		};

		// hadamard can be used to switch back and forth to X basis
		// this gate is for switching to Y basis and back
		template<class MatrixClass = Eigen::MatrixXcd> class HyGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			HyGate()
			{
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
			using OpClass = BaseClass::BaseClass;

			PhaseGate()
			{
				OpClass::operatorMat(0, 0) = 1;
				OpClass::operatorMat(1, 1) = std::complex<double>(0, 1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PhaseShiftGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			PhaseShiftGate(double theta = 0)
			{
				OpClass::operatorMat(0, 0) = 1;
				SetPhaseShift(theta);
			}

			void SetPhaseShift(double theta)
			{
				OpClass::operatorMat(1, 1) = std::polar(1., theta);
			}
		};

		// also known as the flip gate or not gate
		template<class MatrixClass = Eigen::MatrixXcd> class PauliXGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			PauliXGate()
			{
				OpClass::operatorMat(0, 1) = 1;
				OpClass::operatorMat(1, 0) = 1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PauliYGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			PauliYGate()
			{
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -1);
				OpClass::operatorMat(1, 0) = std::complex<double>(0, 1);
			}
		};

		// also known as the phase flip gate
		template<class MatrixClass = Eigen::MatrixXcd> class PauliZGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			PauliZGate()
			{
				OpClass::operatorMat(0, 0) = 1;
				OpClass::operatorMat(1, 1) = -1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SigmaPlusGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			SigmaPlusGate()
			{
				OpClass::operatorMat(0, 1) = 1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SigmaMinusGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			SigmaMinusGate()
			{
				OpClass::operatorMat(1, 0) = 1;
			}
		};


		template<class MatrixClass = Eigen::MatrixXcd> class SquareRootNOTGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			SquareRootNOTGate()
			{	
				OpClass::operatorMat(0, 0) = std::complex<double>(0.5, 0.5);
				OpClass::operatorMat(0, 1) = std::complex<double>(0.5, -0.5);
				OpClass::operatorMat(1, 0) = OpClass::operatorMat(0, 1);
				OpClass::operatorMat(1, 1) = OpClass::operatorMat(0, 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SplitterGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass::BaseClass;

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
		};

		// generic rotation gate
		template<class MatrixClass = Eigen::MatrixXcd> class UGate : public SingleQubitGate<MatrixClass>
		{
		public:
			using BaseClass = SingleQubitGate<MatrixClass>;
			using OpClass = SingleQubitGate<MatrixClass>::BaseClass;

			UGate(double theta = 0, double phi = 0, double lambda = 0)
			{
				SetParams(theta, phi, lambda);
			}

			virtual void SetParams(double theta = 0, double phi = 0, double lambda = 0)
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(0, 0) = std::complex<double>(cos(t2), 0.);
				OpClass::operatorMat(0, 1) = -std::polar(1., lambda) * sin(t2);
				OpClass::operatorMat(1, 0) = std::polar(1., phi) * sin(t2);
				OpClass::operatorMat(1, 1) = std::polar(1., phi + lambda) * cos(t2);
			}
		};


		template<class MatrixClass = Eigen::MatrixXcd> class TwoQubitsGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			using BaseClass = QuantumGateWithOp<MatrixClass>;

			TwoQubitsGate()
				: BaseClass (MatrixClass::Identity(4, 4))
			{
			}

			unsigned int getQubitsNumber() const override
			{
				return 2;
			}

			// controllingQubit2 is ignored, it will be used for three qubit gates only
			// this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
			// even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for details)
			// but I won't bother, I'll keep this different implementation because it might be more clear - construction by tensor product
			// and also could be useful for debugging in the case the optimized version has a problem
			MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const override
			{
				assert(qubit != controllingQubit1);

				const unsigned int nrBasisStates = 1u << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				const unsigned int qubitBit = 1u << qubit;
				const unsigned int ctrlQubitBit = 1u << controllingQubit1;
				const unsigned int mask = qubitBit | ctrlQubitBit;

				// computing the tensor product between the gate matrix and identity operators for the other qubits

				for (unsigned int i = 0; i < nrBasisStates; ++i)
				{
					const unsigned int ind1 = i | mask;
					for (unsigned int j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
							extOperatorMat(i, j) = BaseClass::operatorMat((i & ctrlQubitBit ? 2 : 0) | (i & qubitBit ? 1 : 0), (j & ctrlQubitBit ? 2 : 0) | (j & qubitBit ? 1 : 0));
				}

				return extOperatorMat;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SwapGate : public TwoQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsGate<MatrixClass>;
			using OpClass = TwoQubitsGate<MatrixClass>::BaseClass;

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
			using OpClass = BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass;

			TwoQubitsControlledGate() {};

			bool isControlled() const override
			{
				return true;
			}

			bool isControlQubit(unsigned int qubit) const override
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
			using OpClass = BaseClass::BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass::BaseClass;

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
			using OpClass = BaseClass::BaseClass::BaseClass;

			ControlledZGate()
			{
				OpClass::operatorMat(3, 3) = -1;
			}
		};


		template<class MatrixClass = Eigen::MatrixXcd> class ControlledHadamardGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass::BaseClass;

			ControlledHadamardGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(2, 2) = norm;
				OpClass::operatorMat(2, 3) = norm;
				OpClass::operatorMat(3, 2) = norm;
				OpClass::operatorMat(3, 3) = -norm;
			}
		};

		// controlled-V gate
		template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass::BaseClass;

			ControlledPhaseGate()
			{
				OpClass::operatorMat(3, 3) = std::complex<double>(0, 1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseShiftGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = TwoQubitsControlledGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass::BaseClass;

			ControlledPhaseShiftGate(double theta = 0)
			{
				SetPhaseShift(theta);
			}

			void SetPhaseShift(double theta)
			{
				OpClass::operatorMat(3, 3) = std::polar(1., theta);
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

			unsigned int getQubitsNumber() const override
			{
				return 3;
			}

			// this is not used anymore, instead there is a more efficient implementation in QubitRegister::applyGate
			// even this matrix could be constructed more efficiently in a similar manner (see QubitRegister::applyGate for details)
			// but I won't bother, I'll keep this different implementation because it might be more clear - construction by tensor product
			// and also could be useful for debugging in the case the optimized version has a problem
			MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit1 = 0, unsigned int controllingQubit2 = 0) const override
			{
				assert(qubit != controllingQubit1 && controllingQubit1 != controllingQubit2);

				const unsigned int nrBasisStates = 1u << nrQubits;
				MatrixClass extOperatorMat = MatrixClass::Zero(nrBasisStates, nrBasisStates);

				const unsigned int qubitBit = 1u << qubit;
				const unsigned int qubitBit2 = 1u << controllingQubit1;
				const unsigned int ctrlQubitBit = 1u << controllingQubit2;
				const unsigned int mask = qubitBit | qubitBit2 | ctrlQubitBit;

				// computing the tensor product between the gate matrix and identity operators for the other qubits

				for (unsigned int i = 0; i < nrBasisStates; ++i)
				{
					const unsigned int ind1 = i | mask;
					for (unsigned int j = 0; j < nrBasisStates; ++j)
						if (ind1 == (j | mask)) // the delta 'function'
							extOperatorMat(i, j) = BaseClass::operatorMat((i & ctrlQubitBit ? 4 : 0) | (i & qubitBit2 ? 2 : 0) | (i & qubitBit ? 1 : 0), (j & ctrlQubitBit ? 4 : 0) | (j & qubitBit2 ? 2 : 0) | (j & qubitBit ? 1 : 0));
				}

				return extOperatorMat;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ThreeQubitsControlledGate : public ThreeQubitsGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass;

			ThreeQubitsControlledGate() {};

			bool isControlled() const override
			{
				return true;
			}

			bool isControlQubit(unsigned int qubit) const override
			{
				return qubit == 0;
			}

			void SetOperation(const MatrixClass& U)
			{
				assert(U.rows() == 4 && U.cols() == 4);
				OpClass::operatorMat.block(4, 4, 4, 4) = U;
			}
		};

		// also known as CCNOT
		template<class MatrixClass = Eigen::MatrixXcd> class ToffoliGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass::BaseClass;

			ToffoliGate()
			{
				OpClass::operatorMat(6, 6) = 0;
				OpClass::operatorMat(7, 7) = 0;
				OpClass::operatorMat(6, 7) = 1;
				OpClass::operatorMat(7, 6) = 1;
			}

			bool isControlQubit(unsigned int qubit) const override
			{
				return qubit <= 1;
			}
		};

		// also known as CSWAP
		template<class MatrixClass = Eigen::MatrixXcd> class FredkinGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			using BaseClass = ThreeQubitsControlledGate<MatrixClass>;
			using OpClass = BaseClass::BaseClass::BaseClass;

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


