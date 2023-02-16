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
		};

		template<class MatrixClass = Eigen::MatrixXcd> class QuantumGateWithOp : public QuantumGate<MatrixClass>
		{
		public:
			typedef QuantumGate<MatrixClass> BaseClass;

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

		protected:
			MatrixClass operatorMat;
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SingleQubitGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			typedef QuantumGateWithOp<MatrixClass> BaseClass;

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

			// controllingQubit is ignored, it will be used for two qubit gates
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
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

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
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			HyGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -norm);
				OpClass::operatorMat(1, 0) = std::complex<double>(0, norm);
				OpClass::operatorMat(1, 1) = -norm;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PhaseGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			PhaseGate(double theta = 0)
			{
				OpClass::operatorMat(0, 0) = 1;
				OpClass::operatorMat(1, 1) = std::complex<double>(0, 1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PhaseShiftGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

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

		template<class MatrixClass = Eigen::MatrixXcd> class PauliXGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			PauliXGate()
			{
				OpClass::operatorMat(0, 1) = 1;
				OpClass::operatorMat(1, 0) = 1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PauliYGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			PauliYGate()
			{
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -1);
				OpClass::operatorMat(1, 0) = std::complex<double>(0, 1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class PauliZGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			PauliZGate()
			{
				OpClass::operatorMat(0, 0) = 1;
				OpClass::operatorMat(1, 1) = -1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class NOTSquareRootGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			NOTSquareRootGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = -norm;
				OpClass::operatorMat(1, 0) = norm;
				OpClass::operatorMat(1, 1) = norm;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class SplitterGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			SplitterGate()
			{
				static const double norm = 1. / sqrt(2.);
				OpClass::operatorMat(0, 0) = norm;
				OpClass::operatorMat(0, 1) = std::complex<double>(0, norm);
				OpClass::operatorMat(1, 0) = std::complex<double>(0, norm);
				OpClass::operatorMat(1, 1) = norm;
			}
		};

		// rotation gates

		template<class MatrixClass = Eigen::MatrixXcd> class RotationGate : public SingleQubitGate<MatrixClass>
		{
		public:
			typedef SingleQubitGate<MatrixClass> BaseClass;
			
			virtual void SetTheta(double theta) = 0;
		};

		template<class MatrixClass = Eigen::MatrixXcd> class RxGate : public RotationGate<MatrixClass>
		{
		public:
			typedef RotationGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

			RxGate(double theta = 0)
			{
				SetTheta(theta);
			}

			void SetTheta(double theta) override
			{
				const double t2 = theta * 0.5;

				OpClass::operatorMat(0, 0) = std::complex<double>(cos(t2), 0);
				OpClass::operatorMat(0, 1) = std::complex<double>(0, -sin(t2));
				OpClass::operatorMat(1, 0) = std::complex<double>(0, -sin(t2));
				OpClass::operatorMat(1, 1) = std::complex<double>(cos(t2), 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class RyGate : public RotationGate<MatrixClass>
		{
		public:
			typedef RotationGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

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
				OpClass::operatorMat(1, 1) = std::complex<double>(cos(t2), 0);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class RzGate : public RotationGate<MatrixClass>
		{
		public:
			typedef RotationGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

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
			typedef RotationGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

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
			typedef QuantumGateWithOp<MatrixClass> BaseClass;

			TwoQubitsGate()
				: BaseClass (MatrixClass::Identity(4, 4))
			{
			}

			unsigned int getQubitsNumber() const override
			{
				return 2;
			}

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
			typedef TwoQubitsGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

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
			typedef TwoQubitsGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

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
			typedef TwoQubitsGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

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
			typedef TwoQubitsGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			TwoQubitsControlledGate() {};

			void SetOperation(const MatrixClass& U)
			{
				assert(U.rows() == 2 && U.cols() == 2);
				OpClass::operatorMat.block(2, 2, 2, 2) = U;
			}
		};



		template<class MatrixClass = Eigen::MatrixXcd> class CNOTGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			typedef TwoQubitsControlledGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

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

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			typedef TwoQubitsControlledGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

			ControlledPhaseGate()
			{
				OpClass::operatorMat(3, 3) = std::complex<double>(0, 1);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledPhaseShiftGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			typedef TwoQubitsControlledGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

			ControlledPhaseShiftGate(double theta = 0)
			{
				SetPhaseShift(theta);
			}

			void SetPhaseShift(double theta)
			{
				OpClass::operatorMat(3, 3) = std::polar(1., theta);
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ControlledZGate : public TwoQubitsControlledGate<MatrixClass>
		{
		public:
			typedef TwoQubitsControlledGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

			ControlledZGate()
			{
				OpClass::operatorMat(3, 3) = -1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ThreeQubitsGate : public QuantumGateWithOp<MatrixClass>
		{
		public:
			typedef QuantumGateWithOp<MatrixClass> BaseClass;

			ThreeQubitsGate()
				: BaseClass(MatrixClass::Identity(8, 8))
			{
			}

			unsigned int getQubitsNumber() const override
			{
				return 3;
			}

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
			typedef ThreeQubitsGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass OpClass;

			ThreeQubitsControlledGate() {};

			void SetOperation(const MatrixClass& U)
			{
				assert(U.rows() == 4 && U.cols() == 4);
				OpClass::operatorMat.block(4, 4, 4, 4) = U;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class ToffoliGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			typedef ThreeQubitsControlledGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

			ToffoliGate()
			{
				OpClass::operatorMat(6, 6) = 0;
				OpClass::operatorMat(7, 7) = 0;
				OpClass::operatorMat(6, 7) = 1;
				OpClass::operatorMat(7, 6) = 1;
			}
		};

		template<class MatrixClass = Eigen::MatrixXcd> class FredkinGate : public ThreeQubitsControlledGate<MatrixClass>
		{
		public:
			typedef ThreeQubitsControlledGate<MatrixClass> BaseClass;
			typedef BaseClass::BaseClass::BaseClass OpClass;

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


