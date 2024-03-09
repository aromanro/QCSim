#pragma once

#include "QuantumFourierTransform.h"
#include "Function.h"
#include "NQubitsControlledQuantumGate.h"

namespace QC {

	namespace SubAlgo {

		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PhaseEstimationBase : public QuantumSubAlgorithm<VectorClass, MatrixClass>
		{
		public:
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;

			PhaseEstimationBase(size_t N = 7, size_t L = 3)
				: fRegisterStartQubit(L), nrQubits(N), fourier(N, 0, L - 1)
			{
			}

			size_t getFunctionStartQubit() const
			{
				return fRegisterStartQubit;
			}

			size_t getNrQubits() const
			{
				return nrQubits;
			}

			double getPhase(size_t mostMeasuredState, int secondMostMeasuredState, double estimatedProbability, size_t nrSteps = 10000) const
			{
				const size_t nrStates = 1ULL << fRegisterStartQubit;
				if (secondMostMeasuredState == -1)
					return static_cast<double>(mostMeasuredState) / nrStates;

				double p = 1.;
				double phi = 0.;

				const double pref = 1. / static_cast<double>(1ull << (2 * fRegisterStartQubit));
				const std::complex<double> prefe = 2. * M_PI * std::complex(0., 1.);

				for (size_t i = 1; i < nrSteps; ++i)
				{
					const double curPhi = static_cast<double>(i) / nrSteps;
					const std::complex<double> e = prefe * curPhi;
					const std::complex<double> a = (exp(e) - 1.) / (exp(e / static_cast<double>(nrStates)) - 1.);
					const double curProb = pref * norm(a);

					if (abs(estimatedProbability - curProb) < abs(estimatedProbability - p))
					{
						p = curProb;
						phi = curPhi;
					}
				}

				if (mostMeasuredState != 0)
				{
					if (mostMeasuredState < static_cast<size_t>(secondMostMeasuredState) && static_cast<size_t>(secondMostMeasuredState) != nrStates - 1U)
						return (static_cast<double>(mostMeasuredState) + phi) / nrStates;
					else if (mostMeasuredState > static_cast<size_t>(secondMostMeasuredState))
						return (static_cast<double>(mostMeasuredState) - phi) / nrStates;
				}

				return 1. + (static_cast<double>(mostMeasuredState) - phi) / nrStates;
			}

		protected:
			void Init(RegisterClass& reg)
			{
				//QC::Gates::PauliXGate<MatrixClass> x;
				//reg.setToBasisState(0);
				//reg.ApplyGate(x, fRegisterStartQubit);

				// either the commented above, or this:
				reg.setToQubitState(fRegisterStartQubit);

				// apply hadamard over each qubit from the x-register
				ApplyHadamardOnXRegister(reg);
			}

			void ApplyHadamardOnXRegister(RegisterClass& reg) const
			{
				// apply hadamard over each qubit from the x-register
				// reuse the hadamard gate from the fourier transform base class
				for (size_t i = 0; i < fRegisterStartQubit; ++i)
					reg.ApplyGate(fourier.hadamard, i);
			}

			void IQFT(RegisterClass& reg, bool swap = true)
			{
				fourier.IQFT(reg, swap);
			}

			size_t fRegisterStartQubit;
			size_t nrQubits;

			QuantumFourierTransform<VectorClass, MatrixClass> fourier;
		};


		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class ShorPhaseEstimation : public PhaseEstimationBase<VectorClass, MatrixClass>
		{
		public:
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;
			using BaseClass = PhaseEstimationBase<VectorClass, MatrixClass>;

			ShorPhaseEstimation(QC::Function<VectorClass, MatrixClass>& f, size_t N = 7, size_t L = 3)
				: BaseClass(N, L), func(f)
			{
			}

			size_t Execute(RegisterClass& reg) override
			{
				ExecuteWithoutMeasurement(reg);

				// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
				// the f should be measured separately to find out its content

				//return reg.Measure(0, BaseClass::getFunctionStartQubit() - 1);
				return reg.MeasureAll();
			}

			std::map<size_t, size_t> ExecuteWithMultipleMeasurements(RegisterClass& reg, size_t nrMeasurements = 10000)
			{
				ExecuteWithoutMeasurement(reg);

				return reg.RepeatedMeasure(nrMeasurements);
			}

		protected:
			void ExecuteWithoutMeasurement(RegisterClass& reg)
			{
				BaseClass::Init(reg);

				// now the f(x)
				func.Apply(reg);

				// it doesn't really matter if you measure the qubits from f and when you do after the above
				// or if you measure them several times in a row
				//QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(BaseClass::getFunctionStartQubit(), QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);

				// then perform an inverse fourier transform
				BaseClass::IQFT(reg);
			}

			QC::Function<VectorClass, MatrixClass>& func;
		};



		// this one uses a big matrix made by tensor product, using the NQubitsControlledQuantumGate implementation
		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PhaseEstimation : public PhaseEstimationBase<VectorClass, MatrixClass>
		{
		public:
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;
			using BaseClass = PhaseEstimationBase<VectorClass, MatrixClass>;

			PhaseEstimation(const MatrixClass& op, size_t N = 7, size_t L = 3)
				: BaseClass(N, L), U(op)
			{
				assert(U.rows() == U.cols());
			}

			size_t Execute(RegisterClass& reg) override
			{
				ExecuteWithoutMeasurement(reg);

				// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
				// the f should be measured separately to find out its content

				//return reg.MeasureAll();
				return reg.Measure(0, BaseClass::getFunctionStartQubit() - 1);
			}

			std::map<size_t, size_t> ExecuteWithMultipleMeasurements(RegisterClass& reg, size_t nrMeasurements = 10000)
			{
				ExecuteWithoutMeasurement(reg);

				return reg.RepeatedMeasure(0, BaseClass::getFunctionStartQubit() - 1, nrMeasurements);
			}

			const MatrixClass& getU() const
			{
				return U;
			}

		protected:
			void ExecuteWithoutMeasurement(RegisterClass& reg)
			{
				// TODO: check if things are set up all right: size of U, size of reg, etc.

				BaseClass::Init(reg);

				MatrixClass controlledGate = U;
				const size_t lastQubit = BaseClass::getFunctionStartQubit() - 1;

				// There is another way which I'll let commented out here
				// the current one is more efficient due of avoiding the swap gates when doing the IQFT
				// with the IQFT with the swap, just uncomment the following for loop and comment the one after it
				// do a similar thing for the last U gate application

				//for (size_t ctrlQubit = 0; ctrlQubit < lastQubit; ++ctrlQubit)
				for (size_t ctrlQubit = lastQubit; ctrlQubit > 0; --ctrlQubit)
				{
					NQubitsControlledQuantumGate<VectorClass, MatrixClass> UGate(BaseClass::getNrQubits(), controlledGate, BaseClass::getFunctionStartQubit(), ctrlQubit);

					UGate.Execute(reg);

					// power up the controlled gate
					controlledGate *= controlledGate;
				}

				{
					//NQubitsControlledQuantumGate<VectorClass, MatrixClass> UGate(BaseClass::getNrQubits(), controlledGate, BaseClass::getFunctionStartQubit(), lastQubit);
					NQubitsControlledQuantumGate<VectorClass, MatrixClass> UGate(BaseClass::getNrQubits(), controlledGate, BaseClass::getFunctionStartQubit(), 0);

					UGate.Execute(reg);
				}

				// it doesn't really matter if you measure the qubits from f and when you do after the above
				// or if you measure them several times in a row
				//QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(BaseClass::getFunctionStartQubit(), QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);

				// then perform an inverse fourier transform
				BaseClass::IQFT(reg, false);
			}

			MatrixClass U;
		};


		template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PhaseEstimationWithoutTensorProduct : public PhaseEstimationBase<VectorClass, MatrixClass>
		{
		public:
			using RegisterClass = QubitRegister<VectorClass, MatrixClass>;
			using BaseClass = PhaseEstimationBase<VectorClass, MatrixClass>;

			PhaseEstimationWithoutTensorProduct(const MatrixClass& op, size_t N = 7, size_t L = 3)
				: BaseClass(N, L), U(op)
			{
				assert(U.rows() == 2);
			}

			size_t Execute(RegisterClass& reg) override
			{
				ExecuteWithoutMeasurement(reg);

				// any of those following should do, but if one does not do the f register measurement above and here there is no full register measurement
				// the f should be measured separately to find out its content

				//return reg.MeasureAll();
				return reg.Measure(0, BaseClass::getFunctionStartQubit() - 1);
			}

			std::map<size_t, size_t> ExecuteWithMultipleMeasurements(RegisterClass& reg, size_t nrMeasurements = 10000)
			{
				ExecuteWithoutMeasurement(reg);

				return reg.RepeatedMeasure(0, BaseClass::getFunctionStartQubit() - 1, nrMeasurements);
			}

			const MatrixClass& getU() const
			{
				return U;
			}

		protected:
			void ExecuteWithoutMeasurement(RegisterClass& reg)
			{
				// TODO: check if things are set up all right: size of U, size of reg, etc.

				BaseClass::Init(reg);

				MatrixClass gateOp = U;
				const size_t lastQubit = BaseClass::getFunctionStartQubit() - 1;

				// There is another way which I'll let commented out here
				// the current one is more efficient due of avoiding the swap gates when doing the IQFT
				// with the IQFT with the swap, just uncomment the following for loop and comment the one after it
				// do a similar thing for the last U gate application

				//for (size_t ctrlQubit = 0; ctrlQubit < lastQubit; ++ctrlQubit)
				for (size_t ctrlQubit = lastQubit; ctrlQubit > 0; --ctrlQubit)
				{
					MatrixClass opMat = MatrixClass::Identity(4, 4);
					opMat.block(2, 2, 2, 2) = gateOp;
					ControlledGate.setOperator(opMat);
					ControlledGate.setQubit1(BaseClass::getFunctionStartQubit());
					ControlledGate.setQubit2(ctrlQubit);

					reg.ApplyGate(ControlledGate);

					// power up the controlled gate
					gateOp *= gateOp;
				}

				{
					MatrixClass opMat = MatrixClass::Identity(4, 4);
					opMat.block(2, 2, 2, 2) = gateOp;
					ControlledGate.setOperator(opMat);
					ControlledGate.setQubit1(BaseClass::getFunctionStartQubit());
					ControlledGate.setQubit2(0);

					reg.ApplyGate(ControlledGate);
				}

				// it doesn't really matter if you measure the qubits from f and when you do after the above
				// or if you measure them several times in a row
				//QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure(BaseClass::getFunctionStartQubit(), QC::QuantumAlgorithm<VectorClass, MatrixClass>::getNrQubits() - 1);

				// then perform an inverse fourier transform
				BaseClass::IQFT(reg, false);
			}

			MatrixClass U;
			Gates::AppliedGate<MatrixClass> ControlledGate;
		};
	}

}

