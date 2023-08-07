#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

// for details see for example "Experimenting quantum phenomena on NISQ computers using high level quantum programming" by Duc M. Tran and Hung Q. Nguyen
// https://arxiv.org/abs/2111.02896v2

namespace Paradoxes {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumEraser :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		QuantumEraser(int addseed = 0)
			: BaseClass(2, addseed), eraser(false)
		{
		}

		void setEraser(bool e = true)
		{
			eraser = e;
		}

		bool getEraser() const
		{
			return eraser;
		}

		unsigned int Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure();
		}

		std::map<unsigned int, unsigned int> ExecuteWithMultipleMeasurements(unsigned int nrMeasurements = 10000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(nrMeasurements);
		}

	protected:
		void Init()
		{
			BaseClass::setToBasisState(0);
			// the following has the role of the first beam splitter:
			BaseClass::ApplyGate(hadamard, 0);
		}

		void ExecuteWithoutMeasurement()
		{
			Init();
			// now we're in the state given by the first beam splitter

			// the cnot gate has the role of the spontaneous parametric down convertor
			// from one 'photon' it makes out a pair of them, a 'signal' one and an 'idler' one
			BaseClass::ApplyGate(cnot, 1); //controlling qubit is 0 by default

			// the above two gates (hadamard applied in Init and cnot) are equivalent to the entangling gate (noted E or E2)

			// now the action of the next beam splitter:
			BaseClass::ApplyGate(hadamard, 0);

			// the choice of using or not the eraser could be 'delayed' 
			if (eraser) BaseClass::ApplyGate(hadamard, 1);
		}

		bool eraser;

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::CNOTGate<MatrixClass> cnot;
	};

}


