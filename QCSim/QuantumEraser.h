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
		QuantumEraser(int addseed = 0)
			: QC::QuantumAlgorithm<VectorClass, MatrixClass>(2, addseed), eraser(false)
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
			Init();
			// now we're in the state given by the first beam splitter

			// the cnot gate has the role of the spontaneous parametric down convertor
			// from one 'photon' it makes out a pair of them, a 'signal' one and an 'idler' one
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(cnot, 1); //controlling qubit is 0 by default

			// now the action of the next beam splitter:
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, 0);

			// the choice of using or not the eraser could be 'delayed' 
			if (eraser) QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, 1);

			return QC::QuantumAlgorithm<VectorClass, MatrixClass>::Measure();
		}

	protected:
		void Init()
		{
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::setToBasisState(0);
			// the following has the role of the first beam splitter:
			QC::QuantumAlgorithm<VectorClass, MatrixClass>::ApplyGate(hadamard, 0);
		}

		bool eraser;

		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::CNOTGate<MatrixClass> cnot;
	};

}


