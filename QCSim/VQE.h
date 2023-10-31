#pragma once


#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"
#include "PauliString.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include "PauliString.h"

namespace VQE {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliStringVQE :
		public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		PauliStringVQE()
		{
			const double theta = 0.5 * M_PI;
			rx.SetTheta(theta);
			ry.SetTheta(-theta);
		}

		unsigned int Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			const unsigned int nrQubits = reg.getNrQubits();
			const unsigned int lastQubit = nrQubits - 1;

			for (unsigned int qubit = 0; qubit < nrQubits; ++qubit)
				if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opX)
					reg.ApplyGate(ry, qubit);
				else if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opY)
					reg.ApplyGate(rx, qubit);

			return 0; // not used here
		}

		void setOperatorForQubit(unsigned int qubit, PauliString::PauliString::PauliOp op)
		{
			pauliString.setOperatorForQubit(qubit, op);
		}

		PauliString::PauliString::PauliOp getOperatorForQubit(unsigned int qubit) const
		{
			return pauliString.getOperatorForQubit(qubit);
		}

		double getCoefficient() const
		{
			return pauliString.getCoefficient();
		}

		void setCoefficient(double c)
		{
			pauliString.setCoefficient(c);
		}

	protected:
		PauliString::PauliString pauliString;

		QC::Gates::RxGate<MatrixClass> rx;
		QC::Gates::RyGate<MatrixClass> ry;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class VariationalQuantumEigensolver :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		VariationalQuantumEigensolver(unsigned int N = 2, bool addSeed = false)
			: BaseClass(N, addSeed)
		{
		}

		void AddTerm(double coeff, const PauliStringVQE<VectorClass, MatrixClass>& term)
		{
			terms.push_back(term);
			terms.back().setCoefficient(coeff);
		}

		void Clear()
		{
			terms.clear();
		}

		unsigned int Execute() override
		{
			// TODO: implement it, using Nelder-Mead algorithm
			return 0;
		}

	protected:
		std::vector<PauliStringVQE<VectorClass, MatrixClass>> terms;
	};

}