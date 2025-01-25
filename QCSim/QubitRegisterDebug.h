#pragma once


#include "QubitRegister.h"

namespace QC {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QubitRegisterDebug : public QubitRegister<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QubitRegister<VectorClass, MatrixClass>;

		QubitRegisterDebug(size_t N = 3, unsigned int addseed = 0)
			: BaseClass(N, addseed)
		{
		}

		// allows saving it into a file to look at the data
		// for example to check the results of a quantum simulation
		bool writeToFile(const std::string& name, bool amplitude = true, bool append = false) const
		{
			try {
				std::ofstream thefile;
				thefile.open(name, std::ios::out | (append ? std::ios::app : std::ios::trunc));

				if (!thefile.is_open()) return false;

				if (append) thefile << std::endl << std::endl;

				for (size_t i = 0; i < BaseClass::NrBasisStates; ++i)
				{
					thefile << i << "\t";
					if (amplitude) thefile << std::abs(BaseClass::registerStorage(i));
					else thefile << BaseClass::registerStorage(i);
					thefile << std::endl;
				}

				return true;
			}
			catch (...) {};

			return false;
		}

		void displayState(size_t state) const
		{
			const size_t nQubits = BaseClass::getNrQubits();
			std::cout << "|";

			size_t mask = 1ULL << (nQubits - 1);
			for (size_t qubit = 0; qubit < nQubits; ++qubit)
			{
				std::cout << ((state & mask) ? "1" : "0");
				mask >>= 1;
			}

			std::cout << ">    ";
		}

		void displayRegister() const
		{
			const size_t nQubits = BaseClass::getNrQubits();
			const size_t nStates = BaseClass::getNrBasisStates();

			//std::cout << std::fixed;
			std::cout << std::setprecision(4);

			for (size_t state = 0; state < nStates; ++state)
			{
				const std::complex<double> val = BaseClass::getBasisStateAmplitude(state);
				if (abs(real(val)) < 1E-10 && abs(imag(val)) < 1E-10) continue;

				bool r = false;
				if (abs(real(val)) > 1E-10) {
					std::cout << real(val) << " ";
					r = true;
				}
				if (abs(imag(val)) > 1E-10) {
					if (r && imag(val) > 0) std::cout << "+ ";

					if (imag(val) < 0) {
						std::cout << "-";
						if (r) std::cout << " ";
					}
					std::cout << abs(imag(val)) << "i ";
				}

				displayState(nQubits);
			}
		}
	};

}

