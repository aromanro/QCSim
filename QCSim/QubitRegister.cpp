#include "QubitRegister.h"

#include <chrono>

namespace QC {
	QubitRegister::QubitRegister(int N, int addseed)
		: NrQubits(N), NrBasisStates(1u << NrQubits),
		uniformZeroOne(0,1)
	{
		registerStorage = Eigen::VectorXd::Zero(NrBasisStates);

		uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		timeSeed += addseed;
		std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };

		rng.seed(seed);
	}

	void QubitRegister::setToBasisState(unsigned int State)
	{
		if (State >= NrBasisStates) return;

		registerStorage.setZero();
		registerStorage(State) = 1;
	}

	void QubitRegister::setToCatState()
	{
		registerStorage.setZero();
		static const double OneOverSqrt2 = 1. / sqrt(2.);

		registerStorage(0) = OneOverSqrt2;
		registerStorage(NrBasisStates - 1) = OneOverSqrt2;
	}
	
	void QubitRegister::setToEqualSuperposition()
	{
		registerStorage.setConstant(1. / sqrt(NrBasisStates));
	}

	void QubitRegister::setRawAmplitude(unsigned int State, std::complex<double> val)
	{
		if (State >= NrBasisStates) return;

		registerStorage(State) = val;
	}

	void QubitRegister::Normalize()
	{
		const double accum = (registerStorage.adjoint() * registerStorage)(0).real();
		if (accum < 1E-20) return;

		registerStorage = 1. / sqrt(accum) * registerStorage;
	}

	unsigned int QubitRegister::Measure()
	{
		const double prob = uniformZeroOne(rng);
		double accum = 0;
		for (unsigned int i = 0; i < NrBasisStates; ++i)
		{
			accum += (std::conj(registerStorage(i)) * registerStorage(i)).real();
			if (prob < accum)
			{
				setToBasisState(i); // collapse
				return i;
			}
		}
		const unsigned int state = NrBasisStates - 1;
		setToBasisState(state);

		return state;
	}

	void QubitRegister::ApplyGate(const QuantumGate& gate, unsigned int qubit, unsigned int controllingQubit)
	{
		registerStorage = gate.getOperatorMatrix(NrQubits, qubit, controllingQubit) * registerStorage;
	}

	void QubitRegister::ApplyOperatorMatrix(const Eigen::MatrixXcd& m)
	{
		registerStorage = m * registerStorage;
	}
}
