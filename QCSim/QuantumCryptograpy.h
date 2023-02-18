#pragma once

#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "Utils.h"

#include <bitset>

#define _USE_MATH_DEFINES
#include <math.h>

namespace QuantumCryptograpy {

	// about half of the tranmitted bits will have a matched measurement basis, out of those 50% will be used to check for eavesdropping and discarded afterwards (some other percentage, like 25% or 20% is also possible, but currently it's hardwired to 50%)
	template<int nrBits = 512, class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class BB84Protocol :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		typedef QC::QuantumAlgorithm<VectorClass, MatrixClass> BaseClass;

		BB84Protocol(int addseed = 0)
			: BaseClass(1, addseed),
			eavesdropping(false), randomEavesdropping(false),
			dist_bool(0, 1)
		{
			uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			timeSeed += addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);
		}

		// unlike other Execute functions, returns 0 if a mismatch was found and 1 if the verified bits matched
		unsigned int Execute() override
		{
			Init();

			// key generation & transmission
			for (unsigned int b = 0; b < nrBits; ++b)
			{
				sendBasis[b] = chooseMeasurementBasisForSending(); // 0 - Z basis, 1 - X basis
				send[b] = generateRandomBitToSend();

				// switch back to computational basis if Alice measured in X basis
				if (sendBasis[b])
					measurementBasis.switchToOperatorBasis(BaseClass::reg, X, 0, true);

				if (eavesdropping)
				{
					// Eve
					if (!randomEavesdropping || getRandomBool())
						Receive(b, ereceive, ereceiveBasis, true);
				}

				// Bob
				Receive(b, receive, receiveBasis);
			}

			// at the end of the transmission, Alice announces the basis she used
			// Bob checks against his basis
			checkCommonBasisMeasurements();

			// then the matches/mismatches of the measurements basis are publicly announced
			// the mismatched values are discarded on both Alice and Bob ends

			// share publicly a subset of the bits measured in the same basis, so the other end can check them against what he/she measured
			return checkBitsMismatched() ? 1 : 0;
		}

		bool getEavesdropping() const
		{
			return eavesdropping;
		}

		void setEavesdropping(bool val = true)
		{
			eavesdropping = val;
		}

		bool getRandomEavesdropping() const
		{
			return randomEavesdropping;
		}

		void setRandomEavesdropping(bool val = true)
		{
			randomEavesdropping = val;
		}

		// this is to be used for tests, in 'real world' those remain private to Alice and Bob and they are not to be compared
		bool compareKeys() const
		{
			for (unsigned int b = 0; b < nrBits; ++b)
				if (commonBasis[b] && !compare[b]) // key bits are those that had a common measurement basis and were not used in comparison for eavesdrop detection
					if (receive[b] != send[b]) return false;

			return true;
		}

		std::vector<bool> getReceivedKey() const
		{
			std::vector<bool> key;

			for (unsigned int b = 0; b < nrBits; ++b)
				if (commonBasis[b] && !compare[b])
					key.push_back(receive[b]);

			return key;
		}

		std::vector<bool> getSentKey() const
		{
			std::vector<bool> key;

			for (unsigned int b = 0; b < nrBits; ++b)
				if (commonBasis[b] && !compare[b])
					key.push_back(send[b]);

			return key;
		}

	protected:
		void Init()
		{
			send.reset();
			sendBasis.reset();
			ereceive.reset();
			ereceiveBasis.reset();
			receive.reset();
			receiveBasis.reset();
			commonBasis.reset();
			compare.reset();
		};

		void Receive(unsigned int b, std::bitset<nrBits>& recv, std::bitset<nrBits>& recvBasis, bool switchBackToZ = false)
		{
			// pick a basis for measurement
			const bool basis = getRandomBool();

			// if false, it's Z, the computational basis, so remain in that one, otherwise use X
			recvBasis[b] = basis;
			if (basis)
				measurementBasis.switchToOperatorBasis(BaseClass::reg, X.getRawOperatorMatrix());

			recv[b] = BaseClass::Measure();

			if (basis && switchBackToZ) // switch back to computational basis if needed
				measurementBasis.switchToOperatorBasis(BaseClass::reg, X.getRawOperatorMatrix(), 0, true);
		}

		void checkCommonBasisMeasurements()
		{
			for (unsigned int b = 0; b < nrBits; ++b)
				if (sendBasis[b] == receiveBasis[b])
					commonBasis[b] = true;
		}

		bool checkBitsMismatched()
		{
			// 4 values, if it's 0 then pick it up, so about 25% will be verified
			//std::uniform_int_distribution<> dist_percent(0, 3);

			// about 50% verified
			std::uniform_int_distribution<> dist_percent(0, 1);

			// those will be discarded because they are publicly shared
			// the remaining ones will be the actual key
			bool match = true;
			for (unsigned int b = 0; b < nrBits; ++b)
				if (commonBasis[b] && dist_percent(rng) == 0)
				{
					compare[b] = true;
					if (receive[b] != send[b])
						match = false; // don't bail out here, the whole 'compare' should be set for supplementary checks/tests if needed
				}

			return match;
		}

		unsigned int chooseMeasurementBasisForSending()
		{
			// use this to generate a random value used to pick a measurement basis
			// alternatively you can use a random number generator or a pregenerated sequence in whatever way you want
			// but I like this way more
			BaseClass::setToCatState();
			return BaseClass::Measure();
		}

		unsigned int generateRandomBitToSend()
		{
			// another way would be to simply generate a random value for the state in some other way (0 or 1) and then simply set the register - one qubit - to that state
			// but as for generating the random value above, I like this method more
			BaseClass::setToCatState();
			return BaseClass::Measure();
		}

		bool getRandomBool()
		{
			return dist_bool(rng) == 1;
		};

		// to be used for measurement basis
		QC::Gates::PauliXGate<MatrixClass> X;
		QC::MeasurementBasis<VectorClass, MatrixClass> measurementBasis; // 0, use Z, 1 use X

		// Alice
		std::bitset<nrBits> send;
		std::bitset<nrBits> sendBasis;

		// Eve - eavesdropping
		// the values here do not really matter, but just in case somebody wants to look at them...
		std::bitset<nrBits> ereceive;
		std::bitset<nrBits> ereceiveBasis;

		// Bob
		std::bitset<nrBits> receive;
		std::bitset<nrBits> receiveBasis;

		// for comparison - mark only a subset of the bits with a common basis
		std::bitset<nrBits> commonBasis;
		std::bitset<nrBits> compare;

		bool eavesdropping;
		bool randomEavesdropping;

		std::mt19937_64 rng;
		std::uniform_int_distribution<> dist_bool;
	};
}
