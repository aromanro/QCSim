#include "Tests.h"
#include "ErrorCorrection3Qubits.h"
#include "ShorCode.h"

#include <iostream>
#include <map>


std::pair<std::complex<double>, std::complex <double>> getRandomState()
{
	// generate and normalize
	std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
	std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));

	const double norm = sqrt(std::norm(alpha) + std::norm(beta));
	alpha /= norm;
	beta /= norm;

	return { alpha, beta };
}

bool FlipErrorCorrectionTests()
{
	std::complex<double> alpha;
	std::complex<double> beta;

	std::cout << "Qubit flip:" << std::endl;
	ErrorCorrection::ErrorCorrection3QubitsFlip errorCorrectionFlip;
	for (int i = 0; i < 16; ++i)
	{
		std::tie(alpha, beta) = getRandomState();

		std::cout << "Initial state: " << alpha << "|0> + " << beta << "|1>" << std::endl;

		for (unsigned int q = 0; q <= 3; ++q)
		{
			std::cout << "Flipping qubit: " << ((q == 3) ? "No qubit" : std::to_string(q)) << "...";
			errorCorrectionFlip.SetState(alpha, beta);
			errorCorrectionFlip.SetErrorQubit(q); // q = 3 means 'no qubit flip'

			const unsigned int res = errorCorrectionFlip.Execute(); // return values: 3 means that the first qubit was flipped, 0 - no flip, 1 - first qubit was flipped, 2 - the second one was flipped

			// check against return values that show that the qubit flip was not detected:
			// either some qubit flip was done but not detected (first condition in if)
			// or no qubit flip was done but one was reported (second condition in if)
			// or first qubit flip was done but some other detected
			// or one of the other two qubits was flipped but something else was reported
			if ((res == 0 && q != 3) || (q == 3 && res != 0) ||
				(q == 0 && res != 3) ||
				((res == 1 || res == 2) && res != q))
			{
				std::cout << "\n Qubit flip was not corrected, result: " << res << std::endl;
				return false;
			}

			// now check the fidelity of the wavefunction for the first qubit
			// due of the measurement, the wavefunction collapsed, whence the complication:
			QC::QubitRegister reg(3);
			const unsigned int meas = res << 1;
			reg.setRawAmplitude(meas, alpha);
			reg.setRawAmplitude(meas | 1, beta);

			const double fidelity = reg.stateFidelity(errorCorrectionFlip.getRegisterStorage());
			if (fidelity < 0.99999)
			{
				std::cout << "\n Quit flip was not corrected, result: " << res << " Fidelity is too small: " << fidelity << std::endl;
				return false;
			}

			std::cout << " ok" << std::endl;
		}
	}

	return true;
}

bool SignErrorCorrectionTests()
{
	std::complex<double> alpha;
	std::complex<double> beta;

	std::cout << "\nSign change:" << std::endl;
	ErrorCorrection::ErrorCorrection3QubitsSign errorCorrectionSign;
	for (int i = 0; i < 16; ++i)
	{
		std::tie(alpha, beta) = getRandomState();

		std::cout << "Initial state: " << alpha << "|0> + " << beta << "|1>" << std::endl;

		for (unsigned int q = 0; q <= 3; ++q)
		{
			std::cout << "Changing sign for qubit: " << ((q == 3) ? "No qubit" : std::to_string(q)) << "...";
			errorCorrectionSign.SetState(alpha, beta);
			errorCorrectionSign.SetErrorQubit(q); // q = 3 means 'no error'

			const unsigned int res = errorCorrectionSign.Execute(); // return values: 3 means that the first qubit had a sign change, 0 - no change, 1 - first qubit affected, 2 - the second one affected

			// check against return values that show that the qubit sign change was not detected:
			// either some qubit sign change was done but not detected (first condition in if)
			// or no qubit sign change was done but one was reported (second condition in if)
			// or first qubit sign change was done but some other detected
			// or one of the other two qubits was changed but something else was reported
			if ((res == 0 && q != 3) || (q == 3 && res != 0) ||
				(q == 0 && res != 3) ||
				((res == 1 || res == 2) && res != q))
			{
				std::cout << "\n Sign change was not corrected, result: " << res << std::endl;
				return false;
			}

			// now check the fidelity of the wavefunction for the first qubit
			// due of the measurement, the wavefunction collapsed, whence the complication:
			QC::QubitRegister reg(3);
			const unsigned int meas = res << 1;
			reg.setRawAmplitude(meas, alpha);
			reg.setRawAmplitude(meas | 1, beta);

			const double fidelity = reg.stateFidelity(errorCorrectionSign.getRegisterStorage());
			if (fidelity < 0.99999)
			{
				std::cout << "\n Sign change was not corrected, result: " << res << " Fidelity is too small: " << fidelity << std::endl;
				return false;
			}

			std::cout << " ok" << std::endl;
		}
	}

	return true;
}

std::string getErrorDesc(ErrorCorrection::ShorCode<>::ErrorType errorType)
{
	std::string errorDesc;
	if (errorType == ErrorCorrection::ShorCode<>::ErrorType::Flip)
		errorDesc = "Flipping";
	else if (errorType == ErrorCorrection::ShorCode<>::ErrorType::Sign)
		errorDesc = "Changing sign for";
	else if (errorType == ErrorCorrection::ShorCode<>::ErrorType::Both)
		errorDesc = "Flipping and changing sign for ";
	else
		errorDesc = "Doing nothing for";

	return errorDesc;
}

bool ShorCodeTests()
{
	std::complex<double> alpha;
	std::complex<double> beta;

	std::cout << "\nShor Code:" << std::endl;
	ErrorCorrection::ShorCode errorCorrection;
	for (int i = 0; i < 16; ++i)
	{
		std::tie(alpha, beta) = getRandomState();

		std::cout << "Initial state: " << alpha << "|0> + " << beta << "|1>" << std::endl;

		for (unsigned int q = 0; q <= 9; ++q)
		{
			errorCorrection.SetErrorQubit(q); // q = 9 means 'no error'

			for (unsigned int c = 0; c <= 2; ++c)
			{
				ErrorCorrection::ShorCode<>::ErrorType errorType = static_cast<ErrorCorrection::ShorCode<>::ErrorType>(c);
				if (q == 9) errorType = ErrorCorrection::ShorCode<>::ErrorType::None;
				std::string errorDesc = getErrorDesc(errorType);

				std::cout << errorDesc << " qubit: " << ((q == 9) ? "No qubit" : std::to_string(q)) << "...";
				errorCorrection.SetState(alpha, beta);
				errorCorrection.setErrorType(errorType);

				const unsigned int res = errorCorrection.Execute();

				// limited check of the result, more reliance on checking the fidelity
				// this probably could be extended, but I won't bother:
				if ((res == 0 && q != 9) || (q == 9 && res != 0))
				{
					std::cout << "\n Change was not corrected, result: " << res << std::endl;
					return false;
				}

				// check the fidelity of the wavefunction for the first qubit
				// due of the measurement, the wavefunction collapsed, whence the complication:
				QC::QubitRegister reg(9);
				const unsigned int meas = res << 1;
				reg.setRawAmplitude(meas, alpha);
				reg.setRawAmplitude(meas | 1, beta);

				const double fidelity = reg.stateFidelity(errorCorrection.getRegisterStorage());
				if (fidelity < 0.99999)
				{
					std::cout << "\n Qubit change was not corrected, result: " << res << " Fidelity is too small: " << fidelity << std::endl;
					return false;
				}

				std::cout << " ok" << std::endl;

				if (q == 9) break;
			}
		}
	}

	return true;
}

bool ErrorCorrectionTests()
{
	std::cout << "\nTesting Error Correction for 3 qubits encoding of a qubit..." << std::endl;

	return FlipErrorCorrectionTests() && SignErrorCorrectionTests() && ShorCodeTests();
}
