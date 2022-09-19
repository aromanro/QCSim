#include "Tests.h"

#include "QubitRegister.h"

#include <iostream>
#include <map>

bool approxEqual(double val1, double val2)
{
    return abs(val1 - val2) < 1E-10;
}

bool tests()
{
	std::map<int, int> measurements;
    std::map<int, int> fmeasurements;

	const int nrMeasurements = 10000;


    // testing subregister measurement

    // state is 1/2 |00> - i / 2 |01> + 1/sqrt(2) |11>
    // qubits are in little endian order, that is, the first one is to the right

    QC::QubitRegister reg(2); // only a two qubit register

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0); //overriden below

        reg.setRawAmplitude(0, 0.5);
        reg.setRawAmplitude(1, std::complex<double>(0, -0.5));
        reg.setRawAmplitude(3, 1. / sqrt(2));

        reg.Normalize(); // already normalized, but better to ensure it

        // measure first qubit only

        const unsigned int state = reg.Measure(0, 0);
        ++measurements[state];
    }

    // should get 3/4 1 and 1/4 0
    std::cout << "Subregister measurement, first qubit" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    measurements.clear();
    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0); //overriden below

        reg.setRawAmplitude(0, 0.5);
        reg.setRawAmplitude(1, std::complex<double>(0, -0.5));
        reg.setRawAmplitude(3, 1. / sqrt(2));

        reg.Normalize(); // already normalized, but better to ensure it

        // measure second qubit only

        const unsigned int state = reg.Measure(1, 1);
        ++measurements[state];
    }

    // should have 1/2 probability for 0 and 1
    std::cout << "Subregister measurement, second qubit" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    // now do the measurements one after another

    // TODO: check the register to be left in the expected state after subregister measurements

    measurements.clear();
    fmeasurements.clear();

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0); //overriden below

        reg.setRawAmplitude(0, 0.5);
        reg.setRawAmplitude(1, std::complex<double>(0, -0.5));
        reg.setRawAmplitude(3, 1. / sqrt(2));

        reg.Normalize(); // already normalized, but better to ensure it

        unsigned int state = reg.Measure(0, 0);
        ++measurements[state];

        state = reg.Measure(1, 1);
        ++fmeasurements[state];
    }

    std::cout << "Subregister measurement, first qubit again" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
    std::cout << "Subregister measurement, second qubit again" << std::endl;
    for (auto m : fmeasurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    measurements.clear();
    fmeasurements.clear();

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0); //overriden below

        reg.setRawAmplitude(0, 0.5);
        reg.setRawAmplitude(1, std::complex<double>(0, -0.5));
        reg.setRawAmplitude(3, 1. / sqrt(2));

        reg.Normalize(); // already normalized, but better to ensure it

        unsigned int state = reg.Measure(1, 1);
        ++measurements[state];

        state = reg.Measure(0, 0);
        ++fmeasurements[state];
    }

    std::cout << "Subregister measurement, second qubit now" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
    std::cout << "Subregister measurement, first qubit now" << std::endl;
    for (auto m : fmeasurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;


	return true;
}