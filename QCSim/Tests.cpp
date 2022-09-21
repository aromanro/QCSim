#include "Tests.h"

#include "QubitRegister.h"

#include <iostream>
#include <map>

bool approxEqual(double val1, double val2)
{
    return abs(val1 - val2) < 1E-10;
}

bool approxEqual(std::complex<double> val1, std::complex<double> val2)
{
    return approxEqual(val1.real(), val2.real()) && approxEqual(val1.imag(), val2.imag());
}

bool registerMeasurementsTests()
{
    std::map<int, int> measurements;
    std::map<int, int> fmeasurements;

    const int nrMeasurements = 10000;


    // testing subregister measurement

    // state is 1/2 |00> - i / 2 |01> + 1/sqrt(2) |11>
    // qubits are in little endian order, that is, the first one is to the right

    QC::QubitRegister reg(2); // only a two qubit register

    std::complex<double> c0(0, 0);
    std::complex<double> c1(1, 0);

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

        Eigen::VectorXcd amplitudes = reg.getRegisterStorage();
        if (1 == state)
        {
            if (!approxEqual(amplitudes(0), c0) || !approxEqual(amplitudes(1), std::complex<double>(0, -1. / sqrt(3.))) ||
                !approxEqual(amplitudes(2), c0) || !approxEqual(amplitudes(3), std::complex<double>(sqrt(2.) / sqrt(3.), 0)))
            {
                std::cout << "Wrong state after measuring the first qubit to 1" << std::endl;
                return false;
            }
        }
        else if (0 == state)
        {
            if (!approxEqual(amplitudes(0), c1) || !approxEqual(amplitudes(1), c0) ||
                !approxEqual(amplitudes(2), c0) || !approxEqual(amplitudes(3), c0))
            {
                std::cout << "Wrong state after measuring the first qubit to 0" << std::endl;
                return false;
            }
        }
        else
        {
            std::cout << "Wrong state returned after measuring the first qubit" << std::endl;
            return false;
        }
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

        Eigen::VectorXcd amplitudes = reg.getRegisterStorage();
        if (1 == state)
        {
            if (!approxEqual(amplitudes(0), c0) || !approxEqual(amplitudes(1), c0) ||
                !approxEqual(amplitudes(2), c0) || !approxEqual(amplitudes(3), c1))
            {
                std::cout << "Wrong state after measuring the second qubit to 1" << std::endl;
                return false;
            }

        }
        else if (0 == state)
        {
            if (!approxEqual(amplitudes(0), std::complex<double>(0.5 * sqrt(2), 0)) || !approxEqual(amplitudes(1), std::complex<double>(0, -0.5 * sqrt(2))) ||
                !approxEqual(amplitudes(2), c0) || !approxEqual(amplitudes(3), c0))
            {
                std::cout << "Wrong state after measuring the second qubit to 0" << std::endl;
                return false;
            }
        }
        else
        {
            std::cout << "Wrong state returned after measuring the second qubit" << std::endl;
            return false;
        }
    }

    // should have 1/2 probability for 0 and 1
    std::cout << "Subregister measurement, second qubit" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    // now do the measurements one after another

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

bool tests()
{
    std::cout << "\nTests\n";

    bool res = registerMeasurementsTests();

    return res;
}