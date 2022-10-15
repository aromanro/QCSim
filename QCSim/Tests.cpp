#include "Tests.h"

#include "QubitRegister.h"

#include "GroverAlgorithm.h"
#include "ShorAlgorithm.h"
#include "Teleportation.h"
#include "SuperdenseCoding.h"
#include "CheckCHSHInequality.h"
#include "BernsteinVazirani.h"
#include "QuantumCryptograpy.h"

#include <iostream>
#include <map>

bool approxEqual(double val1, double val2, double err = 1E-10)
{
    return abs(val1 - val2) < err;
}

bool approxEqual(std::complex<double> val1, std::complex<double> val2, double err = 1E-10)
{
    return approxEqual(val1.real(), val2.real(), err) && approxEqual(val1.imag(), val2.imag(), err);
}

void setRegister(QC::QubitRegister<>& reg)
{
    reg.setToBasisState(0); //overriden below

    reg.setRawAmplitude(0, 0.5);
    reg.setRawAmplitude(1, std::complex<double>(0, -0.5));
    reg.setRawAmplitude(3, 1. / sqrt(2));

    reg.Normalize(); // already normalized, but better to ensure it
}

bool checkQubit0()
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
        setRegister(reg);

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

    return true;
}


bool checkQubit1()
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
        setRegister(reg);

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

    return true;
}

bool checkAmplitudesAfterSingleQubitMeasurement()
{
    return checkQubit0() && checkQubit1();
}


bool checkSingleQubitMeasurements()
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
        setRegister(reg);

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
        setRegister(reg);

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

bool registerMeasurementsTests()
{
    std::cout << "\nTesting measurements..." << std::endl;

    return checkAmplitudesAfterSingleQubitMeasurement() && checkSingleQubitMeasurements();
}


bool ShorTests()
{
    std::cout << "\nTesting Shor..." << std::endl;

    for (int i = 0; i < 10;)
    {
        Shor::ShorAlgorithm shorAlgo;//(21, 14, 9)
        unsigned int p1;
        unsigned int p2;
        bool res = shorAlgo.factorize(p1, p2);
        if (p2 > p1) std::swap(p1, p2);

        std::cout << (res ? "Quantum algo: " : "Classical algo: ") << p1 << " " << p2 << std::endl;

        if (p1 != 5 && p2 != 3) return false;

        if (res) ++i;
    }

    return true;
}

bool GroverTests()
{
    std::cout << "\nTesting Grover..." << std::endl;

    bool res = true;

    std::map<int, int> measurements;
    const int nrMeasurements = 100;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 0xf);

    for (int i = 0; i < 5; ++i)
    {
        const int ans = dist(gen);

        std::cout << "Testing for answer: " << ans << std::endl;

        measurements.clear();
        Grover::GroverAlgorithm galgo(8);
        galgo.setCorrectQuestionState(ans);
        for (int i = 0; i < nrMeasurements; ++i)
        {
            unsigned int state = galgo.Execute();
            ++measurements[state];
        }

        bool found = false;
        for (auto m : measurements)
        {
            std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

            if (m.first == ans)
            {
                found = true;
                if (static_cast<double>(m.second) / nrMeasurements < 0.9)
                    res = false;
            }
        }

        if (!found) res = false;
    }

    return res;
}


bool TeleportationTests()
{
    std::cout << "\nTesting teleportation..." << std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_bool(0, 1);

    Teleportation::QuantumTeleportationRealization qt;

    for (int i = 0; i < 30; ++i)
    {
        if (dist_bool(gen))
        {
            qt.SetState(0, 1); // set the qubit to teleport to 'up'
            unsigned int state = qt.Execute();

            std::cout << "Teleported 1, measured: " << state << std::endl;

            if (state != 1) return false;
        }
        else
        {
            qt.SetState(1, 0); // set the qubit to teleport to 'down'
            unsigned int state = qt.Execute();

            std::cout << "Teleported 0, measured: " << state << std::endl;

            if (state != 0) return false;
        }
    }

    std::uniform_real_distribution<> dist_ampl(-1., 1.);
    
    std::cout << "Now testing teleporting a generic state:" << std::endl;

    for (int i = 0; i < 16; ++i)
    {
        // generate and normalize
        std::complex<double> alpha(dist_ampl(gen), dist_ampl(gen));
        std::complex<double> beta(dist_ampl(gen), dist_ampl(gen));

        const double norm = sqrt((alpha * std::conj(alpha) + beta * std::conj(beta)).real());
        alpha /= norm;
        beta /= norm;

        qt.SetState(alpha, beta);
        std::cout << "Teleporting " << alpha << "|0> + " << beta << "|1>";

        unsigned int classicalBits = qt.Teleport(i < 8 ? false : true); // also test sending explicitely the two classical bits for half the tests, although it should not make a difference
        std::cout << " Measured values for the two qubits: " << classicalBits;

        // how is the whole thing looking before Bob's measurement?
        std::complex<double> receivedAlpha = qt.getBasisStateAmplitude(classicalBits);
        std::complex<double> receivedBeta = qt.getBasisStateAmplitude(0x4 | classicalBits);
        std::cout << "... Teleported state: " << receivedAlpha << "|0> + " << receivedBeta << "|1>" << std::endl;

        if (!approxEqual(alpha, receivedAlpha) || !approxEqual(beta, receivedBeta)) return false;
    }

    return true;
}

bool SuperdenseCodingTests()
{
    std::cout << "\nTesting superdense coding..." << std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_bool(0, 1);

    Coding::SuperdenseCoding coding;

    for (int i = 0; i < 30; ++i)
    {
        // bit1 is on position 0, bit2 is on position 1
        const bool bit1 = dist_bool(gen) == 1;
        const bool bit2 = dist_bool(gen) == 1;
        coding.SetBits(bit1, bit2);
        unsigned int classicalBits = coding.Execute();

        const bool recvbit1 = (classicalBits & 1) != 0;
        const bool recvbit2 = (classicalBits & 2) != 0;
        std::cout << "Sent " << (bit1 ? "1" : "0") << (bit2 ? "1" : "0") << " using superdense coding, received: " << (recvbit1 ? "1" : "0") << (recvbit2 ? "1" : "0") << std::endl;

        if (bit1 != recvbit1 || bit2 != recvbit2) return false;
    }

    return true;
}

bool BellInequalitiesTests()
{
    std::cout << "\nTesting CHSH inequality..." << std::endl;

    static const double expected = 2. * sqrt(2.);

    BellInequalities::CheckCHSHInequality test;
    bool separateMeasurements = false;
    for (int i = 0; i < 20; ++i)
    {
        if (i >= 10) separateMeasurements = true;
        test.ResetStatistics();
        for (int i = 0; i < 100000; ++i)
            test.Check(separateMeasurements);

        const double val = test.getValue();
        std::cout << (separateMeasurements ? "Measurements separated" : "Measurements together") << " Value: " << val << (val <= 2. ? " Inequality obeyed (you won the lottery!)" : " Inequality violated") << std::endl;

        if (val <= 2.) return false;
        else if (!approxEqual(val, expected, 0.1))
        {
            std::cout << "Expected value to be about 2*sqrt(2), got a value that's too different" << std::endl;
            return false;
        }
    }

    return true;
}

bool BernsteinVaziraniTests()
{
    std::cout << "\nTesting Bernstein-Vazirani..." << std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist_bool(0, 1);
    
    std::cout << "Three qubits:" << std::endl;
    BernsteinVazirani::BernsteinVaziraniAlgorithm bv;

    unsigned int b0 = 1;
    unsigned int b1 = 2;
    unsigned int b2 = 4;

    for (int i = 0; i < 10; ++i)
    {
        const bool bit0 = dist_bool(gen) == 1;
        const bool bit1 = dist_bool(gen) == 1;
        const bool bit2 = dist_bool(gen) == 1;

        unsigned int str = 0;
        if (bit0) str |= b0;
        if (bit1) str |= b1;
        if (bit2) str |= b2;

        std::cout << "Setting string to: " << str << std::endl;
        bv.setString(str);
        unsigned int state = bv.Execute();

        if (state != str)
        {
            std::cout << "Failed, obtained a different string: " << state << std::endl;
            return false;
        }
    }

    std::cout << "Six qubits:" << std::endl;
    unsigned int b3 = 8;
    unsigned int b4 = 16;
    unsigned int b5 = 32;

    BernsteinVazirani::BernsteinVaziraniAlgorithm bvBig(6);
    for (int i = 0; i < 30; ++i)
    {
        const bool bit0 = dist_bool(gen) == 1;
        const bool bit1 = dist_bool(gen) == 1;
        const bool bit2 = dist_bool(gen) == 1;
        const bool bit3 = dist_bool(gen) == 1;
        const bool bit4 = dist_bool(gen) == 1;
        const bool bit5 = dist_bool(gen) == 1;

        unsigned int str = 0;
        if (bit0) str |= b0;
        if (bit1) str |= b1;
        if (bit2) str |= b2;
        if (bit3) str |= b3;
        if (bit4) str |= b4;
        if (bit5) str |= b5;

        std::cout << "Setting string to: " << str << std::endl;
        bvBig.setString(str);
        unsigned int state = bvBig.Execute();

        if (state != str)
        {
            std::cout << "Failed, obtained a different string: " << state << std::endl;
            return false;
        }
    }

    return true;
}

bool QuantumCryptograpyTests()
{
    std::cout << "\nTesting BB84 protocol..." << std::endl;

    QuantumCryptograpy::BB84Protocol bb84;

    for (int i = 0; i < 30; ++i)
    {
        const bool eavesdrop = i > 9;
        const bool randomEavesdrop = i > 19;
        
        bb84.setEavesdropping(eavesdrop);
        bb84.setRandomEavesdropping(randomEavesdrop);

        std::cout << "Transmission with" << (eavesdrop ? " eavesdropping" : "out eavesdropping") << (eavesdrop ? (randomEavesdrop ? " randomly" : " each time") : "") << std::endl;

        const unsigned int match = bb84.Execute();
        if (match)
        {
            if (eavesdrop)
            {
                std::cout << "Match, but eavesdropped" << std::endl;
                return false;
            }

            // additional check, this is only done for tests, the key remains a secret in 'real life'
            if (!bb84.compareKeys())
            {
                std::cout << "Despite a match for checked bits, the keys do not match" << std::endl;
                return false;
            }
        }
        else
        {
            // mismatch
            if (!eavesdrop)
            {
                std::cout << "Mismatch, but not eavesdropped" << std::endl;
                return false;
            }
        }
        //std::cout << "Key length: " << bb84.getReceivedKey().size() << std::endl;
    }

    return true;
}

bool tests()
{
    std::cout << "\nTests\n";

    bool res = registerMeasurementsTests();
    if (res) res = BernsteinVaziraniTests();
    if (res) res = GroverTests();
    if (res) res = ShorTests();
    if (res) res = TeleportationTests();
    if (res) res = SuperdenseCodingTests();
    if (res) res = QuantumCryptograpyTests();
    if (res) res = BellInequalitiesTests();

    std::cout << "\nTests " << (res ? "succeeded" : "failed") << std::endl;

    return res;
}