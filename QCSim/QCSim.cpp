// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>

#include "QubitRegister.h"
#include "GroverAlgorithm.h"
#include "ShorAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>



int main()
{
    std::cout << "Hello Quantum World!\n";

    //QC::QubitRegister reg;

    std::map<int, int> measurements;
    const int nrMeasurements = 10000;

    //QC::HadamardGate hadamard;
    //QC::PhaseShiftGate phaseShift;
    
    //std::cout << hadamard.getOperatorMatrix(3, 1) << std::endl;

    /*
    for (int i = 0; i < nrMeasurements; ++i)
    {
        //reg.setToCatState();
        //reg.setToBasisState(3);
        reg.setToEqualSuperposition();
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    */

    /*
    // a)
    std::cout << "a)" << std::endl;
    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 1);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    // b)
    measurements.clear();
    std::cout << "b)" << std::endl;
    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 0);
        reg.ApplyGate(hadamard, 1);
        reg.ApplyGate(hadamard, 2);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    // c)
    measurements.clear();
    std::cout << "c)" << std::endl;
    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 2);
        reg.ApplyGate(hadamard, 2);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    // d)
    measurements.clear();
    std::cout << "d)" << std::endl;
    phaseShift.SetPhaseShift(M_PI);
    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 2);
        reg.ApplyGate(phaseShift, 2);
        reg.ApplyGate(hadamard, 2);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
        */

    /*
    Grover::GroverAlgorithm algo;
    algo.setCorrectQuestionState(3);
    for (int i = 0; i < nrMeasurements; ++i)
    {
        unsigned int state = algo.Execute();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
        */

    /*
    QC::CNOTGate cnot;

    Eigen::MatrixXcd m = cnot.getOperatorMatrix(3, 0, 1);
    std::cout << m << std::endl;

    std::cout << "*************************************" << std::endl;
    m = cnot.getOperatorMatrix(3, 2, 1);
    std::cout << m << std::endl;

    measurements.clear();
    std::cout << "a) CNOT: entangled state" << std::endl;

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 1);
        reg.ApplyGate(cnot, 2, 1);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    measurements.clear();
    std::cout << "b) CNOT: cat state" << std::endl;

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 1);
        reg.ApplyGate(cnot, 2, 1);
        reg.ApplyGate(cnot, 0, 1);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    measurements.clear();
    std::cout << "c) unobserved superposition" << std::endl;

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 1);
        reg.ApplyGate(hadamard, 1);
        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    measurements.clear();
    std::cout << "d) CNOT: observed superposition" << std::endl;

    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(0);
        reg.ApplyGate(hadamard, 1);
        reg.ApplyGate(cnot, 2, 1);
        //reg.ApplyGate(cnot, 2, 1);
        reg.ApplyGate(hadamard, 1);

        unsigned int state = reg.Measure();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
        */

/*
    measurements.clear();
    Grover::GroverAlgorithm galgo(8);
    galgo.setCorrectQuestionState(13);
    for (int i = 0; i < nrMeasurements; ++i)
    {
        unsigned int state = galgo.Execute();
        ++measurements[state];
    }
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
*/

    measurements.clear();

    Shor::ShorAlgorithm shorAlgo;

    // with 2 (default) f should be 1, 2, 4, 8 for performing period finding
    shorAlgo.setA(7); // with 7 should be 1, 7, 4, 13
    //shorAlgo.setA(4); // should be 1, 4
    //shorAlgo.setA(8); // 1, 8, 4, 2
    //shorAlgo.setA(11); // 1, 11
    //shorAlgo.setA(13); // 1, 13, 4, 7
    //shorAlgo.setA(14); // 1, 14

    std::map<int, int> fmeasurements;

    for (int i = 0; i < nrMeasurements; ++i)
    {
        unsigned int state = shorAlgo.Execute();
        ++measurements[state & 0x7]; //only the l bits matter

        // this should work only if Execute calls full register measurement, otherwise use the commented part below
        ++fmeasurements[state >> 3]; //save the f register measurements for debugging
        
        //state = shorAlgo.reg.Measure(3, 6);
        //++fmeasurements[state];
    }

    std::cout << "f register" << std::endl;
    for (auto m : fmeasurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    std::cout << "x register" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;

    /*
    measurements.clear();

    QC::QuantumFourierTransform fourier(2);

    for (int i = 0; i < nrMeasurements; ++i)
    {
        // qft from total superposition should go into state 0
        // and from state 0 into total superposition

        // the same for inverse qft

        //fourier.reg.setToBasisState(0);
        fourier.reg.setToEqualSuperposition();
        fourier.QFT();
        fourier.IQFT();
        const unsigned int state = fourier.reg.Measure();
        ++measurements[state];
    }

    std::cout << "Fourier" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
    */

    /*
    measurements.clear();
    QC::SwapGate sg;
  
    for (int i = 0; i < nrMeasurements; ++i)
    {
        reg.setToBasisState(2);
        reg.ApplyGate(sg, 0, 1);
        const unsigned int state = reg.Measure();
        ++measurements[state];
    }

    std::cout << "Swap" << std::endl;
    for (auto m : measurements)
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
    */

    // testing subregister measurement
    measurements.clear();

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
}



 

