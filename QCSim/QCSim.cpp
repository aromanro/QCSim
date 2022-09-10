// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>

#include "QubitRegister.h"
#include "GroverAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>



int main()
{
    std::cout << "Hello Quantum World!\n";

    QC::QubitRegister reg;

    std::map<int, int> measurements;
    const int nrMeasurements = 1000;

    QC::HadamardGate hadamard;
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

    Grover::GroverAlgorithm algo(8);
    algo.setCorrectQuestionState(13);
    for (int i = 0; i < nrMeasurements; ++i)
    {
        unsigned int state = algo.Execute();
        ++measurements[state];
    }
    for (auto m : measurements)
    std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
}



 

