// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>

#include "QubitRegister.h"

#define _USE_MATH_DEFINES
#include <math.h>



int main()
{
    std::cout << "Hello Quantum World!\n";

    QC::QubitRegister reg;

    std::map<int, int> measurements;
    const int nrMeasurements = 1000000;

    QC::HadamardGate hadamard;
    QC::PhaseShiftGate phaseShift;
    std::cout << hadamard.getOperatorMatrix(3, 1) << std::endl;

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
}



 

