// QCSim.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <map>

#include "QubitRegister.h"


int main()
{
    std::cout << "Hello Quantum World!\n";

    QC::QubitRegister reg;
    
    std::map<int, int> measurements;
    const int nrMeasurements = 100000000;
    for (int i = 0; i < nrMeasurements; ++i)
    {
        //reg.setToCatState();
        //reg.setToBasisState(3);
        reg.setToEqualSuperposition();
        unsigned int state = reg.Measure();
        ++measurements[state];
    }

    for (auto m : measurements)
    {
        std::cout << "State: " << m.first << " measured " << m.second << " times, that is " << 100. * m.second / nrMeasurements << "%" << std::endl;
    }

    QC::HadamardGate hadamard;
    std::cout << hadamard.getOperatorMatrix(3, 2) << std::endl;
}

