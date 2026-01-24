# QCSim
Quantum computation simulator

Has a statevector simulator (better than the 'naive' matrix multiplication kind), a Matrix Product State simulator, a Clifford gates simulator using the stabilizer formalism and a Pauli propagation simulator (work in progress, for now it works only with Clifford gates and doesn't allow measurements, only expectation values and sampling).

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6a193db170ab432596079c530fc75c77)](https://app.codacy.com/gh/aromanro/QCSim/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![CodeFactor](https://www.codefactor.io/repository/github/aromanro/qcsim/badge)](https://www.codefactor.io/repository/github/aromanro/qcsim)

[![Build and test on Ubuntu](https://github.com/aromanro/QCSim/actions/workflows/cmake-single-platform.yml/badge.svg)](https://github.com/aromanro/QCSim/actions/workflows/cmake-single-platform.yml)


Blog entry on the Computational Physics Blog: [Quantum Computing Simulator](https://compphys.go.ro/quantum-computing-simulator/)

### More info

Implemented algorithms:

*   Grover
*   Deutsch-Jozsa
*   Simon
*   Quantum Fourier Transform
*   Phase estimation
*   Shor
*   Bernstein-Vazirani
*   Quantum counting
*   Quantum teleportation
*   Entanglement swapping
*   Superdense coding
*   Quantum cryptography: BB84 protocol
*   CHSH inequality violation

Quantum error correction:

*   3-qubit error correcting a qubit-flip
*   3-qubit error correcting a sign-flip
*   Shor Code

Quantum adders:

*   Quantum half-adder for 1-qubit
*   Quantum full-adder for 1-qubit
*   Full adder for two N-qubit numbers
*   Draper adder
*   Draper adder with carry

Simulation of quantum simulation:

*   Evolution in time with a Hamiltonian given as a sum of Pauli products with real coefficients
*   Evolution of a 1D Gaussian packet in time, solving the 1D time-dependent Schrodinger equation using a Trotter decomposition and quantum Fourier transform

Paradoxes (although some of the above might be considered paradoxes as well):
*   Quantum eraser
*   General Elitzur-Vaidman Bomb tester/interaction free measurement/counterfactual computation
*   Hardy's paradox

Quantum games:
*   Coin flipping
*   Quantum pseudo-telepathy: The magic square game

Distributed quantum computing:
*   Distributed CNOT (derived from a class that allows distribution of any 2-qubits controlled gate)
*   Quantum CNOT gate teleportation

Quantum Machine Learning
*   2D Q-Means Clustering

Optimization
*   QAOA (quantum approximate optimization algorithm) on Ising Model
*   VQE (variational quantum eigensolver)

First, a brief description about usage, if somebody wants to use it in a separate project. The minimal thing to do is to add to include directories the Eigen directory and this project directory (the one that contains the source files). That's about it. Here is some minimal code that can be extended to execute a quantum algorithm:

```
#include <iostream>

#include "QubitRegister.h"

int main()
{
    QC::QubitRegister reg(3);
    QC::Gates::PauliXGate x;
    reg.ApplyGate(x, 0);
    unsigned int m = reg.MeasureAll();

    std::cout << "Result should be 1, it's: " << m << std::endl;

    return 0;
}
```
This is the statevector simulator, which is used in all examples so far. 

The project also has a Matrix Product States (the simplest Tensor Network, the 1D variant that has qubit sites on a chain, see the TEBD project for more details) simulator implementation.
There is not much going on in the project with it yet, but some tests I added pass.

Also an implementation of a simulator able to execute Clifford gates - using the stabilizer formalism - is provided. Again, for this there is not much going on yet, only some tests that pass.

Both can be improved (and hopefully will be) by using multithreading.

You should compile it with c++ 17 (at least).
Of course, if one wants to use some other stuff in there, more headers might need to be included.
The 'test' cpp files are intended for both tests and examples, not necessarily for usage in some other app, including the headers should usually be enough.

I would recommend trying to use the cmake file for compiling on linux or mac (or something similar), turning on openmp increases the speed a lot (as in several times faster, depending on the processor) and on top of that turning on AVX2 (not available for apple 'metal') reduces the execution time again (something like 35 seconds instead of over 50 on my computer for all examples/tests).

The most important code to look into are in `QuantumGate.h` (namespace `QC::Gates`), where the quantum gates are implemented. The code is quite straightforward, a little bit more complex being the implementation of `getOperatorMatrix` (obsoleted now, it's not used anymore for 1, 2 or 3 qubits gates, creation and multiplication with the big operator matrices for those gates is optimized out now), which extends the operator matrix from the matrix that involves only the qubits on which the operator is applied, to all qubits from the register. 
Another header to look into is `QubitRegister.h`. The code is again straightforward, a little bit more complex being the one that does the measurement on a part of the register only. Also since the usage of the big operator matrices obtained by tensor product for 1, 2, 3 qubits gates is obsoleted, the code from `QubitRegister::ApplyGate` also got more complex so it might not be so easy to grasp.
One used everywhere is the `QuantumAlgorithm.h` although that one is very simple, mostly a proxy for the qubit register.
Another source file worth looking into is `Utils.h`, since `BellState` and `MeasurementBasis` are used in several algorithms.

For general operators/gates applied on the register, see the `NQubitsQuantumGate` and `NQubitsControlledQuantumGate` "subalgorithm" classes (for implementations with a big matrix) or `NControlledGateWithAncilla.h` and `NControlledNotWithAncilla.h` (for implementations with simple gates and ancilla qubits).

To see various algorithms in action, look into the `Test` files.

I intend to come from time to time to add more algorithms, but for now I think I have more here than I intended when I started the project, so I'll leave it like this for a while (except maybe some improvements of the existing code and bug fixes if necessary).

### Tools

The project compiles on Windows with Visual Studio 2026 (the code can be compiled with older versions starting with VS 2015, but it's currently maintained with VS 2026 and C++ 17 or higher).

### Required libraries

The program requires the typical VC++ runtime libraries.
Dealing with matrices is done with the help of [Eigen](https://eigen.tuxfamily.org/).
The quantum fourier transform is checked against [FFTW](http://fftw.org/) for Schrodinger quantum simulation. I guess using this library could be avoided, I actually provided several methods of solving the equation in the code, maybe I'll switch to some other one in the tests in the future. For now the check is done with FFTW.

### Bibliograpy

> [!TIP]
> I mentioned some things on the corresponding blog page, so please [check it out](https://compphys.go.ro/quantum-computing-simulator/), but I'll also refer here the tutorials papers I mentioned there and the code for them, which I also [provided in another repository](https://github.com/InvictusWingsSRL/QiskitTutorials).
> If you want to learn qiskit, I higly recommend those.
