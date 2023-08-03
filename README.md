# QCSim
Quantum computation simulator

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6a193db170ab432596079c530fc75c77)](https://www.codacy.com/gh/aromanro/QCSim/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=aromanro/QCSim&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/aromanro/qcsim/badge)](https://www.codefactor.io/repository/github/aromanro/qcsim)

Implemented algorithms:

*   Grover
*   Deutsch-Jozsa
*   Simon
*   Quantum Fourier Transform
*   Phase estimation
*   Shor
*   Bernstein–Vazirani
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
*   Distributed CNOT
*   Quantum CNOT gate teleportation

The most important code to look into are in `QuantumGate.h` (namespace `QC::Gates`), where the quantum gates are implemented. The code is quite straightforward, a little bit more complex being the implementation of `getOperatorMatrix` (obsoleted, it's not used anymore for 1, 2 or 3 qubits gates, creation and multiplication with the big operator matrices for those gates is optimized out now), which extends the operator matrix from the matrix that involves only the qubits on which the operator is applied, to all qubits from the register. 
Another header to look into is `QubitRegister.h`. The code is again straightforward, a little bit more complex being the one that does the measurement on a part of the register only. Also since the usage of the big operator matrices obtained by tensor product for 1, 2, 3 qubits gates is obsoleted, the code from `QubitRegister::ApplyGate` also got more complex so it might not be so easy to grasp.
One used everywhere is the `QuantumAlgorithm.h` although that one is very simple, mostly a proxy for the qubit register.
Another source file worth looking into is `Utils.h`, since `BellState` and `MeasurementBasis` are used in several algorithms.

For general operators/gates applied on the register, see the `NQubitsQuantumGate` and `NQubitsControlledQuantumGate` "subalgorithm" classes.

To see various algorithms in action, look into the `Test` files.

I intend to come from time to time to add more algorithms, but for now I think I have more here than I intended when I started the project, so I'll leave it like this for a while (except maybe some improvements of the existing code and bug fixes if necessary).

### Required libraries

Dealing with matrices is done with the help of [Eigen](https://eigen.tuxfamily.org/).
The quantum fourier transform is checked against [FFTW](http://fftw.org/) for Schrodinger quantum simulation. I guess using this library could be avoided, I actually provided several methods of solving the equation in the code, maybe I'll switch to some other one in the tests in the future. For now the check is done with FFTW.