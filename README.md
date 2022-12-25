# QCSim
Quantum computation simulator

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6a193db170ab432596079c530fc75c77)](https://www.codacy.com/gh/aromanro/QCSim/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=aromanro/QCSim&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/aromanro/qcsim/badge)](https://www.codefactor.io/repository/github/aromanro/qcsim)

Work in progress, but the following seem to work:
*   Grover's algorithm
*   Deutsch-Jozsa algorithm
*   Shor's algorithm
*   Bernstein–Vazirani algorithm
*   Quantum adder (half-adder & full-adder for 1-qubits and full adder for two N-qubit numbers)
*   Quantum teleportation
*   Superdense coding
*   Quantum cryptography: BB84 protocol
*   CHSH inequality violation

Simulation of quantum simulation:

*   Evolution in time with a Hamiltonian given as a sum of Pauli products with real coefficients
*   Evolution of a 1D Gaussian packet in time, solving the 1D time-dependent Schrodinger equation using a Trotter decomposition and quantum Fourier transform

Paradoxes (although some of the above might be considered paradoxes as well):
*   Quantum eraser
*   General Elitzur-Vaidman Bomb tester/interaction free measurement/counterfactual computation
*   Hardy's paradox

Maybe I'll add some more things in the future, like something related with quantum random walks and/or quantum tomography, but for now this will take a break for a while. 
