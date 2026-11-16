# QCSim
Quantum computation simulator

Has a statevector simulator (better than the 'naive' matrix multiplication kind), a Matrix Product State simulator, a Clifford gates simulator using the stabilizer formalism, a Pauli propagation simulator and a Path Integral simulator. Those are implemented 'from scratch' with no help from AI whatsoever (maybe some of the tests were generated) - later edit: see the Note following.

There is also an 'extended stabilizer' simulator based on the frames formalism. Its single-frame implementation supports Clifford gates and arbitrary single-qubit X, Y and Z rotations. Multiframe coalescing is still work in progress (postponed indefinitely). 
I started to implement it with no AI help and got it working up to executing Clifford gates, then got bored and resumed the work on it with Chat GPT (5.6 Sol Ultra). 
It worked from the first try, but I had to look over the code, spot issues and prompt for refactoring (there was quite a bit of redundancy)... in the end I've to code almost twice as fast as the initial one.

There is also a Matrix Product Operator simulator (which is a compressed Density Matrix in a similar way as the MPS is a compressed statevector). This one is vibed into existence using Claude Opus 4.8 and GPT-5.5 and although I tried to explain the AI quite detailed my intentions and I used the analogy with MPS, explaining the differences as well, I'm not sure if everything is ok as I didn't check all the code yet. The generated tests do pass, so there are some things that work ok, but I wouldn't still rely on it. We might use it in Maestro and until then I'll check it more carefully.

I also added a DensityMatrix simulator and although I added it with AI as well, it should be good as it's based on the statevector 'calculator' and pretty straightforward (although Claude did some mistakes affecting performance which I had to correct by additional prompting).

> [!NOTE]
> The blog describes the Matrix Product State implementation in the Vidal decomposition. The MPO was originally implemented using mirroring the MPS implementation. The current implementation for both is changed to [Hastings' method](https://arxiv.org/abs/0903.3253).
> This way both precision is increased (some issues caused by the divisions with the singular values which amplify errors were removed) and the speed is increased (more for small bond dimensions, less - in percentage - for big bond dimensions).
> I implemented it with the help of Opus 5.5, reviewed the code with Astra.
>
> The statevector - and the Density Matrix, since it uses the same 'calculator', but this one benefits less - are now optimized better, all gates are applied 'in place' - so the statevector uses half the memory it used before - and I tried to apply open mp and avx2 optimizations that should work on both windows and linux.
> More simulators might be optimized with the AI help in the future. Density Matrix now uses a zero-copy vector view of its matrix storage, so the gate kernels operate on its entries as a 2n-qubit vector.
>
> Other simulators optimized with AI in later commits: Path integral, Pauli propagation, Clifford simulator, Density Matrix.

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
    unsigned int m = static_cast<unsigned int>(reg.MeasureAll());

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

### Gate support and qubit ordering

`QuantumGate.h` includes `RZZGate<>`, `RXXGate<>`, `RYYGate<>`, `RZXGate<>`, `RXYGate<>`,
`XXPlusYYGate<>`, `XXMinusYYGate<>`, `ECRGate<>`, `DCXGate<>`, `CSGate<>`, and `CSDGGate<>`.
The existing `iSwapGate<>` and `iSwapDagGate<>` remain available. Interaction rotations use
`exp(-i theta P / 2)`. The first qubit passed to `ApplyGate` is local matrix bit 0;
for example, `RZXGate<>(theta)` acts with Z on the first qubit and X on the second.
`RXYGate<>` similarly acts with X on the first and Y on the second. The XX+YY and XX-YY
gates accept `(theta, beta = 0)` and mix the odd- and even-parity sectors, respectively.
The asymmetric conventions follow the IBM definitions for
[RZX](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.circuit.library.RZXGate),
[XX+YY](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.circuit.library.XXPlusYYGate),
[XX-YY](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.circuit.library.XXMinusYYGate), and
[ECR](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.circuit.library.ECRGate).

Every single-qubit gate, including a caller-supplied 2 by 2 matrix, can be wrapped with
`ControlledGate<>` or `DoublyControlledGate<>`. These retain the original matrix's global
phase, which becomes a relative phase when controls are added. Named controlled and CC
variants are also provided, including `ControlledTGate<>`, `CCYGate<>`, `CCHGate<>`,
`CCSGate<>`, `CCSDGGate<>`, `CCTGate<>`, `CCTDGGate<>`, `CCSXGate<>`, `CCSXDGGate<>`,
`CCRxGate<>`, `CCRyGate<>`, `CCRzGate<>`, `CCPhaseShiftGate<>`, and `CCUGate<>`.

```cpp
QC::Gates::UGate<> u(0.3, 0.2, -0.1);
QC::Gates::ControlledGate<> cu(u);
QC::Gates::DoublyControlledGate<> ccu(u);
reg.ApplyGate(cu, 0, 4);       // target, control
reg.ApplyGate(ccu, 0, 4, 2);  // target, control1, control2
```

Statevector, density matrix, MPS, MPO, and path-integral simulators accept these local
gate matrices. MPS and MPO take local matrices on at most three qubits (and
compact multi-controlled, diagonal and multiplexed gates on any number of qubits, see
below). For a
three-qubit gate, the median physical
site stays fixed while the outer sites are routed to its neighbours. Two SVDs split
the updated tensor without dividing by Schmidt values. MPO applies U and its complex
conjugate separately to the ket and bra indices. Existing bond limits, truncation
modes, and double-precision arithmetic still apply. On a chain whose canonical form
earlier truncations have left approximate, a unitary three-qubit gate does what a
two-qubit gate does: it weighs its two cuts by the Schmidt values to its left and
stays on the canonical path. `setThreeQubitEnvironmentRefresh(true)` (MPS) asks for
the optimal cuts instead: the outer environments of the three sites are first brought
up to date with QR sweeps, already orthonormal ones are reused, and the chain keeps
explicit environments from there on. That costs the sweeps and makes later updates
slower (measured 10-35% on three-qubit circuits at bond dimensions 32-128) for a
fidelity gain in the third to fourth digit. After a Kraus operator, and for
non-unitary three-qubit operators, the refreshed environments are always used. Extreme
finite gate scales are handled with power-of-two rescaling before contraction.
If routing or contraction fails during a non-unitary three-site MPS operation, the
tensors, cached environments, and logical qubit mapping are restored. Routing
notifications are delivered after a successful commit. MPO also accepts three-qubit
operators through `ApplyOperator`, `ApplyOperators`,
and `ApplyOperatorAndNormalize`; its Kraus-channel interface remains limited to two qubits.

Pauli propagation and extended stabilizer expose generic `ApplyGate` entry points for
one- to three-qubit gates as well as named interaction methods. Their cost depends on
Pauli-term or stabilizer-component growth. Both backends cache compiled matrix operations
by matrix contents and recognize the exact matrices of supported canonical Clifford gates.
Named interaction methods additionally retain analytic structure and Clifford factors.
A compact gate (diagonal, multi-controlled, multiplexed) on at most three qubits is
accepted by both, and by the Clifford simulator when it is a Clifford, through its dense
matrix on its own qubits; on more qubits they reject it with `std::invalid_argument`
instead of building a `2^k` matrix. Pauli propagation records a gate object that has
a native kernel (a Pauli, H, S, a square root of X, a rotation about one axis, CX, CY,
CZ, a swap) as that operation, so `GetOperations()` shows the native operation where the
action is the same; `ApplyControlledGate` and `ApplyDoublyControlledGate` find their
table from the 2 by 2 matrix without building the controlled one.
Typed one-qubit, controlled, CC, and interaction gate objects retain constructor
parameters for analytic Pauli transfers, together with the coefficients the
constructor actually stored. The generic bridge checks the entire current matrix
against those stored coefficients before using the parameters, so the check does not
depend on how a compiler rounds a second evaluation (fused multiply-add, constant
folding), and keys its bounded cache by both matrix contents and the descriptor. Mutated gate objects fall back to matrix action;
constructing an `AppliedGate` from a raw matrix does not assert unitarity. Verified
unitary transfers preserve an untouched identity sector exactly, avoiding growth from
round-off coefficients. Signed zeros share the same matrix cache key.
Euler gates form their entries from separate unit-modulus phase factors, avoiding
unstable angle sums, overflow, and loss of small phases when large angles cancel;
a phase times its exact conjugate is taken as real, so that cancellation does not
depend on whether a compiler fuses a multiply-add.
`AppliedGate<>(typedGate, q0, q1, q2)` retains constructor metadata, including through
`.Adjoint()`. Raw-matrix construction remains available with its supplied coefficients;
metadata is checked against the entire current matrix before analytic execution.
Unrecognized raw matrices are not snapped to unitaries or Cliffords by a tolerance
test, and no entry is rescaled or pinned. What is removed is only what rounding alone
can account for. In Pauli propagation a coefficient of a raw matrix's transfer is
dropped when, and only when, it does not exceed the rounding its own sum can carry: the
rounding of the products and additions that formed it, plus half a unit in the last
place of each stored matrix entry. A coefficient with no cancellation behind it (an
entry of `1e-200`, say) is far above that bound and always stays, and a matrix that is
not unitary keeps its actual action. Without this, the stored entries of a unitary given
as a plain matrix leave coefficients of a few `1e-17` that make every Pauli term branch,
including terms the gate does not touch. A plain matrix that is exactly `a I + b P` for
one Pauli string P that moves a bit (a rotation about a Pauli product: Rx, Ry, RXX,
RYY, RZX, ...) gets its table from `a` and `b` directly, under the same rule. The extended stabilizer applies the same idea
to the amplitudes a generic matrix produces: a real or imaginary part within the
rounding bound of its own sum of products is the residue of a cancellation and is
cleared, each part separately, so components that cancel exactly in exact arithmetic
leave the state instead of staying as ghosts.
`ApplyGate(gate, q0, q1, q2)` uses local low-to-high bit order in both.
The extended stabilizer treats a generic non-unitary matrix as a state filter: its
output is normalized, and a zero output is rejected without changing the state.
In exact mode, verified local unitaries share one component snapshot over a bounded
history of up to 32 gates. Exact descriptor inverses or exact matrix adjoints across
typed and raw gate representations
restore and replay that history without an amplitude threshold, preserving arbitrarily
small input coefficients. This handles nested compute/uncompute, reversed typed
adjoints, exact cross-family aliases such as phase/T and U/Ry, interchangeable CC
controls and Fredkin targets, and disjoint intervening gates. Disjoint
basis-only Clifford operations (such as H, S, and CX) preserve the history. Other native
rotations, overlapping native operations, filters, measurement, reset, restore, and
policy changes invalidate it. Saved states remain independent. Cached adjoint hashes
reject mismatches cheaply; a hash match still requires a full exact matrix comparison.
Crossing the history limit starts a new snapshot, so cancellation outside the retained
window can still leave floating-point residuals: amplitudes that carry rounding
inherited from earlier gates are above the bound of the gate that would cancel them and
are kept. No tolerance-based pruning is introduced.
After prior approximation, a generic filter conservatively raises the accumulated
trace-distance bound to one, as postselection can amplify discarded components.
Verified unitary constructors preserve the existing approximation bound, and so does a
raw matrix that is unitary within `1e-12`; that test decides only this bookkeeping,
never a coefficient.
Pauli propagation retains the unnormalized operator action, including scalar magnitude.
Outcome probabilities assume a trace-preserving circuit: after a matrix gate that is
not unitary (within `1e-12`), `Measure`, `Probability(outcome)`, `Sample` and
`SampleCounts` throw `std::domain_error` instead of returning numbers that are neither
the normalized nor the unnormalized distribution; expectation values stay available.
Pauli's named controlled helpers retain `(control, target)` and `(control1, control2,
target)` ordering; extended-stabilizer helpers use `(target, control)` and `(target,
control1, control2)`. The Clifford simulator supports ECR, DCX, iSWAP, Pauli-product
rotations, and XX+YY/XX-YY at supported Clifford parameter points; it rejects other angles.
Its `ApplyRx`, `ApplyRy`, `ApplyRz`, and `ApplyPhase` methods accept exactly encoded
multiples of pi/2 with at most `2^20` quarter turns in magnitude. Larger angles use
ordinary trigonometry in general gate constructors and are rejected by the Clifford
angle APIs: coarse floating-point spacing must not turn arbitrary large angles into
Cliffords. Its generic `ApplyGate` recognizes supported canonical Clifford matrix
contents, including matrix-backed `AppliedGate` copies, and rejects unrecognized matrices
before changing the state. Recognition includes gate adjoints, all 24 canonical
one-qubit Clifford actions, and canonical one- to three-qubit matrices verified by
exact Pauli intertwining and discrete coefficient checks, including the balanced
real/imaginary coefficients of embedded SX and SX-dagger gates. Both existing
`1/sqrt(d)` and `sqrt(1/d)` coefficient encodings are recognized, with one encoding
required consistently throughout a matrix. Constructor coefficients remain unchanged.
A scaled matrix whose entries are identical to a recognized encoding necessarily has
the same classification; other scalings and mixed encodings are not rounded into it.
Synthesized circuits use
bounded caches with shared ownership. This recognition does not round nearly Clifford
matrices or arbitrary floating-point compositions onto Clifford operations.
Adding controls generally makes a Clifford gate non-Clifford, so a controlled or CC
wrapper does not imply Clifford-backend support.

### Structured gates on more qubits

`QubitRegister<>` and `DensityMatrix<>` share compact APIs for diagonals, multi-controlled
single-qubit operators, and uniformly controlled rotations. They do not construct a
full-register gate matrix or require ancillas.

```cpp
reg.ApplyDiagonalGate(std::vector<std::complex<double>>{1., 1., 1., -1.}, {4, 1});
reg.ApplyMCX(0, {1, 2, 3, 4});
reg.ApplyMCZ(0, {1, 2, 3, 4});
reg.ApplyMCP(0.4, 0, {1, 2, 3, 4});
reg.ApplyMCRX(0.3, 0, {1, 2, 3, 4});
reg.ApplyMultiControlledGate(u, 0, {1, 2, 3, 4});
reg.ApplyMultiplexedRY({0.1, 0.2, 0.3, 0.4}, 0, {4, 1});
```

For diagonals, `qubits[i]` selects bit i of the diagonal index. For multiplexers,
`controls[i]` selects bit i of the angle-table index; there must be exactly
`2^controls.size()` angles. Multi-controlled gates activate when all controls are one
by default; an optional control-state integer selects an arbitrary bit pattern using
the same ordering. Controls must be distinct and cannot include the target. Invalid
indices, lengths, or non-finite coefficients are rejected before applying an operation.
Statevector `Compute`/`Uncompute` records these operations compactly and preserves their
order relative to ordinary gates. The table handed to `ApplyDiagonalGate` or to a
multiplexed rotation is applied where it is: nothing of its size is copied unless the
register is recording its gates, and a recording owns its copy.
`ApplyStructuredGateView` takes a `StructuredGateView` (the fields of a
`StructuredGate` with the tables as pointer and size) for callers that keep their
tables elsewhere. A full-width diagonal (one factor per basis state, as a phase oracle
or a cost layer) is a single pass over the state. When a large table is not in the
order of the state bits (the qubits listed from the most significant one, for
instance), walking the state in order would read the table at scattered places; such a
table is taken through a small buffer instead, filled along the table and applied along
the state, which keeps both in order. That is 2 to 5 times faster than the direct walk
from tables of 2^18 entries on (2^20 when several threads share the work), and is not
used below those sizes or for tables that are already in order.

`AppliedGate<>` can store these operations without expanding a dense matrix. Mix them
with ordinary gates in `ApplyGates`, or pass the circuit to the path simulator's
`SetCircuit`:

```cpp
using Applied = QC::Gates::AppliedGate<>;
std::vector<Applied> circuit{
    Applied(QC::Gates::HadamardGate<>().getRawOperatorMatrix(), 0),
    Applied::MultiControlled(u, 0, {1, 2, 3, 4}),
    Applied::Diagonal({1., 1., 1., -1.}, {4, 1}),
    Applied::MultiplexedRotation(QC::Gates::RotationAxis::Y,
                                {0.1, 0.2, 0.3, 0.4}, 0, {4, 1})
};
reg.ApplyGates(circuit);
```

`getQubits()` and `setQubits()` access the complete ordered qubit list (target first,
then controls, or diagonal table order). The tables of a compact gate are validated
once, when the gate is made, and are immutable afterwards: copies, `Adjoint()`,
`Transposed()` and gates with remapped qubits all share them, together with the values
derived from them (fan-out, rotation coefficients), so none of these copies or rescans a
table. `setQubits()` checks the complete list. `setQubit1/2/3` replace one qubit
without checking, as they do for a dense gate, so that a remap may pass through a list
with a repeated qubit; an index beyond the arity is stored and ignored. A gate left with
a repeated or out-of-range qubit is rejected where it is applied, and applying a compact
`AppliedGate` checks only its qubits, not its table entries again. The statevector and
density-matrix kernels read the shared tables of such a gate under its own qubits: a
remapped, adjoint or transposed gate is applied without a copy of its table (an adjoint
diagonal is conjugated as it is read), `Uncompute` replays the adjoints the same way,
and the cosines and sines of a multiplexed rotation are computed at its first
application and kept with the tables.
`getStructuredGate()` exposes the descriptor; for a remapped, adjoint or transposed gate
it is built on the first call (one copy of the table) and kept, while
`getCompactPayload()` gives the shared tables and the gate's own qubits without
building anything. Compact payloads require a dynamically sized `MatrixClass`.
`getRawOperatorMatrix()` rejects compact gates rather than
allocating an exponentially large matrix. Use the one-argument `ApplyGate(applied)`
to execute its stored qubit mapping.

MPS and MPO share the same compact entry points (`ApplyMCX`, `ApplyMCZ`, `ApplyMCP`,
`ApplyMCRX/Y/Z`, `ApplyMCU`, `ApplyMultiControlledGate`, `ApplyDiagonalGate`, the
multiplexed rotations) and take compact `AppliedGate`s through the one-argument
`ApplyGate`. A multi-controlled one-qubit gate is the identity plus the product of the
control projectors times `U - 1`: an operator of bond dimension two along the chain,
whatever the number of controls and wherever they sit. It is applied in place, with no
routing swaps, no ancillas and no `2^k` matrix: the tensors between the first and the
last of its qubits take the direct sum of the two terms, and one sweep cuts the bonds
again with the configured limits. MPS requires a unitary target there; MPO applies
`O rho O^dagger` for any 2 by 2 target. A diagonal or multiplexed table on at most
three qubits is expanded to its matrix and applied as a local gate. On more qubits it
is applied in place as well, as the chain operator of its table: the table is
written as a product of small matrices, one pair per qubit (a tensor train), whose
sizes say how strongly the gate couples the two sides of each cut. A product of
one-qubit phases has rank one, a sum of couplings between neighbouring table bits rank
two, a multiplexer whose angle is a sum of one term per control rank two; a table
without structure reaches `2^(k/2)`. The tensors between the first and the last qubit
of the gate take those matrices, which multiplies their bonds by the rank (by its
square on MPO, where the operator acts on both sides), and one sweep cuts the bonds
again with the configured limits. The cost therefore follows the rank, not the number
of qubits; tables above rank 64 (MPS) or 16 (MPO) are rejected before anything
changes. Ranks are numerical ranks: a singular value is dropped when it is below the
rounding of the table itself. MPS requires a unitary diagonal there; MPO applies
`O rho O^dagger` for any. When a gate is known as a few two-qubit gates on
neighbouring qubits, applying those is cheaper (measured about 2 times on 12 qubits);
the chain operator is for tables that are only known as tables, and for qubits far
apart, where it needs no routing.

An arbitrary k-qubit diagonal needs only `2^k` coefficients. Multi-controlled kernels
visit the active target pairs, with work proportional to `2^(n-c)` for n statevector
qubits and c controls. Multiplexers precompute rotation coefficients and visit each
target pair once. Density-matrix updates apply the corresponding ket/bra actions;
diagonals use a single pass over the density matrix.

Path-integral simulation can apply a diagonal without adding paths, and a
multi-controlled or multiplexed single-qubit gate has at most two outgoing transitions
per path. Arbitrary diagonal and multiplexed tables can require exponentially large
Pauli expansions.
Wide compact gates of every kind are implemented for statevector, density matrix, and
path-integral simulation; MPS and MPO add wide multi-controlled gates and wide tables of
low rank along the chain. Path execution
includes forward/backward endpoint queries, caching, and circuit repartitioning. Its
cost still depends on the number of live paths. In the path integral a compact gate
costs what it does to the paths, not the size of its table: a step reads the entries its
paths select, compact gates on one to three qubits run on the dense transition tables,
and they use the compiled-gate cache like dense gates. Compact and dense gates throw the
same exception types there (`std::out_of_range` for a qubit outside the register,
`std::invalid_argument` for repeated qubits, `std::length_error` beyond the supported
register width).

To see various algorithms in action, look into the `Test` files.

I intend to come from time to time to add more algorithms, but for now I think I have more here than I intended when I started the project, so I'll leave it like this for a while (except maybe some improvements of the existing code and bug fixes if necessary).

> [!NOTE]
> Qcsim simulators are used in Qoro Quantum's [Maestro](https://github.com/QoroQuantum/maestro), where also [qiskit aer](https://github.com/Qiskit/qiskit-aer) and [cuQuantum](https://github.com/NVIDIA/cuQuantum) based cuda simulators are used, so the simulators here are tested against those.
> In there you may find extensions of qcsim simulators, like the p-blocked simulator (I called it 'composite') - which is useful for simulating distributed quantum computation - or a 'proof of concept' tensor network simulator (there is a serious one implemented with cuQuantum that is also used, but that one is close source).

### Tools

The project compiles on Windows with Visual Studio 2026 (the code can be compiled with older versions starting with VS 2015, but it's currently maintained with VS 2026 and C++ 17 or higher).

### Required libraries

The program requires the typical VC++ runtime libraries.
Dealing with matrices is done with the help of [Eigen](https://eigen.tuxfamily.org/).
The quantum fourier transform is checked against [FFTW](http://fftw.org/) for Schrodinger quantum simulation. I guess using this library could be avoided, I actually provided several methods of solving the equation in the code, maybe I'll switch to some other one in the tests in the future. For now the check is done with FFTW.

### Bibliography

> [!TIP]
> I mentioned some things on the corresponding blog page, so please [check it out](https://compphys.go.ro/quantum-computing-simulator/), but I'll also refer here the tutorials papers I mentioned there and the code for them, which I also [provided in another repository](https://github.com/InvictusWingsSRL/QiskitTutorials).
> If you want to learn qiskit, I highly recommend those.
