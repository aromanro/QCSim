#pragma once
#include "QuantumAlgorithm.h"
#include "QuantumGate.h"

namespace Teleportation
{

    template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumTeleportation :
        public QC::QuantumAlgorithm<VectorClass, MatrixClass>
    {
    public:
        QuantumTeleportation(unsigned int N = 3, int addseed = 0)
            : QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed)
        {
        }
    };


    template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QuantumTeleportationRealization : QuantumTeleportation < VectorClass, MatrixClass)
    {
    public:
        QuantumTeleportationRealization(unsigned int N = 3, int addseed = 0)
            : QuantumTeleportation(N, addseed)
        {
        }

    protected:
        // the initial state to be teleported can be set directly or set up from something simple using some gate (hadamard, phase shift)
        // needed to create the Bell states
        QC::HadamardGate hadamard;
        QC::CNOTGate cnot; // also needed on Alice side
        // needed on receiving side:
        QC::ControlledZGate cz;
    };

}
