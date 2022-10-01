#pragma once
#include "QuantumAlgorithm.h"

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

}
