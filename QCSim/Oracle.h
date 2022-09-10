#pragma once
#include "QuantumGate.h"

namespace Grover {

    template<class MatrixClass = Eigen::MatrixXcd> class Oracle :
        public QC::QuantumGate<MatrixClass>
    {
    public:
        void setCorrectQuestionState(unsigned int state)
        {
            correctQuestionState = state;
        }

        MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override
        {
            const unsigned int nrBasisStates = 1u << nrQubits;
            MatrixClass extOperatorMat = MatrixClass::Identity(nrBasisStates, nrBasisStates);

            if (correctQuestionState < nrBasisStates)
                extOperatorMat(correctQuestionState, correctQuestionState) = -1;

            return extOperatorMat;
        }

    protected:
        unsigned int correctQuestionState = 0;
    };

    template<class MatrixClass = Eigen::MatrixXcd> class J : public QC::QuantumGate<MatrixClass>
    {
    public:
        MatrixClass getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override
        {
            const unsigned int nrBasisStates = 1u << nrQubits;
            MatrixClass extOperatorMat = MatrixClass::Identity(nrBasisStates, nrBasisStates);

            extOperatorMat(0, 0) = -1;

            return extOperatorMat;
        }
    };

}

