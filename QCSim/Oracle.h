#pragma once
#include "QuantumGate.h"

namespace Grover {

    class Oracle :
        public QC::QuantumGate
    {
    public:
        void setCorrectQuestionState(unsigned int state)
        {
            correctQuestionState = state;
        }

        Eigen::MatrixXcd getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override;

    protected:
        unsigned int correctQuestionState = 0;
    };

    class J : public QC::QuantumGate
    {
    public:
        Eigen::MatrixXcd getOperatorMatrix(unsigned int nrQubits, unsigned int qubit = 0, unsigned int controllingQubit = 0) const override;
    };

}

