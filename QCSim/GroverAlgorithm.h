#pragma once

#include "QuantumAlgorithm.h"
#include "Oracle.h"

namespace Grover {

    class GroverAlgorithm :
        public QC::QuantumAlgorithm
    {
    public:
        GroverAlgorithm(int N = 3, int addseed = 0) 
            : QC::QuantumAlgorithm(N, addseed)
        {
            J j;
            JOp = j.getOperatorMatrix(reg.getNrQubits());
            
            setCorrectQuestionState(0);
        }

        void setCorrectQuestionState(unsigned int state)
        {
            Oracle o;
            o.setCorrectQuestionState(state);
            OracleOp = o.getOperatorMatrix(reg.getNrQubits());
        }

        unsigned int Execute() override;
 
    protected:
        void Init()
        {
            //reg.setToBasisState(0);
            //ApplyHadamardOnAllQubits();
            reg.setToEqualSuperposition(); // the same thing as commented above
        }

        void ApplyHadamardOnAllQubits()
        {
            for (unsigned int i = 0; i < reg.getNrQubits(); ++i)
                reg.ApplyGate(hadamard, i);
        }

        void ApplyDiffusionOperator()
        {
            ApplyHadamardOnAllQubits();
            reg.ApplyOperatorMatrix(JOp);
            ApplyHadamardOnAllQubits();
        }

        Eigen::MatrixXcd OracleOp;
        Eigen::MatrixXcd JOp;

        QC::HadamardGate hadamard;
    };

}
