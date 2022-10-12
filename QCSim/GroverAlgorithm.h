#pragma once

#include "QuantumAlgorithm.h"
#include "Oracle.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace Grover {

    template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class GroverAlgorithm :
        public QC::QuantumAlgorithm<VectorClass, MatrixClass>
    {
    public:
        GroverAlgorithm(unsigned int N = 3, int addseed = 0) 
            : QC::QuantumAlgorithm<VectorClass, MatrixClass>(N, addseed)
        {
            J j;
            JOp = j.getOperatorMatrix(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits());
            
            setCorrectQuestionState(0);
        }

        void setCorrectQuestionState(unsigned int state)
        {
            Oracle<MatrixClass> o;
            o.setCorrectQuestionState(state);
            OracleOp = o.getOperatorMatrix(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits());
        }

        unsigned int Execute() override
        {
            Init();

            const unsigned int repeatNo = static_cast<unsigned int>(round(M_PI / 4. * sqrt(QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrBasisStates())));

            for (unsigned int i = 0; i < repeatNo; ++i)
            {
                QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyOperatorMatrix(OracleOp);
                ApplyDiffusionOperator();
            }

            return QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.Measure();
        }
 
    protected:
        void Init()
        {
            //reg.setToBasisState(0);
            //ApplyHadamardOnAllQubits();
            QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.setToEqualSuperposition(); // the same thing as commented above
        }

        void ApplyHadamardOnAllQubits()
        {
            for (unsigned int i = 0; i < QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.getNrQubits(); ++i)
                QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyGate(hadamard, i);
        }

        void ApplyDiffusionOperator()
        {
            ApplyHadamardOnAllQubits();
            QC::QuantumAlgorithm<VectorClass, MatrixClass>::reg.ApplyOperatorMatrix(JOp);
            ApplyHadamardOnAllQubits();
        }

        MatrixClass OracleOp;
        MatrixClass JOp;

        QC::HadamardGate<MatrixClass> hadamard;
    };

}
