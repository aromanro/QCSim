#include "GroverAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Grover
{
    unsigned int GroverAlgorithm::Execute()
    {
        Init();

        const unsigned int repeatNo = static_cast<unsigned int>(round(M_PI / 4. * sqrt(reg.getNrBasisStates())));

        for (unsigned int i = 0; i < repeatNo; ++i)
        {
            reg.ApplyOperatorMatrix(OracleOp);
            ApplyDiffusionOperator();
        }

        return reg.Measure();
    }
}