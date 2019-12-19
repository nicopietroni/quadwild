#ifndef QUADFROMPATCHES_COMMON_H
#define QUADFROMPATCHES_COMMON_H

#include <vector>

#define DEFAULTILPMETHOD qfp::ILPMethod::LEASTSQUARES
#define DEFAULTALPHA 0.5
#define DEFAULTTIMELIMIT 10 * 60
#define DEFAULTGAPLIMIT 0.1
#define DEFAULTCHARTSMOOTHINGITERATIONS 5
#define DEFAULTQUADRANGULATIONSMOOTHINGITERATIONS 5

namespace qfp {

enum ILPMethod { LEASTSQUARES, ABS };

struct Parameters {
    ILPMethod ilpMethod;
    double alpha;
    double timeLimit;
    double gapLimit;
    bool finalSmoothing;
    int chartSmoothingIterations;
    int quadrangulationSmoothingIterations;

    Parameters() {
        ilpMethod = DEFAULTILPMETHOD;
        alpha = DEFAULTALPHA;
        chartSmoothingIterations = DEFAULTCHARTSMOOTHINGITERATIONS;
        quadrangulationSmoothingIterations = DEFAULTQUADRANGULATIONSMOOTHINGITERATIONS;
        timeLimit = DEFAULTTIMELIMIT;
        gapLimit = DEFAULTGAPLIMIT;
    }
};

}

#endif // QUADFROMPATCHES_COMMON_H
