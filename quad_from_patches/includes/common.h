#ifndef QUADFROMPATCHES_COMMON_H
#define QUADFROMPATCHES_COMMON_H

#include <vector>

#define DEFAULTILPMETHOD qfp::ILPMethod::LEASTSQUARES
#define DEFAULTALPHA 0.5
#define DEFAULTTIMELIMIT 60 //1 minute
#define DEFAULTGAPLIMIT 0.1
#define DEFAULTCHARTSMOOTHINGITERATIONS 5
#define DEFAULTQUADRANGULATIONSMOOTHINGITERATIONS 5
#define DEFAULTREGULARITYFORNONQUADRILATERAL true
#define DEFAULTNONQUADRILATERALSIMILARITYFACTOR 1.2

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
    bool regularityForNonQuadrilaterals;
    double nonQuadrilateralSimilarityFactor;

    Parameters() {
        ilpMethod = DEFAULTILPMETHOD;
        alpha = DEFAULTALPHA;
        chartSmoothingIterations = DEFAULTCHARTSMOOTHINGITERATIONS;
        quadrangulationSmoothingIterations = DEFAULTQUADRANGULATIONSMOOTHINGITERATIONS;
        timeLimit = DEFAULTTIMELIMIT;
        gapLimit = DEFAULTGAPLIMIT;
        regularityForNonQuadrilaterals = DEFAULTREGULARITYFORNONQUADRILATERAL;
        nonQuadrilateralSimilarityFactor = DEFAULTNONQUADRILATERALSIMILARITYFACTOR;
    }
};

}

#endif // QUADFROMPATCHES_COMMON_H
