#ifndef QUADFROMPATCHES_COMMON_H
#define QUADFROMPATCHES_COMMON_H

#include <vector>

#define DEFAULTILPMETHOD qfp::ILPMethod::LEASTSQUARES
#define DEFAULTALPHA 0.5

#define DEFAULTISOMETRY true
#define DEFAULTREGULARITYFORQUADRILATERAL true
#define DEFAULTREGULARITYFORNONQUADRILATERAL true
#define DEFAULTregularityNonQuadrilateralWeight 0.9
#define DEFAULTHARDPARITYCONSTRAINT false

#define DEFAULTTIMELIMIT 60 //1 minute
#define DEFAULTGAPLIMIT 0.0 //Optimal
#define DEFAULTMINIMUMGAP 0.2

#define DEFAULTCHARTSMOOTHINGITERATIONS 5
#define DEFAULTQUADRANGULATIONSMOOTHINGITERATIONS 5

namespace qfp {

enum ILPMethod { LEASTSQUARES, ABS };

struct Parameters {
    ILPMethod ilpMethod;
    double alpha;

    bool isometry;
    bool regularityForQuadrilaterals;
    bool regularityForNonQuadrilaterals;
    double regularityNonQuadrilateralWeight;
    bool hardParityConstraint;

    double timeLimit;
    double gapLimit;
    double minimumGap;

    int chartSmoothingIterations;
    int quadrangulationSmoothingIterations;

    Parameters() {
        ilpMethod = DEFAULTILPMETHOD;
        alpha = DEFAULTALPHA;

        isometry = DEFAULTISOMETRY;
        regularityForQuadrilaterals = DEFAULTREGULARITYFORQUADRILATERAL;
        regularityForNonQuadrilaterals = DEFAULTREGULARITYFORNONQUADRILATERAL;
        regularityNonQuadrilateralWeight = DEFAULTregularityNonQuadrilateralWeight;
        hardParityConstraint = DEFAULTHARDPARITYCONSTRAINT;

        timeLimit = DEFAULTTIMELIMIT;
        gapLimit = DEFAULTGAPLIMIT;
        minimumGap = DEFAULTMINIMUMGAP;

        chartSmoothingIterations = DEFAULTCHARTSMOOTHINGITERATIONS;
        quadrangulationSmoothingIterations = DEFAULTQUADRANGULATIONSMOOTHINGITERATIONS;
    }
};

}

#endif // QUADFROMPATCHES_COMMON_H
