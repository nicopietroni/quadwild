#ifndef QR_PARAMETERS_H
#define QR_PARAMETERS_H

#include <vector>

#define DEFAULTINITIALREMESHING true
#define DEFAULTEDGEFACTOR 1
#define DEFAULTREPROJECT true
#define DEFAULTSPLITCONCAVES false
#define DEFAULTFINALSMOOTHING true

#define DEFAULTILPMETHOD QuadRetopology::ILPMethod::LEASTSQUARES
#define DEFAULTALPHA 0.5
#define DEFAULTISOMETRY true
#define DEFAULTREGULARITYFORQUADRILATERAL true
#define DEFAULTREGULARITYFORNONQUADRILATERAL true
#define DEFAULTREGULARITYNONQUADRILATERWEIGHT 0.9
#define DEFAULTHARDPARITYCONSTRAINT true
#define DEFAULTFEASIBILITYFIX false

#define DEFAULTTIMELIMIT 60 //1 minute
#define DEFAULTGAPLIMIT 0.0 //Optimal
#define DEFAULTMINIMUMGAP 0.2

#define DEFAULTCHARTSMOOTHINGITERATIONS 5
#define DEFAULTQUADRANGULATIONFIXEDSMOOTHINGITERATIONS 5
#define DEFAULTQUADRANGULATIONNONFIXEDSMOOTHINGITERATIONS 5
#define DEFAULTDOUBLETREMOVAL true

#define DEFAULTRESULTSMOOTHINGITERATIONS 5
#define DEFAULTRESULTSMOOTHINGNRING 3
#define DEFAULTRESULTSMOOTHINGLAPLACIANITERATIONS 2
#define DEFAULTRESULTSMOOTHINGLAPLACIANNRING 3

namespace QuadRetopology {

enum ILPMethod { LEASTSQUARES, ABS };

struct Parameters {
    bool initialRemeshing;
    double initialRemeshingEdgeFactor;
    bool reproject;
    bool splitConcaves;
    bool finalSmoothing;
    
    ILPMethod ilpMethod;
    double alpha;
    bool isometry;
    bool regularityForQuadrilaterals;
    bool regularityForNonQuadrilaterals;
    double regularityNonQuadrilateralWeight;
    bool feasibilityFix;
    bool hardParityConstraint;

    double timeLimit;
    double gapLimit;
    double minimumGap;

    int chartSmoothingIterations;
    int quadrangulationFixedSmoothingIterations;
    int quadrangulationNonFixedSmoothingIterations;
    bool doubletRemoval;
    
    int resultSmoothingIterations;
    double resultSmoothingNRing;
    int resultSmoothingLaplacianIterations;
    double resultSmoothingLaplacianNRing;    

    Parameters() {
        initialRemeshing = DEFAULTINITIALREMESHING;
        initialRemeshingEdgeFactor = DEFAULTEDGEFACTOR;
        reproject = DEFAULTREPROJECT;
        splitConcaves = DEFAULTSPLITCONCAVES;
        finalSmoothing = DEFAULTFINALSMOOTHING;
                
        ilpMethod = DEFAULTILPMETHOD;
        alpha = DEFAULTALPHA;

        isometry = DEFAULTISOMETRY;
        regularityForQuadrilaterals = DEFAULTREGULARITYFORQUADRILATERAL;
        regularityForNonQuadrilaterals = DEFAULTREGULARITYFORNONQUADRILATERAL;
        regularityNonQuadrilateralWeight = DEFAULTREGULARITYNONQUADRILATERWEIGHT;
        hardParityConstraint = DEFAULTHARDPARITYCONSTRAINT;
        feasibilityFix = DEFAULTFEASIBILITYFIX;

        timeLimit = DEFAULTTIMELIMIT;
        gapLimit = DEFAULTGAPLIMIT;
        minimumGap = DEFAULTMINIMUMGAP;

        chartSmoothingIterations = DEFAULTCHARTSMOOTHINGITERATIONS;
        quadrangulationFixedSmoothingIterations = DEFAULTQUADRANGULATIONFIXEDSMOOTHINGITERATIONS;
        quadrangulationNonFixedSmoothingIterations = DEFAULTQUADRANGULATIONNONFIXEDSMOOTHINGITERATIONS;
        doubletRemoval = DEFAULTDOUBLETREMOVAL;
        
        resultSmoothingIterations = DEFAULTRESULTSMOOTHINGITERATIONS;
        resultSmoothingNRing = DEFAULTRESULTSMOOTHINGNRING;
        resultSmoothingLaplacianIterations = DEFAULTRESULTSMOOTHINGLAPLACIANITERATIONS;
        resultSmoothingLaplacianNRing = DEFAULTRESULTSMOOTHINGLAPLACIANNRING;
    }
};

}

#endif // QR_PARAMETERS_H
