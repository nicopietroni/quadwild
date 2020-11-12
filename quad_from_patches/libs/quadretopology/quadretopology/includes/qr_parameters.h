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
#define DEFAULTREGULARITYQUADRILATERALS true
#define DEFAULTREGULARITYNONQUADRILATERALS true
#define DEFAULTREGULARITYNONQUADRILATERALSWEIGHT 0.9
#define DEFAULTALIGNSINGULARITIES false
#define DEFAULTALIGNSINGULARITIESWEIGHT 0.01
#define DEFAULTREPEATLOSINGCONSTRAINTSITERATIONS 0
#define DEFAULTREPEATLOSINGCONSTRAINTSQUADS false
#define DEFAULTREPEATLOSINGCONSTRAINTSNONQUADS false
#define DEFAULTREPEATLOSINGCONSTRAINTSALIGN false
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
    bool regularityQuadrilaterals;
    bool regularityNonQuadrilaterals;
    double regularityNonQuadrilateralsWeight;
    bool alignSingularities;
    double alignSingularitiesWeight;
    bool repeatLosingConstraintsIterations;
    bool repeatLosingConstraintsQuads;
    bool repeatLosingConstraintsNonQuads;
    bool repeatLosingConstraintsAlign;
    bool feasibilityFix;
    bool hardParityConstraint;

    double timeLimit;
    double gapLimit;
    double minimumGap;
    std::vector<float> callbackTimeLimit;
    std::vector<float> callbackGapLimit;

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
        regularityQuadrilaterals = DEFAULTREGULARITYQUADRILATERALS;
        regularityNonQuadrilaterals = DEFAULTREGULARITYNONQUADRILATERALS;
        regularityNonQuadrilateralsWeight = DEFAULTREGULARITYNONQUADRILATERALSWEIGHT;
        alignSingularities = DEFAULTALIGNSINGULARITIES;
        alignSingularitiesWeight = DEFAULTALIGNSINGULARITIESWEIGHT;
        repeatLosingConstraintsIterations = DEFAULTREPEATLOSINGCONSTRAINTSITERATIONS;
        repeatLosingConstraintsQuads = DEFAULTREPEATLOSINGCONSTRAINTSQUADS;
        repeatLosingConstraintsNonQuads = DEFAULTREPEATLOSINGCONSTRAINTSNONQUADS;
        repeatLosingConstraintsAlign = DEFAULTREPEATLOSINGCONSTRAINTSALIGN;
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
