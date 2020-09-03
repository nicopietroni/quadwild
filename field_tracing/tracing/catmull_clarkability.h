#ifndef CATMULL_CLARKABILITY
#define CATMULL_CLARKABILITY

#include <vector>

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2)
{
    ScalarType AvgL=(L1+L2)/2;
    //ScalarType CC=pow((L1-L2),2)/pow(AvgL,2);
    ScalarType CC=fabs(L1-L2)/AvgL;
    assert(CC>=0);
    //CC=std::min<ScalarType>(CC,1);
    return CC;
}

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2,
                     ScalarType L3)
{
    ScalarType AvgL=(L1+L2+L3)/3;
    ScalarType CC=fabs(L1-(L2+L3))/AvgL;
    assert(CC>=0);
    //CC=std::min<ScalarType>(CC,1);
    return CC;
}

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2,
                     ScalarType L3,
                     ScalarType L4,
                     ScalarType L5)
{
    ScalarType AvgL=(L1+L2+L3+L4+L5)/5;
    ScalarType CC=fabs(L1+L2+L3-(L4+L5))/AvgL;
    assert(CC>=0);
    //CC=std::min<ScalarType>(CC,1);
    return CC;
}


template <class ScalarType>
ScalarType CatmullClarkability(ScalarType L1,
                               ScalarType L2,
                               ScalarType L3)
{
    ScalarType CC0=OrderedCC(L1,L2,L3);
    ScalarType CC1=OrderedCC(L2,L3,L1);
    ScalarType CC2=OrderedCC(L3,L1,L2);
    return (std::max(CC0,std::max(CC1,CC2)));
}

template <class ScalarType>
ScalarType CatmullClarkability(ScalarType L1,
                               ScalarType L2,
                               ScalarType L3,
                               ScalarType L4)
{
    ScalarType CC0=OrderedCC(L1,L3);
    ScalarType CC1=OrderedCC(L2,L4);
    return (std::max(CC0,CC1));
}

template <class ScalarType>
ScalarType CatmullClarkability(ScalarType L1,
                               ScalarType L2,
                               ScalarType L3,
                               ScalarType L4,
                               ScalarType L5)
{
    ScalarType CC0=OrderedCC(L1,L2,L3,L4,L5);
    ScalarType CC1=OrderedCC(L2,L3,L4,L5,L1);
    ScalarType CC2=OrderedCC(L3,L4,L5,L1,L2);
    ScalarType CC3=OrderedCC(L4,L5,L1,L2,L3);
    ScalarType CC4=OrderedCC(L5,L1,L2,L3,L4);
    return (std::max(CC0,std::max(CC1,std::max(CC2,std::max(CC3,CC4)))));
}

template <class ScalarType>
ScalarType CatmullClarkability(ScalarType L1,
                               ScalarType L2,
                               ScalarType L3,
                               ScalarType L4,
                               ScalarType L5,
                               ScalarType L6)
{
    ScalarType CC0=OrderedCC(L1,L3,L5);
    ScalarType CC1=OrderedCC(L2,L4,L6);
    return (std::max(CC0,CC1));
}

template <class ScalarType>
ScalarType CatmullClarkability(const std::vector<ScalarType> &EdgeL)
{
//    assert(EdgeL.size()>=3);
//    assert(EdgeL.size()<=6);

    if (EdgeL.size()==3)
        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2]));
    if (EdgeL.size()==4)
        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3]));
    if (EdgeL.size()==5)
        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4]));
    if (EdgeL.size()==6)
        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4],EdgeL[5]));
    return (std::numeric_limits<ScalarType>::max());
}

#endif
