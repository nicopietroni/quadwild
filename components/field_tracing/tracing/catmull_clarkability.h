/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#ifndef CATMULL_CLARKABILITY
#define CATMULL_CLARKABILITY

#include <vector>

//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2)
//{
//    //ScalarType AvgL=(L1+L2)/2;
//    ScalarType AvgL=std::min(L1,L2);
//    //ScalarType CC=pow((L1-L2),2)/pow(AvgL,2);
//    ScalarType CC=(fabs(L1-L2)-AvgL)/AvgL;
//    CC=std::max(CC,(ScalarType)0);
//    assert(CC>=0);
//    //CC=std::min<ScalarType>(CC,1);
//    return CC;
//}

//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2,
//                     ScalarType L3)
//{
//    //ScalarType AvgL=(L1+L2+L3)/3;
//    ScalarType AvgL=std::min(L3,std::min(L1,L2));
//    if (((L2+L3)-AvgL)>L1)return 0;
//    //ScalarType CC=fabs(L1-(L2+L3))/AvgL;
//    ScalarType CC=((L1-(L2+L3))-AvgL)/AvgL;
//    CC=std::max(CC,(ScalarType)0);
//    assert(CC>=0);
//    //CC=std::min<ScalarType>(CC,1);
//    return CC;
//}

//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2,
//                     ScalarType L3,
//                     ScalarType L4,
//                     ScalarType L5)
//{
//    //ScalarType AvgL=(L1+L2+L3+L4+L5)/5;
//    ScalarType AvgL=std::min(L5,std::min(L1,std::min(L2,std::min(L3,L4))));
//    if ((L1+L2+L3-AvgL)>(L4+L5))return 0;
//    //ScalarType CC=fabs(L1+L2+L3-(L4+L5))/AvgL;
//    ScalarType CC=((L4+L5)-(L1+L2+L3)-AvgL)/AvgL;
//    CC=std::max(CC,(ScalarType)0);
//    assert(CC>=0);
//    //CC=std::min<ScalarType>(CC,1);
//    return CC;
//}


//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3)
//{
//    ScalarType CC0=OrderedCC(L1,L2,L3);
//    ScalarType CC1=OrderedCC(L2,L3,L1);
//    ScalarType CC2=OrderedCC(L3,L1,L2);
//    return (std::max(CC0,std::max(CC1,CC2)));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3,
//                               ScalarType L4)
//{
//    ScalarType CC0=OrderedCC(L1,L3);
//    ScalarType CC1=OrderedCC(L2,L4);
//    //std::cout<<"Value CC0: "<<CC0<<std::endl;
//    //std::cout<<"Value CC1: "<<CC1<<std::endl;
//    return (std::max(CC0,CC1));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3,
//                               ScalarType L4,
//                               ScalarType L5)
//{
//    ScalarType CC0=OrderedCC(L1,L2,L3,L4,L5);
//    ScalarType CC1=OrderedCC(L2,L3,L4,L5,L1);
//    ScalarType CC2=OrderedCC(L3,L4,L5,L1,L2);
//    ScalarType CC3=OrderedCC(L4,L5,L1,L2,L3);
//    ScalarType CC4=OrderedCC(L5,L1,L2,L3,L4);
//    return (std::max(CC0,std::max(CC1,std::max(CC2,std::max(CC3,CC4)))));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3,
//                               ScalarType L4,
//                               ScalarType L5,
//                               ScalarType L6)
//{
//    ScalarType CC0=OrderedCC(L1,L3,L5);
//    ScalarType CC1=OrderedCC(L2,L4,L6);
//    return (std::max(CC0,CC1));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(const std::vector<ScalarType> &EdgeL)
//{
////    assert(EdgeL.size()>=3);
////    assert(EdgeL.size()<=6);

//    if (EdgeL.size()==3)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2]));
//    if (EdgeL.size()==4)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3]));
//    if (EdgeL.size()==5)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4]));
//    if (EdgeL.size()==6)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4],EdgeL[5]));
//    return (std::numeric_limits<ScalarType>::max());
//}


//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2)
//{
//    ScalarType CC=std::max(L1,L2)/std::min(L1,L2);
//    return (CC);
//}

//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2,
//                     ScalarType L3)
//{
//    ScalarType CC=L1/(L2+L3);
//    if (CC<1)return 1;
//    return (CC);
//}

//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2,
//                     ScalarType L3,
//                     ScalarType L4,
//                     ScalarType L5)
//{
//    ScalarType CC=(L1+L2)/(L3+L4+L5);
//    if (CC<1)return 1;
//    return CC;
//}


//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3)
//{
//    ScalarType CC0=OrderedCC(L1,L2,L3);
//    ScalarType CC1=OrderedCC(L2,L3,L1);
//    ScalarType CC2=OrderedCC(L3,L1,L2);
//    return (std::max(CC0,std::max(CC1,CC2)));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3,
//                               ScalarType L4)
//{
//    ScalarType CC0=OrderedCC(L1,L3);
//    ScalarType CC1=OrderedCC(L2,L4);
//    return (std::max(CC0,CC1));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3,
//                               ScalarType L4,
//                               ScalarType L5)
//{
//    ScalarType CC0=OrderedCC(L1,L2,L3,L4,L5);
//    ScalarType CC1=OrderedCC(L2,L3,L4,L5,L1);
//    ScalarType CC2=OrderedCC(L3,L4,L5,L1,L2);
//    ScalarType CC3=OrderedCC(L4,L5,L1,L2,L3);
//    ScalarType CC4=OrderedCC(L5,L1,L2,L3,L4);
//    return (std::max(CC0,std::max(CC1,std::max(CC2,std::max(CC3,CC4)))));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(ScalarType L1,
//                               ScalarType L2,
//                               ScalarType L3,
//                               ScalarType L4,
//                               ScalarType L5,
//                               ScalarType L6)
//{
//    ScalarType CC0=OrderedCC(L1,L3,L5);
//    ScalarType CC1=OrderedCC(L2,L4,L6);
//    return (std::max(CC0,CC1));
//}

//template <class ScalarType>
//ScalarType CatmullClarkability(const std::vector<ScalarType> &EdgeL)
//{
//    //    assert(EdgeL.size()>=3);
//    //    assert(EdgeL.size()<=6);

//    if (EdgeL.size()==3)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2]));
//    if (EdgeL.size()==4)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3]));
//    if (EdgeL.size()==5)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4]));
//    if (EdgeL.size()==6)
//        return (CatmullClarkability(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4],EdgeL[5]));
//    return (std::numeric_limits<ScalarType>::max());
//}

//template <class ScalarType>
//size_t AddedSingularities(const std::vector<ScalarType> &EdgeL,
//                          const ScalarType &RatioThr)
//{
//    assert(RatioThr>1);
//    ScalarType CC=CatmullClarkability(EdgeL);
//    if (CC>RatioThr)
//    {
//        return 1;
//        //if (EdgeL.size()==4)return 2;
//        //else return 3;
//    }else
//    {
//        return 0;
////        if (EdgeL.size()==4)return 0;
////        else return 1;
//    }
//}


//template <class ScalarType>
//ScalarType OrderedCC(ScalarType L1,
//                     ScalarType L2)
//{
//    ScalarType CC=std::max(L1,L2)-std::min(L1,L2);
//    return (CC);
//}

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2,
                     ScalarType L3)
{
    ScalarType CC=L1+L2-L3;
    CC+=std::min(L1,std::min(L2,L3));
    //if (CC<0)return 0;
    return (CC);
}

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2,
                     ScalarType L3,
                     ScalarType L4,
                     ScalarType L5)
{
    ScalarType CC=(L3+L4+L5)-(L1+L2);
    CC+=std::min(L1,std::min(L2,std::min(L3,std::min(L4,L5))));
    //if (CC<0)return 0;
    return CC;
}


template <class ScalarType>
ScalarType CatmullClarkability3(ScalarType L1,
                                ScalarType L2,
                                ScalarType L3)
{
    ScalarType CC0=OrderedCC(L1,L2,L3);
    ScalarType CC1=OrderedCC(L2,L3,L1);
    ScalarType CC2=OrderedCC(L3,L1,L2);
    return (std::min(CC0,std::min(CC1,CC2)));
}

template <class ScalarType>
ScalarType CatmullClarkability4(const ScalarType L1,
                                const ScalarType L2,
                                const ScalarType L3,
                                const ScalarType L4)
{
//    ScalarType CC0=std::min(L1,L3)/std::max(L1,L3);
//    if (CC0<Ratio4)CC0=-1;
//    assert(CC0<=1);
//    if(fabs(L1-L3)<ThR)CC0=1;

//    ScalarType CC1=std::min(L2,L4)/std::max(L2,L4);
//    if (CC1<Ratio4)CC1=-1;
//    assert(CC1<=1);
//    if(fabs(L2-L4)<ThR)CC1=1;

//    return (std::min(CC0,CC1));

    ScalarType CC0=std::min(L1,L3)-std::max(L1,L3);
    ScalarType CC1=std::min(L2,L4)-std::max(L2,L4);
//    CC0+=std::min(L2,L4);
//    CC1+=std::min(L1,L3);
    ScalarType CC=std::min(CC0,CC1);
    CC+=std::min(L1,std::min(L2,std::min(L3,L4)));
    return CC;

}

template <class ScalarType>
ScalarType CatmullClarkability5(ScalarType L1,
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
    return (std::min(CC0,std::min(CC1,std::min(CC2,std::min(CC3,CC4)))));
}

template <class ScalarType>
ScalarType CatmullClarkability6(ScalarType L1,
                                ScalarType L2,
                                ScalarType L3,
                                ScalarType L4,
                                ScalarType L5,
                                ScalarType L6)
{
    ScalarType CC0=OrderedCC(L1,L3,L5);
    ScalarType CC1=OrderedCC(L2,L4,L6);
    return (std::min(CC0,CC1));
}

template <class ScalarType>
bool IsCatmullClarkable(//const int NumF,
                        const std::vector<ScalarType> &EdgeL,
                        const ScalarType &SideThr,
                        bool SkipValence4,
                        bool print_debug=false)
{
    //    assert(EdgeL.size()>=3);
    //    assert(EdgeL.size()<=6);
    //assert(Ratio4<1);
    //if (NumF<=EdgeL.size())return true;

    if ((SkipValence4)&&(EdgeL.size()==4))return true;

    ScalarType CC=0;
    if (EdgeL.size()==3)
        CC=(CatmullClarkability3(EdgeL[0],EdgeL[1],EdgeL[2]));
    if (EdgeL.size()==4)
        CC=(CatmullClarkability4(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3]));
    if (EdgeL.size()==5)
        CC=(CatmullClarkability5(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4]));
    if (EdgeL.size()==6)
        CC=(CatmullClarkability6(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4],EdgeL[5]));


    CC+=SideThr;

    if (print_debug)
    {
        std::cout<<"CC Computation "<<CC<<std::endl;
        std::cout<<"Side "<<EdgeL.size()<<std::endl;
        for (size_t i=0;i<EdgeL.size();i++)
            std::cout<<EdgeL[i]<<",";
        std::cout<<std::endl;

        std::cout<<"SideThr: "<<SideThr<<std::endl;
        std::cout<<"CC: "<<CC<<std::endl;
    }
    //CC=std::max(CC,(ScalarType)0);
    //    if ((EdgeL.size()==4)&&(CC>Ratio4))
    //        return true;
    //    if ((EdgeL.size()==4)&&(CC<=Ratio4))
    //        return false;

    return (CC>0);
}

template <class ScalarType>
size_t AddedSingularities(//const int NumF,
                          const std::vector<ScalarType> &EdgeL,
                          const ScalarType &SideThr,
                          bool SkipValence4,
                          bool print_debug=false)
{
    assert(SideThr>0);
    //bool CC=IsCatmullClarkable(NumF,EdgeL,SideThr,SkipValence4,print_debug);
    bool CC=IsCatmullClarkable(EdgeL,SideThr,SkipValence4,print_debug);
    //std::cout<<"test CC "<<CC<<std::endl;
    if (CC)
        return 0;
    else
        return 1;
}

#endif
