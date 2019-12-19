#ifndef QUADFROMPATCHES_ORIENT_FACES_H
#define QUADFROMPATCHES_ORIENT_FACES_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/geodesic.h>

template <class Mesh>
class OrientFaces
{
    typedef typename Mesh::FaceType FaceType;
    typedef typename Mesh::VertexType VertexType;
    typedef typename Mesh::CoordType CoordType;
    typedef typename Mesh::ScalarType ScalarType;

public:

    static void InvertFace(FaceType &f0)
    {
        std::vector<VertexType*> faceVert;
        for (int i=0;i<f0.VN();i++)
            faceVert.push_back(f0.V(i));
        std::reverse(faceVert.begin(),faceVert.end());
        for (int i=0;i<f0.VN();i++)
            f0.V(i)=faceVert[i];
    }

private:

    static bool IsCoherent(const FaceType &f0,const FaceType &f1,const int IndexE)
    {
        assert(f0.cFFp(IndexE)==&f1);
        assert(&f0!=&f1);
        const VertexType *v0=f0.cV(IndexE);
        const VertexType *v1=f0.cV((IndexE+1)%f0.VN());
        int IndexEopp=f0.cFFi(IndexE);
        assert(f1.cFFp(IndexEopp)==&f0);
        const VertexType *v2=f1.cV(IndexEopp);
        const VertexType *v3=f1.cV((IndexEopp+1)%f1.VN());
        if (v0==v2){assert(v1==v3);return false;}
        assert(v0==v3);
        assert(v1==v2);
        return true;
    }

    static void InvertFaces(Mesh &PolyM,
                            const std::vector<int> &ToInvert)
    {
        for (size_t i=0;i<ToInvert.size();i++)
            InvertFace(PolyM.face[ToInvert[i]]);
    }

    static void PropagateFrom(Mesh &PolyM,int &fI0,
                       std::vector<int> &OrientSet0,
                       std::vector<int> &OrientSet1)
    {
        OrientSet0.clear();
        OrientSet1.clear();

        std::vector<int> CoherentSet(PolyM.face.size(),-1);
        CoherentSet[fI0]=0;

        assert(!PolyM.face[fI0].IsS());

        std::vector<int> exploreStack;
        exploreStack.push_back(fI0);
        do{
            //get from the stack and set as explored
            int currF=exploreStack.back();
            exploreStack.pop_back();
            if(PolyM.face[currF].IsS())continue;//already explored
            PolyM.face[currF].SetS();

            //put in the right coherent set
            int currSet=CoherentSet[currF];
            if (currSet==0)
                OrientSet0.push_back(currF);
            else
                OrientSet1.push_back(currF);
            for (int i=0;i<PolyM.face[currF].VN();i++)
            {
                FaceType *f0=&PolyM.face[currF];
                FaceType *f1=PolyM.face[currF].FFp(i);
                if (f1==f0)continue;//border

                int IndexF1=vcg::tri::Index(PolyM,f1);
                if(PolyM.face[IndexF1].IsS())continue;

                exploreStack.push_back(IndexF1);//add to the stack

                //either is coherent or the opposite
                if (IsCoherent(*f0,*f1,i))
                    CoherentSet[IndexF1]=currSet;
                else
                    CoherentSet[IndexF1]=(currSet+1)%2;
            }
        }while(!exploreStack.empty());
    }

public:

    static void AutoOrientFaces(Mesh &PolyM)
    {
        vcg::tri::UpdateTopology<Mesh>::FaceFace(PolyM);
        vcg::tri::UpdateFlags<Mesh>::FaceClearS(PolyM);
        for (int i=0;i<(int)PolyM.face.size();i++)
        {
            if (PolyM.face[i].IsS())continue;
            std::cout<<"Reoriented face."<<std::endl;
            std::vector<int> OrientSet0,OrientSet1;
            PropagateFrom(PolyM,i,OrientSet0,OrientSet1);
            if (OrientSet0.size()<OrientSet1.size())
                InvertFaces(PolyM,OrientSet0);
            else
                InvertFaces(PolyM,OrientSet1);
        }
    }

};

#endif // QUADFROMPATCHES_ORIENT_FACES_H
