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

#ifndef META_MESH
#define META_MESH

#include <vector>
#include "patch_manager.h"
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include "edge_direction_table.h"

template <class MeshType>
class MetaMesh
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename vcg::face::Pos<FaceType> PosType;

//    enum CEdgeMode{CMNone,CMRemovable,CMLenght,CMPathId};

    std::set<size_t> AvoidPartitions;

    struct MetaVert
    {
        size_t MeshV;
        bool FixedEnd;
        MetaVert(size_t _MeshV,bool _FixedEnd)
        {
            MeshV=_MeshV;
            FixedEnd=_FixedEnd;
        }
    };


    struct MetaFace
    {
        std::vector<size_t> V;
        std::vector<ScalarType> MetaL;
        std::vector<std::pair<int,int> > AdjF;
        std::vector<int> PathId;
        std::vector<bool> RealV;
        std::vector<std::pair<size_t,size_t> > VertDir;
        CoordType BaryF;
        int ExpectedVal;
        int NumSing;

        //std::vector<std::vector<vcg::Color4b> > EdgeCol;
    };

public:
    MeshType *mesh;
    int MaxPath;
    ScalarType MaxL;
    std::vector<MetaVert> MVert;
    std::vector<MetaFace> MFaces;

private:
    std::vector<int> MetaToMeshVert;
    std::vector<int> MeshToMetaVert;

    EdgeDirectionTable EDTableMetaMesh;
    //std::map<std::pair<size_t,size_t>, std::vector<PosType> > MetaEdgeMeshPos;

    bool IsEmpty(const size_t indexMetaFace)
    {
        return (MFaces[indexMetaFace].V.size()==0);
    }

    size_t NumVerts(const size_t indexMetaFace)const
    {
        return (MFaces[indexMetaFace].V.size());
    }

public:

    size_t NumSides(const size_t indexMetaFace)const
    {
        size_t numS=0;
        for (size_t i=0;i<MFaces[indexMetaFace].V.size();i++)
        {
            if (!MFaces[indexMetaFace].RealV[i])continue;
            numS++;
        }
        return numS;
    }

    CoordType MetaVertPos(size_t &indexMetaVert)
    {
        assert(indexMetaVert<MVert.size());
        size_t IndexMeshV=MVert[indexMetaVert].MeshV;
        CoordType Pos=(*mesh).vert[IndexMeshV].P();
        return Pos;
    }


    void GetSideExtremes(const size_t indexMetaFace,
                         const size_t indexSide,
                         size_t &IndexV0,size_t &IndexV1)
    {
        std::vector<size_t> MetaEdgeI;
        GetSideMetaEdges(indexMetaFace,indexSide,MetaEdgeI);
        IndexV0=MetaEdgeI[0];
        size_t numV=NumVerts(indexMetaFace);
        IndexV1=(MetaEdgeI.back()+1)%numV;
    }

private:
    VertexType* MeshV(const size_t indexMetaFace,
                      const size_t &indexMetaVert)const
    {
        assert(indexMetaVert<NumVerts(indexMetaFace));
        int IndexMetaVGlobal=MFaces[indexMetaFace].V[indexMetaVert];
        int IndexMeshV=MVert[IndexMetaVGlobal].MeshV;
        assert(IndexMeshV=MetaToMeshVert[IndexMetaVGlobal]);
        return (&(*mesh).vert[IndexMeshV]);
    }

//    bool IsGenusOK(const size_t indexMetaFace)
//    {
//        vcg::tri::UnMarkAll(*mesh);
//        for (size_t i=0;i<MFaces[indexMetaFace].V.size();i++)
//        {
//            int IndexV=MFaces[indexMetaFace].V[i];
//            assert(IndexV>=0);
//            assert(IndexV<(*mesh).vert.size());
//            if (vcg::tri::IsMarked((*mesh),&(*mesh).vert[IndexV]))return false;
//        }
//        return true;
//    }
    bool IsGenusOK(const size_t indexMetaFace)
    {
        std::set<CoordType> VertPos;
        for (size_t i=0;i<MFaces[indexMetaFace].V.size();i++)
        {
            CoordType Pos=MetaVertPos(MFaces[indexMetaFace].V[i]);
            VertPos.insert(Pos);
            //int IndexV=MFaces[indexMetaFace].V[i];
//            assert(IndexV>=0);
//            assert(IndexV<(*mesh).vert.size());
//            if (vcg::tri::IsMarked((*mesh),&(*mesh).vert[IndexV]))return false;
        }
        return (VertPos.size()==MFaces[indexMetaFace].V.size());
        //return true;
    }

    size_t SideStartingE(const size_t indexMetaFace,
                         const size_t indexSide)const
    {
        size_t numS=0;
        for (size_t i=0;i<MFaces[indexMetaFace].V.size();i++)
        {
            if (!MFaces[indexMetaFace].RealV[i])
                continue;

            if (numS==indexSide)
                return i;

            numS++;
        }
        assert(0);
        return 0;
    }

    void SideLenghts(const size_t indexMetaFace,
                     std::vector<ScalarType> &SideLenghts)const
    {
        assert(MFaces[indexMetaFace].RealV[0]);
        ScalarType CurrL=0;
        size_t NumL=MFaces[indexMetaFace].MetaL.size();
        for (size_t i=0;i<NumL;i++)
        {
            CurrL+=MFaces[indexMetaFace].MetaL[i];
            size_t next_i=(i+1)%NumL;
            if (MFaces[indexMetaFace].RealV[next_i])
            {
                SideLenghts.push_back(CurrL);
                CurrL=0;
            }
        }
    }


    void CheckMetaFaceAdj(const size_t indexMetaFace)
    {
        assert(indexMetaFace<MFaces.size());
        assert(MFaces[indexMetaFace].AdjF.size()==NumVerts(indexMetaFace));
        for (size_t i=0;i<MFaces[indexMetaFace].AdjF.size();i++)
        {
            int IndexF=MFaces[indexMetaFace].AdjF[i].first;
            int IndexE=MFaces[indexMetaFace].AdjF[i].second;
            if (IndexF<0)continue;
            assert(IndexE>=0);
            assert(!IsEmpty(IndexF));
            assert(IndexE<(int)NumVerts(IndexF));
            int IndexF1=MFaces[IndexF].AdjF[IndexE].first;
            int IndexE1=MFaces[IndexF].AdjF[IndexE].second;
            if ((IndexF1!=(int)indexMetaFace)||(IndexE1!=(int)i))
            {
                std::cout<<" WARNING NON COHERENT ADJACENCY ACROSS "<<std::endl;
                std::cout<<"* Face "<<indexMetaFace<<" Edge "<<i<<std::endl;
                MFaces[indexMetaFace].AdjF[i].first=-1;
                MFaces[indexMetaFace].AdjF[i].second=-1;
//                std::cout<<"* Face "<<IndexF<<" Edge "<<IndexE<<std::endl;
//                std::cout<<"* Face "<<IndexF1<<" Edge "<<IndexE1<<std::endl;
            }
//            assert(IndexF1==(int)indexMetaFace);
//            assert(IndexE1==(int)i);
        }
    }

    size_t WhichSide(const size_t indexMetaFace,
                     const size_t indexMetaE)const
    {
        if (indexMetaE>=NumVerts(indexMetaFace))
        {
            std::cout<<"Index Meta E:"<<indexMetaE<<std::endl;
            std::cout<<"Num Meta E:"<<NumVerts(indexMetaFace)<<std::endl;
            assert(0);
        }
        size_t numS=0;
        for (size_t i=0;i<=indexMetaE;i++)
        {
            if (MFaces[indexMetaFace].RealV[i])
                numS++;
        }
        assert(numS>0);
        return (numS-1);
    }

    //    void WhichSideFaceHasPath(const size_t IndexPath,
    //                              std::vector<std::pair<size_t,size_t> > &FaceSides)
    //    {
    //        for (size_t i=0;i<MFaces.size();i++)
    //            for (size_t j=0;j<MFaces[i].PathId.size();j++)
    //            {
    //                if (MFaces[i].PathId[j]!=IndexPath)continue;
    //                size_t SideI=WhichSide(i,j);
    //                FaceSides.push_back(std::pair<size_t,size_t>(i,SideI));
    //            }
    //        std::sort(FaceSides.begin(),FaceSides.end());
    //        auto last = std::unique(FaceSides.begin(),FaceSides.end());
    //        FaceSides.erase(last,FaceSides.end());
    //    }

    int WhichEdgeFaceHasPath(const size_t &IndexMetaF,const size_t &IndexPath)const
    {
        for (size_t j=0;j<MFaces[IndexMetaF].PathId.size();j++)
        {
            if (MFaces[IndexMetaF].PathId[j]!=(int)IndexPath)continue;
            return j;
        }
        return -1;
    }


    int WhichSideFaceHasPath(const size_t &IndexMetaF,const size_t &IndexPath)const
    {
        int IndexE=WhichEdgeFaceHasPath(IndexMetaF,IndexPath);
        if (IndexE==-1)return -1;
        assert(IndexE<(int)NumVerts(IndexMetaF));

        return (WhichSide(IndexMetaF,IndexE));
    }

    void WhichSideFaceHasPath(const size_t IndexPath,std::vector<std::pair<size_t,size_t> > &FaceSides)
    {
        for (size_t i=0;i<MFaces.size();i++)
        {
            int IndexS=WhichSideFaceHasPath(i,IndexPath);
            if (IndexS==-1)continue;
            FaceSides.push_back(std::pair<size_t,size_t>(i,IndexS));
        }
        std::sort(FaceSides.begin(),FaceSides.end());
        auto last = std::unique(FaceSides.begin(),FaceSides.end());
        FaceSides.erase(last,FaceSides.end());
    }


    void PathSide(size_t &IndexMetaF,
                  const size_t IndexPath,
                  std::vector<std::pair<size_t,size_t> > &FaceSides)
    {
        for (size_t j=0;j<MFaces[IndexMetaF].PathId.size();j++)
        {
            if (MFaces[IndexMetaF].PathId[j]!=IndexPath)continue;
            assert(j<NumVerts(IndexMetaF));
            size_t SideI=WhichSide(IndexMetaF,j);
            FaceSides.push_back(std::pair<size_t,size_t>(IndexMetaF,SideI));
        }
        std::sort(FaceSides.begin(),FaceSides.end());
        auto last = std::unique(FaceSides.begin(),FaceSides.end());
        FaceSides.erase(last,FaceSides.end());
    }

    void GetSideMetaEdges(const size_t indexMetaFace,
                          const size_t indexSide,
                          std::vector<size_t> &MetaEdgeI)const
    {
        MetaEdgeI.clear();
        size_t nSides=NumSides(indexMetaFace);
        if(!(indexSide<nSides))
        {
            std::cout<<"I side:"<<indexSide<<std::endl;
            std::cout<<"Num side:"<<nSides<<std::endl;
            assert(0);
        }

        size_t VI0=SideStartingE(indexMetaFace,indexSide);
        size_t VI1=SideStartingE(indexMetaFace,(indexSide+1)%nSides);
        assert(VI0!=VI1);
        size_t numV=NumVerts(indexMetaFace);
        for (size_t i=VI0;i!=VI1;i=(i+1)%numV)
            MetaEdgeI.push_back(i);
    }


    void GetSideAdjacencyEdges(const size_t indexMetaFace,const size_t indexSide,
                               std::vector<std::pair<int,int> > &AdjF)const
    {
        std::vector<size_t> MetaEdgeI;
        GetSideMetaEdges(indexMetaFace,indexSide,MetaEdgeI);
        for (size_t i=0;i<MetaEdgeI.size();i++)
        {
            size_t IndexE=MetaEdgeI[i];
            std::pair<int,int> OppVal=MFaces[indexMetaFace].AdjF[IndexE];
            AdjF.push_back(OppVal);
        }
    }

    void GetSideAdjacencySides(const size_t indexMetaFace,
                               const size_t indexSide,
                               std::vector<std::pair<int,int> >  &FaceSides)const
    {
        std::vector<std::pair<int,int> > AdjF;
        GetSideAdjacencyEdges(indexMetaFace,indexSide,AdjF);
        FaceSides.clear();
        for (size_t i=0;i<AdjF.size();i++)
        {
            int IndexF=AdjF[i].first;
            int IndexE=AdjF[i].second;
            if (IndexF==-1)
            {
                FaceSides.push_back(std::pair<int,int>(IndexF,IndexE));
            }
            else
            {
                assert(IndexE<(int)NumVerts(IndexF));
                size_t IndexSideE=WhichSide(IndexF,IndexE);
                FaceSides.push_back(std::pair<int,int>(IndexF,IndexSideE));
            }
        }
        std::sort(FaceSides.begin(),FaceSides.end());
        auto last = std::unique(FaceSides.begin(),FaceSides.end());
        FaceSides.erase(last,FaceSides.end());
    }

    bool IsBorder(const size_t indexMetaFace,
                  const size_t indexE)const
    {
        return (MFaces[indexMetaFace].AdjF[indexE]==std::pair<int,int>(-1,-1));
    }

    void AddMetaVert(const std::vector<std::vector<size_t> > &FaceCorners,
                     const std::vector<size_t> &FixedV)
    {
        MVert.clear();
        MetaToMeshVert.clear();

        //set map to blocked one
        std::vector<bool> ConcaveNarrowVB=std::vector<bool> ((*mesh).vert.size(),false);
        for (size_t i=0;i<FixedV.size();i++)
        {
            assert(FixedV[i]<ConcaveNarrowVB.size());
            ConcaveNarrowVB[FixedV[i]]=true;
        }

        //allocate vertices
        MeshToMetaVert=std::vector<int>((*mesh).vert.size(),-1);
        for (size_t i=0;i<FaceCorners.size();i++)
            for (size_t j=0;j<FaceCorners[i].size();j++)
            {
                size_t IndexMeshV=FaceCorners[i][j];
                assert(IndexMeshV<(*mesh).vert.size());
                if (MeshToMetaVert[IndexMeshV]!=-1)continue;//already added vertx

                //otherwise add it
                MVert.push_back(MetaVert(IndexMeshV,ConcaveNarrowVB[IndexMeshV]));
                MeshToMetaVert[IndexMeshV]=MVert.size()-1;
                MetaToMeshVert.push_back(IndexMeshV);
            }

        //add the fixed ones in case not already there sometimes narrow might be non corners
        for (size_t i=0;i<FixedV.size();i++)
        {
            size_t IndexMeshV=FixedV[i];
            if (MeshToMetaVert[IndexMeshV]!=-1)continue;//already added vertx

            //otherwise add it
            MVert.push_back(MetaVert(IndexMeshV,ConcaveNarrowVB[IndexMeshV]));
            MeshToMetaVert[IndexMeshV]=MVert.size()-1;
            MetaToMeshVert.push_back(IndexMeshV);
        }

    }


    ScalarType MetaLenght(const std::vector<PosType> &EdgePos,
                          const std::map<std::pair<size_t,size_t>,ScalarType> &EdgeMap)
    {
        ScalarType MetaL=0;
        for (size_t k=0;k<EdgePos.size();k++)
        {
            size_t IndexV0=vcg::tri::Index((*mesh),EdgePos[k].VFlip());
            size_t IndexV1=vcg::tri::Index((*mesh),EdgePos[k].V());
            std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            assert(EdgeMap.count(key)>0);
            ScalarType currL=EdgeMap.at(key);
            MetaL+=currL;
        }
        return MetaL;
    }

    //save for each vert pait the face edge
    std::map<std::pair<size_t,size_t>,std::vector<std::pair<size_t,size_t> > > MeshEdgeMap;

    void UpdateAdjacency()
    {
        for (size_t i=0;i<MFaces.size();i++)
        {
            MFaces[i].AdjF.clear();
            MFaces[i].AdjF.resize(MFaces[i].V.size(),std::pair<int,int> (-1,-1));
        }

        std::map<std::pair<size_t,size_t>,std::vector<std::pair<size_t,size_t> > >::iterator IteEMap;
        //allocate
        for (IteEMap=MeshEdgeMap.begin();IteEMap!=MeshEdgeMap.end();IteEMap++)
        {
            if ((*IteEMap).second.size()<2)continue;//it is on the border
            assert((*IteEMap).second.size()==2);

            size_t IndexF0=(*IteEMap).second[0].first;
            size_t IndexE0=(*IteEMap).second[0].second;

            size_t IndexF1=(*IteEMap).second[1].first;
            size_t IndexE1=(*IteEMap).second[1].second;

            if (MFaces[IndexF0].AdjF[IndexE0]!=std::pair<int,int>(-1,-1))
            {
                std::pair<int,int> TestSide((int)IndexF1,(int)IndexE1) ;
                if(MFaces[IndexF0].AdjF[IndexE0]!=TestSide)
                {
//                    std::pair<int,int> TestAdj=MFaces[IndexF0].AdjF[IndexE0];
                    std::cout<<"- WARNING NON COHERENT ADJACENCY -"<<std::endl;
                    std::cout<<"* Face "<<IndexF0<<" Edge "<<IndexE0<<std::endl;
                    std::cout<<"* Face "<<IndexF1<<" Edge "<<IndexE1<<std::endl;
                    MFaces[IndexF0].AdjF[IndexE0].first=-1;
                    MFaces[IndexF0].AdjF[IndexE0].second=-1;
                    MFaces[IndexF1].AdjF[IndexE1].first=-1;
                    MFaces[IndexF1].AdjF[IndexE1].second=-1;
                    continue;
//                    std::cout<<"*Curr Adj is "<<TestSide.first<<","<<TestSide.second<<std::endl;
//                    std::cout<<"*Assugn Adj is "<<TestAdj.first<<","<<TestAdj.second<<std::endl;
                    //assert(0);
                }
            }
            if (MFaces[IndexF1].AdjF[IndexE1]!=std::pair<int,int>(-1,-1))
            {
                std::pair<int,int> TestSide((int)IndexF1,(int)IndexE1);
                //assert(MFaces[IndexF1].AdjF[IndexE1]==TestSide);
                if(MFaces[IndexF0].AdjF[IndexE0]!=TestSide)
                {
                    //std::pair<int,int> TestAdj=MFaces[IndexF0].AdjF[IndexE0];
                    std::cout<<"- WARNING NON COHERENT ADJACENCY -"<<std::endl;
                    std::cout<<"* Face "<<IndexF0<<" Edge "<<IndexE0<<std::endl;
                    std::cout<<"* Face "<<IndexF1<<" Edge "<<IndexE1<<std::endl;
                    MFaces[IndexF0].AdjF[IndexE0].first=-1;
                    MFaces[IndexF0].AdjF[IndexE0].second=-1;
                    MFaces[IndexF1].AdjF[IndexE1].first=-1;
                    MFaces[IndexF1].AdjF[IndexE1].second=-1;
                    continue;
//                    std::cout<<"*Curr Adj is "<<TestSide.first<<","<<TestSide.second<<std::endl;
//                    std::cout<<"*Assign Adj is "<<TestAdj.first<<","<<TestAdj.second<<std::endl;
//                    assert(0);
                }
            }

            //            if (MFaces[IndexF0].AdjF[IndexE0].first!=-1)
            //            {
            //                assert(MFaces[IndexF0].AdjF[IndexE0].first==(int)IndexF1);
            //                assert(MFaces[IndexF0].AdjF[IndexE0].second==(int)IndexE1);
            //            }
            MFaces[IndexF0].AdjF[IndexE0].first=IndexF1;
            MFaces[IndexF0].AdjF[IndexE0].second=IndexE1;

            //            if (MFaces[IndexF1].AdjF[IndexE1].first!=-1)
            //            {
            //                assert(MFaces[IndexF1].AdjF[IndexE1].first==(int)IndexF0);
            //                assert(MFaces[IndexF1].AdjF[IndexE1].second==(int)IndexE0);
            //            }

            MFaces[IndexF1].AdjF[IndexE1].first=IndexF0;
            MFaces[IndexF1].AdjF[IndexE1].second=IndexE0;
        }

        std::cout<<"Initializing Adjacency"<<std::endl;
        for (size_t i=0;i<MFaces.size();i++)
        {
            //std::cout<<"check "<<i<<" out of "<<MFaces.size()<<std::endl;
            CheckMetaFaceAdj(i);
        }
        std::cout<<"Done"<<std::endl;
    }

    size_t UnsolvedEmitters(const size_t indexMetaFace)
    {
        (void)indexMetaFace;
        return 0;//HAS TO BE IMPLEMENTED
    }

    void GetPatchInfo(const size_t indexMetaFace,
                      PatchInfo<ScalarType> &PInfo,
                      ScalarType Thr,
                      bool SkipValence4)
    {
        PInfo.NumEmitters=UnsolvedEmitters(indexMetaFace);
        PInfo.NumCorners=NumSides(indexMetaFace);
        PInfo.Genus=IsGenusOK(indexMetaFace)?1:-1;
        //        std::vector<ScalarType> SideL;
        //        SideLenghts(indexMetaFace,SideL);

        SideLenghts(indexMetaFace,PInfo.CurvedL);
        PInfo.ExpectedValence=MFaces[indexMetaFace].ExpectedVal;
        PInfo.NumSing=MFaces[indexMetaFace].NumSing;

        if ((PInfo.NumCorners<(int)MIN_ADMITTIBLE)||
                (PInfo.NumCorners>(int)MAX_ADMITTIBLE)||
                (PInfo.Genus!=1)||
                (PInfo.NumEmitters>0)||
                (Thr<=0))
            PInfo.CClarkability=false;
        else
            PInfo.CClarkability=IsCatmullClarkable(PInfo.CurvedL,Thr,SkipValence4);
        //PInfo.CClarkability=IsCatmullClarkable(indexMetaFace,PInfo.CurvedL,Thr,SkipValence4);
    }

    void GetFaceEdgeFromVertPair(const size_t IndexV0,const size_t IndexV1,
                                 std::vector<std::pair<size_t,size_t> > &FaceEdge)
    {
        FaceEdge.clear();
        std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
        for (size_t i=0;i<MeshEdgeMap[Key].size();i++)
            FaceEdge.push_back(MeshEdgeMap[Key][i]);
    }

    void UpdatePathId(const std::vector<std::vector<size_t > > &PathSeq,
                      const std::vector<bool> &IsLoop)
    {
        for (size_t i=0;i<MFaces.size();i++)
        {
            MFaces[i].PathId.clear();
            MFaces[i].PathId.resize(MFaces[i].V.size(),-1);
        }

        for (size_t i=0;i<PathSeq.size();i++)
        {
            size_t IndexPath=i;
            size_t Limit=PathSeq[i].size();
            if (!IsLoop[i])Limit--;
            size_t SizePath=PathSeq[i].size();
            for (size_t j=0;j<Limit;j++)
            {
                size_t Vindex0=PathSeq[i][j];
                size_t Vindex1=PathSeq[i][(j+1)%SizePath];
                std::vector<std::pair<size_t,size_t> > FaceEdge;
                GetFaceEdgeFromVertPair(Vindex0,Vindex1,FaceEdge);
                for (size_t k=0;k<FaceEdge.size();k++)
                {
                    size_t IndexF=FaceEdge[k].first;
                    if (AvoidPartitions.count(IndexF)>0)continue;
                    size_t IndexE=FaceEdge[k].second;
                    assert(MFaces[IndexF].V.size()>0);
                    assert(IndexE<MFaces[IndexF].PathId.size()>0);
                    if (MFaces[IndexF].PathId[IndexE]==-1)
                        MFaces[IndexF].PathId[IndexE]=IndexPath;
                    else
                    {
                        if (MFaces[IndexF].PathId[IndexE]!=(int)IndexPath)
                        {
                        MFaces[IndexF].PathId[IndexE]=(int)IndexPath;
                        std::cout<<"WARNING INCONSISTENT META MESH"<<std::endl;
                        }
                        //std::cout<<"WARNING INCONSISTENT META MESH"<<std::endl;
                        //assert(MFaces[IndexF].PathId[IndexE]==(int)IndexPath);
                    }
                }
            }
        }
    }

    //for each
    //std::vector<std::pair<size_t,size_t>,std::pair<size_t,size_t> > EdgeDir;

    void UpdateEdgeMeshMap(const size_t &MFaceIndex,
                           const size_t &MEdgeIndex,
                           const std::vector<PosType> &SideEdge)
    {
        for (size_t i=0;i<SideEdge.size();i++)
        {
            size_t IndexV0=vcg::tri::Index((*mesh),SideEdge[i].VFlip());
            size_t IndexV1=vcg::tri::Index((*mesh),SideEdge[i].V());
            std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            MeshEdgeMap[key].push_back(std::pair<size_t,size_t> (MFaceIndex,MEdgeIndex));
        }
    }

    void GetEdgeDirs(const std::vector<PosType> &EdgePos,
                     const EdgeDirectionTable &EDTableMesh,
                     std::pair<size_t,size_t> &EdgeDirs)
    {
        //get first vertex
        size_t V0=vcg::tri::Index((*mesh),EdgePos[0].VFlip());
        //and the one next to him in the sequence
        size_t V0e=vcg::tri::Index((*mesh),EdgePos[0].V());
        assert(V0!=V0e);

        //get last vertex
        size_t V1=vcg::tri::Index((*mesh),EdgePos.back().V());
        //and the one next to him in the sequence
        size_t V1e=vcg::tri::Index((*mesh),EdgePos.back().VFlip());
        assert(V1!=V1e);

        assert(V0!=V1);

        //retrieve the directions in the original mesh
        EdgeVert EV0(V0,V0e,V0);
        EdgeVert EV1(V1,V1e,V1);
        size_t Dir0=EDTableMesh.GetEdgeDir(EV0);
        size_t Dir1=EDTableMesh.GetEdgeDir(EV1);
        EdgeDirs.first=Dir0;
        EdgeDirs.second=Dir1;

//        //get the vertex in the metamesh
//        int V0Meta=MeshToMetaVert[V0];
//        int V1Meta=MeshToMetaVert[V1];
//        assert(V0Meta>=0);
//        assert(V1Meta>=0);

//        EdgeVert EV0Meta(V0Meta,V1Meta,V0Meta);
//        EdgeVert EV1Meta(V0Meta,V1Meta,V1Meta);
//        if (EDTableMetaMesh.HasEdgeDir(EV0Meta))
//        {
//            assert(EDTableMetaMesh.HasEdgeDir(EV1Meta));
//            size_t Dir0Test=EDTableMetaMesh.GetEdgeDir(EV0Meta);
//            size_t Dir1Test=EDTableMetaMesh.GetEdgeDir(EV1Meta);
//            if (Dir0Test!=Dir0)
//            {
//                std::cout<<"Error Inconsistent Dir0 Dirs "<<Dir0Test<<"!="<<Dir0<<std::endl;
//                assert(0);
//            }
//            if (Dir1Test!=Dir1)
//            {
//                std::cout<<"Error Inconsistent Dir 1 Dirs "<<Dir1Test<<"!="<<Dir1<<std::endl;
//                assert(0);
//            }
//            //assert(Dir0Test==Dir0);
//            //assert(Dir1Test==Dir1);
//        }
//        else
//        {
//            EDTableMetaMesh.AddEdgeDir(EV0Meta,Dir0);
//            EDTableMetaMesh.AddEdgeDir(EV1Meta,Dir1);
//        }
        //then add to the map
        //EDTableMesh.GetEdgeDir(const EdgeVert &EV)()
    }

    void RetrieveMetaFace(const int &MFaceIndex,
                          const std::vector<size_t> &FaceCorners,
                          const std::vector<size_t> &PatchFaces,
                          const std::map<std::pair<size_t,size_t>,ScalarType> &EdgeMap,
                          const EdgeDirectionTable &EDTableMesh)
    {

        //find the first pos from a corner of that face
        assert(FaceCorners.size()>=2);
        assert(PatchFaces.size()>0);

        //        //select corners
        //        PatchManager<MeshType>::SelectVertices((*mesh),FaceCorners);

        //then retrieve the meta edges
        std::vector<std::vector<PosType> > BorderSeqPos;
        PatchManager<MeshType>::DeriveFauxSeqPos((*mesh),PatchFaces,BorderSeqPos,FaceCorners[0]);

        std::set<size_t> FaceCornerSet(FaceCorners.begin(),FaceCorners.end());

        //then add the vertices to the mesh
        MFaces[MFaceIndex].V.clear();
        MFaces[MFaceIndex].VertDir.clear();
        for (size_t i=0;i<BorderSeqPos.size();i++)
        {
            assert(BorderSeqPos[i].size()>0);
            //get the first corner
            size_t IndexMeshV=vcg::tri::Index((*mesh),BorderSeqPos[i][0].VFlip());
            assert(MeshToMetaVert[IndexMeshV]!=-1);//already added vertx
            MFaces[MFaceIndex].V.push_back(MeshToMetaVert[IndexMeshV]);

            //then compute the meta-lenght
            ScalarType EdgePosL=MetaLenght(BorderSeqPos[i],EdgeMap);
            MFaces[MFaceIndex].MetaL.push_back(EdgePosL);
            MaxL=std::max(MaxL,EdgePosL);

            UpdateEdgeMeshMap(MFaceIndex,i,BorderSeqPos[i]);

            //then save per edge per vertex directions
            //UpdateEdgeDir(BorderSeqPos[i],EDTableMesh);
            std::pair<size_t,size_t> EdgeDir;
            GetEdgeDirs(BorderSeqPos[i],EDTableMesh,EdgeDir);
            MFaces[MFaceIndex].VertDir.push_back(EdgeDir);
        }

        MFaces[MFaceIndex].RealV.resize(MFaces[MFaceIndex].V.size(),false);
        for (size_t j=0;j<MFaces[MFaceIndex].V.size();j++)
        {
            size_t IndexV=MFaces[MFaceIndex].V[j];
            size_t IndexVMesh=MetaToMeshVert[IndexV];
            if (FaceCornerSet.count(IndexVMesh)==0)
                continue;
            else
                MFaces[MFaceIndex].RealV[j]=true;
        }
        assert(MFaces[MFaceIndex].RealV[0]==true);
    }

    void AddMetaFaces(const std::vector<int> &FacePatches,
                      const std::vector<std::vector<size_t> > &FaceCorners,
                      const std::map<std::pair<size_t,size_t>,ScalarType> &EdgeMap,
                      const EdgeDirectionTable &EDTableMesh)
    {

        //select patch borders
        PatchManager<MeshType>::SelectMeshPatchBorders((*mesh),FacePatches,true);

        //select the vertices
        vcg::tri::UpdateSelection<MeshType>::VertexClear((*mesh));
        for (size_t i=0;i<MVert.size();i++)
            (*mesh).vert[MVert[i].MeshV].SetS();

        //PatchManager<MeshType>::SelectVertices((*mesh),FaceCorners);

        //derive pert partition faces
        std::vector<std::vector<size_t> > PartitionFaces;
        PatchManager<MeshType>::DerivePerPartitionFaces(FacePatches,PartitionFaces);

        //allocate faces
        MFaces.resize(PartitionFaces.size());

        //then start retrieving the infos
        MaxL=0;
        //size_t numE=0;
        //MetaEdgeMeshPos.clear();
        for (size_t i=0;i<MFaces.size();i++)
        {
            if (AvoidPartitions.count(i)>0)continue;
            RetrieveMetaFace(i,FaceCorners[i],PartitionFaces[i],EdgeMap,EDTableMesh);
            //numE+=MFaces[i].Edges.size();
        }

        for (size_t i=0;i<MFaces.size();i++)
        {
            MFaces[i].BaryF=CoordType(0,0,0);
            for (size_t j=0;j<MFaces[i].V.size();j++)
            {
                CoordType Pos=MetaVertPos(MFaces[i].V[j]);
                MFaces[i].BaryF+=Pos;
            }
            MFaces[i].BaryF/=MFaces[i].V.size();
        }
    }

    void GetAdjacencySidesBetween(const size_t IndexMetaF0,const size_t IndexMetaF1,
                                  std::vector<size_t> &SideIndexes)const
    {
        std::vector<size_t> EdgeIndex;
        assert(MFaces[IndexMetaF0].RealV[0]);
        assert(MFaces[IndexMetaF0].AdjF.size()==NumVerts(IndexMetaF0));
        for (size_t i=0;i<MFaces[IndexMetaF0].AdjF.size();i++)
        {
            if (MFaces[IndexMetaF0].AdjF[i].first==-1)continue;

            //if there is adjacency and not already added the same side
            if (MFaces[IndexMetaF0].AdjF[i].first==(int)IndexMetaF1)
                EdgeIndex.push_back(i);
        }

        SideIndexes.clear();
        for (size_t i=0;i<EdgeIndex.size();i++)
        {
            assert(EdgeIndex[i]<NumVerts(IndexMetaF0));
            size_t IndexS=WhichSide(IndexMetaF0,EdgeIndex[i]);
            SideIndexes.push_back(IndexS);
        }
        std::sort(SideIndexes.begin(),SideIndexes.end());
        auto last = std::unique(SideIndexes.begin(),SideIndexes.end());
        SideIndexes.erase(last,SideIndexes.end());
    }

    bool IsRemovableSide(const size_t IndexMetaF,
                         const size_t IndexSide)const
    {
        //std::cout<<"0"<<std::endl;
        std::vector<size_t> MetaEdgeI;
        GetSideMetaEdges(IndexMetaF,IndexSide,MetaEdgeI);
        //std::cout<<"1"<<std::endl;
        size_t sizeV=MFaces[IndexMetaF].V.size();
        for (size_t i=0;i<MetaEdgeI.size();i++)
        {
            size_t IndexMetaE=MetaEdgeI[i];
            if (IsBorder(IndexMetaF,IndexMetaE))return false;
            size_t metaVSide0=MFaces[IndexMetaF].V[IndexMetaE];
            size_t metaVSide1=MFaces[IndexMetaF].V[(IndexMetaE+1)%sizeV];
            if (MVert[metaVSide0].FixedEnd)return false;
            if (MVert[metaVSide1].FixedEnd)return false;
        }
        //get adjacency
        std::vector<std::pair<int,int> >  OppFaceSides;
        GetSideAdjacencySides(IndexMetaF,IndexSide,OppFaceSides);
        //std::cout<<"3"<<std::endl;
        if (OppFaceSides.size()>1)return false;

        //see if another side linked to the same metaface
        assert(OppFaceSides.size()==1);
        int IndexMetaF0=IndexMetaF;
        int IndexMetaF1=OppFaceSides[0].first;
        if (IndexMetaF1==IndexMetaF0)return false;// adjacent to itself

        std::vector<size_t> SideIndexes;
        GetAdjacencySidesBetween(IndexMetaF0,IndexMetaF1,SideIndexes);
        if (SideIndexes.size()>1)return false;
        //       }

        //check if the two sides coincide
        std::vector<size_t> MetaEdgeOpp;
        GetSideMetaEdges(OppFaceSides[0].first,OppFaceSides[0].second,MetaEdgeOpp);
        //std::cout<<"4"<<std::endl;
        if (MetaEdgeOpp.size()!=MetaEdgeI.size())return false;
        return true;
    }

    bool IsRemovableSideVisual(const size_t IndexMetaF,
                               const size_t IndexSide)const
    {
        //std::cout<<"0"<<std::endl;
        std::vector<size_t> MetaEdgeI;
        GetSideMetaEdges(IndexMetaF,IndexSide,MetaEdgeI);
        //std::cout<<"1"<<std::endl;
        size_t sizeV=MFaces[IndexMetaF].V.size();
        for (size_t i=0;i<MetaEdgeI.size();i++)
        {
            size_t IndexMetaE=MetaEdgeI[i];
            if (IsBorder(IndexMetaF,IndexMetaE))return false;
            size_t metaVSide0=MFaces[IndexMetaF].V[IndexMetaE];
            size_t metaVSide1=MFaces[IndexMetaF].V[(IndexMetaE+1)%sizeV];
            if (MVert[metaVSide0].FixedEnd)return false;
            if (MVert[metaVSide1].FixedEnd)return false;
        }
        //get adjacency
        std::vector<std::pair<int,int> >  OppFaceSides;
        GetSideAdjacencySides(IndexMetaF,IndexSide,OppFaceSides);
        //std::cout<<"3"<<std::endl;
        if (OppFaceSides.size()>1)return false;

        //see if another side linked to the same metaface
        assert(OppFaceSides.size()==1);
        int IndexMetaF0=IndexMetaF;
        int IndexMetaF1=OppFaceSides[0].first;
        if (IndexMetaF1==IndexMetaF0)return false;// adjacent to itself

        std::vector<size_t> SideIndexes;
        GetAdjacencySidesBetween(IndexMetaF0,IndexMetaF1,SideIndexes);
        if (SideIndexes.size()>1)return false;


        //       }

        //check if the two sides coincide
        std::vector<size_t> MetaEdgeOpp;
        GetSideMetaEdges(OppFaceSides[0].first,OppFaceSides[0].second,MetaEdgeOpp);
        //std::cout<<"4"<<std::endl;
        if (MetaEdgeOpp.size()!=MetaEdgeI.size())return false;
        return true;
    }

public:

    bool IsFieldCornerFaceV(const size_t IndexMetaF,
                            const size_t IndexMetaV)
    {
        assert(IndexMetaV<NumVerts(IndexMetaF));
        std::vector<size_t> EdgeDir;
        GetEdgeDirOnV(IndexMetaF,IndexMetaV,EdgeDir);
        return (HasOrthogonalCross(EdgeDir)||
                HasNarrowCross(EdgeDir));
    }

    bool IsRemovableEdgeVisual(const size_t IndexMetaF,
                               const size_t IndexEdge)const
    {
        assert(IndexEdge<NumVerts(IndexMetaF));
        size_t IndexSide=WhichSide(IndexMetaF,IndexEdge);
        return (IsRemovableSideVisual(IndexMetaF,IndexSide));
        return true;
    }


    void GetEdgeDirOnV(const size_t IndexMetaF,
                       const size_t IndexMetaV,
                       std::vector<CoordType> &EdgeDir)
    {
        std::vector<size_t> IDir;
        GetEdgeDirOnV(IndexMetaF,IndexMetaV,IDir);

        VertexType *v= MeshV(IndexMetaF,IndexMetaV);
        for (size_t i=0;i<IDir.size();i++)
        {
            assert(IDir[i]>=0);
            assert(IDir[i]<4);
            CoordType CrossDir=vcg::tri::CrossField<MeshType>::CrossVector(*v,(int)IDir[i]);
            EdgeDir.push_back(CrossDir);
        }
    }

private:

    bool IsRemovableEdge(const size_t IndexMetaF,
                         const size_t IndexEdge)const
    {
        assert(IndexEdge<NumVerts(IndexMetaF));
        size_t IndexSide=WhichSide(IndexMetaF,IndexEdge);
        return (IsRemovableSide(IndexMetaF,IndexSide));
        return true;
    }

    void GetEdgeDirOnV(const size_t IndexMetaF,
                       const size_t IndexMetaV,
                       std::vector<size_t> &EdgeDir)
    {
        assert(IndexMetaV<NumVerts(IndexMetaF));
        assert(IndexMetaV<MFaces[IndexMetaF].VertDir.size());
        assert(MFaces[IndexMetaF].VertDir.size()==NumVerts(IndexMetaF));
        size_t sizeV=NumVerts(IndexMetaF);
        //get current and prev
        size_t IndexE_curr=IndexMetaV;
        size_t IndexE_prev=(IndexMetaV+sizeV-1)%sizeV;
        size_t Dir0=MFaces[IndexMetaF].VertDir[IndexE_curr].first;
        size_t Dir1=MFaces[IndexMetaF].VertDir[IndexE_prev].second;
        EdgeDir.push_back(Dir0);
        EdgeDir.push_back(Dir1);
    }



//    void GlDrawMetaVert(size_t IndexV)
//    {
//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glDisable(GL_LIGHTING);
//        glDepthRange(0,0.9995);

//        if (MVert[IndexV].FixedEnd)
//        {
//            glPointSize(20);
//            glColor(vcg::Color4b(255,0,0,255));
//        }
//        else
//        {
//            glPointSize(10);
//            glColor(vcg::Color4b(0,255,0,255));
//        }

//        CoordType Pos = MetaVertPos(IndexV);
//        glBegin(GL_POINTS);
//        vcg::glVertex(Pos);
//        glEnd();
//        glPopAttrib();
//    }

//    void GlDrawMetaFaceAdj(size_t IndexF)
//    {
//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glDisable(GL_LIGHTING);
//        glDepthRange(0,0.9995);

//        glLineWidth(5);


//        size_t sizeV=MFaces[IndexF].V.size();

//        CoordType Bary0=MFaces[IndexF].BaryF;
//        for (size_t i=0;i<sizeV;i++)
//        {
//            CoordType PosE0= MetaVertPos(MFaces[IndexF].V[i]);
//            CoordType PosE1= MetaVertPos(MFaces[IndexF].V[(i+1)%sizeV]);
//            CoordType AvgE=(PosE0+PosE1)/2;
//            int AdjF=MFaces[IndexF].AdjF[i].first;

//            CoordType Pos0=Bary0*0.5+AvgE*0.5;
//            CoordType Pos1=AvgE;
//            CoordType Pos2=AvgE;
//            vcg::Color4b adj_col(255,0,255,255);
//            if (AdjF!=-1)
//            {
//                assert(MFaces[IndexF].AdjF[i].second!=-1);
//                int AdjF=MFaces[IndexF].AdjF[i].first;
//                //int AdjE=MFaces[IndexF].AdjF[i].second;
//                CoordType Bary1=MFaces[AdjF].BaryF;
//                Pos2=Bary1*0.5+AvgE*0.5;//Bary1;
//                adj_col=vcg::Color4b(0,0,255,255);
//            }
//            glColor(adj_col);
//            glBegin(GL_LINE_STRIP);
//            vcg::glVertex(Pos0);
//            vcg::glVertex(Pos1);
//            vcg::glVertex(Pos2);
//            glEnd();
//        }
//        glPopAttrib();
//    }

//    void GlDrawMetaFaceVert(size_t IndexF)
//    {
//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glDisable(GL_LIGHTING);
//        glDepthRange(0,0.9995);

//        glPointSize(10);
//        glColor(vcg::Color4b(0,0,0,255));

//        size_t sizeV=MFaces[IndexF].V.size();

//        CoordType bary=MFaces[IndexF].BaryF;

//        glBegin(GL_POINTS);
//        for (size_t i=0;i<sizeV;i++)
//        {
//            if (!MFaces[IndexF].RealV[i])continue;
//            CoordType pos=MetaVertPos(MFaces[IndexF].V[i]);
//            pos=pos*0.8+bary*0.2;
//            vcg::glVertex(pos);
//        }
//        glEnd();


//        glPopAttrib();
//    }

//    void GlDrawMetaFaceVertField(size_t IndexF)
//    {
//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glDisable(GL_LIGHTING);
//        glDepthRange(0,0.9995);


//        size_t sizeV=MFaces[IndexF].V.size();

//        CoordType bary=MFaces[IndexF].BaryF;

//        for (size_t i=0;i<sizeV;i++)
//        {
//            CoordType pos0=MetaVertPos(MFaces[IndexF].V[i]);
//            pos0=pos0*0.8+bary*0.2;

//            std::vector<CoordType> EdgeDir;
//            GetEdgeDirOnV(IndexF,i,EdgeDir);
//            assert(EdgeDir.size()==2);
//            ScalarType sizeD=(*mesh).bbox.Diag()*0.01;
//            CoordType pos1=pos0+EdgeDir[0]*sizeD;
//            CoordType pos2=pos0+EdgeDir[1]*sizeD;
//            if (IsFieldCornerFaceV(IndexF,i))
//            {
//               glLineWidth(10);
//               glColor(vcg::Color4b(255,255,0,255));
//            }
//            else
//            {
//               glLineWidth(5);
//               glColor(vcg::Color4b(200,200,200,255));
//            }
//            glBegin(GL_LINES);
//            vcg::glVertex(pos0);
//            vcg::glVertex(pos1);
//            vcg::glVertex(pos0);
//            vcg::glVertex(pos2);
//            glEnd();
//        }



//        glPopAttrib();
//    }

//    void GlDrawMetaFaceSides(size_t IndexF)
//    {
//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glDisable(GL_LIGHTING);
//        glDepthRange(0,0.9995);
//        glDisable(GL_CULL_FACE);
//        //glLineWidth(5);


//        size_t sizeV=NumSides(IndexF);
//        int ExpVal=MFaces[IndexF].ExpectedVal;
//        if ((sizeV<3)||(sizeV>6))
//            glColor(vcg::Color4b(255,0,0,255));
//        if (sizeV==3)
//            glColor(vcg::Color4b(0,255,255,255));
//        if (sizeV==4)
//            glColor(vcg::Color4b(100,100,100,255));
//        if (sizeV==5)
//            glColor(vcg::Color4b(255,0,255,255));
//        if (sizeV==6)
//            glColor(vcg::Color4b(255,255,0,255));
//        if (ExpVal!=sizeV)
//            glColor(vcg::Color4b(255,0,0,255));

//        CoordType bary=MFaces[IndexF].BaryF;

//        //        for (size_t i=0;i<sizeV;i++)
//        //        {
//        //            size_t IndexV0,IndexV1;
//        //            GetSideExtremes(IndexF,i,IndexV0,IndexV1);

//        //            CoordType Pos0= MetaVertPos(MFaces[IndexF].V[IndexV0]);
//        //            Pos0=bary*0.75+Pos0*0.25;
//        //            CoordType Pos1= MetaVertPos(MFaces[IndexF].V[IndexV1]);
//        //            Pos1=bary*0.75+Pos1*0.25;
//        //            glBegin(GL_LINES);
//        //            vcg::glVertex(Pos0);
//        //            vcg::glVertex(Pos1);
//        //            glEnd();
//        //        }

//        glBegin(GL_POLYGON);
//        for (size_t i=0;i<sizeV;i++)
//        {
//            size_t IndexV0,IndexV1;
//            GetSideExtremes(IndexF,i,IndexV0,IndexV1);

//            CoordType Pos0= MetaVertPos(MFaces[IndexF].V[IndexV0]);
//            Pos0=bary*0.75+Pos0*0.25;

//            vcg::glVertex(Pos0);
//        }
//        glEnd();

//        glPopAttrib();
//    }



//    void GlDrawMetaFace(size_t IndexF,const CEdgeMode EMode)
//    {
//        glPushAttrib(GL_ALL_ATTRIB_BITS);
//        glDisable(GL_LIGHTING);
//        glDepthRange(0,0.9995);

//        glLineWidth(5);


//        size_t sizeV=MFaces[IndexF].V.size();
//        for (size_t i=0;i<sizeV;i++)
//        {
//            if (EMode==CMNone)
//                glColor(vcg::Color4b(200,200,200,255));
//            if (EMode==CMLenght)
//                glColor(vcg::Color4b::ColorRamp(0,MaxL,MFaces[IndexF].MetaL[i]));
//            if (EMode==CMRemovable)
//            {
//                if (IsRemovableEdgeVisual(IndexF,i))
//                    glColor(vcg::Color4b(0,255,0,255));
//                else
//                    glColor(vcg::Color4b(255,0,0,255));
//            }
//            if (EMode==CMPathId)
//            {
//                size_t IndexP=MFaces[IndexF].PathId[i];
//                glColor(vcg::Color4b::Scatter(MaxPath,IndexP));
//            }
//            CoordType Pos0= MetaVertPos(MFaces[IndexF].V[i]);
//            CoordType Pos1= MetaVertPos(MFaces[IndexF].V[(i+1)%sizeV]);
//            glBegin(GL_LINES);
//            vcg::glVertex(Pos0);
//            vcg::glVertex(Pos1);
//            glEnd();
//        }
//        glPopAttrib();
//    }

    void ReSortFace(size_t IndexF)
    {
        size_t StartE=0;
        size_t sizeE=MFaces[IndexF].RealV.size();
        for (size_t i=0;i<sizeE;i++)
        {
            if (!MFaces[IndexF].RealV[i])continue;
            StartE=i;
            break;
        }
        if (StartE==0)return;

        std::vector<size_t> NewV;
        std::vector<ScalarType> NewMetaL;
        std::vector<std::pair<int,int> > NewAdjF;
        std::vector<int> NewPathId;
        std::vector<bool> NewRealV;
        std::vector<std::pair<size_t,size_t> > NewVertDir;

        for (size_t i=0;i<sizeE;i++)
        {
            size_t Index=(i+StartE)%sizeE;
            NewV.push_back(MFaces[IndexF].V[Index]);
            NewMetaL.push_back(MFaces[IndexF].MetaL[Index]);
            NewAdjF.push_back(MFaces[IndexF].AdjF[Index]);
            NewPathId.push_back(MFaces[IndexF].PathId[Index]);
            NewRealV.push_back(MFaces[IndexF].RealV[Index]);
            NewVertDir.push_back(MFaces[IndexF].VertDir[Index]);
        }

        MFaces[IndexF].V=NewV;
        MFaces[IndexF].MetaL=NewMetaL;
        MFaces[IndexF].AdjF=NewAdjF;
        MFaces[IndexF].PathId=NewPathId;
        MFaces[IndexF].RealV=NewRealV;
        MFaces[IndexF].VertDir=NewVertDir;
        assert(MFaces[IndexF].RealV[0]);
    }

    void MakeAdjacencyCoherentFrom(size_t IndexF)
    {
        for (size_t i=0;i<MFaces[IndexF].AdjF.size();i++)
        {
            int IndexF1=MFaces[IndexF].AdjF[i].first;
            int IndexE1=MFaces[IndexF].AdjF[i].second;
            if (IndexF1==-1)continue;
            assert(IndexE1!=-1);
            MFaces[IndexF1].AdjF[IndexE1].first=IndexF;
            MFaces[IndexF1].AdjF[IndexE1].second=i;
        }

        for (size_t i=0;i<MFaces[IndexF].AdjF.size();i++)
        {
            int IndexF1=MFaces[IndexF].AdjF[i].first;
            if (IndexF1==-1)continue;
            CheckMetaFaceAdj(IndexF1);
        }
    }

    void RemoveEmptyFaces()
    {
        std::vector<int> NewMap(MFaces.size(),-1);
        std::vector<MetaFace> NewMFaces;
        for (size_t i=0;i<MFaces.size();i++)
        {
            if (IsEmpty(i))continue;
            NewMFaces.push_back(MFaces[i]);
            NewMap[i]=NewMFaces.size()-1;
        }
        for (size_t i=0;i<NewMFaces.size();i++)
            for (size_t j=0;j<NewMFaces[i].AdjF.size();j++)
            {
                //change to new face index
                int OldIdF=NewMFaces[i].AdjF[j].first;
                if (OldIdF==-1)continue;
                //update face adj if needed
                assert(NewMap[OldIdF]!=-1);
                NewMFaces[i].AdjF[j].first=NewMap[OldIdF];
            }
        MFaces=NewMFaces;
    }

    void PrintAdjacency(size_t IndexF)
    {
        std::cout<<"***Index F ***"<<IndexF<<std::endl;
        std::cout<<"Num edges : "<<NumVerts(IndexF)<<std::endl;
        std::cout<<"Num sides : "<<NumSides(IndexF)<<std::endl;

        size_t numE=MFaces[IndexF].AdjF.size();
        //assert((MFaces[IndexF].RealV[0]));
        if (MFaces[IndexF].RealV[0])
            std::cout<<"#";
        for (size_t i=0;i<MFaces[IndexF].AdjF.size();i++)
        {
            std::cout<<MFaces[IndexF].AdjF[i].first;
            if (!(MFaces[IndexF].RealV[(i+1)%numE]))
                std::cout<<",";
            else
                std::cout<<"#";
        }
        std::cout<<std::endl;
    }

    void MergeSide(size_t IndexF,size_t IndexSide,bool print_debug=false)
    {
        assert(IsRemovableSide(IndexF,IndexSide));

        //std::cout<<"0"<<std::endl;
        //        CheckMetaFaceAdj(IndexF);

        size_t IndexF0=IndexF;
        size_t IndexSide0=IndexSide;
        assert(IndexF0<MFaces.size());
        assert(IndexSide0<NumSides(IndexF0));

        //std::vector
        std::vector<std::pair<int,int> >  FaceSides;
        GetSideAdjacencySides(IndexF0,IndexSide0,FaceSides);

        //this is a condition for removability
        assert(FaceSides.size()==1);
        int IndexF1=FaceSides[0].first;
        int IndexSide1=FaceSides[0].second;

        assert((int)IndexF0!=IndexF1);

        if (print_debug)
            PrintAdjacency(IndexF0);

        if (print_debug)
            PrintAdjacency(IndexF1);

        assert(IndexF1>=0);
        assert(IndexF1<(int)MFaces.size());
        assert(IndexSide1<(int)NumSides(IndexF1));

        assert(!IsEmpty(IndexF0));
        assert(!IsEmpty(IndexF1));
        CheckMetaFaceAdj(IndexF0);
        CheckMetaFaceAdj(IndexF1);


        std::vector<size_t> MetaEdgeI0,MetaEdgeI1;
        GetSideMetaEdges(IndexF0,IndexSide0,MetaEdgeI0);
        GetSideMetaEdges(IndexF1,IndexSide1,MetaEdgeI1);
        assert(MetaEdgeI0.size()==MetaEdgeI1.size());

        size_t sizeE0=NumVerts(IndexF0);
        size_t sizeE1=NumVerts(IndexF1);



        //previous join together
        //size_t IndexPrev0=(MetaEdgeI0[0]+sizeE0-1)%sizeE0;
        size_t IndexE0=(MetaEdgeI0[0]);
        size_t IndexNext0=(MetaEdgeI0.back()+1)%sizeE0;
        //size_t IndexPrev1=(MetaEdgeI1[0]+sizeE1-1)%sizeE1;
        size_t IndexE1=(MetaEdgeI1[0]);
        size_t IndexNext1=(MetaEdgeI1.back()+1)%sizeE1;

        std::vector<size_t> NewV;
        std::vector<ScalarType> NewMetaL;
        std::vector<std::pair<int,int> > NewAdjF;
        std::vector<int> NewPathId;
        std::vector<bool> NewRealV;
        std::vector<std::pair<size_t,size_t> > NewVertDir;

        if (print_debug)
        {
            std::cout<<"IndexE0:"<<IndexE0<<std::endl;
            std::cout<<"IndexNext0:"<<IndexNext0<<std::endl;
            std::cout<<"IndexE1:"<<IndexE1<<std::endl;
            std::cout<<"IndexNext1:"<<IndexNext1<<std::endl;
        }

        for (size_t i=IndexNext0;i!=IndexE0;i=(i+1)%sizeE0)
        {
            NewV.push_back(MFaces[IndexF0].V[i]);
            NewMetaL.push_back(MFaces[IndexF0].MetaL[i]);
            NewAdjF.push_back(MFaces[IndexF0].AdjF[i]);

            //otherwise might create non manifold component
            assert(NewAdjF.back().first!=(int)IndexF0);
            assert(NewAdjF.back().first!=(int)IndexF1);

            //check not linking an empty face
            //assert(!IsEmpty(NewAdjF.back().first));

            NewPathId.push_back(MFaces[IndexF0].PathId[i]);
            NewRealV.push_back(MFaces[IndexF0].RealV[i]);
            NewVertDir.push_back(MFaces[IndexF0].VertDir[i]);
        }
        size_t Interval0=NewV.size();

        for (size_t i=IndexNext1;i!=IndexE1;i=(i+1)%sizeE1)
        {
            NewV.push_back(MFaces[IndexF1].V[i]);
            NewMetaL.push_back(MFaces[IndexF1].MetaL[i]);
            NewAdjF.push_back(MFaces[IndexF1].AdjF[i]);

            //otherwise might create non manifold component
            assert(NewAdjF.back().first!=(int)IndexF0);
            assert(NewAdjF.back().first!=(int)IndexF1);

            //check not linking an empty face
            //assert(!IsEmpty(NewAdjF.back().first));

            NewPathId.push_back(MFaces[IndexF1].PathId[i]);
            NewRealV.push_back(MFaces[IndexF1].RealV[i]);
            NewVertDir.push_back(MFaces[IndexF1].VertDir[i]);
        }

        for (size_t i=0;i<NewAdjF.size();i++)
        {
            if (NewAdjF[i].first==-1)continue;
            assert(NewAdjF[i].first!=(int)IndexF0);
            assert(NewAdjF[i].first!=(int)IndexF1);
            //check not linking an empty face
            assert(!IsEmpty(NewAdjF[i].first));
        }
        //size_t Interval1=NewV.size()-1;


        MFaces[IndexF0].V=NewV;
        MFaces[IndexF0].MetaL=NewMetaL;
        MFaces[IndexF0].AdjF=NewAdjF;
        MFaces[IndexF0].PathId=NewPathId;
        MFaces[IndexF0].RealV=NewRealV;
        MFaces[IndexF0].VertDir=NewVertDir;


//        MFaces[IndexF0].RealV[0]=false;
//        MFaces[IndexF0].RealV[Interval0]=false;
        MFaces[IndexF0].RealV[0]=IsFieldCornerFaceV(IndexF0,0);
        MFaces[IndexF0].RealV[Interval0]=IsFieldCornerFaceV(IndexF0,Interval0);

        MFaces[IndexF0].BaryF=CoordType(0,0,0);
        for (size_t i=0;i<MFaces[IndexF0].V.size();i++)
            MFaces[IndexF0].BaryF+=MetaVertPos(MFaces[IndexF0].V[i]);

        MFaces[IndexF0].BaryF/=MFaces[IndexF0].V.size();

        //        std::cout<<"*New Adjacency 0*"<<std::endl;
        //        PrintAdjacency(IndexF0);

        ReSortFace(IndexF0);

        //then update data on the other side
        MakeAdjacencyCoherentFrom(IndexF0);

        if (print_debug)
        {
            std::cout<<"*New Adjacency *"<<std::endl;
            PrintAdjacency(IndexF0);
        }

        //std::cout<<"C"<<std::endl;
        CheckMetaFaceAdj(IndexF0);
        //std::cout<<"D"<<std::endl;

        int ExpVal0=MFaces[IndexF0].ExpectedVal;
        int ExpVal1=MFaces[IndexF1].ExpectedVal;

        if ((ExpVal0==4)&&(ExpVal1==4))
            MFaces[IndexF0].ExpectedVal=4;
        else
            if ((ExpVal0!=4)&&(ExpVal1==4))
                MFaces[IndexF0].ExpectedVal=ExpVal0;
            else
                if ((ExpVal0==4)&&(ExpVal1!=4))
                    MFaces[IndexF0].ExpectedVal=ExpVal1;
                else
                {

                    assert((ExpVal0!=4)&&(ExpVal1!=4));
                    MFaces[IndexF0].ExpectedVal=-1;
                }

        MFaces[IndexF0].NumSing=MFaces[IndexF0].NumSing+MFaces[IndexF1].NumSing;

        //first make sure that all adj points to the old face
        for (size_t i=0;i<MFaces[IndexF1].AdjF.size();i++)
        {
            int OppF=MFaces[IndexF1].AdjF[i].first;
            if (OppF==-1)continue;
            int OppE=MFaces[IndexF1].AdjF[i].second;
            assert(OppE!=-1);

            //assert(OppF!=IndexF1);

            if (OppF!=(int)IndexF0)
            {
                //it should then be re-directed to the remaining face
                assert(!IsEmpty(OppF));
                assert(OppF<(int)MFaces.size());
                assert(OppE<(int)NumVerts(OppF));
                assert(MFaces[OppF].AdjF[OppE].first!=(int)IndexF1);
                assert(MFaces[OppF].AdjF[OppE].first!=-1);
                assert(!IsEmpty(MFaces[OppF].AdjF[OppE].first));

                if (print_debug)
                {
                    std::cout<<"Face: "<<OppF<<std::endl;
                    std::cout<<"Adjacency is: "<<MFaces[OppF].AdjF[OppE].first<<" instead of "<<IndexF0<<std::endl;
                }
                assert(MFaces[OppF].AdjF[OppE].first==(int)IndexF0);
            }
        }
        //assert(0);
        //clear the other face
        MFaces[IndexF1].V.clear();
        MFaces[IndexF1].MetaL.clear();
        MFaces[IndexF1].AdjF.clear();
        MFaces[IndexF1].PathId.clear();
        MFaces[IndexF1].RealV.clear();
        MFaces[IndexF1].VertDir.clear();

        //then remove the empty faces
        //RemoveEmptyFaces();
    }

    void InitExpectedValence(const std::vector<int> &FacePatches,
                             const std::vector<std::vector<size_t> > &FaceCorners)
    {
        //derive pert partition faces
        std::vector<std::vector<size_t> > PartitionFaces;
        PatchManager<MeshType>::DerivePerPartitionFaces(FacePatches,PartitionFaces);

        for (size_t i=0;i<MFaces.size();i++)
        {
            bool SingOnCorner;
            int ExpVal=PatchManager<MeshType>::ExpectedValence((*mesh),PartitionFaces[i],
                                                               FaceCorners[i],SingOnCorner);

            //tolerant wrt corner singularities
            if (!SingOnCorner)
                MFaces[i].ExpectedVal=ExpVal;
            else
                MFaces[i].ExpectedVal=NumSides(i);

            MFaces[i].NumSing=PatchManager<MeshType>::NumSingularities((*mesh),PartitionFaces[i],FaceCorners[i]);
        }
    }
    enum MergeRes{MRNoPath,MRNoRemovable,MRRemoved};

    MergeRes MergePathStep(const size_t IndexMetaF,const size_t IndexPath)
    {
        if (IsEmpty(IndexMetaF))
            return MRNoPath;

        int SideI=WhichSideFaceHasPath(IndexMetaF,IndexPath);
        if (SideI==-1)
            return MRNoPath;

        if (!IsRemovableSide(IndexMetaF,SideI))
            return MRNoRemovable;

        assert(SideI<(int)NumSides(IndexMetaF));
        MergeSide(IndexMetaF,SideI);
        return MRRemoved;
    }

    MergeRes MergePathIndex(const size_t IndexMetaF,const size_t IndexPath)
    {
        MergeRes currRes=MergePathStep(IndexMetaF,IndexPath);
        if (currRes==MRNoPath)return MRNoPath;
        if (currRes==MRNoRemovable)return MRNoRemovable;
        while (currRes==MRRemoved)
            currRes=MergePathStep(IndexMetaF,IndexPath);

        if (currRes==MRNoRemovable)return MRNoRemovable;
        //if has at least removed one and has not fall into unfeasibility
        return MRRemoved;
    }

    MergeRes MergePathIndex(const size_t IndexPath)
    {
        MergeRes res=MRNoPath;
        for (size_t i=0;i<MFaces.size();i++)
        {
            MergeRes currRes=MergePathIndex(i,IndexPath);
            if (currRes==MRNoRemovable)return MRNoRemovable;
            if (currRes==MRRemoved)
                res=MRRemoved;
        }
        return res;
    }

    bool RemoveIfPossible(size_t IndexPath,ScalarType Thr,
                          size_t MinVal,size_t MaxVal,
                          ScalarType CClarkability,ScalarType avgEdge,
                          bool match_valence)
    {
        //collect all the sides
        std::vector<std::pair<size_t,size_t> > FaceSide;
        //std::cout<<"0"<<std::endl;
        WhichSideFaceHasPath(IndexPath,FaceSide);
        //std::cout<<"1"<<std::endl;
        //if there is no side return false
        if (FaceSide.size()==0)return false;

        //then check if they are removable
        for (size_t i=0;i<FaceSide.size();i++)
        {
            size_t IndexF=FaceSide[i].first;
            size_t IndexS=FaceSide[i].second;
            assert(!IsEmpty(IndexF));
            if (!IsRemovableSide(IndexF,IndexS))return false;
        }
        //std::cout<<"2"<<std::endl;
        //get the old configuration
        std::vector<PatchInfo<ScalarType> > PatchInfos0;
        for (size_t i=0;i<FaceSide.size();i++)
        {
            if (IsEmpty(FaceSide[i].first))continue;
            PatchInfo<ScalarType> PInfo;
            GetPatchInfo(FaceSide[i].first,PInfo,Thr,match_valence);
            PatchInfos0.push_back(PInfo);
        }
        //std::cout<<"3"<<std::endl;
        //save old configuration
        std::vector<MetaFace> SwapMFaces=MFaces;

        //then perform the merge
        MergeRes MRes=MergePathIndex(IndexPath);

        if (MRes!=MRRemoved)
        {
            MFaces=SwapMFaces;
            return false;
        }
        //std::cout<<"4"<<std::endl;
        //std::cout<<"MERGED TEST RESULTS"<<std::endl;

        std::vector<PatchInfo<ScalarType> > PatchInfos1;
        for (size_t i=0;i<FaceSide.size();i++)
        {
            PatchInfo<ScalarType> PInfo;
            if (IsEmpty(FaceSide[i].first))continue;
            GetPatchInfo(FaceSide[i].first,PInfo,Thr,match_valence);
            PatchInfos1.push_back(PInfo);
        }
        //std::cout<<"5"<<std::endl;
        //check if has generated some infeaseable partition
        for (size_t i=0;i<PatchInfos1.size();i++)
        {
            if ((PatchInfos1[i].NumCorners<(int)MIN_ADMITTIBLE)
              ||(PatchInfos1[i].NumCorners>(int)MAX_ADMITTIBLE))
            {
                //restore
                MFaces=SwapMFaces;
                return false;
            }
        }


        bool CanRemove=true;
        CanRemove=PatchManager<MeshType>::BetterConfiguration(PatchInfos0,PatchInfos1,MinVal,
                                                              MaxVal,CClarkability,avgEdge,
                                                              match_valence,true,false);//,false,false);

        if (!CanRemove)
        {
            MFaces=SwapMFaces;
            return false;
        }
        //std::cout<<"6"<<std::endl;
        return true;
    }

    //    void CheckEdgeOrder()
    //    {
    //        for (size_t i=0;i<MFaces.size();i++)
    //                       for (size_t j=0;j<NumSides(i);j++)
    //                       {

    //                       }
    //    }

    size_t MergeLoopRun(ScalarType Thr,size_t MinVal,size_t MaxVal,
                        ScalarType CClarkability,ScalarType avgEdge,
                        bool match_valence)
    {
        //        size_t t0=clock();
        size_t num_removed=0;
        for (int i=MaxPath;i>=0;i--)
        {
            //std::cout<<"Removing Path "<<i<<std::endl;

            bool has_removed=RemoveIfPossible(i,Thr,MinVal,MaxVal,CClarkability,avgEdge,match_valence);
            if (has_removed)
                num_removed++;
        }
        return num_removed;
        //        size_t t1=clock();

        //        std::cout<<"Possible to remove "<<num_removed<<" paths"<<std::endl;
        //        std::cout<<"Time "<<((ScalarType)(t1-t0))/CLOCKS_PER_SEC<<" secs"<<std::endl;
    }

public:

    void GetRemainingPaths(std::vector<size_t> &RemPath)
    {
        RemPath.clear();
        for (size_t i=0;i<MFaces.size();i++)
        {
            if (IsEmpty(i))continue;
            for (size_t j=0;j<MFaces[i].PathId.size();j++)
                RemPath.push_back(MFaces[i].PathId[j]);
        }
        std::sort(RemPath.begin(),RemPath.end());
        auto last = std::unique(RemPath.begin(),RemPath.end());
        RemPath.erase(last,RemPath.end());
    }

    void MergeLoop(ScalarType Thr,size_t MinVal,size_t MaxVal,
                   ScalarType CClarkability,ScalarType avgEdge,
                   bool match_valence)
    {
        //size_t t0=clock();
        size_t num_removed=0;
        size_t curr_removed=0;
        size_t step_rem=0;
        do{
            curr_removed=MergeLoopRun(Thr,MinVal,MaxVal,CClarkability,avgEdge,match_valence);
            num_removed+=curr_removed;
            step_rem++;
        }while (curr_removed>0);
        //        size_t t1=clock();
        std::cout<<"Performed "<<step_rem<<" removal loops"<<std::endl;
        std::cout<<"Removed "<<num_removed<<" paths"<<std::endl;
     }

//    void GLDraw(const CEdgeMode EMode=CMRemovable,
//                bool DrawMetaEdgeAdj=true)
//    {
//        for (size_t i=0;i<MVert.size();i++)
//            GlDrawMetaVert(i);

//        for (size_t i=0;i<MFaces.size();i++)
//        {
//            GlDrawMetaFace(i,EMode);
//            GlDrawMetaFaceSides(i);

//            if (DrawMetaEdgeAdj)
//                GlDrawMetaFaceAdj(i);

//            GlDrawMetaFaceVert(i);
//            GlDrawMetaFaceVertField(i);
//        }
//        //        for (size_t i=0;i<MFaces.size();i++)
//        //            GlDrawMetaFace(i);
//    }

    void Clear()
    {
        MaxL=0;
        MaxPath=0;
        mesh=NULL;
        MVert.clear();
        MFaces.clear();
        MetaToMeshVert.clear();
        MeshToMetaVert.clear();
        MeshEdgeMap.clear();
        EDTableMetaMesh.Clear();
    }

    void Init(MeshType *_mesh,
              const std::vector<std::vector<size_t> > &FaceCorners,
              const std::vector<int> &FacePatches,
              const std::vector<std::vector<size_t > > &PathSeq,
              const std::vector<bool> &IsLoop,
              const std::vector<size_t> &FixedV,
              const std::map<std::pair<size_t,size_t>,ScalarType> &EdgeMap,
              const EdgeDirectionTable &EDTableMesh,
              const std::set<size_t> &_AvoidPartitions)
    {

        AvoidPartitions=_AvoidPartitions;
        std::cout<<"** FACES IGNORED IN METAMESH **:"<<AvoidPartitions.size()<<std::endl;

        Clear();

        mesh=_mesh;

        AddMetaVert(FaceCorners,FixedV);
        AddMetaFaces(FacePatches,FaceCorners,EdgeMap,EDTableMesh);

        UpdateAdjacency();

        MaxPath=PathSeq.size();

        std::cout<<"There are "<<PathSeq.size()<<" paths labels"<<std::endl;

        UpdatePathId(PathSeq,IsLoop);

        InitExpectedValence(FacePatches,FaceCorners);

        //MergeSideTest();
    }

    MetaMesh(){}
};
#endif
