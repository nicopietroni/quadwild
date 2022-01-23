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

#ifndef EDGE_DIRECTION_TABLE
#define EDGE_DIRECTION_TABLE

#include <map>
#include "vert_field_graph.h"
#include "vertex_classifier.h"

//struct EdgeVert
//{
//    size_t EV0;
//    size_t EV1;
//    size_t CurrV;

//    EdgeVert(size_t _EV0,size_t _EV1,size_t _CurrV)
//    {
//        EV0=std::min(_EV0,_EV1);
//        EV1=std::max(_EV0,_EV1);
//        CurrV=_CurrV;
//    }

//    inline bool operator ==(const EdgeVert &left)const
//    {
//        return ((EV0==left.EV0)&&
//                (EV1==left.EV1)&&
//                (CurrV==left.CurrV));
//    }

//    inline bool operator <(const EdgeVert &left)const
//    {
//        if ((EV0==left.EV0)&&
//                (EV1==left.EV1))
//            return (CurrV<left.CurrV);
//        if (EV0==left.EV0)
//            return (EV1<left.EV1);
//        return (EV0<left.EV0);
//    }
//};


//struct EdgeVertKeyHasher
//{
//    std::size_t operator()(const EdgeVert& k) const
//    {
//        //using std::size_t;
//        //using std::hash;
//        //using std::string;
//        const size_t _HASH_P0 = 73856093u;
//        const size_t _HASH_P1 = 19349663u;
//        const size_t _HASH_P2 = 83492791u;
//        return ((std::hash<size_t>()(k.EV0)*_HASH_P0)
//                ^ (std::hash<size_t>()(k.EV1) *_HASH_P1)
//                ^ (std::hash<size_t>()(k.CurrV) *_HASH_P2));
//    }
//};

//      size_t IndexF=Partitions[i][j];
//      for (size_t e=0;e<3;e++)
//      {
//          size_t IndexV=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[IndexF].V(e));
//          if (VertType[IndexV]==TVConvex)//this is convex then is a corner for sure
//              PartitionCorners[i].push_back(IndexV);
//          else
//          {
//              bool PossibleCorner=(HasOrthogonalCross(DirVert[IndexV])||
//                                   HasNarrowCross(DirVert[IndexV]));
//              if (!PossibleCorner)continue;

//              //check incomplete concave
//              if (VertType[IndexV]==TVConcave)
//              {
//                  std::vector<FaceType*> StarF;
//                  std::vector<int> LocalIndexV;
//                  vcg::face::VFStarVF(&VFGraph.Mesh().vert[IndexV],StarF,LocalIndexV);
//                  ScalarType Angle=0;
//                  for (size_t j=0;j<StarF.size();j++)
//                  {
//                      FaceType* f=StarF[j];
//                      int currV=LocalIndexV[j];
//                      if (f->Q()!=i)continue;
//                      Angle+=vcg::face::WedgeAngleRad(*f,currV);
//                  }
//                  if (Angle>(M_PI))continue;
//              }
//              PartitionCorners[i].push_back(IndexV);
//          }
//      }
//  }
//  std::sort(PartitionCorners[i].begin(),PartitionCorners[i].end());
//  std::vector<size_t>::iterator it;
//  it = std::unique (PartitionCorners[i].begin(),PartitionCorners[i].end());
//  PartitionCorners[i].resize( std::distance(PartitionCorners[i].begin(),it) );

bool HasOrthogonalCross(const std::vector<size_t> &Directions)
{
    if (Directions.size()<2)return false;
    for (size_t i=0;i<Directions.size()-1;i++)
        for (size_t j=(i+1);j<Directions.size();j++)
        {
            size_t Dir0=Directions[i];
            size_t Dir1=Directions[j];
            if ((Dir0 % 2)!=(Dir1 % 2))return true;
        }
    return false;
}

bool HasNarrowCross(const std::vector<size_t> &Directions)
{
    if (Directions.size()<2)return false;
    for (size_t i=0;i<Directions.size()-1;i++)
        for (size_t j=(i+1);j<Directions.size();j++)
        {
            size_t Dir0=Directions[i];
            size_t Dir1=Directions[j];
            if (Dir0 == Dir1)return true;
        }
    return false;
}

class EdgeDirectionTable
{


    bool PossibleCross(const std::vector<size_t> &Directions)const
    {
        return (HasOrthogonalCross(Directions)||
                HasNarrowCross(Directions));
    }

    std::vector<TypeVert> VertType;
    std::map<EdgeVert,int> EdgeMapDir;
    //std::vector<std::vector<size_t> > VertDirections;

public:

    //    std::set<size_t> ConvexV;
    //    std::set<size_t> ConcaveV;

    //std::vector<std::vector<size_t> > EdgeDirVert;

    template <class MeshType>
    void FindPerVertDirs(const VertexFieldGraph<MeshType> &VFGraph,
                         const std::vector<size_t> &Partition,
                         std::map<EdgeVert,size_t> &EdgeDirVert,
                         //std::unordered_map<EdgeVert,size_t,EdgeVertKeyHasher> &EdgeDirVert,
                         std::vector<std::vector<size_t> > &DirVert)
    {
        typedef typename MeshType::CoordType CoordType;

        DirVert.clear();
        DirVert.resize(VFGraph.Mesh().vert.size());
        std::set<std::pair<CoordType,CoordType> > AddedEdges;
        for (size_t i=0;i<Partition.size();i++)
        {
            size_t IndexF=Partition[i];
            for (size_t e=0;e<3;e++)
            {
                size_t IndexV0=vcg::tri::Index(VFGraph.Mesh(),
                                               VFGraph.Mesh().face[IndexF].V0(e));

                size_t IndexV1=vcg::tri::Index(VFGraph.Mesh(),
                                               VFGraph.Mesh().face[IndexF].V1(e));

                CoordType Pos0=VFGraph.Mesh().vert[IndexV0].P();
                CoordType Pos1=VFGraph.Mesh().vert[IndexV1].P();
                std::pair<CoordType,CoordType> KeyE(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                if (AddedEdges.count(KeyE)>0)continue;
                AddedEdges.insert(KeyE);


                size_t MinV=std::min(IndexV0,IndexV1);
                size_t MaxV=std::max(IndexV0,IndexV1);

                EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
                if (EdgeDirVert.count(EdgeKey0)>0)//||(VFGraph.EdgeBorderDir.count(EdgeKey0)>0))
                {
                    size_t EdgeDir=EdgeDirVert[EdgeKey0];
                    DirVert[IndexV0].push_back(EdgeDir);
                }


                EdgeVert EdgeKey1(MinV,MaxV,IndexV1);
                if (EdgeDirVert.count(EdgeKey1)>0)//||(VFGraph.EdgeBorderDir.count(EdgeKey1)>0))
                {
                    size_t EdgeDir=EdgeDirVert[EdgeKey1];
                    DirVert[IndexV1].push_back(EdgeDir);
                }
            }
        }
    }

    void FindPerVertDirs(const std::vector<std::pair<size_t,size_t> > &BorderEdges,
                         std::map<size_t,std::vector<size_t> > &DirVert)const
    {
        DirVert.clear();
        for (size_t i=0;i<BorderEdges.size();i++)
        {
            size_t IndexV0=std::min(BorderEdges[i].first,BorderEdges[i].second);
            size_t IndexV1=std::max(BorderEdges[i].first,BorderEdges[i].second);
            assert(IndexV0!=IndexV1);

            //std::pair<size_t,size_t> key(IndexV0,IndexV1);
            EdgeVert EV0(IndexV0,IndexV1,IndexV0);
            EdgeVert EV1(IndexV0,IndexV1,IndexV1);

            //            if (PossibleCross(VertDirections[IndexV0]))
            //            {
            //                if (EdgeMapDir.count(EV0)>0)
            //                     DirVert[IndexV0].push_back(EdgeMapDir.at(EV0));
            //            }
            //            if (PossibleCross(VertDirections[IndexV1]))
            //            {
            //                if (EdgeMapDir.count(EV1)>0)
            //                     DirVert[IndexV1].push_back(EdgeMapDir.at(EV1));
            //            }
            if (!EdgeMapDir.count(EV0))continue;
            assert(EdgeMapDir.count(EV1)>0);

            DirVert[IndexV0].push_back(EdgeMapDir.at(EV0));
            DirVert[IndexV1].push_back(EdgeMapDir.at(EV1));
        }
    }

public:

    bool HasEdgeDir(const EdgeVert &EV)const
    {
        return (EdgeMapDir.count(EV)>0);
    }

    size_t GetEdgeDir(const EdgeVert &EV)const
    {
        assert(EdgeMapDir.count(EV)>0);
        return(EdgeMapDir.at(EV));
    }

    void AddEdgeDir(const EdgeVert &EV,const size_t Dir)
    {
        //assert(EdgeMapDir.count(EV)==0);
        EdgeMapDir[EV]=Dir;
        //VertDirections[EV.CurrV].push_back(Dir);
    }

    void RemoveEdgeDir(const EdgeVert &EV)
    {
        //        size_t Dir=EdgeMapDir[EV];
        EdgeMapDir.erase(EV);
        //        size_t IndexV=EV.CurrV;
        //        for (size_t i=0;i<VertDirections[IndexV].size();i++)
        //            if (VertDirections[IndexV][i]==Dir)
        //            {
        //                VertDirections[IndexV].erase(VertDirections[IndexV].begin()+i);
        //                return;
        //            }
        //should have been found
        //        assert(0);
    }
    //    void RemoveNodeVert()
    //    {

    //    }

    //    void PrintDifferences(const EdgeDirectionTable &ED1)const
    //    {
    //        std::cout<<"**TESTING DIFFERENCE**"<<std::endl;
    //        if (ConvexV!=ED1.ConvexV)
    //        {
    //            std::cout<<"Different ConvexV"<<std::endl;
    //            std::cout<<"size 0:"<<ConvexV.size()<<std::endl;
    //            std::cout<<"size 1:"<<ED1.ConvexV.size()<<std::endl;
    //            assert(0);
    //        }
    //        if (ConcaveV!=ED1.ConcaveV)
    //        {
    //            std::cout<<"Different ConcaveV"<<std::endl;
    //            assert(0);
    //        }
    //        if (EdgeMapDir!=ED1.EdgeMapDir)
    //        {
    //            std::cout<<"Different EdgeDirMap"<<std::endl;
    //            assert(0);
    //        }
    //    }

    void FindCorners(const std::vector<std::pair<size_t,size_t> > &BorderEdges,
                     std::vector<size_t> &PartitionCorners)const
    {
        PartitionCorners.clear();

        //first check if there is a corner
        if (VertType.size()>0)
        {
            for (size_t i=0;i<BorderEdges.size();i++)
            {
                size_t IndexV0=BorderEdges[i].first;
                size_t IndexV1=BorderEdges[i].second;
                if (VertType[IndexV0]==TVConvex)
                    PartitionCorners.push_back(IndexV0);
                if (VertType[IndexV1]==TVConvex)
                    PartitionCorners.push_back(IndexV1);
            }
        }

        //cumulate per vertex direction
        std::map<size_t,std::vector<size_t> >  DirVert;
        FindPerVertDirs(BorderEdges,DirVert);

        //then iterate over all
        std::map<size_t,std::vector<size_t> > ::iterator IteMap;
        for (IteMap=DirVert.begin();IteMap!=DirVert.end();IteMap++)
        {
            bool PossibleCorner0=HasOrthogonalCross((*IteMap).second);
            bool PossibleCorner1=HasNarrowCross((*IteMap).second);
            if (!(PossibleCorner0 || PossibleCorner1))continue;
            size_t IndexV=(*IteMap).first;
            PartitionCorners.push_back(IndexV);
        }
        std::sort(PartitionCorners.begin(),PartitionCorners.end());
        std::vector<size_t>::iterator it;
        it = std::unique (PartitionCorners.begin(),PartitionCorners.end());
        PartitionCorners.resize( std::distance(PartitionCorners.begin(),it) );
    }

    template <class ScalarType>
    void FindCorners(const std::vector<std::pair<size_t,size_t> > &BorderEdges,
                     const std::map<size_t,ScalarType> &AngleBorders,
                     std::vector<size_t> &PartitionCorners)const
    {
        FindCorners(BorderEdges,PartitionCorners);
        std::vector<size_t> PartitionCornersSwap;
        bool filtered=false;
        for (size_t i=0;i<PartitionCorners.size();i++)
        {
            size_t IndexV=PartitionCorners[i];
            //            if (VertType[IndexV]==TVNarrow)
            //            {
            //                filtered=true;
            //                continue;
            //            }
            if (((VertType[IndexV]==TVConcave)||
                 (VertType[IndexV]==TVNarrow))&&
                    (AngleBorders.at(IndexV)>(M_PI)))
            {
                filtered=true;
                continue;
            }

            PartitionCornersSwap.push_back(IndexV);
        }
        if (filtered)
            PartitionCorners=PartitionCornersSwap;

    }

//    template <class ScalarType>
//    void FindCorners(const std::vector<std::pair<size_t,size_t> > &BorderEdges,
//                     const std::map<size_t,ScalarType> &AngleBorders,
//                     std::vector<size_t> &PartitionCorners)const
//    {
//        FindCorners(BorderEdges,PartitionCorners);
//        std::vector<size_t> PartitionCornersSwap;
//        bool filtered=false;
//        for (size_t i=0;i<PartitionCorners.size();i++)
//        {
//            size_t IndexV=PartitionCorners[i];
//            //            if (VertType[IndexV]==TVNarrow)
//            //            {
//            //                filtered=true;
//            //                continue;
//            //            }
//            if (((VertType[IndexV]==TVConcave)||
//                 (VertType[IndexV]==TVNarrow))&&
//                    (AngleBorders.at(IndexV)>(M_PI)))
//            {
//                filtered=true;
//                continue;
//            }

//            PartitionCornersSwap.push_back(IndexV);
//        }
//        if (filtered)
//            PartitionCorners=PartitionCornersSwap;

//    }

    void FindPossibleCorners(std::vector<size_t> &PartitionCorners)
    {
        PartitionCorners.clear();
        std::map<EdgeVert,int>::iterator IteMap;
        //store all possible borders
        std::vector<std::pair<size_t,size_t> > BorderEdges;
        for (IteMap=EdgeMapDir.begin();IteMap!=EdgeMapDir.end();IteMap++)
        {
            size_t IndexV0=(*IteMap).first.EV0;
            size_t IndexV1=(*IteMap).first.EV1;
            std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            BorderEdges.push_back(Key);
        }
        std::sort(BorderEdges.begin(),BorderEdges.end());
        std::vector<std::pair<size_t,size_t> >::iterator it;
        it = std::unique (BorderEdges.begin(),BorderEdges.end());
        BorderEdges.resize( std::distance(BorderEdges.begin(),it) );
        FindCorners(BorderEdges,PartitionCorners);
    }

    void Clear()
    {
        VertType.clear();
        EdgeMapDir.clear();
    }

    void Init(std::vector<TypeVert> &_VertType)
    {
        VertType=_VertType;
        //        VertDirections.clear();
        //        VertDirections.resize(VertType.size());
        //        ConvexV=std::set<size_t>(ConvexIdx.begin(),ConvexIdx.end());
        //        ConcaveV=std::set<size_t>(ConcaveIdx.begin(),ConcaveIdx.end());
        //std::cout<<"There are "<<ConvexV.size()<<" convex vert"<<std::endl;
        EdgeMapDir.clear();
    }
};


//template <class MeshType>
//void FindPartitionEdges(MeshType &mesh,const std::vector<size_t> &Partition,
//                     std::vector<std::pair<size_t,size_t> > &BorderEdges)
//{
//    typedef typename MeshType::FaceType FaceType;
//    typedef typename MeshType::CoordType CoordType;

//    //int t3=clock();
//    //then collect all border patches
//    std::set<std::pair<CoordType,CoordType> > AddedEdges;
//    BorderEdges.clear();

//    //    vcg::tri::UnMarkAll(mesh);
//    //    for (size_t i=0;i<Partition.size();i++)
//    //    {
//    //        size_t IndexF=Partition[i];
//    //        vcg::tri::Mark(mesh,&mesh.face[IndexF]);
//    //    }

//    for (size_t i=0;i<Partition.size();i++)
//    {
//        size_t IndexF=Partition[i];
//        for (size_t j=0;j<mesh.face[IndexF].VN();j++)
//        {
//            //assert(vcg::tri::IsMarked(mesh,&mesh.face[IndexF]));
//            size_t IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].cV0(j));
//            size_t IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].cV1(j));

//            //position check for multiple border across sharp features
//            CoordType Pos0=mesh.vert[IndexV0].P();
//            CoordType Pos1=mesh.vert[IndexV1].P();
//            std::pair<CoordType,CoordType> KeyE(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
//            if (AddedEdges.count(KeyE)>0)continue;
//            AddedEdges.insert(KeyE);

//            //            FaceType *Fopp=mesh.face[IndexF].cFFp(j);
//            //            bool IsBorder=false;
//            //            if (Fopp==(&mesh.face[IndexF]))
//            //                IsBorder=true;
//            //            if (!vcg::tri::IsMarked(mesh,Fopp))
//            //                IsBorder=true;

//            //            if (IsBorder)
//            std::pair<size_t,size_t> pairV(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
//            BorderEdges.push_back(pairV);
//        }
//    }
//    std::sort(BorderEdges.begin(),BorderEdges.end());
//    std::vector<size_t>::iterator it;
//    it = std::unique (BorderEdges.begin(),BorderEdges.end());
//}

template <class MeshType>
void FindPartitionEdges(MeshType &mesh,
                        const std::vector<size_t> &Partition,
                        std::vector<std::pair<size_t,size_t> > &BorderEdges)
{
    typedef typename MeshType::FaceType FaceType;
    //typedef typename MeshType::CoordType CoordType;

    //int t3=clock();
    //then collect all border patches
    BorderEdges.clear();

    //std::set<std::pair<CoordType,CoordType> > AddedEdges;
    std::set<std::pair<size_t,size_t> > AddedEdges;

    for (size_t i=0;i<Partition.size();i++)
    {
        size_t IndexF=Partition[i];
        for (int j=0;j<mesh.face[IndexF].VN();j++)
        {
            FaceType *f=&mesh.face[IndexF];
            FaceType *fOpp=f->FFp(j);
            int IOpp=f->FFi(j);
            bool isEdgeS=false;
            bool onBorder=(vcg::face::IsBorder(*f,j));
            if (!onBorder)
                isEdgeS=(f->IsFaceEdgeS(j)||fOpp->IsFaceEdgeS(IOpp));

            if (onBorder || isEdgeS)
            {
                //assert(vcg::tri::IsMarked(mesh,&mesh.face[IndexF]));
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[IndexF].cV0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[IndexF].cV1(j));

                //                CoordType Pos0=mesh.vert[IndexV0].P();
                //                CoordType Pos1=mesh.vert[IndexV1].P();
                //                std::pair<CoordType,CoordType> KeyE(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                //                if (AddedEdges.count(KeyE)>0)continue;
                //                AddedEdges.insert(KeyE);


                std::pair<size_t,size_t> pairV(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                if (AddedEdges.count(pairV)>0)continue;
                AddedEdges.insert(pairV);

                BorderEdges.push_back(pairV);
            }
        }
    }
    //    std::sort(BorderEdges.begin(),BorderEdges.end());
    //    std::vector<std::pair<size_t,size_t> >::iterator it;
    //    it = std::unique (BorderEdges.begin(),BorderEdges.end());
}

template <class MeshType>
void InitAngleBorders(MeshType &mesh,const std::vector<size_t> &Partition,
                      std::map<size_t,typename MeshType::ScalarType> &AngleBorders)
{
    typedef typename MeshType::ScalarType ScalarType;

    AngleBorders.clear();
    for (size_t i=0;i<Partition.size();i++)
    {
        size_t IndexF=Partition[i];
        for (int j=0;j<mesh.face[IndexF].VN();j++)
        {
            size_t IndexV=vcg::tri::Index(mesh,mesh.face[IndexF].cV(j));
            ScalarType CurrAngle=vcg::face::WedgeAngleRad(mesh.face[IndexF],j);
            AngleBorders[IndexV]+=CurrAngle;
        }
    }
}

template <class MeshType>
void FindCorners(const EdgeDirectionTable &EDirTable,
                 MeshType &mesh,
                 const std::vector<size_t> &Partition,
                 std::vector<size_t> &PartitionCorners)
{
    typedef typename MeshType::ScalarType ScalarType;

    std::vector<std::pair<size_t,size_t> > PartitionEdges;
    std::map<size_t,ScalarType> AngleBorders;
    //assert(0);
    InitAngleBorders(mesh,Partition,AngleBorders);
    FindPartitionEdges(mesh,Partition,PartitionEdges);
    EDirTable.FindCorners<ScalarType>(PartitionEdges,AngleBorders,PartitionCorners);
}

template <class MeshType>
bool HasConcaveEdge(const EdgeDirectionTable &EDirTable,
                    MeshType &mesh,
                    const std::vector<size_t> &Partition)
{
    //TO BE IMPLEMENTED
}

template <class MeshType>
void FindCorners(const EdgeDirectionTable &EDirTable,MeshType &mesh,
                 const std::vector<std::vector<size_t> > &Partitions,
                 std::vector<std::vector<size_t> > &PartitionCorners)
{
    PartitionCorners.resize(Partitions.size());
    for (size_t i=0;i<Partitions.size();i++)
        FindCorners(EDirTable,mesh,Partitions[i],PartitionCorners[i]);
}

template <class MeshType>
void AddEdgeNodes(const std::vector<size_t> &Nodes,
                  const bool IsLoop,
                  EdgeDirectionTable &EDirTable)
{
    size_t Limit=Nodes.size()-1;
    if (IsLoop)
        Limit++;
    size_t size=Nodes.size();

    for (size_t j=0;j<Limit;j++)
    {
        size_t IndexN0=Nodes[j];
        size_t IndexN1=Nodes[(j+1)%size];
        size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
        size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);

        assert(IndexV0!=IndexV1);
        size_t MinV=std::min(IndexV0,IndexV1);
        size_t MaxV=std::max(IndexV0,IndexV1);

        size_t DirV0=VertexFieldGraph<MeshType>::NodeDirI(IndexN0);
        size_t DirV1=VertexFieldGraph<MeshType>::NodeDirI(IndexN1);
        EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
        EdgeVert EdgeKey1(MinV,MaxV,IndexV1);

        EDirTable.AddEdgeDir(EdgeKey0,DirV0);
        EDirTable.AddEdgeDir(EdgeKey1,(DirV1+2)%4);
        //        assert(EDirTable.EdgeMapDir.count(EdgeKey0)==0);
        //        assert(EDirTable.EdgeMapDir.count(EdgeKey1)==0);

        //        EDirTable.EdgeMapDir[EdgeKey0]=DirV0;
        //        EDirTable.EdgeMapDir[EdgeKey1]=((DirV1+2)%4);//put the inverse cause look internally the interval

    }
}

template <class MeshType>
void RemoveEdgeNodes(const std::vector<size_t> &Nodes,
                     const bool IsLoop,
                     EdgeDirectionTable &EDirTable)
{
    size_t Limit=Nodes.size()-1;
    if (IsLoop)
        Limit++;
    size_t size=Nodes.size();

    for (size_t j=0;j<Limit;j++)
    {
        size_t IndexN0=Nodes[j];
        size_t IndexN1=Nodes[(j+1)%size];
        size_t IndexV0=VertexFieldGraph<MeshType>::NodeVertI(IndexN0);
        size_t IndexV1=VertexFieldGraph<MeshType>::NodeVertI(IndexN1);

        assert(IndexV0!=IndexV1);
        size_t MinV=std::min(IndexV0,IndexV1);
        size_t MaxV=std::max(IndexV0,IndexV1);

        //        size_t DirV0=VertexFieldGraph<MeshType>::NodeDirI(IndexN0);
        //        size_t DirV1=VertexFieldGraph<MeshType>::NodeDirI(IndexN1);
        EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
        EdgeVert EdgeKey1(MinV,MaxV,IndexV1);

        //        assert(EDirTable.EdgeMapDir.count(EdgeKey0)==0);
        //        assert(EDirTable.EdgeMapDir.count(EdgeKey1)==0);

        //EDirTable.EdgeMapDir.erase(EdgeKey0);
        //EDirTable.EdgeMapDir.erase(EdgeKey1);

        EDirTable.RemoveEdgeDir(EdgeKey0);
        EDirTable.RemoveEdgeDir(EdgeKey1);

        //        EDirTable.EdgeMapDir[EdgeKey0]=DirV0;
        //        EDirTable.EdgeMapDir[EdgeKey1]=((DirV1+2)%4);//put the inverse cause look internally the interval

    }
}

template <class MeshType>
void AddEdgeNodes(const std::vector<std::vector<size_t> > &Nodes,
                  const std::vector<bool> &IsLoop,
                  EdgeDirectionTable &EDirTable)
{
    for (size_t i=0;i<Nodes.size();i++)
        AddEdgeNodes<MeshType>(Nodes[i],IsLoop[i],EDirTable);
}


template <class MeshType>
void AddBorder(const VertexFieldGraph<MeshType> &VFGraph,
               EdgeDirectionTable &EDirTable)
{
    //do the same for borders
    for (size_t i=0;i<VFGraph.Mesh().face.size();i++)
        for (int j=0;j<VFGraph.Mesh().face[i].VN();j++)
        {
            if (!vcg::face::IsBorder(VFGraph.Mesh().face[i],j))continue;
            size_t IndexV0=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[i].cV0(j));
            size_t IndexV1=vcg::tri::Index(VFGraph.Mesh(),VFGraph.Mesh().face[i].cV1(j));
            size_t DirFlatV0,DirFlatV1;
            VFGraph.GetEdgeDir(IndexV0,IndexV1,DirFlatV0,DirFlatV1);

            assert(IndexV0!=IndexV1);
            size_t MinV=std::min(IndexV0,IndexV1);
            size_t MaxV=std::max(IndexV0,IndexV1);

            EdgeVert EdgeKey0(MinV,MaxV,IndexV0);
            EdgeVert EdgeKey1(MinV,MaxV,IndexV1);


            EDirTable.AddEdgeDir(EdgeKey0,DirFlatV0);
            EDirTable.AddEdgeDir(EdgeKey1,DirFlatV1);
        }
}


#endif
