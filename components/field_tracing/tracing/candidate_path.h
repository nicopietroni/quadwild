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

#ifndef CANDIDATE_PATH
#define CANDIDATE_PATH


#include "vert_field_graph.h"
#include "graph_query.h"
#include "vertex_classifier.h"


enum TraceType{TraceDirect,DijkstraReceivers,TraceLoop};

struct CandidateTrace
{
    TypeVert FromType;
    TypeVert ToType;
    TraceType TracingMethod;
    size_t InitNode;
    std::vector<size_t> PathNodes;
    bool IsLoop;
    bool Updated;
    float Priority;
    bool Unremovable;

    inline bool operator <(const CandidateTrace &C1)const
    {return (Priority<C1.Priority);}

    CandidateTrace(){}

    CandidateTrace(TypeVert _FromType,
                   TypeVert _ToType,
                   TraceType _TracingMethod,
                   size_t _InitNode)
    {
        FromType=_FromType;
        ToType=_ToType;
        TracingMethod=_TracingMethod;
        InitNode=_InitNode;
        IsLoop=false;
        Updated=false;
        Unremovable=false;
        Priority=0;
    }
};

void GetCandidateNodes(const std::vector<CandidateTrace> &TraceSet,
                       std::vector<size_t> &ChoosenNodes)
{
    ChoosenNodes.clear();
    for (size_t i=0;i<TraceSet.size();i++)
        for (size_t j=0;j<TraceSet[i].PathNodes.size();j++)
        {
            size_t IndexN=TraceSet[i].PathNodes[j];
            ChoosenNodes.push_back(IndexN);
        }
}

void GetCandidateNodes(const CandidateTrace &CurrTrace,
                       std::vector<size_t> &ChoosenNodes)
{
    for (size_t j=0;j<CurrTrace.PathNodes.size();j++)
    {
        size_t IndexN=CurrTrace.PathNodes[j];
        ChoosenNodes.push_back(IndexN);
    }
}

template <class MeshType>
void GetCandidateTangentNodes(const CandidateTrace &CurrTrace,
                              std::vector<size_t> &ChoosenTangNodes)
{
    if (CurrTrace.IsLoop)
    {
        for (size_t j=0;j<CurrTrace.PathNodes.size();j++)
        {
            size_t IndexN=VertexFieldGraph<MeshType>::TangentNode(CurrTrace.PathNodes[j]);
            ChoosenTangNodes.push_back(IndexN);
        }
    }
    else
    {
        for (size_t j=1;j<CurrTrace.PathNodes.size()-1;j++)
        {
            size_t IndexN=VertexFieldGraph<MeshType>::TangentNode(CurrTrace.PathNodes[j]);
            ChoosenTangNodes.push_back(IndexN);
        }
    }
}

//template <class MeshType>
//void GetCandidateNodesNodesAndTangent(const std::vector<CandidateTrace> &TraceSet,
//                                      std::vector<size_t> &ChoosenNodes)
//{
//    GetCandidateNodes(TraceSet,ChoosenNodes);
//    std::vector<size_t> TangentNodes=ChoosenNodes;
//    VertexFieldGraph<MeshType>::TangentNodes(TangentNodes);
//    ChoosenNodes.insert(ChoosenNodes.end(),TangentNodes.begin(),TangentNodes.end());
//}

template <class MeshType>
void GetCandidateNodesNodesAndTangent(const std::vector<CandidateTrace> &TraceSet,
                                      std::vector<size_t> &ChoosenNodes)
{
    GetCandidateNodes(TraceSet,ChoosenNodes);
    for (size_t i=0;i<TraceSet.size();i++)
    {
        std::vector<size_t> TangentNodes;
        GetCandidateTangentNodes<MeshType>(TraceSet[i],TangentNodes);
        ChoosenNodes.insert(ChoosenNodes.end(),TangentNodes.begin(),TangentNodes.end());
    }
}

template <class MeshType>
void GetPathPos(VertexFieldGraph<MeshType> &VFGraph,
                const std::vector<CandidateTrace> &TraceSet,
                std::vector<std::vector<vcg::face::Pos<typename MeshType::FaceType> > > &Paths)
{
    //typedef typename MeshType::FaceType FaceType;

    Paths.clear();
    Paths.resize(TraceSet.size());
    for (size_t i=0;i<TraceSet.size();i++)
    {
        VFGraph.GetNodesPos(TraceSet[i].PathNodes,
                            TraceSet[i].IsLoop,
                            Paths[i]);
    }
}

template <class MeshType>
bool CollideCandidates(VertexFieldGraph<MeshType> &VFGraph,
                       const CandidateTrace &CT0,
                       const CandidateTrace &CT1)
{
    return (VertexFieldQuery<MeshType>::CollideTraces(VFGraph,CT0.PathNodes,CT0.IsLoop,
                                                      CT1.PathNodes,CT1.IsLoop));
}

template <class MeshType>
bool CollideWithCandidateSet(VertexFieldGraph<MeshType> &VFGraph,
                             const CandidateTrace &TestTrace,
                             const std::vector<CandidateTrace> &TraceSet)
{
    for (size_t i=0;i<TraceSet.size();i++)
    {
        bool collide=CollideCandidates<MeshType>(VFGraph,TestTrace,TraceSet[i]);
        if (collide)return true;
    }
    return false;
}

template <class MeshType>
bool UpdateCandidate(VertexFieldGraph<MeshType> &VFGraph,
                     CandidateTrace &ToUpdate,
                     const typename MeshType::ScalarType &Drift,
                     const typename MeshType::ScalarType &MaxDijstraDist)//,
                     //bool DebugMsg)
{
    assert(!ToUpdate.Updated);
    ToUpdate.Updated=true;

    size_t IndexN0=ToUpdate.InitNode;
    assert(VFGraph.IsActive(IndexN0));

    if (ToUpdate.TracingMethod==TraceDirect)
    {
        std::vector<size_t> PathN;
        bool hasTraced=TraceDirectPath(VFGraph,IndexN0,PathN);//,DebugMsg);
        if (!hasTraced)return false;
        ToUpdate.PathNodes=PathN;
        ToUpdate.IsLoop=false;
        return true;
    }
    if (ToUpdate.TracingMethod==DijkstraReceivers)
    {
        std::vector<size_t> PathN;
        bool hasTraced=TraceDijkstraPath(VFGraph,IndexN0,Drift,MaxDijstraDist,PathN);
        if (!hasTraced)return false;
        ToUpdate.PathNodes=PathN;
        ToUpdate.IsLoop=false;
        return true;
    }
    if (ToUpdate.TracingMethod==TraceLoop)
    {
        std::vector<size_t> PathN;
        bool hasTraced=TraceLoopPath(VFGraph,IndexN0,Drift,PathN);
        if (!hasTraced)return false;
        ToUpdate.PathNodes=PathN;
        ToUpdate.IsLoop=true;
        return true;
    }
    return false;
}



#endif
