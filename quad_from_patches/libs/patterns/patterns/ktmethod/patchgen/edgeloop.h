//#include "Pattern.h"
#ifndef EDGELOOP_H
#define EDGELOOP_H

#include <vcg/complex/algorithms/create/platonic.h>
#include <set>
//just add topology
template <class MeshType>
void insert_edgeloop(MeshType& mesh, vcg::face::Pos<typename MeshType::FaceType>& h_start)
{
    //store corners
    vector<int> corners;
    for(typename MeshType::VertexIterator vi=mesh.vert.begin();vi!=mesh.vert.end();vi++){
        if(vi->IsS()){
            corners.push_back(vcg::tri::Index(mesh,&*vi));
        }
    }

    UpdateTopology<MeshType>::FaceFace(mesh);
    UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
    UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

    // check if this edge loop forms a loop by walking toward h1 of e
    vcg::face::Pos<typename MeshType::FaceType> h=h_start;
    bool isloop=false;
    int count=0;    
    while(count<1000){
        vcg::face::Pos<typename MeshType::FaceType> h_next_next= h;
        h_next_next.FlipE();
        h_next_next.FlipV();
        h_next_next.FlipE();
        h_next_next.FlipF();
        if(h.IsBorder()){
            // reached a boundary
            h_start = h;
            break;
        }
        else if(h_start==h_next_next){
            // loop detected!
            isloop=true;
            break;
        }
        count++;
    }    
    if(count==1000)
        return;

   // for each edge, insert a vertex at its "middle"
    h=h_start;
    vcg::face::Pos<typename MeshType::FaceType> h_next_next;
    assert(h.V()->IsB());
    vector<std::pair<int,int>> vertices;
   do{
       h_next_next= h;
       //cout<<"index h.v "<<vcg::tri::Index(mesh,h.V())<<endl;
       //cout<<"index h.v from "<<vcg::tri::Index(mesh,h.VFlip())<<endl;
       h_next_next.FlipE();
       //cout<<"after flipE "<<vcg::tri::Index(mesh,h_next_next.VFlip())<<endl;
       h_next_next.FlipV();
       //cout<<"after flipV "<<vcg::tri::Index(mesh,h_next_next.VFlip())<<endl;
       h_next_next.FlipE();
       //cout<<"after flipE "<<vcg::tri::Index(mesh,h_next_next.VFlip())<<endl;
       h_next_next.FlipF();
       //cout<<"after flipF from "<<vcg::tri::Index(mesh,h_next_next.VFlip())<<endl;
       //cout<<"after flipF to "<<vcg::tri::Index(mesh,h_next_next.V())<<endl;
       vertices.push_back(make_pair(vcg::tri::Index(mesh,h.V()),vcg::tri::Index(mesh,h.VFlip())));
       h.F()->SetS(); //to delete after
       if (h.IsBorder() && h!=h_start)
           break;
       if (h_next_next == h_start)
           break;

       h=h_next_next;
       //cout<<"index final h.v "<<vcg::tri::Index(mesh,h.V())<<endl;
       //cout<<"index final h.v from"<<vcg::tri::Index(mesh,h.VFlip())<<endl;
   }while(true);   

    //deleting old faces
    for(typename MeshType::FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
        if(fi->IsS())
            vcg::tri::Allocator<MeshType>::DeleteFace(mesh,*fi);
    }
   vcg::tri::Allocator<MeshType>::CompactFaceVector(mesh);
   //cout<<"number vertices "<<mesh.vert.size();

   // add new vertices
   auto side_indexL=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("LeftSide"));
   auto side_indexR=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("RightSide"));
   int numberVert=vertices.size();
   assert(numberVert>1);
   typename MeshType::VertexIterator vi=vcg::tri::Allocator<MeshType>::AddVertices(mesh,numberVert);
   vector<typename MeshType::VertexPointer> pnewvertices;
   for(int i=0;i<numberVert;i++){
         int toL=side_indexL[vertices[i].first];
         int toR=side_indexR[vertices[i].first];
         int fromL=side_indexL[vertices[i].second];
         int fromR=side_indexR[vertices[i].second];
         set<int> to,from,intersect;
         to.insert(toL);to.insert(toR);
         from.insert(fromL);from.insert(fromR);
         std::set_intersection(to.begin(),to.end(),from.begin(),from.end(),
                           std::inserter(intersect,intersect.begin()));
         if(!intersect.empty()){
             int side=*(intersect.begin());
             if(intersect.size()==2 && numberVert==2){
                 side=1;
             }
             side_indexL[vi]=side;
             side_indexR[vi]=side;
         }
         else{
             side_indexL[vi]=-1;
             side_indexR[vi]=-1;
         }
         pnewvertices.push_back(&*vi);
         if(i<numberVert-1)
            vi++;
   }

   //pnewvertices.push_back(&*vi);

   vcg::tri::Allocator<MeshType>::CompactVertexVector(mesh);
   //add new faces
   int numberFaces=0;
   if(isloop)
        numberFaces=2*numberVert;
   else
        numberFaces=2*(numberVert-1);

   typename MeshType::FaceIterator fi=vcg::tri::Allocator<MeshType>::AddFaces(mesh,numberFaces);
   typename MeshType::FacePointer fstart=&*fi;
   for(int i=0;i<numberVert-1;i++){
       (*fi).Alloc(4);
       (*fi).V(0)=&mesh.vert[vertices[i].first];(*fi).V(1)=pnewvertices[i];(*fi).V(2)=pnewvertices[i+1];(*fi).V(3)=&mesh.vert[vertices[i+1].first];
       /*cout<<"face information "<<endl;
       for(int j=0;j<fi->VN();j++)
           cout<<" "<<vcg::tri::Index(mesh,fi->V(j));
       cout<<endl;*/
       fi++;
       (*fi).Alloc(4);
       (*fi).V(0)=pnewvertices[i];(*fi).V(1)=&mesh.vert[vertices[i].second];(*fi).V(2)=&mesh.vert[vertices[i+1].second];(*fi).V(3)=pnewvertices[i+1];
       /*cout<<"face information "<<endl;
       for(int j=0;j<fi->VN();j++)
           cout<<" "<<vcg::tri::Index(mesh,fi->V(j));
       cout<<endl;*/
       if(i<numberVert-2)
         fi++;
   }
   if(isloop){
       fi++;
       (*fi).Alloc(4);
       (*fi).V(0)=&mesh.vert[vertices.back().first];(*fi).V(1)=pnewvertices.back();(*fi).V(2)=pnewvertices[0];(*fi).V(3)=&mesh.vert[vertices[0].first];
       fi++;
       (*fi).Alloc(4);
       (*fi).V(0)=pnewvertices.back();(*fi).V(1)=&mesh.vert[vertices.back().second];(*fi).V(2)=&mesh.vert[vertices[0].second];(*fi).V(3)=pnewvertices[0];
   }
   //cout<<"added "<<numberFaces<<" faces"<<endl;
   UpdateTopology<MeshType>::FaceFace(mesh);
   UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
   UpdateFlags<MeshType>::VertexBorderFromFaceAdj(mesh);

   // reset h_start to a correct halfedge
   for(int k=0;k<fstart->VN();k++){
       if(fstart->V(k)==&mesh.vert[vertices[0].first] && fstart->V(fi->Next(k))==pnewvertices[0]){
            h_start.Set(fstart,k,fstart->V(k));
            //cout<<"tudo certo"<<endl;
            break;
       }
       if(fstart->V(k)==pnewvertices[0] && fstart->V(fi->Next(k))==&mesh.vert[vertices[0].first]){
            h_start.Set(fstart,k,fstart->V(k));
            h_start.FlipV();
            //cout<<"tudo certo"<<endl;
            break;
       }
   }

   // restore corners
   UpdateFlags<MeshType>::VertexClearS(mesh);
   for(int i=0;i<corners.size();i++){
       mesh.vert[corners[i]].SetS();
   }

}

template <class MeshType>
vcg::face::Pos<typename MeshType::FaceType> getPosFromIndex(MeshType& mesh,int id1,int id2){
    UpdateTopology<MeshType>::FaceFace(mesh);
    vcg::face::Pos<typename MeshType::FaceType> result;
    result.SetNull();
    for(typename MeshType::FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
        if(!fi->IsD()){
            for(int k=0;k<fi->VN();k++){
                if(vcg::tri::Index(mesh,fi->V(k))==id1 ){
                    result.Set(&*fi,k,fi->V(k));
                    break;
                }
            }
            if(!result.IsNull())
                break;
        }
    }
    while(vcg::tri::Index(mesh,result.VFlip())!=id2){        
        result.NextE();
    }
    return result;
}
template <class MeshType>
bool is_on_sideK(MeshType& mesh, typename MeshType::VertexPointer &vp,int k){
    auto side_indexL=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("LeftSide"));
    auto side_indexR=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("RightSide"));
    int L=side_indexL[vp];
    int R=side_indexR[vp];
    if(L==R){
        if(L==k)
            return true;
    }
    return false;
}
template <class MeshType>
vector<int> get_SideK(MeshType& mesh,int side){

    vcg::tri::RequireCompactness(mesh);
    vector<int> result;
    auto side_indexL=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("LeftSide"));
    auto side_indexR=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("RightSide"));
    bool founded=false;
    vcg::face::Pos<typename MeshType::FaceType> start; //pos to a extreme of the side
    set<int> startend; //store begin and end index of the side
    for(typename MeshType::FaceIterator fi=mesh.face.begin();fi!=mesh.face.end();fi++){
        for(int k=0;k<fi->VN();k++){
            if(fi->V(k)->IsS()){
                if(side_indexL[fi->V(k)]==side || side_indexR[fi->V(k)]==side){
                    //if(side_indexL[fi->V(fi->Next(k))]==side){
                        start.Set(&*fi,k,fi->V(k));
                        //start.FlipV();
                        founded=true;
                        break;
                    //}
                }
            }
        }
        if(founded)
            break;
    }
    for(typename MeshType::VertexIterator vi=mesh.vert.begin();vi!=mesh.vert.end();vi++){
        if(vi->IsS()){
            if(side_indexL[vi]==side || side_indexR[vi]==side){
                startend.insert(vcg::tri::Index(mesh,&*vi));
            }
        }
    }
    assert(startend.size()==2);
    int startindex=vcg::tri::Index(mesh,start.V());
    set<int>::iterator it=startend.begin();
    int startset=*it; it++;
    int endset=*it;
    assert(startindex==startset || startindex==endset);
    int bpos=startindex;
    int epos=(startindex==startset)?endset:startset;


    while(side_indexL[start.VFlip()]!=side && vcg::tri::Index(mesh,start.VFlip())!=epos){
        start.NextE();
    }

    start.FlipV();

    result.push_back(bpos);
    do{
        result.push_back(vcg::tri::Index(mesh,start.V()));
        vcg::face::Pos<typename MeshType::FaceType> copy=start;
        do{
            start.NextE();
        }while(side_indexL[start.VFlip()]!=side && side_indexR[start.VFlip()]!=side);
        if(copy==start)
            break;
        start.FlipV();
        if(side_indexL[start.V()]!=side_indexR[start.V()]){
            result.push_back(vcg::tri::Index(mesh,start.V()));           
            break;
        }
    }while(side_indexL[start.V()]==side_indexR[start.V()]);

    /*cout<<"side "<<side<<endl;
    for(int i=0;i<result.size();i++)
        cout<<" "<<result[i];
    cout<<endl;*/

    //store corners
    //vector<int> corners;
    vector<int> corners;
    for(typename MeshType::VertexIterator vi=mesh.vert.begin();vi!=mesh.vert.end();vi++){
        if(vi->IsS()){
            corners.push_back(vcg::tri::Index(mesh,&*vi));
        }
    }
    //cout<<"corners after side K "<<corners.size()<<endl;


    return result;
}
template <class MeshType>
void addQuadStrip_on_sideK(MeshType &mesh,int side){
    UpdateTopology<MeshType>::FaceFace(mesh);
    vector<int> sideVertices=get_SideK(mesh,side);
    /*cout<<"intial side\n";
    for(int i=0;i<sideVertices.size();i++){
        cout<<sideVertices[i]<<" ";
    }
    cout<<endl;*/
    int num_vertices=sideVertices.size();
    assert(num_vertices>1);
    assert(mesh.vert[sideVertices[0]].IsS());    
    assert(mesh.vert[sideVertices.back()].IsS());
    auto side_indexL=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("LeftSide"));
    auto side_indexR=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("RightSide"));
    typename MeshType::VertexIterator vi=vcg::tri::Allocator<MeshType>::AddVertices(mesh,num_vertices);
    typename MeshType::FaceIterator fi=vcg::tri::Allocator<MeshType>::AddFaces(mesh,num_vertices-1);
    vector<typename MeshType::VertexPointer> newvertices;
    for(int i=0;i<num_vertices;i++){
        newvertices.push_back(&*vi);
        side_indexL[vi]=side_indexL[sideVertices[i]];
        side_indexR[vi]=side_indexR[sideVertices[i]];
        side_indexL[sideVertices[i]]=-1;
        side_indexR[sideVertices[i]]=-1;
        if(i<num_vertices-1)
            vi++;
    }

    //side indexes extremes
    int el0=side_indexL[newvertices[0]];
    int er0=side_indexR[newvertices[0]];
    int eln=side_indexL[newvertices.back()];
    int ern=side_indexR[newvertices.back()];
    side_indexL[sideVertices[0]]=(el0!=side)?el0:er0;
    side_indexR[sideVertices[0]]=side_indexL[sideVertices[0]];
    side_indexL[sideVertices.back()]=(eln!=side)?eln:ern;
    side_indexR[sideVertices.back()]=side_indexL[sideVertices.back()];
    /*cout<<"side indexes after"<<endl;
    for(int i=0;i<sideVertices.size();i++){
        cout<<side_indexL[sideVertices[i]]<<" "<<side_indexR[sideVertices[i]]<<endl;
    }
    cout<<"new side indexes "<<endl;
    for(int i=0;i<sideVertices.size();i++){
        cout<<side_indexL[newvertices[i]]<<" "<<side_indexR[newvertices[i]]<<endl;
    }*/
    // update corners
    mesh.vert[sideVertices[0]].ClearS();
    mesh.vert[sideVertices.back()].ClearS();
    newvertices[0]->SetS();
    newvertices.back()->SetS();

    vcg::Point3<typename MeshType::ScalarType> normalBefore=(mesh.face[0].V(1)->P()-mesh.face[0].V(0)->P())^(mesh.face[0].V(3)->P()-mesh.face[0].V(0)->P()); // normal to some face
    for(int i=0;i<num_vertices-1;i++){
        (*fi).Alloc(4);        
        (*fi).V(0)=newvertices[i];(*fi).V(1)=newvertices[i+1];(*fi).V(2)=&mesh.vert[sideVertices[i+1]];(*fi).V(3)=&mesh.vert[sideVertices[i]];
        //(*fi).V(0)=&mesh.vert[sideVertices[i]];(*fi).V(1)=&mesh.vert[sideVertices[i+1]];(*fi).V(2)=newvertices[i+1];(*fi).V(3)=newvertices[i];
        if(i<num_vertices-2)
            fi++;
    }
    /*//correct orientation
    vcg::Point3<typename MeshType::ScalarType> normalAfter=(mesh.face.back().V(1)->P()-mesh.face.back().V(0)->P())^(mesh.face.back().V(3)->P()-mesh.face.back().V(0)->P());
    if(normalAfter.Z()*normalBefore.Z()<0){ // change orientation
        for(int j=1;j<=num_vertices-1;j++){
            typename MeshType::FaceType fii=mesh.face[mesh.FN()-j];
            typename MeshType::VertexPointer temp=fii.V(3);
            fii.V(3)=fii.V(1);
            fii.V(1)=temp;
        }
    }
    UpdateTopology<MeshType>::FaceFace(mesh);
    */
}
#endif // EDGELOOP_H
