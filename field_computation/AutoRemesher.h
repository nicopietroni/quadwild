#ifndef AUTOREMESHER_H
#define AUTOREMESHER_H

#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>

#include <vcg/complex/append.h>

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>

#include <memory>

#define DEBUG
template <class Mesh>
class AutoRemesher {

	typedef typename Mesh::ScalarType ScalarType;
	typedef typename Mesh::CoordType  CoordType;

	typedef typename Mesh::VertexType    VertexType;
	typedef typename Mesh::VertexPointer VertexPointer;

	typedef typename Mesh::FaceType    FaceType;
	typedef typename Mesh::FacePointer FacePointer;

	typedef vcg::GridStaticPtr<FaceType, ScalarType> StaticGrid;

	static ScalarType computeAR(Mesh & m, const double perc = 0.05)
	{
		vcg::tri::ForEachFace(m, [] (FaceType & f) {
			f.Q() = vcg::Quality(f.cP(0), f.cP(1), f.cP(2));
		});

		vcg::Histogram<ScalarType> hist;
		vcg::tri::Stat<Mesh>::ComputePerFaceQualityHistogram(m, hist);

		//		    return hist.MinV();
		return hist.Percentile(perc);
	}

	static void removeColinearFaces(Mesh & m)
	{
		vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

		int count = 0;

		Mesh projectMesh;
		vcg::tri::Append<Mesh, Mesh>::MeshCopy(projectMesh, m);
		vcg::tri::UpdateBounding<Mesh>::Box(projectMesh);
		vcg::tri::UpdateNormal<Mesh>::PerVertexNormalizedPerFace(projectMesh);
		StaticGrid grid;
		grid.Set(projectMesh.face.begin(), projectMesh.face.end());


		do
		{
			vcg::tri::UpdateTopology<Mesh>::FaceFace(m);
			vcg::tri::UnMarkAll(m);

			count = 0;
			for (size_t i = 0; i < size_t(m.FN()); ++i)
			{
				FaceType & f = m.face[i];

				ScalarType quality = vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2));
				ScalarType area = vcg::DoubleArea(f);

				if ((quality <= 0.00000001 /*&& area <= 0.00000001*/) /*|| minEdge <= 0.000001*/)
				{
					//find longest edge
					double edges[3];
					edges[0] = vcg::Distance(f.cP(0), f.cP(1));
					edges[1] = vcg::Distance(f.cP(1), f.cP(2));
					edges[2] = vcg::Distance(f.cP(2), f.cP(0));
					int longestIdx = std::find(edges, edges+3, std::max(std::max(edges[0], edges[1]), edges[2])) - (edges);

					if (vcg::tri::IsMarked(m, f.V2(longestIdx)))
						continue;


					auto f1 = f.cFFp(longestIdx);
					vcg::tri::Mark(m,f.V2(longestIdx));
					if (!vcg::face::IsBorder(f, longestIdx) && vcg::face::IsManifold(f, longestIdx) && vcg::face::checkFlipEdgeNotManifold<FaceType>(f, longestIdx))  {

						// Check if EdgeFlipping improves quality
						FacePointer g = f.FFp(longestIdx); int k = f.FFi(longestIdx);
						vcg::Triangle3<ScalarType> t1(f.P(longestIdx), f.P1(longestIdx), f.P2(longestIdx)), t2(g->P(k), g->P1(k), g->P2(k)),
						        t3(f.P(longestIdx), g->P2(k), f.P2(longestIdx)), t4(g->P(k), f.P2(longestIdx), g->P2(k));

						if ( std::min( QualityFace(t1), QualityFace(t2) ) <= std::min( QualityFace(t3), QualityFace(t4) ))
						{
							ScalarType dist;
							CoordType closest;
							auto fp0 = vcg::tri::GetClosestFaceBase(projectMesh, grid, vcg::Barycenter(t3), 0.000001, dist, closest);
							if (fp0 == NULL)
								continue;

							auto fp1 = vcg::tri::GetClosestFaceBase(projectMesh, grid, vcg::Barycenter(t4), 0.000001, dist, closest);
							if (fp1 == NULL)
								continue;

							vcg::face::FlipEdgeNotManifold<FaceType>(f, longestIdx);
							++count;
						}
					}
				}
			}
			std::cout << count << std::endl;
		} while (count);
	}

public:

	typedef struct Params {

		int iterations = 15;

		ScalarType targetAspect = 0.35;
		int targetDeltaFN = 5000;

		int initialApproximateFN = 20000;

		ScalarType creaseAngle = 25.;
		bool userSelectedCreases = true;
		bool surfDistCheck = true;

	} Params;



	static std::shared_ptr<Mesh> CleanMesh(Mesh & m, const bool & splitNonManifold = false)
	{
		std::shared_ptr<Mesh> ret = std::make_shared<Mesh>();

		vcg::tri::Append<Mesh, Mesh>::MeshCopy(*ret, m);

		vcg::tri::Clean<Mesh>::RemoveDuplicateFace(*ret);
		vcg::tri::Clean<Mesh>::RemoveUnreferencedVertex(*ret);
		vcg::tri::Allocator<Mesh>::CompactEveryVector(*ret);

		vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);

		std::cout << "[AutoRemeshCleaner] Removing colinear faces by flip..." << std::endl;
		removeColinearFaces(*ret);


		if (splitNonManifold)
		{
			vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);
			vcg::tri::Clean<Mesh>::SplitNonManifoldVertex(*ret, 0.0000001);
			vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);
			vcg::tri::ForEachFace(*ret, [] (FaceType & f) {
				for (int i = 0; i < 3; ++i)
				{
					if (!vcg::face::IsManifold(f, i))
						f.SetFaceEdgeS(i);
				}
			});
			vcg::tri::CutMeshAlongSelectedFaceEdges(*ret);
		}

		vcg::tri::UpdateTopology<Mesh>::FaceFace(*ret);

		std::cout << "[AutoRemeshCleaner] Attempting to coherently orient faces..." << std::endl;
		bool oriented = false, orientable = false;
		vcg::tri::Clean<Mesh>::OrientCoherentlyMesh(*ret, oriented, orientable);

		std::cout << "[AutoRemeshCleaner] Orientable:" << orientable << std::endl;
		std::cout << "[AutoRemeshCleaner] Oriented:"   << oriented << std::endl;

		return ret;
	}

	//for big meshes disabling par.surfDistCheck provides big perf improvements, sacrificing result accuracy
	static std::shared_ptr<Mesh> Remesh (Mesh & m, Params & par)
	{
		std::shared_ptr<Mesh> ret = std::make_shared<Mesh>();

		vcg::tri::UpdateBounding<Mesh>::Box(m);
		vcg::tri::UpdateTopology<Mesh>::FaceFace(m);

		vcg::tri::IsotropicRemeshing<Mesh>::Params para;
		para.iter = par.iterations;
		para.SetFeatureAngleDeg(par.creaseAngle);
		para.splitFlag    = true;
		para.swapFlag     = true;
		para.collapseFlag = true;
		para.smoothFlag   = true;
		para.projectFlag  = true;
		para.selectedOnly = false;
		para.adapt=false;

		para.aspectRatioThr = 0;

		para.maxSurfDist = m.bbox.Diag() / 2000.;
		para.surfDistCheck = par.surfDistCheck;
		para.userSelectedCreases = par.userSelectedCreases;

		ScalarType prevFN = m.FN();
		ScalarType deltaFN = m.FN();

		ScalarType aspect = 0;
		ScalarType edgeLow  = 0;
		ScalarType edgeL = std::sqrt(vcg::tri::Stat<Mesh>::ComputeMeshArea(m) * 2 / par.initialApproximateFN);
		ScalarType edgeHigh = edgeL * 2.;

		do {
			vcg::tri::Append<Mesh, Mesh>::MeshCopy(*ret, m);

			para.SetTargetLen(edgeL);
#ifdef DEBUG
			std::cout << edgeLow << "--" << edgeL << "--" << edgeHigh << std::endl;
#endif

			vcg::tri::IsotropicRemeshing<Mesh>::Do(*ret, para);

			aspect = computeAR(*ret);
			ScalarType newFN = ret->FN();

			deltaFN = std::abs(ret->FN() - prevFN);
			prevFN = newFN;

			if (aspect >= par.targetAspect /*&& ret->FN() > par.initialApproximateFN*/)
			{
				edgeLow = edgeL;
			} else {
				edgeHigh = edgeL;
			}

			edgeL = (edgeHigh + edgeLow) / 2.;
#ifdef DEBUG
			std::cout << "Delta FN " << deltaFN << std::endl;
			std::cout << "Mesh FN  " << ret->FN() << std::endl;
			std::cout << "Min AR   " << computeAR(*ret) << std::endl;
			vcg::tri::io::Exporter<Mesh>::Save(
			            *ret,
			            "resultrem.ply",
			            vcg::tri::io::Mask::IOM_FACEQUALITY
			            );
#endif
		}
		while (aspect < par.targetAspect || deltaFN > par.targetDeltaFN);

		std::cerr << "[REMESH] RemeshedFaces:" << ret->FN() << std::endl;
		std::cerr << "[REMESH] RemeshedAspect:" << computeAR(*ret) << std::endl;

		return ret;
	}
};

#endif // AUTOREMESHER_H
