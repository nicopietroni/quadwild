#include <iostream>

#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QDirIterator>

#include <mesh_def.h>

#include <json.hpp>

#define WITHOUT_NUMPY
#include <matplotlibcpp.h>

#define GENERATE_CSV

using namespace std;

typedef vcg::Histogram<ScalarType> Histogram;

/// Dumps histogram h to stream with the following schema:
/// h.sum,h.min,h.max,h.avg,h.p(25),h.p(50),h.p(75),h.p(90),h.p(95);
static void dumpHistogram(Histogram & h, QTextStream & stream)
{
	stream << h.Sum() << "," << h.MinV() << "," << h.MaxV() << "," << h.Avg() << ",";
	stream << h.Percentile(.25) << "," << h.Percentile(.5) << "," << h.Percentile(.75) << "," << h.Percentile(.90) << "," << h.Percentile(.95);
}

static nlohmann::json dumpHistogram(Histogram & h)
{
	nlohmann::json j = {
	    { "total", h.Sum() },
	    { "min",   h.MinV() },
	    { "max",   h.MaxV() },
	    { "mean",  h.Avg() },
	    { "P25",   h.Percentile(.25) },
	    { "P50",   h.Percentile(.50) },
	    { "P75",   h.Percentile(.75) },
	    { "P90",   h.Percentile(.90) },
	    { "P95",   h.Percentile(.95) }
	};

	return j;
}

static void writeDataToFile(const std::vector<ScalarType> & data, const std::string & filename)
{
	QFile txtFile(filename.c_str());

	if(!txtFile.open(QFile::WriteOnly |QFile::Truncate))
	{
		std::cerr << "[PolyMetrics] Error opening text file " << filename << std::endl;
		std::cerr << "[polyMetrics] Continuing..." << std::endl;
		return;
	}

	QTextStream txtStream(&txtFile);

	for (size_t i = 0; i < data.size(); ++i)
		txtStream << data[i] << "\n";
	txtStream.flush();

	txtFile.close();
}

static void dumpVertQualityFile(PolyMesh & m, const std::string filename)
{
	std::vector<ScalarType> quality(size_t(m.VN()));

	for (size_t i = 0; i < size_t(m.VN()); ++i)
		quality[i] = m.vert[i].cQ();

	writeDataToFile(quality, filename);
}

static void dumpFaceQualityFile(PolyMesh & m, const std::string filename)
{
	std::vector<ScalarType> quality(size_t(m.FN()));

	for (size_t i = 0; i < size_t(m.FN()); ++i)
		quality[i] = m.face[i].cQ();

	writeDataToFile(quality, filename);
}

static void createHistPlot (const std::vector<ScalarType> & data, const std::string & title, const std::string & filename)
{
	std::cout << filename << std::endl;
	matplotlibcpp::figure_size(1920,1080);
	std::map<std::string, std::string> gridPref;
	gridPref["linestyle"] = "--";
//	gridPref["alpha"]     = "1";
	matplotlibcpp::grid(true, "both", "y", gridPref);
	matplotlibcpp::hist(data, 50);
	matplotlibcpp::title(title);
	matplotlibcpp::save(filename);
	matplotlibcpp::close();
}

static void dumpVertQualityHistogramPlot(PolyMesh & m, const std::string title, const std::string filename)
{
	std::vector<ScalarType> quality(size_t(m.VN()));

	for (size_t i = 0; i < size_t(m.VN()); ++i)
		quality[i] = m.vert[i].cQ();

	createHistPlot(quality, title, filename);
}
static void dumpFaceQualityHistogramPlot(PolyMesh & m, const std::string title, const std::string  filename)
{
	std::vector<ScalarType> quality(size_t(m.FN()));

	for (size_t i = 0; i < size_t(m.FN()); ++i)
		quality[i] = m.face[i].cQ();

	createHistPlot(quality, title, filename);

}



static void computeModelStatsAndDump(PolyMesh & m, std::string & m_id, nlohmann::json & json)
{
	/* valence stats */
	vcg::tri::UpdateTopology<PolyMesh>::FaceFace(m);
	vcg::Histogram<ScalarType> valenceHist;
	vcg::tri::UpdateQuality<PolyMesh>::VertexValence(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, valenceHist);

	int non_manifold_edges    = vcg::tri::Clean<PolyMesh>::CountNonManifoldEdgeFF(m);
	int non_manifold_vertices = vcg::tri::Clean<PolyMesh>::CountNonManifoldVertexFF(m);

	int num_edges = 0; int num_nonManifold_edges = 0; int num_boundary_edges = 0;
	vcg::tri::Clean<PolyMesh>::CountEdgeNum(m, num_edges, num_boundary_edges, num_nonManifold_edges);
	int num_holes = -1;//vcg::tri::Clean<PolyMesh>::CountHoles(m);

	bool manifold = (non_manifold_edges + non_manifold_vertices) == 0;

	int num_cc = vcg::tri::Clean<PolyMesh>::CountConnectedComponents(m);

	int genus = vcg::tri::Clean<PolyMesh>::MeshGenus(m.VN(), num_edges, m.FN(), num_holes, num_cc);
	int euler = m.VN() - num_edges + m.FN();

	bool watertight = vcg::tri::Clean<PolyMesh>::IsWaterTight(m);

	//	topologyStream << "ID,NUM_VERTS,NUM_EDGES,NUM_FACES,EULER,GENUS,NUM_HOLES,NUM_COMPONENTS,IS_WATERTIGHT,IS_MANIFOLD\n";
	json["id"] = m_id;
	json["num_verts"] = m.VN();
	json["num_edges"] = num_edges;
	json["num_faces"] = m.FN();
	json["euler_num"] = euler;
	json["genus"] = genus;
	json["num_holes"] = num_holes;
	json["num_cc"] = num_cc;
	json["is_watertight"] = watertight;
	json["is_manifold"] = manifold;

	std::cout<<"computing geo stats" << std::endl;
	/* area stats */
	Histogram areaHist;
	for (size_t i = 0; i < m.face.size(); ++i)
		m.face[i].Q() = vcg::PolyArea(m.face[i]);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, areaHist);

	json["area"] = dumpHistogram(areaHist);
	dumpFaceQualityHistogramPlot(m, "Quad Area", "areaHistogram.png");
	dumpFaceQualityFile(m, "faceArea.txt");


	/* */
	Histogram faceAngleDeviationHist;
	vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QAngle);
	for (size_t i = 0; i < size_t(m.FN()); ++i)
		m.face[i].Q() *= 90;

	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceAngleDeviationHist);

	json["angledeviation"] = dumpHistogram(faceAngleDeviationHist);
	dumpFaceQualityHistogramPlot(m, "Angle Deviation", "angleDevHistogram.png");
	dumpFaceQualityFile(m, "faceAngleDeviation.txt");

	/* */
	Histogram faceFlatnessHist;
	vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QPlanar);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceFlatnessHist);

	json["flatness"] = dumpHistogram(faceFlatnessHist);
	dumpFaceQualityHistogramPlot(m, "Flatness", "flatnessHistogram.png");
	dumpFaceQualityFile(m, "faceFlatness.txt");

	/* */
//	Histogram faceAspectHist;
//	vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QTemplate);
//	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceAspectHist);

//	json["template"] = dumpHistogram(faceAspectHist);
//	dumpFaceQualityFile(m, "faceAspectRatio.txt");
//	dumpFaceQualityHistogramPlot(m, "Template", "aspectRatioHistogram.png");

	/* */
	Histogram faceBendingHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityFaceBending(m);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceBendingHist);

	json["bending"] = dumpHistogram(faceBendingHist);
	dumpFaceQualityHistogramPlot(m, "Bending", "bendingHistogram.png");
	dumpFaceQualityFile(m, "faceBending.txt");

	/* */
	Histogram faceTorsionHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityFaceTorsion(m);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceTorsionHist);

	json["torsion"] = dumpHistogram(faceTorsionHist);
	dumpFaceQualityHistogramPlot(m, "Torsion", "torsionHistogram.png");
	dumpFaceQualityFile(m, "faceTorsion.txt");

	/* */
	Histogram vertEdgeLenHist, vertEdgeLenDevHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertEdgeLenght(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertEdgeLenHist);

	double avg = vertEdgeLenHist.Avg();
	for (size_t i = 0; i < size_t(m.VN()); ++i)
		m.vert[i].Q() = (std::fabs(m.vert[i].Q() - avg) / avg) * 100.;

	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertEdgeLenDevHist);

	json["edgelength"] = dumpHistogram(vertEdgeLenDevHist);
	dumpVertQualityHistogramPlot(m, "Mean Edge Length Deviaiton", "edgeLenHistogram.png");
	dumpVertQualityFile(m, "vertEdgeLen.txt");

	/* */
	Histogram vertVoronoiAreaHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertVoronoiArea(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertVoronoiAreaHist);

	json["voroarea"] = dumpHistogram(vertVoronoiAreaHist);
	dumpVertQualityHistogramPlot(m, "Voronoi Area", "voroAreaHistogram.png");
	dumpVertQualityFile(m, "vertVoroArea.txt");

	std::ofstream out("stats.json");
	out << json;
	out.close();
}


static void computeModelStatsAndDump(PolyMesh & m, std::string & m_id, QTextStream &topologyStream, QTextStream &geometryStream, std::string & dataPath)
{
	/* valence stats */
	vcg::tri::UpdateTopology<PolyMesh>::FaceFace(m);
	vcg::Histogram<ScalarType> valenceHist;
	vcg::tri::UpdateQuality<PolyMesh>::VertexValence(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, valenceHist);

	int non_manifold_edges    = vcg::tri::Clean<PolyMesh>::CountNonManifoldEdgeFF(m);
	int non_manifold_vertices = vcg::tri::Clean<PolyMesh>::CountNonManifoldVertexFF(m);

	int num_edges = 0; int num_nonManifold_edges = 0; int num_boundary_edges = 0;
	vcg::tri::Clean<PolyMesh>::CountEdgeNum(m, num_edges, num_boundary_edges, num_nonManifold_edges);

	int num_holes = vcg::tri::Clean<PolyMesh>::CountHoles(m);

	bool manifold = (non_manifold_edges + non_manifold_vertices) == 0;

	int num_cc = vcg::tri::Clean<PolyMesh>::CountConnectedComponents(m);

	int genus = vcg::tri::Clean<PolyMesh>::MeshGenus(m.VN(), num_edges, m.FN(), num_holes, num_cc);
	int euler = m.VN() - num_edges + m.FN();

	bool watertight = vcg::tri::Clean<PolyMesh>::IsWaterTight(m);

	//	topologyStream << "ID,NUM_VERTS,NUM_EDGES,NUM_FACES,EULER,GENUS,NUM_HOLES,NUM_COMPONENTS,IS_WATERTIGHT,IS_MANIFOLD\n";
	topologyStream << m_id.c_str() << "," << m.VN() << "," << num_edges << "," << m.FN() << "," << euler << "," << genus << ",";
	topologyStream << num_holes << "," << num_cc << "," << watertight << "," << manifold << "\n";
	topologyStream.flush();

	std::cout<<"computing geo stats" << std::endl;
	geometryStream << m_id.c_str() << ",";
	/* area stats */
	Histogram areaHist;
	for (size_t i = 0; i < m.face.size(); ++i)
		m.face[i].Q() = vcg::PolyArea(m.face[i]);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, areaHist);

	std::cout << dataPath << std::endl;

	dumpHistogram(areaHist, geometryStream);
	dumpFaceQualityHistogramPlot(m, "Quad Area", dataPath + "areaHistogram.png");
	dumpFaceQualityFile(m, dataPath + "faceArea.txt");


	/* */
	Histogram faceAngleDeviationHist;
	vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QAngle);
	for (size_t i = 0; i < size_t(m.FN()); ++i)
		m.face[i].Q() *= 90;

	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceAngleDeviationHist);


	dumpHistogram(faceAngleDeviationHist, geometryStream);
	dumpFaceQualityHistogramPlot(m, "Angle Deviation", dataPath + "angleDevHistogram.png");
	dumpFaceQualityFile(m, dataPath + "faceAngleDeviation.txt");

	/* */
	Histogram faceFlatnessHist;
	vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QPlanar);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceFlatnessHist);

	dumpHistogram(faceFlatnessHist, geometryStream);
	dumpFaceQualityHistogramPlot(m, "Flatness", dataPath + "flatnessHistogram.png");
	dumpFaceQualityFile(m, dataPath + "faceFlatness.txt");

	/* */
//	Histogram faceAspectHist;
//	vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QTemplate);
//	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceAspectHist);

//	dumpHistogram(faceAspectHist, geometryStream);
//	dumpFaceQualityHistogramPlot(m, "aspectRatioHistogram.png");

	/* */
	Histogram faceBendingHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityFaceBending(m);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceBendingHist);

	dumpHistogram(faceBendingHist, geometryStream);
	dumpFaceQualityHistogramPlot(m, "Bending", dataPath + "bendingHistogram.png");
	dumpFaceQualityFile(m, dataPath + "faceBending.txt");

	/* */
	Histogram faceTorsionHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityFaceTorsion(m);
	vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceTorsionHist);

	dumpHistogram(faceTorsionHist, geometryStream);
	dumpFaceQualityHistogramPlot(m, "Torsion", dataPath + "torsionHistogram.png");
	dumpFaceQualityFile(m, dataPath + "faceTorsion.txt");

	/* */
	Histogram vertEdgeLenHist, vertEdgeLenDevHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertEdgeLenght(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertEdgeLenHist);

	double avg = vertEdgeLenHist.Avg();
	for (size_t i = 0; i < size_t(m.VN()); ++i)
		m.vert[i].Q() = (std::fabs(m.vert[i].Q() - avg) / avg) * 100.;

	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertEdgeLenDevHist);

	dumpVertQualityHistogramPlot(m, "Mean Edge Length Deviaiton", "edgeLenHistogram.png");
	dumpVertQualityFile(m, "vertEdgeLen.txt");
	/* */
	Histogram vertVoronoiAreaHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertVoronoiArea(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertVoronoiAreaHist);

	dumpHistogram(vertVoronoiAreaHist, geometryStream);
	dumpVertQualityHistogramPlot(m, "Voronoi Area", dataPath + "voroAreaHistogram.png");
	dumpVertQualityFile(m, dataPath + "vertVoroArea.txt");
	geometryStream << "\n";
	geometryStream.flush();
}

#ifdef GENERATE_CSV
static void computeModelStatsCSV(PolyMesh & m, std::string & m_id, QTextStream &topologyStream, QTextStream &geometryStream, std::string & dataPath)
{
    /* valence stats */
    vcg::tri::UpdateTopology<PolyMesh>::FaceFace(m);
    vcg::Histogram<ScalarType> valenceHist;
    vcg::tri::UpdateQuality<PolyMesh>::VertexValence(m);
    vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, valenceHist);

    int non_manifold_edges    = vcg::tri::Clean<PolyMesh>::CountNonManifoldEdgeFF(m);
    int non_manifold_vertices = vcg::tri::Clean<PolyMesh>::CountNonManifoldVertexFF(m);

    int num_edges = 0; int num_nonManifold_edges = 0; int num_boundary_edges = 0;
    vcg::tri::Clean<PolyMesh>::CountEdgeNum(m, num_edges, num_boundary_edges, num_nonManifold_edges);

    int num_holes = vcg::tri::Clean<PolyMesh>::CountHoles(m);

    bool manifold = (non_manifold_edges + non_manifold_vertices) == 0;

    int num_cc = vcg::tri::Clean<PolyMesh>::CountConnectedComponents(m);

    int genus = vcg::tri::Clean<PolyMesh>::MeshGenus(m.VN(), num_edges, m.FN(), num_holes, num_cc);
    int euler = m.VN() - num_edges + m.FN();

    bool watertight = vcg::tri::Clean<PolyMesh>::IsWaterTight(m);

    int numSingularities = 0;
    vcg::tri::UpdateQuality<PolyMesh>::VertexValence(m);
    for (size_t i = 0; i < m.vert.size(); i++) {
        if (m.vert[i].Q() != 4) {
            numSingularities++;
        }
    }

    //	topologyStream << "ID,NUM_VERTS,NUM_EDGES,NUM_FACES,EULER,GENUS,NUM_HOLES,NUM_COMPONENTS,IS_WATERTIGHT,IS_MANIFOLD,NUM_SINGULARITIES\n";
    topologyStream << m_id.c_str() << "," << m.VN() << "," << num_edges << "," << m.FN() << "," << euler << "," << genus << ",";
    topologyStream << num_holes << "," << num_cc << "," << watertight << "," << manifold << "," << numSingularities << "\n";
    topologyStream.flush();

    //   geometryStream << "ID,"
    std::cout<<"computing geo stats" << std::endl;
    geometryStream << m_id.c_str() << ",";


    /* area stats */
    //   geometryStream << "TOTAL_AREA,MIN_AREA,MAX_AREA,MEAN_AREA,P25AREA,P50AREA,P75AREA,P90AREA,P95AREA,";
    Histogram areaHist;
    for (size_t i = 0; i < m.face.size(); ++i)
        m.face[i].Q() = vcg::PolyArea(m.face[i]);
    vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, areaHist);
    std::cout << dataPath << std::endl;
    dumpHistogram(areaHist, geometryStream);
    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_ANGLEDEV,MIN_ANGLEDEV,MAX_ANGLEDEV,MEAN_ANGLEDEV,P25ANGLEDEV,P50ANGLEDEV,P75ANGLEDEV,P90ANGLEDEV,P95ANGLEDEV,";
    Histogram faceAngleDeviationHist;
    vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QAngle);

	for (size_t i = 0; i < size_t(m.FN()); ++i)
		m.face[i].Q() *= 90;


    vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceAngleDeviationHist);
    dumpHistogram(faceAngleDeviationHist, geometryStream);
    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_FLATNESS,MIN_FLATNESS,MAX_FLATNESS,MEAN_FLATNESS,P25FLATNESS,P50FLATNESS,P75FLATNESS,P90FLATNESS,P95FLATNESS,";
    Histogram faceFlatnessHist;
    vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QPlanar);
    vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceFlatnessHist);
    dumpHistogram(faceFlatnessHist, geometryStream);
    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_TEMPLATE,MIN_TEMPLATE,MAX_TEMPLATE,MEAN_TEMPLATE,P25TEMPLATE,P50TEMPLATE,P75TEMPLATE,P90TEMPLATE,P95TEMPLATE,";
//    Histogram faceAspectHist;
//    vcg::PolygonalAlgorithm<PolyMesh>::UpdateQuality(m, vcg::PolygonalAlgorithm<PolyMesh>::QTemplate);
//    vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceAspectHist);
//    dumpHistogram(faceAspectHist, geometryStream);
//    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_BENDING,MIN_BENDING,MAX_BENDING,MEAN_BENDING,P25BENDING,P50BENDING,P75BENDING,P90BENDING,P95BENDING,";
    Histogram faceBendingHist;
    vcg::PolygonalAlgorithm<PolyMesh>::InitQualityFaceBending(m);
    vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceBendingHist);
    dumpHistogram(faceBendingHist, geometryStream);
    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_TORSION,MIN_TORSION,MAX_TORSION,MEAN_TORSION,P25TORSION,P50TORSION,P75TORSION,P90TORSION,P95TORSION,";
    Histogram faceTorsionHist;
    vcg::PolygonalAlgorithm<PolyMesh>::InitQualityFaceTorsion(m);
    vcg::tri::Stat<PolyMesh>::ComputePerFaceQualityHistogram(m, faceTorsionHist);
    dumpHistogram(faceTorsionHist, geometryStream);
    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_EDGELEN,MIN_EDGELEN,MAX_EDGELEN,MEAN_EDGELEN,P25EDGELEN,P50EDGELEN,P75EDGELEN,P90EDGELEN,P95EDGELEN,";
	Histogram vertEdgeLenHist, vertEdgeLenDevHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertEdgeLenght(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertEdgeLenHist);

	double avg = vertEdgeLenHist.Avg();
	for (size_t i = 0; i < size_t(m.VN()); ++i)
		m.vert[i].Q() = (std::fabs(m.vert[i].Q() - avg) / avg) * 100.;

	dumpHistogram(vertEdgeLenDevHist, geometryStream);
    geometryStream << ",";

    /* */
    //   geometryStream << "TOTAL_VOROAREA,MIN_VOROAREA,MAX_VOROAREA,MEAN_VOROAREA,P25VOROAREA,P50VOROAREA,P75VOROAREA,P90VOROAREA,P95VOROAREA,";
    Histogram vertVoronoiAreaHist;
    vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertVoronoiArea(m);
    vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertVoronoiAreaHist);
    dumpHistogram(vertVoronoiAreaHist, geometryStream);

    geometryStream << "\n";
    geometryStream.flush();
}
#endif


static int openMesh(PolyMesh & m, std::string & name)
{
	int mask;
	int err = vcg::tri::io::ImporterOBJ<PolyMesh>::Open(m, name.c_str(), mask);
	if (err)
	{
		std::cerr << "[polyMetrics] Import Error: " << vcg::tri::io::ImporterOBJ<PolyMesh>::ErrorMsg(err) <<  std::endl;
		return err;
	}

	return 0;
}

int main(int argc, char * argv[])
{

	if (argc < 3)
	{
        std::cout <<  "[USAGE] polyMetrics baseDirectory globpattern" << std::endl;
        std::cout << "BaseDirectory: the directory that contains all the target subdirectories" << std::endl;
        std::cout << "globpattern: the pattern that matches the filename of the quad mesh to compute the metrics" << std::endl;

		return 1;
	}

	matplotlibcpp::backend("WXAgg");

#ifdef GENERATE_CSV
    QFile topologyFile("./topology.csv");
    QFile geometryFile("./geometry.csv");

    if(!topologyFile.open(QFile::WriteOnly |QFile::Truncate) || !geometryFile.open(QFile::WriteOnly |QFile::Truncate))
    {
        std::cerr << "[PolyMetrics] Error opening files" << std::endl;
        return -1;
    }

    QTextStream topologyStream(&topologyFile);
    QTextStream geometryStream(&geometryFile);

    topologyStream << "ID,NUM_VERTS,NUM_EDGES,NUM_FACES,EULER,GENUS,NUM_HOLES,NUM_COMPONENTS,IS_WATERTIGHT,IS_MANIFOLD,NUM_SINGULARITIES\n";

    geometryStream << "ID,TOTAL_AREA,MIN_AREA,MAX_AREA,MEAN_AREA,P25AREA,P50AREA,P75AREA,P90AREA,P95AREA,";
    geometryStream << "TOTAL_ANGLEDEV,MIN_ANGLEDEV,MAX_ANGLEDEV,MEAN_ANGLEDEV,P25ANGLEDEV,P50ANGLEDEV,P75ANGLEDEV,P90ANGLEDEV,P95ANGLEDEV,";
    geometryStream << "TOTAL_FLATNESS,MIN_FLATNESS,MAX_FLATNESS,MEAN_FLATNESS,P25FLATNESS,P50FLATNESS,P75FLATNESS,P90FLATNESS,P95FLATNESS,";
//    geometryStream << "TOTAL_TEMPLATE,MIN_TEMPLATE,MAX_TEMPLATE,MEAN_TEMPLATE,P25TEMPLATE,P50TEMPLATE,P75TEMPLATE,P90TEMPLATE,P95TEMPLATE,";
    geometryStream << "TOTAL_BENDING,MIN_BENDING,MAX_BENDING,MEAN_BENDING,P25BENDING,P50BENDING,P75BENDING,P90BENDING,P95BENDING,";
    geometryStream << "TOTAL_TORSION,MIN_TORSION,MAX_TORSION,MEAN_TORSION,P25TORSION,P50TORSION,P75TORSION,P90TORSION,P95TORSION,";
    geometryStream << "TOTAL_EDGELEN,MIN_EDGELEN,MAX_EDGELEN,MEAN_EDGELEN,P25EDGELEN,P50EDGELEN,P75EDGELEN,P90EDGELEN,P95EDGELEN,";
    geometryStream << "TOTAL_VOROAREA,MIN_VOROAREA,MAX_VOROAREA,MEAN_VOROAREA,P25VOROAREA,P50VOROAREA,P75VOROAREA,P90VOROAREA,P95VOROAREA\n";
#endif

	QDir dir;
	QString basePath = dir.currentPath();
	QString targetDir = argv[1];
	QString globPattern = argv[2];

    std::cout << "Target Dir: " << targetDir.toStdString() << " glob pattern: " << globPattern.toStdString() << std::endl;

//	QDirIterator it(targetDir, QStringList() << globPattern, QDir::Dirs);
    QDirIterator it(targetDir, QStringList() << globPattern, QDir::Files, QDirIterator::Subdirectories);
	while (it.hasNext())
	{
        it.next();
        dir.setCurrent(it.fileInfo().dir().path());
        std::cout << dir.currentPath().toStdString() << std::endl;

//		if (dir.exists("edgeLenHistogram.png") && dir.exists("flatnessHistogram.png") && dir.exists("voroAreaHistogram.png") && dir.exists("torsionHistogram.png"))
//		{
//			dir.setCurrent(basePath);
//			continue;
//		}

		nlohmann::json json;
//		std::string dataDirName = QDir::toNativeSeparators(dir.absolutePath() + QDir::separator() + outDirName + QDir::separator()).toStdString();

		PolyMesh m;
        std::string mesh = it.fileName().toStdString();// + "_rem_p0_0_quadrangulation_smooth.obj";
		int err = openMesh(m, mesh);

		std::cout << mesh << " " << m.VN() << std::endl;

		if (err)
			continue;

		std::cout << "computing stats..." << std::endl;

		std::string baseName = it.fileInfo().baseName().toStdString();
		computeModelStatsAndDump(m, baseName, json);
		
#ifdef GENERATE_CSV
        std::string dataDirName = targetDir.toStdString();
        computeModelStatsCSV(m, baseName, topologyStream, geometryStream, dataDirName);
#endif

		dir.setCurrent(basePath);
	}


#ifdef GENERATE_CSV
	topologyFile.close();
	geometryFile.close();
#endif


	return 0;
}
