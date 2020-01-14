#include <iostream>

#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QDirIterator>

#include <mesh_def.h>

#define WITHOUT_NUMPY
#include <matplotlibcpp.h>

using namespace std;

typedef vcg::Histogram<ScalarType> Histogram;

/// Dumps histogram h to stream with the following schema:
/// h.sum,h.min,h.max,h.avg,h.p(25),h.p(50),h.p(75),h.p(90),h.p(95);
static void dumpHistogram(Histogram & h, QTextStream & stream)
{
	stream << h.Sum() << "," << h.MinV() << "," << h.MaxV() << "," << h.Avg() << ",";
	stream << h.Percentile(.25) << "," << h.Percentile(.5) << "," << h.Percentile(.75) << "," << h.Percentile(.90) << "," << h.Percentile(.95);
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
	Histogram vertEdgeLenHist;
	vcg::PolygonalAlgorithm<PolyMesh>::InitQualityVertEdgeLenght(m);
	vcg::tri::Stat<PolyMesh>::ComputePerVertexQualityHistogram(m, vertEdgeLenHist);

	dumpHistogram(vertEdgeLenHist, geometryStream);
	dumpVertQualityHistogramPlot(m, "Edge Length", dataPath + "edgeLenHistogram.png");
	dumpVertQualityFile(m, dataPath + "vertEdgeLen.txt");

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
		std::cerr <<  "[USAGE] polyMetrics meshDir globPattern outputdir" << std::endl;
		return 1;
	}

	QFile topologyFile("./topology.csv");
	QFile geometryFile("./geometry.csv");

	if(!topologyFile.open(QFile::WriteOnly |QFile::Truncate) || !geometryFile.open(QFile::WriteOnly |QFile::Truncate))
	{
		std::cerr << "[PolyMetrics] Error opening files" << std::endl;
		return -1;
	}

	QTextStream topologyStream(&topologyFile);
	QTextStream geometryStream(&geometryFile);

	topologyStream << "ID,NUM_VERTS,NUM_EDGES,NUM_FACES,EULER,GENUS,NUM_HOLES,NUM_COMPONENTS,IS_WATERTIGHT,IS_MANIFOLD\n";
	geometryStream << "ID,TOTAL_AREA,MIN_AREA,MAX_AREA,MEAN_AREA,P25AREA,P50AREA,P75AREA,P90AREA,P95AREA,";
	geometryStream << "TOTAL_ANGLEDEV,MIN_ANGLEDEV,MAX_ANGLEDEV,MEAN_ANGLEDEV,P25ANGLEDEV,P50ANGLEDEV,P75ANGLEDEV,P90ANGLEDEV,P95ANGLEDEV,";
	geometryStream << "TOTAL_FLATNESS,MIN_FLATNESS,MAX_FLATNESS,MEAN_FLATNESS,P25FLATNESS,P50FLATNESS,P75FLATNESS,P90FLATNESS,P95FLATNESS,";
	geometryStream << "TOTAL_BENDING,MIN_BENDING,MAX_BENDING,MEAN_BENDING,P25BENDING,P50BENDING,P75BENDING,P90BENDING,P95BENDING,";
	geometryStream << "TOTAL_TORSION,MIN_TORSION,MAX_TORSION,MEAN_TORSION,P25TORSION,P50TORSION,P75TORSION,P90TORSION,P95TORSION,";
	geometryStream << "TOTAL_EDGELEN,MIN_EDGELEN,MAX_EDGELEN,MEAN_EDGELEN,P25EDGELEN,P50EDGELEN,P75EDGELEN,P90EDGELEN,P95EDGELEN,";
	geometryStream << "TOTAL_VOROAREA,MIN_VOROAREA,MAX_VOROAREA,MEAN_VOROAREA,P25VOROAREA,P50VOROAREA,P75VOROAREA,P90VOROAREA,P95VOROAREA,";


	QDir dir;
	QString basePath = dir.currentPath();
	QString targetDir = argv[1];
	QString globPattern = argv[2];
	QString outDir = argv[3];

	QDirIterator it(targetDir, QStringList() << globPattern, QDir::Files, QDirIterator::Subdirectories);
	while (it.hasNext())
	{
		it.next();


		std::string mesh = it.fileName().toStdString();

		std::cout << it.fileInfo().absolutePath().toStdString() << std::endl;


		QString outDirName = outDir + it.fileInfo().baseName();
		bool mkdir = dir.mkdir(outDirName);
		if (!mkdir)
		{
			std::cerr << "[polyMetrics] Error creating data directory for " << it.fileName().toStdString() << std::endl;
			std::cerr << "[polyMetrics] Skipping " << it.fileName().toStdString() << " ..." << std::endl;
			continue;
		}

		std::string dataDirName = QDir::toNativeSeparators(dir.absolutePath() + QDir::separator() + outDirName + QDir::separator()).toStdString();

		dir.setCurrent(it.fileInfo().absolutePath());

		std::cout << mesh << std::endl;
		std::cout << dir.currentPath().toStdString() << std::endl;
		PolyMesh m;
		int err = openMesh(m, mesh);

		if (err)
		{
//			continue;
//					dir.setCurrent(basePath);
		}

		std::cout << "computing stats..." << std::endl;

		std::string baseName = it.fileInfo().baseName().toStdString();
		computeModelStatsAndDump(m, baseName, topologyStream, geometryStream, dataDirName);

		dir.setCurrent(basePath);
	}

	topologyFile.close();
	geometryFile.close();

	return 0;
}
