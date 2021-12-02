#include <iostream>

#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QDirIterator>

#include <mesh_def.h>

using namespace std;

static void writeDataToFile(const std::vector<std::string> & data, const std::string & filename)
{
    QFile txtFile(filename.c_str());

    if(!txtFile.open(QFile::WriteOnly |QFile::Truncate))
    {
        std::cerr << "[orientability_check] Error opening text file " << filename << std::endl;
        std::cerr << "[orientability_check] Continuing..." << std::endl;
        return;
    }

    QTextStream txtStream(&txtFile);

    for (size_t i = 0; i < data.size(); ++i)
        txtStream << data[i].c_str() << "\n";
    txtStream.flush();

    txtFile.close();
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
        std::cout <<  "[USAGE] polyMetrics baseDirectory globpattern" << std::endl;
        std::cout << "BaseDirectory: the directory that contains all the target subdirectories" << std::endl;
        std::cout << "globpattern: the pattern that matches the filename of the quad mesh to compute the metrics" << std::endl;

		return 1;
	}

	QDir dir;
	QString basePath = dir.currentPath();
	QString targetDir = argv[1];
	QString globPattern = argv[2];

    std::cout << "Target Dir: " << targetDir.toStdString() << " glob pattern: " << globPattern.toStdString() << std::endl;

    std::vector<std::string> notOrientables;

    QDirIterator it(targetDir, QStringList() << globPattern, QDir::Files, QDirIterator::Subdirectories);
	while (it.hasNext())
	{
        it.next();
        dir.setCurrent(it.fileInfo().dir().path());
        std::cout << dir.currentPath().toStdString() << std::endl;

		PolyMesh m;
        std::string mesh = it.fileName().toStdString();// + "_rem_p0_0_quadrangulation_smooth.obj";
		int err = openMesh(m, mesh);

		std::cout << mesh << " " << m.VN() << std::endl;

		if (err)
			continue;

        std::cout << "computing orientability..." << std::endl;

        bool isOriented = false;
        bool isOrientable = false;

        vcg::tri::UpdateTopology<PolyMesh>::FaceFace(m);
        vcg::tri::Clean<PolyMesh>::OrientCoherentlyMesh(m, isOriented, isOrientable);

        std::string baseName = dir.current().dirName().toStdString();

        if (!isOrientable)
        {
           notOrientables.push_back(baseName);
        }

		dir.setCurrent(basePath);
	}

    writeDataToFile(notOrientables, "notOrientable.txt");

	return 0;
}
