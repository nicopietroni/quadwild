#include "load_save.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "assert.h"

std::vector<std::vector<size_t>> loadPatches(const std::string& filename)
{
    std::vector<std::vector<size_t>> partitions;

    std::ifstream input;
    input.open(filename.c_str());
    if (!input.is_open())
    {
        std::cout<<"ERROR LOADING PATCH FILE"<<std::endl;
        exit(0);
    }

    size_t numFaces;
    input >> numFaces;

    int maxPartitionId = 0;
    std::vector<int> facePartition(numFaces,-1);

    for (size_t i=0; i < numFaces;i++)
    {
        input >> facePartition[i];
        maxPartitionId = std::max(facePartition[i], maxPartitionId);
    }

    input.close();

    partitions.clear();
    partitions.resize(static_cast<size_t>(maxPartitionId+1));
    for(size_t i = 0; i < facePartition.size(); i++)
    {
        int currentPartition = facePartition[i];
        assert(currentPartition >= 0);
        assert(currentPartition < static_cast<int>(partitions.size()));
        partitions[static_cast<size_t>(currentPartition)].push_back(i);
    }

    return partitions;
}

std::vector<std::vector<size_t>> loadCorners(const std::string& filename)
{
    std::vector<std::vector<size_t>> corners;

    std::ifstream input;
    input.open(filename.c_str());
    if (!input.is_open())
    {
        std::cout<<"ERROR LOADING CORNERS FILE"<<std::endl;
        exit(0);
    }

    size_t numPartitions;
    input >> numPartitions;

    corners.clear();
    corners.resize(numPartitions);
    for (size_t i = 0; i < numPartitions; i++)
    {
        size_t numCorners;
        input >> numCorners;

        corners[i].resize(numCorners);
        for (size_t j = 0; j < numCorners; j++)
            input >> corners[i][j];
    }

    input.close();

    return corners;
}


std::vector<std::pair<size_t,size_t> > LoadFeatures(const std::string &filename)
{
    std::vector<std::pair<size_t,size_t> > features;
    FILE *f=NULL;
    f=fopen(filename.c_str(),"rt");
    if (f==NULL)
    {
        std::cout<<"ERROR LOADING FEATURES FILE"<<std::endl;
        exit(0);
    }

    int numFeatures;
    fscanf(f,"%d\n",&numFeatures);

    for (size_t i=0;i<numFeatures;i++)
    {
        int FIndex,EIndex;
        fscanf(f,"%d,%d\n",&FIndex,&EIndex);

        features.push_back(std::pair<size_t,size_t>(FIndex,EIndex));
    }
    fclose(f);

    return features;
}
