#ifndef _ITERATIVE_CLOSEST_POINT_INCLUDE
#define _ITERATIVE_CLOSEST_POINT_INCLUDE


#include "PointCloud.h"
#include "NearestNeighbors.h"
#include "CommonOperations.h"

class IterativeClosestPoint
{

public:
	void setClouds(PointCloud *pointCloud1, PointCloud *pointCloud2);
	
	void markBorderPoints();
	vector<int> *computeCorrespondence();
	glm::mat4 computeICPStep();
	
	vector<int> *computeFullICP(unsigned int maxSteps = 100);
	std::vector<glm::vec4> cloud1_colors;
	std::vector<uint> getBorderIndexes();
private:
	PointCloud *cloud1, *cloud2;
	NearestNeighbors knn;
	CommonOperations co;
	std::vector<uint> border_indexes;

};


#endif // _ITERATIVE_CLOSEST_POINT_INCLUDE


