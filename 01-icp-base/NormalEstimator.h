#ifndef _NORMAL_ESTIMATOR_INCLUDE
#define _NORMAL_ESTIMATOR_INCLUDE


#include <vector>
#include <glm/glm.hpp>

#include "CommonOperations.h"


#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>


using namespace std;


class NormalEstimator
{
private:
	NearestNeighbors nn;
	CommonOperations co;
public:
	void computePointCloudNormals(const vector<glm::vec3> &points, vector<glm::vec3> &normals);

};


#endif // _NORMAL_ESTIMATOR_INCLUDE


