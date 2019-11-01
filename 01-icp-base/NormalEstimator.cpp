#include "NormalEstimator.h"
#include <iostream>




// This method has to compute a normal per point in the 'points' vector and put it in the 
// 'normals' vector. The 'normals' vector already has the same size as the 'points' vector. 
// There is no need to push_back new elements, only to overwrite them ('normals[i] = ...;')
// In order to compute the normals you should use PCA. The provided 'NearestNeighbors' class
// wraps the nanoflann library that computes K-nearest neighbors effciently. 

void NormalEstimator::computePointCloudNormals(const vector<glm::vec3> &points, vector<glm::vec3> &normals)
{
	// // TODO
	// nn = NearestNeighbors();

	// nn.setPoints(&points);
	co = CommonOperations();
	co.setNN(&(points));
	std::vector<Eigen::Vector3f> neighbors;

	std::vector<Eigen::Vector3f> adj_points;
	Eigen::Vector3f centroid;

	for (uint i = 0; i < points.size(); i++) {
		
		neighbors = co.getNN(points[i], points, 10);
		
		centroid = co.computeCentroidOfPoints(neighbors);

		adj_points = co.adjustPointsByCentroid(centroid, neighbors);

		Eigen::Matrix3f covariance;

		covariance = co.covarianceMatrix(adj_points);

		// TODO Calcuate spectral decompos
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(covariance);
		es.compute(covariance);

		Eigen::Vector3cf eigen_values = es.eigenvalues();
		Eigen::Matrix3cf eigen_vectors = es.eigenvectors();

		// std::cout << es.eigenvectors().row(minindex) << std::endl;
		normals[i] = glm::vec3(eigen_vectors.col(0)(0,0).real(), eigen_vectors.col(0)(1,0).real(), eigen_vectors.col(0)(2,0).real());
		if(normals[i].z < 0.0f) normals[i].z = - normals[i].z;

	}





	
}


