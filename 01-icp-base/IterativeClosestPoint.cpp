#include <iostream>
#include <algorithm>
#include <set>
#include <math.h>
#include "IterativeClosestPoint.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <glm/gtc/matrix_transform.hpp>


void IterativeClosestPoint::setClouds(PointCloud *pointCloud1, PointCloud *pointCloud2)
{
	cloud1 = pointCloud1;
	cloud2 = pointCloud2;
	co = CommonOperations();
	co.setNN(&cloud1->getPoints());
}

std::vector<uint> IterativeClosestPoint::getBorderIndexes() {
	return border_indexes;
}

// This method should mark the border points in cloud 1. It also changes their color (for example to red).
// You will need to add an attribute to this class that stores this property for all points in cloud 1. 

void IterativeClosestPoint::markBorderPoints()
{
	// TODO
	std::vector<glm::vec3> points_1 = cloud1->getPoints();
	std::vector<size_t> neighbors;
	std::vector<float> squared_dists;
	std::vector<Eigen::Vector3f> neigh_points;
	std::vector<Eigen::Vector3f> adj_points;
	Eigen::Vector3f centroid;
	Eigen::Matrix3f covariance;
	for (uint i = 0; i < points_1.size(); i++) {
		neigh_points = co.getNN(points_1[i], points_1, 60);
		Eigen::Vector3f point_in_eigen(points_1[i].x, points_1[i].y, points_1[i].z);
		centroid = co.computeCentroidOfPoints(neigh_points);
		adj_points = co.adjustPointsByCentroid(centroid, neigh_points);
		covariance = co.covarianceMatrix(adj_points);

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(covariance);
		es.compute(covariance);

		Eigen::Vector3cf eigen_values = es.eigenvalues();
		Eigen::Matrix3cf eigen_vectors = es.eigenvectors();

		Eigen::Vector3f v1(eigen_vectors.col(2)(0,0).real(),
											 eigen_vectors.col(2)(1,0).real(),
											 eigen_vectors.col(2)(2,0).real());
		
		Eigen::Vector3f v2(eigen_vectors.col(1)(0,0).real(),
										   eigen_vectors.col(1)(1,0).real(),
											 eigen_vectors.col(1)(2,0).real());

		Eigen::Vector3f v3(eigen_vectors.col(0)(0,0).real(),
											 eigen_vectors.col(0)(1,0).real(),
											 eigen_vectors.col(0)(2,0).real());

		std::vector<Eigen::Vector3f> transformed_neighs;
		Eigen::Vector3f e;
		Eigen::Vector3f n(0,0,1);
		for (uint j = 0; j < neigh_points.size(); j++) {
			e = Eigen::Vector3f(v1.dot(neigh_points[j] - point_in_eigen),
													v2.dot(neigh_points[j] - point_in_eigen),
													v3.dot(neigh_points[j] - point_in_eigen));
			e = e - e.dot(n)*n; // project e into z=0 plane
			transformed_neighs.push_back(e);
		}

		std::vector<float> angles;
		for (uint j = 0; j < neigh_points.size(); j++) {
			angles.push_back(std::atan2(transformed_neighs[j].y(), transformed_neighs[j].x()));
		}

		std::sort(angles.begin(), angles.end());

		float delta;
		float max_delta = std::numeric_limits<float>::min();
		for (uint j = 0; j < angles.size()-1; j++) {
			delta = angles[j+1] - angles[j];
			if (delta > max_delta) {
				max_delta = delta;
			}
		}

		delta = 2*M_PI + angles[0] - angles[angles.size()-1];
		if (delta > max_delta) max_delta = delta;
		if (max_delta > M_PI/2) border_indexes.push_back(i);
	}

	// cloud1_colors = cloud1->getColors();
}


// This method should compute the closest point in cloud 1 for all non border points in cloud 2. 
// This correspondence will be useful to compute the ICP step matrix that will get cloud 2 closer to cloud 1.
// Store the correspondence in this class as the following method is going to need it.
// As it is evident in its signature this method also returns the correspondence. The application draws this if available.

vector<int> *IterativeClosestPoint::computeCorrespondence()
{
	// TODO
	std::set<size_t> closests;
	for (uint i = 0; i < border_indexes.size(); i++) {
		closests.insert(border_indexes[i]);
	}
	std::vector<int>* result = new std::vector<int>();
	std::vector<glm::vec3> points_2 = cloud2->getPoints();
	std::vector<size_t> neighbors;
	for (uint i = 0; i < points_2.size(); i++) {
		neighbors = co.getNN_indexes(points_2[i],1);

		if (closests.find(neighbors[0]) != closests.end()) result->push_back(-1);
		else result->push_back(neighbors[0]);
	}
	return result;
}


// This method should compute the rotation and translation of an ICP step from the correspondence
// information between clouds 1 and 2. Both should be encoded in the returned 4x4 matrix.
// To do this use the SVD algorithm in Eigen.

glm::mat4 IterativeClosestPoint::computeICPStep()
{
	// TODO
	
	return glm::mat4(1.f);
}


// This method should perform the whole ICP algorithm with as many steps as needed.
// It should stop when maxSteps are performed, when the Frobenius norm of the transformation matrix of
// a step is smaller than a small threshold, or when the correspondence does not change from the 
// previous step.

vector<int> *IterativeClosestPoint::computeFullICP(unsigned int maxSteps)
{
	// TODO
	
	return NULL;
}





