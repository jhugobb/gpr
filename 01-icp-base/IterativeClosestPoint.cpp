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
	cout << border_indexes.size() << endl;
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
	correspondence = new std::vector<int>();
	std::vector<glm::vec3> points_2 = cloud2->getPoints();
	std::vector<size_t> neighbors;
	for (uint i = 0; i < points_2.size(); i++) {
		neighbors = co.getNN_indexes(points_2[i],1);

		if (closests.find(neighbors[0]) != closests.end()) correspondence->push_back(-1);
		else correspondence->push_back(neighbors[0]);
	}
	return correspondence;
}


// This method should compute the rotation and translation of an ICP step from the correspondence
// information between clouds 1 and 2. Both should be encoded in the returned 4x4 matrix.
// To do this use the SVD algorithm in Eigen.

glm::mat4 IterativeClosestPoint::computeICPStep()
{
	// TODO
	std::vector<glm::vec3> points_1 = cloud1->getPoints();
	std::vector<glm::vec3> points_2 = cloud2->getPoints();

	std::vector<Eigen::Vector3f> P;
	std::vector<Eigen::Vector3f> Q;

	Eigen::Vector3f p,q;
	int index;
	for (uint i  = 0; i < points_2.size(); i++) {
		index = correspondence->at(i);
		if (index == -1) continue;
		q = Eigen::Vector3f(points_2[i].x, points_2[i].y, points_2[i].z);
		p = Eigen::Vector3f(points_1[correspondence->at(i)].x, 
												points_1[correspondence->at(i)].y, 
												points_1[correspondence->at(i)].z);
		Q.push_back(q);
		P.push_back(p);
	} 
	Eigen::Vector3f Q_centroid = co.computeCentroidOfPoints(Q);
	Eigen::Vector3f P_centroid = co.computeCentroidOfPoints(P);

	Q = co.adjustPointsByCentroid(Q_centroid, Q);
	P = co.adjustPointsByCentroid(P_centroid, P);

	Eigen::Matrix3Xf Q_mat(3, Q.size()), P_mat(3, P.size());

	for (uint i  = 0; i < Q.size(); i++) {
		Q_mat.col(i) = Q[i];
		P_mat.col(i) = P[i];
	}

	Eigen::Matrix3f S_mat = Q_mat * P_mat.transpose();
	Eigen::JacobiSVD<Eigen::Matrix3Xf> svd(S_mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3f R_mat = svd.matrixV() * svd.matrixU().transpose();

	if (R_mat.determinant() == -1) {
		Eigen::Matrix3f I;
		I << 1,  0,  0,  0,
				 0,  1,  0,  0,
				 0,  0,  1,  0,
				 0,  0,  0, -1;
		R_mat = svd.matrixV() * I * svd.matrixU().transpose();
	}	

	Eigen::Vector3f translation = P_centroid - R_mat * Q_centroid;

	glm::mat4 result;
	Eigen::Matrix4f Rtrans;

	Rtrans.col(0) = Eigen::Vector4f(R_mat.col(0)[0],R_mat.col(0)[1],R_mat.col(0)[2], 1);
	Rtrans.col(1) = Eigen::Vector4f(R_mat.col(1)[0],R_mat.col(1)[1],R_mat.col(1)[2], 1);
	Rtrans.col(2) = Eigen::Vector4f(R_mat.col(2)[0],R_mat.col(2)[1],R_mat.col(2)[2], 1);
	Rtrans.col(3) = Eigen::Vector4f(translation.x(), translation.y(), translation.z(),  1);
	frob_norm = Rtrans.norm();
	for (uint i = 0; i < Rtrans.rows(); i++) {
		for (uint j = 0; j < Rtrans.cols(); j++) {
			result[i][j] = Rtrans.row(i)[j];
		}
	}
	return glm::transpose(result);
}

// Returns true if the two arrays are equal, and false otherwise
bool IterativeClosestPoint::checkCorrespondence(std::vector<int>* curr, std::vector<int>* prev) {
	for (uint i = 0; i < curr->size(); i++) {
		if (curr->at(i) != prev->at(i)) return false;
	}
	return true;
}


// This method should perform the whole ICP algorithm with as many steps as needed.
// It should stop when maxSteps are performed, when the Frobenius norm of the transformation matrix of
// a step is smaller than a small threshold, or when the correspondence does not change from the 
// previous step.

vector<int> *IterativeClosestPoint::computeFullICP(unsigned int maxSteps)
{
	// TODO
	glm::mat4 icpTransform;
	std::vector<int>* prev_corr = computeCorrespondence();
	std::vector<int>* curr_corr;

	for (uint iter = 0; iter < maxSteps; iter++) {
		icpTransform = computeICPStep();
		cloud2->transform(icpTransform);
		if (frob_norm < 0.0000001f) break;
		curr_corr = computeCorrespondence();
		if (checkCorrespondence(curr_corr, prev_corr)) break;
		prev_corr = curr_corr;
	}
	
	return curr_corr;
}





