#include "CommonOperations.h"


CommonOperations::CommonOperations() {
  nn = NearestNeighbors();
}

CommonOperations::~CommonOperations() {}

void CommonOperations::setNN(const std::vector<glm::vec3>* points) {
  nn.setPoints(points);
}

// Computes the K-nearest neighbors of p_i in a point set points
std::vector<Eigen::Vector3f> CommonOperations::getNN(glm::vec3 p_i, std::vector<glm::vec3> points, uint k) {
	std::vector<size_t> neighbors;
	std::vector<float> dists_squared;
	std::vector<Eigen::Vector3f> result;
	nn.getKNearestNeighbors(p_i, k, neighbors, dists_squared);

	for (uint i = 0; i < neighbors.size(); i++) {
		glm::vec3 neighbor = points[neighbors[i]];
		Eigen::Vector3f n(neighbor.x, neighbor.y, neighbor.z);
		result.push_back(n);
	}

	return result;
}

std::vector<size_t> CommonOperations::getNN_indexes(glm::vec3 p_i, uint k) {
	std::vector<size_t> neighbors;
	std::vector<float> dists_squared;
	std::vector<Eigen::Vector3f> result;
	nn.getKNearestNeighbors(p_i, k, neighbors, dists_squared);

	return neighbors;
}

// Computes the centroid of a point set points
Eigen::Vector3f CommonOperations::computeCentroidOfPoints(std::vector<Eigen::Vector3f> points) {
	Eigen::Vector3f centroid(0,0,0);
	for (uint i = 0; i < points.size(); i++) {
		centroid += points[i];
	}
	centroid = centroid / (float) points.size();
	return centroid;
}

// Returns a vector of points displaced by a centroid
std::vector<Eigen::Vector3f> CommonOperations::adjustPointsByCentroid(Eigen::Vector3f centroid, std::vector<Eigen::Vector3f> points) {
	std::vector<Eigen::Vector3f> result;

	for (uint i = 0; i < points.size(); i++) {
		Eigen::Vector3f adj_point = points[i] - centroid;
		result.push_back(adj_point);
	}

	return result;
}

Eigen::Matrix3f CommonOperations::covarianceMatrix(std::vector<Eigen::Vector3f> points) {
  Eigen::Matrix3f covariance;
  covariance << 0.0f, 0.0f, 0.0f,
          			0.0f, 0.0f, 0.0f,
          			0.0f, 0.0f, 0.0f;
  for (uint i = 0; i < points.size(); i++) {
    covariance += points[i] * points[i].transpose();
  }

  return covariance;
}