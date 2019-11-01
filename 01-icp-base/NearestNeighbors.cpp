#include "NearestNeighbors.h"


NearestNeighbors::NearestNeighbors()
{
	nnPoints.points = NULL;
	kd_tree = NULL;
}

NearestNeighbors::~NearestNeighbors()
{
	if(kd_tree != NULL)
		delete kd_tree;
}


void NearestNeighbors::setPoints(const vector<glm::vec3> *points)
{
	nnPoints.points = points;
	if(kd_tree != NULL)
		delete kd_tree;
	kd_tree = new kd_tree_t(3, nnPoints, nanoflann::KDTreeSingleIndexAdaptorParams(10) );
	kd_tree->buildIndex();
}

unsigned int NearestNeighbors::getKNearestNeighbors(const glm::vec3 &P, unsigned int K, vector<size_t> &neighbors, vector<float> &dists_squared)
{
	if(nnPoints.points == NULL)
		return 0;

	neighbors.clear();
	neighbors.resize(K);
	dists_squared.clear();
	dists_squared.resize(K);
	size_t num_results = kd_tree->knnSearch(&P[0], K, &neighbors[0], &dists_squared[0]);
	neighbors.resize(num_results);
	dists_squared.resize(num_results);
	
	return num_results;
}

unsigned int NearestNeighbors::getNeighborsInRadius(const glm::vec3 &P, float radius, vector< pair<size_t, float> > &neighbors)
{
	if(nnPoints.points == NULL)
		return 0;
	
	neighbors.clear();
	size_t nMatches = kd_tree->radiusSearch(&P[0], radius, neighbors, nanoflann::SearchParams());
	
	return nMatches;
}

size_t NN_Points::kdtree_get_point_count() const
{
	if(points == NULL)
		return 0;
	
	return points->size();
}

float NN_Points::kdtree_get_pt(const size_t idx, const size_t dim) const
{
	if(points == NULL)
		return 0.f;
	
	return (*points)[idx][dim];
}


