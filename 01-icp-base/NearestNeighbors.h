#ifndef _NEAREST_NEIGHBORS_INCLUDE
#define _NEAREST_NEIGHBORS_INCLUDE


#include <cstdlib>
#include <vector>
#include <glm/glm.hpp>
#include "nanoflann/nanoflann.hpp"


using namespace std;


struct NN_Points
{

public:
	size_t kdtree_get_point_count() const;
	float kdtree_get_pt(const size_t idx, const size_t dim) const;
	template <class BBOX> bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

public:
	const vector<glm::vec3> *points;

};


typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, NN_Points>, NN_Points, 3> kd_tree_t;


class NearestNeighbors
{

public:
	NearestNeighbors();
	~NearestNeighbors();
	
	void setPoints(const vector<glm::vec3> *points);
	unsigned int getKNearestNeighbors(const glm::vec3 &P, unsigned int K, vector<size_t> &neighbors, vector<float> &dists_squared);
	unsigned int getNeighborsInRadius(const glm::vec3 &P, float radius, vector< pair<size_t, float> > &neighbors);

private:
	NN_Points nnPoints;
	kd_tree_t *kd_tree;

};


#endif // _NEAREST_NEIGHBORS_INCLUDE


