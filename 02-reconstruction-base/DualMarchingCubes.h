#ifndef _DUAL_MARCHING_CUBES_INCLUDE
#define _DUAL_MARCHING_CUBES_INCLUDE


#include <map>
#include <set>
#include "PointCloud.h"
#include "ImplicitFunction.h"
#include "TriangleMesh.h"


struct NodePosition
{
	int coords[3];
	
	bool operator<(const NodePosition &pos) const;
};


class DualMarchingCubes
{

public:
	TriangleMesh *extractZeroIsosurface(const ImplicitFunction &function, unsigned int resolution = 64);
	TriangleMesh *extractZeroIsosurfaceFast(const PointCloud *cloud, float radius, const ImplicitFunction &function, unsigned int resolution = 64);
	
private:
	void computeSurfaceNodes(const PointCloud *cloud, float radius, set<NodePosition> &nodes);
	void extractVertices(const ImplicitFunction &function, vector<glm::vec3> &vertices, map<NodePosition, unsigned int> &nodes2Vertices);
	void extractVerticesFast(set<NodePosition> &nodes, const ImplicitFunction &function, vector<glm::vec3> &vertices, map<NodePosition, unsigned int> &nodes2Vertices);
	void extractFaces(const ImplicitFunction &function, map<NodePosition, unsigned int> &nodes2Vertices, vector<unsigned int> &triangles);

	bool isMixedValidNode(const ImplicitFunction &function, const NodePosition &pos, float values[2][2][2]);
	glm::vec3 extractVertex(const NodePosition &pos, float values[2][2][2]);
	
	void printProgressBar(float progress);

private:	
	unsigned int numNodes;
	float nodeSize;
	
};


#endif // _DUAL_MARCHING_CUBES_INCLUDE




