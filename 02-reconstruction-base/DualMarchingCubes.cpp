#include <iostream>
#include "DualMarchingCubes.h"
#include <Eigen/Dense>


#define sqr(x) ((x)*(x))


bool NodePosition::operator<(const NodePosition &pos) const
{
	return (coords[0] < pos.coords[0]) || (coords[0] == pos.coords[0] && coords[1] < pos.coords[1]) ||
	       (coords[0] == pos.coords[0] && coords[1] == pos.coords[1] && coords[2] < pos.coords[2]);
}


TriangleMesh *DualMarchingCubes::extractZeroIsosurface(const ImplicitFunction &function, unsigned int resolution)
{
	vector<glm::vec3> vertices;
	map<NodePosition, unsigned int> nodes2Vertices;
	vector<unsigned int> triangles;
	TriangleMesh *mesh;
	
	nodeSize = 1.f / resolution;
	numNodes = resolution;
	
	extractVertices(function, vertices, nodes2Vertices);
	extractFaces(function, nodes2Vertices, triangles);
	
	mesh = new TriangleMesh();
	mesh->init(vertices, triangles);
	
	return mesh;
}

TriangleMesh *DualMarchingCubes::extractZeroIsosurfaceFast(const PointCloud *cloud, float radius, const ImplicitFunction &function, unsigned int resolution)
{
	set<NodePosition> nodes;
	vector<glm::vec3> vertices;
	map<NodePosition, unsigned int> nodes2Vertices;
	vector<unsigned int> triangles;
	TriangleMesh *mesh;
	
	nodeSize = 1.f / resolution;
	numNodes = resolution;
	
	computeSurfaceNodes(cloud, radius, nodes);
	extractVerticesFast(nodes, function, vertices, nodes2Vertices);
	extractFaces(function, nodes2Vertices, triangles);
	
	mesh = new TriangleMesh();
	mesh->init(vertices, triangles);
	
	return mesh;
}

void DualMarchingCubes::computeSurfaceNodes(const PointCloud *cloud, float radius, set<NodePosition> &nodes)
{
	glm::vec3 box[2];
	NodePosition posBox[2], node;
	
	for(unsigned int i=0; i<cloud->getPoints().size(); i++)
	{
		box[0] = cloud->getPoints()[i] - glm::vec3(radius, radius, radius);
		box[1] = cloud->getPoints()[i] + glm::vec3(radius, radius, radius);
		/*
		box[0] = cloud->getPoints()[i];
		box[0].x -= radius * (1 - sqr(cloud->getNormals()[i].x));
		box[0].y -= radius * (1 - sqr(cloud->getNormals()[i].y));
		box[0].z -= radius * (1 - sqr(cloud->getNormals()[i].z));
		box[1] = cloud->getPoints()[i];
		box[1].x += radius * (1 - sqr(cloud->getNormals()[i].x));
		box[1].y += radius * (1 - sqr(cloud->getNormals()[i].y));
		box[1].z += radius * (1 - sqr(cloud->getNormals()[i].z));
		*/
		
		posBox[0].coords[0] = int((box[0].x + 0.5f) / nodeSize);
		posBox[0].coords[1] = int((box[0].y + 0.5f) / nodeSize);
		posBox[0].coords[2] = int((box[0].z + 0.5f) / nodeSize);
		posBox[1].coords[0] = int((box[1].x + 0.5f) / nodeSize);
		posBox[1].coords[1] = int((box[1].y + 0.5f) / nodeSize);
		posBox[1].coords[2] = int((box[1].z + 0.5f) / nodeSize);
		for(node.coords[0]=posBox[0].coords[0]; node.coords[0]<=posBox[1].coords[0]; node.coords[0]++)
			for(node.coords[1]=posBox[0].coords[1]; node.coords[1]<=posBox[1].coords[1]; node.coords[1]++)
				for(node.coords[2]=posBox[0].coords[2]; node.coords[2]<=posBox[1].coords[2]; node.coords[2]++)
				{
					if(node.coords[0] >= 0 && node.coords[0] < numNodes && node.coords[1] >= 0 && node.coords[1] < numNodes && 
					   node.coords[2] >= 0 && node.coords[2] < numNodes)
						nodes.insert(node);
				}
	}
	
	cout << "#Surface nodes = " << nodes.size() << endl;
	cout << "#Total nodes = " << (numNodes*numNodes*numNodes) << endl;
}

void DualMarchingCubes::extractVertices(const ImplicitFunction &function, vector<glm::vec3> &vertices, map<NodePosition, unsigned int> &nodes2Vertices)
{
	NodePosition pos;
	float values[2][2][2];
	
	cout << "Extracting vertices" << endl;
	printProgressBar(0.f);
	
	for(pos.coords[2]=0; pos.coords[2]<numNodes; pos.coords[2]++)
	{
		for(pos.coords[1]=0; pos.coords[1]<numNodes; pos.coords[1]++)
			for(pos.coords[0]=0; pos.coords[0]<numNodes; pos.coords[0]++)
			{
				if(isMixedValidNode(function, pos, values))
				{
					vertices.push_back(extractVertex(pos, values));
					nodes2Vertices[pos] = vertices.size()-1;
				}
			}
		printProgressBar(float(pos.coords[2]) / numNodes);
	}
	
	cout << endl;
}

void DualMarchingCubes::extractVerticesFast(set<NodePosition> &nodes, const ImplicitFunction &function, vector<glm::vec3> &vertices, map<NodePosition, unsigned int> &nodes2Vertices)
{
	unsigned int i = 0, size = nodes.size();
	NodePosition pos;
	float values[2][2][2];
	
	cout << "Extracting vertices" << endl;
	printProgressBar(0.f);
	
	for(set<NodePosition>::iterator it=nodes.begin(); it!=nodes.end(); it++)
	{
		if(isMixedValidNode(function, *it, values))
		{
			vertices.push_back(extractVertex(*it, values));
			nodes2Vertices[*it] = vertices.size()-1;
		}
		printProgressBar(float(i) / size);
		i++;
	}
	
	cout << endl;
}

void DualMarchingCubes::extractFaces(const ImplicitFunction &function, map<NodePosition, unsigned int> &nodes2Vertices, vector<unsigned int> &triangles)
{
	NodePosition pos;
	unsigned int face[4];
	glm::vec3 P;
	float v1, v2;
	
	cout << "Extracting faces" << endl;
	printProgressBar(0.f);

	// Extract faces in X direction
	for(pos.coords[2]=0; pos.coords[2]<numNodes-1; pos.coords[2]++)
	{
		for(pos.coords[1]=0; pos.coords[1]<numNodes-1; pos.coords[1]++)
			for(pos.coords[0]=0; pos.coords[0]<numNodes; pos.coords[0]++)
			{
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[0] = nodes2Vertices[pos];
				else
					continue;
				pos.coords[1]++;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[1] = nodes2Vertices[pos];
				else
				{
					pos.coords[1]--;
					continue;
				}
				pos.coords[2]++;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[2] = nodes2Vertices[pos];
				else
				{
					pos.coords[1]--;
					pos.coords[2]--;
					continue;
				}
				pos.coords[1]--;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[3] = nodes2Vertices[pos];
				else
				{
					pos.coords[2]--;
					continue;
				}
				pos.coords[2]--;
				
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(pos.coords[0], pos.coords[1] + 1.f, pos.coords[2] + 1.f);
				function(P, v1);
				function(P + glm::vec3(nodeSize, 0.0f, 0.0f), v2);
				if(glm::sign(v1) == glm::sign(v2))
				{
					function(P - glm::vec3(0.15f * nodeSize, 0.0f, 0.0f), v1);
					function(P + glm::vec3(nodeSize + 0.3f * nodeSize, 0.0f, 0.0f), v2);
				}

				if(v1 > 0.f)
				{
					triangles.push_back(face[2]);
					triangles.push_back(face[1]);
					triangles.push_back(face[0]);
						
					triangles.push_back(face[0]);
					triangles.push_back(face[3]);
					triangles.push_back(face[2]);
				}
				else
				{
					triangles.push_back(face[0]);
					triangles.push_back(face[1]);
					triangles.push_back(face[2]);
						
					triangles.push_back(face[2]);
					triangles.push_back(face[3]);
					triangles.push_back(face[0]);
				}
			}
		printProgressBar(0.33f * float(pos.coords[2]) / numNodes);
	}

	// Extract faces in Y direction
	for(pos.coords[2]=0; pos.coords[2]<numNodes-1; pos.coords[2]++)
	{
		for(pos.coords[1]=0; pos.coords[1]<numNodes; pos.coords[1]++)
			for(pos.coords[0]=0; pos.coords[0]<numNodes-1; pos.coords[0]++)
			{
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[0] = nodes2Vertices[pos];
				else
					continue;
				pos.coords[2]++;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[1] = nodes2Vertices[pos];
				else
				{
					pos.coords[2]--;
					continue;
				}
				pos.coords[0]++;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[2] = nodes2Vertices[pos];
				else
				{
					pos.coords[2]--;
					pos.coords[0]--;
					continue;
				}
				pos.coords[2]--;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[3] = nodes2Vertices[pos];
				else
				{
					pos.coords[0]--;
					continue;
				}
				pos.coords[0]--;
				
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(pos.coords[0] + 1.f, pos.coords[1], pos.coords[2] + 1.f);
				function(P, v1);
				function(P + glm::vec3(0.0f, nodeSize, 0.0f), v2);
				if(glm::sign(v1) == glm::sign(v2))
				{
					function(P - glm::vec3(0.0f, 0.15f * nodeSize, 0.0f), v1);
					function(P + glm::vec3(0.0f, nodeSize + 0.3f * nodeSize, 0.0f), v2);
				}

				if(v1 > 0.f)
				{
					triangles.push_back(face[2]);
					triangles.push_back(face[1]);
					triangles.push_back(face[0]);
						
					triangles.push_back(face[0]);
					triangles.push_back(face[3]);
					triangles.push_back(face[2]);
				}
				else
				{
					triangles.push_back(face[0]);
					triangles.push_back(face[1]);
					triangles.push_back(face[2]);
						
					triangles.push_back(face[2]);
					triangles.push_back(face[3]);
					triangles.push_back(face[0]);
				}
			}
		printProgressBar(0.33f + 0.33f * float(pos.coords[2]) / numNodes);
	}

	// Extract faces in Z direction
	for(pos.coords[2]=0; pos.coords[2]<numNodes; pos.coords[2]++)
	{
		for(pos.coords[1]=0; pos.coords[1]<numNodes-1; pos.coords[1]++)
			for(pos.coords[0]=0; pos.coords[0]<numNodes-1; pos.coords[0]++)
			{
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[0] = nodes2Vertices[pos];
				else
					continue;
				pos.coords[0]++;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[1] = nodes2Vertices[pos];
				else
				{
					pos.coords[0]--;
					continue;
				}
				pos.coords[1]++;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[2] = nodes2Vertices[pos];
				else
				{
					pos.coords[0]--;
					pos.coords[1]--;
					continue;
				}
				pos.coords[0]--;
				if(nodes2Vertices.find(pos) != nodes2Vertices.end())
					face[3] = nodes2Vertices[pos];
				else
				{
					pos.coords[1]--;
					continue;
				}
				pos.coords[1]--;
				
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(pos.coords[0] + 1.f, pos.coords[1] + 1.f, pos.coords[2]);
				function(P, v1);
				function(P + glm::vec3(0.0f, 0.0f, nodeSize), v2);
				if(glm::sign(v1) == glm::sign(v2))
				{
					function(P - glm::vec3(0.0f, 0.0f, 0.15f * nodeSize), v1);
					function(P + glm::vec3(0.0f, 0.0f, nodeSize + 0.3f * nodeSize), v2);
				}

				if(v1 > 0.f)
				{
					triangles.push_back(face[2]);
					triangles.push_back(face[1]);
					triangles.push_back(face[0]);
						
					triangles.push_back(face[0]);
					triangles.push_back(face[3]);
					triangles.push_back(face[2]);
				}
				else
				{
					triangles.push_back(face[0]);
					triangles.push_back(face[1]);
					triangles.push_back(face[2]);
						
					triangles.push_back(face[2]);
					triangles.push_back(face[3]);
					triangles.push_back(face[0]);
				}
			}
		printProgressBar(0.67f + 0.33f * float(pos.coords[2]) / numNodes);
	}
			
	cout << endl;
}

bool DualMarchingCubes::isMixedValidNode(const ImplicitFunction &function, const NodePosition &pos, float values[2][2][2])
{
	glm::vec3 P;
	glm::ivec3 glmPos = glm::ivec3(pos.coords[0], pos.coords[1], pos.coords[2]), delta;
	
	for(delta.z=0; delta.z<2; delta.z++)
		for(delta.y=0; delta.y<2; delta.y++)
			for(delta.x=0; delta.x<2; delta.x++)
			{
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				if(!function(P, values[delta.x][delta.y][delta.z]))
					return false;
			}

	for(delta.z=0; delta.z<2; delta.z++)
		for(delta.y=0; delta.y<2; delta.y++)
			for(delta.x=0; delta.x<2; delta.x++)
				if(glm::sign(values[delta.x][delta.y][delta.z]) != glm::sign(values[0][0][0]))
					return true;
	
	return false;
}

glm::vec3 DualMarchingCubes::extractVertex(const NodePosition &pos, float values[2][2][2])
{
	glm::vec3 vertex, P, Q;
	glm::ivec3 glmPos = glm::ivec3(pos.coords[0], pos.coords[1], pos.coords[2]), delta;
	vector<glm::vec3> points;
	float v1, v2;
	
	// Compute intersections with sticks of this node
	// Edges parallel to X
	delta.x = 0;
	for(delta.y=0; delta.y<2; delta.y++)
		for(delta.z=0; delta.z<2; delta.z++)
		{
			v1 = values[delta.x][delta.y][delta.z];
			v2 = values[delta.x+1][delta.y][delta.z];
			if(glm::sign(v1) != glm::sign(v2))
			{
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				delta.x = 1;
				Q = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				delta.x = 0;
				if(v1 < 0.f)
					v1 = -v1;
				if(v2 < 0.f)
					v2 = -v2;
				if(v1 + v2 < 1e-6)
					vertex = (P + Q) / 2.f;
				else
					vertex = (v2 * P + v1 * Q) / (v2 + v1);
				points.push_back(vertex);
			}
		}

	// Edges parallel to Y
	delta.y = 0;
	for(delta.x=0; delta.x<2; delta.x++)
		for(delta.z=0; delta.z<2; delta.z++)
		{
			v1 = values[delta.x][delta.y][delta.z];
			v2 = values[delta.x][delta.y+1][delta.z];
			if(glm::sign(v1) != glm::sign(v2))
			{
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				delta.y = 1;
				Q = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				delta.y = 0;
				if(v1 < 0.f)
					v1 = -v1;
				if(v2 < 0.f)
					v2 = -v2;
				if(v1 + v2 < 1e-6)
					vertex = (P + Q) / 2.f;
				else
					vertex = (v2 * P + v1 * Q) / (v2 + v1);
				points.push_back(vertex);
			}
		}

	// Edges parallel to Z
	delta.z = 0;
	for(delta.x=0; delta.x<2; delta.x++)
		for(delta.y=0; delta.y<2; delta.y++)
		{
			v1 = values[delta.x][delta.y][delta.z];
			v2 = values[delta.x][delta.y][delta.z+1];
			if(glm::sign(v1) != glm::sign(v2))
			{
				P = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				delta.z = 1;
				Q = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * glm::vec3(glmPos + delta);
				delta.z = 0;
				if(v1 < 0.f)
					v1 = -v1;
				if(v2 < 0.f)
					v2 = -v2;
				if(v1 + v2 < 1e-6)
					vertex = (P + Q) / 2.f;
				else
					vertex = (v2 * P + v1 * Q) / (v2 + v1);
				points.push_back(vertex);
			}
		}
	
	// Estimate plane through points
	glm::vec3 centroid;
	Eigen::Matrix3f covariance;
	Eigen::Vector3f V;
	glm::vec3 normal;
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver;
	
	centroid = glm::vec3(0.f);
	for(unsigned int j=0; j<points.size(); j++)
		centroid += points[j];
	centroid /= points.size();
	covariance << 0,0,0,0,0,0,0,0,0;
	for(unsigned int j=0; j<points.size(); j++)
	{
		V << points[j].x-centroid.x, points[j].y-centroid.y, points[j].z-centroid.z;
		covariance += V * V.transpose();
	}
	if(covariance.norm() < 1e-6)
		vertex = centroid;
	else
	{
		eigensolver.compute(covariance);
		if (eigensolver.info() == Eigen::Success)
		{
			V = eigensolver.eigenvectors().col(0);
			normal = glm::normalize(glm::vec3(V(0), V(1), V(2)));
		}
		
		// Project node center to estimated plane
		vertex = glm::vec3(-0.5f, -0.5f, -0.5f) + nodeSize * (glm::vec3(glmPos) + glm::vec3(0.5f, 0.5f, 0.5f));
		vertex = vertex - (glm::dot(normal, vertex) - glm::dot(normal, centroid)) * normal;
	}
	
	return vertex;
}

void DualMarchingCubes::printProgressBar(float progress)
{
	int numAsterisks = int(50.99f * progress);
	
	for(unsigned int i=0; i<60; i++)
		cout << "\b";
	cout << "[";
	for(unsigned int i=0; i<numAsterisks; i++)
		cout << "*";
	for(unsigned int i=0; i<(50-numAsterisks); i++)
		cout << ".";
	cout << "]";
	cout.flush();
}








