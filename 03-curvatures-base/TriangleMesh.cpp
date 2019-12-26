#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include "TriangleMesh.h"


TriangleMesh::TriangleMesh()
{
	bGLObjsInit = false;
}


void TriangleMesh::init(const vector<glm::vec3> &newVertices, const vector<unsigned int> &newTriangles)
{
	copy(newVertices.begin(), newVertices.end(), back_inserter(vertices));
	copy(newTriangles.begin(), newTriangles.end(), back_inserter(triangles));
	colors.resize(newVertices.size(), glm::vec4(0.5f, 0.5f, 0.5f, 1.0f));
	computePerVertexNormals();
}

bool TriangleMesh::load(const string &filename)
{
	ifstream fin;
	char line[100];
	int nPoints, nFaces, nVertsInFace, face[4];
	glm::vec3 P;
	
	fin.open(filename.c_str());
	if(!fin.is_open())
		return false;
	
	// Check file format
	fin.getline(line, 100);
	if(strncmp(line, "ply", strlen("ply")))
		return false;
	fin.getline(line, 100);
	if(strncmp(line, "format ascii 1.0", strlen("format ascii 1.0")))
		return false;
	
	// Read # vertices, # faces & skip the rest of the header
	fin.getline(line, 100);
	while(strncmp(line, "end_header", strlen("end_header")) != 0)
	{
		if(strncmp(line, "element vertex", strlen("element vertex")) == 0)
			nPoints = atoi(&line[strlen("element vertex")]);
		else if(strncmp(line, "element face", strlen("element face")) == 0)
			nFaces = atoi(&line[strlen("element face")]);
		fin.getline(line, 100);
	}
	
	// Read vertices
	vertices.clear();
	for(unsigned int i=0; i<nPoints; i++)
	{
		fin >> P.x >> P.y >> P.z;
		vertices.push_back(P);
	}
	
	// Read faces
	triangles.clear();
	for(unsigned int i=0; i<nFaces; i++)
	{
		fin >> nVertsInFace;
		for(unsigned int j=0; j<nVertsInFace; j++)
			fin >> face[j];
		for(unsigned int j=1; j<nVertsInFace-1; j++)
		{
			triangles.push_back(face[0]);
			triangles.push_back(face[j]);
			triangles.push_back(face[j+1]);
		}
	}
	colors.clear();
	colors.resize(vertices.size(), glm::vec4(0.5f, 0.5f, 0.5f, 1.0f));
	computePerVertexNormals();
	
	return true;
}

void TriangleMesh::getBBox(glm::vec3 box[2]) const
{
	box[0] = vertices[0];
	box[1] = vertices[0];
	for(unsigned int i=1; i<vertices.size(); i++)
	{
		box[0] = glm::min(box[0], vertices[i]);
		box[1] = glm::max(box[1], vertices[i]);
	}
}

void TriangleMesh::transform(const glm::mat4 &matrix)
{
	for(unsigned int i=0; i<vertices.size(); i++)
		vertices[i] = glm::vec3(matrix * glm::vec4(vertices[i], 1.f));
}

void TriangleMesh::sendToOpenGL(ShaderProgram &program)
{
	vector<glm::vec3> *glVertices = new vector<glm::vec3>;
	vector<glm::vec3> *glNormals = new vector<glm::vec3>;
	vector<glm::vec4> *glColors = new vector<glm::vec4>;
	glm::vec3 normal;
	
	for(unsigned int i=0; i<triangles.size(); i+=3)
	{
		glVertices->push_back(vertices[triangles[i]]);
		glVertices->push_back(vertices[triangles[i+1]]);
		glVertices->push_back(vertices[triangles[i+2]]);
		
		normal = glm::cross(vertices[triangles[i+1]] - vertices[triangles[i]], vertices[triangles[i+2]] - vertices[triangles[i]]);
		normal = glm::normalize(normal);
		
		glNormals->push_back(normal);
		glNormals->push_back(normal);
		glNormals->push_back(normal);

		glColors->push_back(colors[triangles[i]]);
		glColors->push_back(colors[triangles[i+1]]);
		glColors->push_back(colors[triangles[i+2]]);
	}
	
	// Send data to OpenGL
	if(!bGLObjsInit)
		glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	if(!bGLObjsInit)
	{
		glGenBuffers(1, &vboPosition);
		glGenBuffers(1, &vboNormal);
		glGenBuffers(1, &vboColor);
	}
	bGLObjsInit = true;
	// VBO for vertex positions
	glBindBuffer(GL_ARRAY_BUFFER, vboPosition);
	glBufferData(GL_ARRAY_BUFFER, 3 * glVertices->size() * sizeof(float), &(*glVertices)[0].x, GL_STATIC_DRAW);
	posLocation = program.bindVertexAttribute("position", 3, 0, 0);
	// VBO for vertex normals
	glBindBuffer(GL_ARRAY_BUFFER, vboNormal);
	glBufferData(GL_ARRAY_BUFFER, 3 * glNormals->size() * sizeof(float), &(*glNormals)[0].x, GL_STATIC_DRAW);
	normalLocation = program.bindVertexAttribute("normal", 3, 0, 0);
	// VBO for vertex colors
	glBindBuffer(GL_ARRAY_BUFFER, vboColor);
	glBufferData(GL_ARRAY_BUFFER, 4 * glColors->size() * sizeof(float), &(*glColors)[0].r, GL_STATIC_DRAW);
	colorLocation = program.bindVertexAttribute("color", 4, 0, 0);

	delete glVertices;
	delete glNormals;
	delete glColors;
}

void TriangleMesh::render() const
{
	glBindVertexArray(vao);
	glEnableVertexAttribArray(posLocation);
	glEnableVertexAttribArray(normalLocation);
	glEnableVertexAttribArray(colorLocation);
	glDrawArrays(GL_TRIANGLES, 0, triangles.size());
}

void TriangleMesh::free()
{
	glDeleteBuffers(1, &vboPosition);
	glDeleteBuffers(1, &vboNormal);
	glDeleteBuffers(1, &vboColor);
	glDeleteVertexArrays(1, &vao);
	
	vertices.clear();
	triangles.clear();
	colors.clear();
	
	bGLObjsInit = false;
}

void TriangleMesh::computePerVertexNormals()
{
	glm::vec3 u, v, normal;

	perVertexNormals.resize(vertices.size(), glm::vec3(0.f));
	for(unsigned int i=0; i<triangles.size(); i+=3)
	{
		u = vertices[triangles[i+1]] - vertices[triangles[i]];
		v = vertices[triangles[i+2]] - vertices[triangles[i]];
		normal = glm::normalize(glm::cross(u, v));
		perVertexNormals[triangles[i]] += normal * acos(glm::dot(glm::normalize(u), glm::normalize(v)));
		u = -u;
		v = u + v;
		perVertexNormals[triangles[i+1]] += normal * acos(glm::dot(glm::normalize(u), glm::normalize(v)));
		v = -v;
		u = v + u;
		perVertexNormals[triangles[i+2]] += normal * acos(glm::dot(glm::normalize(u), glm::normalize(v)));
	}
	for(unsigned int i=0; i<vertices.size(); i+=3)
	{
		perVertexNormals[i] = glm::normalize(perVertexNormals[i]);
	}
}







