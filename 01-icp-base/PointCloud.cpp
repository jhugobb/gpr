#include <iostream>
#include <fstream>
#include <cstring>
#include <glm/gtc/matrix_inverse.hpp>
#include "PointCloud.h"
#include "NormalEstimator.h"


PointCloud::PointCloud()
{
	bGLObjsInit = false;
}


void PointCloud::addPoint(const glm::vec3 &P, const glm::vec3 &N)
{
	points.push_back(P);
	normals.push_back(N);
}

bool PointCloud::load(const string &filename)
{
	ifstream fin;
	char line[100];
	int nPoints;
	glm::vec3 P, N = glm::vec3(0.f, 0.f, 1.f);
	bool bNormals = false;
	
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
	
	// Read # points & skip the rest of the header
	fin.getline(line, 100);
	while(strncmp(line, "end_header", strlen("end_header")) != 0)
	{
		if(strncmp(line, "element vertex", strlen("element vertex")) == 0)
			nPoints = atoi(&line[strlen("element vertex")]);
		else if(strncmp(line, "property float nx", strlen("property float nx")) == 0)
			bNormals = true;
		fin.getline(line, 100);
	}
	
	// Read point cloud
	points.clear();
	for(unsigned int i=0; i<nPoints; i++)
	{
		fin >> P.x >> P.y >> P.z;
		if(bNormals)
			fin >> N.x >> N.y >> N.z;
		addPoint(P, N);
	}
	if(!bNormals)
	{
		NormalEstimator nEstimator;
		
		nEstimator.computePointCloudNormals(points, normals);
	}
	colors.resize(points.size(), glm::vec4(0.5f, 0.5f, 0.5f, 1.0f));
	
	return true;
}

bool PointCloud::save(const string &filename) const
{
	ofstream fout;
	
	fout.open(filename.c_str());
	if(!fout.is_open())
		return false;
	
	// Write header
	fout << "ply" << endl;
	fout << "format ascii 1.0" << endl;
	fout << "element vertex " << points.size() << endl;
	fout << "property float x" << endl;
	fout << "property float y" << endl;
	fout << "property float z" << endl;
	fout << "property float nx" << endl;
	fout << "property float ny" << endl;
	fout << "property float nz" << endl;
	fout << "end_header" << endl;
	
	// Write data
	for(unsigned int pointId=0; pointId<points.size(); pointId++)
	{
		fout << points[pointId].x << " ";
		fout << points[pointId].y << " ";
		fout << points[pointId].z << endl;

		fout << normals[pointId].x << " ";
		fout << normals[pointId].y << " ";
		fout << normals[pointId].z << endl;
	}
	fout.close();

	return true;
}

unsigned int PointCloud::size() const
{
	return points.size();
}

void PointCloud::getBBox(glm::vec3 box[2]) const
{
	box[0] = points[0];
	box[1] = points[0];
	for(unsigned int i=0; i<points.size(); i++)
	{
		box[0] = glm::min(box[0], points[i]);
		box[1] = glm::max(box[1], points[i]);
	}
}

void PointCloud::sendToOpenGL(ShaderProgram &program)
{
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
	glBufferData(GL_ARRAY_BUFFER, 3 * points.size() * sizeof(float), &points[0].x, GL_STATIC_DRAW);
	posLocation = program.bindVertexAttribute("position", 3, 0, 0);
	// VBO for vertex normals
	glBindBuffer(GL_ARRAY_BUFFER, vboNormal);
	glBufferData(GL_ARRAY_BUFFER, 3 * normals.size() * sizeof(float), &normals[0].x, GL_STATIC_DRAW);
	normalLocation = program.bindVertexAttribute("normal", 3, 0, 0);
	// VBO for vertex normals
	glBindBuffer(GL_ARRAY_BUFFER, vboColor);
	glBufferData(GL_ARRAY_BUFFER, 4 * colors.size() * sizeof(float), &colors[0].x, GL_STATIC_DRAW);
	colorLocation = program.bindVertexAttribute("color", 4, 0, 0);
}

void PointCloud::render() const
{
	glBindVertexArray(vao);
	glEnableVertexAttribArray(posLocation);
	glEnableVertexAttribArray(normalLocation);
	glEnableVertexAttribArray(colorLocation);
	glDrawArrays(GL_POINTS, 0, points.size());
}

void PointCloud::free()
{
	glDeleteBuffers(1, &vboPosition);
	glDeleteBuffers(1, &vboNormal);
	glDeleteBuffers(1, &vboColor);
	glDeleteVertexArrays(1, &vao);
	
	points.clear();
}

void PointCloud::transform(const glm::mat4 &matrix)
{
	for(unsigned int i=0; i<points.size(); i++)
	{
		points[i] = glm::vec3(matrix * glm::vec4(points[i], 1.0f));
		normals[i] = glm::inverseTranspose(glm::mat3(matrix)) * normals[i];
	}
}

void PointCloud::changeAllColors(const glm::vec4 &color, ShaderProgram &program)
{
	for(unsigned int i=0; i<colors.size(); i++)
		colors[i] = color;
	sendToOpenGL(program);
}

void PointCloud::changeColorsOfIndexes(const glm::vec4 &color, ShaderProgram &program, std::vector<uint> indexes)
{
	for(unsigned int i=0; i<indexes.size(); i++)
		colors[indexes[i]] = color;
	sendToOpenGL(program);
}







