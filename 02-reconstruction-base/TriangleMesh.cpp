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
	colors.resize(newTriangles.size(), glm::vec4(0.5f, 0.5f, 0.5f, 1.0f));
}

void TriangleMesh::sendToOpenGL(ShaderProgram &program)
{
	vector<glm::vec3> *glVertices = new vector<glm::vec3>;
	vector<glm::vec3> *glNormals = new vector<glm::vec3>;
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
	// VBO for vertex normals
	glBindBuffer(GL_ARRAY_BUFFER, vboColor);
	glBufferData(GL_ARRAY_BUFFER, 4 * colors.size() * sizeof(float), &colors[0].x, GL_STATIC_DRAW);
	colorLocation = program.bindVertexAttribute("color", 4, 0, 0);

	delete glVertices;
	delete glNormals;
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





