#include "SegmentCloud.h"
#include <iostream>

SegmentCloud::SegmentCloud()
{
	bGLObjsInit = false;
}


void SegmentCloud::buildFromCorrespondence(const PointCloud *cloud1, const PointCloud *cloud2, const vector<int> &correspondence)
{
	points.clear();
	
	for(unsigned int i=0; i<correspondence.size(); i++)
	{
		if (correspondence[i] == -1) std::cout << "hello" << endl;
		if(correspondence[i] != -1)
		{
			points.push_back(cloud2->getPoints()[i]);
			points.push_back(cloud1->getPoints()[correspondence[i]]);
		}
	}
}

void SegmentCloud::sendToOpenGL(ShaderProgram &program)
{
	// Send data to OpenGL
	if(!bGLObjsInit)
		glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	if(!bGLObjsInit)
		glGenBuffers(1, &vboPosition);
	bGLObjsInit = true;
	// VBO for vertex positions
	glBindBuffer(GL_ARRAY_BUFFER, vboPosition);
	glBufferData(GL_ARRAY_BUFFER, 3 * points.size() * sizeof(float), &points[0].x, GL_STATIC_DRAW);
	posLocation = program.bindVertexAttribute("position", 3, 0, 0);
}

void SegmentCloud::render() const
{
	if(bGLObjsInit)
	{
		glBindVertexArray(vao);
		glEnableVertexAttribArray(posLocation);
		glDrawArrays(GL_LINES, 0, points.size());
	}
}

void SegmentCloud::free()
{
	glDeleteBuffers(1, &vboPosition);
	glDeleteVertexArrays(1, &vao);
	
	points.clear();
}




