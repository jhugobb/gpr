#ifndef _SEGMENT_CLOUD_INCLUDE
#define _SEGMENT_CLOUD_INCLUDE


#include <vector>
#include <glm/glm.hpp>
#include "ShaderProgram.h"
#include "PointCloud.h"


using namespace std;


class SegmentCloud
{

public:
	SegmentCloud();
	
	void buildFromCorrespondence(const PointCloud *cloud1, const PointCloud *cloud2, const vector<int> &correspondence);
	
	void sendToOpenGL(ShaderProgram &program);
	void render() const;
	void free();
	
private:
	vector<glm::vec3> points;

	bool bGLObjsInit;
	GLuint vao;
	GLuint vboPosition;
	GLint posLocation;

};


#endif // _SEGMENT_CLOUD_INCLUDE



