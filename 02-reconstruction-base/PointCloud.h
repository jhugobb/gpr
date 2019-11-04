#ifndef _POINT_CLOUD_INCLUDE
#define _POINT_CLOUD_INCLUDE


#include <vector>
#include <glm/glm.hpp>
#include "ShaderProgram.h"


using namespace std;


class PointCloud
{
public:
	PointCloud();

	void addPoint(const glm::vec3 &P, const glm::vec3 &N);

	bool load(const string &filename);
	bool save(const string &filename) const;
	
	void changeAllColors(const glm::vec4 &color, ShaderProgram &program);

	unsigned int size() const;
	void getBBox(glm::vec3 box[2]) const;
	void transform(const glm::mat4 &matrix);
	
	void sendToOpenGL(ShaderProgram &program);
	void render() const;
	void free();
	
	const vector<glm::vec3> &getPoints() const { return points; }
	vector<glm::vec3> &getPoints() { return points; }
	const vector<glm::vec3> &getNormals() const { return normals; }
	vector<glm::vec3> &getNormals() { return normals; }
	const vector<glm::vec4> &getColors() const { return colors; }
	vector<glm::vec4> &getColors() { return colors; }
	
	float computeDistanceBetweenPoints();
	
private:
	vector<glm::vec3> points, normals;
	vector<glm::vec4> colors;

	bool bGLObjsInit;
	GLuint vao;
	GLuint vboPosition, vboNormal, vboColor;
	GLint posLocation, normalLocation, colorLocation;

};


#endif // _POINT_CLOUD_INCLUDE


