#ifndef _SCENE_INCLUDE
#define _SCENE_INCLUDE


#include <glm/glm.hpp>
#include "Camera.h"
#include "ShaderProgram.h"
#include "PointCloud.h"
#include "SegmentCloud.h"
#include "IterativeClosestPoint.h"


// Scene contains all the entities of our game.
// It is responsible for updating and render them.


class Scene
{

public:
	Scene();
	~Scene();

	void init();
	void update(int deltaTime);
	void render();
	
	bool loadScans(const char *filename1, const char *filename2);
	
	void render_gui();
	
	void getBBox(glm::vec3 bbox[2]) const;
	void addCloud(PointCloud *cloud);
	void transform(const glm::mat4 &matrix);

	Camera &getCamera();
	
	void switchPolygonMode();

private:
	void initShaders();
	void computeModelViewMatrix();
	void applyICPStep();
	void applyFullICP();

private:
	Camera camera;
	vector<PointCloud *> clouds;
	SegmentCloud segments;
	ShaderProgram pointsProgram, segmentsProgram;
	float currentTime;
	
	IterativeClosestPoint icp;
	
	bool bLighting;
	bool bDraw1, bDraw2, bDrawLines;
	int cloudToDraw, pointSize;
	glm::mat4 icpTransform;

};


#endif // _SCENE_INCLUDE

