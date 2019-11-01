#include <iostream>
#include <fstream>
#include <cmath>
#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include "imgui/imgui.h"
#include "Scene.h"


Scene::Scene()
{
}

Scene::~Scene()
{
	for(vector<PointCloud *>::iterator it=clouds.begin(); it!=clouds.end(); it++)
		delete *it;
}


void Scene::init()
{
	initShaders();

	currentTime = 0.0f;
	
	camera.init(2.0f);
	
	bDraw1 = true;
	bDraw2 = true;
	bDrawLines = true;
	pointSize = 5;
}

bool Scene::loadScans(const char *filename1, const char *filename2)
{
	ifstream fin;
	
	// Check if files exist
	fin.open(filename1);
	if(!fin.is_open())
		return false;
	fin.open(filename2);
	if(!fin.is_open())
		return false;
	
	// Loading point clouds
	PointCloud *cloud;
	
	cloud = new PointCloud();
	if(cloud->load(filename1))
	{
		cloud->changeAllColors(glm::vec4(0.15f, 0.15f, 0.85f, 1.0f), pointsProgram);
		addCloud(cloud);
		cout << "Point cloud " << filename1 << " loaded. " << cloud->size() << " points." << endl;
	}
	else
		delete cloud;

	cloud = new PointCloud();
	if(cloud->load(filename2))
	{
		cloud->changeAllColors(glm::vec4(0.15f, 0.85f, 0.85f, 1.0f), pointsProgram);
		addCloud(cloud);
		cout << "Point cloud " << filename2 << " loaded. " << cloud->size() << " points." << endl;
	}
	else
		delete cloud;
		
	// Compute global bounding box
	glm::vec3 bbox[2], center, size;
	glm::mat4 matrix;
	
	getBBox(bbox);
	center = (bbox[0] + bbox[1]) / 2.0f;
	size = bbox[1] - bbox[0];
	size = glm::vec3(glm::max(size.x, glm::max(size.y, size.z)));
	matrix = glm::mat4(1.0f);
	matrix = glm::scale(matrix, glm::vec3(1.0f, 1.0f, 1.0f) / size);
	matrix = glm::translate(matrix, -center);
	transform(matrix);

	icp.setClouds(clouds[0], clouds[1]);
	icp.markBorderPoints();
	clouds[0]->changeColorsOfIndexes(glm::vec4(1,0,0,1), pointsProgram, icp.getBorderIndexes());
	clouds[0]->sendToOpenGL(pointsProgram);
	vector<int> *correspondence = icp.computeCorrespondence();
	if(correspondence != NULL)
	{
		segments.buildFromCorrespondence(clouds[0], clouds[1], *correspondence);
		segments.sendToOpenGL(segmentsProgram);
	}
	icpTransform = icp.computeICPStep();
	
	return true;
}

void Scene::update(int deltaTime)
{
	currentTime += deltaTime;
}

void Scene::render()
{
	pointsProgram.use();
	pointsProgram.setUniform1i("bLighting", int(bLighting));
	pointsProgram.setUniformMatrix4f("projection", camera.getProjectionMatrix());
	pointsProgram.setUniformMatrix4f("modelview", camera.getModelViewMatrix());
	pointsProgram.setUniformMatrix3f("normalMatrix", glm::inverseTranspose(glm::mat3(camera.getModelViewMatrix())));
	glPointSize(float(pointSize));
	
	if(bDraw1)
		clouds[0]->render();
	if(bDraw2)
		clouds[1]->render();

	segmentsProgram.use();
	segmentsProgram.setUniform4f("color", 0.f, 0.f, 0.f, 1.f);
	segmentsProgram.setUniformMatrix4f("projection", camera.getProjectionMatrix());
	segmentsProgram.setUniformMatrix4f("modelview", camera.getModelViewMatrix());
	//glLineWidth(1.f);
	
	if(bDrawLines)
		segments.render();
}

void Scene::render_gui()
{
	ImGui::Begin("Options", 0, ImGuiWindowFlags_AlwaysAutoResize);
	ImGui::Text("Render");
	ImGui::SliderInt("Point size", &pointSize, 1, 7);
	ImGui::Checkbox("Lighting", &bLighting);
	ImGui::Checkbox("Draw cloud 1", &bDraw1);
	ImGui::Checkbox("Draw cloud 2", &bDraw2);
	ImGui::Checkbox("Draw correspondence", &bDrawLines);
	if(ImGui::Button("ICP Step"))
		applyICPStep();
	if(ImGui::Button("Full ICP"))
		applyFullICP();
	ImGui::End();
}

Camera &Scene::getCamera()
{
	return camera;
}

void Scene::initShaders()
{
	Shader vShader, fShader;

	// Load shader for point rendering
	vShader.initFromFile(VERTEX_SHADER, "shaders/points.vert");
	if(!vShader.isCompiled())
	{
		cout << "Vertex Shader Error" << endl;
		cout << "" << vShader.log() << endl << endl;
	}
	fShader.initFromFile(FRAGMENT_SHADER, "shaders/points.frag");
	if(!fShader.isCompiled())
	{
		cout << "Fragment Shader Error" << endl;
		cout << "" << fShader.log() << endl << endl;
	}
	pointsProgram.init();
	pointsProgram.addShader(vShader);
	pointsProgram.addShader(fShader);
	pointsProgram.link();
	if(!pointsProgram.isLinked())
	{
		cout << "Shader Linking Error" << endl;
		cout << "" << pointsProgram.log() << endl << endl;
	}
	pointsProgram.bindFragmentOutput("outColor");
	vShader.free();
	fShader.free();

	// Load shader for segment rendering
	vShader.initFromFile(VERTEX_SHADER, "shaders/unlit.vert");
	if(!vShader.isCompiled())
	{
		cout << "Vertex Shader Error" << endl;
		cout << "" << vShader.log() << endl << endl;
	}
	fShader.initFromFile(FRAGMENT_SHADER, "shaders/unlit.frag");
	if(!fShader.isCompiled())
	{
		cout << "Fragment Shader Error" << endl;
		cout << "" << fShader.log() << endl << endl;
	}
	segmentsProgram.init();
	segmentsProgram.addShader(vShader);
	segmentsProgram.addShader(fShader);
	segmentsProgram.link();
	if(!segmentsProgram.isLinked())
	{
		cout << "Shader Linking Error" << endl;
		cout << "" << segmentsProgram.log() << endl << endl;
	}
	segmentsProgram.bindFragmentOutput("outColor");
	vShader.free();
	fShader.free();
}

void Scene::addCloud(PointCloud *cloud)
{
	clouds.push_back(cloud);
	cloud->sendToOpenGL(pointsProgram);
}

void Scene::getBBox(glm::vec3 bbox[2]) const
{
	if(clouds.size() == 0)
		return;

	glm::vec3 cloudBBox[2];
	
	clouds[0]->getBBox(bbox);
	for(unsigned int i=1; i<clouds.size(); i++)
	{
		clouds[i]->getBBox(cloudBBox);
		bbox[0] = glm::min(bbox[0], cloudBBox[0]);
		bbox[1] = glm::max(bbox[1], cloudBBox[1]);
	}
}

void Scene::transform(const glm::mat4 &matrix)
{
	for(unsigned int i=0; i<clouds.size(); i++)
	{
		clouds[i]->transform(matrix);
		clouds[i]->sendToOpenGL(pointsProgram);
	}
}

void Scene::applyICPStep()
{
	clouds[1]->transform(icpTransform);
	clouds[1]->sendToOpenGL(pointsProgram);
	vector<int> *correspondence = icp.computeCorrespondence();
	if(correspondence != NULL)
	{
		segments.buildFromCorrespondence(clouds[0], clouds[1], *correspondence);
		segments.sendToOpenGL(segmentsProgram);
	}
	icpTransform = icp.computeICPStep();
}

void Scene::applyFullICP()
{
	vector<int> *correspondence = icp.computeFullICP();
	clouds[1]->sendToOpenGL(pointsProgram);
	if(correspondence != NULL)
	{
		segments.buildFromCorrespondence(clouds[0], clouds[1], *correspondence);
		segments.sendToOpenGL(segmentsProgram);
	}
}









