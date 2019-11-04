#include <iostream>
#include <fstream>
#include <cmath>
#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include "imgui/imgui.h"
#include "Scene.h"
#include "SphereFunction.h"
#include "SimpleDistance.h"
#include "DualMarchingCubes.h"
#include "RBFFunction.h"


Scene::Scene()
{
	mesh = NULL;
	selectedFunction = 0;
	dmcResolution = 1;
}

Scene::~Scene()
{
	if(mesh != NULL)
	{
		mesh->free();
		delete mesh;
	}
}


void Scene::init()
{
	initShaders();

	currentTime = 0.0f;
	
	camera.init(2.0f);
	
	bLighting = true;
	bDrawPointCloud = true;
	bDrawMesh = true;
	pointSize = 5;
}

bool Scene::loadScan(const char *filename)
{
	ifstream fin;
	
	// Check if file exists
	fin.open(filename);
	if(!fin.is_open())
		return false;
	
	if(cloud.load(filename))
		cout << "Point cloud " << filename << " loaded. " << cloud.size() << " points." << endl;
	else
		return false;

	// Compute global bounding box
	glm::vec3 bbox[2], center, size;
	glm::mat4 matrix;
	
	getBBox(bbox);
	center = (bbox[0] + bbox[1]) / 2.0f;
	size = bbox[1] - bbox[0];
	size = glm::vec3(glm::max(size.x, glm::max(size.y, size.z)));
	matrix = glm::mat4(1.0f);
	matrix = glm::scale(matrix, glm::vec3(1.0f - 0.1f, 1.0f - 0.1f, 1.0f - 0.1f) / size);
	matrix = glm::translate(matrix, -center);
	transform(matrix);
	
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
	
	if(bDrawPointCloud)
		cloud.render();
	
	if(bDrawMesh && mesh != NULL)
		mesh->render();
}

void Scene::render_gui()
{
	ImGui::Begin("Options", 0, ImGuiWindowFlags_AlwaysAutoResize);

	ImGui::Text("Render");
	ImGui::Spacing();
	ImGui::SliderInt("Point size", &pointSize, 1, 7);
	ImGui::Checkbox("Lighting", &bLighting);
	ImGui::Checkbox("Draw point cloud", &bDrawPointCloud);
	ImGui::Separator();
	
	ImGui::Text("Dual Marching Cubes");
	ImGui::Spacing();
	ImGui::RadioButton("Sphere", &selectedFunction, 0); ImGui::SameLine();
	ImGui::RadioButton("Simple", &selectedFunction, 1); ImGui::SameLine();
	ImGui::RadioButton("RBF", &selectedFunction, 2);
	ImGui::SliderInt("", &dmcResolution, 0, 3, ""); ImGui::SameLine();
	ImGui::Text("Resolution: %d", 1 << (dmcResolution+5));
	ImGui::Checkbox("Draw extracted mesh", &bDrawMesh);
	ImGui::Separator();

	ImGui::Spacing();
	if(ImGui::Button("Extract mesh"))
		extractMesh();

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
}

void Scene::getBBox(glm::vec3 bbox[2]) const
{
	cloud.getBBox(bbox);
}

void Scene::transform(const glm::mat4 &matrix)
{
	cloud.transform(matrix);
	cloud.sendToOpenGL(pointsProgram);
}

void Scene::extractMesh()
{
	SphereFunction sphereDistanceFunction;
	SimpleDistance simpleDistanceFunction;
	DualMarchingCubes mc;
	RBFFunction rbfFunction;

	if(mesh != NULL)
	{
		mesh->free();
		delete mesh;
		mesh = NULL;
	}
	
	float cloudRadius = cloud.computeDistanceBetweenPoints();

	switch(selectedFunction)
	{
	case 0:
		sphereDistanceFunction.init(glm::vec3(0.f), 0.25f);
		mesh = mc.extractZeroIsosurface(sphereDistanceFunction, 1 << (dmcResolution+5));
		break;
	case 1:
		simpleDistanceFunction.init(&cloud, 0.1f);
		mesh = mc.extractZeroIsosurfaceFast(&cloud, 0.5f * cloudRadius, simpleDistanceFunction, 1 << (dmcResolution+5));
		break;
	case 2:
		rbfFunction.init(&cloud, cloudRadius, 3.0f * cloudRadius);
		mesh = mc.extractZeroIsosurfaceFast(&cloud, 0.5f * cloudRadius, rbfFunction, 1 << (dmcResolution+5));
		break;
	}
	if(mesh != NULL)
		mesh->sendToOpenGL(pointsProgram);
}









