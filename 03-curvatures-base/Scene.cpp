#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include "imgui/imgui.h"
#include "Scene.h"
#include "MongePatch.h"
#include "NearestNeighbors.h"


Scene::Scene()
{
	mesh = NULL;
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
	bCurvature = false;
	selectedCurvature = 2;
}

bool Scene::loadScan(const char *filename)
{
	ifstream fin;
	
	// Check if file exists
	fin.open(filename);
	if(!fin.is_open())
		return false;
	
	if(mesh != NULL)
		delete mesh;
	mesh = new TriangleMesh();
	if(mesh->load(filename))
	{
		cout << "Triangle mesh " << filename << " loaded. ";
		cout << mesh->getVertices().size() << " vertices. ";
		cout << (mesh->getTriangles().size()/3) << " faces." << endl;
	}
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
	
	computeCurvatures();
	greyColors();
	
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
	
	if(mesh != NULL)
		mesh->render();
}

void Scene::render_gui()
{
	ImGui::Begin("Options", 0, ImGuiWindowFlags_AlwaysAutoResize);

	ImGui::Text("Render");
	ImGui::Spacing();
	ImGui::Checkbox("Lighting", &bLighting);
	ImGui::Spacing();
	ImGui::Separator();

	bool prevUseCurvature = bCurvature;
	int prevCurvature = selectedCurvature;
	
	ImGui::Spacing();
	ImGui::Checkbox("Curvatures", &bCurvature);
	ImGui::Spacing();
	ImGui::RadioButton("Min", &selectedCurvature, 0); ImGui::SameLine();
	ImGui::RadioButton("Max", &selectedCurvature, 1); ImGui::SameLine();
	ImGui::RadioButton("Mean", &selectedCurvature, 2); ImGui::SameLine();
	ImGui::RadioButton("Gaussian", &selectedCurvature, 3);
	if(prevUseCurvature != bCurvature)
	{
		if(bCurvature)
			changeCurvatureColors();
		else
			greyColors();
	}
	else if(bCurvature && (prevCurvature != selectedCurvature))
		changeCurvatureColors();

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
	mesh->getBBox(bbox);
}

void Scene::transform(const glm::mat4 &matrix)
{
	mesh->transform(matrix);
	mesh->sendToOpenGL(pointsProgram);
}

void Scene::computeCurvatures()
{
	const vector<glm::vec3> &vertices = mesh->getVertices();
	const vector<glm::vec3> &normals = mesh->getPerVertexNormals();
	NearestNeighbors knn;
	vector<size_t> neighbors;
	vector<float> dists_squared;
	vector<glm::vec3> neighborsCoords;
	vector<float> curvatures[2];

	kmin.resize(vertices.size());
	kmax.resize(vertices.size());
	knn.setPoints(&vertices);
	for(unsigned int i=0; i<vertices.size(); i++)
	{
		MongePatch patch;
		float curv1, curv2;

		neighbors.clear();
		dists_squared.clear();
		neighborsCoords.clear();
		knn.getKNearestNeighbors(vertices[i], 21, neighbors, dists_squared);
		for(unsigned int j=0; j<neighbors.size(); j++)
			neighborsCoords.push_back(vertices[neighbors[j]]);
		patch.init(vertices[i], normals[i], neighborsCoords);
		patch.principalCurvatures(curv1, curv2);
		kmin[i] = curv1;
		kmax[i] = curv2;
		curvatures[0].push_back(fabsf(curv1));
		curvatures[1].push_back(fabsf(curv2));
	}
	sort(curvatures[0].begin(), curvatures[0].end());
	maxCurvature[0] = curvatures[0][int(0.9f * curvatures[0].size())];
	sort(curvatures[1].begin(), curvatures[1].end());
	maxCurvature[1] = curvatures[1][int(0.9f * curvatures[1].size())];
}

void Scene::changeCurvatureColors()
{
	const vector<glm::vec3> &vertices = mesh->getVertices();
	vector<glm::vec4> &colors = mesh->getColors();
	float curvature, scale;
	glm::vec4 color;

	for(unsigned int i=0; i<vertices.size(); i++)
	{
		switch(selectedCurvature)
		{
		case 0:
			curvature = kmin[i];
			scale = maxCurvature[0];
			break;
		case 1:
			curvature = kmax[i];
			scale = maxCurvature[1];
			break;
		case 2:
			curvature = (kmin[i] + kmax[i]) / 2.f;
			scale = (maxCurvature[0] + maxCurvature[1]) / 2.f;
			break;
		case 3:
			curvature = kmin[i] * kmax[i];
			scale = maxCurvature[0] * maxCurvature[1];
			break;
		}
		curvature = glm::clamp((curvature + scale) / (2.f * scale), 0.f, 1.f);
		if(curvature < 0.25f)
			color = glm::mix(glm::vec4(0.f, 0.f, 1.f, 1.f), glm::vec4(0.f, 1.f, 1.f, 1.f), 4.f * curvature);
		else if(curvature < 0.5f)
			color = glm::mix(glm::vec4(0.f, 1.f, 1.f, 1.f), glm::vec4(0.f, 1.f, 0.f, 1.f), 4.f * curvature - 1.f);
		else if(curvature < 0.75f)
			color = glm::mix(glm::vec4(0.f, 1.f, 0.f, 1.f), glm::vec4(1.f, 1.f, 0.f, 1.f), 4.f * curvature - 2.f);
		else
			color = glm::mix(glm::vec4(1.f, 1.f, 0.f, 1.f), glm::vec4(1.f, 0.f, 0.f, 1.f), 4.f * curvature - 3.f);
		colors[i] = color;
	}
	mesh->sendToOpenGL(pointsProgram);
}

void Scene::greyColors()
{
	const vector<glm::vec3> &vertices = mesh->getVertices();
	vector<glm::vec4> &colors = mesh->getColors();

	colors.clear();
	colors.resize(vertices.size(), glm::vec4(0.5f, 0.5f, 0.5f, 1.0f));
	mesh->sendToOpenGL(pointsProgram);
}










