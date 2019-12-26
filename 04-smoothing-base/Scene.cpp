#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include "imgui/imgui.h"
#include "Scene.h"
#include "LaplacianSmoothing.h"


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
	bWireframe = false;
	lambda = 1.f;
	selectedIterativeFunction = 0;
	selectedGlobalFunction = 0;
	percentConstraints = 100;
	constraintWeight = 1.0f;
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
	
	// Save original vertices for reload
	originalVertices.clear();
	copy(mesh->getVertices().begin(), mesh->getVertices().end(), std::back_inserter(originalVertices));
	
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
	pointsProgram.setUniform1i("bBlack", false);
	pointsProgram.setUniformMatrix4f("projection", camera.getProjectionMatrix());
	pointsProgram.setUniformMatrix4f("modelview", camera.getModelViewMatrix());
	pointsProgram.setUniformMatrix3f("normalMatrix", glm::inverseTranspose(glm::mat3(camera.getModelViewMatrix())));
	
	if(bWireframe)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}
	if(mesh != NULL)
		mesh->render();
	if(bWireframe)
	{
		glDisable(GL_POLYGON_OFFSET_FILL);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		pointsProgram.setUniform1i("bBlack", true);
		mesh->render();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

void Scene::render_gui()
{
	ImGui::Begin("Options", 0, ImGuiWindowFlags_AlwaysAutoResize);
	
	// Model options
	ImGui::Text("Model");
	ImGui::Spacing();
	if(ImGui::Button("Reload"))
		reloadMesh();
	ImGui::Spacing();

	// Render options
	ImGui::Text("Render");
	ImGui::Spacing();
	ImGui::Checkbox("Lighting", &bLighting);
	ImGui::Spacing();
	ImGui::Checkbox("Wireframe", &bWireframe);
	ImGui::Spacing();

	ImGui::Separator();
	ImGui::Spacing();

	// Iterative smoothing
	ImGui::Text("Iterative Smoothing");
	ImGui::Spacing();
	
	ImGui::Text("Lambda: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputFloat("##Lambda", &lambda, 0.01f, 0.1f, "%.2f");
	lambda = glm::max(0.f, glm::min(lambda, 1.f));
	ImGui::Spacing();

	ImGui::Text("#Iterations: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(75);
	ImGui::InputInt("##Iterations", &nIterations, 1, 10);
	nIterations = glm::max(1, glm::min(nIterations, 100));
	ImGui::Spacing();
	ImGui::RadioButton("Laplacian", &selectedIterativeFunction, 0); ImGui::SameLine();
	ImGui::RadioButton("Bilaplacian", &selectedIterativeFunction, 1); ImGui::SameLine();
	ImGui::RadioButton("Lambda-Nu", &selectedIterativeFunction, 2);
	ImGui::Spacing();
	if(ImGui::Button("Apply"))
	{
		LaplacianSmoothing smoothing;
		
		smoothing.setMesh(mesh);
		switch(selectedIterativeFunction)
		{
		case 0:
			smoothing.iterativeLaplacian(nIterations, lambda);
			break;
		case 1:
			smoothing.iterativeBilaplacian(nIterations, lambda);
			break;
		case 2:
			smoothing.iterativeLambdaNu(nIterations, lambda);
			break;
		}
		mesh->sendToOpenGL(pointsProgram);
	}
	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();

	// Global smoothing
	ImGui::Text("Global Smoothing");
	ImGui::Spacing();

	ImGui::Text("Constraints: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputInt("##Percentage", &percentConstraints, 1, 10);
	percentConstraints = glm::max(1, glm::min(percentConstraints, 100));
	ImGui::Spacing();

	ImGui::Text("Weight: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputFloat("##Weight", &constraintWeight, 0.01f, 0.1f, "%.2f");
	constraintWeight = glm::max(0.f, glm::min(constraintWeight, 1.f));
	ImGui::Spacing();

	ImGui::RadioButton("Laplacian##2", &selectedGlobalFunction, 0); ImGui::SameLine();
	ImGui::RadioButton("Bilaplacian##2", &selectedGlobalFunction, 1);
	ImGui::Spacing();
	if(ImGui::Button("Apply##2"))
	{
		vector<unsigned int> selected;
		vector<bool> constraints;
		LaplacianSmoothing smoothing;
		
		srand(time(NULL));
		selected.resize(mesh->getVertices().size());
		for(unsigned int i=0; i<mesh->getVertices().size(); i++)
			selected[i] = i;
		random_shuffle(selected.begin(), selected.end());
		selected.resize(selected.size() * percentConstraints / 100);
		constraints.resize(mesh->getVertices().size(), false);
		for(unsigned int i=0; i<selected.size(); i++)
			constraints[selected[i]] = true;
		smoothing.setMesh(mesh);
		switch(selectedGlobalFunction)
		{
		case 0:
			smoothing.globalLaplacian(constraints);
			break;
		case 1:
			smoothing.globalBilaplacian(constraints, constraintWeight);
			break;
		}
		mesh->sendToOpenGL(pointsProgram);
	}

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

void Scene::reloadMesh()
{
	mesh->getVertices().clear();
	copy(originalVertices.begin(), originalVertices.end(), std::back_inserter(mesh->getVertices()));
	mesh->sendToOpenGL(pointsProgram);
}











