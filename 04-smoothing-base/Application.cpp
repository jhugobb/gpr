#include <iostream>
#include <cstdlib>
#include <cstring>
#include "Application.h"
#include <GL/freeglut.h>
#include <glm/gtc/matrix_transform.hpp>


void Application::init()
{
	bPlay = true;
	glClearColor(1.f, 1.f, 1.f, 1.0f);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	scene.init();
	
	for(unsigned int i=0; i<256; i++)
	{
		keys[i] = false;
		specialKeys[i] = false;
	}
	mouseButtons[0] = false;
	mouseButtons[1] = false;
	lastMousePos = glm::ivec2(-1, -1);
}

bool Application::loadScan(const char *filename)
{
	return scene.loadScan(filename);
}

bool Application::update(int deltaTime)
{
	scene.update(deltaTime);
	
	return bPlay;
}

void Application::render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	scene.render();
}

void Application::render_gui()
{
	scene.render_gui();
}

void Application::resize(int width, int height)
{
	glViewport(0, 0, width, height);
	scene.getCamera().resizeCameraViewport(width, height);
}

void Application::keyPressed(int key)
{
	if(key == 27) // Escape code
		bPlay = false;
	keys[key] = true;
}

void Application::keyReleased(int key)
{
	keys[key] = false;
}

void Application::specialKeyPressed(int key)
{
	specialKeys[key] = true;
}

void Application::specialKeyReleased(int key)
{
	specialKeys[key] = false;
}

void Application::mouseMove(int x, int y)
{
	// Rotation
	if(mouseButtons[0] && lastMousePos.x != -1)
		scene.getCamera().rotateCamera(0.5f * (y - lastMousePos.y), 0.5f * (x - lastMousePos.x));

	// Zoom
	if(mouseButtons[1] && lastMousePos.x != -1)
		scene.getCamera().zoomCamera(0.01f * (y - lastMousePos.y));

 	lastMousePos = glm::ivec2(x, y);
}

void Application::mousePress(int button)
{
	mouseButtons[button] = true;
}

void Application::mouseRelease(int button)
{
	mouseButtons[button] = false;
	if(!mouseButtons[0] && !mouseButtons[1])
		lastMousePos = glm::ivec2(-1, -1);
}

bool Application::getKey(int key) const
{
	return keys[key];
}

bool Application::getSpecialKey(int key) const
{
	return specialKeys[key];
}





