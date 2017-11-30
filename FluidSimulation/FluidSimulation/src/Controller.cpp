#include "Controller.h"
#include "Camera.h"
#include "Configuration.h"
#include <glm/glm.hpp>

CController::CController()
{
	mScreenWidth = Screen::WIDTH;
	mScreenHeight = Screen::HEIGHT;

	glm::vec3 pos = glm::vec3(14.0, 15.0, 33.0);
	glm::vec3 dir = glm::normalize(glm::vec3(-pos.x, -pos.y, -pos.z));
	mCamera = CCamera(pos, dir, mScreenWidth, mScreenHeight, 60.0, 0.5, 100.0);
}

CController::~CController()
{

}

void CController::InitializeGL()
{
	static const GLfloat lightPos[4] = { 20.0f, 20.0f, 20.0f, 1.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void CController::PaintGL()
{
	mCamera.set();

	float scale = 2.0;

	glColor3f(0.0, 0.0, 0.0);
	glPointSize(6.0);
	
	mCamera.unset();
}

