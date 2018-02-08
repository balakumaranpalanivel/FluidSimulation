//#include "Controller.h"
//#include "Camera.h"
//#include "Configuration.h"
//#include <glm/glm.hpp>
//
//CController::CController()
//{
//	mScreenWidth = Screen::WIDTH;
//	mScreenHeight = Screen::HEIGHT;
//
//	glm::vec3 pos = glm::vec3(14.0, 15.0, 33.0);
//	glm::vec3 dir = glm::normalize(glm::vec3(-pos.x, -pos.y, -pos.z));
//	mCamera = CCamera(pos, dir, mScreenWidth, mScreenHeight, 60.0, 0.5, 100.0);
//}
//
//CController::~CController()
//{
//
//}
//
//void CController::InitializeGL()
//{
//	static const GLfloat lightPos[4] = { 20.0f, 20.0f, 20.0f, 1.0f };
//	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
//	glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT0);
//	glEnable(GL_DEPTH_TEST);
//	glEnable(GL_COLOR_MATERIAL);
//	glEnable(GL_NORMALIZE);
//	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
//}
//
//void CController::PaintGL()
//{
//	mCamera.set();
//
//	float scale = 2.0;
//
//	glColor3f(0.0, 0.0, 0.0);
//	glPointSize(6.0);
//	
//	mCamera.unset();
//}

/****************************************************************************
**
** Copyright (C) 2013 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/
#include <GL\glew.h>
#include <string>
#include "BMPLoader.h"

//#include <QtWidgets\qopenglwidget.h>
//#include <QtGui\qcursor.h>
//#include <QtGui\QKeyEvent>
//#include <QtGui\qopengltexture.h>
//#include <QtCore\qtime>
//#include <QtCore\qtimer.h>
#include <algorithm>
#include "Controller.h"
#include "glm/glm.hpp"
#include "quaternion.h"
#include "Gradient.h"
#include <iostream>

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//! [0]
GLWidget::GLWidget()
{

	screenWidth = 1280;
	screenHeight = 720;

	//// update timer
	//float updatesPerSecond = 60;
	//updateTimer = new QTimer(this);
	//connect(updateTimer, SIGNAL(timeout()), this, SLOT(updateSimulation()));
	//updateTimer->start(1000.0 / updatesPerSecond);

	//deltaTimer = new QTime();
	//deltaTimer->start();

	// Initialize camera
	glm::vec3 pos = glm::vec3(14.0, 15.0, 33.0);
	glm::vec3 dir = glm::normalize(glm::vec3(-pos.x, -pos.y, -pos.z));
	camera = CCamera(pos, dir, screenWidth, screenHeight,
		60.0, 0.5, 100.0);

	// for updating camera movement
	isMovingForward = false;
	isMovingBackward = false;
	isMovingRight = false;
	isMovingLeft = false;
	isMovingUp = false;
	isMovingDown = false;
	isRotatingRight = false;
	isRotatingLeft = false;
	isRotatingUp = false;
	isRotatingDown = false;
	isRollingRight = false;
	isRollingLeft = false;

	// simulation system
	minDeltaTimeModifier = 0.125;
	maxDeltaTimeModifier = 1.0;
	deltaTimeModifier = maxDeltaTimeModifier;
	runningTime = 0.0;
	currentFrame = 0;

	initializeSimulation();
	fluidSim.SetCamera(&camera);

}

GLWidget::~GLWidget()
{
}

//QSize GLWidget::minimumSizeHint() const
//{
//	return QSize(50, 50);
//}

//QSize GLWidget::sizeHint() const
//{
//	return QSize(screenWidth, screenHeight);
//}

void GLWidget::initializeGL()
{
	//QOpenGLTexture textureBall(QImage(QString("images/blankball.png")));
	//texture[0] = textureBall.textureId();

	//glBindTexture(texture[0], GL_TEXTURE_2D);
	//fluidSim.SetTexture(&texture[0]);

	//GLuint temp = 0;
	//fluidSim.SetTexture(&temp);

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


// move and rotate camera
void GLWidget::updateCameraMovement(float dt) {
	float camSpeed = 10.0;
	float camRotSpeed = 2.0;

	if (isMovingForward) { camera.moveForward(camSpeed * dt); }
	if (isMovingBackward) { camera.moveBackward(camSpeed * dt); }
	if (isMovingRight) { camera.moveRight(camSpeed * dt); }
	if (isMovingLeft) { camera.moveLeft(camSpeed * dt); }
	if (isMovingUp) { camera.moveUp(camSpeed * dt); }
	if (isMovingDown) { camera.moveDown(camSpeed * dt); }
	if (isRotatingRight) { camera.rotateRight(camRotSpeed * dt); }
	if (isRotatingLeft) { camera.rotateLeft(camRotSpeed * dt); }
	if (isRotatingUp) { camera.rotateUp(camRotSpeed * dt); }
	if (isRotatingDown) { camera.rotateDown(camRotSpeed * dt); }
	if (isRollingRight) { camera.rollRight(camRotSpeed * dt); }
	if (isRollingLeft) { camera.rollLeft(camRotSpeed * dt); }
}

void GLWidget::initializeSimulation()
{

	simulationFPS = mSimulationConfig.FPS;
	isSimulationPaused = mSimulationConfig.IS_SIMULATION_PAUSED;

	double radius = mSimulationConfig.SMOOTHING_RADIUS;
	fluidSim = CSPHFluidSimulation(radius);

	double minx = 0.0;
	double maxx = 14.0;
	double miny = 0.0;
	double maxy = 8.0;
	double minz = 0.0;
	double maxz = 6.0;
	double rx = 0.3;
	double ry = 0.8;
	double rz = 1.0;

	fluidSim.SetBounds(minx, maxx, miny, maxy, minz, maxz);
	std::vector<glm::vec3> points;
	double n = mSimulationConfig.NUM_PARTICLES;
	for (int i = 0; i<n; i++) {
		float x = minx + rx*((float)rand() / RAND_MAX) * (maxx - minx);
		float y = miny + ry*((float)rand() / RAND_MAX) * (maxy - miny);
		float z = minz + rz*((float)rand() / RAND_MAX) * (maxz - minz);
		glm::vec3 p = glm::vec3(x, y, z);
		points.push_back(p);
	}
	fluidSim.AddFluidParticles(points);

	double damp = mSimulationConfig.INITIAL_DAMPING_CONSTANT;
	fluidSim.SetDampingConstant(damp);

	fluidGradient = Gradients::GetSkyblueGradient();
	minColorDensity = mGraphicsConfig.MIN_COLOR_DENSITY;
	maxColorDensity = mGraphicsConfig.MAX_COLOR_DENSITY;

	// create obstacles
	/*
	float widthPad = 0.2;
	float heightPad = 0.4;
	float width = (maxz - minz) / 3 + widthPad;
	float height = (maxy - miny) + heightPad;
	float depth = 4.2;
	float r = 0.6*radius;
	int layers = 6;
	bool isStaggered = true;

	glm::vec3 p1 = glm::vec3(depth, height / 2, ((maxz - minz) / 3) - ((maxz - minz) / 6));
	glm::vec3 p2 = glm::vec3(depth, height / 2, (2 * (maxz - minz) / 3) - ((maxz - minz) / 6));
	glm::vec3 p3 = glm::vec3(depth, height / 2, (3 * (maxz - minz) / 3) - ((maxz - minz) / 6));

	std::vector<glm::vec3> o1, o2, o3;
	o1 = Utils::CreatePointPanel(width, height, r, layers,
		glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 1.0, 0.0),
		isStaggered);
	o2 = Utils::CreatePointPanel(width, height + 1.0, r, layers,
		glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 1.0, 0.0),
		isStaggered);
	o3 = Utils::CreatePointPanel(width, height, r, layers,
		glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 1.0, 0.0),
		isStaggered);

	o1 = Utils::TranslatePoints(o1, p1);
	//o2 = utils::translatePoints(o2, p2);
	o3 = Utils::TranslatePoints(o3, p3);

	int id1 = fluidSim.AddObstacleParticles(o1);
	int id2 = fluidSim.AddObstacleParticles(o2);
	int id3 = fluidSim.AddObstacleParticles(o3);
	//fluidSim.translateObstacle(id1, p1);
	fluidSim.TranslateObstacle(id2, p2);
	//fluidSim.translateObstacle(id3, p3);
	*/

}

void GLWidget::activateSimulation() {

	double minx = 0.0;
	double maxx = 14.0;
	double miny = 0.0;
	double maxy = 8.0;
	double minz = 0.0;
	double maxz = 6.0;
	fluidSim.SetBounds(minx, maxx, miny, maxy, minz, maxz);

	double damp = mSimulationConfig.FINAL_DAMPING_CONSTANT;
	fluidSim.SetDampingConstant(damp);
	isRendering = true;
}

void GLWidget::stopSimulation() {
}

void GLWidget::updateSimulationSettings() {

	minColorDensity = mGraphicsConfig.MIN_COLOR_DENSITY;
	maxColorDensity = mGraphicsConfig.MAX_COLOR_DENSITY;
	simulationFPS = mSimulationConfig.FPS;
	isSimulationPaused = mSimulationConfig.IS_SIMULATION_PAUSED;
	isSimulationDrawn = mSimulationConfig.IS_SIMULATION_DRAWN;
}

double last_time1 = 0.0f;

void GLWidget::updateSimulation(float dt) {

	// find delta time
	double  curr_time = timeGetTime();
	float  delta = (curr_time - last_time1) * 0.001f;

	last_time1 = curr_time;
	updateCameraMovement(delta);
	fluidSim.SetCamera(&camera);

	dt = 1.0 / simulationFPS;
	if (isSimulationPaused) {
		dt = 0.0;
	}
	runningTime = runningTime + dt;

	bool isRenderingEnabled = mSimulationConfig.IS_RENDERING_ENABLED;

	if (!isSimulationPaused) {
		if (isRendering) {
			float speed = -1.33;
			fluidSim.TranslateObstacle(1, glm::vec3(0.0, speed*dt, 0.0));
		}
		fluidSim.Update(dt);
	}


	if (isRendering && isRenderingEnabled && !isSimulationPaused) {
		writeFrame();
	}
	else {
		writeFrame();
	}
}

void GLWidget::writeFrame() {
	paintGL();

	// image file name
	std::string s = std::to_string(currentFrame);
	s.insert(s.begin(), 6 - s.size(), '0');
	//s = "test_render/" + s + ".png";
	//bool r = saveFrameToFile(s + ".bmp");

	//qDebug() << r << QString::fromStdString(s);

	//// log timings to file
	//using namespace luabridge;
	//lua_State* L = luaL_newstate();
	//luaL_openlibs(L);

	//if (luaL_dofile(L, "scripts/fluid_config.lua") != 0) {
	//	qDebug() << "Error loading script";
	//	exit(1);
	//}
	//LuaRef t = getGlobal(L, "settings");
	//std::string filename = t["logfile"].cast<std::string>();

	//QString f = QString::fromStdString(filename);
	//QString filepath = "logs/" + f;
	//QFile file(filepath);
	//if (file.open(QIODevice::Append)) {
	//	QTextStream stream(&file);
	//	stream << currentFrame << " " << fluidSim.getTimingData() << endl;
	//	file.close();
	//}
	//else {
	//	qDebug() << "error opening file: ";
	//}

	currentFrame += 1;
}

bool GLWidget::saveFrameToFile(std::string fileName) {
	//GLubyte *data = (GLubyte*)malloc(4 * (int)screenWidth*(int)screenHeight);
	//if (data) {
	//	glReadPixels(0, 0, screenWidth, screenHeight,
	//		GL_RGBA, GL_UNSIGNED_BYTE, data);
	//}

	FILE *fScreenshot;
	int nSize = screenWidth * screenHeight * 3;

	GLubyte *pixels = new GLubyte[nSize];
	if (pixels == NULL)
		return false;
	fScreenshot = fopen(fileName.c_str(), "w+");
	if (fScreenshot == NULL)
		return false;
	else 
		fclose(fScreenshot);

	fScreenshot = fopen(fileName.c_str(), "r+");

	glReadPixels(0, 0, screenWidth, screenHeight, GL_RGB,
		GL_UNSIGNED_BYTE, pixels);

	BMPLoader SaveBmp;
	SaveBmp.SaveBMP(fileName.c_str(), screenWidth, screenHeight, pixels);
	////convert to BGR format    
	//GLubyte temp;
	//int i = 0;
	//while (i < nSize)
	//{
	//	temp = pixels;       //grab blue
	//	pixels = pixels[i + 2];//assign red to blue
	//	pixels[i + 2] = temp;     //assign blue to red

	//	i += 3;     //skip to next blue byte
	//}

	//unsigned char TGAheader[12] = { 0,0,2,0,0,0,0,0,0,0,0,0 };
	//unsigned char header[6] = { int(screenWidth) % 256,screenWidth / 256,
	//	int(screenHeight) % 256,screenHeight / 256,24,0 };

	//fwrite(TGAheader, sizeof(unsigned char), 12, fScreenshot);
	//fwrite(header, sizeof(unsigned char), 6, fScreenshot);
	//fwrite(pixels, sizeof(GLubyte), nSize, fScreenshot);
	//fclose(fScreenshot);

	delete[] pixels;

	return true;
	//convert to BGR format    
	//GLubyte temp;
	//int i = 0;
	//while (i < nSize)
	//{
	//	temp = pixels;       //grab blue
	//	pixels = pixels[i + 2];//assign red to blue
	//	pixels[i + 2] = temp;     //assign blue to red

	//	i += 3;     //skip to next blue byte
	//}

	//QImage image((int)screenWidth, (int)screenHeight, QImage::Format_RGB32);
	//for (int j = 0; j<screenHeight; j++) {
	//	for (int i = 0; i<screenWidth; i++) {
	//		int idx = 4 * (j*screenWidth + i);
	//		char r = data[idx + 0];
	//		char g = data[idx + 1];
	//		char b = data[idx + 2];

	//		// sets 32 bit pixel at (x,y) to yellow.
	//		//uchar *p = image.scanLine(j) + i;
	//		//*p = qRgb(255, 0, 0);
	//		QRgb value = qRgb(r, g, b);
	//		image.setPixel(i, screenHeight - j - 1, value);
	//	}
	//}
	//bool r = image.save(fileName);
	//free(data);

	//return r;
	return true;
}

void GLWidget::paintGL()
{
	camera.set();
	Utils::DrawGrid();

	float scale = 2.0;
	glm::mat4 scaleMat = glm::transpose(glm::mat4(scale, 0.0, 0.0, 0.0,
		0.0, scale, 0.0, 0.0,
		0.0, 0.0, scale, 0.0,
		0.0, 0.0, 0.0, 1.0));
	glPushMatrix();
	glMultMatrixf((GLfloat*)&scaleMat);

	if (isSimulationDrawn) {
		fluidSim.Draw();
	}
	else {
		fluidSim.DrawBounds();
	}

	glColor3f(0.0, 0.0, 0.0);
	glPointSize(6.0);

	glPopMatrix();

	camera.unset();
}

void GLWidget::resizeGL(int width, int height)
{
	screenWidth = width;
	screenHeight = height;
	camera.resize((float)width, (float)height);
	paintGL();
}

//void GLWidget::keyPressEvent(QKeyEvent *event)
//{
//	// turn on camera movement
//	isMovingForward = event->key() == Qt::Key_W;
//	isMovingBackward = event->key() == Qt::Key_S;
//	isMovingRight = event->key() == Qt::Key_D;
//	isMovingLeft = event->key() == Qt::Key_A;
//	isMovingUp = event->key() == Qt::Key_T;
//	isMovingDown = event->key() == Qt::Key_G;
//
//	isRotatingRight = event->key() == Qt::Key_E;
//	isRotatingLeft = event->key() == Qt::Key_Q;
//	isRotatingUp = event->key() == Qt::Key_F;
//	isRotatingDown = event->key() == Qt::Key_R;
//	isRollingRight = event->key() == Qt::Key_X;
//	isRollingLeft = event->key() == Qt::Key_Z;
//
//	// slow down simulation
//	if (event->key() == Qt::Key_C) {
//		deltaTimeModifier = minDeltaTimeModifier;
//	}
//}

//void GLWidget::keyReleaseEvent(QKeyEvent *event)
//{
//	// turn off camera movement
//	isMovingForward = event->key() == Qt::Key_W ? false : isMovingForward;
//	isMovingBackward = event->key() == Qt::Key_S ? false : isMovingBackward;
//	isMovingRight = event->key() == Qt::Key_D ? false : isMovingRight;
//	isMovingLeft = event->key() == Qt::Key_A ? false : isMovingLeft;
//	isMovingUp = event->key() == Qt::Key_T ? false : isMovingUp;
//	isMovingDown = event->key() == Qt::Key_G ? false : isMovingDown;
//
//	isRotatingRight = event->key() == Qt::Key_E ? false : isRotatingRight;
//	isRotatingLeft = event->key() == Qt::Key_Q ? false : isRotatingLeft;
//	isRotatingDown = event->key() == Qt::Key_R ? false : isRotatingDown;
//	isRotatingUp = event->key() == Qt::Key_F ? false : isRotatingUp;
//	isRollingRight = event->key() == Qt::Key_X ? false : isRollingRight;
//	isRollingLeft = event->key() == Qt::Key_Z ? false : isRollingLeft;
//
//	if (event->key() == Qt::Key_C) {
//		deltaTimeModifier = maxDeltaTimeModifier;
//	}
//
//	if (event->key() == Qt::Key_H) {
//		activateSimulation();
//	}
//
//}










