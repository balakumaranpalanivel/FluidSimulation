/*#pragma once
#include "Camera.h"

class CController
{
public:
	CController();
	~CController();

	void InitializeGL();
	void PaintGL();

private:
	float mScreenWidth;
	float mScreenHeight;

	CCamera mCamera;
};*/	

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

#pragma once
//#include <QtWidgets\qopenglwidget.h>
#include <vector>
#include <tuple>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
//#include <QtGui\qimage.h>
#include <cmath>
#include "glm/glm.hpp"
#include "Camera.h"
#include "Utils.h"
#include "SpatialGrid.h"
#include "SPHFluidSimulation.h"
#include "Configuration.h"

/**
* @class  GLWidget
* @brief  our OpenGL view derived from QGLWidget.
* We have to override several functions for our
* application-specific OpenGL drawing functionality
*/
class GLWidget
{
public:
	GLWidget();
	~GLWidget();

	//void keyPressEvent(QKeyEvent *event);
	//void keyReleaseEvent(QKeyEvent *event);

	//QSize minimumSizeHint() const;
	//QSize sizeHint() const;

public:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void updateSimulation(float dt);

private:
	float screenWidth;
	float screenHeight;

	//void updatePerspective();
	void updateCameraMovement(float dt);
	//void drawGrid();
	//void drawAnimation();
	//bool compareByZDistance(const SPHParticle &p1, const SPHParticle &p2);
	//void drawBillboard(GLuint *tex, glm::vec3 p, float width);
	void initializeSimulation();
	void activateSimulation();
	void updateSimulationSettings();
	void stopSimulation();
	void writeFrame();
	bool saveFrameToFile();
	bool isRendering = false;

	//// update/draw tiemrs
	//QTimer *drawTimer;
	//QTimer *updateTimer;
	//QTime *deltaTimer;

	// camera and movement
	CCamera camera;
	bool isMovingForward;
	bool isMovingBackward;
	bool isMovingRight;
	bool isMovingLeft;
	bool isMovingUp;
	bool isMovingDown;
	bool isRotatingRight;
	bool isRotatingLeft;
	bool isRotatingUp;
	bool isRotatingDown;
	bool isRollingRight;
	bool isRollingLeft;

	// simulation system
	float minDeltaTimeModifier;
	float maxDeltaTimeModifier;
	float deltaTimeModifier;
	float runningTime;
	int currentFrame;
	float simulationFPS = 30.0;
	bool isSimulationPaused = false;
	bool isSimulationDrawn = true;

	// billboard test
	GLuint texture[1];

	CSPHFluidSimulation fluidSim;
	std::vector<std::array<double, 3>> fluidGradient;
	double minColorDensity = 0.0;
	double maxColorDensity = 100.0;

	// Configuration 
	SimulationConfig mSimulationConfig;
	Graphics mGraphicsConfig;
	SPHFluidSimulation mFluidConfig;
};







