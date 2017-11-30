#pragma once
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
};	