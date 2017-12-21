#pragma once

#include <GL\glew.h>
#include <glm\glm.hpp>
#include <vector>
#include <cmath>
#include "Camera.h"
#include "Quaternion.h"

namespace Utils
{
	void DrawBillboard(CCamera *camera, GLuint *tex, glm::vec3 p, float width);
	void DrawGrid();
	void DrawWireFrameCube(glm::vec3 pos, float size);
	void DrawWireFrameCube(glm::vec3 pos, float width, float height, float depth);

	float Lerp(float x1, float x2, float t);
	float SmoothStep(float t);

	std::vector<glm::vec3> CreatePointPanel(float width, float height,
		float spacing, int numLayers,
		glm::vec3 w, glm::vec3 h, bool isStaggered);
	std::vector<glm::vec3> TranslatePoints(std::vector<glm::vec3> points, glm::vec3 trans);

	std::vector<glm::vec3> MergePoints(std::vector<glm::vec3> points1,
		std::vector<glm::vec3> points2);

}