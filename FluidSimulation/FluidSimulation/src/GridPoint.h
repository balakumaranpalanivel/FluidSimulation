#pragma once

#include <vector>
#include <glm\glm.hpp>

struct GridPoint 
{
	glm::vec3 position;
	int id;
	double tx, ty, tz;
	int i, j, k;
	bool isInGridCell = false;
	bool isMarkedForRemoval = false;
};
