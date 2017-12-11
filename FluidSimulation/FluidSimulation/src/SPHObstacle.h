#pragma once

#include <vector>
#include "glm/glm.hpp"
#include "SPHPartcile.h"

struct SPHObstacle
{
	glm::vec3 position;
	std::vector<SPHParticle*> particles;
	bool isVisible = true;
	int id;
};
