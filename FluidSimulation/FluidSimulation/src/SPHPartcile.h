#pragma once
#include <vector>
#include <glm/glm.hpp>

struct SPHParticle
{
	glm::vec3 position;
	glm::vec3 prevPosition;
	glm::vec3 velocity;
	glm::vec3 acceleration;

	double mass;
	double denstiy;
	double pressure;

	std::vector<SPHParticle*> neighbours;
	int gridID;		// for lookup ??
	bool isObstacle;
	bool isMarkedForRemoval = false;
	bool isVisible = true;
	double zDistance = 0.0;

	// graphics
	glm::vec3 color;
	double colorDensity;
	double colorVelocity;
	bool isStuckInBoundary = false;
	double boundaryAlphaValue = 1.0;
	double alpha = 1.0;
};