#pragma once
#include "SpatialGrid.h"
#include "Gradient.h"

class CSPHFluidSimulation
{
public:
	CSPHFluidSimulation();
	CSPHFluidSimulation(double smoothingRadius);
	~CSPHFluidSimulation();

	//void Update(float dt);
	//void Draw();
	//void DrawBounds();

private:
	// Initialise
	void InitSimulationConstants();
	void InitKernelConstants();

	// Simulation
	void InitBoundaryParticles();

	// kernel constants
	double mPoly6Coefficient;
	double mSpikeyGradCoefficient;
	double mViscocityLaplacianCoefficient;

	// graphics
	std::vector<std::array<double, 3>> mFluidGradient;

	double mSmoothingRadius;
	CSpatialGrid mGrid;
	glm::vec3 mCameraPosition;
};
