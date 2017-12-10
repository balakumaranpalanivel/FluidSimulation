#include "SPHFluidSimulation.h"

CSPHFluidSimulation::CSPHFluidSimulation()
{
	mSmoothingRadius = 1.0;
}

CSPHFluidSimulation::~CSPHFluidSimulation()
{

}

CSPHFluidSimulation::CSPHFluidSimulation(double smoothingRadius)
{
	mSmoothingRadius = smoothingRadius;
	mGrid = CSpatialGrid(mSmoothingRadius);

	InitSimulationConstants();
	InitKernelConstants();
	InitBoundaryParticles();

	mCameraPosition = glm::vec3(0.0f, 0.0f, 0.0f);
	mFluidGradient = Gradients::GetSkyblueGradient();
}

void CSPHFluidSimulation::InitSimulationConstants()
{
	// TODO - load the constants
}

void CSPHFluidSimulation::InitBoundaryParticles()
{
	// TODO - 
}

void CSPHFluidSimulation::InitKernelConstants()
{
	double pi = 3.1415926535897;

	mPoly6Coefficient = 315.0 / (64.0*pi*powf(mSmoothingRadius, 9.0));
	mSpikeyGradCoefficient = -45.0 / (pi*powf(mSmoothingRadius, 6.0));
	mViscocityLaplacianCoefficient = 45.0 / (pi*powf(mSmoothingRadius, 6.0f));
}