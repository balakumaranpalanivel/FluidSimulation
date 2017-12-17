#pragma once
#include "SpatialGrid.h"
#include "Gradient.h"
#include "SPHObstacle.h"

class CSPHFluidSimulation
{
public:
	CSPHFluidSimulation();
	CSPHFluidSimulation(double smoothingRadius);
	~CSPHFluidSimulation();

	//void Update(float dt);
	//void Draw();
	//void DrawBounds();

	int AddObstacleParticles(std::vector<glm::vec3> points);
	void RemoveObstacle(int id);

private:
	// Initialise
	void InitSimulationConstants();
	void InitKernelConstants();

	SPHParticle* CreateSPHParticle(glm::vec3 pos, glm::vec3 velocity);
	SPHParticle* CreateSPHObstacleParticle(glm::vec3 pos);
	SPHParticle* AddObstacleParticle(glm::vec3 pos);
	int GetUniqueObstacleID();
	int mCurrentObstacleID = 0;

	// Simulation
	void InitBoundaryParticles();
	inline double EvaluateSpeedOfSound(SPHParticle *sp);
	inline double EvaluateSpeedOfSoundSquared(SPHParticle *sp);
	void RemoveSPHParticlesMarkedForRemoval();
	void UpdateFluidConstants();
	void UpdateObstacleVelocity(double dt);
	void UpdateGrid();
	double CalculateTimeStep();
	void UpdateNearestNeighbours();
	void UpdateFluidDensityAndPressure();
	void UpdateFluidAcceleration();
	glm::vec3 CalculateBoundaryAcceleration(SPHParticle *sp);

	// graphics
	void UpdateZSortingDistance();
	void UpdateFluidColor(double dt);
	bool IsFluidParticleStuckToBoundary(SPHParticle *sp);
	void UpdateFluidParticleColorDensity(double dt, SPHParticle *sp);
	glm::vec3 CalculateFluidParticleColor(SPHParticle *sp);
	void UpdateFluidParticleAlpha(double dt, SPHParticle *sp);
	bool mIsCameraInitialised = false;
	bool mIsTextureInitialised = false;
	CCamera *mCamera;

	// kernel constants
	double mPoly6Coefficient;
	double mSpikeyGradCoefficient;
	double mViscocityLaplacianCoefficient;

	// graphics
	std::vector<std::array<double, 3>> mFluidGradient;

	// simulation constants
	double mSmoothingRadius;
	double mSmoothingRadiusSquared;
	glm::vec3 mGravityForce;
	double mCourantSafetyFactor = 1.0;
	double mMinTimeStep = 1.0 / 240.0;
	double mGravityMagnitude;
	double mInitialDensity;
	double mPressureCoefficient;
	double mParticleMass;
	bool mIsMotionDampingEnabled = false;
	bool mDisplayConsoleOutput = false;
	double mMotionDampingCoefficient;
	double mBoundaryDampingCoefficient;
	double mRatioOfSpecificHeats;
	double mViscosityCoefficient;
	double mMaximumVelocity;
	double mMaximumAcceleration;

	// boundary constraints
	double mBoundaryForceRadius = 0.1;
	double mMinBoundaryForce = 0.0;
	double mMaxBoundaryForce = 0.0;
	double mXmin = 0.0;
	double mXmax = 1.0;
	double mYmin = 0.0;
	double mYmax = 1.0;
	double mZmin = 0.0;
	double mZmax = 1.0;
	int mBoundaryObstacleID;
	bool mIsBoundaryObstacleInitialized = false;
	bool mIsHiddenBoundaryParticlesEnabled = true;
	bool mIsBoundaryParticlesEnabled = false;

	// graphics
	double mMaxColorVelocity = 1.0;
	double mMaxColorAcceleration = 1.0;
	double mMinColorDensity = 0.0;
	double mMaxColorDensity = 100.0;
	double mColorArrivalRadius = 0.5;
	double mStuckToBoundaryRadius = 0.01;
	double mStuckToBoundaryAlphaVelocity = 1.0;
	bool mIsTextureInitialized = false;

	CSpatialGrid mGrid;

	glm::vec3 mCameraPosition;
	std::unordered_map<int, SPHObstacle*> mObstaclesByID;
	std::unordered_map<int, SPHParticle*> mParticlesByGridID;
	std::vector<SPHParticle*> mObstacleParticles;
	std::vector<SPHParticle*> mAllParticles;
	std::vector<SPHObstacle*> mObstacles;
	std::vector<SPHParticle*> mFluidParticles;
	bool mIsSPHParticleRemoved = false;

};
