#include "SPHFluidSimulation.h"
#include "Configuration.h"

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
	SPHFluidSimulation fluidSim;
	Graphics graphicsConfig;
	SimulationConfig simConfig;

	mSmoothingRadiusSquared = mSmoothingRadius*mSmoothingRadius;
	mRatioOfSpecificHeats = fluidSim.RATIO_OF_SPECIFIC_HEAT;
	mPressureCoefficient = fluidSim.PRESSURE_COEFFICIENT;
	mInitialDensity = fluidSim.INITIAL_DENSITY;
	mViscosityCoefficient = fluidSim.VISCOSITY_COEFFICIENT;
	mParticleMass = fluidSim.PARTICLE_MASS;
	mMaximumVelocity = fluidSim.MAXIMUM_VELOCITY;
	mMaximumAcceleration = fluidSim.MAXIMUM_ACCELERATION;
	mMotionDampingCoefficient = fluidSim.MOTION_DAMPING_COEFFICIENT;
	mBoundaryDampingCoefficient = fluidSim.BOUNDARY_DAMPING_COEFFICIENT;
	mGravityMagnitude = fluidSim.GRAVITY_MAGNITUDE;
	mIsMotionDampingEnabled = fluidSim.IS_MOTION_DAMPING_ENABLED;
	mIsBoundaryParticlesEnabled = fluidSim.IS_BOUNDARY_PARTICLES_ENABLED;
	mIsHiddenBoundaryParticlesEnabled = graphicsConfig.IS_HIDDEN_BOUNDARY_PARTICLES_ENABLED;
	mDisplayConsoleOutput = fluidSim.DISPLAY_SIMULATION_CONSOLE_OUTPUT;

	mMinColorDensity = graphicsConfig.MIN_COLOR_DENSITY;
	mMaxColorDensity = graphicsConfig.MAX_COLOR_DENSITY;
	mMaxColorVelocity = graphicsConfig.MAX_COLOR_VELOCITY;
	mMaxColorAcceleration = graphicsConfig.MAX_COLOR_ACCELERATION;
	mColorArrivalRadius = graphicsConfig.COLOR_ARRIVAL_RADIUS;
	mStuckToBoundaryRadius = graphicsConfig.STUCK_TO_BOUNDARY_RADIUS;
	mStuckToBoundaryAlphaVelocity = graphicsConfig.STUCK_TO_BOUNDARY_ALPHA_VELOCITY;

	mGravityForce = glm::vec3(0.0, -mGravityMagnitude, 0.0);
}

void CSPHFluidSimulation::InitBoundaryParticles()
{
	// CONTINUE
	if (!mIsBoundaryParticlesEnabled)
	{
		return;
	}

	if (mIsBoundaryObstacleInitialized)
	{
		RemoveObstacle(mBoundaryObstacleID);
	}

	std::vector<glm::vec3> obsPoints;
	glm::vec3 xDir = glm::vec3(1.0, 0.0, 0.0);
	glm::vec3 yDir = glm::vec3(0.0, 1.0, 0.0);
	glm::vec3 zDir = glm::vec3(0.0, 0.0, 1.0);

	float xWidth = mXmax - mXmin;
	float yWidth = mYmax - mYmin;
	float zWidth = mZmax - mZmin;
	float pad = mSmoothingRadius;
	int layers = 1;
	bool isStaggered = true;

	// x-z (ymin) plane
	glm::vec3 o = glm::vec3(
		mXmin + 0.5*xWidth, 
		0.0,
		mZmin + 0.5*zWidth);
	std::vector<glm::vec3> points = Utils::CreatePointPanel(xWidth, zWidth, pad, 
		layers, xDir, zDir, isStaggered);
	points = Utils::TranslatePoints(points, o);
	obsPoints = Utils::MergePoints(obsPoints, points);

	// x-z (ymax) plane
	o = glm::vec3(
		mXmin + 0.5*xWidth,
		mYmin + yWidth,
		mZmin + 0.5*zWidth);
	points = Utils::CreatePointPanel(xWidth, zWidth, pad, layers, xDir, zDir, isStaggered);
	points = Utils::TranslatePoints(points, o);
	obsPoints = Utils::MergePoints(obsPoints, points);

	// x-y (zmin) plane
	o = glm::vec3(
		mXmin + 0.5*xWidth,
		mYmin + 0.5*xWidth,
		0.0);
	points = Utils::CreatePointPanel(xWidth, yWidth, pad, layers, xDir, yDir, isStaggered);
	points = Utils::TranslatePoints(points, o);
	obsPoints = Utils::MergePoints(obsPoints, points);

	// x-y (zmax) plane
	o = glm::vec3(
		mXmin + 0.5*xWidth,
		mYmin + 0.5*xWidth,
		mZmin + zWidth
	);
	points = Utils::CreatePointPanel(xWidth, yWidth, pad, layers, xDir, yDir, isStaggered);
	points = Utils::TranslatePoints(points, o);
	obsPoints = Utils::MergePoints(obsPoints, points);

	// y-z (xmin) plane
	o = glm::vec3(
		0.0,
		mYmin + 0.5*xWidth,
		mZmin + 0.5*zWidth);
	points = Utils::CreatePointPanel(yWidth, zWidth, pad, layers, yDir, zDir, isStaggered);
	points = Utils::TranslatePoints(points, o);
	obsPoints = Utils::MergePoints(obsPoints, points);

	// y-z (xmax) plane
	o = glm::vec3(
		mXmin + xWidth,
		mYmin + 0.5*xWidth,
		mZmin + 0.5*zWidth);
	points = Utils::CreatePointPanel(yWidth, zWidth, pad, layers, yDir, zDir, isStaggered);
	points = Utils::TranslatePoints(points, o);
	obsPoints = Utils::MergePoints(obsPoints, points);

	mBoundaryObstacleID = AddObstacleParticles(obsPoints);

	SPHObstacle *obs = mObstaclesByID[mBoundaryObstacleID];
	obs->isVisible = false;

	for (unsigned int i=0; i<obs->particles.size(); i++)
	{
		obs->particles[i]->isVisible = false;
	}

	mIsBoundaryObstacleInitialized = true;

}

void CSPHFluidSimulation::RemoveObstacle(int obstacleID)
{
	if (mObstaclesByID.find(obstacleID) == mObstaclesByID.end())
	{
		return;
	}
	mIsSPHParticleRemoved = true;

	SPHObstacle* o = mObstaclesByID[obstacleID];
	mObstaclesByID.erase(obstacleID);

	SPHParticle* p;
	for (unsigned int i = 0; i < o->particles.size(); i++)
	{
		p = o->particles[i];
		mGrid.RemovePoint(p->gridID);
		mParticlesByGridID.erase(p->gridID);
		p->isMarkedForRemoval = true;
	}

	for (unsigned int i = 0; i < mObstacles.size(); i++)
	{
		SPHObstacle* op = mObstacles[i];
		if (op->id == o->id)
		{
			mObstacles.erase(mObstacles.begin() + i);
			break;
		}
	}

	delete o;
}

int CSPHFluidSimulation::AddObstacleParticles(std::vector<glm::vec3> points)
{
	SPHObstacle *obs = new SPHObstacle();
	obs->id = GetUniqueObstacleID();

	SPHParticle *p;
	for (unsigned int i = 0; i < points.size(); i++)
	{
		p = AddObstacleParticle(points[i]);
		obs->particles.push_back(p);
	}
	mObstacles.push_back(obs);

	std::pair<int, SPHObstacle*> pair(obs->id, obs);
	mObstaclesByID.insert(pair);

	return obs->id;
}

void CSPHFluidSimulation::InitKernelConstants()
{
	double pi = 3.1415926535897;

	mPoly6Coefficient = 315.0 / (64.0*pi*powf(mSmoothingRadius, 9.0));
	mSpikeyGradCoefficient = -45.0 / (pi*powf(mSmoothingRadius, 6.0));
	mViscocityLaplacianCoefficient = 45.0 / (pi*powf(mSmoothingRadius, 6.0f));
}

SPHParticle* CSPHFluidSimulation::CreateSPHParticle(glm::vec3 pos, glm::vec3 velocity)
{
	SPHParticle *s = new SPHParticle();

	s->position = pos;
	s->velocity = velocity;
	s->acceleration = glm::vec3(0.0, 0.0, 0.0);

	// create pressure offset from initial pressure
	// uniform densties will cause uniform pressure of 0, meaning no acceleration
	s->denstiy = mInitialDensity;
	
	// mass of sphere
	s->mass = mParticleMass;

	// initial pressure will be calculated once all particles are in place
	s->pressure = 0.0;

	// graphics
	s->color = glm::vec3(1.0, 1.0, 1.0);
	s->colorDensity = s->denstiy;

	return s;
}

SPHParticle* CSPHFluidSimulation::CreateSPHObstacleParticle(glm::vec3 pos)
{
	SPHParticle *p = CreateSPHParticle(pos, glm::vec3(0.0, 0.0, 0.0));
	return p;
}

int CSPHFluidSimulation::GetUniqueObstacleID()
{
	int id = mCurrentObstacleID;
	mCurrentObstacleID++;
	return id;
}

SPHParticle* CSPHFluidSimulation::AddObstacleParticle(glm::vec3 pos)
{
	SPHParticle *sp = CreateSPHObstacleParticle(pos);
	sp->isObstacle = true;

	sp->gridID = mGrid.InsertPoint(sp->position);
	std::pair<int, SPHParticle*> pair(sp->gridID, sp);
	mParticlesByGridID.insert(pair);

	mObstacleParticles.push_back(sp);
	mAllParticles.push_back(sp);

	return sp;
}

inline double CSPHFluidSimulation::EvaluateSpeedOfSound(SPHParticle *sp)
{
	double sqr = mRatioOfSpecificHeats * (sp->pressure) / sp->denstiy;
	if (sqr < 0)
	{
		sqr = -sqr;
	}
	return sqrt(sqr);
}

inline double CSPHFluidSimulation::EvaluateSpeedOfSoundSquared(SPHParticle *sp)
{
	if (sp->denstiy < 0.00001)
	{
		return 0.0;
	}
	return mRatioOfSpecificHeats * (sp->pressure) / sp->denstiy;
}

void CSPHFluidSimulation::RemoveSPHParticlesMarkedForRemoval()
{
	if (mAllParticles.size() == 0 || !mIsSPHParticleRemoved)
	{
		return;
	}

	SPHParticle *p;
	for (int i = (int)mFluidParticles.size() - 1; i >= 0; i--)
	{
		if (mFluidParticles[i]->isMarkedForRemoval)
		{
			p = mFluidParticles[i];
			mFluidParticles.erase(mFluidParticles.begin() + i);
		}
	}

	for (int i = (int)mObstacleParticles.size() - 1; i >= 0; i--)
	{
		if (mObstacleParticles[i]->isMarkedForRemoval)
		{
			p = mObstacleParticles[i];
			mObstacleParticles.erase(mObstacleParticles.begin() + i);
		}
	}

	for (int i = (int)mAllParticles.size() - 1; i >= 0; i--)
	{
		if (mAllParticles[i]->isMarkedForRemoval)
		{
			p = mAllParticles[i];
			mAllParticles.erase(mAllParticles.begin() + i);
			delete p;
		}
	}

	mIsSPHParticleRemoved = false;
}