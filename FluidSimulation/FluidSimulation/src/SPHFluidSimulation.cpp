#include "SPHFluidSimulation.h"
#include "Configuration.h"
#include "StopWatch.h"
#include <iostream>

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
	s->prevPosition = pos;
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

void CSPHFluidSimulation::UpdateFluidConstants()
{
	InitSimulationConstants();
}

void CSPHFluidSimulation::UpdateObstacleVelocity(double dt)
{
	SPHObstacle *obs;
	SPHParticle *sp;
	glm::vec3 trans;
	
	for (unsigned int i = 0; i < mObstacles.size(); i++)
	{
		obs = mObstacles[i];
		for (unsigned int j = 0; i < obs->particles.size(); j++)
		{
			sp = obs->particles[j];
			if (sp->position == sp->prevPosition)
			{
				sp->velocity = glm::vec3(0.0, 0.0, 0.0);
			}

			trans = sp->position - sp->prevPosition;
			double dist = glm::length(trans);
			double eps = 0.00000001;
			if (dist > eps)
			{
				float speed = fmin((dist / dt), mMaximumVelocity);
				sp->velocity = (trans / (float)dist)*speed;
			}
			else
			{
				sp->velocity = sp->position;
			}
		}
	}
}

void CSPHFluidSimulation::UpdateGrid()
{
	SPHParticle *sp;
	for (unsigned int i = 0; i < mAllParticles.size(); i++)
	{
		sp = mAllParticles[i];
		mGrid.MovePoint(sp->gridID, sp->position);
	}
	mGrid.Update();
}

double CSPHFluidSimulation::CalculateTimeStep()
{
	double maxvsq = 0.0;	// max velocity squared
	double maxcsq = 0.0;	// max speed of sound squard
	double maxasq = 0.0;	// max acceleration squared

	SPHParticle *sp;
	for (unsigned int i = 0; i < mFluidParticles.size(); i++)
	{
		sp = mFluidParticles[i];
		double vsq = glm::dot(sp->velocity, sp->velocity);
		double asq = glm::dot(sp->acceleration, sp->acceleration);
		double csq = EvaluateSpeedOfSoundSquared(sp);
		
		if (vsq > maxvsq)
		{
			maxvsq = vsq;
		}

		if (csq > maxcsq)
		{
			maxcsq = csq;
		}

		if (asq > maxasq)
		{
			maxasq = asq;
		}
	}

	double maxv = sqrt(maxvsq);
	double maxc = sqrt(maxcsq);
	double maxa = sqrt(maxasq);

	double vStep = mCourantSafetyFactor*mSmoothingRadius / fmax(1.0, maxv);
	double cStep = mCourantSafetyFactor*mSmoothingRadius / maxc;
	double aStep = sqrt(mSmoothingRadius / maxa);
	double tempMin = fmin(vStep, cStep);

	return fmax(mMinTimeStep, fmin(tempMin, aStep));
}

void CSPHFluidSimulation::UpdateNearestNeighbours()
{
	SPHParticle *sp;
	for (unsigned int i = 0; i < mAllParticles.size(); i++)
	{
		sp = mAllParticles[i];
		sp->neighbours.clear();
		std::vector<int> refs = mGrid.GetIDsInRadiusOfPoint(sp->gridID, mSmoothingRadius);
		for (unsigned int j = 0; i < refs.size(); j++)
		{
			sp->neighbours.push_back(mParticlesByGridID[refs[j]]);
		}
	}
}

void CSPHFluidSimulation::UpdateFluidDensityAndPressure()
{
	// once we find particles density, we can find its pressure
	SPHParticle *pi, *pj;
	glm::vec3 r;
	for (unsigned int i = 0; i < mAllParticles.size(); i++)
	{
		pi = mAllParticles[i];
		double density = 0.0;
		for (unsigned int j = 0; j < pi->neighbours.size(); j++)
		{
			pj = pi->neighbours[j];
			r = pi->position - pj->position;
			double distsq = glm::dot(r, r);
			double diff = mSmoothingRadiusSquared - distsq;
			density += pj->mass*mPoly6Coefficient*diff*diff*diff;
		}
		pi->denstiy = fmax(density, mInitialDensity); // less than initial density
		// produces negative pressures
		pi->pressure = mPressureCoefficient*(pi->denstiy - mInitialDensity);
	}
}

void CSPHFluidSimulation::UpdateFluidAcceleration()
{
	SPHParticle *pi, *pj;
	glm::vec3 acc;
	glm::vec3 r;
	glm::vec3 vdiff;
	for (unsigned int i = 0; i < mFluidParticles.size(); i++)
	{
		pi = mFluidParticles[i];
		acc = glm::vec3(0.0, 0.0, 0.0);

		for (unsigned int j = 0; i < pi->neighbours.size(); j++)
		{
			pj = pi->neighbours[j];
			r = pi->position - pj->position;
			double dist = glm::length(r);

			if (dist == 0.0)
			{
				continue;
			}

			float inv = 1 / dist;
			r = inv*r;

			// acceleration due to pressure
			float diff = mSmoothingRadius - dist;
			float spikey = mSpikeyGradCoefficient*diff*diff;
			float massRatio = pj->mass / pi->mass;
			float pterm = (pi->pressure + pj->pressure) / (2 * pi->denstiy*pj->denstiy);
			acc -= (float)(massRatio*pterm*spikey)*r;

			// acceleration due to viscosity
			if (!pj->isObstacle)
			{
				float lap = mViscocityLaplacianCoefficient*diff;
				vdiff = pj->velocity - pi->velocity;
				acc += (float)(mViscosityCoefficient*massRatio*(1 / pj->denstiy)*lap)*vdiff;
			}
		}

		// acceleration due to gravity
		acc += mGravityForce;

		// acceleration due to simulation bounds
		acc += CalculateBoundaryAcceleration(pi);

		// motion damping
		double mag = glm::length(acc);
		if (mIsMotionDampingEnabled)
		{
			glm::vec3 damp = pi->velocity * (float)mMotionDampingCoefficient;
			if (glm::length(damp) > mag)
			{
				acc = glm::vec3(0.0, 0.0, 0.0);
			}
			else
			{
				acc -= damp;
			}
		}

		if (mag > mMaximumAcceleration)
		{
			acc = (acc / (float)mag)*(float)mMaximumAcceleration;
		}

		pi->acceleration = acc;
	}
}

glm::vec3 CSPHFluidSimulation::CalculateBoundaryAcceleration(SPHParticle *sp)
{
	double r = mBoundaryForceRadius;
	double minf = mMinBoundaryForce;
	double maxf = mMaxBoundaryForce;

	glm::vec3 p = sp->position;
	glm::vec3 acceleration = glm::vec3(0.0, 0.0, 0.0);

	// X
	if (p.x < mXmin + r)
	{
		double dist = fmax(0.0, p.x - mXmin);
		double force = Utils::Lerp(maxf, minf, dist / r);
		acceleration += glm::vec3(force / sp->mass, 0.0, 0.0);
	}
	else if (p.x < mXmax + r)
	{
		double dist = fmax(0.0, mXmax - p.x);
		double force = Utils::Lerp(maxf, minf, dist / r);
		acceleration += glm::vec3(-force / sp->mass, 0.0, 0.0);
	}

	// Y
	if (p.y < mYmin + r)
	{
		double dist = fmax(0.0, p.y - mYmin);
		double force = Utils::Lerp(maxf, minf, dist / r);
		acceleration += glm::vec3(0.0, force / sp->mass, 0.0);
	}
	else if (p.y > mYmax - r)
	{
		double dist = fmax(0.0, mYmax - p.y);
		double force = Utils::Lerp(maxf, minf, dist / r);
		acceleration += glm::vec3(0.0, -force / sp->mass, 0.0);
	}

	// Z
	if (p.z < mZmin + r)
	{
		double dist = fmax(0.0, p.z - mZmin);
		double force = Utils::Lerp(maxf, minf, dist / r);
		acceleration += glm::vec3(0.0, 0.0, force / sp->mass);
	}
	else if (p.z > mZmax - r)
	{
		double dist = fmax(0.0, mZmax - p.z);
		double force = Utils::Lerp(maxf, minf, dist / r);
		acceleration += glm::vec3(0.0, 0.0, -force / sp->mass);
	}
	
	return acceleration;
}

void CSPHFluidSimulation::UpdateZSortingDistance()
{
	if (!mIsCameraInitialised || !mIsTextureInitialised)
	{
		return;
	}

	glm::vec3 r;
	glm::vec3 cpos = mCamera->position;
	for (unsigned int i = 0; i < mAllParticles.size(); i++)
	{
		r = cpos - mAllParticles[i]->position;
		mAllParticles[i]->zDistance = glm::dot(r, r);
	}
}

bool CSPHFluidSimulation::IsFluidParticleStuckToBoundary(SPHParticle *sp)
{
	double r = mStuckToBoundaryRadius;
	bool isStuck = false;
	glm::vec3 p = sp->position;

	if (p.x < mXmin + r || p.x > mXmax - r ||
		p.y < mYmin + r || p.y > mYmax - r ||
		p.z < mZmin + r || p.z > mZmax - r)
	{
		isStuck = true;
	}

	return isStuck;
}

void CSPHFluidSimulation::UpdateFluidParticleColorDensity(double dt,
	SPHParticle *sp)
{
	// update color density
	// find color density target
	double targetDensity = sp->denstiy;
	double currentDensity = sp->colorDensity;
	double desired = targetDensity - currentDensity;
	if (fabs(desired) < mColorArrivalRadius)
	{
		double r = fabs(desired) / mColorArrivalRadius;
		double newMag = Utils::Lerp(0, mMaxColorVelocity, r);
		if (desired > 0)
		{
			desired = newMag;
		}
		else
		{
			desired = -newMag;
		}
	}
	else
	{
		if (desired > 0)
		{
			desired = mMaxColorVelocity;
		}
		else
		{
			desired = -mMaxColorVelocity;
		}
	}

	// color densiry acceleration
	double acc = desired - sp->colorVelocity;
	if (fabs(acc) > mMaxColorAcceleration)
	{
		if (acc > 0)
		{
			acc = mMaxColorAcceleration;
		}
		else
		{
			acc = -mMaxColorAcceleration;
		}
	}

	// update color density velocity and position from acceleration
	sp->colorVelocity += acc*dt;
	sp->colorDensity += sp->colorVelocity*dt;
	sp->colorDensity = fmax(0.0, sp->colorDensity);
}

glm::vec3 CSPHFluidSimulation::CalculateFluidParticleColor(SPHParticle *sp)
{
	int idxOffset = 50;
	int maxIdx = mFluidGradient.size() - 1;
	int minIdx = fmin(0 + idxOffset, maxIdx);
	double minD = mMinColorDensity;
	double maxD = mMaxColorDensity;
	double invDiff = 1 / (maxD - minD);
	glm::vec3 targetColor;
	std::array<double, 3> cmin;
	std::array<double, 3> cmax;

	// update color based upon color density
	double d = sp->colorDensity;
	double r = 1 - (d - minD)*invDiff;
	r = fmax(0.0, r);
	r = fmin(1.0, r);

	// find target color
	double lerpVal = Utils::Lerp(minIdx, maxIdx, r);
	int minCIdx = floor(lerpVal);
	int maxCIdx = ceil(lerpVal);
	if (minCIdx == maxCIdx)
	{
		cmin = mFluidGradient[minCIdx];
		targetColor = glm::vec3(cmin[0], cmin[1], cmin[2]);
	}
	else
	{
		cmin = mFluidGradient[minCIdx];
		cmax = mFluidGradient[maxCIdx];
		double f = lerpVal - minCIdx;

		targetColor = glm::vec3(Utils::Lerp(cmin[0], cmax[0], f),
			Utils::Lerp(cmin[1], cmax[1], f),
			Utils::Lerp(cmin[2], cmax[2], f));
	}

	return targetColor;
}

void CSPHFluidSimulation::UpdateFluidParticleAlpha(double dt, SPHParticle *sp)
{
	sp->isStuckInBoundary = IsFluidParticleStuckToBoundary(sp);

	if (sp->isStuckInBoundary && sp->boundaryAlphaValue > 0.0)
	{
		sp->boundaryAlphaValue -= mStuckToBoundaryAlphaVelocity*dt;
		sp->boundaryAlphaValue = fmax(0.0, sp->boundaryAlphaValue);
		sp->alpha = Utils::SmoothStep(sp->boundaryAlphaValue);
	}
	else if (!sp->isStuckInBoundary && sp->boundaryAlphaValue < 1.0)
	{
		sp->boundaryAlphaValue += mStuckToBoundaryAlphaVelocity*dt;
		sp->boundaryAlphaValue = fmin(1.0, sp->boundaryAlphaValue);
		sp->alpha = Utils::SmoothStep(sp->boundaryAlphaValue);
	}
}

void CSPHFluidSimulation::UpdateFluidColor(double dt)
{
	if (dt == 0) { return; }

	SPHParticle *sp;

	// find target color index in gradient. Lerp between color indices
	// Move color to target based upon smoothed value of particle density.
	// Particle density is smoothed by keeping track of acceleration and velocity
	// of smoothed value changes
	for (unsigned int i = 0; i < mFluidParticles.size(); i++)
	{
		sp = mFluidParticles[i];

		UpdateFluidParticleColorDensity(dt, sp);
		UpdateFluidParticleAlpha(dt, sp);
		sp->color = CalculateFluidParticleColor(sp);

		// update fluid alpha
	}
}

bool CompareByZDistance(const SPHParticle *p1, const SPHParticle *p2)
{
	return p1->zDistance > p2->zDistance;
}

void CSPHFluidSimulation::UpdateGraphics(double dt)
{
	UpdateZSortingDistance();
	UpdateFluidColor(dt);
	std::sort(mAllParticles.begin(), mAllParticles.end(), CompareByZDistance);
}

void CSPHFluidSimulation::Update(float dt)
{
	UpdateFluidConstants();
	RemoveSPHParticlesMarkedForRemoval();

	CStopWatch neighbourTimer = CStopWatch();
	CStopWatch simulationTimer = CStopWatch();
	CStopWatch graphicsTimer = CStopWatch();

	int numSteps = 0;
	double timeLeft = dt;
	while (timeLeft > 0.0)
	{
		neighbourTimer.Start();
		UpdateGrid();
		UpdateNearestNeighbours();
		neighbourTimer.Stop();

		simulationTimer.Start();
		UpdateFluidDensityAndPressure();
		UpdateFluidAcceleration();

		// calculate next time step
		double timeStep = CalculateTimeStep();
		timeLeft -= timeStep;
		if (timeLeft < 0.0)
		{
			timeStep = timeStep + timeLeft;
			timeLeft = 0.0;
		}
		numSteps += 1;

		UpdateFluidPosition((double)timeStep);
		UpdateObstacleVelocity((double)timeStep);
		simulationTimer.Stop();
	}
	graphicsTimer.Start();
	UpdateGraphics(dt);
	graphicsTimer.Stop();

	mNeighbourSearchTime = neighbourTimer.GetTime();
	mSimulationTime = simulationTimer.GetTime();
	mGraphicsUpdateTime = graphicsTimer.GetTime();

	double total = mNeighbourSearchTime + mSimulationTime;
	std::cout << "Total: " << total << " PCT Neighbour: " << (100 * (mNeighbourSearchTime / total)) << " Particles: " << mAllParticles.size() << std::endl;
}

void CSPHFluidSimulation::UpdateFluidPosition(double dt)
{
	SPHParticle *p;
	for (unsigned int i = 0; i < mFluidParticles.size(); i++)
	{
		p = mFluidParticles[i];

		// calculate velocity at half timestep interval for leapfrog integeration
		if (p->isHalfTimeStepVelocityInitialized)
		{
			p->velocityAtHalfTimeStep += (float)dt*p->acceleration;
		}
		else
		{
			p->velocityAtHalfTimeStep = p->velocity + (float)(0.5*dt)*p->acceleration;
			p->isHalfTimeStepVelocityInitialized = true;
		}

		// new position calculated with half time step for leap frog intergeration
		p->position += (float)dt*p->velocityAtHalfTimeStep;

		// update sph velocity by advancing half time step velocity by
		// 1/2 interval
		p->velocity = p->velocityAtHalfTimeStep + (float)(0.5*dt)*p->acceleration;
		if (glm::length(p->velocity) > mMaximumVelocity)
		{
			glm::vec3 unit = glm::normalize(p->velocity);
			p->velocity = (float)mMaximumVelocity * unit;
		}

		EnforceFluidParticlePositionBounds(p);
	}
}

void CSPHFluidSimulation::EnforceFluidParticlePositionBounds(SPHParticle *p)
{
	if (!mIsEnforcingFluidParticlePositionBoundsThisTimeStep)
	{
		mIsEnforcingFluidParticlePositionBoundsThisTimeStep = true;
		return;
	}

	double eps = 0.001;
	float d = mBoundaryDampingCoefficient;
	if (p->position.x < mXmin)
	{
		p->position = glm::vec3(mXmin + eps, p->position.y, p->position.z);
		p->velocity = glm::vec3(-d*p->velocity.x, p->velocity.y, p->velocity.z);
	}
	else if (p->position.x > mXmax)
	{
		p->position = glm::vec3(mXmax - eps, p->position.y, p->position.z);
		p->velocity = glm::vec3(-d*p->velocity.x, p->velocity.y, p->velocity.z);
	}

	if (p->position.y < mYmin)
	{
		p->position = glm::vec3(p->position.x, mYmin + eps, p->position.z);
		p->velocity = glm::vec3(p->velocity.x, -d*p->velocity.y, p->velocity.z);
	}
	else if (p->position.y > mYmax)
	{
		p->position = glm::vec3(p->position.x, mYmax - eps, p->position.z);
		p->velocity = glm::vec3(p->velocity.x, -d*p->velocity.y, p->velocity.z);
	}

	if (p->position.z < mZmin)
	{
		p->position = glm::vec3(p->position.x, p->position.y, mZmin + eps);
		p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -d*p->velocity.z);
	}
	else if (p->position.z > mZmax)
	{
		p->position = glm::vec3(p->position.x, p->position.y, mZmax - eps);
		p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -d*p->velocity.z);
	}
}

void CSPHFluidSimulation::DrawBounds()
{
	double w = mXmax - mXmin;
	double h = mYmax - mYmin;
	double d = mZmax - mZmin;
	Utils::DrawWireFrameCube(glm::vec3(mXmin+w/2, mYmin+h/2, mZmin+d/2), 
		w, h, d);
}

void CSPHFluidSimulation::Draw()
{
	CStopWatch drawTimer = CStopWatch();

	SPHParticle *sp;
	float size = 0.5*mSmoothingRadius;
	for (unsigned int i = 0; i < mAllParticles.size(); i++)
	{
		sp = mAllParticles[i];
		glm::vec3 p = sp->position;

		if (sp->isObstacle)
		{
			glColor3f(0.5, 0.5, 0.5);
		}
		else
		{
			if (mIsHiddenBoundaryParticlesEnabled)
			{
				glColor4f(sp->color.x, sp->color.y, sp->color.z, sp->alpha);
			}
			else
			{
				glColor3f(sp->color.x, sp->color.y, sp->color.z);
			}
		}

		if (sp->isVisible)
		{
			Utils::DrawBillboard(mCamera, mTexture, p, size);
		}
	}

	DrawBounds();

	drawTimer.Stop();
	mGraphicsDrawTime = drawTimer.GetTime();
}

void CSPHFluidSimulation::AddFluidParticles(std::vector<glm::vec3> points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		AddFluidParticle(points[i], glm::vec3(0.0, 0.0, 0.0));
	}
}

void CSPHFluidSimulation::AddFluidParticle(glm::vec3 pos, glm::vec3 velocity)
{
	SPHParticle *sp = CreateSPHParticle(pos, velocity);
	sp->isObstacle = false;

	sp->gridID = mGrid.InsertPoint(sp->position);
	std::pair<int, SPHParticle*> pair(sp->gridID, sp);
	mParticlesByGridID.insert(pair);

	mFluidParticles.push_back(sp);
	mAllParticles.push_back(sp);
}

void CSPHFluidSimulation::SetBounds(double xmin, double xmax,
	double ymin, double ymax,
	double zmin, double zmax)
{
	mXmin = xmin; mXmax = xmax;
	mYmin = ymin; mYmax = ymax;
	mZmin = zmin; mZmax = zmax;

	InitBoundaryParticles();

	mIsEnforcingFluidParticlePositionBoundsThisTimeStep = false;
}

void CSPHFluidSimulation::SetDampingConstant(double c)
{
	mMotionDampingCoefficient = c;
}

void CSPHFluidSimulation::SetTexture(GLuint *tex)
{
	mTexture = tex;
	mIsTextureInitialized = true;
}

void CSPHFluidSimulation::SetCamera(CCamera *cam)
{
	mCamera = cam;
	mIsCameraInitialised = true;
}

std::vector<SPHParticle*> CSPHFluidSimulation::GetFluidParticles()
{
	return mFluidParticles;
}

std::vector<SPHParticle*> CSPHFluidSimulation::GetObstacleParticles()
{
	return mObstacleParticles;
}

std::vector<SPHParticle*> CSPHFluidSimulation::GetAllParticles()
{
	return mAllParticles;
}

float CSPHFluidSimulation::GetParticleSize()
{
	return mSmoothingRadius;
}

float CSPHFluidSimulation::GetInitialDensity()
{
	return mInitialDensity;
}

void CSPHFluidSimulation::SetObstaclePosition(int id, glm::vec3 pos)
{
	if (mObstaclesByID.find(id) == mObstaclesByID.end()) {
		return;
	}

	SPHObstacle *o = mObstaclesByID[id];
	TranslateObstacle(id, pos - o->position);
}

void CSPHFluidSimulation::TranslateObstacle(int id, glm::vec3 trans)
{
	if (mObstaclesByID.find(id) == mObstaclesByID.end())
	{
		return;
	}
	SPHObstacle *o = mObstaclesByID[id];

	SPHParticle *sp;
	for (unsigned int i = 0; i < o->particles.size(); i++)
	{
		sp = o->particles[i];
		sp->position += trans;
	}

	o->position += trans;
}

void CSPHFluidSimulation::RotateObstacle(int id, Quaternion q)
{
	if (mObstaclesByID.find(id) == mObstaclesByID.end())
	{
		return;
	}

	SPHObstacle *o = mObstaclesByID[id];
	glm::mat4 rot = q.getRotationMatrix();

	SPHParticle *sp;
	glm::vec3 op = o->position;
	glm::vec4 p;
	for (unsigned int i = 0; i<o->particles.size(); i++)
	{
		sp = o->particles[i];
		p = glm::vec4(
			sp->position.x - op.x,
			sp->position.y - op.y,
			sp->position.z - op.z, 1.0);
		p = rot*p;
		sp->position = glm::vec3(p.x + op.x, p.y + op.y, p.z + op.z);
	}

}
