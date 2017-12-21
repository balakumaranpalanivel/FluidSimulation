#pragma once
enum Screen
{
	WIDTH = 1280,
	HEIGHT = 720
};

class SPHFluidSimulation
{
public:
	const float RATIO_OF_SPECIFIC_HEAT	= 1.0;
	const float PRESSURE_COEFFICIENT	= 20.0;
	const float INITIAL_DENSITY = 20.0;
	const float VISCOSITY_COEFFICIENT = 0.018;
	const float PARTICLE_MASS = 1.0;
	const float MAXIMUM_VELOCITY = 75.0;
	const float MAXIMUM_ACCELERATION = 75.0;
	const float MOTION_DAMPING_COEFFICIENT = 0.0;
	const float BOUNDARY_DAMPING_COEFFICIENT = 0.6;
	const float GRAVITY_MAGNITUDE = 9.8;
	const bool IS_MOTION_DAMPING_ENABLED = true;
	const bool IS_BOUNDARY_PARTICLES_ENABLED = false;
	const bool DISPLAY_SIMULATION_CONSOLE_OUTPUT = true;
};

static class Graphics
{
public:
	const float MIN_COLOR_DENSITY = 0.0;
	const float MAX_COLOR_DENSITY = 100;
	const float MAX_COLOR_VELOCITY = 100.0;
	const float MAX_COLOR_ACCELERATION = 100.0;
	const float COLOR_ARRIVAL_RADIUS = 0.1;
	const float STUCK_TO_BOUNDARY_RADIUS = 0.01;
	const float STUCK_TO_BOUNDARY_ALPHA_VELOCITY = 3.0;
	const bool IS_HIDDEN_BOUNDARY_PARTICLES_ENABLED = true;
};

static class SimulationConfig
{
public:
	const float FPS = 120;
	const float SMOOTHING_RADIUS = 0.2;
	const int NUM_PARTICLES = 120000;
	const float INITIAL_DAMPING_CONSTANT = 2.0;
	const float FINAL_DAMPING_CONSTANT = 0.0;
	const bool IS_RENDERING_ENABLED = true;
	const bool IS_SIMULATION_PAUSED = false;
	const bool IS_SIMULATION_DRAWN = true;
};