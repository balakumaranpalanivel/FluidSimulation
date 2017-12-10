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

