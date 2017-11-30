#pragma once
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
	double mSmoothingRadius;
};
