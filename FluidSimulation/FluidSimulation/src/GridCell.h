#pragma once
#include <vector>
#include <glm\glm.hpp>
#include "GridPoint.h"

class CGridCell
{
public:
	CGridCell();
	void Reset();
	void Initialize(int i, int j, int k);
	void InsertGridPoint(GridPoint* p);
	void RemoveGridPoint(GridPoint *p);
	bool IsEmpty();
	std::vector<GridPoint*> GetGridPoints();

	std::vector<CGridCell*> mNeighbours;
	int i, j, k;

private:
	std::vector<GridPoint*> mPoints;
};