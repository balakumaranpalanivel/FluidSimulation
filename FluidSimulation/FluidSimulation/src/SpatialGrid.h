#pragma once
#include <vector>
#include <unordered_map>
#include <GL\glew.h>
#include <cmath>
#include <time.h>
#include <glm\glm.hpp>
#include "CellHash.h"
#include "Utils.h"
#include "GridCell.h"

class CSpatialGrid
{
public:
	CSpatialGrid();
	CSpatialGrid(double cellSize);
	int InsertPoint(glm::vec3 point);
	void MovePoint(int id, glm::vec3 position);
	void RemovePoint(int id);
	std::vector<glm::vec3> GetObjectsInRadiusOfPoint(int ref, double radius);
	std::vector<int> GetIDsInRadiusOfPoint(int ref, double radius);
	void Update();
	void Draw();

private:
	int GenerateUniqueGridPointID();
	void InitializeFreeCells();
	void InsertGridPointIntoGrid(GridPoint *p);
	void PositionToIJK(glm::vec3 p, int *i, int *j, int *k);
	glm::vec3 IJKToPosition(int i, int j, int k);
	CGridCell* GetNewGridCell(int i, int j, int k);
	void UpdateGridPointCellOffset(GridPoint *gp, int i, int j, int k);
	std::vector<int> FastIDNeighbourSearch(int ref, double r, GridPoint *gp);
	void RemoveGridPointsMarkedForRemoval();

	double mSize;
	int mCurrentGridPointID;
	std::vector<GridPoint*> mPoints;
	std::unordered_map<int, GridPoint*> mGridPointsById;
	std::vector<CGridCell*> mFreeCells;
	int mNumInitialFreeCells;
	CCellHash mCellHashTable;
	bool mIsCellRemoved = false;
};