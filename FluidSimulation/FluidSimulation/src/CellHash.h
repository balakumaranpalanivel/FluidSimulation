#pragma once
#include <unordered_map>
#include <vector>
#include "GridCell.h"

class CCellHash
{
public:
	CCellHash();
	bool IsGridCellInHash(int i, int j, int k);
	void InsertGridCell(CGridCell *cell);
	void RemoveGridCell(CGridCell *cell);
	CGridCell* GetGridCell(int i, int j, int k);
	CGridCell* FindGridCell(int i, int j, int k, bool *IsGridCellFound);
	void GetGridCells(std::vector<CGridCell*> *cells);

private:
	inline long ComputeHash(int i, int j, int k);

	long mMaxNumHashValues;
	std::unordered_map<long, std::vector<CGridCell*>> mCellMap;
};