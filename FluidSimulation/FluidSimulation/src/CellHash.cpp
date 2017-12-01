#include "CellHash.h"
#include <iostream>

CCellHash::CCellHash()
{
	mMaxNumHashValues = 10000;
}

inline long CCellHash::ComputeHash(int i, int j, int k)
{
	return (abs(541 * (long)i + 79 * (long)j + 31 * (long)k) % mMaxNumHashValues);
}

void CCellHash::InsertGridCell(CGridCell *cell)
{
	long h = ComputeHash(cell->i, cell->j, cell->k);

	if (mCellMap.find(h) == mCellMap.end())
	{
		std::vector<CGridCell*> newChain;
		std::pair<long, std::vector<CGridCell*>> pair(h, newChain);
		mCellMap.insert(pair);
	}
	
	mCellMap[h].push_back(cell);
}

void CCellHash::RemoveGridCell(CGridCell *cell)
{
	int i = cell->i;
	int j = cell->j;
	int k = cell->k;
	long h = ComputeHash(i, j, k);

	if (mCellMap.find(h) == mCellMap.end())
	{
		std::cout << "Could Not Find Cell" << std::endl;
		return;
	}

	// remove from the chain
	bool isRemoved = false;
	std::vector<CGridCell*> chain = mCellMap[h];
	for (int idx = 0; idx<(int)chain.size(); idx++)
	{
		CGridCell *cell = (mCellMap[h])[idx];
		if (cell->i == i && cell->j == j && cell->k == k)
		{
			mCellMap[h].erase(mCellMap[h].begin() + idx);
			isRemoved = true;
			break;
		}
	}

	if (!isRemoved) {
		std::cout << "Could not find/remove gridcell" << i << j << k;
	}

	// remove chain from map if empty
	if (chain.size() == 0) {
		mCellMap.erase(h);
	}
}

CGridCell* CCellHash::FindGridCell(int i, int j, int k, bool *IsGridCellFound)
{
	long h = ComputeHash(i, j, k);

	CGridCell *returnCell = NULL;
	std::vector<CGridCell*> chain = mCellMap[h];

	for (int idx = 0; idx < (int)chain.size(); idx++)
	{
		returnCell = chain[idx];
		if (returnCell->i == i && returnCell->j == j && returnCell->k == k)
		{
			*IsGridCellFound = true;
			return returnCell;
		}
	}

	*IsGridCellFound = false;
	return returnCell;
}

bool CCellHash::IsGridCellInHash(int i, int j, int k)
{
	// TODO - Use Find Grid Cell function itself
	long h = ComputeHash(i, j, k);

	if (mCellMap.find(h) == mCellMap.end()) {
		return false;
	}

	CGridCell *cell;
	std::vector<CGridCell*> chain = mCellMap[h];
	for (int idx = 0; idx<(int)chain.size(); idx++) {
		cell = chain[idx];
		if (cell->i == i && cell->j == j && cell->k == k) {
			return true;
		}
	}
	return false;
}

void CCellHash::GetGridCells(std::vector<CGridCell*> *cells)
{
	for (std::pair<int, std::vector<CGridCell*>> pair : mCellMap)
	{
		for (int i = 0; i < (int)pair.second.size(); i++)
		{
			cells->push_back(pair.second[i]);
		}
	}
}

CGridCell* CCellHash::GetGridCell(int i, int j, int k) 
{
	long h = ComputeHash(i, j, k);

	CGridCell *cell = NULL;
	std::vector<CGridCell*> chain = mCellMap[h];
	for (int idx = 0; idx<(int)chain.size(); idx++)
	{
		cell = chain[idx];
		if (cell->i == i && cell->j == j && cell->k == k)
		{
			return cell;
		}
	}
	return cell;
}