#include "GridCell.h"
#include <iostream>

CGridCell::CGridCell()
{
	i = j = k = 0;
}

void CGridCell::Initialize(int ii, int jj, int kk)
{
	i = ii;
	j = jj;
	k = kk;
}

void CGridCell::InsertGridPoint(GridPoint* gp)
{
	gp->i = i;
	gp->j = j;
	gp->k = k;
	gp->isInGridCell = true;
	mPoints.push_back(gp);
}

void CGridCell::RemoveGridPoint(GridPoint* gp)
{
	for (int i = 0; i < (int)mPoints.size(); i++)
	{
		if (mPoints[i]->id == gp->id)
		{
			gp->isInGridCell = false;
			mPoints.erase(mPoints.begin() + i);
			return;
		}
	}

	std::cout << "Cannot Find the GridPoint " << gp->i << " "
		<< gp->j << " "
		<< gp->k << " "
		<< gp->id << " "
		<< std::endl;
}

std::vector<GridPoint*> CGridCell::GetGridPoints()
{
	return mPoints;
}

bool CGridCell::IsEmpty()
{
	return (mPoints.size() == 0);
}

void CGridCell::Reset()
{
	for (int i = 0; i < mPoints.size(); i++)
	{
		mPoints[i]->isInGridCell = false;
	}
	mPoints.clear();
	i = j = k = 0;
}