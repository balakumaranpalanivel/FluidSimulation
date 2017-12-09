#include "SpatialGrid.h"

CSpatialGrid::CSpatialGrid()
{

}

CSpatialGrid::CSpatialGrid(double cellSize)
{
	mSize = cellSize;
	mCurrentGridPointID = 0;
	mNumInitialFreeCells = 10000;
	
	InitializeFreeCells();

}

void CSpatialGrid::InitializeFreeCells()
{
	for (int i = 0; i < mNumInitialFreeCells; i++)
	{
		CGridCell* cell = new CGridCell();
		mFreeCells.push_back(cell);
	}
}

int CSpatialGrid::InsertPoint(glm::vec3 p)
{
	GridPoint *point = new GridPoint();
	point->position = p;
	point->id = GenerateUniqueGridPointID();
	point->isMarkedForRemoval = false;
	mPoints.push_back(point);

	std::pair<int, GridPoint*> pair(point->id, point);
	mGridPointsById.insert(pair);

	InsertGridPointIntoGrid(point);
	return point->id;
}

void CSpatialGrid::PositionToIJK(glm::vec3 p, int *i, int *j, int *k)
{
	double inv = 1 / mSize;
	*i = ceil(p.x*inv) - 1;
	*j = ceil(p.y*inv) - 1;
	*k = ceil(p.z*inv) - 1;

	double eps = 0.00000000001;
	if (fabs(fmod(p.x, mSize)) < eps) {
		*i = *i + 1;
	}
	if (fabs(fmod(p.y, mSize)) < eps) {
		*j = *j + 1;
	}
	if (fabs(fmod(p.z, mSize)) < eps) {
		*k = *k + 1;
	}
}

glm::vec3 CSpatialGrid::IJKToPosition(int i, int j, int k)
{
	return glm::vec3(i*mSize, j*mSize, k*mSize);
}

int CSpatialGrid::GenerateUniqueGridPointID()
{
	int id = mCurrentGridPointID;
	mCurrentGridPointID++;
	return id;
}

void CSpatialGrid::InsertGridPointIntoGrid(GridPoint *p)
{
	int i, j, k;
	PositionToIJK(p->position, &i, &j, &k);

	bool isCellInTable = false;
	CGridCell *cell = mCellHashTable.FindGridCell(i, j, k, &isCellInTable);

	if (isCellInTable)
	{
		cell->InsertGridPoint(p);
	}
	else
	{
		cell = GetNewGridCell(i, j, k);
		cell->InsertGridPoint(p);
		mCellHashTable.InsertGridCell(cell);
	}

	UpdateGridPointCellOffset(p, i, j, k);
}

void CSpatialGrid::UpdateGridPointCellOffset(GridPoint *gp, int i, int j, int k)
{
	glm::vec3 cp = IJKToPosition(i, j, k);
	gp->tx = gp->position.x - cp.x;
	gp->ty = gp->position.y - cp.y;
	gp->tz = gp->position.z - cp.z;
}

void CSpatialGrid::MovePoint(int id, glm::vec3 newPosition)
{
	if (mGridPointsById.find(id) == mGridPointsById.end())
	{
		return;
	}

	GridPoint *point = mGridPointsById[id];
	int i = point->i;
	int j = point->j;
	int k = point->k;

	glm::vec3 trans = newPosition - point->position;
	point->tx += trans.x;
	point->ty += trans.y;
	point->tz += trans.z;
	point->position = newPosition;

	// point has moved to new cell
	if (point->tx >= mSize || point->ty >= mSize || point->tz >= mSize ||
		point->tx < 0 || point->ty < 0 || point->tz < 0)
	{
		int nexti, nextj, nextk;
		PositionToIJK(point->position, &nexti, &nextj, &nextk);

		// remove grid point from old cell
		CGridCell* oldCell = mCellHashTable.GetGridCell(i, j, k);
		oldCell->RemoveGridPoint(point);

		// remove the cell itself if its empty
		if (oldCell->IsEmpty())
		{
			mCellHashTable.RemoveGridCell(oldCell);
			oldCell->Reset();
			mFreeCells.push_back(oldCell);
		}

		// Insert into new cell
		bool isCellInTable = false;
		CGridCell *cell = mCellHashTable.FindGridCell(nexti, nextj, nextk, &isCellInTable);
		if (isCellInTable)
		{
			cell->InsertGridPoint(point);
		}
		else
		{
			CGridCell *cell = GetNewGridCell(nexti, nextj, nextk);
			cell->InsertGridPoint(point);
			mCellHashTable.InsertGridCell(cell);
		}

		UpdateGridPointCellOffset(point, nexti, nextj, nextk);
	}
}

void CSpatialGrid::RemovePoint(int id)
{
	if (mGridPointsById.find(id) == mGridPointsById.end())
	{
		return;
	}
	mIsCellRemoved = true;

	GridPoint *point = mGridPointsById[id];
	mGridPointsById.erase(id);

	int i = point->i;
	int j = point->j;
	int k = point->k;

	point->isMarkedForRemoval = true;
	bool isCellInTable = false;
	CGridCell *cell = mCellHashTable.FindGridCell(i, j, k, &isCellInTable);

	if (!isCellInTable)
	{
		return;
	}

	cell->RemoveGridPoint(point);
	if (cell->IsEmpty())
	{
		mCellHashTable.RemoveGridCell(cell);
		cell->Reset();
		mFreeCells.push_back(cell);
	}
}

std::vector<glm::vec3> CSpatialGrid::GetObjectsInRadiusOfPoint(int ref, double r)
{
	std::vector<glm::vec3> objects;

	if (mGridPointsById.find(ref) == mGridPointsById.end())
	{
		return objects;
	}

	GridPoint *p = mGridPointsById[ref];
	double tx = p->tx;
	double ty = p->ty;
	double tz = p->tz;
	int i, j, k;
	PositionToIJK(p->position, &i, &j, &k);
	double inv = 1 / mSize;
	double rsq = r*r;

	int imin = i - fmax(0, ceil((r - tx)*inv));
	int jmin = j - fmax(0, ceil((r - ty)*inv));
	int kmin = k - fmax(0, ceil((r - tz)*inv));
	int imax = i + fmax(0, ceil((r - mSize + tx)*inv));
	int jmax = j + fmax(0, ceil((r - mSize + ty)*inv));
	int kmax = k + fmax(0, ceil((r - mSize + tz)*inv));

	CGridCell *cell;
	GridPoint *gp;
	glm::vec3 v;
	std::vector<GridPoint*> points;
	for (int ii = imin; ii <= imax; ii++)
	{
		for (int jj = jmin; jj <= jmax; jj++)
		{
			for (int kk = kmin; kk <= kmax; kk++)
			{

				bool isInHash = false;
				cell = mCellHashTable.FindGridCell(ii, jj, kk, &isInHash);
				if (isInHash)
				{
					points = cell->GetGridPoints();
					for (int idx = 0; idx<(int)points.size(); idx++)
					{
						gp = points[idx];
						if (gp->id != ref)
						{
							v = p->position - gp->position;
							if (glm::dot(v, v) < rsq)
							{
								objects.push_back(gp->position);
							}
						}
					}
				}
			}
		}
	}
	return objects;
}

std::vector<int> CSpatialGrid::FastIDNeighbourSearch(int ref, double r, GridPoint *p) 
{
	std::vector<int> objects;

	bool isInHash = false;
	CGridCell *cell = mCellHashTable.FindGridCell(p->i, p->j, p->k, &isInHash);
	if (!isInHash) {
		return objects;
	}

	std::vector<GridPoint*> points = cell->GetGridPoints();
	GridPoint *gp;
	glm::vec3 v;
	double rsq = r*r;
	for (unsigned int i = 0; i<points.size(); i++)
	{
		gp = points[i];
		if (gp->id != ref)
		{
			v = p->position - gp->position;
			if (glm::dot(v, v) < rsq)
			{
				objects.push_back(gp->id);
			}
		}
	}

	std::vector<CGridCell*> neighbours = cell->mNeighbours;
	for (unsigned int i = 0; i<neighbours.size(); i++)
	{
		points = neighbours[i]->GetGridPoints();
		for (unsigned int j = 0; j<points.size(); j++)
		{
			gp = points[j];
			v = p->position - gp->position;
			if (glm::dot(v, v) < rsq)
			{
				objects.push_back(gp->id);
			}
		}
	}

	return objects;
}

std::vector<int> CSpatialGrid::GetIDsInRadiusOfPoint(int ref, double r)
{

	if (mGridPointsById.find(ref) == mGridPointsById.end())
	{
		std::vector<int> objects;
		return objects;
	}

	GridPoint *p = mGridPointsById[ref];
	double tx = p->tx;
	double ty = p->ty;
	double tz = p->tz;
	int i, j, k;
	PositionToIJK(p->position, &i, &j, &k);
	double inv = 1 / mSize;
	double rsq = r*r;

	int imin = i - fmax(0, ceil((r - tx)*inv));
	int jmin = j - fmax(0, ceil((r - ty)*inv));
	int kmin = k - fmax(0, ceil((r - tz)*inv));
	int imax = i + fmax(0, ceil((r - mSize + tx)*inv));
	int jmax = j + fmax(0, ceil((r - mSize + ty)*inv));
	int kmax = k + fmax(0, ceil((r - mSize + tz)*inv));

	if (imax - imin <= 3 && imax - imin >= 1)
	{
		return FastIDNeighbourSearch(ref, r, p);
	}

	std::vector<int> objects;
	GridPoint *gp;
	CGridCell *cell;
	glm::vec3 v;
	std::vector<GridPoint*> points;
	for (int ii = imin; ii <= imax; ii++) 
	{
		for (int jj = jmin; jj <= jmax; jj++)
		{
			for (int kk = kmin; kk <= kmax; kk++)
			{

				bool isInHash = false;
				cell = mCellHashTable.FindGridCell(ii, jj, kk, &isInHash);
				if (isInHash)
				{
					points = cell->GetGridPoints();
					for (int idx = 0; idx<(int)points.size(); idx++)
					{
						gp = points[idx];
						if (gp->id != ref)
						{
							v = p->position - gp->position;
							if (glm::dot(v, v) < rsq)
							{
								objects.push_back(gp->id);
							}
						}
					}
				}

			}
		}
	}

	return objects;
}

CGridCell* CSpatialGrid::GetNewGridCell(int i, int j, int k)
{
	if (mFreeCells.size() == 0)
	{
		int n = 200;
		for (int i = 0; i<n; i++)
		{
			mFreeCells.push_back(new CGridCell());
		}
	}

	CGridCell *cell = mFreeCells.back();
	mFreeCells.pop_back();

	cell->Initialize(i, j, k);
	return cell;
}

void CSpatialGrid::RemoveGridPointsMarkedForRemoval() 
{
	if (mPoints.size() == 0 || !mIsCellRemoved)
	{
		return;
	}

	GridPoint *p;
	for (int i = (int)mPoints.size() - 1; i >= 0; i--)
	{
		if (mPoints[i]->isMarkedForRemoval)
		{
			p = mPoints[i];
			mPoints.erase(mPoints.begin() + i);

			delete p;
		}
	}

	mIsCellRemoved = false;
}

void CSpatialGrid::Update()
{
	RemoveGridPointsMarkedForRemoval();

	// update each cell's cell neighbours
	std::vector<CGridCell*> cells;
	mCellHashTable.GetGridCells(&cells);

	CGridCell* cell;
	CGridCell* gc;
	for (unsigned int idx = 0; idx<cells.size(); idx++)
	{
		cell = cells[idx];
		cell->mNeighbours.clear();

		int ii = cell->i;
		int jj = cell->j;
		int kk = cell->k;

		for (int k = kk - 1; k <= kk + 1; k++)
		{
			for (int j = jj - 1; j <= jj + 1; j++)
			{
				for (int i = ii - 1; i <= ii + 1; i++)
				{
					if (!(i == ii && j == jj && k == kk))
					{
						bool isInTable = false;
						gc = mCellHashTable.FindGridCell(i, j, k, &isInTable);
						if (isInTable)
						{
							cell->mNeighbours.push_back(gc);
						}
					}
				}
			}
		}

	}
}

void CSpatialGrid::Draw()
{
	if (mPoints.size() == 0)
	{ 
		return; 
	}

	glColor3f(1.0, 0.4, 0.0);
	glPointSize(6.0);
	glBegin(GL_POINTS);
	for (int i = 0; i<(int)mPoints.size(); i++)
	{
		glm::vec3 p = mPoints[i]->position;
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();

	std::vector<CGridCell*> cells;
	mCellHashTable.GetGridCells(&cells);

	glLineWidth(1.0);
	glColor4f(0.0, 0.0, 1.0, 0.4);
	for (int i = 0; i<(int)cells.size(); i++)
	{
		CGridCell *c = cells[i];
		glm::vec3 pos = IJKToPosition(c->i, c->j, c->k);
		pos = pos + glm::vec3(0.5*mSize, 0.5*mSize, 0.5*mSize);
		Utils::DrawWireFrameCube(pos, mSize);
	}
}
