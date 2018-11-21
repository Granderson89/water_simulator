#include "MarchingCubes.h"

// Initialise the Marching Cubes grid
void MarchingCubes::initMCGrid()
{
	// Set the cell size
	m_cellSize = 3.0f;
	// Set MC grid dimensions
	unsigned int width = 5;
	unsigned int height = 6;
	unsigned int breadth = 3;
	// Create GRIDCELLs and set positions
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
			for (int k = 0; k < breadth; k++)
			{
				GRIDCELL cell;
				cell.p[0] = glm::vec3(i, j, k) * m_cellSize + m_origin;
				cell.p[1] = cell.p[0] + glm::vec3(m_cellSize, 0.0f, 0.0f);
				cell.p[2] = cell.p[0] + glm::vec3(m_cellSize, 0.0f, m_cellSize);
				cell.p[3] = cell.p[0] + glm::vec3(0.0f, 0.0f, m_cellSize);
				cell.p[4] = cell.p[0] + glm::vec3(0.0f, m_cellSize, 0.0f);
				cell.p[5] = cell.p[0] + glm::vec3(m_cellSize, m_cellSize, 0.0f);
				cell.p[6] = cell.p[0] + glm::vec3(m_cellSize, m_cellSize, m_cellSize);
				cell.p[7] = cell.p[0] + glm::vec3(0.0f, m_cellSize, m_cellSize);
				m_mcGrid.push_back(cell);
			}
	// Initialise the neighbour lists for each vertex
	m_neighbours.resize(width * height * breadth);
	for (int i = 0; i < m_neighbours.size(); i++)
		m_neighbours.at(i).resize(8);
}

// Update the neighbouring particles for every vertex of every GRIDCELL
void MarchingCubes::updateMCNeighbours(Model & model, ParticleData & particles, std::map<std::vector<unsigned int>, std::vector<unsigned int>> &searchGrid, float searchGridCellSize)
{
	// Cell identification
	int cellNumber = 0;
	// Update the neighbours
	for (GRIDCELL cell : m_mcGrid)
	{
		// For every every vertex
		for (int vertex = 0; vertex < 8; vertex++)
		{
			// Begin a list of neighbours
			std::vector<unsigned int> neighbours;
			// Get cell key that this cell occupies within the SearchGrid
			glm::vec3 cell_key = cell.p[0] / searchGridCellSize;
			unsigned int cell1 = cell_key.x;
			unsigned int cell2 = cell_key.y;
			unsigned int cell3 = cell_key.z;
			// Loop over neighbouring cells
			for (int i = -1; i < 2; i++)
				for (int j = -1; j < 2; j++)
					for (int k = -1; k < 2; k++)
					{
						// Get neighbouring cell key
						std::vector<unsigned int> neighbourCell;
						neighbourCell.push_back(cell1 + i);
						neighbourCell.push_back(cell2 + j);
						neighbourCell.push_back(cell3 + k);
						// If there is an entry in the neighbouring cell...
						if (searchGrid.find(neighbourCell) != searchGrid.end())
						{
							// Get occupants in this neighbouring cell
							std::vector<unsigned int> occupants = searchGrid[neighbourCell];
							for (int n = 0; n < occupants.size(); n++)
							{
								// If occupant is close enough, add to the list of neighbours
								float distance = glm::length(cell.p[vertex] - particles.getProj(occupants.at(n)));
								if (distance < 5.0f)
									neighbours.push_back(occupants.at(n));

							}
						}
					}
			// Set the neighbours for this vertex
			m_neighbours.at(cellNumber).at(vertex) = neighbours;
		}
		// Move on to the next cell
		cellNumber++;
	}
}

// Update the scalar field value of the surface for every vertex of every GRIDCELL
void MarchingCubes::updateScalarValues(glm::vec3 x[])
{
	// Cell identification
	int cellNumber = 0;
	// Update the scalar values
	for (GRIDCELL &cell : m_mcGrid)
	{
		// For every vertex
		for (int vertex = 0; vertex < 8; vertex++)
		{
			float scalarFieldValue = 0.0f;
			// Smoothing kernel radius
			float h = 2.0f;
			// Accumulate value for every neighbour of this vertex
			for (int neighbour = 0; neighbour < m_neighbours.at(cellNumber).at(vertex).size(); neighbour++)
			{
				int neighbourIndex = m_neighbours.at(cellNumber).at(vertex)[neighbour];
				scalarFieldValue += Kernels::poly6(cell.p[vertex] - x[neighbourIndex], h);
			}
			// Set the value for this vertex
			cell.val[vertex] = scalarFieldValue;
		}
		// Move on to the next cell
		cellNumber++;
	}
}

// Polygonise a GRIDCELL
// Calculates the triangles needed to represent the surface through a GRIDCELL
// The number of triangles is returned (0 if the cell is entirely above or below the surface)
int MarchingCubes::polygonise(GRIDCELL grid, std::vector<TRIANGLE> &triangles)
{
	// Build the index to look up the edge table to get those vertices
	// that are inside the surface
	int cubeIndex = 0;
	if (grid.val[0] > 0.0f) cubeIndex |= 1;
	if (grid.val[1] > 0.0f) cubeIndex |= 2;
	if (grid.val[2] > 0.0f) cubeIndex |= 4;
	if (grid.val[3] > 0.0f) cubeIndex |= 8;
	if (grid.val[4] > 0.0f) cubeIndex |= 16;
	if (grid.val[5] > 0.0f) cubeIndex |= 32;
	if (grid.val[6] > 0.0f) cubeIndex |= 64;
	if (grid.val[7] > 0.0f) cubeIndex |= 128;

	// If all vertices are in or outside of the surface then no polygonisation
	// is required
	if (edgeTable[cubeIndex] == 0)
		return(0);

	// Find the vertices where the surface intersects the cube
	glm::vec3 vertices[12];
	if (edgeTable[cubeIndex] & 1)
		vertices[0] = vertexInterp(grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
	if (edgeTable[cubeIndex] & 2)
		vertices[1] = vertexInterp(grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
	if (edgeTable[cubeIndex] & 4)
		vertices[2] = vertexInterp(grid.p[2], grid.p[3], grid.val[2], grid.val[3]);
	if (edgeTable[cubeIndex] & 8)
		vertices[3] = vertexInterp(grid.p[3], grid.p[0], grid.val[3], grid.val[0]);
	if (edgeTable[cubeIndex] & 16)
		vertices[4] = vertexInterp(grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
	if (edgeTable[cubeIndex] & 32)
		vertices[5] = vertexInterp(grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
	if (edgeTable[cubeIndex] & 64)
		vertices[6] = vertexInterp(grid.p[6], grid.p[7], grid.val[6], grid.val[7]);
	if (edgeTable[cubeIndex] & 128)
		vertices[7] = vertexInterp(grid.p[7], grid.p[4], grid.val[7], grid.val[4]);
	if (edgeTable[cubeIndex] & 256)
		vertices[8] = vertexInterp(grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
	if (edgeTable[cubeIndex] & 512)
		vertices[9] = vertexInterp(grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
	if (edgeTable[cubeIndex] & 1024)
		vertices[10] = vertexInterp(grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
	if (edgeTable[cubeIndex] & 2048)
		vertices[11] = vertexInterp(grid.p[3], grid.p[7], grid.val[3], grid.val[7]);

	// Create the triangle and add to the list
	int num_triangles = 0;
	for (int i = 0; triTable[cubeIndex][i] != -1; i += 3) {
		TRIANGLE t;
		t.p[0] = vertices[triTable[cubeIndex][i]];
		t.p[1] = vertices[triTable[cubeIndex][i + 1]];
		t.p[2] = vertices[triTable[cubeIndex][i + 2]];
		triangles.push_back(t);
		num_triangles++;
	}
	// Return the number of triangles that were required
	return(num_triangles);
}

// Interpolate the position that the surface cuts an edge between two vertices
// (based on the scalar value at each vertex)
glm::vec3 MarchingCubes::vertexInterp(glm::vec3 p1, glm::vec3 p2, float val1, float val2)
{
	glm::vec3 interpolatedPos = p1 + (-val1 / (val2 - val1)) * (p2 - p1);
	return interpolatedPos;
}