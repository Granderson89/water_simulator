#include "SearchGrid.h"



void SearchGrid::initSearchGrid(const float tankWidth, const float tankDepth, const float tankHeight, const float numParticles)
{
	m_gridMin = glm::vec3(-0.5f * tankWidth, 0, -0.5f * tankDepth);
	m_numParticles = numParticles;
}

void SearchGrid::updateSearchGrid(Model &model)
{
	// Reset the grid
	m_grid.clear();
	// Get the particle data
	ParticleData particles = model.getParticles();
	// Record cell position of each particle
	for (int i = 0; i < particles.getSize(); i++)
	{
		// Get the particle's position
		glm::vec3 position = particles.getProj(i);
		// Check which cell it is in, if the particle is not already recorded
		// as being present in the cell then add it
		int col = floor((position.x - m_gridMin.x) / m_gridCellSize);
		int row = floor((position.y - m_gridMin.y) / m_gridCellSize);
		int cell = floor((position.z - m_gridMin.z) / m_gridCellSize);
		std::string key_string = std::to_string(col) + std::to_string(row) + std::to_string(cell);
		int key = std::stoi(key_string);

		if (std::find(m_grid[key].begin(), m_grid[key].end(), i) == m_grid[key].end())
		{
			// Record in the grid map
			m_grid[key].push_back(i);
			// Record occupied cell for the particle
			model.getCell(i) = key;
		}
	}
}

std::vector<unsigned int> SearchGrid::getNeighbours(Model &model, const unsigned int index)
{
	std::vector<unsigned int> neighbours;
	// Get other occupants in this particle's cell
	std::vector<unsigned int> occupants = m_grid[model.getCell(index)];
	// Add to these to the list of neighbours
	for (int i = 0; i < occupants.size(); i++)
	{
		if (occupants[i] != index)
		{
			neighbours.push_back(occupants[i]);
		}
	}
	return neighbours;
}
