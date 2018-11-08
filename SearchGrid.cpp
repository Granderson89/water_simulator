#include "SearchGrid.h"
#include <chrono>

using namespace std::chrono;

void SearchGrid::initSearchGrid(const float tankWidth, const float tankDepth, const float tankHeight, const float numParticles)
{
	m_gridMin = glm::vec3(-0.5f * tankWidth, 0, -0.5f * tankDepth);
	m_numParticles = numParticles;
	m_neighbours.resize(m_numParticles);
}

void SearchGrid::updateSearchGrid(Model &model, ParticleData &particles)
{
	// Reset the grid
	m_grid.clear();
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
		std::vector<unsigned int> key;
		key.push_back(col);
		key.push_back(row);
		key.push_back(cell);
		{
			if (std::find(m_grid[key].begin(), m_grid[key].end(), i) == m_grid[key].end())
			{
				// Record in the grid map
				m_grid[key].push_back(i);
				// Record occupied cell for the particle
				model.getCell(i) = key;
			}
		}
	}
	updateNeighbours(model, particles);
}

void SearchGrid::updateNeighbours(Model & model, ParticleData & particles)
{
	for (int p = 0; p < particles.getSize(); p++)
	{
		// Begin a list of neighbours
		std::vector<unsigned int> neighbours;
		// Get cell key that this particle occupies
		std::vector<unsigned int> cell = model.getCell(p);
		unsigned int cell1 = cell.at(0);
		unsigned int cell2 = cell.at(1);
		unsigned int cell3 = cell.at(2);
		// Loop over neighbouring cells
		for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 2; j++)
			{
				for (int k = -1; k < 2; k++)
				{
					std::vector<unsigned int> neighbourCell;
					neighbourCell.push_back(cell1 + i);
					neighbourCell.push_back(cell2 + j);
					neighbourCell.push_back(cell3 + k);
					if (m_grid.find(neighbourCell) != m_grid.end())
					{
						// Get occupants in this neighbouring cell
						std::vector<unsigned int> occupants = m_grid[neighbourCell];
						for (int n = 0; n < occupants.size(); n++)
						{
							if (occupants.at(n) != p)
							{
								float distance = glm::length(particles.getProj(p) - particles.getProj(occupants.at(n)));
								if (distance < 5.0f)
									neighbours.push_back(occupants.at(n));
							}
						}
					}
				}

			}
		}
		m_neighbours.at(p) = neighbours;
	}
}
