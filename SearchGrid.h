#pragma once
#include <map>
#include <vector>
#include "glm/glm.hpp"
#include "Model.h"
#include "ParticleData.h"

// SearchGrid splits space into 3D cells and groups particles who occupy
// the same cell.
// Contains:
//  - Map of cell (row, column, depth) to occupants (list of particle indices)
//  - Grid cell size
//  - Min point of the grid
//  - The number of particles in the model
//  - List of neighbours for each particle

class SearchGrid
{
public:
	SearchGrid() : m_gridCellSize(3.0f) {}
	~SearchGrid() {}

	void initSearchGrid(const float tankWidth, const float tankDepth, const float tankHeight, 
		const float numParticles);
	void updateSearchGrid(Model &model, ParticleData &particles);
	std::vector<unsigned int> getNeighbours(const unsigned int index) { return m_neighbours.at(index); };
	void updateNeighbours(Model &model, ParticleData &particles);
	std::map<std::vector<unsigned int>, std::vector<unsigned int>> getGrid() { return m_grid; }

private:
	// Map of cell (row, column, depth) to occupants (list of particle indices)
	std::map<std::vector<unsigned int>, std::vector<unsigned int>> m_grid;
	// Grid cell size
	const float m_gridCellSize;
	// Min point of the grid
	glm::vec3 m_gridMin;
	// The number of particles in the model
	float m_numParticles;
	// List of neighbours for each particle
	std::vector<std::vector<unsigned int>> m_neighbours;
};