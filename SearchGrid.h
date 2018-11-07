#pragma once
#include <map>
#include <vector>
#include "glm/glm.hpp"
#include "Model.h"
#include "ParticleData.h"

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

private:
	std::map<std::vector<unsigned int>, std::vector<unsigned int>> m_grid;
	const float m_gridCellSize;
	glm::vec3 m_gridMin;
	float m_numParticles;
	std::vector<std::vector<unsigned int>> m_neighbours;
};