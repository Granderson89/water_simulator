#pragma once
#include "glm/glm.hpp"
#include <vector>

class PBFluids {
public:
	static bool calculateDensity(
		int particleIndex,
		int numParticles,
		glm::vec3 x[],
		float mass[],
		int numNeighbours,
		int neighbours[],
		float density0,
		float &density);
	static bool calculateLambda(
		int particleIndex,
		int numParticles,
		glm::vec3 x[],
		float mass[],
		float density,
		int numNeighbours,
		int neighbours[],
		float density0,
		float &lambda);
	static bool solveDensityConstraint(
		int particleIndex,
		int numParticles,
		glm::vec3 x[],
		float mass[],
		int numNeighbours,
		int neighbours[],
		float density0,
		float lambda[],
		glm::vec3 &corr);
};