#pragma once
#include <vector>
#include <algorithm>
#include "glm/glm.hpp"
#include "Kernels.h"

// PBFluids takes the place of a density constraint.
// Contains functions to calculate the necessary information to solve
// a density constraint

class PBFluids {
public:
	static bool calculateDensity(
		int particleIndex,
		int numParticles,
		glm::vec3 x[],
		float mass[],
		unsigned int numNeighbours,
		unsigned int neighbours[],
		float density0,
		float &density);
	static bool calculateLambda(
		int particleIndex,
		int numParticles,
		glm::vec3 x[],
		float mass[],
		float density,
		unsigned int numNeighbours,
		unsigned int neighbours[],
		float density0,
		float &lambda);
	static bool solveDensityConstraint(
		int particleIndex,
		int numParticles,
		glm::vec3 x[],
		float mass[],
		unsigned int numNeighbours,
		unsigned int neighbours[],
		float density0,
		float lambda[],
		glm::vec3 &corr);
};