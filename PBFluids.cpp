#include "PBFluids.h"
#include "Kernels.h"
#include <algorithm>
#include <iostream>

// Calculate the density for particle i
bool PBFluids::calculateDensity(int particleIndex, int numParticles, glm::vec3 x[], float mass[], int numNeighbours, int neighbours[], float density0, float &density)
{
	density = 0.0f;
	float h = 2.0f;
	for (int j = 0; j < numNeighbours; j++) {
		int neighbourIndex = neighbours[j];
		if (neighbourIndex < numParticles) {
			auto y = x[neighbourIndex];
			density += Kernels::poly6(x[particleIndex] - x[neighbourIndex], h);
		}
	}

	return true;
}

// Calculate the lambda multiplier for particle i
bool PBFluids::calculateLambda(int particleIndex, int numParticles, glm::vec3 x[], float mass[], float density, int numNeighbours, int neighbours[], float density0, float &lambda)
{
	float h = 2.0f;

	// Set relaxation parameter
	float epsilon = 1.0e-6;

	// Evaluate the constraint function
	float c = std::max(density / density0 - 1.0f, 0.0f);

	// If constraint is not satisfied...
	if (c != 0.0f) {
		// Calculate gradients
		float sum_grad_c = 0.0f;
		glm::vec3 gradC_i(0.0f);

		for (int j = 0; j < numNeighbours; j++) {
			int neighbourIndex = neighbours[j];
			glm::vec3 gradC_j = -1.0f / density0 * Kernels::spiky_grad_vec(x[particleIndex] - x[neighbourIndex], h);
			sum_grad_c += pow(gradC_j.length(), 2);
			gradC_i -= gradC_j;
		}

		sum_grad_c += pow(gradC_i.length(), 2);

		// Calculate lambda
		lambda = -c / (sum_grad_c + epsilon);
	}
	else
		lambda = 0.0f;

	return true;
}
// Calculate position correction
bool PBFluids::solveDensityConstraint(int particleIndex, int numParticles, glm::vec3 x[], float mass[], int numNeighbours, int neighbours[], float density0, float lambda[], glm::vec3 & corr)
{
	float h = 2.0f;
	corr = glm::vec3(0.0f);
	for (int j = 0; j < numNeighbours; j++) {
		int neighbourIndex = neighbours[j];
		glm::vec3 gradC_j = -mass[neighbourIndex] / density0 * Kernels::spiky_grad_vec(x[particleIndex] - x[neighbourIndex], h);
		corr -= (lambda[particleIndex] + lambda[neighbourIndex]) * gradC_j;
	}

	return true;
}
