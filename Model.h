#pragma once
#include "ParticleData.h"

// The Model contains information on the properties of the system including:
//	- water particles present
//  - the rest density of the water
//  - The densities at each particle
//  - The lambda values at each particle
//  - The position correction for each particle
//  - The cell occupied by each particle

class Model
{
public:
	Model() {};
	~Model() { m_particles.~ParticleData(); };

	// Gets and sets
	ParticleData &getParticles() { return m_particles; }

	const float &getRestDensity() { return m_restDensity; }
	void setRestDensity(float density) { m_restDensity = density; }

	float& getDensity(const unsigned int i) { return m_density.at(i); }
	void setDensity(const unsigned int i, const float &density) { m_density.at(i) = density; }

	float& getLambda(const unsigned int i) { return m_lambda.at(i); }
	void setLambda(const unsigned int i, const float &lambda) { m_lambda.at(i) = lambda; }

	glm::vec3& getDp(const unsigned int i) { return m_dp.at(i); }
	void setDp(const unsigned int i, const glm::vec3 &dp) { m_dp.at(i) = dp; }

	std::vector<unsigned int> & getCell(const unsigned int i) { return m_cell.at(i); }
	void setCell(const unsigned int i, const std::vector<unsigned int>  &cell) { m_cell.at(i) = cell; }

	// Initialise model
	void initModel(const unsigned int numParticles, glm::vec3* particles, Mesh* meshes, const float restDensity);

	// Resize the particle data
	void resizeParticleData(const unsigned int size);

protected:
	// The particles in the model
	ParticleData m_particles;
	// Rest density of the water
	float m_restDensity;
	// Densities for every particle
	std::vector<float> m_density;
	// Lambdas for every particle
	std::vector<float> m_lambda;
	// Position corrections for every particle
	std::vector<glm::vec3> m_dp;
	// Cell occupied by each particle
	std::vector<std::vector<unsigned int>> m_cell;
};