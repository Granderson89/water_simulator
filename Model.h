#pragma once
#include "ParticleData.h"

// The Model contains information on the properties of the system including:
//	- water particles present

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

	// Initialise model
	void initModel(const unsigned int numParticles, glm::vec3* particles, Mesh* meshes, const float restDensity);

	// Resize the particle data
	void resizeParticleData(const unsigned int size);

protected:
	ParticleData m_particles;
	float m_restDensity;
	std::vector<float> m_density;
	std::vector<float> m_lambda;
	std::vector<glm::vec3> m_dp;
};