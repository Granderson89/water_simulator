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

	// Initialise model
	void initModel(const unsigned int numParticles, glm::vec3* particles, Mesh* meshes);

	// Resize the particle data
	void resizeParticleData(const unsigned int size);

protected:
	ParticleData m_particles;
};