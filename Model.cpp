#include "Model.h"

void Model::initModel(const unsigned int numParticles, glm::vec3* particles, Mesh* meshes)
{
	resizeParticleData(numParticles);
	for (int i = 0; i < numParticles; i++)
	{
		m_particles.setMesh(i, meshes[i]);
		m_particles.setPos(i, particles[i]);
	}
}

void Model::resizeParticleData(const unsigned int size)
{
	m_particles.resize(size);
}
