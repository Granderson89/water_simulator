#include "Model.h"

void Model::initModel(const unsigned int numParticles, glm::vec3* particles, Mesh* meshes, const float restDensity)
{
	resizeParticleData(numParticles);
	setRestDensity(restDensity);
	for (int i = 0; i < numParticles; i++)
	{
		m_particles.setMesh(i, meshes[i]);
		m_particles.setPos(i, particles[i]);
		m_particles.setMass(i, 1.0f);
	}
	
}

void Model::resizeParticleData(const unsigned int size)
{
	m_particles.resize(size);
	m_density.resize(size);
	m_lambda.resize(size);
	m_dp.resize(size);
	m_cell.resize(size);
}
