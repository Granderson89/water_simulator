#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Mesh.h"

// The ParticleData class contains information on all particles including:
//	- current positions
//	- projected positions
//	- current velocities
//	- inverse masses
//	- meshes for rendering

class ParticleData
{
public:
	ParticleData() : m_pos(),
		m_proj(),
		m_vel(),
		m_invMasses(),
		m_meshes()
	{}

	~ParticleData() {
		m_pos.clear();
		m_proj.clear();
		m_vel.clear();
		m_invMasses.clear();
		m_meshes.clear();
	}

	// Gets and sets
	glm::vec3 &getPos(const unsigned int i) { return m_pos[i]; }
	void setPos(const unsigned int i, const glm::vec3 pos)
	{
		m_pos[i] = pos;
		m_meshes[i].setPos(pos);
	}
	glm::vec3 &getProj(const unsigned int i) { return m_proj[i]; }
	void setProj(const unsigned int i, const glm::vec3 proj)
	{
		m_proj[i] = proj;
	}
	glm::vec3 &getVel(const unsigned int i) { return m_vel[i]; }
	void setVel(const unsigned int i, const glm::vec3 vel)
	{
		m_vel[i] = vel;
	}
	float &getInvMass(const unsigned int i) { return m_invMasses[i]; }
	void setMass(const unsigned int i, const float mass)
	{
		if (mass != 0.0f)
			m_invMasses[i] = 1.0f / mass;
		else
			m_invMasses[i] = 0.0f;
	}
	Mesh &getMesh(const unsigned int i) { return m_meshes[i]; }
	void setMesh(const unsigned int i, const Mesh mesh)
	{
		m_meshes[i] = mesh;
	}
	// Get number of particles
	unsigned int getSize() { return m_pos.size(); }

	// Resize data lists
	void resize(const unsigned int size)
	{
		m_pos.resize(size);
		m_proj.resize(size);
		m_vel.resize(size);
		m_invMasses.resize(size);
		m_meshes.resize(size);
	}

	

private:
	// Positions
	std::vector<glm::vec3> m_pos;
	// Projected positions
	std::vector<glm::vec3> m_proj;
	// Velocities
	std::vector<glm::vec3> m_vel;
	// Inverse masses
	std::vector<float> m_invMasses;
	// Meshes
	std::vector<Mesh> m_meshes;
};