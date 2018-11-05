#pragma once
#include <glm/glm.hpp>
#include <vector>

class Constraint {
public:
	Constraint() {}
	~Constraint() {}

	struct Particle
	{
		Particle() {}
		unsigned int index;
		glm::vec3 position;
		float invMass;
	};

protected:
	int m_n; // caridnality
	std::vector<Particle*> m_particles; // list of particles in the constraint
	float m_k; // stiffness
	bool m_type;
	float m_restLength;
	glm::vec3 m_normal;
};

/*
**  DISTANCE CONSTRAINT
*/
class BoundaryConstraint : public Constraint {
public:
	// Constructor
	BoundaryConstraint(std::vector<Particle*> particles, float stiffness, float restLength, glm::vec3 normal) {
		m_particles = particles; m_k = stiffness; m_restLength = restLength; m_normal = normal;
	}
	// Solve constraint (distance)
	bool solveDistanceConstraint(glm::vec3 &corr1);

	unsigned int getIndex() { return m_particles.at(0)->index; }
};