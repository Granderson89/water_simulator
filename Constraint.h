#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <glm/glm.hpp>
#include <vector>

// A Constraint acts on a selection of particles and contains:
//  - A Particle struct to keep necessary information on the particles involved
//  - A cardinality (number of particles it affects)
//  - The list of Particles that are in the constraint
//  - The stiffness of the constraint
//  - The rest length of the constraint

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
	// Caridnality
	int m_n; 
	// List of particles in the constraint
	std::vector<Particle*> m_particles;
	// Stiffness
	float m_k;
	// Rest length
	float m_restLength;
};

/*
**  BOUNDARY CONSTRAINT
*/
// A Boundary Constraint acts on one particle and contains:
//  - The normal of the boundary
//  - A function to solve the constraint and update the position correction
class BoundaryConstraint : public Constraint {
public:
	// Constructor
	BoundaryConstraint(std::vector<Particle*> particles, float stiffness, float restLength, glm::vec3 normal) {
		m_particles = particles; m_k = stiffness; m_restLength = restLength; m_normal = normal;
	}
	// Solve constraint (distance)
	bool solveDistanceConstraint(glm::vec3 &corr);

	unsigned int getIndex() { return m_particles.at(0)->index; }

protected:
	// Normal of the boundary
	glm::vec3 m_normal;
};