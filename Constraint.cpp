#include "Constraint.h"
#define _USE_MATH_DEFINES
#include <math.h>

bool BoundaryConstraint::solveDistanceConstraint(glm::vec3 &corr1) {
	float sumInvMass = 0;
	for (int i = 0; i < m_particles.size(); i++) {
		sumInvMass += m_particles.at(i)->invMass;
	}
	if (sumInvMass == 0.0f) {
		return false;
	}

	glm::vec3 n = m_particles.at(1)->position - m_particles.at(0)->position;
	float d = glm::length(n);

	glm::vec3 corr;
	corr1 = m_k * m_normal * (d - m_restLength) / sumInvMass;
	return true;
}
