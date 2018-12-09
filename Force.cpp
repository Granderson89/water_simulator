#include <iostream>
#include <cmath>
#include "Force.h"
#include "Body.h"
#include "glm/ext.hpp"

glm::vec3 Force::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	return glm::vec3(0.0f);
}

/*
** GRAVITY
*/
glm::vec3 Gravity::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	return getGravity() / mass;
}

/*
** DRAG
*/
glm::vec3 Drag::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	return getDrag() / mass;
}

glm::vec3 Drag::getDrag() {
	if (m_b->getVel() == glm::vec3(0.0f))
		return glm::vec3(0.0f);
	// Direction
	glm::vec3 e = -m_b->getVel() / glm::length(-m_b->getVel());
	// Calculate drag
	glm::vec3 force_drag = 0.5f * m_density * pow((glm::length(m_b->getVel())), 2) * m_drag_coeff * m_area * e;
	return force_drag;
}

/*
** AERO
*/
glm::vec3 Aero::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	return getAero() / mass;
}

glm::vec3 Aero::getAero() {
	// Get the velocity
	glm::vec3 vel_surface = (m_b1->getVel() + m_b2->getVel() + m_b3->getVel()) / 3;
	// Take wind into account
	glm::vec3 vel = vel_surface - *m_wind;
	// Get the normal
	glm::vec3 norm = glm::cross((m_b1->getPos() - m_b2->getPos()), (m_b3->getPos() - m_b1->getPos()) / glm::length(glm::cross((m_b1->getPos() - m_b2->getPos()), (m_b3->getPos() - m_b1->getPos()))));
	// Get the area of the triangle
	float tri_area = 0.5f * glm::length(glm::cross((m_b1->getPos() - m_b2->getPos()), (m_b3->getPos() - m_b1->getPos())));
	// Area exposed to air flow
	if (vel == glm::vec3(0.0f))
	{
		return glm::vec3(0.0f);
	}
	float area = tri_area * glm::dot(vel, norm) / glm::length(vel);
	// Calculate aerodynamic force on triangle
	glm::vec3 force_aero = 0.5f * m_density * pow(glm::length(vel), 2) * m_drag_coeff * area * norm;
	// Return force on one particle
	return force_aero / 3;
}
	

/*
** HOOKE'S LAW
*/
glm::vec3 Hooke::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	return getHooke() / mass;
}

glm::vec3 Hooke::getHooke() {
	// Compute distance and direction of displacement
	float distance = glm::length(m_b2->getPos() - m_b1->getPos());
	glm::vec3 e = (m_b2->getPos() - m_b1->getPos()) / distance;
	// Compute 1D velocities
	float v1 = glm::dot(e, m_b1->getVel());
	float v2 = glm::dot(e, m_b2->getVel());
	// Compute 1D forces and map back into 3D
	float force_sd = -m_ks * (m_rest - distance) - m_kd * (v1 - v2);
	return -force_sd *e;
}