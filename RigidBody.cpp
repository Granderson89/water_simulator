// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include "RigidBody.h"
#include <glm/glm.hpp>
#include "glm/ext.hpp"

RigidBody::RigidBody()
{
	setMesh(Mesh::Mesh());
	Body::scale(glm::vec3(1.0f));

	// set dynamic values
	setAcc(glm::vec3(0.0f, 0.0f, 0.0f));
	setVel(glm::vec3(0.0f, 0.0f, 0.0f));

	// physical properties
	Body::setMass(1.0f);
	setCor(1.0f);

	m_invInertia = updateInvInertia();

	m_obb = new OBB();
}


RigidBody::~RigidBody()
{
}

glm::mat3 RigidBody::updateInvInertia()
{
	// Width
	float w = getScale()[0][0] * 2.0f;
	// Height
	float h = getScale()[1][1] * 2.0f;
	// Depth
	float d = getScale()[2][2] * 2.0f;
	// Mass
	float m = getMass();

	// Inertia x
	float Ix = 1.0f / 12.0f * m * (pow(h, 2) + pow(d, 2));
	// Inertia y
	float Iy = 1.0f / 12.0f * m * (pow(w, 2) + pow(d, 2));
	// Inertia z
	float Iz = 1.0f / 12.0f * m * (pow(w, 2) + pow(h, 2));

	// Build inverse inertia tensor
	glm::mat3 I = glm::scale(glm::mat4(1.0f), glm::vec3(Ix, Iy, Iz));
	return glm::inverse(I);
}

void RigidBody::scale(glm::vec3 vect)
{
	Body::scale(vect);
	setInvInertia(updateInvInertia());
}

void RigidBody::setMass(float mass)
{
	Body::setMass(mass);
	setInvInertia(updateInvInertia());
}

void RigidBody::updateObb()
{
	m_obb->c = getPos();
	m_obb->u[0] = glm::normalize(glm::vec3(getMesh().getRotate()[0]));
	m_obb->u[1] = glm::normalize(glm::vec3(getMesh().getRotate()[1]));
	m_obb->u[2] = glm::normalize(glm::vec3(getMesh().getRotate()[2]));
	m_obb->e[0] = getScale()[0][0];
	m_obb->e[1] = getScale()[1][1];
	m_obb->e[2] = getScale()[2][2];
}

void RigidBody::sleepTest()
{
	if (m_canSleep == false)
		return;

	auto motion = glm::dot(getVel(), getVel()) +
		glm::dot(getAngVel(), getAngVel());
	if (motion < getMotionThreshold())
	{
		setIsAwake(false);
	}
	else
	{
		setIsAwake(true);
	}
}

