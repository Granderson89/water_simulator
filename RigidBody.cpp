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

void RigidBody::updatePlanes()
{
	m_planes.clear();
	// Get the vertices
	auto vertices = getMesh().getVertices();
	// Get the transform matrix
	auto transform = getMesh().getModel();
	// Transform to world space
	std::vector<glm::vec3> ws_vertices;
	for each (Vertex v in vertices)
	{
		glm::vec3 ws = transform * glm::vec4(v.getCoord(), 1.0f);
		ws_vertices.push_back(ws);
	}
	// Update the planes
	Plane p0 = computePlane(ws_vertices[0], ws_vertices[3], ws_vertices[1]);
	m_planes.push_back(p0);

	Plane p1 = computePlane(ws_vertices[1], ws_vertices[2], ws_vertices[5]);
	m_planes.push_back(p1);

	Plane p2 = computePlane(ws_vertices[5], ws_vertices[6], ws_vertices[4]);
	m_planes.push_back(p2);

	Plane p3 = computePlane(ws_vertices[4], ws_vertices[7], ws_vertices[0]);
	m_planes.push_back(p3);

	Plane p4 = computePlane(ws_vertices[3], ws_vertices[7], ws_vertices[2]);
	m_planes.push_back(p4);

	Plane p5 = computePlane(ws_vertices[4], ws_vertices[0], ws_vertices[5]);
	m_planes.push_back(p5);
}

Plane RigidBody::computePlane(glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
	Plane p;
	p.n = glm::normalize(glm::cross(b - a, c - a));
	p.d = glm::dot(p.n, a);
	return p;
}

