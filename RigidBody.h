#pragma once
#include "Body.h"

// Region R = { x | x = c+r*u[0]+s*u[1]+t*u[2] }, |r|<=e[0], |s|<=e[1], |t|<=e[2]
struct OBB {
	glm::vec3 c; // OBB centre point
	glm::vec3 u[3]; // Local x, y and z axes
	float e[3]; // Positive halfwidth extents of OBB along each axis
};

struct Plane {
	glm::vec3 n; // Plane normal. Point X on the plane satisfies Dot(n, X) = d
	float d;  // d = dot(n, p) for a given point on the plane
};

class RigidBody :
	public Body
{
public:
	RigidBody();
	~RigidBody();

	// set and get methods
	void setAngVel(const glm::vec3 &omega) { m_angVel = omega; }
	void setAngAcc(const glm::vec3 &alpha) { m_angAcc = alpha; }
	void setMass(float mass);
	void setInvInertia(const glm::mat3 &invInertia) { m_invInertia = invInertia; }

	glm::vec3 getAngVel() { return m_angVel; }
	glm::vec3 getAngAcc() { return m_angAcc; }
	glm::mat3 getInvInertia() { return m_invInertia; };
	OBB getObb() { return *m_obb; }
	std::vector<int> getCells() { return m_cells; }
	std::vector<Plane> getPlanes() { return m_planes; }

	void scale(glm::vec3 vect);
	void addCell(int cell) { m_cells.push_back(cell); }
	void clearCells() { m_cells.clear(); }
	glm::mat3 updateInvInertia();
	void updateObb();

	void updatePlanes();
	Plane computePlane(glm::vec3 a, glm::vec3 b, glm::vec3 c);

private:
	float m_density;
	glm::mat3 m_invInertia; // Inverse Inertia
	glm::vec3 m_angVel; // angular velocity
	glm::vec3 m_angAcc; // angular acceleration

	OBB *m_obb;
	std::vector<int> m_cells;
	std::vector<Plane> m_planes;
};