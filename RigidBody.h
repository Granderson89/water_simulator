#pragma once
#include "Body.h"

// Region R = { x | x = c+r*u[0]+s*u[1]+t*u[2] }, |r|<=e[0], |s|<=e[1], |t|<=e[2]
struct OBB {
	glm::vec3 c; // OBB centre point
	glm::vec3 u[3]; // Local x, y and z axes
	float e[3]; // Positive halfwidth extents of OBB along each axis
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
	void setIsAwake(bool isAwake) { m_isAwake = isAwake; }
	void setMotionThreshold(float motionThreshold) { m_motionThreshold = motionThreshold; }
	void setCanSleep(bool canSleep) { m_canSleep = canSleep; }

	glm::vec3 getAngVel() { return m_angVel; }
	glm::vec3 getAngAcc() { return m_angAcc; }
	glm::mat3 getInvInertia() { return m_invInertia; };
	OBB getObb() { return *m_obb; }
	std::vector<int> getCells() { return m_cells; }
	bool getIsAwake() { return m_isAwake; }
	float getMotionThreshold() { return m_motionThreshold; }
	std::vector<RigidBody*> getTested() { return m_tested; }

	void scale(glm::vec3 vect);
	void addCell(int cell) { m_cells.push_back(cell); }
	void clearCells() { m_cells.clear(); }
	void addTested(RigidBody* rb) { m_tested.push_back(rb); }
	void clearTested() { m_tested.clear(); }
	glm::mat3 updateInvInertia();
	void updateObb();
	void sleepTest();

private:
	float m_density;
	glm::mat3 m_invInertia; // Inverse Inertia
	glm::vec3 m_angVel; // angular velocity
	glm::vec3 m_angAcc; // angular acceleration

	OBB *m_obb;
	std::vector<int> m_cells;
	std::vector<RigidBody*> m_tested;

	// Sleep variables
	bool m_isAwake;
	bool m_canSleep;
	float m_motionThreshold;
};