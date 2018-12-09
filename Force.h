#pragma once
#include <glm/glm.hpp>

class Body; // forward declaration to avoid circular dependencies
class RigidBody;

class Force
{
public:
	Force() {}
	~Force() {}

	virtual glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);
};

/*
** GRAVITY CLASS
*/
class Gravity : public Force {

public:
	// constructors
	Gravity() {}
	Gravity(const glm::vec3 &gravity) { m_gravity = gravity; }

	// get and set methods
	glm::vec3 getGravity() const { return m_gravity; }
	void setGravity(glm::vec3 gravity) { m_gravity = gravity; }

	// physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	glm::vec3 m_gravity = glm::vec3(0.0f, -9.8f, 0.0f);
};

/*
** DRAG CLASS
*/
class Drag : public Force {
public:
	Drag() {}
	Drag(Body* b, float density, float drag_coeff, float area) { m_b = b;  
		m_density = density; m_drag_coeff = drag_coeff; m_area = area; 
	}

	// get and set methods
	glm::vec3 getDrag();
	void setDrag(glm::vec3 drag) { m_drag = drag; }

	// physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	glm::vec3 m_drag;
	float m_density; // density of medium
	float m_drag_coeff; // coefficient of drag
	float m_area; // cross sectional area

	Body* m_b; // pointer to the body affected by drag
};

/*
** AERO CLASS
*/
class Aero : public Force {
public:
	Aero() {}
	Aero(Body* b_1, Body* b_2, Body* b_3, glm::vec3* wind,  float density, float drag_coeff) {
		m_b1 = b_1; m_b2 = b_2; m_b3 = b_3; m_wind = wind;
		m_density = density; m_drag_coeff = drag_coeff; }
	
	// get and set methods
	glm::vec3 getAero();

	// physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	glm::vec3 m_aero;
	float m_density; // density of medium
	float m_drag_coeff; // coefficient of drag

	// Pointers to each corner of the triangle
	Body* m_b1;
	Body* m_b2;
	Body* m_b3;
	// Pointer to the wind
	glm::vec3* m_wind;
};

/*
** HOOKE'S LAW
*/
class Hooke : public Force {
public:
	Hooke() {}
	Hooke(Body* b1, Body* b2, float ks, float kd, float rest) {
		m_ks = ks; m_kd = kd; m_rest = rest; m_b1 = b1; m_b2 = b2;
	}

	// get and set methods
	glm::vec3 getHooke();
	
	// physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	float m_ks; // spring stiffness
	float m_kd; // damping coefficient
	float m_rest; // spring rest length

	Body* m_b1; // pointer to the body connected to one extremity of the spring
	Body* m_b2; // pointer to the body connected to the other extremity
};
