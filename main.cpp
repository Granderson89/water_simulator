#pragma once
// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include <random>
#include <math.h>
#include <map>
#include <omp.h>

// Std. Includes
#include <string>
#include <time.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "glm/ext.hpp"

// Other Libs
#include "SOIL2/SOIL2.h"

// project includes
#include "Application.h"
#include "Shader.h"
#include "Mesh.h"
#include "Model.h"
#include "PBFluids.h"
#include "Constraint.h"
#include "SearchGrid.h"
#include "MarchingCubes.h"
#include "RigidBody.h"

#define _PARTICLES

void initialiseParticles();
void boundaryCollisionDetection(int i, glm::vec3 particle);
void createCollisionConstraint(glm::vec3 particle, glm::vec3 ghost, unsigned int index, glm::vec3 collisionNormal);
void ClosestPtPointOBB(glm::vec3 p, OBB b, glm::vec3 &q);
void rigidBodyCollisionDetection(int i, glm::vec3 particle, RigidBody rigidBody);
int TestPointPolyhedron(glm::vec3 p, Plane *h, int n);
float DistPointPlane(glm::vec3 q, Plane p);

using namespace std;

Model model;
// Single water particle properties
const float radius = 0.25f;
const float diam = 2.0f * radius;
// Volume of water
const float width = 10.0f;
const float depth = 2.0f;
const float height = 15.0f;
const unsigned int num_particles = width * depth * height;
const float restDensity = 0.5f;
// Volume of tank
const float tankWidth = (width + 1) * diam;
const float tankDepth = (depth + 1) * diam;
const float tankHeight = 11.0f;
// List of collision constraints
std::vector<BoundaryConstraint*> collisionConstraints;


// Force plane
glm::vec3 force(-tankWidth / 2.0f - 1.0f, 0.0f, 0.0f);

// Solver iterations
const unsigned int solverIterations = 10;

// Timestep
const GLfloat dt = 0.2f;

// Particle shader
Shader particleShader;

int main() {
	std::cout << "DAM BREAK SCENARIO" << std::endl;
	std::cout << "Volume of water falls into tank" << std::endl;
	/*std::cout << "Use cursor to push water from the left" << std::endl;
	std::cout << "Position cursor at the far left of the screen before starting" << std::endl;
	*/
	std::cout << "Push any button to continue..." << std::endl;
	//cin.get();
	// Create application
	Application app = Application::Application();
	app.initRender();
	Application::camera.setCameraPosition(glm::vec3(0.0f, 5.0f, 20.0f));
	particleShader = Shader("resources/shaders/core.vert", "resources/shaders/core_green.frag");

	// Set up particles
	initialiseParticles();
	ParticleData particles = model.getParticles();

	// Set up a cubic rigid body
	RigidBody rb = RigidBody::RigidBody();
	Mesh m = Mesh::Mesh(Mesh::CUBE);
	rb.setMesh(m);
	Shader rbShader = Shader("resources/shaders/core.vert", "resources/shaders/core_blue.frag");
	rb.getMesh().setShader(rbShader);
	rb.scale(glm::vec3(0.5f, 1.0f, 5.0f));
	rb.setMass(2.0f);
	rb.setPos(glm::vec3(0.0f, 0.0f, 0.0f));
	rb.updateObb();

	int i = 0;
	for each (Vertex v in rb.getMesh().getVertices())
	{
		glm::vec3 ws = rb.getMesh().getModel() * glm::vec4(v.getCoord(), 1.0f);
		std::cout << "Vertex " << i << " ("
			<< ws.x << ", "
			<< ws.y << ", "
			<< ws.z << ")"
			<< std::endl;
			i++;
	}
	rb.updatePlanes();

	// Create a search grid for the neighbour search
	SearchGrid searchGrid;
	searchGrid.initSearchGrid(tankWidth, tankDepth, tankHeight, particles.getSize());

	MarchingCubes mC;
	mC.setOrigin(glm::vec3(-0.5f * tankWidth + diam - 3.0f, diam - 3.0f, -0.5f * tankDepth + diam - 3.0f));
	mC.initMCGrid(tankWidth, tankDepth, tankHeight);

	// Rigid body variables
	GLfloat t = 0.0f;
	GLfloat accumulator = 0.0f;
	GLfloat currentTime = (GLfloat)glfwGetTime();

	// Game loop
	while (!glfwWindowShouldClose(app.getWindow())) {
		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(dt);
		// Get the mouse position
		double mouse_x;
		double mouse_y;
		glfwGetCursorPos(app.getWindow(), &mouse_x, &mouse_y);
		double xx = 2 * mouse_x / app.SCREEN_WIDTH - 1.0f;
		rb.setPos(0, 5* xx);
		/*
		** POSITION BASED FLUIDS
		*/
		// Apply forces and predict new position
		for (unsigned int p = 0; p < particles.getSize(); p++)
		{
			// Apply gravity
			particles.getVel(p) += dt * glm::vec3(0.0f, -9.8f, 0.0f);
			// Predict new position
			particles.getProj(p) = particles.getPos(p) + dt * particles.getVel(p);
		}
		// Neighbour search
		// Update the search grid
		searchGrid.updateSearchGrid(model, particles);
		mC.updateMCNeighbours(model, particles, searchGrid.getGrid(), searchGrid.getGridCellSize());
		mC.updateScalarValues(&particles.getProj(0));

		rb.updateObb();
		rb.updatePlanes();
		// Clear the collision constraints from the last run
		collisionConstraints.clear();
		// Create collision constraints
		for (unsigned int p = 0; p < particles.getSize(); p++)
		{
			boundaryCollisionDetection(p, particles.getProj(p));
			rigidBodyCollisionDetection(p, particles.getProj(p), rb);
		}
		// Solver loop
		unsigned int iters = 0;
		while (iters < solverIterations)
		{	
			
			// DENSITY CONSTRAINT
			// Calculate lambda
			for (int p = 0; p < particles.getSize(); p++)
			{
				// Get particle p's neighbours
				std::vector<unsigned int> neighbours = searchGrid.getNeighbours(p);
				unsigned int numNeighbours = neighbours.size();

				if (numNeighbours > 0)
				{
					// Calculate density
					PBFluids::calculateDensity(p, particles.getSize(), &particles.getProj(0), &particles.getMass(0),
						numNeighbours, &neighbours.at(0), model.getRestDensity(), model.getDensity(p));
					// Calculate lambda
					PBFluids::calculateLambda(p, particles.getSize(), &particles.getProj(0), &particles.getMass(0),
						model.getDensity(p), numNeighbours, &neighbours.at(0), model.getRestDensity(), model.getLambda(p));
				}
			}
			// Calculate dp
			for (int p = 0; p < particles.getSize(); p++)
			{
				// Get particle p's neighbours
				std::vector<unsigned int> neighbours = searchGrid.getNeighbours(p);
				unsigned int numNeighbours = neighbours.size();
				// Calculate dp
				if (numNeighbours > 0)
				{
					glm::vec3 dp;
					PBFluids::solveDensityConstraint(p, particles.getSize(), &particles.getProj(0), &particles.getMass(0),
						numNeighbours, &neighbours.at(0), model.getRestDensity(), &model.getLambda(0), dp);
					model.getDp(p) = dp;
				}				
			}
			// Update projected positions
			for (int p = 0; p < particles.getSize(); p++)
			{
				particles.getProj(p) = particles.getProj(p) + model.getDp(p);
			}

			// COLLISION CONSTRAINTS
			// Solve collision constraints
			for (int i = 0; i < collisionConstraints.size(); i++)
			{
				glm::vec3 dp;
				collisionConstraints.at(i)->solveDistanceConstraint(dp);
				int p = collisionConstraints.at(i)->getIndex();
				model.getDp(p) = dp;
				particles.getProj(p) = particles.getProj(p) + model.getDp(p);
			}
			iters++;
		}

		// Commit the position change
		for (int i = 0; i < particles.getSize(); i++)
		{
			// Update velocity
			particles.getVel(i) = (particles.getProj(i) - particles.getPos(i)) / dt;
			// Update position
			particles.setPos(i, particles.getProj(i));
		}


		// Set frame time
		GLfloat newTime = (GLfloat)glfwGetTime();
		GLfloat frameTime = newTime - currentTime;
		currentTime = newTime;
		accumulator += frameTime;
		while (accumulator >= dt)
		{
			// Simulate each particle's motion
			// Calculate particle's acceleration from resultant force and particle's mass
			rb.setAcc(rb.applyForces(rb.getPos(), rb.getVel(), t, dt));

			// Calculate new velocity and position
			rb.setVel(rb.getVel() + dt * rb.getAcc());
			rb.setPos(glm::vec3(rb.getPos()) + dt * rb.getVel());

			// Calcuate new rotation
			rb.setAngVel(rb.getAngVel() + dt * rb.getAngAcc());
			glm::vec3 dRot = rb.getAngVel() * dt;
			if (glm::dot(dRot, dRot) > 0)
			{
				rb.rotate(sqrt(glm::dot(dRot, dRot)), dRot);
			}
			// Update accumulator and time
			accumulator -= dt;
			t += dt;
		}

		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
#ifdef _PARTICLES
		// draw particles
		for (int i = 0; i < particles.getSize(); i++)
		{
			app.draw(particles.getMesh(i));
		}
		app.draw(rb.getMesh());
#else
		// MARCHING CUBES
		// Generate a list of triangles
		vector<MarchingCubes::TRIANGLE> triangles;
		for (MarchingCubes::GRIDCELL &cell : mC.getCells())
			mC.polygonise(cell, triangles);
		if (triangles.size() != 0)
		{
			// Convert to a list of vertices
			vector<Vertex> vertices;
			for (MarchingCubes::TRIANGLE t : triangles)
			{
				vertices.push_back(t.p[0]);
				vertices.push_back(t.p[1]);
				vertices.push_back(t.p[2]);
			}
			// Create the surface mesh from the vertices and draw
			Mesh surface = Mesh::Mesh();
			surface.initMesh(&vertices.at(0), vertices.size());
			surface.setShader(particleShader);
			app.draw(surface);
		}
#endif
		app.display();
	}

	app.terminate();
	return EXIT_SUCCESS;
}

// Deals with the initial creation and placement of particles and initialises
// the model
void initialiseParticles()
{
	// Water particle initialisation
	// Origin
	glm::vec3 origin(-0.5f * tankWidth + diam, diam, -0.5f * tankDepth + diam);
	const float y_offset = sqrt(3.0f) * radius;
	// Particle start positions
	std::vector<glm::vec3> waterParticles;
	waterParticles.resize(width * depth * height);
	// Particle meshes for rendering
	std::vector<Mesh> waterMeshes;
	waterMeshes.resize(width * depth * height);
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
			for (int k = 0; k < depth; k++)
			{
				int index = i * height * depth + j * depth + k;
				// Position
				waterParticles[index] = diam * glm::vec3(i, j, k) + origin;
				// Mesh
				waterMeshes[index] = Mesh::Mesh();
				waterMeshes[index].scale(glm::vec3(.1f, .1f, .1f));
				waterMeshes[index].rotate((GLfloat)M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
				waterMeshes[index].setShader(particleShader);
			}
	// Initialise the model
	model.initModel(waterParticles.size(), waterParticles.data(), waterMeshes.data(), restDensity);
	std::cout << "Model initialised" << std::endl;
}

// Detect all collisions with the boundaries and create a constraint as a response
void boundaryCollisionDetection(int i, glm::vec3 particle)
{
	glm::vec3 ghost;
	if (particle.y < 0.0f)
	{
		// Create 'ghost' particle to act as boundary
		ghost = glm::vec3(particle.x, 0.0f, particle.z);
		glm::vec3 collisionNormal(0.0f, 1.0f, 0.0f);
		createCollisionConstraint(particle, ghost, i, collisionNormal);
	}
	if (particle.x < -tankWidth / 2.0f)
	//if (particle.x < 0.0f)
	{
		// Create 'ghost' particle to act as boundary
		glm::vec3 ghost(-tankWidth / 2.0f, particle.y, particle.z);
		glm::vec3 collisionNormal(1.0f, 0.0f, 0.0f);
		createCollisionConstraint(particle, ghost, i, collisionNormal);
	}
	else if (particle.x > tankWidth / 2.0f)
	{
		// Create 'ghost' particle to act as boundary
		glm::vec3 ghost(tankWidth / 2.0f, particle.y, particle.z);
		glm::vec3 collisionNormal(-1.0f, 0.0f, 0.0f);
		createCollisionConstraint(particle, ghost, i, collisionNormal);
	}
	if (particle.z < -tankDepth / 2.0f)
	{
		// Create 'ghost' particle to act as boundary
		glm::vec3 ghost(particle.x, particle.y, -tankDepth / 2.0f);
		glm::vec3 collisionNormal(0.0f, 0.0f, 1.0f);
		createCollisionConstraint(particle, ghost, i, collisionNormal);
	}
	else if (particle.z > tankDepth / 2.0f)
	{
		// Create 'ghost' particle to act as boundary
		glm::vec3 ghost(particle.x, particle.y, tankDepth / 2.0f);
		glm::vec3 collisionNormal(0.0f, 0.0f, -1.0f);
		createCollisionConstraint(particle, ghost, i, collisionNormal);
	}
}

// Create a collision constraint, given a particle, it's boundary ghost particle and the boundary normal
void createCollisionConstraint(glm::vec3 particle, glm::vec3 ghost, unsigned int index, glm::vec3 collisionNormal)
{
	// Set up water particle info for the constraint 
	Constraint::Particle *p1 = new Constraint::Particle();
	p1->position = particle;
	p1->invMass = model.getParticles().getInvMass(index);
	p1->index = index;
	// Create a ghost particle at the boundary
	// Bigger inverse mass -> smaller response
	Constraint::Particle *p2 = new Constraint::Particle();
	p2->position = ghost;
	p2->invMass = 20.0f;
	std::vector<Constraint::Particle*> particles;
	particles.push_back(p1);
	particles.push_back(p2);
	// Create a boundary constraint between water and ghost particle
	BoundaryConstraint *constraint = new BoundaryConstraint(particles, 1.0f, 0.1f, collisionNormal);
	collisionConstraints.push_back(constraint);
}

// Detect all collisions with the rigid body and create a constraint as a response
void rigidBodyCollisionDetection(int i, glm::vec3 particle, RigidBody rigidBody)
{
	if (TestPointPolyhedron(particle, &rigidBody.getPlanes().at(0), 6))
	{
		std::cout << "Collision Detected" << std::endl;
	}
}

float DistPointPlane(glm::vec3 q, Plane p)
{
	// return Dot(q, p.n) - p.d; if plane equation normalized (||p.n||==1)
	return (glm::dot(p.n, q) - p.d) / glm::dot(p.n, p.n);
}

// Test if point p inside polyhedron given as the intersection volume of n halfspaces
int TestPointPolyhedron(glm::vec3 p, Plane *h, int n) {
	for (int i = 0; i < n; i++) {
		// Exit with ‘no containment’ if p ever found outside a halfspace
		if (DistPointPlane(p, h[i]) > 0.0f) return 0;
	}
	// p inside all halfspaces, so p must be inside intersection volume
	return 1;
}

// Given point p, return point q on (or in) OBB b, closest to p
void ClosestPtPointOBB(glm::vec3 p, OBB b, glm::vec3 &q)
{
	glm::vec3 d = p - b.c;
	// Start result at center of box; make steps from there
	q = b.c;
	// For each OBB axis...
	for (int i = 0; i < 3; i++) {
		// ...project d onto that axis to get the distance
		// along the axis of d from the box center
		float dist = glm::dot(d, b.u[i]);
		// If distance farther than the box extents, clamp to the box
		if (dist > b.e[i]) dist = b.e[i];
		if (dist < -b.e[i]) dist = -b.e[i];
		// Step that distance along the axis to get world coordinate
		q += dist * b.u[i];
	}
}
