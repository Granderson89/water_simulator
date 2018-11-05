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

void initialiseParticles();
void boundaryCollisionDetection(ParticleData &particleData);
void createCollisionConstraint(glm::vec3 particle, glm::vec3 ghost, unsigned int index, glm::vec3 collisionNormal);

using namespace std;

Model model;
// Single water particle properties
const float radius = 1.0f;
const float diam = 2.0f * radius;
// Volume of water
const float width = 10.0f;
const float depth = 10.0f;
const float height = 5.0f;
const float restDensity = 0.01f;
// Volume of tank
const float tankWidth = (width + 1) * diam;
const float tankDepth = (depth + 1) * diam;
const float tankHeight = 6.0f;
// List of collision constraints
std::vector<BoundaryConstraint*> collisionConstraints;

// Solver iterations
const unsigned int solverIterations = 4;

// Timestep
const GLfloat dt = 0.1f;

// Particle shader
Shader particleShader;

int main() {
	// Create application
	Application app = Application::Application();
	app.initRender();
	Application::camera.setCameraPosition(glm::vec3(0.0f, 5.0f, 20.0f));
	particleShader = Shader("resources/shaders/core.vert", "resources/shaders/core_blue.frag");
	
	// Set up particles
	initialiseParticles();
	ParticleData particles = model.getParticles();

	// Create a search grid
	SearchGrid searchGrid;
	searchGrid.initSearchGrid(tankWidth, tankDepth, tankHeight, particles.getSize());
	
	// Game loop
	while (!glfwWindowShouldClose(app.getWindow())) {
		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(dt);
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
		searchGrid.updateSearchGrid(model);
		// Create boundary collision constraints
		boundaryCollisionDetection(particles);

		// Solver loop
		unsigned int iters = 0;
		while (iters < solverIterations)
		{	
			// DENSITY CONSTRAINT
			// Calculate lambda
			for (int p = 0; p < particles.getSize(); p++)
			{
				// Get particle p's neighbours
				std::vector<unsigned int> neighbours = searchGrid.getNeighbours(model, p);
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
				std::vector<unsigned int> neighbours = searchGrid.getNeighbours(model, p);
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

		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw particles
		for (int i = 0; i < particles.getSize(); i++)
		{
			app.draw(particles.getMesh(i));
		}
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
				if (j % 2 == 0)
					waterParticles[index] += glm::vec3(diam / 2.0f, 0.0f, 0.0f);
				if (k % 2 == 0)
					waterParticles[index] += glm::vec3(0.0f, 0.0f, diam / 2.0f);
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
void boundaryCollisionDetection(ParticleData &particles)
{
	collisionConstraints.clear();
	for (int i = 0; i < model.getParticles().getSize(); i++)
	{
		auto particle = particles.getProj(i);
		glm::vec3 ghost;
		if (particle.y < 0.0f)
		{
			// Create 'ghost' particle to act as boundary
			ghost = glm::vec3(particle.x, 0.0f, particle.z);
			glm::vec3 collisionNormal(0.0f, 1.0f, 0.0f);
			createCollisionConstraint(particle, ghost, i, collisionNormal);
		}
		if (particle.x < -tankWidth / 2.0f)
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
	BoundaryConstraint *constraint = new BoundaryConstraint(particles, 1.0f, 0.2f, collisionNormal);
	collisionConstraints.push_back(constraint);
}