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

void initialiseParticles();
void updateSearchGrid();
std::vector<int> getNeighbours(const unsigned int index);

using namespace std;

Model model;
// Single water particle properties
const float radius = 1.0f;
const float diam = 2.0f * radius;
// Volume of water
const float width = 10.0f;
const float depth = 5.0f;
const float height = 5.0f;
const float restDensity = 0.01f;
// Volume of tank
const float tankWidth = (width + 1) * diam;
const float tankDepth = (depth + 1) * diam;
const float tankHeight = 6.0f;
// Search grid
std::map<int, std::vector<int>> grid;
const float gridCellSize = 3.0f;
const glm::vec3 gridMin(-0.5f * tankWidth, 0, -0.5f * tankDepth);;

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

		// For all particles do:
		//	Find neighbouring particles based on proj
		// End for
		// Neighbour search
		// Update the search grid
		updateSearchGrid();

		// While iter < solverIterations do:
		//	For all particles do:
		//		Calculate lambda
		//	End for
		//	For all particles do:
		//		Calculate dp
		//		Perform collision detection and response
		//	End for
		//	For all particles do:
		//		Update proj (proj = proj + dp)
		//	End for
		// End while
		// Solver loop
		unsigned int iters = 0;
		while (iters < solverIterations)
		{			
#pragma omp parallel for num_threads(8)
			// Calculate lambda
			for (int p = 0; p < particles.getSize(); p++)
			{
				// Get particle p's neighbours
				std::vector<int> neighbours = getNeighbours(p);
				// Calculate density
				if (neighbours.size() > 0)
				{
					PBFluids::calculateDensity(p, particles.getSize(), &particles.getProj(0), &particles.getMass(0),
						neighbours.size(), &neighbours.at(0), model.getRestDensity(), model.getDensity(p));
					// Calculate lambda
					PBFluids::calculateLambda(p, particles.getSize(), &particles.getProj(0), &particles.getMass(0),
						model.getDensity(p), neighbours.size(), &neighbours.at(0), model.getRestDensity(), model.getLambda(p));
				}

			}
			

#pragma omp parallel for num_threads(8)
			// Calculate dp
			for (int p = 0; p < particles.getSize(); p++)
			{
				// Get particle p's neighbours
				std::vector<int> neighbours = getNeighbours(p);
				// Calculate dp
				if (neighbours.size() > 0)
				{
					glm::vec3 dp;
					PBFluids::solveDensityConstraint(p, particles.getSize(), &particles.getProj(0), &particles.getMass(0),
						neighbours.size(), &neighbours.at(0), model.getRestDensity(), &model.getLambda(0), dp);
					model.getDp(p) = dp;
				}
				// Perform collision detection and response
				glm::vec3 preCollision = particles.getProj(p) + model.getDp(p);
				float penetration;
				if (preCollision.y < 0.0f) {
					penetration = preCollision.y;
					model.getDp(p).y -= 1.5f * penetration;
				}
				if (preCollision.x < -tankWidth / 2) {
					penetration = preCollision.x - (-tankWidth / 2);
					model.getDp(p).x -= 1.5f * penetration;
				}
				if (preCollision.x > tankWidth / 2) {
					penetration = preCollision.x - (tankWidth / 2);
					model.getDp(p).x -= 1.5f * penetration;
				}
				if (preCollision.z < -tankDepth / 2) {
					penetration = preCollision.z - (-tankDepth / 2);
					model.getDp(p).z -= 1.5f * penetration;
				}
				if (preCollision.z > tankDepth / 2) {
					penetration = preCollision.z - (tankDepth / 2);
					model.getDp(p).z -= 1.5f * penetration;
				}
			}
#pragma omp parallel for num_threads(8)
			// Update projected positions
			for (int p = 0; p < particles.getSize(); p++)
			{
				particles.getProj(p) = particles.getProj(p) + model.getDp(p);
			}
			iters++;
		}

		// For all particles do:
		//	Update velocity (v = (proj - pos) / dt)
		//	Apply velocity confinement and XSPH viscosity
		//	Update position (pos = proj)
		// End for
#pragma omp parallel for num_threads(8)
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

/*
*  BROAD PHASE GRID UPDATE
*/
// Update the grid for each particle. Record which cells they occupy
void updateSearchGrid()
{
	// Reset the grid
	grid.clear();
	int numParticles = width * depth * height;
#pragma omp parallel for num_threads(8)
	for (int i = 0; i < numParticles; i++)
	{
		// Get the particle's position
		glm::vec3 position = model.getParticles().getProj(i);
		// Check which cell it is in, if the particle is not already recorded
		// as being present in the cell then add it
		int col = floor((position.x - gridMin.x) / gridCellSize);
		int row = floor((position.y - gridMin.y) / gridCellSize);
		int cell = floor((position.z - gridMin.z) / gridCellSize);
		std::string key_string = std::to_string(col) + std::to_string(row) + std::to_string(cell);
		int key = std::stoi(key_string);
#pragma omp critical 
		{
			if (std::find(grid[key].begin(), grid[key].end(), i) == grid[key].end())
			{
				grid[key].push_back(i);
				model.getCell(i) = key;
			}
		}
	}
}

std::vector<int> getNeighbours(const unsigned int index)
{
	std::vector<int> neighbours;
	std::vector<int> occupants = grid[model.getCell(index)];
	for (int i = 0; i < occupants.size(); i++)
	{
		if (occupants[i] != index)
		{
			neighbours.push_back(occupants[i]);
		}
	}
	return neighbours;
}