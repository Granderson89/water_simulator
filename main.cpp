#pragma once
// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include <random>
#include <math.h>
#include <map>

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

void initialiseParticles();
void updateSearchGrid(glm::vec3* positions);

using namespace std;

Model model;
// Single water particle properties
const float radius = 1.0f;
const float diam = 2.0f * radius;
// Volume of water
const float width = 3.0f;
const float depth = 2.0f;
const float height = 3.0f;
// Volume of tank
const float tankWidth = (width + 1) * diam;
const float tankDepth = (depth + 1) * diam;
const float tankHeight = 6.0f;
// Search grid
std::map<int, std::vector<int>> grid;
const float gridCellSize = 4.0f;
const glm::vec3 gridMin(-0.5f * tankWidth, 0, -0.5f * tankDepth);;

// Time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;
// Timestep
const GLfloat dt = 0.001f;

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
			//particles.getVel(p) = particles.getVel(p) + dt * glm::vec3(0.0f, -9.8f, 0.0f);
			// Predict new position
			particles.getProj(p) = particles.getPos(p) + dt * particles.getVel(p);
		}

		// For all particles do:
		//	Find neighbouring particles based on proj
		// End for
		// Neighbour search
		// Update the search grid
		updateSearchGrid(particles.getAllProj().data());

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

		// For all particles do:
		//	Update velocity (v = (proj - pos) / dt)
		//	Apply velocity confinement and XSPH viscosity
		//	Update position (pos = proj)
		// End for
		// Commit the position change
		for (unsigned int i = 0; i < particles.getSize(); i++)
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
		for (unsigned int i = 0; i < particles.getSize(); i++)
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
				// Mesh
				waterMeshes[index] = Mesh::Mesh();
				waterMeshes[index].scale(glm::vec3(.1f, .1f, .1f));
				waterMeshes[index].rotate((GLfloat)M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
				waterMeshes[index].setShader(particleShader);
			}
	// Initialise the model
	model.initModel(waterParticles.size(), waterParticles.data(), waterMeshes.data());
	std::cout << "Model initialised" << std::endl;
}

/*
*  BROAD PHASE GRID UPDATE
*/
// Update the grid for each particle. Record which cells they occupy
void updateSearchGrid(glm::vec3* positions)
{
	// Reset the grid
	grid.clear();
	for (int i = 0; i < width * depth * height; i++)
	{
		// Get the particle's position
		glm::vec3 position = positions[i];
		// Check which cell it is in, if the particle is not already recorded
		// as being present in the cell then add it
		int col = floor((position.x - gridMin.x) / gridCellSize);
		int row = floor((position.y - gridMin.y) / gridCellSize);
		int cell = floor((position.z - gridMin.z) / gridCellSize);
		std::string key_string = std::to_string(col) + std::to_string(row) + std::to_string(cell);
		int key = std::stoi(key_string);
		if (std::find(grid[key].begin(), grid[key].end(), i) == grid[key].end())
		{
			grid[key].push_back(i);
		}
	}
}