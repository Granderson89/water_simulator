#pragma once
// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include <random>
#include <math.h>

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

using namespace std;

Model model;
// Single water particle properties
const float radius = 1.0f;
const float diam = 2.0f * radius;
// Volume of water
const float width = 10.0f;
const float depth = 5.0f;
const float height = 5.0f;
// Volume of tank
const float tankWidth = (width + 1) * diam;
const float tankDepth = (depth + 1) * diam;
const float tankHeight = 6.0f;

// Time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;
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
		ParticleData particles = model.getParticles();
		// For all particles do:
		//	Apply forces to velocity (v = v + dt * f)
		//	predict new position based on velocity (proj = pos + dt * v)
		// End for

		// For all particles do:
		//	Find neighbouring particles based on proj
		// End for

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