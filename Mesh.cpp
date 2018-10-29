#include "Mesh.h"
#include <errno.h>

/*
**	MESH 
*/

// default constructor creates a horizontal plane or dimensions 1 x 1 centered on the origin
Mesh::Mesh()
{
	// Create vertices
	Vertex vertices[] = { 
		Vertex(glm::vec3(-0.5f ,0.0f ,-0.5f)),
		Vertex(glm::vec3(0.5f, 0.0f, -0.5f)),
		Vertex(glm::vec3(-0.5f ,0.0f ,0.5f)),
		Vertex(glm::vec3(0.5f, 0.0f, -0.5f)),
		Vertex(glm::vec3(-0.5f ,0.0f ,0.5f)),
		Vertex(glm::vec3(0.5f ,0.0f ,0.5f))
	};

	//create mesh
	initMesh(vertices, sizeof(vertices) / sizeof(vertices[0]));

	// initialise tansform matrices (identify)
	initTransform();
}

Mesh::~Mesh()
{
}

/* 
** INIT METHODS 
*/

// initialise transform matrices to identity
void Mesh::initTransform() {
	m_translate = glm::mat4(1.0f);
	m_rotate = glm::mat4(1.0f);
	m_scale = glm::mat4(1.0f);
}

// create mesh from vertices
void Mesh::initMesh(Vertex* vertices, unsigned int numVertices) {
	m_numIndices = numVertices;

	glGenVertexArrays(1, &m_vertexArrayObject);
	glBindVertexArray(m_vertexArrayObject);
	glGenBuffers(NUM_BUFFERS, m_vertexArrayBuffers);
	glBindBuffer(GL_ARRAY_BUFFER, m_vertexArrayBuffers[POSITION_VB]);
	glBufferData(GL_ARRAY_BUFFER, numVertices * sizeof(vertices[0]), vertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glBindVertexArray(0);
}

/*
** TRANSFORMATION METHODS
*/

// translate
void Mesh::translate(const glm::vec3 &vect) {
	m_translate = glm::translate(m_translate, vect);
}

// rotate
void Mesh::rotate(const float &angle, const glm::vec3 &vect) {
	m_rotate = glm::rotate(m_rotate, angle, vect);
}

// scale
void Mesh::scale(const glm::vec3 &vect) {
	m_scale = glm::scale(m_scale, vect);
}


