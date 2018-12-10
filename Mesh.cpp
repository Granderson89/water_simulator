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

// construct a mesh from a type
Mesh::Mesh(MeshType type)
{
	Vertex vertices[36];

	switch (type)
	{
	case TRIANGLE:
		// Create triangle
		vertices[0] = Vertex(glm::vec3(-1.0f, -1.0f, 0.0f));
		vertices[1] = Vertex(glm::vec3(0.0f, 1.0f, 0.0f));
		vertices[2] = Vertex(glm::vec3(1.0f, -1.0f, 0.0f));
		break;
	case QUAD:
		// create quad
		vertices[0] = Vertex(glm::vec3(-1.0f, 0.0f, -1.0f));
		vertices[1] = Vertex(glm::vec3(1.0f, 0.0f, -1.0f));
		vertices[2] = Vertex(glm::vec3(-1.0f, 0.0f, 1.0f));
		vertices[3] = Vertex(glm::vec3(1.0f, 0.0f, -1.0f));
		vertices[4] = Vertex(glm::vec3(-1.0f, 0.0f, 1.0f));
		vertices[5] = Vertex(glm::vec3(1.0f, 0.0f, 1.0f));
		break;
	case CUBE:
		// create cube
		vertices[0] = Vertex(glm::vec3(-1.0f, -1.0f, -1.0f));
		vertices[1] = Vertex(glm::vec3(1.0f, -1.0f, -1.0f));
		vertices[2] = Vertex(glm::vec3(1.0f, 1.0f, -1.0f));
		vertices[3] = Vertex(glm::vec3(-1.0f, -1.0f, -1.0f));
		vertices[4] = Vertex(glm::vec3(1.0f, 1.0f, -1.0f));
		vertices[5] = Vertex(glm::vec3(-1.0f, 1.0f, -1.0f));
		vertices[6] = Vertex(glm::vec3(-1.0f, -1.0f, 1.0f));
		vertices[7] = Vertex(glm::vec3(1.0f, -1.0f, 1.0f));
		vertices[8] = Vertex(glm::vec3(1.0f, 1.0f, 1.0f));
		vertices[9] = Vertex(glm::vec3(-1.0f, -1.0f, 1.0f));
		vertices[10] = Vertex(glm::vec3(1.0f, 1.0f, 1.0f));
		vertices[11] = Vertex(glm::vec3(-1.0f, 1.0f, 1.0f));
		vertices[12] = Vertex(glm::vec3(-1.0f, -1.0f, -1.0f));
		vertices[13] = Vertex(glm::vec3(1.0f, -1.0f, -1.0f));
		vertices[14] = Vertex(glm::vec3(1.0f, -1.0f, 1.0f));
		vertices[15] = Vertex(glm::vec3(-1.0f, -1.0f, -1.0f));
		vertices[16] = Vertex(glm::vec3(1.0f, -1.0f, 1.0f));
		vertices[17] = Vertex(glm::vec3(-1.0f, -1.0f, 1.0f));
		vertices[18] = Vertex(glm::vec3(-1.0f, 1.0f, -1.0f));
		vertices[19] = Vertex(glm::vec3(1.0f, 1.0f, -1.0f));
		vertices[20] = Vertex(glm::vec3(1.0f, 1.0f, 1.0f));
		vertices[21] = Vertex(glm::vec3(-1.0f, 1.0f, -1.0f));
		vertices[22] = Vertex(glm::vec3(1.0f, 1.0f, 1.0f));
		vertices[23] = Vertex(glm::vec3(-1.0f, 1.0f, 1.0f));
		vertices[24] = Vertex(glm::vec3(-1.0f, -1.0f, -1.0f));
		vertices[25] = Vertex(glm::vec3(-1.0f, 1.0f, -1.0f));
		vertices[26] = Vertex(glm::vec3(-1.0f, 1.0f, 1.0f));
		vertices[27] = Vertex(glm::vec3(-1.0f, -1.0f, -1.0f));
		vertices[28] = Vertex(glm::vec3(-1.0f, 1.0f, 1.0f));
		vertices[29] = Vertex(glm::vec3(-1.0f, -1.0f, 1.0f));
		vertices[30] = Vertex(glm::vec3(1.0f, -1.0f, -1.0f));
		vertices[31] = Vertex(glm::vec3(1.0f, 1.0f, -1.0f));
		vertices[32] = Vertex(glm::vec3(1.0f, 1.0f, 1.0f));
		vertices[33] = Vertex(glm::vec3(1.0f, -1.0f, -1.0f));
		vertices[34] = Vertex(glm::vec3(1.0f, 1.0f, 1.0f));
		vertices[35] = Vertex(glm::vec3(1.0f, -1.0f, 1.0f));
		break;
	}

	// generate unique vertex vector (no duplicates)
	m_vertices = std::vector<Vertex>(std::begin(vertices), std::end(vertices));
	for (int i = 0; i < m_vertices.size(); i++)
	{
		for (int j = 0; j < m_vertices.size(); j++)
		{
			if (i == j)
				continue;
			if (m_vertices[i].getCoord() == m_vertices[j].getCoord())
			{
				m_vertices[j] = Vertex(glm::vec3(0.0f));
			}
		}
	}
	for (int i = 0; i < m_vertices.size(); i++)
	{
		if (m_vertices[i].getCoord() == glm::vec3(0.0f))
		{
			m_vertices.erase(m_vertices.begin() + i);
			i--;
		}
	}

	// create mesh
	initMesh(vertices, sizeof(vertices) / sizeof(vertices[0]));

	// create model matrix (identity)
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


