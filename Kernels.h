#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "glm/glm.hpp"
#include "Kernels.h"

class Kernels {
public:
	static float poly6(glm::vec3 r, float h);

	static float spiky_grad(glm::vec3 r, float h);

	static glm::vec3 spiky_grad_vec(glm::vec3 r, float h);
	
};