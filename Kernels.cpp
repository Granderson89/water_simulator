#include "Kernels.h"

float Kernels::poly6(glm::vec3 r, float h) {
	float mag_r = glm::length(r);
	if (mag_r >= 0.0f && mag_r <= h) {
		return 315.0f * pow((pow(h, 2) - pow(mag_r, 2)), 3) / (64.0f * M_PI * pow(h, 9));
	}
	else {
		return 0.0f;
	}
}

float Kernels::spiky_grad(glm::vec3 r, float h) {
	float mag_r = glm::length(r);
	if (mag_r >= 0.0f && mag_r <= h) {
		return -45.0f * pow(h - mag_r, 2) / (M_PI * pow(h, 6));
	}
	else {
		return 0.0f;
	}
}

glm::vec3 Kernels::spiky_grad_vec(glm::vec3 r, float h) {
	float mag_r = glm::length(r);
	if (mag_r >= 0.0f && mag_r <= h) {
		return (float)(-45.0f * pow(h - mag_r, 2) / (M_PI * pow(h, 6))) * glm::normalize(r);
	}
	else {
		return glm::vec3(0.0f);
	}
}
