#include "SphereFunction.h"


void SphereFunction::init(const glm::vec3 &center, float radius)
{
	C = center;
	R = radius;
}

bool SphereFunction::operator()(const glm::vec3 &P, float &value) const
{
	value = glm::distance(P, C) - R;
	
	return true;
}




