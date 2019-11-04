#ifndef _SPHERE_FUNCTION_INCLUDE
#define _SPHERE_FUNCTION_INCLUDE


#include <glm/glm.hpp>
#include "ImplicitFunction.h"


class SphereFunction : public ImplicitFunction
{

public:
	void init(const glm::vec3 &center, float radius);

	bool operator()(const glm::vec3 &P, float &value) const;
	
private:
	glm::vec3 C;
	float R;

};


#endif // _SPHERE_FUNCTION_INCLUDE



