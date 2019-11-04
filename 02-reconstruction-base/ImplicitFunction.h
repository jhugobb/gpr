#ifndef _IMPLICIT_FUNCTION_INCLUDE
#define _IMPLICIT_FUNCTION_INCLUDE


#include <glm/glm.hpp>


class ImplicitFunction
{

public:
	virtual bool operator()(const glm::vec3 &P, float &value) const = 0;

};


#endif // _IMPLICIT_FUNCTION_INCLUDE




