#ifndef _SIMPLE_DISTANCE_INCLUDE
#define _SIMPLE_DISTANCE_INCLUDE


#include <glm/glm.hpp>
#include "ImplicitFunction.h"
#include "PointCloud.h"
#include "CommonOperations.h"

class SimpleDistance : public ImplicitFunction
{

public:
	void init(const PointCloud *pointCloud, float samplingRadius);

	// If point is too far away (larger than sampling radius), return false
	bool operator()(const glm::vec3 &P, float &value) const;
	
private:

	float radius;
	CommonOperations* co;
	const PointCloud* points;

};


#endif // _SIMPLE_DISTANCE_INCLUDE



