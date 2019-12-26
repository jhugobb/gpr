#ifndef _RBF_FUNCTION_INCLUDE
#define _RBF_FUNCTION_INCLUDE


#include "ImplicitFunction.h"
#include "PointCloud.h"
#include "CommonOperations.h"


class RBFFunction : public ImplicitFunction
{

public:
	void init(const PointCloud *pointCloud, float standardDeviation, float supportRadius);

	bool operator()(const glm::vec3 &P, float &value) const;
	
private:
	float gaussian(float r) const;
	const PointCloud* points;
	float std_dev;
	float radius;
	CommonOperations* co;
	std::vector<double> c_is;
	std::vector<glm::vec3>* ps;
};


#endif // _RBF_FUNCTION_INCLUDE


