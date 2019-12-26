#include "SimpleDistance.h"

// TODO
/* Initialize everything to be able to compute the implicit distance of [Hoppe92] 
   at arbitrary points that are close enough to the point cloud.
 */

void SimpleDistance::init(const PointCloud *pointCloud, float samplingRadius)
{
  co = new CommonOperations();
  co->setNN(&pointCloud->getPoints());
  radius = samplingRadius;
  points = pointCloud;
}


/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the implicit function of [Hoppe92]
   at point P.
 */

bool SimpleDistance::operator()(const glm::vec3 &P, float &value) const
{
  std::vector<size_t> neighbors;
  neighbors = co->getNN_indexes(P, 1);
  size_t index = neighbors[0];

  glm::vec3 p_i = points->getPoints()[index];
  glm::vec3 n_i = points->getNormals()[index];

  float possible = glm::dot(P - p_i, n_i);

  glm::vec3 z = p_i - (possible) * n_i;

  if (glm::distance(z, P) < radius) {
    value = possible;
    return true;
  }
	return false;
}






