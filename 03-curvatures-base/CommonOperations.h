#ifndef _COMMON_OPERATIONS_INCLUDE
#define _COMMON_OPERATIONS_INCLUDE

#include "NearestNeighbors.h"
#include <Eigen/Dense>

class CommonOperations {

  public:
    CommonOperations();
    ~CommonOperations();
    void setNN(const std::vector<glm::vec3>* points);
    std::vector<Eigen::Vector3f> getNN(glm::vec3 p_i, std::vector<glm::vec3> points, uint k);
    std::vector<size_t> getNN_indexes(glm::vec3 p_i, uint k);
    std::vector<size_t> getNN_in_radius(glm::vec3 p_i, float radius);
    Eigen::Vector3f computeCentroidOfPoints(std::vector<Eigen::Vector3f> points);
    std::vector<Eigen::Vector3f> adjustPointsByCentroid(Eigen::Vector3f centroid, std::vector<Eigen::Vector3f> points);
    Eigen::Matrix3f covarianceMatrix(std::vector<Eigen::Vector3f> points);
  private:
    NearestNeighbors nn;
};
#endif //_COMMON_OPERATIONS_INCLUDE
