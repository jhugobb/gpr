#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "RBFFunction.h"

//TODO


/* Initialize everything to be able to compute the implicit distance to the reconstructed
   point cloud at arbitrary points that are close enough to the point cloud. As should be
   obvious by the name of the class, the distance has to be computed using RBFs.
 */

void RBFFunction::init(const PointCloud *pointCloud, float standardDeviation, float supportRadius)
{
  std_dev = standardDeviation;
  radius = supportRadius;
  points = pointCloud;

  co = new CommonOperations();

  // float d = supportRadius/2.0f; 
  Eigen::initParallel();

  std::vector<glm::vec3> p_coords = points->getPoints();
  std::vector<glm::vec3> p_normals = points->getNormals();
  float min_x, min_y,min_z;
  min_x = min_y = min_z = std::numeric_limits<float>::max();
  float max_x, max_y, max_z;
  max_x = max_y = max_z = std::numeric_limits<float>::min();


  // Calculate optimal d as 1% of the diagonal of the bounding box of the model
  for (uint i = 0; i < p_coords.size(); i++) {
    max_x = max(max_x, p_coords[i].x);
    max_y = max(max_y, p_coords[i].y);
    max_z = max(max_z, p_coords[i].z);

    min_x = min(min_x, p_coords[i].x);
    min_y = min(min_y, p_coords[i].y);
    min_z = min(min_z, p_coords[i].z);
  }

  glm::vec3 max_v = glm::vec3(max_x, max_y, max_z);
  glm::vec3 min_v = glm::vec3(min_x, min_y, min_z);
  float d = glm::distance(max_v, min_v) * 0.01f;

  // Create point cloud+ with +d and -d points
   ps = new std::vector<glm::vec3>();
  std::vector<float> vs;
  for (uint i = 0; i < p_coords.size(); i++) {

    ps->push_back(p_coords[i]);
    vs.push_back(0);

    ps->push_back(p_coords[i] + d * p_normals[i]);
    vs.push_back(d);

    ps->push_back(p_coords[i] - d * p_normals[i]);
    vs.push_back(-d);
  
  } 

  co->setNN(ps);
  Eigen::SparseMatrix<double> gaussians(ps->size(), ps->size());
  std::vector<size_t> neighs;
  glm::vec3 p1;
  glm::vec3 p2;
  float g;


  cout << "Starting Matrix creation" << endl;
  // Create Gaussian matrix
  std::vector<Eigen::Triplet<float>> trips;
  for (uint i = 0; i < ps->size(); i++) {
    p1 = ps->at(i);
    neighs = co->getNN_in_radius(p1, radius);
    for (uint j = 0; j < neighs.size(); j++) {
      p2 = ps->at(neighs[j]);
      g = gaussian(glm::distance(p1,p2));

      trips.push_back(Eigen::Triplet<float>(i, neighs[j], g));
    }
  }
  gaussians.setFromTriplets(trips.begin(), trips.end());

  // Add lambda I to make the system more stable
  Eigen::SparseMatrix<double> iden_lamb(ps->size(), ps->size());
  iden_lamb.setIdentity();
  iden_lamb *= 0.1f;

  gaussians = gaussians + iden_lamb;


  cout << "ended creating matrix" << endl;
  Eigen::VectorXd vs_eigen(vs.size());

  for (uint i = 0; i < vs.size(); i++) {
    vs_eigen[i] = vs[i];
  }

  // Solve system 
  Eigen::VectorXd cs;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(gaussians);
  cout << "Starting to solve..." << endl;
  cs = solver.solve(vs_eigen);

  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;
  // Save C values
  for (uint i = 0; i < cs.size(); i++) {
    c_is.push_back(cs[i]);
  }
}

/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the RBF implicit distance at point P.
 */

bool RBFFunction::operator()(const glm::vec3 &P, float &value) const
{
  // Once we have the c_i's, sum all gaussians of points inside support radius time c_i's and return that 
  std::vector<size_t> neighs;
  neighs = co->getNN_in_radius(P, radius);
  float result = 0;
  glm::vec3 p_notconst = P;
  if (neighs.size() != 0) {
    for (uint i = 0; i < neighs.size(); i++) {
      result += gaussian(glm::distance(p_notconst, ps->at(neighs[i]))) * c_is[neighs[i]];
    }
    value = result;
    return true;
  } 

	return false;
}

float RBFFunction::gaussian(float r) const {
  float r_2 = pow(r, 2);
  float bottom = 2 * pow(std_dev, 2);

  return exp(-r_2/bottom);
}




