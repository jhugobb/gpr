#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include "LaplacianSmoothing.h"


void LaplacianSmoothing::setMesh(TriangleMesh *newMesh)
{
	mesh = newMesh;
}

/* This method should apply nIterations iterations of the laplacian vector multiplied by lambda 
   to each of the vertices. */

void LaplacianSmoothing::iterativeLaplacian(int nIterations, float lambda)
{
  std::vector<uint> neighbors;
  std::vector<glm::vec3> new_points;
  for (uint i = 0; i < nIterations; i++) {
    std::vector<glm::vec3> points = mesh->getVertices();
    new_points.clear();
    new_points.resize(points.size());
    for (uint k = 0; k < points.size(); k++) {
      glm::vec3 p_i = points[k];
      neighbors.clear();
      mesh->getNeighbors(k, neighbors);
      glm::vec3 laplacian = glm::vec3(0,0,0);
      for (uint j = 0; j < neighbors.size(); j++) {
        glm::vec3 p_j = points[neighbors[j]];
        laplacian += p_j - p_i;
      }
      laplacian /= neighbors.size();
      new_points[k] = p_i + lambda * laplacian;
    }

    for (uint k = 0; k < points.size(); k++) {
      mesh->getVertices()[k] = new_points[k];
    }
  }
}

/* This method should apply nIterations iterations of the bilaplacian operator using lambda 
   as a scaling factor. */

void LaplacianSmoothing::iterativeBilaplacian(int nIterations, float lambda)
{
  std::vector<uint> neighbors;
  std::vector<glm::vec3> new_points;
  for (uint i = 0; i < nIterations; i++) {
    std::vector<glm::vec3> points = mesh->getVertices();
    new_points.clear();
    new_points.resize(points.size());
    for (uint k = 0; k < points.size(); k++) {
      glm::vec3 p_i = points[k];
      neighbors.clear();
      mesh->getNeighbors(k, neighbors);
      glm::vec3 laplacian = glm::vec3(0,0,0);
      for (uint j = 0; j < neighbors.size(); j++) {
        glm::vec3 p_j = points[neighbors[j]];
        laplacian += p_j - p_i;
      }
      laplacian /= neighbors.size();
      new_points[k] = p_i + lambda * laplacian;
    }

    for (uint k = 0; k < new_points.size(); k++) {
      glm::vec3 p_i = new_points[k];
      neighbors.clear();
      mesh->getNeighbors(k, neighbors);
      glm::vec3 laplacian = glm::vec3(0,0,0);
      for (uint j = 0; j < neighbors.size(); j++) {
        glm::vec3 p_j = new_points[neighbors[j]];
        laplacian += p_j - p_i;
      }
      laplacian /= neighbors.size();
      new_points[k] = p_i - lambda * laplacian;
    }

    for (uint k = 0; k < points.size(); k++) {
      mesh->getVertices()[k] = new_points[k];
    }
  }
}

/* This method should apply nIterations iterations of Taubin's operator using lambda 
   as a scaling factor, and computing the corresponding nu value. */

void LaplacianSmoothing::iterativeLambdaNu(int nIterations, float lambda)
{
  float nu = 1.f/(1.f/10.f - 1.f/lambda);
  std::vector<uint> neighbors;
  std::vector<glm::vec3> new_points;
  for (uint i = 0; i < nIterations; i++) {
    std::vector<glm::vec3> points = mesh->getVertices();
    new_points.clear();
    new_points.resize(points.size());
    for (uint k = 0; k < points.size(); k++) {
      glm::vec3 p_i = points[k];
      neighbors.clear();
      mesh->getNeighbors(k, neighbors);
      glm::vec3 laplacian = glm::vec3(0,0,0);
      for (uint j = 0; j < neighbors.size(); j++) {
        glm::vec3 p_j = points[neighbors[j]];
        laplacian += p_j - p_i;
      }
      laplacian /= neighbors.size();
      new_points[k] = p_i + lambda * laplacian;
    }

    for (uint k = 0; k < new_points.size(); k++) {
      glm::vec3 p_i = new_points[k];
      neighbors.clear();
      mesh->getNeighbors(k, neighbors);
      glm::vec3 laplacian = glm::vec3(0,0,0);
      for (uint j = 0; j < neighbors.size(); j++) {
        glm::vec3 p_j = new_points[neighbors[j]];
        laplacian += p_j - p_i;
      }
      laplacian /= neighbors.size();
      new_points[k] = p_i + nu * laplacian;
    }

    for (uint k = 0; k < points.size(); k++) {
      mesh->getVertices()[k] = new_points[k];
    }
  }
}

/* This method should compute new vertices positions by making the laplacian zero, while
   maintaing the vertices marked as constraints fixed. */

void LaplacianSmoothing::globalLaplacian(const vector<bool> &constraints)
{
  // For every point 1/neighs on the position of the neighs, 0 otherwise (sparse), -1 on the index of the row we are creating
  // b will be 0 if no point is constrained, 1/neighs * x/y/z of the constrained points
  std::vector<glm::vec3> points = mesh->getVertices();
  
  Eigen::SparseMatrix<float> A(3*points.size(), 3*points.size());
  std::vector<uint> neighs;
  glm::vec3 p1;
  glm::vec3 p2;
  Eigen::VectorXf b(3*points.size());

  cout << "Starting Matrix creation" << endl;
  // Create Gaussian matrix
  std::vector<Eigen::Triplet<float>> trips;
  for (uint i = 0; i < points.size(); i++) {
    p1 = points[i];
    mesh->getNeighbors(i, neighs);
    b[3*i + 0] = 0;
    b[3*i + 1] = 0;
    b[3*i + 2] = 0;
    if (constraints[i]) {
      trips.push_back(Eigen::Triplet<float>(3*i + 0, 3*i + 0, -1.0f));
      trips.push_back(Eigen::Triplet<float>(3*i + 1, 3*i + 1, -1.0f));
      trips.push_back(Eigen::Triplet<float>(3*i + 2, 3*i + 2, -1.0f));
      b[3*i + 0] = p1.x;
      b[3*i + 1] = p1.y;
      b[3*i + 2] = p1.z;
    } else {
      for (uint j = 0; j < neighs.size(); j++) {
        p2 = points[neighs[j]];
        if (constraints[neighs[j]]) {
          b[3*i + 0] -= p2.x;
          b[3*i + 1] -= p2.y;
          b[3*i + 2] -= p2.z;
        } else {
          // x
          trips.push_back(Eigen::Triplet<float>(3*i + 0, 3*neighs[j] + 0, 1.0f/neighs.size()));

          // y
          trips.push_back(Eigen::Triplet<float>(3*i + 1, 3*neighs[j] + 1, 1.0f/neighs.size()));

          // z
          trips.push_back(Eigen::Triplet<float>(3*i + 2, 3*neighs[j] + 2, 1.0f/neighs.size()));
        }

      }
      b[3*i + 0] /= neighs.size();
      b[3*i + 1] /= neighs.size();
      b[3*i + 2] /= neighs.size();
      trips.push_back(Eigen::Triplet<float>(3*i + 0, 3*i + 0, -1));
      trips.push_back(Eigen::Triplet<float>(3*i + 1, 3*i + 1, -1));
      trips.push_back(Eigen::Triplet<float>(3*i + 2, 3*i + 2, -1));
    }
  }
  A.setFromTriplets(trips.begin(), trips.end());
  // A.makeCompressed();
  Eigen::VectorXf x(3*points.size());
  
  cout << "Starting to compute..." << endl;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<float>>  solver;
  solver.compute(A);

  cout << "Starting to solve..." << endl;
  x = solver.solve(b);

  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

  for (uint i = 0; i < points.size(); i++) {
    if (constraints[i]) continue;
    mesh->getVertices()[i] = glm::vec3(x[3*i + 0], x[3*i + 1], x[3*i + 2]);
  }
}

/* This method has to optimize the vertices' positions in the least squares sense, 
   so that the laplacian is close to zero and the vertices remain close to their 
   original locations. The constraintWeight parameter is used to control how close 
   the vertices have to be to their original positions. */

void LaplacianSmoothing::globalBilaplacian(const vector<bool> &constraints, float constraintWeight)
{
  std::vector<glm::vec3> points = mesh->getVertices();
  uint constraints_size = 0;
  for (uint i = 0; i < points.size(); i++) {
    if (constraints[i]) constraints_size++;
  }
  Eigen::SparseMatrix<float> A(3*(points.size()+constraints_size), 3*points.size());
  std::vector<uint> neighs;
  glm::vec3 p1;
  glm::vec3 p2;
  Eigen::VectorXf b(3*(points.size()+constraints_size));

  cout << "Starting Matrix creation" << endl;
  // Create Gaussian matrix
  std::vector<Eigen::Triplet<float>> trips;
  for (uint i = 0; i < points.size(); i++) {
    p1 = points[i];
    mesh->getNeighbors(i, neighs);
    b[3*i + 0] = 0;
    b[3*i + 1] = 0;
    b[3*i + 2] = 0;
    for (uint j = 0; j < neighs.size(); j++) {
      p2 = points[neighs[j]];
      if (constraints[neighs[j]]) {
        b[3*i + 0] -= p2.x;
        b[3*i + 1] -= p2.y;
        b[3*i + 2] -= p2.z;
      } else {
        // x
        trips.push_back(Eigen::Triplet<float>(3*i + 0, 3*neighs[j] + 0, 1.0f/neighs.size()));

        // y
        trips.push_back(Eigen::Triplet<float>(3*i + 1, 3*neighs[j] + 1, 1.0f/neighs.size()));

        // z
        trips.push_back(Eigen::Triplet<float>(3*i + 2, 3*neighs[j] + 2, 1.0f/neighs.size()));
      }

    }
    b[3*i + 0] /= neighs.size();
    b[3*i + 1] /= neighs.size();
    b[3*i + 2] /= neighs.size();
    trips.push_back(Eigen::Triplet<float>(3*i + 0, 3*i + 0, -1));
    trips.push_back(Eigen::Triplet<float>(3*i + 1, 3*i + 1, -1));
    trips.push_back(Eigen::Triplet<float>(3*i + 2, 3*i + 2, -1));
  }
  uint count = 0;
  uint s = points.size();
  for (uint i = 0; i < points.size(); i++) {
    if (constraints[i]) {
      glm::vec3 p1 = points[i];
      b[s + 3*count + 0] = p1.x * constraintWeight;
      b[s + 3*count + 1] = p1.y * constraintWeight;
      b[s + 3*count + 2] = p1.z * constraintWeight;
      trips.push_back(Eigen::Triplet<float>(s + 3*count + 0, s + 3*count + 0, constraintWeight));
      trips.push_back(Eigen::Triplet<float>(s + 3*count + 1, s + 3*count + 1, constraintWeight));
      trips.push_back(Eigen::Triplet<float>(s + 3*count + 2, s + 3*count + 2, constraintWeight));
    }
  }


  A.setFromTriplets(trips.begin(), trips.end());
  // A.makeCompressed();
  Eigen::VectorXf x(3*(points.size()+constraints_size));
  
  cout << "Starting to compute..." << endl;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<float>>  solver;
  solver.compute(A.transpose() * A);

  cout << "Starting to solve..." << endl;
  x = solver.solve(A.transpose() * b);

  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

  for (uint i = 0; i < points.size(); i++) {
    if (constraints[i]) continue;
    mesh->getVertices()[i] = glm::vec3(x[3*i + 0], x[3*i + 1], x[3*i + 2]);
  }
  // Don't delete columns, simply add equations like p'_i = p_i where i is a constraint
}








