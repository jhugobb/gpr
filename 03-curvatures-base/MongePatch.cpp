#include <iostream>
#include "MongePatch.h"
#include <Eigen/Dense>

// Given a point P, its normal, and its closest neighbors (including itself) 
// compute a quadratic Monge patch that approximates the neighborhood of P.
// The resulting patch will be used to compute the principal curvatures of the 
// surface at point P.

void MongePatch::init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest)
{
	// N = (nx, ny, nz)
	// calc max{nx, ny, nz}
	// for example, if nx is greater,
	// 							use (0,1,0)
	//              if ny is greater,
	//              use (0,0,1)
	// 							else use (1,0,0)
	glm::vec3 w = -normal;
	w = glm::normalize(w);
	glm::vec3 g;
	glm::vec3 u;
	float m = max(max(abs(w.x), abs(w.y)), abs(w.z));

	if (m == abs(w.x)) {
		g = glm::vec3(0,1,0);
	} else if (m == abs(w.y)) {
		g = glm::vec3(0,0,1);
	} else {
		g = glm::vec3(1,0,0);
	}
	
	u = glm::cross(g, w);
	u = glm::normalize(u);
	
	glm::vec3 v = glm::cross(w, u);
	v = glm::normalize(v);

	// Transform points to the new coordinate system (u,v,w)
	std::vector<glm::vec3> trans_points;
	glm::vec3 t;
	for (glm::vec3 q : closest) {
		t = glm::vec3(glm::dot(q-P, u), glm::dot(q-P, v), glm::dot(q-P, w));
		trans_points.push_back(t);
	}

  // Create the linear system Ax = b
	// A = w(u,v)
	// x = S
	// b = sum(w_i * q_i)
	Eigen::MatrixXd wuv(6,6);
	wuv = Eigen::MatrixXd::Zero(6,6);
	Eigen::VectorXd S(6);
	S = Eigen::VectorXd::Zero(6);
	Eigen::VectorXd wqs(6);
	wqs = Eigen::VectorXd::Zero(6);

	for (glm::vec3 p : trans_points) {
		Eigen::VectorXd q(6);
		q = Eigen::VectorXd::Zero(6);
		q << pow(p.x,2), p.x*p.y, pow(p.y,2), p.x, p.y, 1;
		wuv += q * q.transpose();
		wqs += p.z * q;
	}

	// Solve system
	S = wuv.colPivHouseholderQr().solve(wqs);

	// Define Hessian
	Eigen::Matrix2d H;

	double a = S(0);
	double b = S(1);
	double c = S(2);

	H << 2*a, b,
			 b, 2*c;

	// Get eigenvalues
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(H);
	solver.compute(H);

	Eigen::Vector2d vals = solver.eigenvalues();
	
	// Store result for later
	kmin_ = vals(0);
	kmax_ = vals(1);

}

// Return the values of the two principal curvatures for this patch

void MongePatch::principalCurvatures(float &kmin, float &kmax) const
{
	kmin = kmin_;
	kmax = kmax_;
}


