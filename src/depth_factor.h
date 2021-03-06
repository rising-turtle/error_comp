/*
	Oct. 25 2019, He Zhang, hzhang8@vcu.edu 

	depth based factor 

*/

#pragma once 

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "utility.h"


class SingleInvDepthFactor : public ceres::SizedCostFunction<1, 1>
{
 public:
 	SingleInvDepthFactor(double inv_i); 
 	virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

 	double inv_depth_i;
 	static double sqrt_info; 
};


class ProjectionDepthFactor : public ceres::SizedCostFunction<3, 7, 7, 7, 1>
{
public:
	ProjectionDepthFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, double inv_j); 

	virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
	void check(double **parameters); 

	// double inv_depth_i; 
	double inv_depth_j; 
	Eigen::Vector3d pts_i; 
	Eigen::Vector3d pts_j; 
	// static double sqrt_info; 
	static Eigen::Matrix3d sqrt_info; 
};

