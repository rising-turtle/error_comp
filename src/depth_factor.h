/*
	Oct. 25 2019, He Zhang, hzhang8@vcu.edu 

	depth based factor 

*/

#pragma once 

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "utility.h"


class SingleDepthFactor : public ceres::SizedCostFunction<1, 7, 1>
{
public:
	SingleDepthFactor(const double T, )
};


class ProjectionDepthFactor : public ceres::SizedCostFunction<1, 7, 7, 7, 1>
{
public:
	ProjectionDepthFactor(double inv_depth_i, double inv_depth_j); 

	virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
	void check(double **parameters); 

	double inv_depth_i; 
	double inv_depth_j; 
	static double sqrt_info; 
};

