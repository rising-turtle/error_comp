/*
	Oct. 25 2019, He Zhang, hzhang8@vcu.edu 

	depth based factor 

*/

#include "depth_factor.h"


double ProjectionDepthFactor::sqrt_info = 1.; // need to set before use this factor 

ProjectionDepthFactor::ProjectionDepthFactor(double inv_i, double inv_j):
inv_depth_i(inv_i),
inv_depth_j(inv_j)
{}

bool ProjectionDepthFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    Eigen::Vector3d pts_camera_i = pts_i_td / inv_depth_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic2.inverse() * (pts_imu_j - tic2);

    residuals[0] = (dep_j) - 1./inv_depth_j;

    
}