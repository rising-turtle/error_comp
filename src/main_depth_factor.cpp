/*

	Oct. 24 2019, He Zhang, hzhang8@vcu.edu 

	test depth factor, depth are available in multiple frames 

*/

#include "estimator_depth_factor.h"
#include "cases.h"
#include "mean_std.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 
#include "depth_factor.h"

double focal_length = 240.; 

double run_once(int n_feature, double noise,  double& err_dis, double& err_ori ); 
void test_depth_factor(); 


int main(int argc, char* argv[])
{

	double e1, e2;
    run_once(50, 2./focal_length,e1, e2 ); 
	// test_depth_factor(); 
	return 1; 

}


double run_once(int n_feature, double noise, double& err_dis, double& err_ori )
{

	Case* pc = new Case_forward(); 
	pc->init_pose_features(n_feature);
	pc->gen_observations();
	pc->add_noise(noise);

	double err[2]; 
	{
		EstimatorDepth est(noise); 
		est.optimize(pc, err); 
	}
	cout<<endl<<"run it again!"<<endl; 
	{
		EstimatorDepth est(noise); 
		est.m_depth_factor_type = RGBD; 
		est.optimize(pc, err); 
	}
	err_dis = err[0]; 
	err_ori = err[1]; 

	delete pc;

	return err_dis;
}


void test_depth_factor()
{
    Eigen::Vector3d Pi(0, 0, 0);
    Eigen::Quaterniond Qi(1, 0, 0, 0);

    Eigen::Vector3d Pj(0, 0, 1);
    // Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);
    Eigen::Quaterniond Qj(1, 0, 0, 0);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 0.159081; // 0.5

    double ** para = new double *[4];//*[4]; 

    para[0] = new double[7]; // pose_i 
    para[0][0] = 0; para[0][1] = 0; para[0][2] = 0; 
    para[0][3] = 0; para[0][4] = 0; para[0][5] = 0; para[0][6] = 1; 
    
    para[1] = new double[7]; // pose_j
    para[1][0] = 0; para[1][1] = 0; para[1][2] = 0.2; 
    para[1][3] = 0; para[1][4] = 0; para[1][5] = 0.199; para[1][6] = 0.98; 
    
    para[2] = new double[7]; // T_I2C
    para[2][0] = tic(0); para[2][1] = tic(1); para[2][2] = tic(2); 
    para[2][3] = qic.x(); para[2][4] = qic.y(); para[2][5] = qic.z(); para[2][6] = qic.w(); 
    
    para[3] = new double[1]; 
    para[3][0] = 1; // depth 

    Eigen::Vector3d pts_i(-0.378468, -0.227461, 1.); 
    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
    double dep_j = pts_camera_j.z();
    double inv_dep_j = 1./dep_j;
    // Eigen::Vector3d pts_j(pts_camera_j.x()/dep_j - 0.1, pts_camera_j.y()/dep_j -0.3, 1.);
    // Eigen::Vector3d pts_j(-0.310314, 0.332055, 1.0);

    ProjectionDepthFactor * f = new ProjectionDepthFactor(pts_i, inv_dep_j);
    f->check(para);
    return ; 
}
