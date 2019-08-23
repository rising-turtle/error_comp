/*
	Aug. 21, 2019, He Zhang, hzhang8@vcu.edu 

	test jacobians in sampson 

*/

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "utility.h"
#include "sampson_factor.h"

using namespace std; 

void test_J1(); 

void test_d_res_d_pose(); 


int main(int argc, char* argv[])
{

	// test_J1(); //

	test_d_res_d_pose(); 

	return 0; 
}

void test_d_res_d_pose()
{
	Eigen::Vector3d Pi(0, 0, 0);
    Eigen::Quaterniond Qi(1, 0, 0, 0);

    Eigen::Vector3d Pj(0, 0, 0.2);
    Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 1;

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
    para[3][0] = 1;

    Eigen::Vector3d pts_i(0.5, 0.5, 1.); 
    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
	double dep_j = pts_camera_j.z();
    Eigen::Vector3d pts_j(pts_camera_j.x()/dep_j - 0.1, pts_camera_j.y()/dep_j -0.3, 1.);
 
    SampsonFactor * f = new SampsonFactor(pts_i, pts_j);
    f->check(para);
    return ; 
}


void test_J1()
{

	Eigen::Vector3d Pi(0, 0, 0);
    Eigen::Quaterniond Qi(1, 0, 0, 0);

    Eigen::Vector3d Pj(0, 0, 0.2);
    Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 1;

    Eigen::Vector3d pts_i(0.5, 0.5, 1.); 
    Eigen::Vector3d pts_i_gt(0.5, 0.5, 1.); 

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    

	double dep_j = pts_camera_j.z();
    Eigen::Vector3d pts_j(pts_camera_j.x()/dep_j - 0.1, pts_camera_j.y()/dep_j -0.3, 1.);
    Eigen::Vector2d epsilon = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
    Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
    d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
    						0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
    Eigen::Matrix<double, 3, 3> d_pts_cam_j_d_pts_cam_i; 
    d_pts_cam_j_d_pts_cam_i = qic.inverse()*Qj.inverse()*Qi*qic; 
    
    Eigen::Matrix<double, 3, 2> d_pts_cam_i_d_pts_i; 
    d_pts_cam_i_d_pts_i << 1./inv_dep_i, 0, 
    					   0, 1./inv_dep_i,
    					   0, 0;

    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
   	Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_j = Eigen::Matrix<double,2,2>::Identity()*-1; 
    Eigen::Matrix<double, 2, 4> d_epsilon_d_X; 
    d_epsilon_d_X.block<2,2>(0, 0) = d_epsilon_d_pts_i; 
    d_epsilon_d_X.block<2,2>(0, 2) = d_epsilon_d_pts_j;

    const double eps = 1e-6; 
    Eigen::Matrix<double, 2, 4> num_jacobians; 

    Eigen::Vector4d X(pts_i(0), pts_i(1), pts_j(0), pts_j(1)); 

    for(int k=0; k<4; k++){

    	Eigen::Vector4d delta = Eigen::Vector4d(k==0, k==1, k==2, k==3)*eps; 
    	Eigen::Vector4d tmpX = X+delta; 
    	pts_i << tmpX(0), tmpX(1), 1.; 
    	pts_j << tmpX(2), tmpX(3), 1.;

    	pts_camera_i = pts_i / inv_dep_i;
    	pts_imu_i = qic * pts_camera_i + tic;
    	pts_w = Qi * pts_imu_i + Pi;
     	pts_imu_j = Qj.inverse() * (pts_w - Pj);
    	pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    	dep_j = pts_camera_j.z();

    	Eigen::Vector2d tmp_epsilon = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    	num_jacobians.col(k) = (tmp_epsilon - epsilon)/eps; 

    }

    cout <<"analytal jacobians: "<<endl; 

    cout <<d_epsilon_d_X<<endl;

    cout <<"numeric jacobians: "<<endl; 

    cout <<num_jacobians<<endl; 


    // test d_J_d_qi 
    Eigen::Matrix<double, 2, 2> J = d_epsilon_d_X.block<2,2>(0,0); 

    // J = d_epsilon_d_pts_cam_j *qic.inverse()*Qj.inverse()*Qi*qic * d_pts_cam_i_d_pts_i; 

    Eigen::Matrix<double, 2, 3> d_Jc1_d_qi; 
    Eigen::Matrix<double, 3, 1> t1 = qic * d_pts_cam_i_d_pts_i.col(0); 
    // d_Jc1_d_qi = d_epsilon_d_pts_cam_j * qic.inverse()*Qj.inverse() * Qi * Utility::skewSymmetric(-t1);
    Eigen::Matrix3d tmp ; 
    tmp = qic.inverse()*Qj.inverse() * Qi.toRotationMatrix()*Utility::skewSymmetric(-t1);
    d_Jc1_d_qi = d_epsilon_d_pts_cam_j * tmp; 
	// d_Jc1_d_qi = d_epsilon_d_pts_cam_j * qic.inverse()*Qj.inverse() * Qi.toRotationMatrix(); //Utility::skewSymmetric(-t1);
    Eigen::Matrix<double, 2, 3> num_jacobians_Jc1; 

    for(int k=0; k<3; k++){

    	Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
    	Eigen::Quaterniond tmp_Qi = Qi * Utility::deltaQ(delta); 
   		
   		Eigen::Matrix<double, 3, 3> tmp_m = qic.inverse()*Qj.inverse()*tmp_Qi.toRotationMatrix()*qic.toRotationMatrix();
   		Eigen::Matrix<double, 2, 2> tmp_J = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_pts_i;

    	Eigen::Vector2d tmp_diff_J = tmp_J.col(0)  - J.col(0);

    	num_jacobians_Jc1.col(k) = tmp_diff_J/eps; 
    }

    cout <<"analytal d_Jc1_d_qi: "<<endl;

    cout << d_Jc1_d_qi <<endl; 

    cout <<"numeric d_Jc1_d_qi: "<<endl;

    cout <<num_jacobians_Jc1<<endl;

    Eigen::Matrix<double, 2, 3> d_Jc2_d_qi; 
    Eigen::Matrix<double, 3, 1> t2 = qic * d_pts_cam_i_d_pts_i.col(1); 
    // d_Jc1_d_qi = d_epsilon_d_pts_cam_j * qic.inverse()*Qj.inverse() * Qi * Utility::skewSymmetric(-t1);
    tmp = qic.inverse()*Qj.inverse() * Qi.toRotationMatrix()*Utility::skewSymmetric(-t2);
    d_Jc2_d_qi = d_epsilon_d_pts_cam_j * tmp; 
	// d_Jc1_d_qi = d_epsilon_d_pts_cam_j * qic.inverse()*Qj.inverse() * Qi.toRotationMatrix(); //Utility::skewSymmetric(-t1);
    Eigen::Matrix<double, 2, 3> num_jacobians_Jc2; 

    for(int k=0; k<3; k++){

    	Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
    	Eigen::Quaterniond tmp_Qi = Qi * Utility::deltaQ(delta); 
   		
   		Eigen::Matrix<double, 3, 3> tmp_m = qic.inverse()*Qj.inverse()*tmp_Qi.toRotationMatrix()*qic.toRotationMatrix();
   		Eigen::Matrix<double, 2, 2> tmp_J = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_pts_i;

    	Eigen::Vector2d tmp_diff_J = tmp_J.col(1)  - J.col(1);

    	num_jacobians_Jc2.col(k) = tmp_diff_J/eps; 
    }

    cout <<"analytal d_Jc2_d_qi: "<<endl;

    cout << d_Jc2_d_qi <<endl; 

    cout <<"numeric d_Jc2_d_qi: "<<endl;

    cout <<num_jacobians_Jc2<<endl;


    Eigen::Matrix<double, 2, 2> JtJ = d_epsilon_d_X*d_epsilon_d_X.transpose(); 
    Eigen::Vector4d res = -d_epsilon_d_X.transpose()*JtJ.inverse()*epsilon;

    cout <<"res = "<<res.transpose()<<endl;
    cout <<"J.inverse() = "<< d_epsilon_d_X.inverse()<<endl;
    // Eigen::Vector4d simple_res = -d_epsilon_d_X.inverse()*epsilon; 

    // cout <<"simple_res = "<<simple_res.transpose()<<endl;


}