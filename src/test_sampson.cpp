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

#define SQ(x) (((x)*(x)))

void test_J1(); 

void test_d_res_d_pose(); 

void test_J1_lambda();

void test_d_res_d_pose_lambda();

void test_sampson_jacobian(); 

int main(int argc, char* argv[])
{

    cout.precision(8);

	// test_J1(); //

	// test_d_res_d_pose(); 

    test_sampson_jacobian(); 

    // test_J1_lambda(); 

    // test_d_res_d_pose_lambda();

	return 0; 
}


void test_sampson_jacobian()
{
    Eigen::Vector3d Pi(1, 2, 3);
    Eigen::Quaterniond Qi(sqrt(1-SQ(0.1)-SQ(0.3)-SQ(0.5)), 0.1, -0.3, 0.5);

    Eigen::Vector3d Pj(0.1, 0.2, 0.3);
    // Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);
    Eigen::Quaterniond Qj(sqrt(1-SQ(0.4)-SQ(0.1)-SQ(0.2)), 0.4, 0.1, -0.2);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 1/3.; // 0.5

    double ** para = new double *[4];//*[4]; 

    para[0] = new double[7]; // pose_i 
    para[0][0] = Pi(0); para[0][1] = Pi(1); para[0][2] = Pi(2); 
    para[0][3] = Qi.x(); para[0][4] = Qi.y(); para[0][5] = Qi.z(); para[0][6] = Qi.w(); 
    
    para[1] = new double[7]; // pose_j
    para[1][0] = Pj(0); para[1][1] = Pj(1); para[1][2] = Pj(2); 
    para[1][3] = Qj.x(); para[1][4] = Qj.y(); para[1][5] = Qj.z(); para[1][6] = Qj.w(); 
    
    para[2] = new double[7]; // T_I2C
    para[2][0] = tic(0); para[2][1] = tic(1); para[2][2] = tic(2); 
    para[2][3] = qic.x(); para[2][4] = qic.y(); para[2][5] = qic.z(); para[2][6] = qic.w(); 
    
    para[3] = new double[1]; 
    para[3][0] = inv_dep_i;
    Eigen::Vector3d pts_i(0.1, 0.2, 1.);
    Eigen::Vector3d pts_j(-0.3, 1., 1.0);

   {
        SampsonFactorCross * f = new SampsonFactorCross(pts_i, pts_j);
        // f->compute_Jacobian_pose_i(para, NULL); 
        f->check(para);
    }

}
void test_d_res_d_pose()
{
	/*Eigen::Vector3d Pi(0, 0, 0);
    Eigen::Quaterniond Qi(1, 0, 0, 0);

    Eigen::Vector3d Pj(0, 0, 0.2);
    Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 0.5;
*/

    Eigen::Vector3d Pi(0, 0, 0);
    Eigen::Quaterniond Qi(1, 0, 0, 0);

    Eigen::Vector3d Pj(0.353549, -0.118023, 2.12716);
    // Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);
    Eigen::Quaterniond Qj(0.9174208456426309, 0.3488320, -0.187413, 0.0391356);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 0.470548; // 0.5

    double ** para = new double *[4];//*[4]; 

  	para[0] = new double[7]; // pose_i 
    para[0][0] = 0; para[0][1] = 0; para[0][2] = 0; 
    para[0][3] = 0; para[0][4] = 0; para[0][5] = 0; para[0][6] = 1; 
    
    para[1] = new double[7]; // pose_j
    para[1][0] = Pj(0); para[1][1] = Pj(1); para[1][2] = Pj(2); 
    para[1][3] = Qj.x(); para[1][4] = Qj.y(); para[1][5] = Qj.z(); para[1][6] = Qj.w(); 
  	
  	para[2] = new double[7]; // T_I2C
    para[2][0] = tic(0); para[2][1] = tic(1); para[2][2] = tic(2); 
    para[2][3] = qic.x(); para[2][4] = qic.y(); para[2][5] = qic.z(); para[2][6] = qic.w(); 
    
    para[3] = new double[1]; 
    para[3][0] = inv_dep_i;

    // Eigen::Vector3d pts_i(0.5, 0.5, 1.); 
    Eigen::Vector3d pts_i(-1.3163, 0.660394, 1.); 
    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
	double dep_j = pts_camera_j.z();
    // Eigen::Vector3d pts_j(pts_camera_j.x()/dep_j - 0.1, pts_camera_j.y()/dep_j -0.3, 1.);
    Eigen::Vector3d pts_j(-0.418893, 0.0679902, 1.0);

    {
        // cout <<" check SampsonFactor: "<<endl; 
        // SampsonFactorEssential * f = new SampsonFactorEssential(pts_i, pts_j);
        // f->check(para);
    }
    {
        cout<<" check SampsonFactorCross: "<<endl; 
        SampsonFactorCross * f = new SampsonFactorCross(pts_i, pts_j);
        f->check(para); 
    }

    {
        // cout <<" check SampsonFactorCrossWithLambda: "<<endl;
        // SampsonFactorCrossWithLambda * f = new SampsonFactorCrossWithLambda(pts_i, pts_j);
        // f->check(para);
    }

    return ; 
}

void test_d_res_d_pose_lambda()
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
    para[3][0] = 1;

    Eigen::Vector3d pts_i(-0.378468, -0.227461, 1.); 
    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
    double dep_j = pts_camera_j.z();
    // Eigen::Vector3d pts_j(pts_camera_j.x()/dep_j - 0.1, pts_camera_j.y()/dep_j -0.3, 1.);
    Eigen::Vector3d pts_j(-0.310314, 0.332055, 1.0);

    // SampsonFactorWithLambda * f = new SampsonFactorWithLambda(pts_i, pts_j);
    SampsonFactorCrossWithLambda * f = new SampsonFactorCrossWithLambda(pts_i, pts_j);
    f->check(para);
    return ; 
}


void test_J1_lambda()
{
    Eigen::Vector3d Pi(0, 0, 0);
    Eigen::Quaterniond Qi(1, 0, 0, 0);

    Eigen::Vector3d Pj(0, 0, 0.2);
    Eigen::Quaterniond Qj(0.98, 0, 0, 0.199);

    Eigen::Vector3d tic(0, 0, 0);
    Eigen::Quaterniond qic(1, 0, 0, 0);

    double inv_dep_i = 0.5;

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

    Eigen::Matrix<double, 3, 1> d_pts_cam_i_d_lambda; 
    double lambda_isq = 1./(inv_dep_i*inv_dep_i); 
    d_pts_cam_i_d_lambda = pts_i *(-lambda_isq); 


    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_j = Eigen::Matrix<double,2,2>::Identity()*-1; 
    Eigen::Matrix<double, 2, 1> d_epsilon_d_lambda = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_lambda;
    Eigen::Matrix<double, 2, 5> d_epsilon_d_X; 
    d_epsilon_d_X.block<2,2>(0, 0) = d_epsilon_d_pts_i; 
    d_epsilon_d_X.block<2,2>(0, 2) = d_epsilon_d_pts_j;
    d_epsilon_d_X.block<2,1>(0, 4) = d_epsilon_d_lambda;

    const double eps = 1e-7; 

    {
        Eigen::Vector3d delta_pts_i = d_pts_cam_i_d_lambda * eps; 
        Eigen::Vector3d new_pi = pts_camera_i + delta_pts_i; 
        cout <<" new pts_camera_i: "<<pts_camera_i<<endl; 

        double new_ind = inv_dep_i + eps; 
        Eigen::Vector3d new_nu_pi = pts_i / new_ind; 
        cout <<" new numeric pts_camera_i "<<new_nu_pi<<endl;
    }

    Eigen::Matrix<double, 2, 5> num_jacobians; 

    Eigen::Matrix<double, 5, 1> X; X << pts_i(0), pts_i(1), pts_j(0), pts_j(1), inv_dep_i; 

    for(int k=0; k<5; k++){

        Eigen::Matrix<double, 5, 1> delta = Eigen::Matrix<double, 5, 1>::Zero(); 
        delta(k) = eps; // = Eigen::Vector4d(k==0, k==1, k==2, k==3, k==4)*eps; 
        Eigen::Matrix<double, 5, 1> tmpX = X+delta; 
        pts_i << tmpX(0), tmpX(1), 1.; 
        pts_j << tmpX(2), tmpX(3), 1.;
        inv_dep_i = tmpX(4); 

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

/*
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
    cout <<"J.inverse() = "<< d_epsilon_d_X.inverse()<<endl;*/
    // Eigen::Vector4d simple_res = -d_epsilon_d_X.inverse()*epsilon; 

    // cout <<"simple_res = "<<simple_res.transpose()<<endl;


}