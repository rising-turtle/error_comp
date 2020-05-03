/*
	May. 2, 2020, He Zhang, hzhang8@vcu.edu 

	analytical computation of jacobians w.r.t Pose_i and Pose_j 

*/

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "utility.h"
#include "sampson_factor.h"

using namespace std; 
using namespace Eigen;
namespace{

	Eigen::Matrix<double, 2, 4> J_func(double const *const *parameters, const Vector3d& pts_i, const Vector3d& pts_j)
	{
		Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
	    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

	    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
	    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

	    double inv_dep_i = parameters[3][0];

	    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
	    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
	    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
	    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
	    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
    	double dep_j = pts_camera_j.z();

	    // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
	    Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
	    // d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
	    //                        0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
	    d_epsilon_d_pts_cam_j << 0., 1., -pts_j.y(), 
	                             -1, 0, pts_j.x(); 

	    Eigen::Matrix<double, 3, 3> d_pts_cam_j_d_pts_cam_i; 
	    d_pts_cam_j_d_pts_cam_i = qic.inverse()*Qj.inverse()*Qi*qic; 
	    Eigen::Matrix<double, 3, 2> d_pts_cam_i_d_pts_i; 
	    d_pts_cam_i_d_pts_i << 1./inv_dep_i, 0, 
	                           0, 1./inv_dep_i,
	                           0, 0;

	    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
	    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_j = Eigen::Matrix<double,2,2>::Identity()*-1; 
	    d_epsilon_d_pts_j << 0, -dep_j,
	                        dep_j, 0; 

	    Eigen::Matrix<double, 2, 4> d_epsilon_d_X; 
	    d_epsilon_d_X.block<2,2>(0, 0) = d_epsilon_d_pts_i; 
	    d_epsilon_d_X.block<2,2>(0, 2) = d_epsilon_d_pts_j;

	    return d_epsilon_d_X; 

	}

	void jacobian_J(double const *const *parameters, const Vector3d& pts_i, const Vector3d& pts_j, Eigen::Matrix<double, 8, 6>& dJ_dpose_i, Eigen::Matrix<double, 8, 6>& dJ_dpose_j)
	{
		Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
	    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

	    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
	    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

	    double inv_dep_i = parameters[3][0];

	    // Eigen::Matrix3d Rji_cam = qic.inverse()*Qj.inverse()*Qi*qic; 

	    Eigen::Matrix<double, 2, 3> A; 
	    A << 0, 1, -pts_j(1), 
	    	-1, 0, pts_j(0);

	    Eigen::Vector3d d1; d1 << 1./inv_dep_i,  0, 0; 
	    Eigen::Vector3d d2; d2 << 0., 1./inv_dep_i, 0; 

	   	Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();

	    // J = [J1 J2 J3 J4]
	    Eigen::Matrix<double, 2, 3> dJ1_dtheta_i = -A*ric.transpose()*Rj.transpose()*Ri*Utility::skewSymmetric(qic*d1); 
	    Eigen::Matrix<double, 2, 3> dJ1_dtheta_j = A*ric.transpose()*Utility::skewSymmetric(Qj.inverse()*Qi*qic*d1); 
	    Eigen::Matrix<double, 2, 3> dJ1_dti = Eigen::Matrix<double, 2, 3>::Zero(); 
	    Eigen::Matrix<double, 2, 3> dJ1_dtj = Eigen::Matrix<double, 2, 3>::Zero(); 

	   	Eigen::Matrix<double, 2, 3> dJ2_dtheta_i = -A*ric.transpose()*Rj.transpose()*Ri*Utility::skewSymmetric(qic*d2); 
	    Eigen::Matrix<double, 2, 3> dJ2_dtheta_j = A*ric.transpose()*Utility::skewSymmetric(Qj.inverse()*Qi*qic*d2); 
	    Eigen::Matrix<double, 2, 3> dJ2_dti = Eigen::Matrix<double, 2, 3>::Zero(); 
	    Eigen::Matrix<double, 2, 3> dJ2_dtj = Eigen::Matrix<double, 2, 3>::Zero(); 

	    // cout<<"dJ1_dtheta_j: "<<endl<<dJ1_dtheta_j<<endl; 
	    // cout<<"dJ2_dtheta_j: "<<endl<<dJ2_dtheta_j<<endl;

	    // d_epsilon_d_pj
	    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    	Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    	Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    	Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    	Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
    	// double dep_j = pts_camera_j.z();
    	// Eigen::Vector2d epsilon = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    	// using cross product as epsilon 
    	// Eigen::Vector2d epsilon(pts_camera_j.y() - dep_j*pts_j.y(), -pts_camera_j.x() + dep_j * pts_j.x()); 

    	// d_pj_d_pose_i 
    	Eigen::Matrix3d d_pj_d_theta_i = - ric.transpose() * Rj.transpose() * Ri * Utility::skewSymmetric(pts_imu_i);
    	Eigen::Matrix3d d_pj_d_ti = ric.transpose()  * Rj.transpose() ; 
    	Eigen::Matrix3d d_pj_d_theta_j = ric.transpose() *Utility::skewSymmetric(pts_imu_j); 
    	Eigen::Matrix3d d_pj_d_tj = -ric.transpose()  * Rj.transpose(); 

    	// Vector3d tji = Rj.transpose()*(Pi - Pj); 
    	// cout<<"Pi: "<<Pi<<endl;
    	// cout <<"Pj: "<<Pj<<endl;
    	// cout <<"Ri:" <<endl<< Ri<<endl; 
    	// cout <<"Rj': "<<Rj.transpose()<<endl;
    	// cout <<"tji: "<<tji<<endl;

    	// cout<<"pts_imu_j: "<<endl<<pts_imu_j<<endl;
    	// cout<<"d_pj_d_theta_j: "<<endl<<d_pj_d_theta_j<<endl; 

    	Eigen::Matrix<double, 1, 3> r3(0, 0, 1); 
    	Eigen::Matrix<double, 1, 3> d_zj_d_theta_i = r3 * d_pj_d_theta_i; 
    	Eigen::Matrix<double, 1, 3> d_zj_d_ti = r3 * d_pj_d_ti; 

    	Eigen::Matrix<double, 4, 6> dJ34_dpose_i = Eigen::Matrix<double, 4, 6>::Zero(); 
    	dJ34_dpose_i.block<1,3>(1, 0) = -d_zj_d_ti; dJ34_dpose_i.block<1,3>(1,3) = -d_zj_d_theta_i;
    	dJ34_dpose_i.block<1,3>(2, 0) =  d_zj_d_ti; dJ34_dpose_i.block<1,3>(2,3) =  d_zj_d_theta_i;

    	// d_pj_d_pose_j
    	Eigen::Matrix<double, 1, 3> d_zj_d_theta_j = r3 * d_pj_d_theta_j; 
    	Eigen::Matrix<double, 1, 3> d_zj_d_tj = r3 * d_pj_d_tj; 

    	Eigen::Matrix<double, 4, 6> dJ34_dpose_j = Eigen::Matrix<double, 4, 6>::Zero(); 
    	dJ34_dpose_j.block<1,3>(1, 0) = -d_zj_d_tj; dJ34_dpose_j.block<1,3>(1,3) = -d_zj_d_theta_j;
    	dJ34_dpose_j.block<1,3>(2, 0) =  d_zj_d_tj; dJ34_dpose_j.block<1,3>(2,3) =  d_zj_d_theta_j;

    	// cout<<"d_zj_d_theta_j: "<<endl<<d_zj_d_theta_j<<endl; 
    	// cout<<"d_zj_d_tj: "<<endl<<d_zj_d_tj<<endl;
    	// cout<<"dJ34_dpose_j: "<<endl<<dJ34_dpose_j<<endl; 

    	// finally combine the result 
    	dJ_dpose_i.block<1,3>(0,0) = dJ1_dti.block<1,3>(0,0); dJ_dpose_i.block<1,3>(0,3) = dJ1_dtheta_i.block<1,3>(0,0); 
    	dJ_dpose_i.block<1,3>(1,0) = dJ2_dti.block<1,3>(0,0); dJ_dpose_i.block<1,3>(1,3) = dJ2_dtheta_i.block<1,3>(0,0);
    	dJ_dpose_i.block<2,6>(2,0) = dJ34_dpose_i.block<2,6>(0,0);
    	dJ_dpose_i.block<1,3>(4,0) = dJ1_dti.block<1,3>(1,0); dJ_dpose_i.block<1,3>(4,3) = dJ1_dtheta_i.block<1,3>(1,0);
    	dJ_dpose_i.block<1,3>(5,0) = dJ2_dti.block<1,3>(1,0); dJ_dpose_i.block<1,3>(5,3) = dJ2_dtheta_i.block<1,3>(1,0);
    	dJ_dpose_i.block<2,6>(6,0) = dJ34_dpose_i.block<2,6>(2,0);

    	dJ_dpose_j.block<1,3>(0,0) = dJ1_dtj.block<1,3>(0,0); dJ_dpose_j.block<1,3>(0,3) = dJ1_dtheta_j.block<1,3>(0,0); 
    	dJ_dpose_j.block<1,3>(1,0) = dJ2_dtj.block<1,3>(0,0); dJ_dpose_j.block<1,3>(1,3) = dJ2_dtheta_j.block<1,3>(0,0);
    	dJ_dpose_j.block<2,6>(2,0) = dJ34_dpose_j.block<2,6>(0,0);
    	dJ_dpose_j.block<1,3>(4,0) = dJ1_dtj.block<1,3>(1,0); dJ_dpose_j.block<1,3>(4,3) = dJ1_dtheta_j.block<1,3>(1,0);
    	dJ_dpose_j.block<1,3>(5,0) = dJ2_dtj.block<1,3>(1,0); dJ_dpose_j.block<1,3>(5,3) = dJ2_dtheta_j.block<1,3>(1,0);
    	dJ_dpose_j.block<2,6>(6,0) = dJ34_dpose_j.block<2,6>(2,0);

    	return ; 
	}

	void jacobian_J_transpose(Eigen::Matrix<double, 8, 6>& dtmp_dpose_i, Eigen::Matrix<double, 8, 6>& dJ_dpose_i)
	{
		dJ_dpose_i.block<1,6>(0,0) = dtmp_dpose_i.block<1,6>(0,0); 
		dJ_dpose_i.block<1,6>(1,0) = dtmp_dpose_i.block<1,6>(4,0); 
		dJ_dpose_i.block<1,6>(2,0) = dtmp_dpose_i.block<1,6>(1,0);
		dJ_dpose_i.block<1,6>(3,0) = dtmp_dpose_i.block<1,6>(5,0);
		dJ_dpose_i.block<1,6>(4,0) = dtmp_dpose_i.block<1,6>(2,0);
		dJ_dpose_i.block<1,6>(5,0) = dtmp_dpose_i.block<1,6>(6,0);
		dJ_dpose_i.block<1,6>(6,0) = dtmp_dpose_i.block<1,6>(3,0);     
		dJ_dpose_i.block<1,6>(7,0) = dtmp_dpose_i.block<1,6>(7,0); 
	}

	void jacobian_Jt(double const *const *parameters, const Vector3d& pts_i, const Vector3d& pts_j, Eigen::Matrix<double, 8, 6>& dJ_dpose_i, Eigen::Matrix<double, 8, 6>& dJ_dpose_j)
	{
		Eigen::Matrix<double, 8, 6> dtmp_dpose_i, dtmp_dpose_j; 
		jacobian_J(parameters, pts_i, pts_j, dtmp_dpose_i, dtmp_dpose_j); 

		jacobian_J_transpose(dtmp_dpose_i, dJ_dpose_i); 
		jacobian_J_transpose(dtmp_dpose_j, dJ_dpose_j);
		// dJ_dpose_i.block<1,6>(0,0) = dtmp_dpose_i.block<1,6>(0,0); 
		// dJ_dpose_i.block<1,6>(1,0) = dtmp_dpose_i.block<1,6>(4,0); 
		// dJ_dpose_i.block<1,6>(2,0) = dtmp_dpose_i.block<1,6>(1,0);
		// dJ_dpose_i.block<1,6>(3,0) = dtmp_dpose_i.block<1,6>(5,0);
		// dJ_dpose_i.block<1,6>(4,0) = dtmp_dpose_i.block<1,6>(2,0);
		// dJ_dpose_i.block<1,6>(5,0) = dtmp_dpose_i.block<1,6>(6,0);
		// dJ_dpose_i.block<1,6>(6,0) = dtmp_dpose_i.block<1,6>(3,0);     
		// dJ_dpose_i.block<1,6>(7,0) = dtmp_dpose_i.block<1,6>(7,0); 

		// dJ_dpose_j.block<1,6>(0,0) = dtmp_dpose_j.block<1,6>(0,0); 
		// dJ_dpose_j.block<1,6>(1,0) = dtmp_dpose_j.block<1,6>(4,0); 
		// dJ_dpose_j.block<1,6>(2,0) = dtmp_dpose_j.block<1,6>(1,0);
		// dJ_dpose_j.block<1,6>(3,0) = dtmp_dpose_j.block<1,6>(5,0);
		// dJ_dpose_j.block<1,6>(4,0) = dtmp_dpose_j.block<1,6>(2,0);
		// dJ_dpose_j.block<1,6>(5,0) = dtmp_dpose_j.block<1,6>(6,0);
		// dJ_dpose_j.block<1,6>(6,0) = dtmp_dpose_j.block<1,6>(3,0);     
		// dJ_dpose_j.block<1,6>(7,0) = dtmp_dpose_j.block<1,6>(7,0); 
	}

	void jacobian_JJt(double const *const *parameters, const Vector3d& pts_i, const Vector3d& pts_j, Eigen::Matrix<double, 4, 6>& dJJt_dpose_i, Eigen::Matrix<double, 4, 6>& dJJt_dpose_j)
	{
		Eigen::Matrix<double, 2, 4> J = J_func(parameters, pts_i, pts_j); 
		Eigen::Matrix<double, 4, 2> Jt = J.transpose(); 

		Eigen::Matrix<double, 8, 6> dJ_dpose_i, dJ_dpose_j, dJt_dpose_i, dJt_dpose_j; 
		jacobian_J(parameters, pts_i, pts_j, dJ_dpose_i, dJ_dpose_j); 
		jacobian_J_transpose(dJ_dpose_i, dJt_dpose_i); 
		jacobian_J_transpose(dJ_dpose_j, dJt_dpose_j);

		Eigen::Matrix<double, 2, 2> JJt = J*Jt; 
		dJJt_dpose_j.setZero(); 
		dJJt_dpose_i.setZero(); 

		for(int t=0; t<4; t++){
			 int i = std::floor((t)/2) ; 
			 int j = t - (i)*2;
			 // cout<<" t: "<<t<<" i: "<<i<<" j: "<<j<<endl; 
			 for (int k=0; k<4; k++){ 
				 dJJt_dpose_i.row(t) = dJJt_dpose_i.row(t) + Jt(k,j)*dJ_dpose_i.row(4*i+k) + J(i,k)*dJt_dpose_i.row(j+k*2); 
				 dJJt_dpose_j.row(t) = dJJt_dpose_j.row(t) + Jt(k,j)*dJ_dpose_j.row(4*i+k) + J(i,k)*dJt_dpose_j.row(j+k*2); 
			}
		}
	}

	void jacobian_inv_JJt(double const *const *parameters, const Vector3d& pts_i, const Vector3d& pts_j, Eigen::Matrix<double, 4, 6>& diJJt_dpose_i, Eigen::Matrix<double, 4, 6>& diJJt_dpose_j)
	{
		Eigen::Matrix<double, 2, 4> J = J_func(parameters, pts_i, pts_j); 
		Eigen::Matrix<double, 4, 2> Jt = J.transpose(); 
	
		Eigen::Matrix<double, 2, 2> JJt = J*Jt; 
		Eigen::Matrix<double, 2, 2> iJJt = JJt.inverse(); 

		diJJt_dpose_i.setZero(); 
		diJJt_dpose_j.setZero(); 

		Eigen::Matrix<double, 4, 6> dJJt_dpose_i, dJJt_dpose_j; 
    	jacobian_JJt(parameters, pts_i, pts_j, dJJt_dpose_i, dJJt_dpose_j); 

		for(int c=0; c<6; c++){
			{
				Eigen::Matrix<double, 4, 1> dA_dx = dJJt_dpose_i.col(c); 
				Eigen::Map<Matrix<double, 2, 2, Eigen::RowMajor> > dA_dc(dA_dx.data()); 
				Matrix<double, 2, 2> dy_dc = - iJJt * dA_dc * iJJt; 
				Eigen::Map<Matrix<double, 4, 1> > dy_dc_i(dy_dc.data()); 
				diJJt_dpose_i.col(c) = dy_dc_i; 
			}
			{
				Eigen::Matrix<double, 4, 1> dA_dx = dJJt_dpose_j.col(c); 
				Eigen::Map<Matrix<double, 2, 2, Eigen::RowMajor> > dA_dc(dA_dx.data()); 
				Matrix<double, 2, 2> dy_dc = - iJJt * dA_dc * iJJt; 
				Eigen::Map<Matrix<double, 4, 1> > dy_dc_j(dy_dc.data()); 
				diJJt_dpose_j.col(c) = dy_dc_j; 
			}
		}
		return ; 
	}

	void jacobian_algebra_error(double const *const *parameters, const Vector3d& pts_i, const Vector3d& pts_j, Eigen::Matrix<double, 2, 6>& de_dpose_i, Eigen::Matrix<double, 2, 6>& de_dpose_j)
	{
		Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
	    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

	    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
	    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

	    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
	    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

	    double inv_dep_i = parameters[3][0];
		

		Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    	Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    	Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    	Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);

		Eigen::Matrix3d Rj = Qj.toRotationMatrix(); 
		Eigen::Matrix3d Ri = Qi.toRotationMatrix(); 
		Eigen::Matrix3d Ric = qic.toRotationMatrix(); 

		Eigen::Matrix3d de_dpj = -Utility::skewSymmetric(pts_j); 
		Eigen::Matrix3d dpj_dti = Ric.transpose() * Rj.transpose(); 
		Eigen::Matrix3d dpj_dtheta_i = -Ric.transpose()*Rj.transpose() *Ri*Utility::skewSymmetric(pts_imu_i); 
		Eigen::Matrix3d dpj_dtj = -Ric.transpose() * Rj.transpose(); 
		Eigen::Matrix3d dpj_dtheta_j = Ric.transpose()*Utility::skewSymmetric(pts_imu_j);

		Eigen::Matrix<double, 2, 3> ext_m; 
		ext_m << 1, 0, 0, 
				 0, 1, 0;
		Eigen::Matrix<double, 3, 6> dpj_dpose_i, dpj_dpose_j; 
		dpj_dpose_i.block<3,3>(0,0) = dpj_dti; dpj_dpose_i.block<3,3>(0,3) = dpj_dtheta_i; 
		dpj_dpose_j.block<3,3>(0,0) = dpj_dtj; dpj_dpose_j.block<3,3>(0,3) = dpj_dtheta_j;
		de_dpose_i = ext_m * de_dpj * dpj_dpose_i; 
		de_dpose_j = ext_m * de_dpj * dpj_dpose_j; 
		return ;
	}

}

void SampsonFactorCross::compute_Jacobian_pose_i(double const *const *parameters, double** jacobians) const
{
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    double inv_dep_i = parameters[3][0];

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    
    double dep_j = pts_camera_j.z();
    // Eigen::Vector2d epsilon = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    // using cross product as epsilon 
    Eigen::Vector2d epsilon(pts_camera_j.y() - dep_j*pts_j.y(), -pts_camera_j.x() + dep_j * pts_j.x()); 

    // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
    Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
    // d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
    //                        0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
    d_epsilon_d_pts_cam_j << 0., 1., -pts_j.y(), 
                             -1, 0, pts_j.x(); 

    Eigen::Matrix<double, 3, 3> d_pts_cam_j_d_pts_cam_i; 
    d_pts_cam_j_d_pts_cam_i = qic.inverse()*Qj.inverse()*Qi*qic; 
    Eigen::Matrix<double, 3, 2> d_pts_cam_i_d_pts_i; 
    d_pts_cam_i_d_pts_i << 1./inv_dep_i, 0, 
                           0, 1./inv_dep_i,
                           0, 0;

    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
    Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_j = Eigen::Matrix<double,2,2>::Identity()*-1; 
    d_epsilon_d_pts_j << 0, -dep_j,
                        dep_j, 0; 

    Eigen::Matrix<double, 2, 4> d_epsilon_d_X; 
    d_epsilon_d_X.block<2,2>(0, 0) = d_epsilon_d_pts_i; 
    d_epsilon_d_X.block<2,2>(0, 2) = d_epsilon_d_pts_j;

    // sampson approximation 
    // Eigen::Map<Eigen::Matrix<double, 4, 1>> residual(residuals); 
    Eigen::Matrix<double, 2, 4> J = d_epsilon_d_X; 
    Eigen::Matrix<double, 4, 2> Jt = J.transpose();
    Eigen::Matrix<double, 2, 2> JtJ = J*Jt; 
    Eigen::Matrix<double, 2, 2> A = JtJ.inverse();
    Eigen::Matrix<double, 2, 1> ze = A * epsilon; 
    // Eigen::Matrix<double, 4, 2> JJ = J.transpose()*JtJ.inverse(); 
    // residual = -JJ*epsilon;
    // residual = sqrt_info * residual; 

    // Eigen::Matrix<double, 8, 6> dJ_dpose_i, dJ_dpose_j; 
    // jacobian_J(parameters, pts_i, pts_j, dJ_dpose_i, dJ_dpose_j); 
    // cout<<"dJ_dpose_i: "<<endl<<dJ_dpose_i<<endl; 
    // cout<<"dJ_dpose_j: "<<endl<<dJ_dpose_j<<endl;

    Eigen::Matrix<double, 8, 6> dJt_dpose_i, dJt_dpose_j;
    jacobian_Jt(parameters, pts_i, pts_j, dJt_dpose_i, dJt_dpose_j);

    // cout << "dJt_dpose_i: "<<endl<<dJt_dpose_i<<endl; 
    // cout << "dJt_dpose_j: "<<endl<<dJt_dpose_j<<endl;

    // Eigen::Matrix<double, 4, 6> dJJt_dpose_i, dJJt_dpose_j; 
    // jacobian_JJt(parameters, pts_i, pts_j, dJJt_dpose_i, dJJt_dpose_j); 
    // cout << "dJJt_dpose_i: "<<endl<<dJJt_dpose_i<<endl; 
    // cout << "dJJt_dpose_j: "<<endl<<dJJt_dpose_j<<endl; 

  	Eigen::Matrix<double, 4, 6> diJJt_dpose_i, diJJt_dpose_j; 
    jacobian_inv_JJt(parameters, pts_i, pts_j, diJJt_dpose_i, diJJt_dpose_j); 
    // cout << "diJJt_dpose_i: "<<endl<<diJJt_dpose_i<<endl; 
    // cout << "diJJt_dpose_j: "<<endl<<diJJt_dpose_j<<endl; 

  	Eigen::Matrix<double, 2, 6> de_dpose_i, de_dpose_j; 
    jacobian_algebra_error(parameters, pts_i, pts_j, de_dpose_i, de_dpose_j); 
    // cout << "de_dpose_i: "<<endl<<de_dpose_i<<endl; 
    // cout << "de_dpose_j: "<<endl<<de_dpose_j<<endl; 

    // dz_dpose_i, z = inv(JJ')*epsilon 
    Eigen::Matrix<double, 2, 6> dz_dpose_i = Eigen::Matrix<double, 2, 6>::Zero(); 
    Eigen::Matrix<double, 2, 6> dz_dpose_j = Eigen::Matrix<double, 2, 6>::Zero(); 

    for(int i=0; i<2; i++){
    	for(int j=0; j<2; j++){
    		dz_dpose_i.row(i) = dz_dpose_i.row(i) + epsilon(j)*diJJt_dpose_i.row(2*i+j) + A(i,j)*de_dpose_i.row(j); 
    		dz_dpose_j.row(i) = dz_dpose_j.row(i) + epsilon(j)*diJJt_dpose_j.row(2*i+j) + A(i,j)*de_dpose_j.row(j);
    	}
    }
    // cout << "dz_dpose_i: "<<endl<<dz_dpose_i<<endl; 
    // cout << "dz_dpose_j: "<<endl<<dz_dpose_j<<endl; 

      // dz_dpose_i, z = inv(JJ')*epsilon 
    Eigen::Matrix<double, 4, 6> dse_dpose_i = Eigen::Matrix<double, 4, 6>::Zero(); 
    Eigen::Matrix<double, 4, 6> dse_dpose_j = Eigen::Matrix<double, 4, 6>::Zero(); 

    for(int i=0; i<4; i++){
    	for(int j=0; j<2; j++){
    		dse_dpose_i.row(i) = dse_dpose_i.row(i) + ze(j)*dJt_dpose_i.row(2*i+j) + Jt(i,j) * dz_dpose_i.row(j); 
    		dse_dpose_j.row(i) = dse_dpose_j.row(i) + ze(j)*dJt_dpose_j.row(2*i+j) + Jt(i,j) * dz_dpose_j.row(j);
    	}
    }
    cout << "dse_dpose_i: "<<endl<<dse_dpose_i<<endl; 
    cout << "dse_dpose_j: "<<endl<<dse_dpose_j<<endl; 

    Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]); 
    Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]); 
    jacobian_pose_i.setZero();
    jacobian_pose_j.setZero(); 
    jacobian_pose_i.block<4,6>(0,0) = -sqrt_info*dse_dpose_i;
    jacobian_pose_j.block<4,6>(0,0) = -sqrt_info*dse_dpose_j;


    return ; 
}


void SampsonFactorCross::compute_Jacobian_pose_j(double const *const *parameters, double** jacobians) const
{



}
