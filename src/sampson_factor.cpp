/*
	Aug. 21, 2019, He Zhang, hzhang8@vcu.edu 

	sampson approximation to the geometric distance

*/

#include "sampson_factor.h"
#include <iomanip>      // std::setprecision

using namespace std;


Eigen::Matrix<double, 4,4> SampsonFactor::sqrt_info;

SampsonFactor::SampsonFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j)
// const SampsonCostFunctor* functor)
 : pts_i(_pts_i), pts_j(_pts_j)//, functor_(functor)
{
	SampsonFactor::sqrt_info = 24* Eigen::Matrix<double, 4, 4>::Identity(); 
}


bool SampsonFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
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

    // sampson approximation 
    Eigen::Map<Eigen::Matrix<double, 4, 1>> residual(residuals); 
    Eigen::Matrix<double, 2, 4> J = d_epsilon_d_X; 
    Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
    Eigen::Matrix<double, 4, 2> JJ = J.transpose()*JtJ.inverse(); 
    residual = -JJ*epsilon;
 
    residual = sqrt_info * residual; 

    // cout<<" in factor: epsilon: "<<endl << epsilon<<endl; 
    // cout <<" JJ = "<<endl<<JJ<<endl;

    double eps = 1e-6; 

    if(jacobians){

    	Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();
        Eigen::Matrix<double, 2, 3> reduced_epsilon_d_pts_cam_j(2, 3);

        // reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),
        //    0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j);

        // reduce = sqrt_info * reduce;

        Eigen::Matrix<double, 2, 4> J_new; 
        Eigen::Matrix<double, 2, 2> JtJ_new;  
    	Eigen::Matrix<double, 4, 2> JJ_new; 
    	Eigen::Matrix<double, 2,1> epsilon_new; 
    	Eigen::Matrix<double, 4,1> residual_new; 

    	if(jacobians[0]){

    		Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]); 

    		Eigen::Matrix<double, 3, 6> jaco_i; // d_epsilon_d_pose_i
            jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();
            jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);

            // d_res_d_pi
            for(int k=0; k<3; k++){
            	Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
            	Eigen::Vector3d tmp_Pi = Pi + delta; 

    			Eigen::Vector3d pts_w = Qi * pts_imu_i + tmp_Pi;
    			Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    			Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

            	// epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
            	double dep_j = pts_camera_j.z();
    			epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    			// compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
    			Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
    			d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
    						0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
    			Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
            	J_new.block<2,2>(0,0) = d_epsilon_d_pts_i;
            	J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 

            	JtJ_new = J_new*J_new.transpose(); 
   				JJ_new = J_new.transpose()*JtJ_new.inverse(); 
    			residual_new = -sqrt_info*JJ_new*epsilon_new;
    			
            	// jacobian_pose_i.leftCols<3>() = -sqrt_info * JJ * reduce * jaco_i.leftCols<3>();
            	jacobian_pose_i.col(k) = (residual_new - residual)/eps; 
        	}

            // d_res_d_qi 
            for(int k=0; k<3; k++){

            	Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
    			Eigen::Quaterniond tmp_Qi = Qi * Utility::deltaQ(delta);

    			Eigen::Vector3d pts_w = tmp_Qi * pts_imu_i + Pi;
    			Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    			Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

            	// epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
            	double dep_j = pts_camera_j.z();
    			epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    			// compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
    			Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
    			d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
    						0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
   				Eigen::Matrix<double, 3, 3> tmp_m = qic.inverse()*Qj.inverse()*tmp_Qi.toRotationMatrix()*qic.toRotationMatrix();
   				Eigen::Matrix<double, 2, 2> tmp_J = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_pts_i;

   				J_new.block<2,2>(0,0) = tmp_J; 
   				J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 
   				JtJ_new = J_new*J_new.transpose(); 
   				JJ_new = J_new.transpose()*JtJ_new.inverse(); 

   				// Eigen::Matrix<double, 2, 1> tt = reduce * jaco_i.rightCols<3>() * (delta); 
   				// epsilon_new = epsilon + tt; // jaco_i.rightCols<3>() * (delta);

   				residual_new = -sqrt_info*JJ_new*epsilon_new;

   				jacobian_pose_i.col(k+3) = (residual_new - residual)/eps; 
            }
            jacobian_pose_i.rightCols<1>().setZero();

    	}
    	if(jacobians[1]){

    		Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]); 

    		Eigen::Matrix<double, 3, 6> jaco_j; // d_epsilon_d_pose_j
            jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
            jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);

            // d_res_d_pj
            // jacobian_pose_j.leftCols<3>() = -sqrt_info * JJ * reduce * jaco_j.leftCols<3>();

            // d_res_d_pj
            for(int k=0; k<3; k++){
            	Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
            	Eigen::Vector3d tmp_Pj = Pj + delta; 

    			Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - tmp_Pj);
    			Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

            	// epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
            	double dep_j = pts_camera_j.z();
    			epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

    			// compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
    			Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
    			d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
    						0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
    			Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
            	J_new.block<2,2>(0,0) = d_epsilon_d_pts_i;
            	J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 

            	JtJ_new = J_new*J_new.transpose(); 
   				JJ_new = J_new.transpose()*JtJ_new.inverse(); 
    			residual_new = -sqrt_info*JJ_new*epsilon_new;
    			
            	// jacobian_pose_i.leftCols<3>() = -sqrt_info * JJ * reduce * jaco_i.leftCols<3>();
            	jacobian_pose_j.col(k) = (residual_new - residual)/eps; 
        	}


            // d_res_d_qj
            for(int k=0; k<3; k++){

            	Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
    			Eigen::Quaterniond tmp_Qj = Qj * Utility::deltaQ(delta); 

    			Eigen::Vector3d pts_imu_j = tmp_Qj.inverse() * (pts_w - Pj);
    			Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

            	// epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
            	double dep_j = pts_camera_j.z();
    			epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
				// compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
    			Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
    			d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
    						0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
   				Eigen::Matrix<double, 3, 3> tmp_m = qic.inverse()*tmp_Qj.inverse()*Qi.toRotationMatrix()*qic.toRotationMatrix();
   				Eigen::Matrix<double, 2, 2> tmp_J = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_pts_i;

   				J_new.block<2,2>(0,0) = tmp_J; 
   				J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 
   				JtJ_new = J_new*J_new.transpose(); 
   				JJ_new = J_new.transpose()*JtJ_new.inverse(); 

   				// Eigen::Matrix<double, 2, 1> tt = reduce * jaco_j.rightCols<3>() * (delta); 
   				// epsilon_new = epsilon + tt;

   				residual_new = -sqrt_info*JJ_new*epsilon_new;

   				jacobian_pose_j.col(k+3) = (residual_new - residual)/eps; 
            }
            jacobian_pose_j.rightCols<1>().setZero();
    	}
    	if(jacobians[2]){
    		//TODO: d_res_d_pose_ic 

    	}
    	if(jacobians[3]){
    		// TODO: d_res_d_lambda
    	}

    }
    return true; 
}


void SampsonFactor::check(double ** parameters)
{
	double *res = new double[4];
    double **jaco = new double *[2];
    jaco[0] = new double[4 * 7];
    jaco[1] = new double[4 * 7];
    cout.precision(8);
    // jaco[2] = new double[2 * 7];
    // jaco[3] = new double[2 * 1];
    Evaluate(parameters, res, jaco);
    puts("check begins");

    puts("my");

    Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> J_pi(jaco[0]); 
    Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> J_pj(jaco[1]); 
    

    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;

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

    // sampson approximation 
    // Eigen::Map<Eigen::Matrix<double, 4, 1>> residual(residuals); 

    Eigen::Matrix<double, 4, 1> residual;
    Eigen::Matrix<double, 2, 4> J = d_epsilon_d_X; 
    Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
    Eigen::Matrix<double, 4, 2> JJ = J.transpose()*JtJ.inverse(); 
    residual = -JJ*epsilon;
 
    residual = sqrt_info * residual;
    Eigen::Matrix<double, 4, 12> num_jacobian;

    double eps = 1e-6; 
    Eigen::Matrix<double, 2, 3> reduce(2, 3);
    reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),
        0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j);

    Eigen::Matrix3d Ri = Qi.toRotationMatrix();
    Eigen::Matrix3d Rj = Qj.toRotationMatrix();
    Eigen::Matrix3d ric = qic.toRotationMatrix();

 	Eigen::Matrix<double, 3, 6> jaco_i; // d_epsilon_d_pose_i
        jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();
        jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);

	Eigen::Matrix<double, 3, 6> jaco_j; // d_epsilon_d_pose_j
        jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
        jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);

    
//    cout <<"recompute epsilon: "<<endl<<epsilon<<endl;
//	cout <<" JJ = "<<endl << JJ<<endl; 

/*
	{
    	Eigen::Vector3d delta = Eigen::Vector3d(0, 0, 1) * eps; 
    	// Eigen::Vector3d Pii = Pi + delta;
    	Pi += delta; 

    	Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
	    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
	    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
	    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
	    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

	    double dep_j = pts_camera_j.z();
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

	    // sampson approximation 
	   	Eigen::Matrix<double, 4, 1> residual_new; 
	    Eigen::Matrix<double, 2, 4> J = d_epsilon_d_X; 

	    Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
	    Eigen::Matrix<double, 4, 2> JJ = J.transpose()*JtJ.inverse(); 
	    residual_new = -JJ*epsilon;
	    residual_new = sqrt_info * residual_new;

	    cout <<"epsilon: "<<endl<<epsilon<<endl;
	    cout <<" JJ = "<<endl << JJ<<endl; 

	}
*/

    for(int k=0; k<12; k++){

    	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    	Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    	Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    	Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    	int a = k/3; int b = k%3; 
    	Eigen::Vector3d delta = Eigen::Vector3d(b==0, b==1, b==2) * eps; 

    	if(a==0) {
    		Pi += delta; 
    		// check epsilon
    		// print expected J

    		// Eigen::Matrix<double, 4, 3> tmp_J_pi = -sqrt_info * JJ * reduce * jaco_i.leftCols<3>(); 
    		// Eigen::Matrix<double, 2, 1> tmp_epsilon = reduce*jaco_i.leftCols<3>() * delta;
    		// tmp_epsilon = epsilon + tmp_epsilon;
    		// Eigen::Vector4d tt = tmp_J_pi * delta; 
    		// Eigen::Vector4d tmp_e1 = residual + tt;
    		// tmp_e1 = sqrt_info * tmp_e1; 
    		// if(k == 2)
    		// {	
    		// 	cout <<"k = "<<k<<" expected residual: "<<endl<<tmp_e1.transpose()<<endl;
    		// 	cout <<" residual = "<<endl<<residual.transpose()<<endl;
    		// 	cout <<" tmp_J_pi = "<<endl<<tmp_J_pi<<endl; 

    		// 	cout <<" JJ = "<<endl<<JJ<<endl; 
	    	// 	cout <<" epsilon: "<<endl<<tmp_epsilon<<endl;
    		// 	cout <<" J_pi = "<<endl<<J_pi.block<4,3>(0,0)<<endl;
    		// }
    	}
    	else if(a == 1)
    		Qi = Qi * Utility::deltaQ(delta); 
    	else if(a == 2){
    		Pj += delta;
			// print expected J
			// Eigen::Matrix<double, 4, 3> tmp_J_pj = -sqrt_info * JJ * reduce * jaco_j.leftCols<3>(); 
   //  		Eigen::Vector4d tt = tmp_J_pj * delta;
   //  		Eigen::Vector4d tmp_e1 = residual + tt; 
   //  		if(k == 2){
   //  			cout << "a = "<<a<<" expected: "<<endl<<tmp_e1.transpose()<<endl;
   //  			cout <<" tmp_J_pj = "<<endl<<tmp_J_pj<<endl; 
   //  			cout <<" J_pj = "<<endl<<J_pj.block<4,3>(0,0)<<endl;
   //  		}
    	}
    	else if(a == 3)
    		Qj = Qj * Utility::deltaQ(delta); 

    	Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
	    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
	    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
	    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
	    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

	    double dep_j = pts_camera_j.z();
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

	    // sampson approximation 
	   	Eigen::Matrix<double, 4, 1> residual_new; 
	    Eigen::Matrix<double, 2, 4> J = d_epsilon_d_X; 

	    Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
	    Eigen::Matrix<double, 4, 2> JJ = J.transpose()*JtJ.inverse(); 
	    residual_new = -JJ*epsilon;
	    residual_new = sqrt_info * residual_new;

	    num_jacobian.col(k) = (residual_new - residual)/eps; 

	    //  if(a == 0 || a == 2){
	   	// 	if(k == 2){
	    // 		cout << "k = "<<k<<" residual new: "<<endl<<residual_new.transpose()<<endl;
	    // 		cout <<" residual = "<<endl<<residual.transpose()<<endl;
	    // 		cout <<" JJ = "<<endl<<JJ<<endl; 
	    // 		cout <<" epsilon: "<<endl<<epsilon<<endl;
	    // 		cout <<" num_jacobian(2) = "<<endl<<num_jacobian.col(k)<<endl;
	   	// 	}
	    // }

    }

    std::cout<<"sampson_factor.cpp: in check num_jacobian: "<<std::endl; 

    std::cout <<num_jacobian<<std::endl;

}


////////////////////////////////////////////////////////////////////////

// below is sampson factor with lambda  

////////////////////////////////////////////////////////////////////////


Eigen::Matrix<double, 5,5> SampsonFactorWithLambda::sqrt_info;

SampsonFactorWithLambda::SampsonFactorWithLambda(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j)
// const SampsonCostFunctor* functor)
 : pts_i(_pts_i), pts_j(_pts_j)//, functor_(functor)
{
    SampsonFactorWithLambda::sqrt_info = 24*Eigen::Matrix<double, 5, 5>::Identity(); 
}


bool SampsonFactorWithLambda::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
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

    // sampson approximation 
    Eigen::Map<Eigen::Matrix<double, 5, 1>> residual(residuals); 
    Eigen::Matrix<double, 2, 5> J = d_epsilon_d_X; 
    Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
    Eigen::Matrix<double, 5, 2> JJ = J.transpose()*JtJ.inverse(); 
    residual = -JJ*epsilon;
 
    residual = sqrt_info * residual; 

    // cout<<" in factor: epsilon: "<<endl << epsilon<<endl; 
    // cout <<" JJ = "<<endl<<JJ<<endl;

    double eps = 1e-6; 

    if(jacobians){

        Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();

        Eigen::Matrix<double, 2, 5> J_new; 
        Eigen::Matrix<double, 2, 2> JtJ_new;  
        Eigen::Matrix<double, 5, 2> JJ_new; 
        Eigen::Matrix<double, 2,1> epsilon_new; 
        Eigen::Matrix<double, 5,1> residual_new; 

        if(jacobians[0]){

            Eigen::Map<Eigen::Matrix<double, 5, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]); 

            Eigen::Matrix<double, 3, 6> jaco_i; // d_epsilon_d_pose_i
            jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();
            jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);

            // d_res_d_pi
            for(int k=0; k<3; k++){
                Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
                Eigen::Vector3d tmp_Pi = Pi + delta; 

                Eigen::Vector3d pts_w = Qi * pts_imu_i + tmp_Pi;
                Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
                Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

                // epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
                double dep_j = pts_camera_j.z();
                epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

                // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
                Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
                d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
                            0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
                Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
                Eigen::Matrix<double, 2, 1> d_epsilon_d_lambda = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_lambda; 
                J_new.block<2,2>(0,0) = d_epsilon_d_pts_i;
                J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 
                J_new.block<2,1>(0,4) = d_epsilon_d_lambda;

                JtJ_new = J_new*J_new.transpose(); 
                JJ_new = J_new.transpose()*JtJ_new.inverse(); 
                residual_new = -sqrt_info*JJ_new*epsilon_new;
                
                // jacobian_pose_i.leftCols<3>() = -sqrt_info * JJ * reduce * jaco_i.leftCols<3>();
                jacobian_pose_i.col(k) = (residual_new - residual)/eps; 
            }

            // d_res_d_qi 
            for(int k=0; k<3; k++){

                Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
                Eigen::Quaterniond tmp_Qi = Qi * Utility::deltaQ(delta);

                Eigen::Vector3d pts_w = tmp_Qi * pts_imu_i + Pi;
                Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
                Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

                // epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
                double dep_j = pts_camera_j.z();
                epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

                // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
                Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
                d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
                            0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
                Eigen::Matrix<double, 3, 3> tmp_m = qic.inverse()*Qj.inverse()*tmp_Qi.toRotationMatrix()*qic.toRotationMatrix();
                Eigen::Matrix<double, 2, 2> tmp_J = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_pts_i;
                Eigen::Matrix<double, 2, 1> d_epsilon_d_lambda = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_lambda; 

                J_new.block<2,2>(0,0) = tmp_J; 
                J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 
                J_new.block<2,1>(0,4) = d_epsilon_d_lambda;
                JtJ_new = J_new*J_new.transpose(); 
                JJ_new = J_new.transpose()*JtJ_new.inverse(); 

                // Eigen::Matrix<double, 2, 1> tt = reduce * jaco_i.rightCols<3>() * (delta); 
                // epsilon_new = epsilon + tt; // jaco_i.rightCols<3>() * (delta);

                residual_new = -sqrt_info*JJ_new*epsilon_new;

                jacobian_pose_i.col(k+3) = (residual_new - residual)/eps; 
            }
            jacobian_pose_i.rightCols<1>().setZero();

        }
        if(jacobians[1]){

            Eigen::Map<Eigen::Matrix<double, 5, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]); 

            Eigen::Matrix<double, 3, 6> jaco_j; // d_epsilon_d_pose_j
            jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
            jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);

            // d_res_d_pi
            for(int k=0; k<3; k++){
                Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
                Eigen::Vector3d tmp_Pj = Pj + delta; 

                Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - tmp_Pj);
                Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

                // epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
                double dep_j = pts_camera_j.z();
                epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();

                // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
                Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
                d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
                            0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
                Eigen::Matrix<double, 2, 2> d_epsilon_d_pts_i = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_pts_i; 
                Eigen::Matrix<double, 2, 1> d_epsilon_d_lambda = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_lambda; 

                J_new.block<2,2>(0,0) = d_epsilon_d_pts_i;
                J_new.block<2,2>(0,2) = d_epsilon_d_pts_j; 
                J_new.block<2,1>(0,4) = d_epsilon_d_lambda;

                JtJ_new = J_new*J_new.transpose(); 
                JJ_new = J_new.transpose()*JtJ_new.inverse(); 
                residual_new = -sqrt_info*JJ_new*epsilon_new;
                
                // jacobian_pose_i.leftCols<3>() = -sqrt_info * JJ * reduce * jaco_i.leftCols<3>();
                jacobian_pose_j.col(k) = (residual_new - residual)/eps; 
            }


            // d_res_d_qi 
            for(int k=0; k<3; k++){

                Eigen::Vector3d delta = Eigen::Vector3d(k==0, k==1, k==2)*eps; 
                Eigen::Quaterniond tmp_Qj = Qj * Utility::deltaQ(delta); 

                Eigen::Vector3d pts_imu_j = tmp_Qj.inverse() * (pts_w - Pj);
                Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

                // epsilon_new = reduce * jaco_i.leftCols<3>() * delta; 
                double dep_j = pts_camera_j.z();
                epsilon_new = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
                // compute J = d_epsilon_d_X, [2x4], where X = (pi_x, pi_y, pj_x, pj_y) 
                Eigen::Matrix<double, 2, 3> d_epsilon_d_pts_cam_j;
                d_epsilon_d_pts_cam_j << 1./dep_j, 0, -pts_camera_j.x()/(dep_j*dep_j), 
                            0, 1./dep_j, -pts_camera_j.y()/(dep_j*dep_j); 
                Eigen::Matrix<double, 3, 3> tmp_m = qic.inverse()*tmp_Qj.inverse()*Qi.toRotationMatrix()*qic.toRotationMatrix();
                Eigen::Matrix<double, 2, 2> tmp_J = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_pts_i;
                Eigen::Matrix<double, 2, 1> d_epsilon_d_lambda = d_epsilon_d_pts_cam_j * tmp_m * d_pts_cam_i_d_lambda; 

                J_new.block<2,2>(0,0) = tmp_J; 
                J_new.block<2,2>(0,2) = d_epsilon_d_pts_j;
                J_new.block<2,1>(0,4) = d_epsilon_d_lambda;
                JtJ_new = J_new*J_new.transpose(); 
                JJ_new = J_new.transpose()*JtJ_new.inverse(); 

                // Eigen::Matrix<double, 2, 1> tt = reduce * jaco_j.rightCols<3>() * (delta); 
                // epsilon_new = epsilon + tt;

                residual_new = -sqrt_info*JJ_new*epsilon_new;

                jacobian_pose_j.col(k+3) = (residual_new - residual)/eps; 
            }
            jacobian_pose_j.rightCols<1>().setZero();
        }
        if(jacobians[2]){
            //TODO: d_res_d_pose_ic 

        }
        if(jacobians[3]){
            // TODO: d_res_d_lambda
        }

    }
    return true; 
}


void SampsonFactorWithLambda::check(double ** parameters)
{
    double *res = new double[5];
    double **jaco = new double *[2];
    jaco[0] = new double[5 * 7];
    jaco[1] = new double[5 * 7];
    cout.precision(8);
    Evaluate(parameters, res, jaco);
    puts("check begins");

    puts("my");

    Eigen::Map<Eigen::Matrix<double, 5, 7, Eigen::RowMajor>> J_pi(jaco[0]); 
    Eigen::Map<Eigen::Matrix<double, 5, 7, Eigen::RowMajor>> J_pj(jaco[1]); 
    

    std::cout << Eigen::Map<Eigen::Matrix<double, 5, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 5, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 5, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;

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

    // sampson approximation 
    // Eigen::Map<Eigen::Matrix<double, 4, 1>> residual(residuals); 

    Eigen::Matrix<double, 5, 1> residual;
    Eigen::Matrix<double, 2, 5> J = d_epsilon_d_X; 
    Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
    Eigen::Matrix<double, 5, 2> JJ = J.transpose()*JtJ.inverse(); 
    residual = -JJ*epsilon;
 
    residual = sqrt_info * residual;
    Eigen::Matrix<double, 5, 12> num_jacobian;

    double eps = 1e-7; 
    Eigen::Matrix3d Ri = Qi.toRotationMatrix();
    Eigen::Matrix3d Rj = Qj.toRotationMatrix();
    Eigen::Matrix3d ric = qic.toRotationMatrix();

    Eigen::Matrix<double, 3, 6> jaco_i; // d_epsilon_d_pose_i
        jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();
        jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);

    Eigen::Matrix<double, 3, 6> jaco_j; // d_epsilon_d_pose_j
        jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
        jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);


    for(int k=0; k<12; k++){

        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        int a = k/3; int b = k%3; 
        Eigen::Vector3d delta = Eigen::Vector3d(b==0, b==1, b==2) * eps; 

        if(a==0) 
            Pi += delta; 
        else if(a == 1)
            Qi = Qi * Utility::deltaQ(delta); 
        else if(a == 2)
            Pj += delta;
        else if(a == 3)
            Qj = Qj * Utility::deltaQ(delta); 

        Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
        Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
        Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
        Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
        Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

        double dep_j = pts_camera_j.z();
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
        Eigen::Matrix<double, 2, 1> d_epsilon_d_lambda = d_epsilon_d_pts_cam_j * d_pts_cam_j_d_pts_cam_i * d_pts_cam_i_d_lambda; 
        Eigen::Matrix<double, 2, 5> d_epsilon_d_X; 
        d_epsilon_d_X.block<2,2>(0, 0) = d_epsilon_d_pts_i; 
        d_epsilon_d_X.block<2,2>(0, 2) = d_epsilon_d_pts_j;
        d_epsilon_d_X.block<2,1>(0, 4) = d_epsilon_d_lambda;

        // sampson approximation 
        Eigen::Matrix<double, 5, 1> residual_new; 
        Eigen::Matrix<double, 2, 5> J = d_epsilon_d_X; 

        Eigen::Matrix<double, 2, 2> JtJ = J*J.transpose(); 
        Eigen::Matrix<double, 5, 2> JJ = J.transpose()*JtJ.inverse(); 
        residual_new = -JJ*epsilon;
        residual_new = sqrt_info * residual_new;

        num_jacobian.col(k) = (residual_new - residual)/eps; 

    }

    std::cout<<"sampson_factor.cpp: in check num_jacobian: "<<std::endl; 

    std::cout <<num_jacobian<<std::endl;

}
