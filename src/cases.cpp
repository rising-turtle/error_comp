/*

	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 


	cases to be tested 

*/

#include "cases.h"
#include <cstdlib>
#include <random>
#include "Eigen/Core"
#include "Eigen/Dense"

Case::Case():
m_thre_cnt_feats(100){}

Case::~Case(){}

void Case::init_pose_features(int N_feats)
{
	cout<<"cases.cpp: in Case init_pose_features(), should be implemented in sub classes"<<endl; 
}

void Case::gen_observations()
{
	cout<<"cases.cpp: in Case gen_observations(), should be implemented in sub classes"<<endl; 
}

void Case::add_noise(double pix_std)
{

	double focal_length = 400; 

	static std::random_device rd{};
    static std::mt19937 gen{rd()};
 
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> noise{0,pix_std};

    std::uniform_real_distribution<> rangle(0. , 360.);

    double dis, angle; 

	for(int i=0; i<mv_obs.size(); i++){
		for(int j=0; j<mv_obs[i].size(); j++){
			FeatureMeasurement& fm = mv_obs[i][j]; 

			dis = noise(gen); 
			angle = rangle(gen)*M_PI/180.; 

			// add noise for pixel location 
			fm.xy[0] += (dis*cos(angle)); // noise(gen); 
			fm.xy[1] += (dis*sin(angle)); // noise(gen); 

			// add noise to depth measurement 
			// depth_std = 1.45e-3*(depth)^2 
			double dpt_std = 1.45*0.001*fm.depth*fm.depth; 
			std::normal_distribution<> depth_noise{0,dpt_std};
			fm.depth += depth_noise(gen); 
		}
	}
}

double Case::rmse(double p[][7], double* perr){
	
	double sum_err_pos = 0; 
	double sum_err_ori = 0; 

	for(int i=1; i<mv_gt_poses.size(); i++){
		Eigen::Vector3d gt_x(mv_gt_poses[i][0], mv_gt_poses[i][1], mv_gt_poses[i][2]); 
		Eigen::Quaterniond gt_q(mv_gt_poses[i][6], mv_gt_poses[i][3], mv_gt_poses[i][4], mv_gt_poses[i][5]);
	
		Eigen::Vector3d x(p[i][0], p[i][1], p[i][2]); 
		Eigen::Quaterniond q(p[i][6], p[i][3], p[i][4], p[i][5]); 

		Eigen::Vector3d dx = gt_x - x; 
		sum_err_pos += dx.dot(dx); 

		Eigen::Quaterniond dq = gt_q.inverse()*q; 
		Eigen::Matrix<double, 3, 3> R = dq.toRotationMatrix(); 
		double theta = acos(((R(0,0)+R(1,1) + R(2,2))-1)/2);
		sum_err_ori += theta*theta;
	}

	if(mv_gt_poses.size() > 1){
		sum_err_pos /= (mv_gt_poses.size()-1);
		sum_err_ori /= (mv_gt_poses.size()-1);
		sum_err_pos = sqrt(sum_err_pos);
		sum_err_ori = sqrt(sum_err_ori); 
	}

	perr[0] = sum_err_pos; 
	perr[1] = sum_err_ori; 

	return sum_err_pos;

}

Case_forward::Case_forward(){}
Case_forward::~Case_forward(){}

void Case_forward::init_pose_features(int N_feats)
{

	// features in front of the camera lie between x ~ [-2,2] y ~[-2, 2] z ~ [2, 10]
	// int N_feats = 300; 

	m_thre_cnt_feats = N_feats; 
	int cnt_gen_feats = N_feats*3; 
	if(cnt_gen_feats < 200) cnt_gen_feats = 200;

	mv_feats.resize(cnt_gen_feats); // usually generate more features  

	double x, y, z; 

	std::random_device rd; 
	std::mt19937 gen(rd()); 
 	std::uniform_real_distribution<> xy_dis(-2.0, 2.0);
 	std::uniform_real_distribution<> z_dis(2.0, 7.0); // 2-10

	for(int i=0; i<cnt_gen_feats; i++){
		Feature& ft = mv_feats[i]; 
		ft.m_id = i+1; 
		ft.m_loc[0] = xy_dis(gen); 
		ft.m_loc[1] = xy_dis(gen); 
		ft.m_loc[2] = z_dis(gen);
	}

	// generate poses, move forward one step  

	vector<double> p1{0, 0, 0, 0, 0, 0, 1.};
	vector<double> p2{0, 0, 1.0, 0, 0, 0, 1.}; 

	mv_gt_poses.push_back(p1); 
	mv_gt_poses.push_back(p2);  

	mv_poses = mv_gt_poses; 

}


void Case_forward::gen_observations()
{
	// mv_poses = mv_gt_poses;
	 // mv_poses[1][0] += 0.1; 
	 // mv_poses[1][1] += 0.1; 
	 // mv_poses[1][2] -= 0.1;

	for(int j=0; j<mv_gt_poses.size(); j++){

		Eigen::Vector3d Pj(mv_gt_poses[j][0], mv_gt_poses[j][1], mv_gt_poses[j][2]);
    	Eigen::Quaterniond Qj(mv_gt_poses[j][6], mv_gt_poses[j][3], mv_gt_poses[j][4], mv_gt_poses[j][5]);

    	vector<FeatureMeasurement> obs; 

    	for(int i=0; i<mv_feats.size(); i++){

    		Feature& fi = mv_feats[i]; 

    		Eigen::Vector3d pts_w(fi.m_loc[0], fi.m_loc[1], fi.m_loc[2]); 
    		Eigen::Vector3d pts_j = Qj.inverse()*(pts_w - Pj); 

    		if(pts_j.z() > 0){
    			double px = pts_j.x()/pts_j.z(); 
    			double py = pts_j.y()/pts_j.z(); 
    			// cout <<"feat i = "<<i+1<<" "<<fi.m_loc[0]<<" "<<fi.m_loc[1]<<" "<<fi.m_loc[2]<<" project to "<<px<<" "<<py<<endl; 

    			// simulate a camera frame [640, 480] with f = 400 
    			if((fabs(px) < (320./400)) && (fabs(py) < (240./400)) ){ // visible 
    				FeatureMeasurement m; 
    				m.feat_id = fi.m_id; 
    				m.pose_id = j; 
    				m.feat_index = obs.size();
    				m.xy[0] = px; m.xy[1] = py; 
    				m.depth = pts_j.z(); 

    				if(fi.mv_pose_id.size() == 0){
    					fi.m_ref_id = 0; // use the first observation 
    				}

    				fi.mv_pose_id.push_back(j); 
    				fi.mv_feat_idx.push_back(m.feat_index); 
    				obs.push_back(m);
    				// cout <<"feat i = "<<i+1<<" is visible at pose = "<<j<<endl;
    			}
    		}
    	}
    	mv_obs.push_back(obs);

	}

}
