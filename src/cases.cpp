/*

	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 


	cases to be tested 

*/

#include "case.h"
#include <cstdlib>
#include <random>

Case::Case(){}

Case::~Case(){}

void Case::init_pose_features()
{
	cout<<"cases.cpp: in Case init_pose_features(), should be implemented in sub classes"<<endl; 
}

void Case::gen_observations()
{
	cout<<"cases.cpp: in Case gen_observations(), should be implemented in sub classes"<<endl; 
}

Case_forward::Case_forward(){}
Case_forward::~Case_forward(){}

void Case_forward::init_pose_features()
{

	// features in front of the camera lie between x ~ [-2,2] y ~[-2, 2] z ~ [2, 10]
	int N_feats = 200; 

	mv_feats.resize(N_feats); 

	double x, y, z; 

	std::random_device rd; 
	std::mt19937 gen(rd()); 
	std::uniform_int_distribution<> xy_dis(-200, 200); 
	std::uniform_int_distribution<> z_dis(200, 1000); 

	for(int i=0; i<N_feats; i++){
		Feature& ft = mv_feats[i]; 
		ft.m_id = i+1; 
		ft.loc[0] = gen(xy_dis)*0.01; 
		ft.loc[1] = gen(xy_dis)*0.01; 
		ft.loc[2] = gen(z_dis)*0.01;
	}

	// generate poses, move forward one step  

	vector<double> p1{0, 0, 0, 0, 0, 0, 1.};
	vector<double> p2{0, 0, 0.5, 0, 0, 0, 1.}; 

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

    		Eigen::Vector3d pts_w(fi.loc[0], fi.loc[1], fi.loc[2]); 
    		Eigen::Vector3d pts_j = Qj.inverse()*(pts_w - Pj); 

    		if(pts_j.z() > 0){
    			double px = pts_j.x()/pts_j.z(); 
    			double py = pts_j.y()/pts_j.z(); 

    			// simulate a camera frame [640, 480] with f = 400 
    			if(fabs(px) < (320/400) && fabs(py) < (240/400) ){ // visible 
    				FeatureMeasurement m; 
    				m.feat_id = fi.m_id; 
    				m.pose_id = j; 
    				m.feat_index = obs.size(); 
    				obs.push_back(m); 
    				m.xy[0] = px; m.xy[1] = py; 
    				m.depth = pts_j.z(); 

    				if(fi.mv_pose_id.size() == 0){
    					fi.m_ref_id = 0; // use the first observation 
    				}

    				fi.mv_pose_id.push_back(j); 
    				fi.mv_feat_idx.push_back(m.feat_index); 
    			}
    		}
    	}

	}

}
