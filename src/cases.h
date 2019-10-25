/*

	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 

	cases to be tested 

*/


#pragma once 

#include <vector>
#include <iostream>

using namespace std; 

class Feature{

public:
	int m_id; // feature id 
	double m_loc[3];  // localation in world coordinate 

	vector<int> mv_pose_id; // has been seen in different frames 
	vector<int> mv_feat_idx; // feature index in the frame at different pose 
	int m_ref_id; // where depth is referenced 
};

class FeatureMeasurement{	

public:
	int feat_id; 
	int pose_id; 	// observed at pose
	int feat_index; // feature index in the pose 
	double xy[2]; // feature coordinate in the normalized plane, x/z, y/z
	// float disparity; // disparity 
	double depth; 

	// to simulate a stereo vision observation 
	double r_xy[2]; // right coordinate 
	static double rig_len; // rig lenth 0.1 meter
};

class Case{

public:
	Case(); 
	virtual ~Case(); 

	vector<vector<double>> mv_poses;  // each pose has seven parameters x,y,z,qx,qy,qz,qw 
	vector<vector<FeatureMeasurement>> mv_obs; // observations at each pose 

	virtual void init_pose_features(int N_feats=300); // init features and poses 
	virtual void gen_observations();   // generate observations 
	virtual void add_noise(double pix_std); // pixel std 
	virtual double rmse(double p[][7], double* perr); // compute rmse 
	vector<Feature> mv_feats; 			// ground truth features 
	vector<vector<double>> mv_gt_poses; // ground truth poses
	int m_thre_cnt_feats; // maximum used features 
};

class Case_forward : public Case
{
public:
	Case_forward(); 
	virtual ~Case_forward(); 
	virtual void init_pose_features(int N_feats = 300);
	virtual void gen_observations(); 
};

