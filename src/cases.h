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

	vector<int> m_pose_id; // has been seen in different frames 
	vector<int> m_feat_idx; // feature index in the frame at different pose 
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
};

class Case{

public:
	Case(); 
	virtual ~Case(); 

	vector<vector<double>> mv_poses;  // each pose has seven parameters x,y,z,qx,qy,qz,qw 
	vector<vector<FeatureMeasurement>> mv_obs; // observations at each pose 

	virtual void init_pose_features(); // init features and poses 
	virtual void gen_observations();   // generate observations 
	vector<Feature> mv_feats; 			// ground truth features 
	vector<vector<double>> mv_gt_poses; // ground truth poses
};

class Case_forward : public Case
{
public:
	Case_forward(); 
	virtual ~Case_forward(); 
	virtual void init_pose_features();
	virtual void gen_observations(); 
};

