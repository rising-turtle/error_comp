/*

	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 

	estimator, given a ca build optimizaton structure 

*/


#include "estimator.h"
#include "cases.h"
#include "projection_factor.h"
#include "pose_local_parameterization.h"
#include <iostream>

using namespace std; 

Estimator::Estimator():m_cnt_pose(0),
m_cnt_feat(0){}
Estimator::~Estimator(){}

void Estimator::before_optimize(Case& ca)
{
	m_cnt_pose = ca.mv_poses.size(); 

	if(m_cnt_pose > WINDOW_SIZE) {
		cerr << "estimator.cpp: m_cnt_pose = "<<m_cnt_pose<<" > WINDOW_SIZE: " <<WINDOW_SIZE<<endl; 
		return ; 
	}

	// assign poses 
	for(int i=0; i<m_cnt_pose; i++){
		for(int j=0; j<7; j++){
			para_Pose[i][j] = ca.mv_poses[i][j];
		}
	}

	// assign feature depth 
	m_cnt_feat = ca.mv_feats.size(); 
	if(m_cnt_feat >= NUM_OF_F){
		cerr << "estimator.cpp: m_cnt_feat = "<<m_cnt_feat<<" > NUM_OF_F: "<<NUM_OF_F<<endl; 
		return ; 
	}

	int cnt_valid = 0; 
	for(int i=0; i<m_cnt_feat; i++){
		if(ca.mv_feats[i].mv_pose_id.size() > 0){
			int node_id = ca.mv_feats[i].mv_pose_id[ca.mv_feats[i].m_ref_id]; // first observation 
			int obs_idx = ca.mv_feats[i].mv_feat_idx[ca.mv_feats[i].m_ref_id];
			para_Feature[i][0] = 1.0/(ca.mv_obs[node_id][obs_idx].depth);
			++cnt_valid; 
		}
		else{
			para_Feature[i][0] = -1; 
		}
	}

	// assign ric 
	for(int i=0; i<6; i++)
		para_Ex_Pose[0][i] = 0; 
	para_Ex_Pose[0][6] = 1;
	para_Tr[0][0] = 0; 
	para_Td[0][0] = 0; 
	cout <<"estimator.cpp: m_cnt_pose = "<<m_cnt_pose<<" m_cnt_feat = "<<m_cnt_feat<<endl;
	cout <<"estimator.cpp: valid feat: "<<cnt_valid<<endl;
	return ; 
}


void Estimator::optimize(Case* ca)
{
	// build structure 
	before_optimize(*ca); 

	// build optimization factors
	ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = NULL;
    loss_function = new ceres::HuberLoss(1.0);

    // add pose 
    for (int i = 0; i < m_cnt_pose; i++){
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
    }
    problem.SetParameterBlockConstant(para_Pose[0]);

    // time offset 
    // problem.AddParameterBlock(para_Td[0], 1);
    // problem.SetParameterBlockConstant(para_Td[0]); 

    // ex pose
    ceres::LocalParameterization *local_param = new PoseLocalParameterization; 
    problem.AddParameterBlock(para_Ex_Pose[0], 7, local_param); 
    problem.SetParameterBlockConstant(para_Ex_Pose[0]);


    int feat_fac_cnt = 0; 
    // add features 
    for (int i=0; i<m_cnt_feat; i++){
    	if(para_Feature[i][0] < 0) continue; 

    	// 
    	int ref_id = ca->mv_feats[i].m_ref_id;
    	int ref_node_id = ca->mv_feats[i].mv_pose_id[ref_id]; 
    	int ref_feat_idx = ca->mv_feats[i].mv_feat_idx[ref_id]; 
    	Eigen::Vector3d pi(ca->mv_obs[ref_node_id][ref_feat_idx].xy[0], ca->mv_obs[ref_node_id][ref_feat_idx].xy[1], 1); 
    	for(int j=0; j<ca->mv_feats[i].mv_pose_id.size(); j++){
    		int node_id = ca->mv_feats[i].mv_pose_id[j];
    		if(node_id == ref_node_id) continue; 
    		int feat_id = ca->mv_feats[i].mv_feat_idx[j]; 
    		Eigen::Vector3d pj(ca->mv_obs[node_id][feat_id].xy[0], ca->mv_obs[node_id][feat_id].xy[1], 1.); 
    		
    		ProjectionFactor *f = new ProjectionFactor(pi, pj); 
    		f->sqrt_info = 24 * Eigen::Matrix2d::Identity(); 

    		problem.AddResidualBlock(f, loss_function, para_Pose[ref_node_id], para_Pose[node_id], 
    			para_Ex_Pose[0], para_Feature[i]);
    		problem.SetParameterBlockConstant(para_Feature[i]); 
    		++feat_fac_cnt; 

    	}
    }

    cout <<"estimator.cpp: feature factor number = "<<feat_fac_cnt<<endl;

    // solve the problem 
    ceres::Solver::Options options; 
    options.linear_solver_type = ceres::DENSE_SCHUR; 
    options.trust_region_strategy_type = ceres::DOGLEG; 
    options.max_num_iterations = 100; 
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary; 
    ceres::Solve(options, &problem, &summary); 
	cout << summary.BriefReport() << endl;

	for(int i=1; i<m_cnt_pose; i++){
		for(int j=0; j<7; j++)
			cout<<para_Pose[i][j]<<" ";
		cout<<endl;
	}

}