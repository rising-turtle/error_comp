/*
	
	Oct. 24 2019, He Zhang, hzhang8@vcu.edu 

	depth factor based estimator

*/

#include "cases.h"
#include "estimator_depth_factor.h"
#include "projection_factor.h"
#include "projectionOneFrameTwoCamFactor.h"
#include "projectionTwoFrameOneCamFactor.h"
#include "projectionTwoFrameTwoCamFactor.h"
#include "pose_local_parameterization.h"
#include "depth_factor.h"

using namespace Eigen;

EstimatorDepth::EstimatorDepth(double sigma_dis):
Estimator(sigma_dis),
m_depth_factor_type(STEREO)
{

	ProjectionTwoFrameOneCamFactor::sqrt_info = sigma_dis * Matrix2d::Identity(); 
	ProjectionOneFrameTwoCamFactor::sqrt_info = sigma_dis * Matrix2d::Identity(); 
	ProjectionTwoFrameTwoCamFactor::sqrt_info = sigma_dis * Matrix2d::Identity();

	// ur = ul - rig_len * lambda; lambda = 1./depth
	SingleInvDepthFactor::sqrt_info = sigma_dis/FeatureMeasurement::rig_len; 
	ProjectionDepthFactor::sqrt_info = (sigma_dis/FeatureMeasurement::rig_len)*Eigen::Matrix3d::Identity();

}

EstimatorDepth::~EstimatorDepth(){}


void EstimatorDepth::optimize(Case* ca, double* perr)
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

	{
    	// ex pose
    	ceres::LocalParameterization *local_param = new PoseLocalParameterization; 
    	problem.AddParameterBlock(para_Ex_Pose[0], 7, local_param); 
    	problem.SetParameterBlockConstant(para_Ex_Pose[0]);
	}
	{
		ceres::LocalParameterization *local_param = new PoseLocalParameterization; 
    	problem.AddParameterBlock(para_Ex_Pose[1], 7, local_param); 
    	problem.SetParameterBlockConstant(para_Ex_Pose[1]);
	}

    // int feat_fac_cnt = 0; 

    double sum_proj_dis = 0; 
    int num_proj = 0; 

    int used_features = 0; 

    // add features 
    for (int i=0; i<m_cnt_feat; i++){
    	if(para_Feature[i][0] < 0) continue; 

        if(ca->mv_feats[i].mv_pose_id.size() >= 2){
            if(used_features++ >= ca->m_thre_cnt_feats){
                break; 
            }
        }

    	// 
    	int ref_id = ca->mv_feats[i].m_ref_id;
    	int ref_node_id = ca->mv_feats[i].mv_pose_id[ref_id]; 
    	int ref_feat_idx = ca->mv_feats[i].mv_feat_idx[ref_id]; 
    	Eigen::Vector3d pi(ca->mv_obs[ref_node_id][ref_feat_idx].xy[0], ca->mv_obs[ref_node_id][ref_feat_idx].xy[1], 1); 

    	Vector2d velocity_i(0, 0); 
    	Vector2d velocity_j(0, 0); 
    	double cur_td = 0;

    	for(int j=0; j<ca->mv_feats[i].mv_pose_id.size(); j++){
    		int node_id = ca->mv_feats[i].mv_pose_id[j];
    		int feat_id = ca->mv_feats[i].mv_feat_idx[j]; 
    		Eigen::Vector3d pj(ca->mv_obs[node_id][feat_id].xy[0], ca->mv_obs[node_id][feat_id].xy[1], 1.); 
    		Eigen::Vector3d right_pj(ca->mv_obs[node_id][feat_id].r_xy[0], ca->mv_obs[node_id][feat_id].r_xy[1], 1.); 
            double inv_depth = 1./ca->mv_obs[node_id][feat_id].depth; 

            double proj_dis = projDis(pi, pj, para_Feature[i][0], para_Pose[ref_node_id], para_Pose[node_id]);

            sum_proj_dis += proj_dis; 
            num_proj ++;

    		if(node_id != ref_node_id){


                if(m_depth_factor_type == STEREO){
                    // 1. add projection factor 
                    ProjectionTwoFrameOneCamFactor* f = new ProjectionTwoFrameOneCamFactor(pi, pj, velocity_i, velocity_j, cur_td, cur_td); 
                    problem.AddResidualBlock(f, loss_function, para_Pose[ref_node_id], para_Pose[node_id], para_Ex_Pose[0], para_Feature[i], para_Td[0]); 

    			    // 2. add stereo factor 
                    ProjectionTwoFrameTwoCamFactor* fs = new ProjectionTwoFrameTwoCamFactor(pi, right_pj, velocity_i, velocity_j, cur_td, cur_td); 
                    problem.AddResidualBlock(fs, loss_function, para_Pose[ref_node_id], para_Pose[node_id], para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[i], para_Td[0]);
                }else if(m_depth_factor_type == RGBD){
                    // 1. add projection depth factor 
                    ProjectionDepthFactor *fs = new ProjectionDepthFactor(pi, pj, inv_depth); 
                    problem.AddResidualBlock(fs, loss_function, para_Pose[ref_node_id], para_Pose[node_id], para_Ex_Pose[0], para_Feature[i]);
                }
    		}else{

                if(m_depth_factor_type == STEREO){
    			// 3. add stereo factor 
                    ProjectionOneFrameTwoCamFactor* fs = new ProjectionOneFrameTwoCamFactor(pi, right_pj, velocity_i, velocity_j, cur_td, cur_td); 
                    problem.AddResidualBlock(fs, loss_function, para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[i], para_Td[0]); 
                }else if(m_depth_factor_type == RGBD){
                    SingleInvDepthFactor* fs = new SingleInvDepthFactor(inv_depth); 
                    problem.AddResidualBlock(fs, loss_function, para_Feature[i]);
                }

    		}

            // ++feat_fac_cnt; 
    	}
    }
    if(m_depth_factor_type == STEREO)
        problem.SetParameterBlockConstant(para_Td[0]);

    cout <<"estimator.cpp: used_features: "<<used_features-1<<endl;
    cout <<"estimator.cpp: num_proj = "<<num_proj<<" average proj_dis: "<<(sum_proj_dis/(num_proj+1e-6))*240<<endl; 
    // cout <<"estimator.cpp: valid feature factor number = "<<feat_fac_cnt<<endl;

    // solve the problem 
    ceres::Solver::Options options; 
    options.linear_solver_type = ceres::DENSE_SCHUR; 
    options.trust_region_strategy_type = ceres::DOGLEG; 
    options.max_num_iterations = 100; 
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary; 
    ceres::Solve(options, &problem, &summary); 
	cout << summary.BriefReport() << endl;

    /*if(summary.termination_type == ceres::NUMERICAL_FAILURE){

        cout <<"this case is failed !"<<endl;
    }*/

	for(int i=1; i<m_cnt_pose; i++){
		for(int j=0; j<7; j++)
			cout<<para_Pose[i][j]<<" ";
		cout<<endl;
	}

    if(perr != NULL){
        ca->rmse(para_Pose, perr); 
    }

    return ;

}