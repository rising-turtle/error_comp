/*
	
	Oct. 24 2019, He Zhang, hzhang8@vcu.edu 

	depth factor based estimator

*/

#include "estimator_depth_factor.h"
#include "projection_factor.h"
#include "projectionOneFrameTwoCamFactor.h"
#include "projectionTwoFrameOneCamFactor.h"
#include "projectionTwoFrameTwoCamFactor.h"

EstiamatorDepth::EstiamatorDepth(double sigma_dis):
Estimator(sigma_dis)
{}

EstiamatorDepth::~EstiamatorDepth(){}



void EstiamatorDepth::optimize(Case* ca, double* perr)
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
    	for(int j=0; j<ca->mv_feats[i].mv_pose_id.size(); j++){
    		int node_id = ca->mv_feats[i].mv_pose_id[j];
    		if(node_id == ref_node_id) continue; 
    		int feat_id = ca->mv_feats[i].mv_feat_idx[j]; 
    		Eigen::Vector3d pj(ca->mv_obs[node_id][feat_id].xy[0], ca->mv_obs[node_id][feat_id].xy[1], 1.); 
    		
            // detect whether this pair of matched features is an inlier 

            double proj_dis = projDis(pi, pj, para_Feature[i][0], para_Pose[ref_node_id], para_Pose[node_id]);

            sum_proj_dis += proj_dis; 
            num_proj ++;

            // if(proj_dis > m_std_dis*3){
                // continue; 
            // }
       

            {
    		    ProjectionFactor *f = new ProjectionFactor(pi, pj); 
                problem.AddResidualBlock(f, loss_function, para_Pose[ref_node_id], para_Pose[node_id], 
                para_Ex_Pose[0], para_Feature[i]);
            }
            problem.SetParameterBlockConstant(para_Feature[i]); 
            ++feat_fac_cnt; 
            
            // 

    		// Eigen::Matrix<double, 4, 4> sqrt_info = 24 * Eigen::Matrix<double, 4, 4>::Identity();
    		// ceres::CostFunction* f =
      			// new ceres::AutoDiffCostFunction<SampsonCostFunctor, 4, 7, 7, 7, 1>(
        			// new SampsonCostFunctor(pi, pj, sqrt_info));
    		
    	}
    }

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