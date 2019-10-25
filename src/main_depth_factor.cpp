/*

	Oct. 24 2019, He Zhang, hzhang8@vcu.edu 

	test depth factor, depth are available in multiple frames 

*/

#include "estimator_depth_factor.h"
#include "cases.h"
#include "mean_std.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

double focal_length = 240.; 

double run_once(int n_feature, double noise,  double& err_dis, double& err_ori ); 

int main(int argc, char* argv[])
{

	double e1, e2;
    run_once(50, 2./focal_length,e1, e2 ); 

	return 1; 

}


double run_once(int n_feature, double noise, double& err_dis, double& err_ori )
{

	Case* pc = new Case_forward(); 
	pc->init_pose_features(n_feature);
	pc->gen_observations();
	pc->add_noise(noise);

	double err[2]; 
	EstimatorDepth est(noise); 
	est.optimize(pc, err); 

	err_dis = err[0]; 
	err_ori = err[1]; 

	delete pc;

	return err_dis;
}
