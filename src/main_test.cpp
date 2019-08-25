/*
	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 

	main test file 

*/

#include "projection_factor.h"
#include "estimator.h"
#include "cases.h"
#include "mean_std.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string> 

using namespace std ; 

struct Result
{
public:
	Result(int n, double noise_, double m1, double std1, double m2, double std2):
	nfeats(n), noise(noise_), mean_rmse_dis(m1), std_rmse_dis(std1), mean_rmse_ori(m2), 
	std_rmse_ori(std2)
	{}
	Result(){}
	int nfeats; 
	double noise; 
	double mean_rmse_dis; 
	double std_rmse_dis;
	double mean_rmse_ori;
	double std_rmse_ori; 
};

void save_in_file(string f_name, std::vector<Result>& v);

double run_once(int n_feature, double noise,  double& err_dis, double& err_ori ); 

void test_noise(int N_feats = 300); 

void test_feat_num(double noise = 0.03); 

FACTOR_TYPE g_fac_type = FACTOR_TYPE::TRANSFER_E; 

int main(int argc, char* argv[])
{

	test_noise();

	return 1; 

}


void test_noise(int N_feats)
{
	int N = 20; 
	vector<double> noise{0.01, 0.04, 0.07}; // , 0.1, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0}; 
	vector<double> ve_dis(N); 
	vector<double> ve_ori(N); 
	vector<Result> v_rec(noise.size()); 
	double err_dis, err_ori;

	double mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori; 

	for(int k=0; k<noise.size(); k++)
	{
		for(int i=0; i<N; i++)
		{
			run_once(N_feats, noise[k], err_dis, err_ori); 
			ve_dis[i] = err_dis; 
			ve_ori[i] = err_ori; 
		}

		compute_mu_sigma(ve_dis.data(), ve_dis.size(), mean_rmse_dis, std_rmse_dis);
		compute_mu_sigma(ve_ori.data(), ve_ori.size(), mean_rmse_ori, std_rmse_ori); 

		Result r(N_feats, noise[k], mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori);
		v_rec[k] = r;
	}
	// cout << "err_dis = "<<err_dis<<" err_ori = "<<err_ori<<endl;

	if(g_fac_type == FACTOR_TYPE::TRANSFER_E){
		save_in_file("transfer_err.log", v_rec);
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON){
		save_in_file("sampson_err.log", v_rec); 
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_D){
		save_in_file("sampson_lambda_err.log", v_rec); 
	}
	return; 
}

void test_feat_num(double noise)
{
	int N = 20; 
	vector<double> v_n_feats{10, 20, 30}; // , 0.1, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0}; 
	vector<Result> v_rec(v_n_feats.size()); 
	vector<double> ve_dis(N); 
	vector<double> ve_ori(N); 
	double err_dis, err_ori;

	double mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori; 

	for(int k=0; k<v_n_feats.size(); k++)
	{
		int N_feats = v_n_feats[k];
		
		for(int i=0; i<N; i++)
		{
			run_once(N_feats, noise, err_dis, err_ori); 
			ve_dis[i] = err_dis; 
			ve_ori[i] = err_ori; 
		}

		compute_mu_sigma(ve_dis.data(), ve_dis.size(), mean_rmse_dis, std_rmse_dis);
		compute_mu_sigma(ve_ori.data(), ve_ori.size(), mean_rmse_ori, std_rmse_ori); 

		Result r(N_feats, noise, mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori);
		v_rec[k] = r;
	}
	// cout << "err_dis = "<<err_dis<<" err_ori = "<<err_ori<<endl;

	if(g_fac_type == FACTOR_TYPE::TRANSFER_E){
		save_in_file("transfer_feat_err.log", v_rec);
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON){
		save_in_file("sampson_feat_err.log", v_rec); 
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_D){
		save_in_file("sampson_lambda_feat_err.log", v_rec); 
	}
	return; 
}

double run_once(int n_feature, double noise, double& err_dis, double& err_ori )
{

	Case* pc = new Case_forward(); 
	pc->init_pose_features(n_feature);
	pc->gen_observations();
	pc->add_noise(noise);

	double err[2]; 
	Estimator est; 
	est.m_factor_type = g_fac_type; // FACTOR_TYPE::TRANSFER_E; 
	est.optimize(pc, err); 

	err_dis = err[0]; 
	err_ori = err[1]; 

	delete pc;

	return err_dis;
}

void save_in_file(string f_name, std::vector<Result>& v)
{
	ofstream ouf(f_name.c_str()); 

	for(int i=0; i<v.size(); i++){
		Result& r = v[i]; 
		ouf << r.nfeats<<"\t "<<r.noise<<"\t "<<r.mean_rmse_dis
		<<"\t "<<r.std_rmse_dis<<"\t "<<r.mean_rmse_ori<<"\t "<<r.std_rmse_ori<<endl; 
	}
	ouf.close();
}