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

void print_result(std::vector<Result>& v);

void save_in_file(string f_name, std::vector<Result>& v);

double run_once(int n_feature, double noise,  double& err_dis, double& err_ori ); 

void test_noise(int N_feats = 300); 

void test_feat_num(double noise = 0.03); 

FACTOR_TYPE g_fac_type = FACTOR_TYPE::SAMPSON_CD; // SAMPSON_C SAMPSON_D SAMPSON; //TRANSFER_E; 

int main(int argc, char* argv[])
{

	vector<FACTOR_TYPE> st{FACTOR_TYPE::SAMPSON_C, FACTOR_TYPE::TRANSFER_E}; 

	for(auto& a : st){
		g_fac_type = a;
		test_feat_num(5.); 
		test_noise(30);
	}

	return 1; 

}

void test_noise(int N_feats)
{
	int N = 5000; // number of repeats
	vector<double> noise{1., 1.5, 2.0, 2.5, 3.0, 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10., 10.5, 11., 11.5, 12.}; // , 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1}; // , 0.1, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0}; 
	vector<double> ve_dis(N); 
	vector<double> ve_ori(N); 
	vector<Result> v_rec(noise.size()); 
	double focal_length = 240; // 300; 
	double err_dis, err_ori;

	double mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori; 

	for(int k=0; k<noise.size(); k++)
	{
		for(int i=0; i<N; i++)
		{
			run_once(N_feats, noise[k]/focal_length, err_dis, err_ori); 
			ve_dis[i] = err_dis; 
			ve_ori[i] = err_ori; 
		}

		compute_mu_sigma(ve_dis.data(), ve_dis.size(), mean_rmse_dis, std_rmse_dis);
		compute_mu_sigma(ve_ori.data(), ve_ori.size(), mean_rmse_ori, std_rmse_ori); 

		Result r(N_feats, noise[k], mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori);
		v_rec[k] = r;
	}
	// cout << "err_dis = "<<err_dis<<" err_ori = "<<err_ori<<endl;

	print_result(v_rec);

	if(g_fac_type == FACTOR_TYPE::TRANSFER_E){
		save_in_file("transfer_err.log", v_rec);
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON){
		save_in_file("sampson_err.log", v_rec); 
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_C){
		save_in_file("sampson_c_err.log", v_rec);
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_D){
		save_in_file("sampson_d_err.log", v_rec); 
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_CD){
		save_in_file("sampson_cd_err.log", v_rec); 
	}
	return; 
}

void test_feat_num(double noise)
{
	int N = 5000; 
	vector<double> v_n_feats{10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120}; // , 0.1, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0}; 
	vector<Result> v_rec(v_n_feats.size()); 
	vector<double> ve_dis(N); 
	vector<double> ve_ori(N); 

	double focal_length = 240; // 300; 

	double err_dis, err_ori;

	double mean_rmse_dis, std_rmse_dis, mean_rmse_ori, std_rmse_ori; 

	for(int k=0; k<v_n_feats.size(); k++)
	{
		int N_feats = v_n_feats[k];
		
		for(int i=0; i<N; i++)
		{
			run_once(N_feats, noise/focal_length, err_dis, err_ori); 
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
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_C){
		save_in_file("sampson_c_feat_err.log", v_rec);
	}else if(g_fac_type == FACTOR_TYPE::SAMPSON_D){
		save_in_file("sampson_d_feat_err.log", v_rec); 
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
	Estimator est(noise); 
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

void print_result(std::vector<Result>& v)
{	
	for(int i=0; i<v.size(); i++){
		Result& r = v[i]; 
		cout << r.nfeats<<"\t "<<r.noise<<"\t "<<r.mean_rmse_dis
		<<"\t "<<r.std_rmse_dis<<"\t "<<r.mean_rmse_ori<<"\t "<<r.std_rmse_ori<<endl; 
	}
}