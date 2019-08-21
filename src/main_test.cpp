/*
	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 

	main test file 

*/

#include "projection_factor.h"
#include "estimator.h"
#include "cases.h"

int main(int argc, char* argv[])
{

	Case* pc = new Case_forward(); 
	pc->init_pose_features();
	pc->gen_observations();
	pc->add_noise(0.03);
	Estimator est; 
	est.optimize(pc); 

	return 1; 

}