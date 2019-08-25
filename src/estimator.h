/*

	Aug. 20, 2019, He Zhang, hzhang8@vcu.edu 

	estimator, given a case build optimizaton structure 

*/

#pragma once 

class Case;

const int WINDOW_SIZE = 10;
const int NUM_OF_F = 1000;

enum SIZE_PARAMETERIZATION
{
    SIZE_POSE = 7,
    SIZE_SPEEDBIAS = 9,
    SIZE_FEATURE = 1
};


enum FACTOR_TYPE
{
    TRANSFER_E = 1, 
    SAMPSON, 
    SAMPSON_D 
};

class Estimator
{
public:
	Estimator(); 
	~Estimator(); 

    void before_optimize(Case& ca);

    void optimize(Case* , double* perr = 0); 

	int m_cnt_pose; 
	int m_cnt_feat;

	double para_Pose[WINDOW_SIZE + 1][SIZE_POSE];
    double para_SpeedBias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS];
    double para_Feature[NUM_OF_F][SIZE_FEATURE];
    double para_Ex_Pose[2][SIZE_POSE];
    double para_Retrive_Pose[SIZE_POSE];
    double para_Td[1][1];
    double para_Tr[1][1];

    FACTOR_TYPE m_factor_type; 

    // Matrix3d ric[2];
    // Vector3d tic[2];

    // Vector3d        Ps[(WINDOW_SIZE + 1)];
    // Vector3d        Vs[(WINDOW_SIZE + 1)];
    // Matrix3d        Rs[(WINDOW_SIZE + 1)];
    // Vector3d        Bas[(WINDOW_SIZE + 1)];
    // Vector3d        Bgs[(WINDOW_SIZE + 1)];

};