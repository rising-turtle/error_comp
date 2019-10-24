/*
	
	Oct. 24 2019, He Zhang, hzhang8@vcu.edu 

	depth factor based estimator

*/

#pragma once

#include "estimator.h"

enum DEPTH_FACTOR_TYPE
{
    STEREO = 1, 
    RGBD, 
};


class EstimatorDepth : public Estimator 
{
public:
	EstimatorDepth(double sigma_dis = 1.); 
	virtual ~EstimatorDepth(); 
	
    virtual void optimize(Case* , double* perr = 0); 

    DEPTH_FACTOR_TYPE m_depth_factor_type; 
};