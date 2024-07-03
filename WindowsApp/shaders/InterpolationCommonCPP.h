#ifndef _INTERPOLATION_COMMON_CPP_H
#define _INTERPOLATION_COMMON_CPP_H

struct 
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
GConstants
{
	uint2 NumControlPoints;
	uint2 NumOuputPoints;
	float TCBTension;
	float TCBContinuity;
	float TCBBias;
};

#define INTERPOLATION_METHOD_CR 0
#define INTERPOLATION_METHOD_TCB 1
#define INTERPOLATION_METHOD_LAGRANGE 2
//...
#define INTERPOLATION_METHOD_MAX 3

#endif//INTERPOLATION_COMMON_CPP_H