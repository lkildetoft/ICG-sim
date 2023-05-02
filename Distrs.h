#ifndef Distrs_h
#define Distrs_h
/*
Author: Love Kildetoft
*/
extern const float pi;
float phaseFunc(float theta, float g);
float flowGrad(float idx, float max_idx, float max_conc, float max_low, float grad_low, float high_factor);
float expFunc(float x, float mue);

#endif /* Distrs_h */
