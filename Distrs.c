#include "Distrs.h"
#include <math.h>
#include <stdio.h>

/*
Author: Love Kildetoft
*/
//Everyones favourite numerical constant
const float pi = 3.141592653589793;

float phaseFunc(float theta, float g) {
    /*
    Henvey-Greenstein scattering phase function. 
    
    Arguments:
    theta: angle in radians. Range depends on coordinate 
    system. 
    g: anisotropy factor of medium, between -1 and 1.

    Returns:
    Value of the scattering phase function at theta, g. 
    */
    float p = (float)((1/(4*pi))*((1-pow(g, 2))/(pow(1+pow(g, 2)-2*g*cos(theta), 3/2))));
    return p;
};

float flowGrad(float idx, float max_idx, float max_conc, float max_low, float grad_low, float high_factor) {
    float grad;
    if (idx > max_idx) {
        grad = high_factor;
    } else if (idx < max_idx) {
        if ((max_low - grad_low*(idx/max_idx) >= 0)) {
            grad = max_low - grad_low*(idx/max_idx);
        } else {
            grad = 0;
        }
    };
    return grad;
}

float expFunc(float x, float mue) {
    /*
    Exponentially decaying distribution (see Beer-Lambert law).
    */
    float p = (float)(mue*exp(-mue*x));
    return p;
};
