/*
Author: Love Kildetoft
*/
#include <math.h>
#include <stdlib.h>

#include "Distrs.h"
#include "RandChoice.h"
#include "Simulation.h"

void genPhaseVals(float *plist, int size, float lim, float g) {
    /*
    Modify an empty float array to contain scattering phase function
    values. 
    */
    int n = 0;
    for (float i = -1 * lim; i < lim; i += (2 * lim / size)) {
        plist[n++] = phaseFunc(i, g);
    };
};

void genUniform(float *ulist, int nvals) {
    /* 
    Just the same value over and over. In statistics we call it a 
    uniform distribution :D
    */
    for (int i = 0; i < nvals; i++) {
        ulist[i] = (float)1/nvals;
    };
};

void interact(photon *p, esophagus *eso, indocyn *icg, float conc, float *t_angles, float *t_phasevals, int len) {
    /*
    Simulates the interaction of photons in a medium with appropriate parameters.
    A "photon" is absorbed and scattered at each call.
    */
    float w = p->w;
    //Account for a certain concentration of ICG in the tissue
    float choices[2] = {0, 1};
    float choice_prob[2] = {conc, 1.0-conc};

    float mua;
    float mus;
    //Find the appropriate coefficients
    int choice = randChoice(choices, choice_prob, 2);
    if (choice == 0) {
        mua = icg->mua;
        mus = icg->mus;
        p->curr_abs = 1;
    } else {
        mua = eso->mua_ex;
        mus = eso->mus_ex;
        p->curr_abs = 0;
    };

    float mue = mua + mus;
    float l;
    //Calculate new path from angles and randomly generated path length
    float theta = randChoice(t_angles, t_phasevals, len);
    float phi = 2*pi*((float)rand()/(float)RAND_MAX);
    float path = (float)rand()/(float)RAND_MAX;
    if (path == 0) {
        l = 0;
    } else {
        l = -1*log(path)/(mue);
    }
    float xs = l*sin(theta)*cos(phi);
    float ys = l*sin(theta)*sin(phi);
    float zs = l*cos(theta);
    //Update the weight and coordinates of the photon
    p->w -= (mua/mue)*w;
    p->x += xs;
    p->y += ys;
    p->z += zs;
};

void fluoress(photon *p, esophagus *eso, float *t_angles, float *t_phasevals, int len) {
    /*
    Simulates the interaction of photons in a medium with appropriate parameters.
    A "photon" is absorbed and scattered at each call.
    */
    float w = p->w;
    float mua = eso->mua_fl;
    float mus = eso->mus_fl;

    float mue = mua + mus;
    
    float theta = randChoice(t_angles, t_phasevals, len);
    float phi = 2*pi*((float)rand()/(float)RAND_MAX);
    float path = (float)rand()/(float)RAND_MAX;
    float l;
    if (path == 0) {
        l = 0;
    } else {
        l = -1*log(path)/(mue);
    }    
    float xs = l*sin(theta)*cos(phi);
    float ys = l*sin(theta)*sin(phi);
    float zs = l*cos(theta);

    p->w -= (mua/mue)*w;
    p->x += xs;
    p->y += ys;
    p->z += zs;
};
