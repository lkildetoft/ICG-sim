//
//  Simulation.h
//  PDT-SIM
//
//  Created by Love Kildetoft on 2022-11-23.
//

#ifndef Simulation_h
#define Simulation_h

typedef struct photon {
    float x;
    float y;
    float z;
    float w;
    int outside;
    int curr_abs;
    int fluor;
} photon;

typedef struct esophagus {
    float xmax;
    float ymax;
    float zmax;
    float mus_ex;
    float mua_ex;
    float mus_fl;
    float mua_fl;
} esophagus;

typedef struct lightsource {
    float x;
    float y;
    float z;
    float spread_r;
} lightsource;

typedef struct indocyn {
    float mus;
    float mua;
} indocyn;

void genPhaseVals(float *plist, int size, float lim, float g);
void genExpVals(float *glist, int size, float lim, float mu);
void genUniform(float *ulist, int nvals);
void interact(photon *p, esophagus *eso, indocyn *icg, float conc, float *t_angles, float *t_phasevals, int len);
void fluoress(photon *p, esophagus *eso, float *t_angles, float *t_phasevals, int len);
#endif /* Simulation_h */
