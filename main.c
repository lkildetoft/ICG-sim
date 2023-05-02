/*
                ██╗ ██████╗ ██████╗       ███████╗██╗███╗   ███╗                   
                ██║██╔════╝██╔════╝       ██╔════╝██║████╗ ████║
                ██║██║     ██║  ███╗█████╗███████╗██║██╔████╔██║
                ██║██║     ██║   ██║╚════╝╚════██║██║██║╚██╔╝██║
                ██║╚██████╗╚██████╔╝      ███████║██║██║ ╚═╝ ██║
                ╚═╝ ╚═════╝ ╚═════╝       ╚══════╝╚═╝╚═╝     ╚═╝

ICG-SIM: A Monte-Carlo simulation of light propagation and fluorescence from 
indocyanine-green (ICG) in tissues. Can easily be modified for other fluorophores.
Specifically developed to simulate mal-perfusion in the anastomosis after esophagectomy,
as such allows for three distinct areas with ICG-concentrations growing at different speeds.

Author: Love Kildetoft

Based on MCML by Steven Jaques et. al. and its extension to simulate 
fluorescence, but written and adapted from scratch. 
Saves the current "face" and depth image of the absorbed and fluorescence 
photons at each iteration. These files may become very large in the end 
depending on the amount of photons/iterations/binning, so tread 
carefully! 

To run a simulation: 
Right now, the parameters have to modified directly in the main function. If
you stricly follow my guide here, nothing should break :) All are set to 
appropriate defaults corresponding to the human esophagus. 
After the parameters have been set, compile and run using your favourite 
compiler! 

Parameters:
    int niters: total number of iterations
    *******************************************************************************
    int nvals: the number of values (granularity/resolution) for the sim,
    eg. the number of paths/angles/total grid on which the photons may 
    propagate.
    *******************************************************************************
    int nphotons: total number of "photons" used in the simulation, usually
    a very large number. 
    *******************************************************************************
    float ls_x: x-position of the light source
    float ls_y: y-position of the light source
    float ls_z: z-position of the light source
    float ls_sp_r: radial spread, a point source can be simulated by making 
    this variable very small.
    *******************************************************************************
    float e_xmax: x-length of the tissue.
    float e_ymax: y-length of the tissue.
    float e_zmax: z-length of the tissue.
    float e_mus_ex: scattering coefficient of the tissue, for excitation photons.
    float e_mua_ex: absorption coefficient of the tissue, for excitation photons.
    float e_mus_fl: scattering coefficient of the tissue, for fluorescence photons.
    float e_mua_fl: absorption coefficient of the tissue, for fluorescence photons.
    *******************************************************************************
    float i_mus: scattering coefficient of indocyanine green.
    float i_mua: absorption coefficient of indocyanine green.
    *******************************************************************************
    float g_ex: overall anisotropy parameter for excitation photons.
    float g_fl: overall anisotropy parameter for fluorescence photons.
    *******************************************************************************
    float injected_conc: total "injected"/maximum ICG-concentration
    float mal_conc: initial concentration for mal-perfusive region.
    float icg_conc_upper: initial concentration above mal-perfusive region.
    float icg_conc_lower: initial concentration below mal-perfusive region.
    *******************************************************************************
    float mal_factor: ICG-concentration increase at each iteration, mal-perfusive 
    region.
    float icg_upper_factor: ICG-concentration increase at each iteration, upper
    region.
    float icg_lower_factor: ICG-concentration increase at each iteration, lower
    region.
    *******************************************************************************
    float l_mal_area: lower z-bound for mal-perfusive region. 
    float h_mal_area: upper z-bound for mal-perfusive region. 
    *******************************************************************************
    float fl_ratio: 1 - probability for photon to undergo fluorescence by ICG.
    *******************************************************************************
    int nrows: number of rows in the output histograms/images.
    int ncols: number of columns in the output histograms/images.
    *******************************************************************************                                                   
    Good luck!
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "Distrs.h"
#include "Simulation.h"
#include "RandChoice.h"
#include "ddHist.h"

int main(int argc, const char * argv[]) {
    //Sim parameters, see above.
    //clock_t begin = clock();
    int niters = 2000; 
    int nvals = 300;
    int nphotons = (int)1e8;
    
    float ls_x = 1;
    float ls_y = 0;
    float ls_z = 12.5;
    float ls_sp_r = 12.5;

    float e_xmax = 2;
    float e_ymax = 2;
    float e_zmax = 5;
    float e_mus_ex = 29.5;
    float e_mua_ex = 3.3;
    float e_mus_fl = 33;
    float e_mua_fl = 3.3;

    float i_mus = 2.1;
    float i_mua = 209;

    float g_ex = 0.85;
    float g_fl = 0.8;

    float inj_conc = 0.25;
    float icg_upper_factor = 0.1;
    float icg_lower_factor = 0.001;

    int ana_idx = 24;

    int fl_ratio = 20;

    int nrows = 32;
    int ncols = 24;

    float *conc_arr = calloc(nrows, sizeof(float));

    /*
    Here, the actual simulation starts and no more parameters should
    require modification. 
    */
    lightsource *source = &(lightsource){ls_x, ls_y, ls_z, ls_sp_r};
    esophagus *eso = &(esophagus){e_xmax, e_ymax, e_zmax, e_mus_ex, e_mua_ex, e_mus_fl, e_mua_fl};
    indocyn *icg = &(indocyn){i_mus, i_mua};
    int npackets = nphotons/niters;
    float fl_choice[2] = {1, 0};
    float fl_prob[2] = {(float)(1 - 1/fl_ratio), (float)(1/fl_ratio)};
    //Seed random generator
    srand((unsigned)time(NULL));
    //Create output files and their corresponding file-pointers.
    FILE *frontFluor = fopen("front_fl_frames.txt", "w");
    FILE *frontAbs = fopen("front_abs_frames.txt", "w");
    FILE *depthFluor = fopen("depth_fl_frames.txt", "w");
    FILE *depthAbs = fopen("depth_abs_frames.txt", "w");
    //Allocate space in memory for everything needed for the sim
    photon *photons = calloc(nphotons, sizeof(photon));

    float *theta = calloc(nvals, sizeof(float));

    float *t_plist_ex = calloc(nvals, sizeof(float));
    float *t_plist_fl = calloc(nvals, sizeof(float));

    float *rl_edges = calloc(nrows, sizeof(float));
    float *rh_edges = calloc(nrows, sizeof(float));
    genBinEdges(rl_edges, rh_edges, 0, eso->zmax, nrows);
    printf("Anastomosis located at %f cm", (rl_edges[ana_idx] + rh_edges[ana_idx])/2);
    float *cl_edges = calloc(nrows, sizeof(float));
    float *ch_edges = calloc(nrows, sizeof(float));
    genBinEdges(cl_edges, ch_edges, 0, eso->xmax, ncols);

    int **front_hist_fl = calloc(nrows, sizeof(int*));
    int **front_hist_abs = calloc(nrows, sizeof(int*));
    int **depth_hist_fl = calloc(nrows, sizeof(int*));
    int **depth_hist_abs = calloc(nrows, sizeof(int*));

    for (int i = 0; i < nrows; i++) {
        front_hist_fl[i] = calloc(ncols, sizeof(int));
        front_hist_abs[i] = calloc(ncols, sizeof(int));
        depth_hist_fl[i] = calloc(ncols, sizeof(int));
        depth_hist_abs[i] = calloc(ncols, sizeof(int));
    };
    //Fill scattering phase function arrays
    genPhaseVals(t_plist_ex, nvals, pi, g_ex);   
    genPhaseVals(t_plist_fl, nvals, pi, g_fl);  
    //Fill polar angle array
    int k = 0;
    for (float i = -pi; i < pi; i+=(pi)/(float)nvals) {
        theta[k++] = i;
    };
    
    printf("Initializing %d photons, please wait \n", nphotons);
    //Fill photon array with appropriate homogenous distribution over the tissue surface.
    #pragma omp parallel for
    for (int i = 0; i < nphotons; i++) {
        float ang = 2*pi*((float)rand()/(float)RAND_MAX);
        float sp = (source->spread_r)*((float)rand()/(float)RAND_MAX);
        photons[i] = (photon){((source->x)+(sp)*cos(ang)), source->y, ((source->z) + (sp)*sin(ang)), 1.0, 0, 0, 0};
    };    

    printf("Simulation started with %d photons, please wait \n", nphotons);
    //Main simulation loop
    int curr_iter_max = npackets;
    for (int iter = 0; iter < niters; iter++) {
        //Keep track of the number of absorbed/fluorescense photons each iteration.
        printf("Running simulation, %f %% done \n", ((float)iter/niters)*100);
        printf("Working on photons 0 to %d \n", curr_iter_max);
        //Update the ICG-concentration in each area
        for (int i = 0; i < nrows; i++) {
            if (conc_arr[i] < inj_conc) {
                conc_arr[i] += flowGrad(i, ana_idx, inj_conc, icg_lower_factor, icg_upper_factor);
            } else {
            }
        }
        //Loop over all photons
        #pragma omp parallel for
        for (int q = 0; q < curr_iter_max; q++) {
            photon *p;
            p = &photons[q];
            //Check fluorescence
            if ((p->w <= 0.1) && (p->outside == 0) && (p->curr_abs == 1)) {
                float fl_event = randChoice(fl_choice, fl_prob, 2);
                if (fl_event == 1) {
                    p->w = 1.0;
                    p->fluor = 1;
                }
            }
            //Let the current photon interact/fluoress appropriately depending on its current parameters
            if ((p->w > 0) && (p->outside == 0) && p->fluor == 0) {
                if (((p->x)*(p->x - eso->xmax) > 0) || ((p->y)*(p->y - eso->ymax) > 0) || ((p->z)*(p->z - eso->zmax) > 0)) {
                    p->outside = 1;
                } else {
                    for (int i = 0; i < nrows; i++) {
                        float curr_bin_l = rl_edges[i];
                        float curr_bin_h = rh_edges[i];
                        if ((p->z - curr_bin_l)*(p->z - curr_bin_h) <= 0) {
                            interact(p, eso, icg, conc_arr[i], theta, t_plist_ex, nvals);
                        }
                    }
                }
            } else if ((p->w > 0) && (p->outside == 0) && p->fluor == 1) {
                if (((p->x)*(p->x - eso->xmax) > 0) || ((p->y)*(p->y - eso->ymax) > 0) || ((p->z)*(p->z - eso->zmax) > 0)) {
                    p->outside = 1;
                } else {
                    for (int i = 0; i < nrows; i++) {
                        float curr_bin_l = rl_edges[i];
                        float curr_bin_h = rh_edges[i];
                        if ((p->z - curr_bin_l)*(p->z - curr_bin_h) <= 0) {
                            fluoress(p, eso, icg, conc_arr[i], theta, t_plist_fl, nvals);
                        }
                    }
                }
            }
            
            //Generate histograms
            if ((p->w <= 0.1) && (p->outside == 0) && (p->curr_abs == 1) && (p->fluor == 0)) {
                genHist(front_hist_abs, rl_edges, rh_edges, cl_edges, ch_edges, nrows, ncols, p->x, p->z);
                genHist(depth_hist_abs, rl_edges, rh_edges, cl_edges, ch_edges, nrows, ncols, p->y, p->z);
            } 
    
            if ((p->w > 0) && (p->fluor == 1) && (p->outside == 0) && (p->y <= 0.1)) {
                genHist(front_hist_fl, rl_edges, rh_edges, cl_edges, ch_edges, nrows, ncols, p->x, p->z);
                genHist(depth_hist_fl, rl_edges, rh_edges, cl_edges, ch_edges, nrows, ncols, p->y, p->z);
            }
        }
        //Write each histogram at each iteration to the appropriate file
        for (int row = 0; row < nrows; row++) {
            for (int col = 0; col < ncols; col++) {
                fprintf(frontFluor, "%d", front_hist_fl[row][col]);
                fprintf(depthFluor, "%d", depth_hist_fl[row][col]);
                fprintf(frontAbs, "%d", front_hist_abs[row][col]);
                fprintf(depthAbs, "%d", depth_hist_abs[row][col]);
                if (col < (ncols - 1)) {
                    fprintf(frontFluor, ",");
                    fprintf(depthFluor, ",");
                    fprintf(frontAbs, ",");
                    fprintf(depthAbs, ",");
                } 
            }
            fprintf(frontFluor, "\n");
            fprintf(depthFluor, "\n");
            fprintf(frontAbs, "\n");
            fprintf(depthAbs, "\n");
        }
        curr_iter_max += npackets;
    }
    //Finished! Nice
    printf("Done \n");
    //clock_t end = clock();
    //double benchmark = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("%lf s", benchmark);
    //Free up memory
    fclose(frontFluor);
    fclose(frontAbs);
    fclose(depthFluor);
    fclose(depthAbs);

    free(photons);
    free(theta);
    free(t_plist_ex);
    free(t_plist_fl);
    for (int i = 0; i < nrows; i++) {
        free(front_hist_fl[i]);
        free(front_hist_abs[i]);
        free(depth_hist_fl[i]);
        free(depth_hist_abs[i]);
    };
    free(front_hist_fl);
    free(front_hist_abs);
    free(depth_hist_fl);
    free(depth_hist_abs);
    //Close all file-pointers


    return 0;
}
