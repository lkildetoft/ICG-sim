#include "ddHist.h"
/*
Author: Love Kildetoft
*/
void genBinEdges(float *lbin_edge, float *hbin_edge, float low, float high, int nbins) {
    float width = (high - low)/nbins;
    float curr_bin = low;
    for (int i = 0; i < nbins; i++) {
        lbin_edge[i] = curr_bin;
        hbin_edge[i] = curr_bin + width;
        curr_bin += width;
    }
};

void genHist(int **hist, float *rl_edges, float *rh_edges, float *cl_edges, float *ch_edges, int nrows, int ncols, float x, float y) {
    int row_idx, col_idx;

    for (int i = 0; i < nrows; i++) {
        float curr_row_l = rl_edges[i];
        float curr_row_h = rh_edges[i];

        for (int j = 0; j < ncols; j++) {
            float curr_col_l = cl_edges[j];
            float curr_col_h = ch_edges[j];
            int row_cond = 0;
            int col_cond = 0;
            
            if ((y > curr_row_l) && (y < curr_row_h)) {
                row_idx = i;
                row_cond = 1;
            }

            if ((x > curr_col_l) && (x < curr_col_h)) {
                col_idx = j;
                col_cond = 1;
            } 
            
            if ((y == curr_row_l) && (i != 0)) {
                row_idx = i - 1;
                row_cond = 1;
            } 
            
            if ((y == curr_row_h) && (j != (nrows - 1))) {
                row_idx = i + 1;
                row_cond = 1;
            } 
            
            if ((x == curr_col_l) && (j != 0)) {
                col_idx = j - 1;
                col_cond = 1;
            } 
            
            if ((x == curr_col_h) && (j != (ncols - 1))) {
                col_idx = j + 1;
                col_cond = 1;
            }

            if ((row_cond == 1) && (col_cond == 1)) {
                hist[row_idx][col_idx] += 1;
            } else {
                continue;
            }
        }
    }
};
