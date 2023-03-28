#ifndef ddHist_h
#define ddHist_h

void genBinEdges(float *lbin_edge, float *hbin_edge, float low, float high, int nbins);
void genHist(int **hist, float *rl_edges, float *rh_edges, float *cl_edges, float *ch_edges, int nrows, int ncols, float x, float y);

#endif /* ddHist_h */