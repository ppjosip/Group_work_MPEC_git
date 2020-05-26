#include "mex.h"
#include "qp_W.h"
#include "qp_gradJ.h"
#include "qp_gradhT.h"
#include "qp_h.h"

/* The computational routine */
int* getArrSizes(int *sz_qp_A_arg, int *sz_qp_A_res, int *sz_qp_A_iw, int *sz_qp_A_w, double y[]){
	
	qp_W_work(sz_qp_A_arg, sz_qp_A_res, sz_qp_A_iw, sz_qp_A_w);
	y[0] = (double) *sz_qp_A_arg;
	y[1] = (double) *sz_qp_A_res;
	y[2] = (double) *sz_qp_A_iw;
	y[3] = (double) *sz_qp_A_w;	
	
	qp_gradJ_work(sz_qp_A_arg, sz_qp_A_res, sz_qp_A_iw, sz_qp_A_w);
	y[4] = (double) *sz_qp_A_arg;
	y[5] = (double) *sz_qp_A_res;
	y[6] = (double) *sz_qp_A_iw;
	y[7] = (double) *sz_qp_A_w;
	
	qp_gradhT_work(sz_qp_A_arg, sz_qp_A_res, sz_qp_A_iw, sz_qp_A_w);
	y[8] = (double) *sz_qp_A_arg;
	y[9] = (double) *sz_qp_A_res;
	y[10] = (double) *sz_qp_A_iw;
	y[11] = (double) *sz_qp_A_w;
	
	qp_h_work(sz_qp_A_arg, sz_qp_A_res, sz_qp_A_iw, sz_qp_A_w);
	y[12] = (double) *sz_qp_A_arg;
	y[13] = (double) *sz_qp_A_res;
	y[14] = (double) *sz_qp_A_iw;
	y[15] = (double) *sz_qp_A_w;
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */
	int sz_qp_A_arg;
	int sz_qp_A_res;
	int sz_qp_A_iw;
	int sz_qp_A_w;
    
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
	
    /* call the computational routine */
    getArrSizes(&sz_qp_A_arg, &sz_qp_A_res, &sz_qp_A_iw, &sz_qp_A_w, outMatrix);
}
