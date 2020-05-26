#include "mex.h"
#include "qp_W.h"
#include "qp_gradJ.h"
#include "qp_gradhT.h"
#include "qp_h.h"

/* The computational routine */
void getSparsity(double y[]){
	const int* test_out;
	int n_row; 
	int n_col; 
	int nnz;
	
	test_out = qp_W_sparsity_out(0);	
	n_row = *(test_out+0); 
	n_col = *(test_out+1); 	
	nnz   = *(test_out+2+n_col);
	y[0] = (double) n_row;
	y[1] = (double) n_col;
	y[2] = (double) nnz;
	
	test_out = qp_gradJ_sparsity_out(0);	
	n_row = *(test_out+0); 
	n_col = *(test_out+1); 	
	nnz   = *(test_out+2+n_col);
	y[3] = (double) n_row;
	y[4] = (double) n_col;
	y[5] = (double) nnz;
	
	test_out = qp_gradhT_sparsity_out(0);	
	n_row = *(test_out+0); 
	n_col = *(test_out+1); 	
	nnz   = *(test_out+2+n_col);
	y[6] = (double) n_row;
	y[7] = (double) n_col;
	y[8] = (double) nnz;
	
	test_out = qp_h_sparsity_out(0);	
	n_row = *(test_out+0); 
	n_col = *(test_out+1); 	
	nnz   = *(test_out+2+n_col);
	y[9] = (double) n_row;
	y[10] = (double) n_col;
	y[11] = (double) nnz;
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
    plhs[0] = mxCreateDoubleMatrix(3,4,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
	
    /* call the computational routine */
    getSparsity(outMatrix);
}
