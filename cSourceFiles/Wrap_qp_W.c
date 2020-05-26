#include <stdlib.h>
#include <stdio.h>
#include "Wrap_qp_W.h"
#include "qp_W.h"

void Wrap_qp_W(double xu[], double p[], double y[]){
	// Declare all variables for ANSI C
const double* argument[42]; 
double buff[19600]; 
double* result[41]; 
int iw[141]; 
double w[6958]; 
int n_row = 140; 
int n_col = 140; 
int nnz = 954; 
int arr_off[140]; 
int arr_nnz[954]; 
	
    int mem;
    int n_out = 0;
    const int* test_out;
    int ii;
    int nnz_count = 0;
    int kk;
	
	// Define and declare working arrays
	argument[0] = xu;
    argument[1] = p;
	result[0] = buff;
	mem = 0;

	// Calc numerical result
	qp_W(argument, result, iw, w, mem);

	// Decompose CCS
	test_out = qp_W_sparsity_out(n_out);

	for(ii=0; ii<n_col; ii++){
		arr_off[ii] = *(test_out+3+ii)-*(test_out+2+ii);
	}

	for(ii=0; ii<nnz; ii++){
		arr_nnz[ii] = *(test_out+2+n_col+1+ii);
	}

	for(ii = 0; ii < n_col; ii++){
		if(arr_off[ii]){
			for(kk = 0; kk<arr_off[ii]; kk++){
				y[ii*n_row+arr_nnz[nnz_count]] = *(result[0]+nnz_count);
				nnz_count++;
			}
		}
	}
}



