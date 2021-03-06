/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

int qp_gradJ(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem);
void qp_gradJ_incref(void);
void qp_gradJ_decref(void);
int qp_gradJ_n_out(void);
int qp_gradJ_n_in(void);
const char* qp_gradJ_name_in(int i);
const char* qp_gradJ_name_out(int i);
const int* qp_gradJ_sparsity_in(int i);
const int* qp_gradJ_sparsity_out(int i);
int qp_gradJ_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif
