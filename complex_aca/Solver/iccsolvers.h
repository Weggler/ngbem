int solve_sym_system_cg_mv(long nEq, void (*mv)(double*,double*,void*), void* data, double* rhs, double* solution, double eps, long maxIter) ;
int solve_system_gmres_mv(long nEq, void (*mat_vec)(double*,double*,void*), void* data, void (*pre_con)(double*,double*,void*), void* pdata, double* rhs, double* solution, double eps, long maxIter, void (*prn_fnc)(int,double)) ;

