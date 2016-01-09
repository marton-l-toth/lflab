#ifndef __qwe_nz_h__
#define __qwe_nz_h__

#define NZ_N_FUN0 10u
#define NZ_N_MAP0 10u

typedef void (*nz_fun0_t)(double*,int);
extern void nz_init();

int nz_bv_c(unsigned int * to, double v, int n);
int rndbits(int * st2, int n); // st2[0]:v st2[1]:n
void mk_gwn(double * to, int n);

extern nz_fun0_t nz_fun0[NZ_N_FUN0], nz_map0[NZ_N_MAP0];

#endif // __qwe_nz_h__
