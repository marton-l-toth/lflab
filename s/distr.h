typedef double (*dfun1_t)(double);
static double df_gau(double x) { return exp(-16.0*x*x); }
static double df_exp(double x) { return exp(-16.0*fabs(x)); }
static double df_crc(double x) { return sqrt(1.0-x*x); }
static double df_hat(double x) { double y=x*x-1.0; return y*y; }
static dfun1_t df_tab[16] = { 0,0,0,df_crc, 0,df_exp,0,df_gau, df_hat,0,0,0, 0,0,0,0 };
#define DISTR_FUN(C) (df_tab[(C)&15])
#define DISTR_NF(C)  dfun1_t f = DISTR_FUN(C); int n = 256 - 4*((C)&32)
typedef struct { double ar; int ty, ix ; } rnd_dsc_t;
