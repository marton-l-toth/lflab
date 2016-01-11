#ifndef __qwe_mx_h__
#define __qwe_mx_h__

#define MX_L_WIN(X) (0xe00007 + 16*(X))
#define MXLF_WIN 1

class BoxInst; class BoxGen; class BoxModel;
struct au16w_t; struct fa_writer;

int mx_mkroot();
int mx_mkctl(BoxGen * bx);
int mx_mklive(BoxInst * bxi);
int mx_clear(int ix);
int mx_del(int ix);
// arg: tlim vlim x[1] ... x[ni-1]
int mx_add_filter(int trgi, BoxModel * mdl, int ni, double * arg, int osel);
int mx_add_box(int trgi, BoxInst * bxi, const char * updnnino, const double * arg, int osel, int delay = 0);
int mx_calc(int ix, double *to1, double *to2, int n, int f);
int mx_r_isemp(int ix);
int mx_c_add(int ci, int bi, int ky);
int mx_c_stop(int ci, int ky, int f);
int mx_c_dump_keys(char * to, int ci);
int mx_c_bpm_ugly_hack(int ci, int bp10m);
int mx_c_unlink(int ci);
int mx_tr_add(int bi, BoxGen * bx);
void mx_tr_rm(int tri);
unsigned char * mx_l_dat(int li); // uchar[224]
int mx_l_op(int li, int ix, int val); // (ix/-1,ff,1)-closewin (ix,ff,2)-boxrm (ix,ff,3)-getflg
int mx_au16_cfg(au16w_t * to, int nch, const char * s);
int mx_calc_int(int ix, short * to, au16w_t * cfg, fa_writer * fa, int n);

#endif // __qwe_mx_h__
