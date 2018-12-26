#ifndef __qwe_pt_h__
#define __qwe_pt_h__

#define QENV(c)  pt_qenv_v[(int)c&31]
#define QENVL(c) pt_qenv_l[(int)c&31]

#define ACV_WIN(J) (0xf00007 + 16*(J&65535))

extern volatile int pt_sig_flg;
extern int pt_nullfd, 	      pt_qenv_l[32];
extern const char *pt_hello, *pt_qenv_v[32];
extern int pt_errtemp, pt_xapp_bv[];

int cfg_write(int lg);
struct cfg_ent;
void cfg_export(cfg_ent *p);
void cfg_setint(cfg_ent *p, int k);
int cfg_setstr(cfg_ent *p, const char *s);

const char ** pt_init(int ac, char ** av, int *pfd_io, int *pfd_con, int *pfd_wrk);
void pt_sig_act();
int pt_reg_prc(int pid, const char** av, int rprt);
int pt_iw_cmd_sn(int trg, const char *s, int n);
inline int pt_iocmd_sn(const char *s, int n) { return pt_iw_cmd_sn(0, s, n); }
inline int pt_wrk_cmd (const char *s, int n) { return pt_iw_cmd_sn(1, s, n); }
int pt_iocmd(char *s);
int pt_con_op(const char * arg); // -1:start -2:stop -3:restart -4:killed /:prclink
int pt_acv_op(int id, int op, const char *a1, const char *a2);
int pt_kill_pa(int flg);
int pt_show_lic();
void pt_calc_xbv();
double * pt_samp_shm(int bits);
int pt_wrk_start(int re);
int pt_io_dead(), pt_wrk_dead(), pt_acv_dead(const char *);
void errtemp_cond(const char* s);
int pt_plot(int nid, int flg, int bits, int len);
void pt_debug();

#endif // __qwe_pt_h__
