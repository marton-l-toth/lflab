#ifndef __qwe_pt_h__
#define __qwe_pt_h__

#define PT_IOP   0
#define PT_GUI   1
#define PT_DOT   2
#define PT_GPLT  3
#define PT_ACV   4
#define PT_COUNT 5
#define PT_STR "ioprc\0  gui\0    dot\0    gnuplot\0auconv\0 "

#define QENV(c)  pt_qenv_v[(int)c&31]
#define QENVL(c) pt_qenv_l[(int)c&31]

#define ACV_WIN(J) (0xf00007 + 16*(J&65535))

typedef int (*pt_wfun)(int,int,int);
extern volatile int pt_chld_flg;
extern int pt_nullfd, 	      pt_qenv_l[32];
extern const char *pt_hello, *pt_qenv_v[32];

int cfg_write(int lg); 
struct cfg_ent; 
void cfg_setint(cfg_ent *p, int k); 
int cfg_setstr(cfg_ent *p, const char *s); 

const char ** pt_init(int ac, char ** av, int *pfd_io, int *pfd_con);
void pt_reg(int ix, int pid, pt_wfun fun);
void pt_wait();
void pt_chld_act();
int pt_iocmd_sn(const char *s, int n);
int pt_iocmd(char *s);
int pt_con_op(const char * arg); // -1:start -2:stop -3:restart -4:killed /:prclink
int pt_acv_op(int id, int op, const char *a1, const char *a2);
int pt_kill_pa(int flg);
int pt_show_lic();
void pt_debug();

#endif // __qwe_pt_h__
