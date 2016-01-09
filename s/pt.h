#ifndef __qwe_pt_h__
#define __qwe_pt_h__

#define PT_IOP   0
#define PT_GUI   1
#define PT_DOT   2
#define PT_GPLT  3
#define PT_ACV   4
#define PT_COUNT 5
#define PT_STR "ioprc\0  gui\0    dot\0    gnuplot\0auconv\0 "

typedef int (*pt_wfun)(int,int,int);
extern volatile int pt_chld_flg;
extern int pt_cp_i2m;

#define ACV_WIN(J) (0xf00007 + 16*(J&65535))
void pt_reg(int ix, int pid, pt_wfun fun);
void pt_wait();
void pt_chld_act();
int pt_iocmd(char *s);
int pt_con_op(int x); // -1: start >0: pid
int pt_acv_op(int id, int op);
int pt_show_lic();

#endif // __qwe_pt_h__
