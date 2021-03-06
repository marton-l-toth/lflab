#ifndef __qwe_util_h__
#define __qwe_util_h__

#include <cstdarg>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>

#include "errtab.inc"
#include "lwarr.h"
#include "uc0.h"
#include "uc1.h"

#define BVFOR_JM(X) for (unsigned int j, m = (X); j=__builtin_ffs(m), j-- > 0; m &= ~(1u<<j))

const double c0_freq = 16.35159783128737;

void log(const char * fmt, ...);
void bug(const char * fmt, ...) __attribute__ ((noreturn));

char * bigblk(int n);
char * alloc_32k();
inline char * nf_alloc(int x) { extern int nfa_siz; extern char *nfa_cur, *nf_alloc_2(int); x+=7,x&=~7;
	char*p=nfa_cur; int n=nfa_siz-x; return (n<0) ? nf_alloc_2(x) : (nfa_cur+=x, nfa_siz=n, p); }

extern int sample_rate; extern double sample_length;
extern const double smalldiv128[130], fib7s_tab[128];

inline int bitfill(int x) { return x|=(x>>16), x|=(x>>8), x|=(x>>4), x|=(x>>2), x|(x>>1); }
int next_prime17(int k);
inline double fib7s(int n) { return fib7s_tab[n&127]; }
int approx_cmp(double x, double y);
int approx_cmp_v(const double *x, const double *y, int n);
double ipow(double x, int k); // k>=0
inline double ipows(double x, int k) { return k<0 ? ipow(1.0/x, -k) : ipow(x, k); }
void log_n(const char * fmt, ...);
void log_sn(const char * pref, const char * s, int l, int exnl = 0);
int u_sleep(int x);

bool save_to_file(const char * fname);
void bye(int x);
void shuffle(int * p, int n);
int namelen(const char * s);
inline void filld(double *p, double x, int n) { 
	for (int i=0; i<n; i++) p[i] = x; }
inline int strchr_i(const char * s, int c) {
	const char * p = strchr(s, c);
	return p ? p - s : -1;
}
int lorem(char * to, int l);
int text_wid16(const char * s);

static inline double cubic(double x_orig, double x0, double v0, double d0, double x1, double v1, double d1) {
	double x = x_orig - x0;
	double l = x1 - x0;
	double xx = x * x;
	double l_1 = 1.0 / l;
	double ll_1 = l_1 * l_1;
	double vd = v1 - v0;
	return v0 + d0*x + ll_1 * (3.0*vd - l*(d0+d0+d1)) * xx +
		           l_1 * ll_1 * (-2.0*vd + l*(d1+d0)) * xx * x;
}

void wav_mono_head(int fd, int len);
void ssort(char ** pp, int len, int k);
void wrnan_seq(void * to, const int * p, int n);
void wrnan_lbl(void * to, const char * s, int n);
inline static int parse_sep(const char *s, int c) { for (int i=0; 1; i++) if (s[i]!=32) return i+(s[i]==c); }
int parse_num(void * to, const char * s);
inline static double at0f(const char *s) { double r = 0.0; parse_num(&r, s); return r; }
const char * dbl2str_s(int ix, double x);
inline void bump_dec2(char *s) { if (s[1]=='9') ++*s, s[1]='0'; else ++s[1]; }
inline void p_close(int *pfd) { if (*pfd>=0) close(*pfd), *pfd=-1; }

struct blk2k { unsigned short xp16, lo[818]; unsigned char rsrv, hi[409]; };
struct fa_writer { blk2k blk; unsigned short toc[1024]; int fd,cnt,cur,nch,id; unsigned int buf[818]; };
struct au16w_t { int ovf, rs4; unsigned char nch, oflg, vol, trg[5], cp[16]; };

int fa_start(fa_writer * fa, int nch);
int fa_add12(fa_writer * fa, const double * x, const double * y, int n);
int fa_end(fa_writer * fa);
int is_asv_name(const char *s);  // /.../__asv.lf /.../__asv--x.lf
int coward(const char * fn);

class AReader { public: virtual int line(char * s) = 0; virtual ~AReader() {} }; // ret: 0:done <0:err
int qstat_op(int op);
int qstat_cfg(int op, const char * arg);
int qstat_chk0(AReader ** ppr, const char *s);
int nd_path_uf(char *to, int id, int max); // def: node.c

class Scale01 {
	public:
		static double f0(double v0, double v1, int ty, double x);
		static void vec(double *q, double x0, double x1, int n, int ty);
		Scale01() : m_ty(1), m_x(0), m_v0(0.0), m_v1(1.0), m_t0(0.0), m_t1(1.0) {}
		int ty() { return m_ty; }   
		int x()  { return m_x; }
		int ty_c() { return "CQhl-qc"[m_ty+3]; }
		double v0() { return m_v0; }
		double v1() { return m_v1; }
		void set_ty(int ty) { set_ty_2(ty); upd_t01(); }
		void set_x(int xx) { m_x = xx; }
		void set_v0(double v0) { m_v0 = v0; upd_t01(); }
		void set_v1(double v0) { m_v0 = v0; upd_t01(); }
		void set_all(double v0, double v1, int ty) { set_ty_2(ty); m_v0=v0; m_v1=v1; upd_t01(); }
		void set_all_x(double v0, double v1, int ty, int xx) { set_all(v0,v1,ty); m_x = xx; }
		double f(double x);
	protected:
		void upd_t01();
		void set_ty_2(int ty);
		int m_ty, m_x;
		double m_v0, m_v1;
		double m_t0, m_t1, m_sg;
};

#endif // __qwe_util_h__
