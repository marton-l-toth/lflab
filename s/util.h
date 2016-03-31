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
void bug(const char * fmt, ...);

char * bigblk(int n);
char * alloc_32k();

extern int sample_rate; extern double sample_length;
extern const double smalldiv63[64];

int approx_cmp(double x, double y);
int approx_cmp_v(const double *x, const double *y, int n);
double ipow(double x, int k); // k>=0
inline double ipows(double x, int k) { return k<0 ? ipow(1.0/x, -k) : ipow(x, k); }
void bug(const char * fmt, ...);
void log(const char * fmt, ...);
void log_n(const char * fmt, ...);
void log_sn(const char * pref, const char * s, int l);
void log_sortedblk(short *p, int n, bool absf, const char * pref = 0);
int u_sleep(int x);
int cfg_write();
struct cfg_ent;
void cfg_setint(cfg_ent *p, int k);
int cfg_setstr(cfg_ent *p, const char *s);

bool save_to_file(const char * fname);
void bye(int x);
unsigned int bit_rev(unsigned int n, int bits);
double * fft(double * re, double * im, int bits, bool reverse);
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
void set_fd(int *to, int from);

struct blk2k { unsigned short xp16, lo[818]; unsigned char rsrv, hi[409]; };
struct fa_writer { blk2k blk; unsigned short toc[1024]; int fd,cnt,cur,nch,id; unsigned int buf[818]; };
struct au16w_t { int ovf, rs4; unsigned char nch, oflg, vol, trg[5], cp[16]; };

const char * au_file_name(const char *dir, int dlen, int id, const char *a1, const char *a2, const char *ext);
int fa_start(fa_writer * fa, int nch);
int fa_add12(fa_writer * fa, const double * x, const double * y, int n);
int fa_end(fa_writer * fa);
void samp_stat(const double *p, int n, int k, bool dB, double dBy, double *pmin, double *pmax, double *pavg);
int is_asv_name(const char *s);  // /.../__asv.lf /.../__asv--x.lf
int coward(const char * fn);

class AReader { public: virtual int line(char * s) = 0; virtual ~AReader() {} }; // ret: 0:done <0:err

class QuickStat {
	public:
		QuickStat() : m_siz(63) {}
		int size() const { return m_siz; }
		void store(const double * p, int n);
		int cmd(const char *s);
		AReader * chk0(const char *s);
		int chk1(const double *p);
	protected:
		int m_siz;
		int m_pos[64];
		double m_val[64];
};
extern QuickStat qstat;

class Scale01 {
	public:
		static double f0(double v0, double v1, int ty, double x);
		Scale01() : m_ty(1), m_v0(0.0), m_v1(1.0), m_t0(0.0), m_t1(1.0) {}
		int ty() { return m_ty; }
		int ty_c() { return "CQhl-qc"[m_ty+3]; }
		double v0() { return m_v0; }
		double v1() { return m_v1; }
		void set_ty(int ty) { set_ty_2(ty); upd_t01(); }
		void set_v0(double v0) { m_v0 = v0; upd_t01(); }
		void set_v1(double v0) { m_v0 = v0; upd_t01(); }
		void set_all(double v0, double v1, int ty) { set_ty_2(ty); m_v0=v0; m_v1=v1; upd_t01(); }
		double f(double x);
	protected:
		void upd_t01();
		void set_ty_2(int ty);
		int m_ty;
		double m_v0, m_v1;
		double m_t0, m_t1, m_sg;
};

#endif // __qwe_util_h__
