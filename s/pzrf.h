#ifndef __qwe_pzrf_h__
#define __qwe_pzrf_h__

#include "box0.h"

#define RECF_SB 0
#define RECF_CHEB 1
#define RECF_BSTOP 2
#define RECF_BIQUAD 3
#define RECF_HARM 4
#define RECF_ARC 5

#define QUICK_PZFILT(NM) struct NM : PZFInst { static scf_t sc_ini; int mk_filter(double**); }; \
	BX_SCALC(NM::sc_ini) { SCALC_BXI(NM); return ini2(abxi,inflg,inb,outb,n,bxi->mk_filter(inb+1)); } \
	int NM::mk_filter(double **inb)

class RECF_PZItem {
	        public: double eval(double re, double im) { return c * eval1(re, im); }
			double eval1(double re, double im);
			int eql(double wfq, double wid, double amp);
			void spp(double v, double a) { a *= M_PI_2; pr = -v*sin(a); pi = v*cos(a); }
			void spz(double v, double a) { a *= M_PI_2; zr = -v*sin(a); zi = v*cos(a); }
			double c, zr, zi, pr, pi;
};
 
struct RECF_ABXY { double x21[2], y21[2], rsrv, a0, a1, a2, b1, b2; };

void pzrf_show_last();
double fq_warp(double fq);
double fq_warp2(double fq2);
void rfpz_transform(RECF_ABXY * to, RECF_PZItem * pz, int n, int flg = 1);

struct PZFInst : BoxInst_B1 {
	static scf_t sc_one, sc_seq, sc_par;
	static int ini2(BoxInst * abxi, int inflg, double** inb, double** outb, int n, int r);
	void set_n_d(int n, double d=0.0);
	void dcy();
	void rf_debug(int n);

	int m_n_bq, m_t;
	double m_dcy;
	RECF_ABXY * m_ab;
};

#endif // __qwe_pzrf_h__
