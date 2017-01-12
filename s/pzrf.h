#ifndef __qwe_pzrf_h__
#define __qwe_pzrf_h__

#include "box0.h"

#define RECF_SB 0
#define RECF_CHEB 1
#define RECF_BSTOP 2
#define RECF_BIQUAD 3
#define RECF_HARM 4
#define RECF_ARC 5

#define QUICK_APZFILT(C,NM) class NM : public C##PZFInst { \
		public: virtual int mk_filter(double ** inb); }; \
	int NM::mk_filter(double **inb)

#define QUICK_APZFILT_1(C,NM) class NM : public C##PZFInst { \
	protected: int m_arg; \
	public: NM(int arg) : m_arg(arg) {} \
		virtual int mk_filter(double ** inb); }; \
	int NM::mk_filter(double **inb)

#define QUICK_PZFILT(NM)    QUICK_APZFILT(S,NM)
#define QUICK_PZFILT_1(NM)  QUICK_APZFILT_1(S,NM)
#define QUICK_PPZFILT(NM)   QUICK_APZFILT(P,NM)
#define QUICK_PPZFILT_1(NM) QUICK_APZFILT_1(P,NM)

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

class PZFInst : public BoxInst {
	public: 
		PZFInst() : m_n_bq(-1), m_t(0), m_damp(0.0), m_blk(0) {}
		virtual ~PZFInst();
	protected:
		virtual int mk_filter(double ** inb) = 0;
		void set_n(int n);
		int ini(double ** inb);
		void damp();
		void rf_debug(int n);

		int m_n_bq, m_t;
		double m_damp;
		char * m_blk;
		RECF_ABXY * m_ab;
};

class SPZFInst: public PZFInst { public: virtual int calc(int inflg, double** inb, double** outb, int n); };
class PPZFInst: public PZFInst { public: virtual int calc(int inflg, double** inb, double** outb, int n); };

#endif // __qwe_pzrf_h__
