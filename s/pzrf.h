#ifndef __qwe_pzrf_h__
#define __qwe_pzrf_h__

#include "box0.h"

#define RECF_SB 0
#define RECF_CHEB 1
#define RECF_BSTOP 2
#define RECF_BIQUAD 3
#define RECF_HARM 4
#define RECF_ARC 5

#define QUICK_PZFILT(NM) class NM : public PZFInst { \
		public: virtual int mk_filter(double ** inb); }; \
	int NM::mk_filter(double **inb)

#define QUICK_PZFILT_1(NM) class NM : public PZFInst { \
	protected: int m_arg; \
	public: NM(int arg) : m_arg(arg) {} \
		virtual int mk_filter(double ** inb); }; \
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

class PZFInst : public BoxInst {
	public:
		static scf_t sc_ini, sc_one, sc_seq, sc_par;
		PZFInst() : BoxInst(sc_ini), m_t(0), m_damp(0.0), m_blk(0) {}
		virtual ~PZFInst();
	protected:
		virtual int mk_filter(double ** inb) = 0;
		void set_n(int n);
		void damp();
		void rf_debug(int n);

		int m_n_bq, m_t;
		double m_damp;
		char * m_blk;
		RECF_ABXY * m_ab;
};

#endif // __qwe_pzrf_h__
