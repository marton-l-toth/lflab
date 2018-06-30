#include "util.h"
#include "box0.h"

static const unsigned char minifilt_ix[66] = {
124,9, 168,9, 116,8, 160,8, 109,7, 153,7, 103,6, 147,6, 98,5, 142,5, 94,4, 138,4, 91,3, 135,3, 89,2, 133,2,
0,1, 45,2, 1,2, 47,3, 3,3, 50,4, 6,4, 54,5, 10,5, 59,6, 15,6, 65,7, 21,7, 72,8, 28,8, 80,9, 36,9 };

static const double minifilt_td[33] = { 4.0,2.0, 3.5,1.75, 3.0,1.5, 2.5,1.25, 2.0,1.0, 1.5,.75,
	1.0,.5, .5,.25, 0.0, .25,.5, .5,1.0, .75,1.5, 1.0,2.0, 1.25,2.5, 1.5,3.0, 1.75,3.5, 2.0,4.0 };

static const double minifilt9[177] = { 1.0,
        0.5,0.5,
        0.25,0.5,0.25,
        0.125,0.375,0.375,0.125,
        0.0625,0.25,0.375,0.25,0.0625,
        0.03125,0.15625,0.3125,0.3125,0.15625,
                0.03125,
        0.015625,0.09375,0.234375,0.3125,0.234375,
                0.09375,0.015625,
        0.0078125,0.0546875,0.1640625,0.2734375,0.2734375,
                0.1640625,0.0546875,0.0078125,
        0.00390625,0.03125,0.109375,0.21875,0.2734375,
                0.21875,0.109375,0.03125,0.00390625,
        0.75,0.25,
        0.5625,0.375,0.0625,
        0.421875,0.421875,0.140625,0.015625,
        0.31640625,0.421875,0.2109375,0.046875,0.00390625,
        0.2373046875,0.3955078125,0.263671875,0.087890625,0.0146484375,
                0.0009765625,
        0.177978515625,0.35595703125,0.296630859375,0.1318359375,0.032958984375,
                0.00439453125,0.000244140625,
        0.13348388671875,0.31146240234375,0.31146240234375,0.17303466796875,0.05767822265625,
                0.01153564453125,0.00128173828125,6.103515625e-05,
        0.100112915039062,0.2669677734375,0.31146240234375,0.2076416015625,0.086517333984375,
                0.0230712890625,0.00384521484375,0.0003662109375,1.52587890625e-05,
        -0.5,0.5,
        0.25,-0.5,0.25,
        -0.125,0.375,-0.375,0.125,
        0.0625,-0.25,0.375,-0.25,0.0625,
        -0.03125,0.15625,-0.3125,0.3125,-0.15625,
                0.03125,
        0.015625,-0.09375,0.234375,-0.3125,0.234375,
                -0.09375,0.015625,
        -0.0078125,0.0546875,-0.1640625,0.2734375,-0.2734375,
                0.1640625,-0.0546875,0.0078125,
        0.00390625,-0.03125,0.109375,-0.21875,0.2734375,
                -0.21875,0.109375,-0.03125,0.00390625,
        0.75,-0.25,
        0.5625,-0.375,0.0625,
        0.421875,-0.421875,0.140625,-0.015625,
        0.31640625,-0.421875,0.2109375,-0.046875,0.00390625,
        0.2373046875,-0.3955078125,0.263671875,-0.087890625,0.0146484375,
                -0.0009765625,
        0.177978515625,-0.35595703125,0.296630859375,-0.1318359375,0.032958984375,
                -0.00439453125,0.000244140625,
        0.13348388671875,-0.31146240234375,0.31146240234375,-0.17303466796875,0.05767822265625,
                -0.01153564453125,0.00128173828125,-6.103515625e-05,
        0.100112915039062,-0.2669677734375,0.31146240234375,-0.2076416015625,0.086517333984375,
                -0.0230712890625,0.00384521484375,-0.0003662109375,1.52587890625e-05 };

//? {{{!._spF}}}
//? echo filter (sparse recursive filter) / comb filter
//? in: filter input
//? fq1: base freq. (*)
//? [fi]: list of FIR mini-filters (-15...15 each)
//? 0:none +:low-p -:high-p. even: symmetric odd: assym.
//? t<i>, c<i>: time, constant
//? t>=0: input-t, t<0: output-t
//? ---
//? (*) if fq1 is negative, times are in samples, and multiplied
//? by -fq1 (fq1 = -1.0 is raw samples mode)
#define SPFUN_LS c10,c11,c12,c13,c19,c20,c21,c22,c23,c29,c30,c31,c32,c33,c39,c90,c91,c92,c93,c99
class SparseRF; typedef int (comb_fun_dt)(SparseRF*, int, double*, double*, int),
      			    (*comb_fun_t)(SparseRF*, int, double*, double*, int);
class SparseRF : public BoxInst {
	public:
		static scf_t sc_f0, sc_f1;
		SparseRF(int nxy) : BoxInst(sc_f0), m_nxy(nxy), m_px(0), m_xc(0) {}
		virtual ~SparseRF() { free(m_px); free(m_xc); }
	protected:
		static comb_fun_dt c00, SPFUN_LS;
		void ini(double ** ap);
		int m_nxy, m_nx, m_ny, m_siz, m_t, *m_xt, *m_yt;
		double *m_px, *m_py, *m_xc, *m_yc;
		comb_fun_t m_fun;
};

// in fq1 cfg t1 c1 ... tn cn dmp
void SparseRF::ini(double** p) {
	static comb_fun_t ftab[25] = {c00,c00,c00,c00,c00, SPFUN_LS};
	int dmpf, nxy = m_nxy, n_xy1 = 0, siz;
	double f1=**p, tmul = f1<0 ? -f1 : (double)sample_rate / f1;
	NAN_UNPK_8(mf, p[1], 0);
	double dmpc = p[2*nxy+2][0], 
	       dmp1 = (dmpf = (fabs(dmpc)>1e-100)) ? exp(dmpc *= -sample_length) : 1.0;
	for (int j, i=0; i<nxy; i++) j = mf_x[i]+16, j = mf_x[i] = (unsigned int)j<33u ? j : 16,
				     n_xy1 += minifilt_ix[2*j+1];
	char *blk = (char*)malloc(12*n_xy1);
	double *xc = m_xc = (double*)blk, *yc = xc + n_xy1;
	int *xt = m_xt = (int*)(blk+8*n_xy1), *yt = xt + n_xy1, tlim = 0;
	for (int i=0; i<nxy; i++) {
		int *qi, yf, df2, k = mf_x[i], n2 = minifilt_ix[2*k+1];
		double *qd, t0 = p[2+2*i][0], c = p[3+2*i][0]; 
		const double *mfv = minifilt9 + (int)minifilt_ix[2*k];
		if (t0<-1e-10) t0 = -t0, yf = 1, df2 = dmpf, qi = (yt-=n2), qd = (yc-=n2);
		else yf = df2 = 0, qi = xt, qd = xc, xt += n2, xc += n2;
		int t = (int)lround(t0 * tmul + minifilt_td[k]);
		if (t<yf) t=yf; else if (t>441000) t = 441000;
		int tlimc = t + n2; if (tlimc>tlim) tlim = tlimc;
		/*if (df2)*/ for (int j=(c*=exp(dmpc*t),00); j<n2; j++) qi[j]=t+j, qd[j]=c*mfv[j], c*=dmp1;
		//else     for (int j=0; 			 j<n2; j++) qi[j]=t+j, qd[j]=c*mfv[j];
	}
	if (xc!=yc || xt!=yt) bug("sparse/ini"); 
	m_nx = xc - m_xc, m_ny = n_xy1 - m_nx; m_yc = xc; m_yt = xt;
	for (siz = 16; siz<tlim; siz += siz);
	m_px = (double*)malloc(16*siz); m_py = m_px + siz; m_t = 0; m_siz = siz;
	for (int i=siz-tlim; i<siz; i++) m_px[i] = 0.0;
	for (int i=siz-tlim; i<siz; i++) m_py[i] = 0.0;
	m_fun = ftab[5*min_i(m_nx, 4)+min_i(m_ny,4)];
}

BX_SCALC(SparseRF::sc_f0) { SCALC_BXI(SparseRF); bxi->ini(inb+1); CALC_FW(sc_f1); }

BX_SCALC(SparseRF::sc_f1) {
	static double x; SCALC_BXI(SparseRF); double *in = inb[0], *out = outb[0];
	int iif; if ( !(iif = -(inflg&1)) ) x = *in, in = &x;
	return (*bxi->m_fun)(bxi, iif, in, out, n); }

int SparseRF::c00(SparseRF *p, int iim, double *in, double *out, int n) { *out=0.0; return 0; }

#define SPCI(Z,J) Z##t##J = p->m_##Z##t[J]
#define SPCD(Z,J) Z##c##J = p->m_##Z##c[J]
#define SPC0Y0
#define SPC0X1 int SPCI(x,0);			  double SPCD(x,0);
#define SPC0Y1 int SPCI(y,0);			  double SPCD(y,0);
#define SPC0X2 int SPCI(x,0),SPCI(x,1);		  double SPCD(x,0),SPCD(x,1);
#define SPC0Y2 int SPCI(y,0),SPCI(y,1);		  double SPCD(y,0),SPCD(y,1);
#define SPC0X3 int SPCI(x,0),SPCI(x,1),SPCI(x,2); double SPCD(x,0),SPCD(x,1),SPCD(x,2);
#define SPC0Y3 int SPCI(y,0),SPCI(y,1),SPCI(y,2); double SPCD(y,0),SPCD(y,1),SPCD(y,2);
#define SPC0X9 int nx = p->m_nx, *xt = p->m_xt; double *xc = p->m_xc;
#define SPC0Y9 int ny = p->m_ny, *yt = p->m_yt; double *yc = p->m_yc;
#define SPCP(Z,J) Z##c##J * p##Z[(t-Z##t##J)&m]
#define SPC1Y0
#define SPC1X1 v =  SPCP(x,0);
#define SPC1Y1 v += SPCP(y,0);
#define SPC1X2 v =  SPCP(x,0) + SPCP(x,1);
#define SPC1Y2 v += SPCP(y,0) + SPCP(y,1);
#define SPC1X3 v =  SPCP(x,0) + SPCP(x,1) + SPCP(x,2);
#define SPC1Y3 v += SPCP(y,0) + SPCP(y,1) + SPCP(y,2);
#define SPC1X9 v =  0.0; for (int j=0; j<nx; j++) v += xc[j] * px[(t-xt[j]) & m];
#define SPC1Y9		 for (int j=0; j<ny; j++) v += yc[j] * py[(t-yt[j]) & m];

#define SPCDEF(I,J) int SparseRF::c##I##J(SparseRF *p, int iim, double*in, double*q, int n) { \
	int t = p->m_t, m = p->m_siz-1;  p->m_t += n; double *px = p->m_px, *py = p->m_py;  SPC0X##I SPC0Y##J \
	for (int k,i=0; i<n; i++, t++) { px[k=t&m] = in[i&iim]; \
					 double SPC1X##I SPC1Y##J q[i]=py[k]=CUT300(v); } return 1; }

SPCDEF(1,0) SPCDEF(1,1) SPCDEF(1,2) SPCDEF(1,3) SPCDEF(1,9)
SPCDEF(2,0) SPCDEF(2,1) SPCDEF(2,2) SPCDEF(2,3) SPCDEF(2,9)
SPCDEF(3,0) SPCDEF(3,1) SPCDEF(3,2) SPCDEF(3,3) SPCDEF(3,9)
SPCDEF(9,0) SPCDEF(9,1) SPCDEF(9,2) SPCDEF(9,3) SPCDEF(9,9)

//? {{{!._nacF}}}
//? normalised accumulator filter
//? in:filter input
//? md: mode (0:lp, 1:lp*, 2:hp, 3:hp*) *:1-sample-delay
//? dmp: dampening 
//? rpt: repeat count
//? sprd: dmp spread: (1.0-spr)*dmp ... (1.0+spr)*dmp
//? dmp can be non-constant

#define ACCLOOP(X) for (int i=0; i<n; i++) (X); break
#define ACCLOOPZ(X) for (int i=0; i<n; i++) (X), q[i] = z; break;
#define ACC_CALC(NM, MUL) static int NM(int md, double *st, double *q, double *p0, double *p1, int n) { \
	double x, y, y1, z, v = *st; switch(md) { 					   \
		case  0: y = MUL(*p1), y1 = *p0 * (1.0-y); ACCLOOP(q[i] = v = (y*v + y1));  \
		case  1: y = MUL(*p1), y1 = 1.0 - y; ACCLOOP(q[i] = v = (y*v + y1*p0[i]));   \
		case  2: x = *p0; ACCLOOP((y = MUL(p1[i]), q[i] = v = y*v + (1.0-y)*x));      \
		case  3: ACCLOOP((y = MUL(p1[i]), q[i] = v = y*v + (1.0-y)*p0[i]));  	       \
		case  4: y = MUL(*p1), y1 = *p0 * (1.0-y); ACCLOOPZ((z = v, v = (y*v + y1)));   \
		case  5: y = MUL(*p1), y1 = 1.0 - y; ACCLOOPZ((z = v, v = (y*v + y1*p0[i])));    \
		case  6: x = *p0; ACCLOOPZ((y = MUL(p1[i]), z = v, v = y*v + (1.0-y)*x));	  \
		case  7: ACCLOOP((y = MUL(p1[i]), q[i] = v, v = y*v + (1.0-y)*p0[i]));	           \
		case  8: y = MUL(*p1), y1 = (x=*p0) * (1.0-y); ACCLOOP(q[i] = x - (v = (y*v + y1)));\
		case  9: y = MUL(*p1), y1 = 1.0 - y; ACCLOOP(q[i] = p0[i] - (v = (y*v + y1*p0[i])));\
		case 10: x = *p0; ACCLOOP((y = MUL(p1[i]), q[i] = x - (v = y*v + (1.0-y)*x)));      \
		case 11: ACCLOOP((y = MUL(p1[i]), q[i] = p0[i] - (v = y*v + (1.0-y)*p0[i])));	    \
		case 12: y = MUL(*p1), y1 = (x=*p0) * (1.0-y); ACCLOOPZ((z = x-v, v = (y*v + y1))); \
		case 13: y = MUL(*p1), y1 = 1.0 - y; ACCLOOPZ((z = p0[i]-v, v = (y*v + y1*p0[i]))); \
		case 14: x = *p0; ACCLOOPZ((y = MUL(p1[i]), z = x-v, v = y*v + (1.0-y)*x));         \
		case 15: ACCLOOPZ((y = MUL(p1[i]), z = p0[i]-v, v = y*v + (1.0-y)*p0[i]));          \
		default: memset(q, 0, 8*n); return 1; } *st = fabs(v)>1e-250 ? v : 0.0; return 1;   }

ACC_CALC(acc_calc_d, att2mul)
ACC_CALC(acc_calc_m, )

class NAccFilt : public BoxInst {
	public:
		static scf_t sc_f0, sc_f1, sc_fn;
		NAccFilt() : BoxInst(sc_f0), m_rpt(0) {}
		virtual ~NAccFilt() { if(m_rpt>1) free(m_u.p); }
	protected:
		int m_rpt;
		union { double v; double *p; } m_u;
};

BX_SCALC(NAccFilt::sc_f0) {
	SCALC_BXI(NAccFilt); int m = bxi->m_rpt = ivlim((int)lround(inb[3][0]),0,256);
	switch(m){ case 0:					CALC_FW(sc_cp0);
		   case 1:  bxi->m_u.v = 0.0;			CALC_FW(sc_f1);
		   default: bxi->m_u.p = (double*)calloc(m,8);	CALC_FW(sc_fn); }}

BX_SCALC(NAccFilt::sc_f1) { SCALC_BXI(NAccFilt); int f1 = inflg&1, md = (int)lround(inb[1][0]) & 3;
			    return acc_calc_d(f1+((inflg&4)>>1)+4*md, &bxi->m_u.v, *outb, *inb, inb[2], n); }

BX_SCALC(NAccFilt::sc_fn) {
	SCALC_BXI(NAccFilt); double *out = outb[0], *in = inb[0], *pv = bxi->m_u.p;
	int md = (int)lround(inb[1][0]) & 3, m = bxi->m_rpt, flg1 = inflg&1;
	if (inflg&16) {
		double x, y, mul[n], mul2[n],
		       *pm = inb[2], *ps = inb[4], mm = -2.0*sample_length / (double)(m-1);
		int dmsk = -((inflg>>2)&1), f43 = 3+4*md;
		for (int i=0; i<n; i++) x = pm[i&dmsk], y = x*ps[i], 
			mul[i] = att2mul(x-y), mul2[i] = exp(x*mm);
		acc_calc_m(flg1+2+4*md, pv, out, in, mul, n);
		for (int i=1; i<m; i++) { for (int j=0; j<n; j++) mul[j] *= mul2[j];
					  acc_calc_m(f43, pv+i, out,out, mul,n); }
	} else {
		double spr = inb[4][0];
		if (spr<1e-14) {
			if (inflg&4) {
				double *pm = inb[2], mul[n]; for (int i=0; i<n; i++) mul[i] = att2mul(pm[i]);
				acc_calc_m(flg1 + 2 + 4*md, pv, out, in, mul, n);
				for (int i=1,f43=3+4*md; i<m; i++) acc_calc_m(f43, pv+i, out, out, mul, n);
			} else {
				double mul = att2mul(inb[2][0]);
				acc_calc_m(flg1 + 4*md, pv, out, in, &mul, n);
				for (int i=1,f41=1+4*md; i<m; i++) acc_calc_m(f41, pv+i, out, out, &mul, n);
		}} else {
			double div = .5 * (double)(m-1), sp2 = spr / div;
			if (inflg&4) {
				double *pm = inb[2], mul[n], spr1 = 1.0-spr;
				for (int i=0; i<n; i++) mul[i] = att2mul(pm[i]*spr1);
				acc_calc_m(flg1 + 2 + 4*md, pv, out, in, mul, n);
				for (int i=1,f43=3+4*md; i<m; i++) acc_calc_m(f43, pv+i, out, out, mul, n);
			} else {
				double x = inb[2][0], mul = att2mul(x*(1.0-spr)), mul2 = att2mul(x*sp2);
				// log("nacc: %.15g %.15g", mul, mul2);
				acc_calc_m(flg1 + 4*md, pv, out, in, &mul, n);
				for (int i=1,f41=1+4*md; i<m; i++) mul *= mul2,
					acc_calc_m(f41, pv+i, out,out, &mul,n);
	}}}
	return 1;
}

//? {{{!._accF}}}
//? accumulator filter
//? out[n] = in[n] + exp(-sample_length*dmp[n])*out[n-1]
//? dmp can be non-constant
class AccFilt : public BoxInst {
	public:
		static scf_t sc_f0;
		AccFilt() : BoxInst(sc_f0), m_v(0.0) {}
	protected:
		double m_v;
};

BX_SCALC(AccFilt::sc_f0) {
	SCALC_BXI(AccFilt); double x, y, v=bxi->m_v, *p = *outb, *q0, *q1;
	switch (inflg & 3) {
		case 0 : x = sample_length * **inb, y = exp(-sample_length * inb[1][0]); 
			 for (int i=0; i<n; i++) p[i] = ((v+=x)*=y); break;
		case 1 : q0 = *inb, y = exp(-sample_length * inb[1][0]); 
			 for (int i=0; i<n; i++) p[i] = ((v+=sample_length*q0[i])*=y); break;
		case 2 : x = sample_length * **inb, q1 = inb[1]; 
			 for (int i=0; i<n; i++) p[i] = ((v+=x)*=exp(-sample_length*q1[i])); break;
		case 3 : q0 = *inb, q1 = inb[1]; 
			 for (int i=0;i<n;i++) p[i]= ((v+=sample_length*q0[i])*=exp(-sample_length*q1[i])); }
	bxi->m_v = v; return 1;
}

static void qlh_ini(double *ed, double f, double a) {
	double cf = cos(2*M_PI*sample_length*f), cf_1 = cf-1.0, acf_1 = a*cf_1, e2 = 1.0 - cf;
	ed[0] = e2 + e2; ed[1] = (M_SQRT2 * sqrt(-cf_1*cf_1*cf_1) + acf_1) / acf_1;  }

static int qlh_calc(double *q, double *p, double *yv, const double *ed, int n, int md) {
	double x, y = yv[0], v = yv[1], e = ed[0], d = ed[1];
	switch(md) {
		case 0: x=*p; for (int i=0; i<n; i++) v*=d, y+=(v+=e*(  x -y)), q[i] = y; break;
		case 1:       for (int i=0; i<n; i++) v*=d, y+=(v+=e*(p[i]-y)), q[i] = y; break;
		case 2: x=*p; for (int i=0; i<n; i++) v*=d, y+=(v+=e*(  x -y)), q[i] = x-y; break;
		case 3:       for (int i=0; i<n; i++) v*=d, y+=(v+=e*(p[i]-y)), q[i] = p[i]-y; break;
	}
	yv[0] = CUT300(y); yv[1] = CUT300(v); return 1;
}

//? {{{!._qh1}}}
//? quick & simple highpass filter
//? in: filter input
//? par: parameter (0<=par<=.5), can be non-constant
//? a=1-par, b=1-2*par, out[j] = a*(in[j]-in[j-1]) + b*out[j-1]
//? normalized: amp=1.0 at (sample_rate/2.0)
class QHiFilt : public BoxInst {
	public:
		static scf_t sc_f0;
		QHiFilt() : BoxInst(sc_f0), m_x1(0.0), m_y1(0.0), m_flg(4) {}
	protected:
		double m_x1, m_y1;
		int m_flg; // 4:cp
};

#define QHC1(J) (z = pz[J], a = 1.0-z, b = a-z)
#define QHC2(J) (x = q[J], y1 = to[J] = b*y1 + a*(x-x1), x1 = x)
#define QHC0(J) (y1 = to[J] = b*y1)
BX_SCALC(QHiFilt::sc_f0) {
	SCALC_BXI(QHiFilt);
	double x, z, a, b, x1 = bxi->m_x1, y1 = bxi->m_y1, *to = outb[0], *q = inb[0], *pz = inb[1];
	switch((inflg|bxi->m_flg)&7) {
		case 4: if (fabs(*pz)<1e-99) return BOX_CP0; else bxi->m_flg &= ~4;
		case 0: QHC1(0); QHC2(0); x1 = *q; for (int i=1; i<n; i++) QHC0(i); break;
		case 5: if (fabs(*pz)<1e-99) return BOX_CP0; else bxi->m_flg &= ~4;
		case 1: QHC1(0); for (int i=0; i<n; i++) QHC2(i);  break;
		case 6: bxi->m_flg &= ~4;
		case 2: QHC1(0); QHC2(0); x1 = *q; for (int i=1; i<n; i++) QHC1(i), QHC0(i); break;
		case 7: bxi->m_flg &= ~4;
		case 3: for (int i=0; i<n; i++) QHC1(i), QHC2(i); break;
	}
	bxi->m_x1 = x1; bxi->m_y1 = y1; return 1;
}

//? {{{!._q1p}}}
//? Single-pole filter: q1p* / q1p+ / q1p-
//? in: filter input
//? a: parameter, -1<=a<=1 (can be non-constant)
//? diff.eq.: out[j] = b*in[j] + a*out[j-1]
//? ---
//? q1p*: b=1-abs(a) (amp=1.0 @ max)
//? q1p+: b=1-a (amp=1.0 @ DC (0Hz))
//? q1p-: b=1+a (amp=1.0 @ samp.rate/2 (22050Hz))
class Q1pFilt : public BoxInst {
	public:
		static scf_t sc_f0;
		Q1pFilt(int md) : BoxInst(sc_f0), m_y1(0.0), m_flg(16+4*md) {}
	protected:
		double m_y1;
		int m_flg; // 12:md 16:cp
};

#define Q1PF(L) if (fabs(*pz)<1e-99) goto cp##L; else bxi->m_flg &= ~16
#define Q1P0(B) x=*q; for (int i=0; i<n; i++) a=pz[i],b=(B),y1=to[i]=a*y1+b* x  ; goto bye
#define Q1P1(B)       for (int i=0; i<n; i++) a=pz[i],b=(B),y1=to[i]=a*y1+b*q[i]; goto bye
BX_SCALC(Q1pFilt::sc_f0) {
	SCALC_BXI(Q1pFilt);
	double x, a, b, y1 = bxi->m_y1, *to = outb[0], *q = inb[0], *pz = inb[1];
	switch((inflg|bxi->m_flg)&31) {
		case 16: Q1PF(0);	  case  0: a = *pz; b = 1.0-fabs(a); goto cc0;
		case 17: Q1PF(1);	  case  1: a = *pz; b = 1.0-fabs(a); goto cc1;
		case 20: Q1PF(0);	  case  4: a = *pz; b = 1.0-   a   ; goto cc0;
		case 21: Q1PF(1);	  case  5: a = *pz; b = 1.0-   a   ; goto cc1;
		case 24: Q1PF(0);	  case  8: a = *pz; b = 1.0+   a   ; goto cc0;
		case 25: Q1PF(1);	  case  9: a = *pz; b = 1.0+   a   ; goto cc1;
		case 18: bxi->m_flg&=~16; case  2: Q1P0(1.0-fabs(a));
		case 19: bxi->m_flg&=~16; case  3: Q1P1(1.0-fabs(a));
		case 22: bxi->m_flg&=~16; case  6: Q1P0(1.0-   a   );
		case 23: bxi->m_flg&=~16; case  7: Q1P1(1.0-   a   );
		case 26: bxi->m_flg&=~16; case 10: Q1P0(1.0+   a   );
		case 27: bxi->m_flg&=~16; case 11: Q1P1(1.0+   a   );
		default: return RTE_BUG;
	}
cp0:	*to = *q; return 0;
cp1:	if (to!=q) memcpy(to, q, 8*n);   return 1;
cc0:	x=*q; for (int i=0; i<n; i++) to[i] = y1 = a*y1+b* x  ;  goto bye;
cc1:	      for (int i=0; i<n; i++) to[i] = y1 = a*y1+b*q[i];  goto bye;
bye:	bxi->m_y1 = y1; return 1;
}

//? {{{!._ffoa}}}
//? first order allpass - constant fractional delay
//? in: filter input
//? dly: delay (in samples, recommended: 0.1 ... 1.1)
class FOAFilt : public BoxInst {
	public:
		static scf_t sc_f0, sc_f1;
		FOAFilt() : BoxInst(sc_f0) {}
	protected:
		double m_x1, m_y1, m_z;
};

BX_SCALC(FOAFilt::sc_f0) { SCALC_BXI(FOAFilt); double x = inb[1][0];
			   bxi->m_x1 = bxi->m_y1 = 0.0; bxi->m_z = (1.0-x)/(1.0+x); CALC_FW(sc_f1); }

BX_SCALC(FOAFilt::sc_f1) {
	SCALC_BXI(FOAFilt); double x, x1=bxi->m_x1, y1=bxi->m_y1, z=bxi->m_z, *q = inb[0], *to = outb[0];
	if (inflg&1) { for (int i=0; i<n; i++) x=q[i], to[i]=y1=z*(x-y1)+x1, x1=x; }
	else { x=*q; to[0]=y1=z*(x-y1)+x1; x1=x; for (int i=1; i<n; i++) to[i]=y1=z*(x-y1)+x1; }
	bxi->m_x1 = x1; bxi->m_y1 = y1; return 1;  }

//? {{{!._qlh}}}
//? quick & simple lowpass/highpass filter
//? in: filter input
//? fq: cutoff frq
//? amp: (approx.) amplify at fq
class QLoHiFilt : public BoxInst {
	public:
		static scf_t sc_f0, sc_f1;
		QLoHiFilt(int f) : BoxInst(sc_f0), m_flg(f) {}
	protected:
		int m_flg;
		double m_ed[2], m_yv[2];
};

BX_SCALC(QLoHiFilt::sc_f0) { SCALC_BXI(QLoHiFilt); qlh_ini(bxi->m_ed, inb[1][0], inb[2][0]);
			     bxi->m_yv[0] = bxi->m_yv[1] = 0.0; CALC_FW(sc_f1); }

BX_SCALC(QLoHiFilt::sc_f1) { SCALC_BXI(QLoHiFilt); 
			     return qlh_calc(*outb,*inb, bxi->m_yv,bxi->m_ed, n, (inflg&1)+(bxi->m_flg&2)); }

//? {{{!._qlh*}}}
//? sequence of quick low/high pass filters
//? in: filter input
//? fq0. fq1: first & last cutoff frq
//? n: number of simple filters (n>0:lowpass, n<0:highpass)
//? scl: freq. scale type (see map01 for description)
//? ==> .!b.misc.map01
class QLHSeqFilt : public BoxInst {
	public:
		static scf_t sc_f0, sc_f1, sc_fn;
		QLHSeqFilt() : BoxInst(sc_f0) {}
		~QLHSeqFilt() { if (m_nf>0) free(m_edyv); }
	protected:
		int m_nf, m_flg;
		double * m_edyv;
};

BX_SCALC(QLHSeqFilt::sc_f0) {
	SCALC_BXI(QLHSeqFilt); int nf = (int)lround(inb[3][0]), flg = 0;
	if (nf<0) nf = -nf, flg = 2;  if (nf>999) nf=999;
	bxi->m_nf = nf; bxi->m_flg = flg;
	double *edyv = bxi->m_edyv = (double*)calloc(nf, 32);
	Scale01 sc; sc.set_all(inb[1][0], inb[2][0], (int)lround(inb[4][0]));
	double t, ts; sc_t psc;
	if (nf==1) t = 	    ts = 0.5, 			psc = sc_f1;
	else 	   t = 0.0, ts = 1.0/(double)(nf-1),	psc = sc_fn;
	for (int i=0; i<nf; i++, edyv+=4, t+=ts) qlh_ini(edyv, sc.f(t), inb[5][0]);
	CALC_FW(psc);
}

BX_SCALC(QLHSeqFilt::sc_f1) { SCALC_BXI(QLHSeqFilt); double *edyv = bxi->m_edyv;
			      return qlh_calc(*outb, *inb, edyv+2, edyv, n, bxi->m_flg+(inflg&1)); }

BX_SCALC(QLHSeqFilt::sc_fn) {
	SCALC_BXI(QLHSeqFilt); int nf = bxi->m_nf, flg = bxi->m_flg;
	double *p = *inb, *q = *outb, *edyv = bxi->m_edyv;
	qlh_calc(q, p, edyv+2, edyv, n, flg+(inflg&1)); ++flg; edyv += 4;
	for (int i=1; i<nf; i++, edyv += 4) qlh_calc(q, q, edyv+2, edyv, n, flg);
	return 1;
}

void b_filt_v_init(ANode *rn) {
	qmk_box(rn, "=qhi1", QMB_ARG0(QHiFilt), 0, 2, 33, "qh1", "i*r", "in$par", "uuu3%a");
	qmk_box(rn, "=q1p*", QMB_ARG1(Q1pFilt), 0, 2, 33, "q1p", "i*R*1", "in$a", "uuu3%a");
	qmk_box(rn, "=q1p+", QMB_ARG1(Q1pFilt), 1, 2, 33, "q1p", "1");
	qmk_box(rn, "=q1p-", QMB_ARG1(Q1pFilt), 2, 2, 33, "q1p", "1");
}

void b_filt_misc_init(ANode * rn) {
	qmk_box(rn, "=acc", QMB_ARG0(AccFilt), 0, 2, 1, "accF", "i*r", "in$dmp", "uuu3%a");
	qmk_box(rn, "=accN", QMB_ARG0(NAccFilt), 0, 5, 1, "nacF", "i*r", "in$md$dmp$rpt$sprd", "uuu3%a");
	qmk_box(rn, "=qlow",  QMB_ARG1(QLoHiFilt), 0, 3, 33, "qlh", "i*R*1", "in$fq$amp", "uuu3%a");
	qmk_box(rn, "=qhigh", QMB_ARG1(QLoHiFilt), 2, 3, 33, "qlh", "1");
	qmk_box(rn, "=qlh*", QMB_ARG0(QLHSeqFilt), 0, 6, 33, "qlh*", "1i*", "in$fq0$fq1$n$scl$amp");
	qmk_box(rn, "=f/dly", QMB_ARG0(FOAFilt), 0, 2, 33, "ffoa", "1i*", "in$dly");
}

void b_filt_echo_init(ANode * rn) {
	char nm[8]; memcpy(nm, "=echo01", 8); 
	qmb_arg_t qa = QMB_ARG1(SparseRF);
	qmk_box(rn, nm, qa, 1, 6, 1, "spF", "i-i:o*R*1i-", 3, "in$fq1$[fi]", 
			2, 1, "c$t", "out", "uuK3%F", 0x501, "dmp");
	for (int i=2; i<14; i++) i==10 ? (nm[5]=49, nm[6]=48) : ++nm[6],
		qmk_box(rn, nm, qa, i, 2*i+4, 1, "spF", "r1i-", 512*i+769, "dmp");
}
