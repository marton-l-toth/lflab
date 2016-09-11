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
//? fq1: base freq.
//? [fi]: list of FIR mini-filters (-15...15 each)
//? 0:none +:low-p -:high-p. even: symmetric odd: assym.
//? t<i>, c<i>: time, constant
//? t>=0: input-t, t<0: output-t
class SparseRF : public BoxInst {
	public:
		SparseRF(int nxy) : m_nxy(nxy), m_px(0), m_xc(0) {}
		virtual ~SparseRF() { free(m_px); free(m_xc); }
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		void ini(double ** ap);
		int m_nxy, m_nx, m_ny, m_siz, m_t, *m_xt, *m_yt;
		double *m_px, *m_py, *m_xc, *m_yc;
};

// in fq1 cfg t1 c1 ... tn cn dmp
void SparseRF::ini(double** p) {
	int dmpf, nxy = m_nxy, n_xy1 = 0, siz;
	double tmul = (double)sample_rate / p[0][0];
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
}

int SparseRF::calc(int inflg, double** inb, double** outb, int n) {
	if (!m_px) { if (n) ini(inb+1); else return 1; }
	double x, *in, *out = *outb; int iim;
	if (inflg&1) in = inb[0], iim = -1; else x = inb[0][0], in = &x, iim = 0;
	int m = m_siz-1, t = m_t;
	for (int k,i=0; i<n; i++, t++) {
		m_px[k = t&m] = in[i&iim];
		double v = 0.0;
		for (int j=0; j<m_nx; j++) v += m_xc[j] * m_px[(t-m_xt[j]) & m];
		for (int j=0; j<m_ny; j++) v += m_yc[j] * m_py[(t-m_yt[j]) & m];
		m_py[k] = out[i] = CUT300(v);
	}
	m_t = t & m; return 1;
}

#define ACCLOOP(X) for (int i=0; i<n; i++) (X); break
#define ACCLOOPZ(X) for (int i=0; i<n; i++) (X), q[i] = z; break;
#define ACC_CALC(NM, MUL) static void NM(int md, double *st, double *q, double *p0, double *p1, int n) { \
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
		default: memset(q, 0, 8*n); return; } *st = fabs(v)>1e-250 ? v : 0.0; return;       }

ACC_CALC(acc_calc_d, att2mul)
ACC_CALC(acc_calc_m, )

//? {{{!._nacF}}}
//? normalised accumulator filter
//? in:filter input
//? md: mode (0:lp, 1:lp*, 2:hp, 3:hp*) *:1-sample-delay
//? dmp: dampening 
//? rpt: repeat count
//? sprd: dmp spread: (1.0-spr)*dmp ... (1.0+spr)*dmp
//? dmp can be non-constant

class NAccFilt : public BoxInst {
	public:
		NAccFilt() : m_rpt(-1) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int m_rpt;
		double * m_v;
};

int NAccFilt::calc(int inflg, double** inb, double** outb, int n) {
	double *out = outb[0], *in = inb[0];
	int md = (int)lround(inb[1][0]) & 3, m = m_rpt, flg1 = inflg&1;
	if (m<0 && ((m = m_rpt = ivlim((int)lround(inb[3][0]),0,256)))) m_v = (double*)calloc(m_rpt,8);
	if (m<2) return (m<1) ? BOX_CP0 : (acc_calc_d(flg1+((inflg&4)>>1)+4*md, m_v, out, in, inb[2], n), 1);
	if (inflg&16) {
		double x, y, mul[n], mul2[n],
		       *pm = inb[2], *ps = inb[4], mm = -2.0*sample_length / (double)(m-1);
		int dmsk = -((inflg>>2)&1), f43 = 3+4*md;
		for (int i=0; i<n; i++) x = pm[i&dmsk], y = x*ps[i], 
			mul[i] = att2mul(x-y), mul2[i] = exp(x*mm);
		acc_calc_m(flg1+2+4*md, m_v, out, in, mul, n);
		for (int i=1; i<m; i++) { for (int j=0; j<n; j++) mul[j] *= mul2[j];
					  acc_calc_m(f43, m_v+i, out,out, mul,n); }
	} else {
		double spr = inb[4][0];
		if (spr<1e-14) {
			if (inflg&4) {
				double *pm = inb[2], mul[n]; for (int i=0; i<n; i++) mul[i] = att2mul(pm[i]);
				acc_calc_m(flg1 + 2 + 4*md, m_v, out, in, mul, n);
				for (int i=1,f43=3+4*md; i<m; i++) acc_calc_m(f43, m_v+i, out, out, mul, n);
			} else {
				double mul = att2mul(inb[2][0]);
				acc_calc_m(flg1 + 4*md, m_v, out, in, &mul, n);
				for (int i=1,f41=1+4*md; i<m; i++) acc_calc_m(f41, m_v+i, out, out, &mul, n);
		}} else {
			double div = .5 * (double)(m-1), sp2 = spr / div;
			if (inflg&4) {
				double *pm = inb[2], mul[n], spr1 = 1.0-spr;
				for (int i=0; i<n; i++) mul[i] = att2mul(pm[i]*spr1);
				acc_calc_m(flg1 + 2 + 4*md, m_v, out, in, mul, n);
				for (int i=1,f43=3+4*md; i<m; i++) acc_calc_m(f43, m_v+i, out, out, mul, n);
			} else {
				double x = inb[2][0], mul = att2mul(x*(1.0-spr)), mul2 = att2mul(x*sp2);
				// log("nacc: %.15g %.15g", mul, mul2);
				acc_calc_m(flg1 + 4*md, m_v, out, in, &mul, n);
				for (int i=1,f41=1+4*md; i<m; i++) mul *= mul2,
					acc_calc_m(f41, m_v+i, out,out, &mul,n);
	}}}
	return 1;
}

//? {{{!._accF}}}
//? accumulator filter
//? out[n] = in[n] + exp(-sample_length*dmp[n])*out[n-1]
//? dmp can be non-constant
class AccFilt : public BoxInst {
	public:
		AccFilt() : m_v(0.0) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		double m_v;
};

int AccFilt::calc(int inflg, double** inb, double** outb, int n) { 
	double x, y, *p = *outb, *q0, *q1; switch (inflg & 3) {
	case 0 : x = sample_length * **inb, y = exp(-sample_length * inb[1][0]); 
		 for (int i=0; i<n; i++) p[i] = ((m_v+=x)*=y); return 1;
	case 1 : q0 = *inb, y = exp(-sample_length * inb[1][0]); 
		 for (int i=0; i<n; i++) p[i] = ((m_v+=sample_length*q0[i])*=y); return 1;
	case 2 : x = sample_length * **inb, q1 = inb[1]; 
		 for (int i=0; i<n; i++) p[i] = ((m_v+=x)*=exp(-sample_length*q1[i])); return 1;
	case 3 : q0 = *inb, q1 = inb[1]; 
		 for (int i=0;i<n;i++) p[i]= ((m_v+=sample_length*q0[i])*=exp(-sample_length*q1[i])); return 1;
	}} // no it doesn't

static void qlh_ini(double *ed, double f, double a) {
	double cf = cos(2*M_PI*sample_length*f), cf_1 = cf-1.0, acf_1 = a*cf_1, e2 = 1.0 - cf;
	ed[0] = e2 + e2; ed[1] = (M_SQRT2 * sqrt(-cf_1*cf_1*cf_1) + acf_1) / acf_1;  }

static void qlh_calc(double *q, double *p, double *yv, const double *ed, int n, int md) {
	double x, y = yv[0], v = yv[1], e = ed[0], d = ed[1];
	switch(md) {
		case 0: x=*p; for (int i=0; i<n; i++) v*=d, y+=(v+=e*(  x -y)), q[i] = y; break;
		case 1:       for (int i=0; i<n; i++) v*=d, y+=(v+=e*(p[i]-y)), q[i] = y; break;
		case 2: x=*p; for (int i=0; i<n; i++) v*=d, y+=(v+=e*(  x -y)), q[i] = x-y; break;
		case 3:       for (int i=0; i<n; i++) v*=d, y+=(v+=e*(p[i]-y)), q[i] = p[i]-y; break;
	}
	yv[0] = CUT300(y); yv[1] = CUT300(v);
}

//? {{{!._qh1}}}
//? quick & simple highpass filter
//? in: filter input
//? par: parameter (0<=par<=.5)
//? a=1-par, b=1-2*par, out[j] = a*(in[j]-in[j-1]) + b*out[j-1]
//? normalized: amp=1.0 at (sample_rate/2.0)
class QHiFilt : public BoxInst {
	public:
		QHiFilt() : m_x1(0.0), m_y1(0.0), m_flg(4) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		double m_x1, m_y1;
		int m_flg;
};

#define QHC1(J) (z = pz[J], a = 1.0-z, b = a-z)
#define QHC2(J) (x = q[J], y1 = to[J] = b*y1 + a*(x-x1), x1 = x)
#define QHC0(J) (y1 = to[J] = b*y1)
int QHiFilt::calc(int inflg, double** inb, double** outb, int n) {
	if (!n) return 0; 
	double x, z, a, b, x1 = m_x1, y1 = m_y1, *to = outb[0], *q = inb[0], *pz = inb[1];
	switch((inflg|m_flg)&7) {
		case 4: if (fabs(*pz)<1e-99) return BOX_CP0; else m_flg &= ~4;
		case 0: QHC1(0); QHC2(0); x1 = *q; for (int i=1; i<n; i++) QHC0(i); break;
		case 5: if (fabs(*pz)<1e-99) return BOX_CP0; else m_flg &= ~4;
		case 1: QHC1(0); for (int i=0; i<n; i++) QHC2(i);  break;
		case 6: m_flg &= ~4;
		case 2: QHC1(0); QHC2(0); x1 = *q; for (int i=1; i<n; i++) QHC1(i), QHC0(i); break;
		case 7: m_flg &= ~4;
		case 3: for (int i=0; i<n; i++) QHC1(i), QHC2(i); break;
	}
	m_x1 = x1; m_y1 = y1; return 1;
}

//? {{{!._qlh}}}
//? quick & simple lowpass/highpass filter
//? in: filter input
//? fq: cutoff frq
//? amp: (approx.) amplify at fq
class QLoHiFilt : public BoxInst {
	public:
		QLoHiFilt(int f) : m_flg(f|64) { m_yv[0] = m_yv[1] = 0.0; }
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int m_flg;
		double m_ed[2], m_yv[2];
};

int QLoHiFilt::calc(int inflg, double** inb, double** outb, int n) {
	if (m_flg & 64) m_flg &= 63, qlh_ini(m_ed, inb[1][0], inb[2][0]);
	return qlh_calc(*outb, *inb, m_yv, m_ed, n, (inflg&1)+(m_flg&2)), 1; }

//? {{{!._qlh*}}}
//? sequence of quick low/high pass filters
//? in: filter input
//? fq0. fq1: first & last cutoff frq
//? n: number of simple filters (n>0:lowpass, n<0:highpass)
//? scl: freq. scale type (see map01 for description)
//? ==> .!b.misc.map01
class QLHSeqFilt : public BoxInst {
	public:
		QLHSeqFilt() : m_nf(-1) {}
		~QLHSeqFilt() { if (m_nf>0) log("qlh_free: %p", m_edyv), free(m_edyv), m_nf=-2; else if (m_nf==-2) bug("qlh* double free"); }
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int m_nf, m_flg;
		double * m_edyv;
};

int QLHSeqFilt::calc(int inflg, double** inb, double** outb, int n) {
	int nf, flg; 
	double *edyv, *p = *inb, *q = *outb;
	if (m_nf <= 0) {
		if (!m_nf || !(nf = m_nf = (int)lround(inb[3][0]))) 
			return (inflg&1) ? (p==q || (memcpy(q, p, 8*n), 1)) : (*q=*p, 0);
		flg = m_flg = (nf<0) ? (nf=-nf, 2) : 0;  m_edyv = edyv = (double*)calloc(m_nf=nf, 32);
		log("qlh_alloc: %p", m_edyv);
		Scale01 sc; sc.set_all(inb[1][0], inb[2][0], (int)lround(inb[4][0]));
		double t, ts; if (nf==1) t = ts = .5; else t = 0.0, ts = 1.0/(double)(nf-1);
		for (int i=0; i<nf; i++, edyv+=4, t+=ts) qlh_ini(edyv, sc.f(t), inb[5][0]); edyv = m_edyv;
	} else { edyv = m_edyv; flg = m_flg; nf = m_nf; }
	qlh_calc(q, p, edyv+2, edyv, n, flg+(inflg&1)); ++flg; edyv += 4;
	for (int i=1; i<nf; i++, edyv += 4) qlh_calc(q, q, edyv+2, edyv, n, flg);
	return 1;
}

void b_filt_v_init(ANode *rn) {
	qmk_box(rn, "=qhi1", QMB_ARG0(QHiFilt), 0, 2, 33, "qh1", "i*r", "in$par", "uuu3%a");
}

void b_filt_misc_init(ANode * rn) {
	qmk_box(rn, "=acc", QMB_ARG0(AccFilt), 0, 2, 1, "accF", "i*r", "in$dmp", "uuu3%a");
	qmk_box(rn, "=accN", QMB_ARG0(NAccFilt), 0, 5, 1, "nacF", "i*r", "in$md$dmp$rpt$sprd", "uuu3%a");
	qmk_box(rn, "=qlow",  QMB_ARG1(QLoHiFilt), 0, 3, 33, "qlh", "i*R*1", "in$fq$amp", "uuu3%a");
	qmk_box(rn, "=qhigh", QMB_ARG1(QLoHiFilt), 2, 3, 33, "qlh", "1");
	qmk_box(rn, "=qlh*", QMB_ARG0(QLHSeqFilt), 0, 6, 33, "qlh*", "1i*", "in$fq0$fq1$n$scl$amp");
}

void b_filt_echo_init(ANode * rn) {
	char nm[8]; memcpy(nm, "=echo01", 8); 
	qmb_arg_t qa = QMB_ARG1(SparseRF);
	qmk_box(rn, nm, qa, 1, 6, 1, "spF", "i-i:o*R*1i-", 3, "in$fq1$[fi]", 
			2, 1, "c$t", "out", "uuK3%F", 0x501, "dmp");
	for (int i=2; i<14; i++) i==10 ? (nm[5]=49, nm[6]=48) : ++nm[6],
		qmk_box(rn, nm, qa, i, 2*i+4, 1, "spF", "r1i-", 512*i+769, "dmp");
}
