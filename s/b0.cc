#include "util.h"
#include "box0.h"
#include "glob.h"
#include "nz.h"
#include "cfgtab.inc"

#define IMPBOX(NM,NA,NT) class NM : public BoxInst { \
	public: NM() : m_t(0) {} \
		virtual int calc(int inflg, double ** inb, double ** outb, int n); \
	        double m_t, m_arg[NA]; int m_tx[NT]; }; \
int NM::calc(int inflg, double ** inb, double ** outb, int n)

#define WVBOX(NM) class NM : public BoxInst { \
	public: NM() : m_phs(NAN) {} \
		virtual int calc(int inflg, double ** inb, double ** outb, int n); \
	protected: double m_phs; }; \
	int NM::calc(int inflg, double ** inb, double ** outb, int n)

#define WVBOX1(NM, X, MUL) WVBOX(NM) { \
	double *p = *inb, *q = *outb, x = (m_phs==m_phs)?m_phs:(MUL*inb[1][0]), y, mul = (MUL)*sample_length;\
	if (inflg&1) { for (int i=0; i<n; i++) y=(X), x += mul*p[i], q[i] = y, x>MUL && (x-=MUL); } \
	else { double d=mul*p[0]; for (int i=0; i<n; i++) q[i]=(X), x+=d, x>MUL && (x-=MUL); }\
	m_phs = x; return 1; }

#define WVBOX2(NM, X, MUL) WVBOX(NM) { \
	double *p = *inb, *q = *outb, x = (m_phs==m_phs)?m_phs:(MUL*inb[2][0]), y, mul = (MUL)*sample_length;\
	PREP_INPUT(y, 1); \
	if (inflg&1) { for (int i=0; i<n; i++) y=(X), x += mul*p[i], q[i]=y, x>MUL && (x-=MUL); }   \
	else { double d=mul*p[0]; for (int i=0; i<n; i++) q[i]=(X), x+=d, x>MUL && (x-=MUL); }\
	m_phs = x; return 1; }

//? {{{!._0}}}
//? this box has no inputs, and outputs constant zero.
//? {{{!._2}}}
//? this box simply copies its input to its output.
STATELESS_BOX_0(ZeroBox) { **outb = 0.0; return 0; }
STATELESS_BOX_0(CopyBox) { return BOX_CP0; }

//? {{{!._m01}}}
//? Interval mapping: maps interval [0..1] to [y0..y1]
//? x - input
//? y0, y1 - interval endpoints (output for x=0 and x=1)
//? scl - scale type
//? ---
//? scl = -3: 1 1/8 1/27 1/64 (rev. cub.)
//? scl = -2: 1 1/4 1/9 1/16 (rev. quad.) (*)
//? scl = -1: 1 1/2 1/3 1/4 (harmonic)
//? scl = 0: 1, 2, 4, 8 (logarithmic) (*)
//? scl = 1: 1, 2, 3, 4 (linear)
//? scl = 2: 1, 4, 9, 16  (quadratic) (*)
//? scl = 3: 1, 8, 27, 64 (cubic)
//? (*): y0 and y1 must have the same sign
//? ---
//? Note: this "from-to-scale" scheme is used in several
//? builtin boxes, and in instrument configuration
//? For non-linear mapping (scl!=1) having non-constant
//? y0 and/or y1 input increases CPU usage.
STATELESS_BOX_0(Map01Box) {
	double *to = outb[0];
	if (!(inflg&15)) return *to = Scale01::f0(inb[1][0], inb[2][0], (int)lround(inb[3][0]), **inb), 0;
	if (inflg & 14) {
		PREP_INPUT(x, 0); PREP_INPUT(v0, 1); PREP_INPUT(v1, 2);
		if (inflg&8) for (int i=0; i<n; i++) to[i] = Scale01::f0(v0p[i&v0msk], v1p[i&v1msk],
								         (int)lround(inb[3][i]), xp[i&xmsk]);
		else for (int i=0, ty = (int)lround(inb[3][0]); i<n; i++) to[i] = Scale01::f0(v0p[i&v0msk],
					v1p[i&v1msk], ty, xp[i&xmsk]);
	} else {
		Scale01 sc; sc.set_all(inb[1][0], inb[2][0], (int)lround(inb[3][0]));
		double * px = inb[0]; for (int i=0; i<n; i++) to[i] = sc.f(px[i]);
	}
	return 1;
}

//? {{{!._m01x}}}
//? Multiple interval mapping
//? [sc] - list of scales (-3 ... 3) as for simple map01
//? x<i>  - mapping input (exp. 0 <= x<i> <= 1)
//? fr<i> - from (output for x<i>=0)
//? to<i> - to   (output for x<i>=1)
//? for each <i>, the 0..1 interval (x<i>) is mapped to
//? fr<i>..to<i>, with the <i>th element of [sc] selecting
//? the scale type.
//? if [sc] is shorted than expected or non-list, the default
//? scale type is 0 (logarithmic)
//? ---
//? Unless you want non-constant scale types, you can replace
//? multiple "map01" boxes with one instance of this box, 
//? simplifying graph box layout.
STATELESS_BOX_1(Map01VBox) {
	NAN_UNPK_8(scl, inb[0], 0);
	int oflg = 0;  inflg >>= 1;  inb++;
	for (int i=0, nmap=m_arg; i<nmap; i++, inflg>>=3, inb+=3) {
		int ty = scl_x[i], flg = inflg & 7;
		double * to = outb[i];
		if (!flg) { *to = Scale01::f0(inb[1][0], inb[2][0], ty, **inb); }
		else if (oflg |= (1<<i), flg==1) {
			Scale01 sc; sc.set_all(inb[1][0], inb[2][0], ty);
			double * px = inb[0]; for (int i=0; i<n; i++) to[i] = sc.f(px[i]); }
		else {  PREP_INPUT(x, 0); PREP_INPUT(f, 1); PREP_INPUT(t, 2);
			for (int i=0;i<n;i++) to[i] = Scale01::f0(fp[i&fmsk], tp[i&tmsk], ty, xp[i&xmsk]); }}
	return oflg;
}

STATELESS_BOX_0(VersionBox) { **outb = (double)v_major + 0.01 * (double)v_minor; return 0; }

//? {{{!._1}}}
//? simple impulse (first sample is 1.0, the rest are all zero)
//? HINT: since a simple impulse has a flat frequency graph, it is
//? ideal for testing the frequency response of filters.
IMPBOX(Impulse1, 0, 0) { return m_t ? (**outb=0.0, 0) 
			   	 : (**outb=1.0, memset(*outb+1, 0, 8*n-8), m_t=1, 1); }
//? {{{!._trI}}}
//? triangle impulse (up and dn are in seconds)
//? ---
//? HINT: Because up, down and amplitude are cleary visible
//? on a plot display, triangle impulses are ideal for checking
//? if scales in wrap/shadow wrap are working as expected.
IMPBOX(TriangImp, 3, 2) {
	int i, n2, t = m_t;
	if (!t) m_tx[0] = sec2samp(inb[0][0]), m_arg[0] = 1.0 / (double)m_tx[0], m_arg[2] = 0.0,
		m_tx[1] = sec2samp(inb[1][0]), m_arg[1] = 1.0 / (double)m_tx[1], m_tx[1] += m_tx[0];
	if (t>=m_tx[1]) return **outb=0.0, 0;
	double d, v = m_arg[2], *q = outb[0];
	if (t<m_tx[0]) {d=m_arg[0]; for(i=0,n2=min_i(n,m_tx[0]-t); i<n2; i++) q[i]=(v+=d); t+=n2;n-=n2;q+=n2;}
	if (t<m_tx[1]) {d=m_arg[1]; for(i=0,n2=min_i(n,m_tx[1]-t); i<n2; i++) q[i]=(v-=d); t+=n2;n-=n2;q+=n2;}
	if (n) memset(q, 0, 8*n);
	m_t = t; m_arg[2] = v; return 1;
}

//? {{{!._gsI}}}
//? Gaussian impulse
//? wid - total width, out@(t=wid/2) = 1.0
//? bits - precision, out@(t=0) = out@(t=wid) = 2^(-bits-1)
IMPBOX(GaussImp, 3, 1) {
	double y,d,dd,*q=*outb;   int t = m_t; m_t += n;
	if(!t){ if (!n) return 0;
		int hw = sec2samp(.5 * **inb);
		if (hw<1) return m_tx[0]=0, *q=1.0, memset(q+1,0,8*n-8), m_t = 1;
		double a = (double)hw, c1 = -M_LN2*inb[1][0], c = c1/(a*a);
		y  = m_arg[0] = exp(c1);
		d  = m_arg[1] = exp(c*(1.0-a-a));
		dd = m_arg[2] = exp(c+c);
		m_tx[0] = 2*hw+1;
	} else {
		int n2 = m_tx[0] - t;
		if (n2<1) return *q=0.0, 0;
		if (n2<n) memset(q+n2, 0, 8*(n-n2)), n = n2;
		y = m_arg[0]; d = m_arg[1]; dd = m_arg[2];
	}
	for (int i=0; i<n; i++) q[i] = y, y *= d, d *= dd;
	m_arg[0] = y; m_arg[1] = d; return 1;
}

//? {{{!._w}}}
//? simple waves (pt - pulse train, sin, sawtooth)
//? fq: frequency, phs: start phase (0...1)
//? Note: for sound sin and saw is recommended (out: -1...1)
//? for control (with map01) use sin01 and saw01 (out: 0...1)
#define PTLOOP(X) for (int i=0; i<n; i++) y = (phs>=1.0) ? (phs-=floor(phs), 1.0) : 0.0, phs+=(X), q[i] = y;
WVBOX(PulseTrain) {
	double *q = outb[0], *fqp = inb[0], y, phs = (m_phs==m_phs)?m_phs:inb[1][0];
	if (inflg&1) { PTLOOP(fqp[i]*sample_length); }
	else { double x = fqp[0]*sample_length; PTLOOP(x); }
	m_phs = phs; return 1;
}

double gpt_f(double x, double y) {
	return x+=.5, x-=floor(x), 0.8*(exp(-16*y*x*x)+exp(-16*y*(x-1)*(x-1))); }

//? {{{!._w1}}}
//? waves with a parameter (gpt - Gaussian pulse train, square)
//? fq: frequency, phs: start phase (0...1), k: shape parameter
WVBOX1(SineWave, sin(x), 2.0*M_PI);
WVBOX1(Sine01Wave, .5+.5*sin(x), 2.0*M_PI);
WVBOX1(Saw01Wave, x, 1.0);
WVBOX1(SawTWave, x+x-1.0, 1.0);
WVBOX2(GPlsTrain, gpt_f(x, yp[i&ymsk]), 1.0);
WVBOX2(SqWave, (x<yp[i&ymsk]) ? 1.0 : 0.0, 1.0);

//? {{{!._nz}}}
//? Simple noise generator
//? ty - noise type
//? 0:normal(Gaussian) 1:linear 2:2*lin 3:3*lin 4:lin^2 5:lin^3
//? 6:lin^4 7:lin^6 8:exponential (signed) 9: exp. (unsigned)
//? {{{!._nz2}}}
//? Two simple noise generators, multiplied
//? ty1, ty2 - noise types
//? 0:normal(Gaussian) 1:linear 2:2*lin 3:3*lin 4:lin^2 5:lin^3
//? 6:lin^4 7:lin^6 8:exponential (signed) 9: exp. (unsigned)

class NZ0Box : public BoxInst {
	public: NZ0Box() : m_fun(0) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected: nz_fun0_t m_fun; 
};

class NZ0m2Box : public BoxInst {
	public: NZ0m2Box() : m_fun1(0), m_fun2(0) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected: nz_fun0_t m_fun1, m_fun2; 
};

int NZ0Box::calc(int inflg, double** inb, double** outb, int n) {
	if (!m_fun) { int k = (int)lround(**inb);
		      m_fun = nz_fun0[(unsigned int)k < NZ_N_FUN0 ? k : 0]; }
	(*m_fun)(*outb, n); return 1; }

int NZ0m2Box::calc(int inflg, double** inb, double** outb, int n) {
	if (!m_fun1) {
		int k = (int)lround(**inb), m = (int)lround(*inb[1]);
		m_fun1 = nz_fun0[(unsigned int)k < NZ_N_FUN0 ? k : 0];
		m_fun2 = nz_fun0[(unsigned int)m < NZ_N_FUN0 ? m : 0];
	}
	double tbuf[n], *q = *outb;
	(*m_fun1)(q, n); (*m_fun2)(tbuf, n);
	for (int i=0; i<n; i++) q[i] *= tbuf[i]; return 1;
}

//? {{{!._rpt}}}
//? Randomized pulse train
//? fq1 - base fq
//? fnz - fq noise -- fq is fq1 * (1.0 + fnz*gwn())
//? vnz - vol noise -- imp vol is (1.0 - vnz) + vnz*gwn()
//? (gwn() is approx. Gaussian white noise (-1.0 ... 1.0) )

class RndPt : public BoxInst {
	public: RndPt() : m_acc(1.00000001) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected: double m_acc, m_fmul;
};

int RndPt::calc(int inflg, double** inb, double** outb, int n) {
	int i; double d, fq, v, *fqp = *inb, *q = *outb, rnd[2], fmul = m_fmul, acc = m_acc;
	PREP_INPUT(fnz, 1); PREP_INPUT(vnz, 2);
	if (inflg&1) for (i=0; i<n; i++) q[i] = (acc<1.0) ? 0.0 : (mk_gwn(rnd, 2), v = vnzp[i&vnzmsk],
			   v = 1.0-v+v*rnd[0], fmul = sample_length*(1.0+rnd[1]*fnzp[i&fnzmsk]), acc-=1.0, v),
			acc += fqp[i] * fmul;
	else for (i=0,fq=*fqp,d=fq*m_fmul; i<n; i++) q[i] = (acc<1.0) ? 0.0 : (mk_gwn(rnd,2),
			    v = vnzp[i&vnzmsk], v = 1.0-v+v*rnd[0],
			    d = fq*(fmul=sample_length*(1.0+rnd[1]*fnzp[i&fnzmsk])), acc-=1.0, v), acc += d;
	m_fmul = fmul; m_acc = acc;
	return 1;
}

STATELESS_BOX_1(DebugBox) { 
	log_n("debug_box[%d]: n=%d", m_arg, n);
	for (int i=0,m=1; i<m_arg; i++,m+=m) {
		log_n(", i%d=", i);
		if (inflg&m) {  memcpy(outb[i], inb[i], 8*n); log_n("(");
				for (int n2 = min_i(n, 3), j=0; j<n2; j++) log_n(" %.15g"+!j, inb[i][j]);
				if (n>3) log_n(" ... %.15g"+4*(n==4), inb[i][n-1]);  log_n(")");
		} else       {  log_n("%.15g", outb[i][0] = inb[i][0]); }}
	log(""); return inflg;
}

//? {{{!._rqo}}}
//? 1/r^2 oscillator [-1 ... 1]
//? f0, f1: force coeff ("charge") at -1 and 1
//? dmp: dampening
//? vl: speed limit (1/sample)
//? samp: oversampling (0: none 1:2x 2:3x 3:5x)
//? (more oversampling means higher CPU usage,
//? but also better quality and a greater range of
//? parameters where the oscillator is stable)
//? samp is expected to be constant
class RQOsc : public BoxInst {
	public:
		RQOsc() : m_y(0.0), m_v(0.0) {}
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		double m_y, m_v;
};

#define RQO_STEP(T) d0=1.0+y; d1=1.0-y; f = T*(f0/(d0*d0) - f1/(d1*d1)); \
	if ((v+=f)>0.0) { if (v> vl) v = vl; if ((y+=T*v)> .99999) y = .99999; }       \
	else            { if (v<-vl) v =-vl; if ((y+=T*v)<-.99999) y =-.99999; } 

#define RQO_LH(EL) for (int i=0; vl=vlp[i&vlmsk], f0=f0p[i&f0msk], f1=f1p[i&f1msk], i<n; i++, v *= (EL))

#define THIRD .3333333333333333
#define RQO0_CL(EL) RQO_LH(EL) { RQO_STEP(1.0) q[i] = y; }
#define RQO1_CL(EL) RQO_LH(EL) { RQO_STEP(.5) ya = y; RQO_STEP(.5); q[i] = .5*(ya+y); }
#define RQO2_CL(EL) RQO_LH(EL) { RQO_STEP(THIRD) ya = y; RQO_STEP(THIRD) ya += y; \
				 RQO_STEP(THIRD) q[i] = THIRD*(ya+y); }
#define RQO3_CL(EL) RQO_LH(EL) { RQO_STEP(.2) ya = y; RQO_STEP(.2) ya+= y; RQO_STEP(.2) ya+= y; \
				 RQO_STEP(.2) ya+= y; RQO_STEP(.2) q[i] = .2*(ya+y); }

int RQOsc::calc(int inflg, double** inb, double** outb, int n) {
	PREP_INPUT(f0, 0); PREP_INPUT(f1, 1); PREP_INPUT(vl, 3);
	int ty = (inflg&4) + ((int)lround(inb[4][0]) & 3);
	double ya, d0, d1, f, f0, f1, vl, el, *pel = inb[2], y = m_y, v = m_v, *q = outb[0];
	switch(ty) {
		case 0: el = exp(-sample_length*pel[0]); RQO0_CL(el); break;
		case 1: el = exp(-sample_length*pel[0]); RQO1_CL(el); break;
		case 2: el = exp(-sample_length*pel[0]); RQO2_CL(el); break;
		case 3: el = exp(-sample_length*pel[0]); RQO3_CL(el); break;
		case 4: RQO0_CL(exp(-sample_length*pel[i])); break;
		case 5: RQO1_CL(exp(-sample_length*pel[i])); break;
		case 6: RQO2_CL(exp(-sample_length*pel[i])); break;
		case 7: RQO3_CL(exp(-sample_length*pel[i])); break;
	}
	m_y = CUT300(y); m_v = CUT300(v); return 1; 
}

inline int ixround(double x, int n) { int r = (int)lround(x);
	return ((unsigned int)r<(unsigned int)n) ? r : (r<0?0:n-1); }

//? {{{!._sel}}}
//? this box selects the sel-th of the inputs in0...in<n-1>
//? (sel is rounded to nearest integer and limited to 0...n-1)
STATELESS_BOX_1(Sel1Box) {
	double *pv, *q = *outb, *pi = inb[m_arg];
	int ixm = 1<<m_arg; 
	if (!(inflg&ixm)) {
		int j = ixround(*pi, m_arg), r = (inflg>>j)&1;
		if (!r) *q = *inb[j]; else if ((pv=inb[j])!=q) for (int i=0;i<n;i++) q[i]=pv[i];
		return r;
	}
	double con[m_arg]; BVFOR_JM((ixm-1)&~inflg) con[j] = inb[j][0];
	for (int j,i=0; i<n; i++) j    = ixround(pi[i], m_arg), 
				  q[i] = (inflg&(1<<j)) ? inb[j][i] : con[j];
	return 1;
}

//? {{{!._sel2}}}
//? this box selects the sel-th group of inputs, copying
//? a<sel>, b<sel>, ... z<sel> to a, b, ... z
//? (sel is rounded to nearest integer and limited to 0...n-1,
//? where n is number of input groups)
STATELESS_BOX_1(Sel2Box) { // TODO: io-ali
	int igl = m_arg>>4, nig = m_arg&15, itot = nig * igl, ixm = 1<<itot;
	double *pv, *po, *pi = inb[itot];
	if (!(inflg&ixm)) {
		int j = ixround(*pi, nig), j2 = j*igl, of = (inflg>>j2)&(ixm-1);
		for (int k=0; k<igl; k++) {
			if (!(of&(1<<k))) *outb[k] = *inb[j2+k]; else if ((pv=inb[j2+k])!=(po=outb[k]))
				for (int i=0;i<n;i++) po[i] = pv[i]; }
		return of;
	}
	double con[itot]; BVFOR_JM((ixm-1)&~inflg) con[j] = inb[j][0];
	for (int i=0; i<n; i++)
		for (int j = ixround(pi[i], nig), j2 = j*igl, f2 = inflg>>j2,  k = 0; k<igl; k++)
			outb[k][i] = (f2&(1<<k)) ? inb[j2+k][i] : con[j2+k];
	return (1<<igl)-1;
}


//? {{{!._sHsg}}}
//? simple histogram calculator
//? in : input
//? t : period (in seconds)
//? in every t-sized time interval, the output is the histogram
//? (the histogram is updated every t seconds, so for the first
//? t seconds the output is zero)
//? interval -1...1 is divided to (t*sample_rate) equal intvals
//? output is normalized (~= const 1.0 for equal distribution)
//? ---
//? While it is not impossible to play the output as sound, the
//? recommended usage of this box is to see stats (e.g. for noise
//? distribution) with gnuplot.
//? TODO: examples
class HistGInst : public BoxInst {
	public:
		HistGInst() : m_len(0) {}
		virtual ~HistGInst() { if (m_len) free(m_oblk); }
		virtual int calc(int inflg, double** inb, double** outb, int n);
		inline void ou(double *to, double *q, int n) { memcpy(to, q+m_pos, 8*n); m_pos+=n; }
		inline int ix(double x, int m) { 
			int j = (int)lround((x+m_add)*m_mul);
			return ((unsigned int)j < (unsigned int)m) ? j : (j<0 ? 0 : m-1); }
		void upd_blk();
	protected:
		int m_len, m_half, m_pos, m_total, *m_pcnt;
		double m_add, m_mul, *m_oblk;
};

void HistGInst::upd_blk() {
	int m = m_len; if (m_pos!=m) bug("hist/updblk: %d != %d", m_pos, m); m_pos = 0;
	double mul = (double)m / (double)(m_total+=m_len);
	for (int i=0; i<m; i++) m_oblk[i] = mul * (double)m_pcnt[i];
}

int HistGInst::calc(int inflg, double** inb, double** outb, int n) {
	int m, *q, sr = sample_rate, xf = inflg&1;  double *ob, *to = outb[0], *p = inb[0];
	if (!m_len) { if (!n) return 0;
		      m = m_len = ivlim((int)lround(inb[1][0]*(double)sr), 2, sr<<CFG_STATBUF_SIZ.i);
		      m_mul = .5 * (double)m, m_add = (double)(m-1) / (double)(m);
		      m_pos = m_total = 0;
		      char * buf = (char*)calloc((3*m+1)&~1, 4);
		      ob = m_oblk = (double*)buf; q = m_pcnt = (int*)(buf+8*m); }
	else 	   {  m = m_len; ob = m_oblk; q = m_pcnt; }
	if(xf){ while(1){ int n1 = m - m_pos;
			  if (n<n1) { for (int i=0; i<n; i++) ++q[ix(p[i],m)]; ou(to, ob, n); return 1; }
			  for (int i=0; i<n1; i++) ++q[ix(p[i],m)];     ou(to, ob, n1); upd_blk();
			  if (!(n-=n1)) return 1; else to += n1, p += n1; }}
	else {  int * q1 = q + ix(*p,m);
		while(1){ int n1 = m - m_pos;
			  if (n<n1) return *q1+=n, ou(to,ob,n), 1;
			  *q1 += n1; ou(to,ob,n1); upd_blk(); 
			  if (!(n-=n1)) return 1; else to += n1; }}}

static void sel_ini(ANode *rn) {
	ANode * sd[7]; char nm[16]; nm[1] = 0;
	for (int i=0; i<7; i++) *nm = 49+i, sd[i] = qmk_dir(rn, nm);
	memcpy(nm, "sel02", 6); qmb_arg_t qa = QMB_ARG1(Sel1Box);
	qmk_box(sd[0], nm, qa, 2, 3, 33, "sel", "R*1i-", "uu%99%", 513, "sel");
	for (int i=3; i<30; i++) nm[3] = 48+i/10, nm[4] = 48+i%10,
		qmk_box(sd[0],nm,qa, i, i+1, 33 ,"sel", "1i-", 1+256*i, "sel");
	memcpy(nm, "sel02x2", 8); qa = QMB_ARG1(Sel2Box);
	qmk_box(sd[1], nm, qa, 0x22, 5, 2, "sel2", "i:o*R*1i-", 0, 0, "a$b", "a$b", "uu%99%", 1025, "sel");
	for (int i=3; i<15; i++) nm[3] = 48+i/10, nm[4] = 48+i%10,
		qmk_box(sd[1],nm,qa, 32+i, 2*i+1, 2, "sel2", "1i-", 1+512*i, "sel");
	nm[4] = 'x'; nm[6] = 0; char inm[48], onm[16];
	for (int i=3; i<8; i++) {
		nm[3] = 50; nm[5] = 48+i;
		for (int j=0; j<i; j++) inm[3*j] = onm[2*j] = inm[24+3*j] = 97+j, 
			inm[3*j+1]=49, inm[3*j+25]=48,  inm[3*j+2]=onm[2*j+1]=inm[3*j+26]=36;
		qmk_box(sd[i-1], nm, qa, 16*i+2, 2*i+1, i, "sel2", "i-i-i-o*R*1", i, inm+24, 257*i, inm,
				512*i+1, "sel", onm, "uu%99%");
		for (int j=3, k=3*i; k<30; j++, k+=i) {
			++nm[3]; for (int v=0; v<i; v++) ++inm[3*v+1];
			qmk_box(sd[i-1], nm, qa, 16*i+j, k+1, i, "sel2", "1i-i-R1",
					256*(k-i)+i, inm, 256*k+1, "sel"); 
		}}}

void b_help_init(ANode * rn) {
	extern const char *hlpn_dir[], *hlpn_box[]; // hlpn.cc (generated)
	qmb_arg_t qa = QMB_ARG0(VersionBox);   const char *s, *s2;   ANode *dir[16]; dir[0] = rn; 
	for (int i=1; (s=hlpn_dir[i]); i++) dir[i] = qmk_dir(dir[hxd2i(*s)], s+1);
	for (int i=0; (s=hlpn_box[i]); i++) { for (s2=s+1; *s2!='.'; s2++);
					      qmk_box(dir[hxd2i(*s)], s2+1, qa, 0, 0, 1, s+1, 0); }}

#define XFT(J) "$x" #J "$fr" #J "$to" #J
void b_map_init(ANode * rn) {
	qmk_box(rn, "map01", QMB_ARG0(Map01Box), 0, 4, 33, "m01", "*i*o*", "kkk%%%", "x$y0$y1$scl", "y");
	char nm[16]; memcpy(nm, "map01*2", 8); qmb_arg_t qa = QMB_ARG1(Map01VBox);
	qmk_box(rn, nm, qa, 2, 7, 2, "m01x", "i-o.R*1", 28, "[sc]" XFT(0) XFT(1) XFT(2) XFT(3) XFT(4)
		   XFT(5) XFT(6) XFT(7) XFT(8), 0, "y", "kkk%%%");
	for (int i=3; i<=9; i++) ++nm[6], qmk_box(rn, nm, qa, i, 3*i+1, i, "m01x", "1");
}

void b_b0_init(ANode * rn) {
	ANode *mc = qmk_dir(rn, "misc"),  *wv = qmk_dir(rn, "wave"),  *im = qmk_dir(rn, "imp"), 
	      *nz = qmk_dir(rn, "noise"), *db = qmk_dir(rn, "debug"), *st = qmk_dir(rn, "stat");
	qmk_box(mc, "zero", QMB_ARG0(ZeroBox), 0, 0, 33, "_0", "*o*", "HHH%%%", "0.0");
	qmk_box(mc, "copy", QMB_ARG0(CopyBox), 0, 1, 33, "_3", "*i*o*", "kkk%%%", "x", "x");
	qmk_box(mc, "map01", QMB_ARG0(Map01Box), 0, 4, 33, "m01", "*i*o*", "kkk%%%", "x$y0$y1$scl", "y");//old
	qmk_box(mc, "rqosc", QMB_ARG0(RQOsc), 0, 5, 33, "rqo", "*i*", "k%k0%0", "f0$f1$dmp$vlim$samp");
	qmk_box(nz, "nz0", QMB_ARG0(NZ0Box), 0, 1, 33, "nz", "i*o*R*1", "ty", "out", "zzzOOO");
	qmk_box(nz, "nz0*2", QMB_ARG0(NZ0m2Box), 0, 2, 33, "nz2", "1i*", "ty1$ty2"); 
	qmk_box(nz, "rnd-pt", QMB_ARG0(RndPt), 0, 3, 33, "rpt", "i*o*R*", "fq1$fnz$vnz", "out", "z%%OOO");
	qmk_box(wv, "~pt",   QMB_ARG0(PulseTrain), 0, 2, 33, "w", "*i*R1", "z%%%%d", "fq$phs");
	qmk_box(wv, "~saw",  QMB_ARG0(SawTWave), 0, 2, 33, "w", "1+R2", "zz%");
	qmk_box(wv, "~saw01",QMB_ARG0(Saw01Wave), 0, 2, 33, "w", "2");
	qmk_box(wv, "~sin",  QMB_ARG0(SineWave), 0, 2, 33, "w", "1+", "uuu");
	qmk_box(wv, "~sin01",  QMB_ARG0(Sine01Wave), 0, 2, 33, "w", "1+", "uuu");
	qmk_box(wv, "~sq",   QMB_ARG0(SqWave), 0, 3, 33, "w1", "1+i*R2", "z%z", "fq$k$phs");
	qmk_box(wv, "~gpt",  QMB_ARG0(GPlsTrain), 0, 3, 33, "w1", "2+", "z%z");
	qmk_box(im, "^1", QMB_ARG0(Impulse1), 0, 0, 33, "_1", "*", "z%%O%%"); 
	qmk_box(im, "^triang", QMB_ARG0(TriangImp), 0, 2, 33, "_trI", "*i*", "zz%O%%", "up$dn");
	qmk_box(im, "^gauss", QMB_ARG0(GaussImp), 0, 2, 33, "_gsI", "*i*", "%zzO%%", "wid$bits");
	qmk_box(st, "histg", QMB_ARG0(HistGInst), 0, 2, 33, "sHsg", "*i*", "%%%kkk", "in$t");
	char nm[16]; memcpy(nm, "debug01", 8); qmb_arg_t qa = QMB_ARG1(DebugBox);
	qmk_box(db, nm, qa, 1, 1, 33, "dbg", "R*1", "zz%z%%");
	for (int i=2; i<31; i++) nm[5] = 48+i/10, nm[6] = 48+i%10, qmk_box(db,nm,qa,i,i,i,"dbg","1");
	sel_ini(qmk_dir(mc, "sel"));
}
