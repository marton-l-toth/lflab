#include "util.h"
#include "box0.h"
#include "glob.h"
#include "cfgtab.inc"

#define IMPBOX(NM,NA,NT) struct NM : BoxInst_B0 { \
	static scf_t sc_ini, sc_f1;   double m_arg[NA]; int m_t, m_tx[NT]; }; \
	BX_SCALC(NM::sc_ini)

#define WVBOX(NM,PJ,MUL) struct NM : BoxInst_B0 { \
		static scf_t sc_ini, sc_f1;   double m_phs; }; \
	BX_SCALC(NM::sc_ini) { SCALC_BXI(NM); bxi->m_phs=MUL*inb[PJ][0]; CALC_FW(sc_f1); } \
	BX_SCALC(NM::sc_f1)

#define PHLOOP(MUL) (x>MUL && (x-=MUL))
#define WVBOX1(NM, X, MUL) WVBOX(NM,1,MUL) { \
	SCALC_BXI(NM); double v, *p = *inb, *q = *outb, x = bxi->m_phs, mul = (MUL)*sample_length;  \
	if (inflg&1) {     	   for(int i=0;i<n;i++)    v=(X), x+=mul*p[i], q[i]=v,  PHLOOP(MUL); } \
	else	     { v=mul*p[0]; for(int i=0;i<n;i++) q[i]=(X), x+=v,			PHLOOP(MUL); } \
	bxi->m_phs = x; return 1; }

#define W2CASE(J,IX,CX,MUL) case J: IX; for(int i=0;i<n;i++) CX, PHLOOP(MUL); break
#define WVBOX2(NM, X, MUL) WVBOX(NM,2,MUL) { \
	SCALC_BXI(NM); double *p=*inb, *z=inb[1], *q=*outb, x=bxi->m_phs, mul=(MUL)*sample_length, y,v; \
	switch(inflg&3){ W2CASE(0, (v=mul*p[0], y=z[0]), (	  q[i]=(X), x+=v		), MUL); \
			 W2CASE(1, (            y=z[0]), (	     v=(X), x+=mul*p[i], q[i]=v ), MUL); \
			 W2CASE(2, (v=mul*p[0]	      ), (y=z[i], q[i]=(X), x+=v		), MUL); \
			 W2CASE(3, (void)0,		 (y=z[i],    v=(X), x+=mul*p[i], q[i]=v ), MUL); }\
	bxi->m_phs = x; return 1; }

struct ALSBox : BoxInst_B1 {
	typedef void (*c2_t) (double *q, double **inb, int n);
	static scf_t sc_fo, sc_ft;
	static int ini2(BoxInst * abxi, int inflg, double** inb, double** outb, int n, c2_t f2);
	int m_y;
};

int ALSBox::ini2(BoxInst * abxi, int inflg, double** inb, double** outb, int n, c2_t f2) {
	SCALC_BXI(ALSBox); int m = bxi->m_arg;
	if (m) { (*f2)((double*)bxi->alloc0(8*m), inb, m); CALC_FW(sc_fo); }
	if ((m=ivlim((int)lround(**inb),0,1024)) > n) {
		(*f2)((double*)bxi->alloc0(8*m), inb+1, m); bxi->m_y=m; CALC_FW(sc_ft); }
	else {  double *to = outb[0]; (*f2)(to, inb+1, m);
		int k = n-m; if(k) memset(to+m,0,8*k);   bxi->m_psc = sc_zero; return 1; }}

BX_SCALC(ALSBox::sc_fo) { SCALC_BXI(ALSBox); double *p = (double*)bxi->m_p0;
			  for (int i=0,m=bxi->m_arg; i<m; i++) outb[i][0] = p[i];   return 0; }

BX_SCALC(ALSBox::sc_ft) {
	SCALC_BXI(ALSBox);  int x = bxi->m_arg, k = bxi->m_y-x; 
	double *to = *outb, *p0 = (double*)bxi->m_p0, *p = p0+x;
	if (n<k) return memcpy(to,p,8*n), bxi->m_arg = x+n, 1;
	memcpy(to,p,8*k); if((n-=k)) memset(to+k,0,8*n);
	free(bxi->m_p0); bxi->m_p0=0; bxi->m_psc=sc_zero; return 1; }

#define LSBOX(NM) struct NM : ALSBox { static scf_t sc_ini; static void calc2(double*,double**,int); }; \
		  BX_SCALC(NM::sc_ini) { return ini2(abxi, inflg, inb, outb, n, &calc2); } \
		  void NM::calc2(double *q, double **inb, int n)

//? {{{!._0}}}
//? this box has no inputs, and outputs constant zero.
//? {{{!._2}}}
//? this box simply copies its input to its output.
//? {{{!._1}}}
//? simple impulse (first sample is 1.0, the rest are all zero)
//? HINT: since a simple impulse has a flat frequency graph, it is
//? ideal for testing the frequency response of filters.

BX_SCALC(BoxInst::sc_zero) { return **outb=0.0, 0; }
BX_SCALC(BoxInst::sc_cp0)  { return BOX_CP0; }
BX_SCALC(sc_imp1) { double *q=*outb; *q=1.0; if(--n)memset(q+1,0,8*n); abxi->m_psc=sc_zero; return 1; }

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
STATELESS_BOX(Map01Box) {
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
//? if [sc] is shorter than expected or non-list, the default
//? scale type is 0 (logarithmic)
//? ---
//? Unless you want non-constant scale types, you can replace
//? multiple "map01" boxes with one instance of this box, 
//? simplifying graph box layout.
STATELESS_BOX(Map01VBox) {
	NAN_UNPK_8(scl, inb[0], 0);
	int oflg = 0;  inflg >>= 1;  inb++;
	for (int i=0, nmap=abxi->m_arg; i<nmap; i++, inflg>>=3, inb+=3) {
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

STATELESS_BOX(VersionBox) { **outb = (double)v_major + 0.01 * (double)v_minor; return 0; }

//? {{{!._trI}}}
//? triangle impulse (up and dn are in seconds)
//? ---
//? HINT: Because up, down and amplitude are cleary visible
//? on a plot display, triangle impulses are ideal for checking
//? if scales in wrap/shadow wrap are working as expected.
IMPBOX(TriangImp, 3, 2) {
	SCALC_BXI(TriangImp); int t0 = sec2samp(inb[0][0]), t1 = sec2samp(inb[1][0]), t2 = t0+t1;
	if (!t2) CALC_FW(sc_imp1); bxi->m_t = 0;
	bxi->m_tx[0]=t0; bxi->m_tx[1]=t2; double *q = bxi->m_arg;
	q[0] = 1.0/(double)t0; q[1] = -1.0/(double)t1; q[2] = t0>0 ? 0.0 : 1.0; CALC_FW(sc_f1); }

#define TRLOOP(A,Z) for (int i=(A); i<(Z); i++) q[i]=(v+=d)
BX_SCALC(TriangImp::sc_f1) {
	SCALC_BXI(TriangImp); int k,j=0, t = bxi->m_t, t0 = bxi->m_tx[0], t1 = bxi->m_tx[1];
	double d, v = bxi->m_arg[2], *q = *outb;
	if ((k=t0-t)>0){ d=bxi->m_arg[0]; if(k>=n) { TRLOOP(0,n); goto done; } else { TRLOOP(0,k); j=k; } }
	if ((k=t1-t)>j){ d=bxi->m_arg[1]; if(k>=n) { TRLOOP(j,n); goto done; } else { TRLOOP(j,k);      } }
	bxi->m_psc = sc_zero; return (k>0) ? (memset(q+k, 0, 8*(n-k)), 1) : (*q=0.0, 0);
done:   bxi->m_t = t+n; bxi->m_arg[2] = v; return 1;
}

//? {{{!._gsI}}}
//? Gaussian impulse
//? wid - total width, out@(t=wid/2) = 1.0
//? bits - precision, out@(t=0) = out@(t=wid) = 2^(-bits-1)
IMPBOX(GaussImp, 3, 1) {
	SCALC_BXI(GaussImp); int hw = sec2samp(.5 * **inb); if (hw<1) CALC_FW(sc_imp1);
	double a = (double)hw, c1 = -M_LN2*inb[1][0], c = c1/(a*a);  bxi->m_t = 0;
	bxi->m_arg[0] = exp(c1); bxi->m_arg[1] = exp(c*(1.0-a-a)); bxi->m_arg[2]=exp(c+c);
	bxi->m_tx[0] = 2*hw+1; CALC_FW(sc_f1); }

BX_SCALC(GaussImp::sc_f1) {
	SCALC_BXI(GaussImp); double y=bxi->m_arg[0], d=bxi->m_arg[1], dd=bxi->m_arg[2], *q=*outb;
	int t = bxi->m_t, n2 = bxi->m_tx[0] - t;
	if (n2<1) return bxi->m_psc=sc_zero, *q=0.0, 0;
	if (n2<n) memset(q+n2, 0, 8*(n-n2)), n=n2, bxi->m_psc=sc_zero;
	for (int i=0; i<n; i++) q[i] = y, y *= d, d *= dd;
	bxi->m_arg[0] = y; bxi->m_arg[1] = d; bxi->m_t = t+n; return 1;
}

//? {{{!._w}}}
//? simple waves (pt - pulse train, sin, sawtooth)
//? fq: frequency, phs: start phase (0...1)
//? Note: for sound sin and saw is recommended (out: -1...1)
//? for control (with map01) use sin01 and saw01 (out: 0...1)
#define PTLOOP(X) for (int i=0; i<n; i++) y = (phs>=1.0) ? (phs-=floor(phs), 1.0) : 0.0, phs+=(X), q[i] = y;
WVBOX(PulseTrain,1,1) {
	SCALC_BXI(PulseTrain); double *q = outb[0], *fqp = inb[0], y, phs = bxi->m_phs;
	if (inflg&1) { PTLOOP(fqp[i]*sample_length); }
	else { double x = fqp[0]*sample_length; PTLOOP(x); }
	bxi->m_phs = phs; return 1;
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
WVBOX2(GPlsTrain, gpt_f(x, y), 1.0);
WVBOX2(SqWave, (x<y) ? 1.0 : 0.0, 1.0);

//? {{{!._Fib}}}
//? Fibonacci series generator
//? i0: index of first (out0) Fibonacci-number (rounded)
//? out0 ... out<n-1>: fib(round(in0)) ... fib(round(in0)+n-1)
//? ---
//? input is expected to be constant
LSBOX(FibBox) { for (int j=0,k=(int)lround(**inb); j<n; j++) q[j] = fib7s(j+k); }

//? {{{!._Pri}}}
//? Outputs prime numbers
//? p0: start (out0 will be next prime)
//? mdif: minimal diff. between primes (at least 1)
//? mmul: min. multipler between primes
//? out0 ... out<n-1>: prime numbers
//? ---
//? input is expected to be constant
//? primes found only up to 131101
LSBOX(PriBox) { int i = next_prime17((int)lround(**inb)), id = max_i(1, (int)lround(inb[1][0]));
		double x, xm = inb[2][0];
		for (int j=0; j<n; j++) q[j] = x = (double)i,
					i = next_prime17(max_i(i+id, (int)lround(x*xm))); }

//? {{{!._Sxy}}}
//? scale/list box
//? Generates a (lin/log/quad etc.) scale from x0 to x1 
//? scl: scale type (-3...3, rev.cub,rev.sq,hrm,log,lin,sq,cub)
//? ==> .!b.map.map01 -- scl explained here
//? inputs are expected to be constant
//? {{{!._Sxm}}}
//? scale/list box
//? Generates a (lin/log/quad etc.) scale from x0 to x0*xm
//? scl: scale type (-3...3, rev.cub,rev.sq,hrm,log,lin,sq,cub)
//? ==> .!b.map.map01 -- scl explained here
//? inputs are expected to be constant
LSBOX(SxyBox) { Scale01::vec(q, inb[0][0], inb[1][0],		n, ivlim((int)lround(inb[2][0]),-3,3)); }
LSBOX(SxmBox) { Scale01::vec(q, inb[0][0], inb[0][0]*inb[1][0], n, ivlim((int)lround(inb[2][0]),-3,3)); }

STATELESS_BOX(DebugBox) {
	if0(glob_flg & GLF_INI0) return 0;
	int arg = abxi->m_arg; log_n("debug_box[%d]: n=%d", arg, n);
	for (int i=0,m=1; i<arg; i++,m+=m) {
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
struct RQOsc : BoxInst_B0 { static scf_t sc_ini, sc_1;   double m_y, m_v; };

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

BX_SCALC(RQOsc::sc_ini) { SCALC_BXI(RQOsc); bxi->m_y = bxi->m_v = 0.0; CALC_FW(sc_1); }
BX_SCALC(RQOsc::sc_1) {
	SCALC_BXI(RQOsc); PREP_INPUT(f0, 0); PREP_INPUT(f1, 1); PREP_INPUT(vl, 3);
	int ty = (inflg&4) + ((int)lround(inb[4][0]) & 3);
	double ya, d0, d1, f, f0, f1, vl, el, *pel = inb[2], y = bxi->m_y, v = bxi->m_v, *q = outb[0];
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
	bxi->m_y = CUT300(y); bxi->m_v = CUT300(v); return 1; 
}

inline int ixround(double x, int n) { int r = (int)lround(x);
	return ((unsigned int)r<(unsigned int)n) ? r : (r<0?0:n-1); }

//? {{{!._sel}}}
//? this box selects the sel-th of the inputs in0...in<n-1>
//? (sel is rounded to nearest integer and limited to 0...n-1)
STATELESS_BOX(Sel1Box) {
	int arg = abxi->m_arg, ixm = 1<<arg; 
	double *pv, *q = *outb, *pi = inb[arg];
	if (!(inflg&ixm)) {
		int j = ixround(*pi, arg), r = (inflg>>j)&1;
		if (!r) *q = *inb[j]; else if ((pv=inb[j])!=q) for (int i=0;i<n;i++) q[i]=pv[i];
		return r;
	}
	double con[arg]; BVFOR_JM((ixm-1)&~inflg) con[j] = inb[j][0];
	for (int j,i=0; i<n; i++) j    = ixround(pi[i], arg), 
				  q[i] = (inflg&(1<<j)) ? inb[j][i] : con[j];
	return 1;
}

//? {{{!._sel2}}}
//? this box selects the sel-th group of inputs, copying
//? a<sel>, b<sel>, ... z<sel> to a, b, ... z
//? (sel is rounded to nearest integer and limited to 0...n-1,
//? where n is number of input groups)
STATELESS_BOX(Sel2Box) { // TODO: io-ali
	int arg = abxi->m_arg, igl = arg>>4, nig = arg&15, itot = nig * igl, ixm = 1<<itot;
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
struct HistGInst : BoxInst_B1 {
	static scf_t sc_ini, sc_f1;
	inline void ou(double *to, double *q, int n) { memcpy(to, q+m_pos, 8*n); m_pos+=n; }
	inline int ix(double x, int m) { 
		int j = (int)lround((x+m_add)*m_mul);
		return ((unsigned int)j < (unsigned int)m) ? j : (j<0 ? 0 : m-1); }
	void upd_blk();
	int m_len, m_half, m_pos, m_total, *m_pcnt;
	double m_add, m_mul;
};

void HistGInst::upd_blk() {
	int m = m_len; if (m_pos!=m) bug("hist/updblk: %d != %d", m_pos, m); m_pos = 0;
	double mul = (double)m / (double)(m_total+=m_len), *q=(double*)m_p0;
	for (int i=0; i<m; i++) q[i] = mul * (double)m_pcnt[i];
}

BX_SCALC(HistGInst::sc_ini) {
	SCALC_BXI(HistGInst); int sr = sample_rate; 
	int m = bxi->m_len = ivlim((int)lround(inb[1][0]*(double)sr), 2, sr<<CFG_STATBUF_SIZ.i);
	bxi->m_mul = .5*(double)m, bxi->m_add = (double)(m-1) / (double)(m);
	int as = 4*((3*m+1)&~1); char *buf = (char*)bxi->alloc0(as); memset(buf,0,as);
	bxi->m_pos = bxi->m_total = 0; bxi->m_pcnt = (int*)(buf+8*m);
	CALC_FW(sc_f1);
}

BX_SCALC(HistGInst::sc_f1) {
	SCALC_BXI(HistGInst); int m=bxi->m_len, *q=bxi->m_pcnt, xf=inflg&1;
	double *to=outb[0], *p=inb[0], *ob = (double*)bxi->m_p0;
	if(xf){ while(1){ int n1 = m - bxi->m_pos;
			  if (n<n1) { for(int i=0;i<n;i++) ++q[bxi->ix(p[i],m)]; bxi->ou(to,ob,n); return 1; }
			  for (int i=0; i<n1; i++) ++q[bxi->ix(p[i],m)]; bxi->ou(to, ob, n1); bxi->upd_blk();
			  if (!(n-=n1)) return 1; else to += n1, p += n1; }}
	else {  int * q1 = q + bxi->ix(*p,m);
		while(1){ int n1 = m - bxi->m_pos;
			  if (n<n1) return *q1+=n, bxi->ou(to,ob,n), 1;
			  *q1 += n1; bxi->ou(to,ob,n1); bxi->upd_blk(); 
			  if (!(n-=n1)) return 1; else to += n1; }}}

static void sel_ini(ANode *rn) {
	ANode * sd[7]; char nm[16]; nm[1] = 0;
	for (int i=0; i<7; i++) *nm = 49+i, sd[i] = qmk_dir(rn, nm);
	memcpy(nm, "sel02", 6); qmb_arg_t qa = QMB_A_SL(Sel1Box);
	qmk_box(sd[0], nm, qa, 2, 3, 33, "sel", "R*1i-", "uu%99%", 513, "sel");
	for (int i=3; i<30; i++) nm[3] = 48+i/10, nm[4] = 48+i%10,
		qmk_box(sd[0],nm,qa, i, i+1, 33 ,"sel", "1i-", 1+256*i, "sel");
	memcpy(nm, "sel02x2", 8); qa = QMB_A_SL(Sel2Box);
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
	qmb_arg_t qa = QMB_A_SL(VersionBox);   const char *s, *s2;   ANode *dir[16]; dir[0] = rn; 
	for (int i=1; (s=hlpn_dir[i]); i++) dir[i] = qmk_dir(dir[hxd2i(*s)], s+1);
	for (int i=0; (s=hlpn_box[i]); i++) { for (s2=s+1; *s2!='.'; s2++);
					      qmk_box(dir[hxd2i(*s)], s2+1, qa, 0, 0, 1, s+1, 0); }}

#define XFT(J) "$x" #J "$fr" #J "$to" #J
void b_map_init(ANode * rn) {
	qmk_box(rn, "map01", QMB_A_SL(Map01Box), 0, 4, 33, "m01", "*i*o*", "kkk%%%", "x$y0$y1$scl", "y");
	char nm[16]; memcpy(nm, "map01*2", 8); qmb_arg_t qa = QMB_A_SL(Map01VBox);
	qmk_box(rn, nm, qa, 2, 7, 2, "m01x", "i-o.R*1", 28, "[sc]" XFT(0) XFT(1) XFT(2) XFT(3) XFT(4)
		   XFT(5) XFT(6) XFT(7) XFT(8), 0, "y", "kkk%%%");
	for (int i=3; i<=9; i++) ++nm[6], qmk_box(rn, nm, qa, i, 3*i+1, i, "m01x", "1");
}

#define LS_DIR(R,NM,NI,IN) ls_dir(R, QMB_A_BX(NM##Box), #NM, NI, IN)
static ANode * ls_dir(ANode *rn, qmb_arg_t qa, const char *nm0, int ni, const char *inm) {
	rn = qmk_dir(rn, nm0);  int l = strlen(nm0); char nm[24]; memcpy(nm, nm0, l);
	memcpy(nm+l, "01", 3); qmk_box(rn, nm, qa, 1, ni,   33, nm0, "i*R1", inm+4);
 	for (int i=2; i<26; i++) nm[l] = 48+i/10, nm[l+1] = 48+i%10, qmk_box(rn,nm,qa,i,ni,i+32,nm0,"1"); 
	return memcpy(nm+l, "*",  2), qmk_box(rn, nm, qa, 0, ni+1, 33, nm0, "i*",   inm); }

static void ls_ini(ANode *rn) {
	LS_DIR(rn, Fib, 1, "siz$i0");		reg_bn(LS_DIR(rn, Sxy, 3, "siz$x0$x1$scl"), 8);
	LS_DIR(rn, Pri, 3, "siz$p0$mdif$mmul"); reg_bn(LS_DIR(rn, Sxm, 3, "siz$x0$xm$scl"), 9); }

void b_b0_init(ANode * rn) {
	ANode *mc = qmk_dir(rn, "misc"),  *wv = qmk_dir(rn, "wave"),  *im = qmk_dir(rn, "imp"), 
	      *db = qmk_dir(rn, "debug"), *st = qmk_dir(rn, "stat");
	qmk_box(mc, "zero", BoxInst::sc_zero, 0, 0, 33, "_0", "*o*", "HHH%%%", "0.0");
	qmk_box(mc, "copy", BoxInst::sc_cp0,  0, 1, 33, "_3", "*i*o*", "kkk%%%", "x", "x");
	qmk_box(mc, "map01", QMB_A_SL(Map01Box), 0, 4, 33, "m01", "*i*o*", "kkk%%%", "x$y0$y1$scl", "y");//old
	qmk_box(mc, "rqosc", QMB_A_BX(RQOsc), 0, 5, 33, "rqo", "*i*", "k%k0%0", "f0$f1$dmp$vlim$samp");
	qmk_box(wv, "~pt",   QMB_A_BX(PulseTrain), 0, 2, 33, "w", "*i*R1", "z%%%%d", "fq$phs");
	qmk_box(wv, "~saw",  QMB_A_BX(SawTWave), 0, 2, 33, "w", "1+R2", "zz%");
	qmk_box(wv, "~saw01",QMB_A_BX(Saw01Wave), 0, 2, 33, "w", "2");
	qmk_box(wv, "~sin",  QMB_A_BX(SineWave), 0, 2, 33, "w", "1+", "uuu");
	qmk_box(wv, "~sin01",  QMB_A_BX(Sine01Wave), 0, 2, 33, "w", "1+", "uuu");
	qmk_box(wv, "~sq",   QMB_A_BX(SqWave), 0, 3, 33, "w1", "1+i*R2", "z%z", "fq$k$phs");
	qmk_box(wv, "~gpt",  QMB_A_BX(GPlsTrain), 0, 3, 33, "w1", "2+", "z%z");
	qmk_box(im, "^1", QMB_A_BX(Impulse1), 0, 0, 33, "_1", "*", "z%%O%%"); 
	qmk_box(im, "^triang", QMB_A_BX(TriangImp), 0, 2, 33, "_trI", "*i*", "zz%O%%", "up$dn");
	qmk_box(im, "^gauss", QMB_A_BX(GaussImp), 0, 2, 33, "_gsI", "*i*", "%zzO%%", "wid$bits");
	qmk_box(st, "histg", QMB_A_BX(HistGInst), 0, 2, 33, "sHsg", "*i*", "%%%kkk", "in$t");
	char nm[16]; memcpy(nm, "debug01", 8); qmb_arg_t qa = QMB_A_SL(DebugBox);
	qmk_box(db, nm, qa, 1, 1, 33, "dbg", "R*1", "zz%z%%");
	for (int i=2; i<31; i++) nm[5] = 48+i/10, nm[6] = 48+i%10, qmk_box(db,nm,qa,i,i,i,"dbg","1");
	ls_ini(qmk_dir(rn, "ls")); sel_ini(qmk_dir(mc, "sel"));
}
