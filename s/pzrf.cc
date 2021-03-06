#include <ctype.h>
#include <math.h>

#include "util.h"
#include "pzrf.h"
#include "util2.h"
#include "glob.h"
#include "guistub.h"
#include "cmd.h"
#include "cfgtab.inc"

double RECF_PZItem::eval1(double re, double im) {
	double v = 1.0;
	if (zr==zr) v  = ( (re-zr)*(re-zr) + (im-zi)*(im-zi) ) 
		       * ( (re-zr)*(re-zr) + (im+zi)*(im+zi) );
	if (pr==pr) v /= ( (re-pr)*(re-pr) + (im-pi)*(im-pi) )
		       * ( (re-pr)*(re-pr) + (im+pi)*(im+pi) );
	return sqrt(v);
}

int RECF_PZItem::eql(double wfq, double wid, double amp) {
	int r = approx_cmp(amp, 1.0); if (!r) return 0;
	int sg = (amp<0.0) && (amp=-amp, 1);
	double pa, za;    if (r<0) pa = wid, za = wid*amp; else za = wid, pa = wid/amp;
	pr = -wfq*sin(pa); pi = wfq*cos(pa); c = 1.0;
	zr =  wfq*sin(za); zi = wfq*cos(za); if (!sg) zr = -zr;   return 1;
}

static const double tanhalf = tan(0.5);
static const double tanhalf_r = 1.0 / tanhalf;
static const double rft_T = 2*tanhalf;
static const double rft_T2 = rft_T * rft_T;

static RECF_PZItem * lpz_p = 0;
static int lpz_n = 0;

double fq_warp(double fq) { return tan(fq * sample_length * M_PI) * tanhalf_r; }
double fq_warp2(double fq2) { return tan(fq2) * tanhalf_r; }
void fq_warp_ir2(double *pr, double *pi, double *p2) {
	if (*p2<1e-10) return;
	double fq = sqrt(*p2), mul = fq_warp(fq)/fq;
	*pr *= mul; *pi *= mul; *p2 *= mul*mul; }

static void rfpz_store(RECF_PZItem * pz, int n) {
	if (lpz_p) delete[](lpz_p); lpz_p = pz; lpz_n = n; }

void rfpz_transform(RECF_ABXY * to, RECF_PZItem * pz, int n, int flg) {
	int stf = flg&(n>0); if (stf) rfpz_store(pz,n);
	for (int i=0; i<n; i++) {
		memset(to+i, 0, 32);
		double pr=pz[i].pr, pi=pz[i].pi;
		double zr=pz[i].zr, zi=pz[i].zi;
		double c = pz[i].c;
		if (pr == pr) {
			double p2 = pr*pr + pi*pi;
			if (flg&2) fq_warp_ir2(&pr, &pi, &p2);
			double mul = 1.0 / (4.0 - 4.0*rft_T*pr + rft_T2*p2);
			to[i].b1 = mul * (8.0 - 2*rft_T2*p2);
			to[i].b2 = mul * (-4.0 - 4.0*rft_T*pr - p2*rft_T2);
			if (zr == zr) {
				double z2 = zr*zr + zi*zi;
				if (flg&2) fq_warp_ir2(&zr, &zi, &z2);
				to[i].a0 = c*mul * (4.0 - 4.0*rft_T*zr + z2*rft_T2);
				to[i].a1 = c*mul * (-8.0 + 2*rft_T2*z2);
				to[i].a2 = c*mul * (4.0 + 4.0*rft_T*zr + z2*rft_T2);
			} else {
				double v = c*mul * rft_T2;
				to[i].a0 = to[i].a2 = v;
				to[i].a1 = v+v;
			}
		} else if (zr == zr) {
			double z2 = zr*zr + zi*zi;
			if (flg&2) fq_warp_ir2(&zr, &zi, &z2);
			double mul = pz[i].c / rft_T2;
			to[i].a0 = c*mul * (4.0 - 4.0*rft_T*zr + z2*rft_T2);
			to[i].a1 = c*mul * (-8.0 + 2*rft_T2*z2);
			to[i].a2 = c*mul * (4.0 + 4.0*rft_T*zr + z2*rft_T2);
			to[i].b1 = -2.0; 
			to[i].b2 = -1.0;
		} else {
			to[i].a0 = c;
			to[i].a1 = to[i].a2 = 0.0;
			to[i].b1 = to[i].b2 = 0.0;
		}
	}
	if (!stf) delete[](pz);
}

static double get_recrad() {
	double r = 0.0, x;
	for (int i=0; i<lpz_n; i++) {
		RECF_PZItem * p = lpz_p + i;
		if (p->pr==p->pr) { if ((x=fabs(p->pr))>r) r=x; if ((x=fabs(p->pi))>r) r=x; }
		if (p->zr==p->zr) { if ((x=fabs(p->zr))>r) r=x; if ((x=fabs(p->zi))>r) r=x; }
	}
	return r;
}

static void paraconj(RECF_PZItem*p,int n) {for(RECF_PZItem*q=p+n;p<q;p++) p->zr=-p->pr, p->zi=p->pi, p->c=1.0;}

void pzrf_show_last() {
	double scl = 372.48 / get_recrad();
	char buf[8*lpz_n];
	int xi ,yi, len = 0;
	for (int i=0; i<lpz_n; i++) {
		RECF_PZItem * p = lpz_p + i;
		if (p->pr==p->pr) xi = 383+(int)lround(scl*p->pr), yi = 383+(int)lround(scl*p->pi),
				  buf[len++] = ((xi>>6)&31) + 48, buf[len++] = (xi&63) + 48,
				  buf[len++] = ((yi>>6)&31) + 48, buf[len++] = (yi&63) + 48;
		if (p->zr==p->zr) xi = 383+(int)lround(scl*p->zr), yi = 383+(int)lround(scl*p->zi),
				  buf[len++] = ((xi>>6)&31) + 80, buf[len++] = (xi&63) + 48,
				  buf[len++] = ((yi>>6)&31) + 48, buf[len++] = (yi&63) + 48;
	}
	gui2.cre(0x77, 'P'); gui2.wupd('P'); gui2.sn(buf, len);
}

void PZFInst::set_n_d(int n, double d) { char *p = (char*)(m_p0 = malloc(n*sizeof(RECF_ABXY)+15));
					 m_n_bq = n; m_ab = (RECF_ABXY*)(p + (-((long)p) & 15)); m_dcy=d; }
void PZFInst::dcy() {
	if (!m_n_bq || m_dcy<1e-14 || m_dcy>1e20) return;
	double d = exp(-m_dcy * sample_length), dd = d*d;
	for (int i=0; i<m_n_bq; i++) m_ab[i].b1*=d, m_ab[i].b2*=dd, m_ab[i].a1*=d, m_ab[i].a2*=dd;
}


void PZFInst::rf_debug(int n) {
	m_t += n; if (m_t < sample_rate) return; else m_t -= sample_rate;
	log_n("rf_debug:"); int m = m_n_bq;
	for (int i=0; i<m; i++) log_n(" (%g %g %g %g)", m_ab[i].x21[1], m_ab[i].x21[0], 
							m_ab[i].y21[1], m_ab[i].y21[0]);
	log("");
}

#define PZF_CI(X){double x=in[0]; if(fabs(x)<1e-270){in=zeroblkD+(X);} else{in=p;for(int i=0;i<n;i++) p[i]=x;}}
#define PZF0(ZX) int m=m_n_bq; if (m<1) {do{if(m){ int r=ini(inb); if(r<0)return r; if ((m=m_n_bq)) break; \
					  	   return (ZX); }}while(0); }
#define PZF_II0 { RECF_ABXY *s=bxi->m_ab, *sl=s+m; do
#define PZF_II1 while(++s < sl); }
#define PZF_COL(Q,P) { double z = *(P); PZF_II0 { \
	double x1=s->x21[1], y1=s->y21[1], y = z*s->a0 + x1*s->a1 + *s->x21*s->a2 + y1*s->b1 + *s->y21*s->b2; \
	*s->x21=x1, s->x21[1]=z, *s->y21=y1, z=s->y21[1]=y; } PZF_II1 *(Q) = z; }

#define PZF_1(I,X) \
        RECF_ABXY * vv = bxi->m_ab + (I); \
        double x1 = vv->x21[1], x2 = vv->x21[0],                y1 = vv->y21[1], y2 = vv->y21[0], \
               a0 = vv->a0,     a1 = vv->a1,     a2 = vv->a2,   b1 = vv->b1,     b2 = vv->b2;     \
        for (int j=0; j<n; j++) {                                  \
                double y = a0*in[j] + a1*x1 + a2*x2 + b1*y1 + b2*y2;\
                x2 = x1; x1 = in[j]; y2 = y1; y1 = y; (X);        } \
        vv->x21[0]=CUT300(x2); vv->x21[1]=CUT300(x1); vv->y21[0]=CUT300(y2); vv->y21[1]=CUT300(y1);

#ifdef __SSE2__
#include <emmintrin.h>
#define PZF_LDC(Z,J) Z##J##J = _mm_load_pd1(&s->Z##J)
#define ADD2 _mm_add_pd
#define MUL2 _mm_mul_pd

int PZFInst::ini2(BoxInst * abxi, int inflg, double** inb, double** outb, int n, int r) {
	static sc_t calctab[6]={sc_cp0, sc_zero, sc_one, sc_one, sc_seq, sc_par};
	if0(r<0) return r; SCALC_BXI(PZFInst);
	bxi->m_t=0; bxi->dcy(); int m = bxi->m_n_bq; CALC_FW(calctab[ 2*(!!m+(m>1)) + (r&1) ]); }

BX_SCALC(PZFInst::sc_one) { SCALC_BXI(PZFInst); double *in=inb[0], *p=outb[0];
			    if (!(inflg&1)) PZF_CI(0)  PZF_1(0, (p[j]=y)); return 1; }

BX_SCALC(PZFInst::sc_seq) {
	SCALC_BXI(PZFInst); int m = bxi->m_n_bq;
	double *in0, *in = inb[0], *p = outb[0];
	switch ( (((long)p>>1)&4) + 2*(inflg&1) + ((n>10)&CFG_DEVEL.i)) {
		case 0: case 4: PZF_CI(0) 
		case 2: case 6: for (int i=0; i<m; i++) { PZF_1(i, (p[j]=y)); in = p; } return 1;
		case 1: PZF_CI(0); goto even; 
		case 5: PZF_CI(1); goto odd;
		case 3: if (  (long)in&8 ) memcpy(p,in,8*n), in=p;goto even;
		case 7: if (!((long)in&8)) memcpy(p,in,8*n), in=p; goto odd;
		default:break;
	}
odd:	PZF_COL(p,in); --n; ++p; ++in;
even:   --n; in0 = in; PZF_II0 {
		__m128d PZF_LDC(a,0), PZF_LDC(a,1), PZF_LDC(a,2), PZF_LDC(b,1), PZF_LDC(b,2),
			x21 = _mm_load_pd(s->x21), y21 = _mm_load_pd(s->y21);
		int j=0; do {
			__m128d p01=_mm_load_pd (in+j), 
				t21=ADD2(ADD2(MUL2(a00,p01), MUL2(a22,x21)),
				    	 ADD2(MUL2(b22,y21), MUL2(a11,_mm_shuffle_pd(x21,p01,1))));
			t21 = ADD2(t21, MUL2(b11, _mm_shuffle_pd(y21, _mm_setzero_pd(), 1))); x21 = p01; \
			_mm_store_pd(p+j, y21=ADD2(t21, MUL2(b11, _mm_shuffle_pd(_mm_setzero_pd(), t21, 0))));
		} while ((j+=2)<n);
		_mm_store_pd(s->x21, x21); _mm_store_pd(s->y21, y21); in = p;
	} PZF_II1
	if (!(n&1)) PZF_COL(p+n, in0+n);
	return 1;
}
#else // __SSE2__
BX_SCALC(PZFInst::sc_seq) {
	SCALC_BXI(PZFInst); int m = bxi->m_n_bq;
	double *in = inb[0], *p = outb[0]; if (!(inflg&1)) PZF_CI(0)
	for (int i=0; i<m; i++) { PZF_1(i, (p[j] = y)); in = p; }
	IFDBGX(RFINST) bxi->rf_debug(n);    		return 1;
}
#endif // __SSE2__

BX_SCALC(PZFInst::sc_par) {
	SCALC_BXI(PZFInst); int m = bxi->m_n_bq; inflg &= 1;
	double *in = inb[0], *p = outb[0];
	int cif = (!inflg && fabs(*in)>1e-270), xtf = (in==p);
	double tbuf[(cif|xtf)?n:1];
	if (cif) { double x = *in; in = tbuf; for (int j=0; j<n; j++) tbuf[j] = x; goto notmp; }
	if (!inflg) { in = zeroblkD; goto notmp; }
	if (!xtf) goto notmp;
	{PZF_1( 0 , tbuf[j]=y);}; for (int i=1; i<m-1; i++) {PZF_1(i, tbuf[j]+=y);} 
	{PZF_1(m-1, p[j]=tbuf[j]+y);} return 1; goto done;
notmp:  {PZF_1(0, p[j]=y);} for (int i=1; i<m; i++) { PZF_1(i, p[j]+=y); }
done:   IFDBGX(RFINST) bxi->rf_debug(n);    		return 1;
}

static void half_zero(RECF_ABXY * to, double fq0, double bw1, double c = 1.0) {
	memset(to, 0, 32);
	double fq = fq_warp(fq0), bw = bw1 * fq, cf = cos(fq),
	       r  = 1.0 - 3.0*bw, rr = r*r, rcf = r*cf, k = (1.0 - 2.0*rcf + rr) / (2.0 - 2.0*cf);
	to->a0 = c*(1.0-k);	to->a1 = (c+c)*(k-r)*cf;	to->a2 = c*(rr-k);
				to->b1 = rcf+rcf;		to->b2 = -rr;
}

QUICK_PZFILT(SBandF) { set_n_d(1); half_zero(m_ab, inb[0][0], inb[1][0]); return 0; }
//? {{{!._bF}}}
//? Simple band-pass/high-pass filter (single pole, half-zero)
//? in - filter input
//? fq - frequency
//? Q - sharpness
//? Q=0.01: narrow band, with some high freq let through
//? Q=0.1: wider band, lot of high freq let through
//? Q=0.3: high-pass with a little bump at "fq"
//? Q>.5: high-pass

//? {{{!._eqF}}}
//? Equalizer filter
//? in - filter input
//? fq1 - base frequency
//? f*<i> - frequency (multplied by fq1)
//? wid<i> - width (0...1)
//? amp<i> - amplification for fq1 * f*<i>
//? amp may be negative, which makes phase different
//? amp=-1 -> (almost) no frq. amplitude change, only phase
QUICK_PZFILT(EqlzF) {
	double fq1 = inb[0][0]; 
	int n = 0, n0 = m_arg;
	RECF_PZItem * pzn = new RECF_PZItem[n0];
	for (int i=0; i<n0; i++) n += pzn[n].eql(fq_warp(fq1*inb[3*i+1][0]), inb[3*i+2][0], inb[3*i+3][0]);
	set_n_d(n); rfpz_transform(m_ab, pzn, n); return 0;
}

//? {{{!._eqLF}}}
//? Equalizer filter (list for harmonic frqs.)
//? in - filter input
//? fq1 - base frequency
//? step - frq. step (*fq1)
//? wid - width (0..1) for base freq.
//? 2div - amp logscale divider
//? [a] - amp list, j-th amp is 2^(a<j>/2div)
QUICK_PZFILT(EqlzLsF) {
	double fq = inb[0][0], step = fq*inb[1][0], wid1 = inb[2][0], 
	       amp1 = exp(M_LN2 / inb[3][0]), wid2 = wid1 * fq;
	NAN_UNPK_32(amp, inb[4], 0);
	int i, j, k, np = amp_n;   set_n_d(np);
	RECF_PZItem * pzn = new RECF_PZItem[np];
	for (i=j=0; i<np; i++, fq+=step) if ((k=amp_x[i])) pzn[j].eql(fq_warp(fq),wid2/fq,ipows(amp1,k)), j++;
	rfpz_transform(m_ab, pzn, m_n_bq = j); return 0;
}

//? {{{!._hppF}}}
//? parallel (added) sequence of single-pole-pair filters
//? in - filter input
//? fq - base frequency
//? Q - angle (0..1) - lower Q->sharper freq.resp., longer snd.
//? np - # of poles
//? z/2p - number of zeroes for each pole pair (0, 1, or 2)
//? stp - freq. step (multiplied by fq)
//? cfac - todo....
//? dcy - decay
QUICK_PZFILT(HPPF) { //"in$fq$Q$np$z/2p$stp$cfac$dcy"
	double fq1 = inb[0][0], bw = inb[1][0], stp1 = inb[4][0], cfac = inb[5][0],
	       fq = fq1, stp = fq1*stp1, c = 1.0;
	int np = (int)lround(inb[2][0]), zp2p = (int)lround(inb[3][0]);
	set_n_d(np, inb[6][0]);
	if (zp2p == 1) { for (int i=0; i<np; i++, fq+=stp, c*=cfac) half_zero(m_ab+i, fq, bw, c); return 1; }
	RECF_PZItem * pz = new RECF_PZItem[np];
	double fqw, bw2 = .5*M_PI*bw, cr = -sin(bw2), ci = cos(bw2), zri = zp2p ? 0.0 : NAN;
	for (int i=0; i<np; i++, fq+=stp, c*=cfac) fqw = fq_warp(fq), pz[i].c = c,
		pz[i].pr = cr*fqw, pz[i].pi = ci*fqw, pz[i].zr = pz[i].zi = zri;
	rfpz_transform(m_ab, pz, np); return 1;
}

//? {{{!._hpsF}}}
//? parallel (added) sequence of single-pole-pair filters
//? in - filter input
//? fq - base frequency
//? Q0, Q1, Qs - angle (0..1) - from, to scale (see map01)
//? c0, c1, cs - const factor - from, to scale (see map01)
//? np - # of poles
//? z/2p - number of zeroes for each pole pair (0, 1, or 2)
//? stp - freq. step (multiplied by fq)
//? dcy - decay
//? ==> .!b.misc.map01
QUICK_PZFILT(HPPscF) { //"in$fq$Q$np$z/2p$stp$Q0$Q1$Qs$c0$c1$cs$dcy"
	int np = (int)lround(inb[1][0]), zp2p = (int)lround(inb[2][0]);
	double fq1 = inb[0][0], stp1 = inb[3][0],
	       fq = fq1, stp = fq1*stp1, t = 0.0, tstp = 1.0/(double)(np-1);
	Scale01 qsc, csc;
	qsc.set_all(inb[4][0], inb[5][0], (int)lround(inb[6][0]));
	csc.set_all(inb[7][0], inb[8][0], (int)lround(inb[9][0]));
	set_n_d(np, inb[10][0]); 
	if (zp2p == 1) { for (int i=0; i<np; i++, fq+=stp, t+=tstp) 
		half_zero(m_ab+i, fq, qsc.f(t), csc.f(t)); return 1; }
	RECF_PZItem * pz = new RECF_PZItem[np];
	double fqw, zri = zp2p ? 0.0 : NAN, ang;
	for (int i=0; i<np; i++, fq+=stp, t+=tstp) fqw = fq_warp(fq), pz[i].c = csc.f(t),
		ang = qsc.f(t) * .5 * M_PI, pz[i].pr = -sin(ang)*fqw, pz[i].pi = cos(ang)*fqw,
		pz[i].zr = pz[i].zi = zri;
	rfpz_transform(m_ab, pz, np); return 1;
}

#define CHECK_NPZ(nm) if (npz[1]>npz[0]) return log("mkf/%s: cowardly refusing"\
	" to create filter with #z(%d) > #p(%d)", #nm, npz[1], npz[0]), RTE_FZZZ; \
	set_n_d(*npz); RECF_PZItem *p, *pzi = new RECF_PZItem[*npz]

#define CHECK_NPZ_2(nm) \
	if (np!=npz[0]) return log("mkf/%s: BUG: wrong num of poles(exp:%d got:%d)", #nm, npz[0],np),RTE_BUG;\
	if (nz!=npz[1]) return log("mkf/%s: BUG: wrong num of zeros(exp:%d got:%d)", #nm, npz[1],nz),RTE_BUG;\
	for (int i=nz; i<np; i++) pzi[i].zr = pzi[i].zi = NAN;

#define FQ_NORM(x) double nfq = fq_warp2(fq1*inb[x][0]); \
	for (int i=0; i<np; i++) pzi[i].c = 1.0 / pzi[i].eval1(0.0, nfq);

#define PZ_END(nm, x) CHECK_NPZ_2(nm); FQ_NORM(x); rfpz_transform(m_ab, pzi, m_n_bq = np); return 0;

//? {{{!._bqF}}}
//? generic biquad sequence filter
//? in - filter input
//? fq1 - base frq
//? [#] - count for each (r<j>, a<j>): >0: pole <0: zero
//? nfq - freq. resp. will be normalized for fq1*nfq
//? f<j> - freq (multiplied by fq1)
//? q<j> - angle (0..1 for pole, -1..1 for zero, mult.by -pi/2)
//? CPU usage is approx. linear with total number of poles
//? HINT: use main menu/filter disp. to visualize poles/zeroes
QUICK_PZFILT(BiQuad1F) {
	log("hello/bq1");
	double fq1 = **inb * M_PI * sample_length;
	NAN_UNPK_8(dsc, inb[1], 0);
	int n = min_i(dsc_n, m_arg); if (!n) return 0;
	int np=0, nz=0, npz[2]; npz[0] = npz[1] = 0;
	for (int k,i=0; i<n; i++) k = dsc_x[i], dsc_x[i] = k = (k>=0) ? (k?k:-127) : -128-k,
				  npz[(k>>7)&1] += k&127;
	CHECK_NPZ(bq);
	for (int i=0; i<n; i++) {
		int k = dsc_x[i] & 127, zf = dsc_x[i] <= 0;
		if (k<0) k = -k; else if (!k) k = 1;
		double v = fq_warp2(inb[2*i+3][0] * fq1),
		       ang = .5 * M_PI * inb[2*i+4][0],
		       re = -v*sin(ang), im = v*cos(ang);
		if (zf) for(int j=0;j<k;j++) p=pzi+(nz++), p->zr=re, p->zi=im;
		else    for(int j=0;j<k;j++) p=pzi+(np++), p->pr=re, p->pi=im;
	}
	PZ_END(bq, 2);
}

//? {{{!._bqsF}}}
//? generic biquad sequence filter with scales
//? in - filter input
//? fq1 - base frq
//? [sc] - [rsc asc] for each (n, ra, rz, aa, az)   (*)
//? [md] - [p/z ri ai] for each (n, ra, rz, aa, az) (**)
//? nfq - freq. resp. will be normalized for fq1*nfq
//? n(j) - number of poles/zeroes
//? f(j)<, >f(j) - interval for rel. frq (mult. by fq1)
//? q(j)<, >q(j) - interval for angles (mult. by -pi/2)
//? (*) scale (-3...3) for rel.frq and angle
//? (**) p/z: 0:zero 1:pole 
//? ri,ai: interval type (0: [] 1: ]] 2: [[ 3:][ )
//? CPU usage is approx. linear with total number of poles
//? HINT: use main menu/filter disp. to visualize poles/zeroes
QUICK_PZFILT(BiQuadSeqF) {
	log("hello/bqs");
	static const double div1[4] = { .5, 1.0, 0.0, .5 };
	double fq1 = **inb * M_PI * sample_length;
	NAN_UNPK_8(scl, inb[1], 1);
	int n = min_i(scl_n/2, m_arg); if (!n) return 0;
	NAN_UNPK_8(mode, inb[2], 0);
	int np=0, nz=0, n2[n], md[2], npz[2]; npz[0] = npz[1] = 0;
	for (int i=0; i<n; i++) npz[~mode_x[3*i]&1] += (n2[i] = (int)lround(inb[5*i+4][0]));
	CHECK_NPZ(bqs);
	Scale01 pz01[2];
	double t0[2], td[2];
	log("bqs: n=%d", n);
	for (int i=0; i<n; i++) {
		int m = n2[i]; if (!m) continue;
		double ** pp = inb + 5*i; char *ps = scl_x + 2*i, *pmd = mode_x + 3*i;
		for (int j=0; j<2; j++) {
			pz01[j].set_all(pp[2*j+5][0], pp[2*j+6][0], ps[j]);
			md[j] = pmd[j+1];
			if (m==1) { t0[j] = div1[md[j]&3]; td[j] = 0.0; continue; }
			double t = 1.0 / (double) ((md[j]&1) + ((md[j]>>1)&1) + m + m - 2);
			t0[j] = (md[j]&1) ? t : 0; td[j] = t + t;
		}
		int zf = !*pmd; 
		log("bqs: adding %d %ses", m, zf?"zero":"pol");
		for (int j=0; j<m; j++,t0[0]+=td[0],t0[1]+=td[1]) {
			double v = fq_warp2(fq1 * pz01[0].f(t0[0])),
			       ang = 0.5 * M_PI * pz01[1].f(t0[1]),
			       re = -v*sin(ang), im = v*cos(ang);
			if (zf) p = pzi + (nz++), p->zr = re, p->zi = im;
			else    p = pzi + (np++), p->pr = re, p->pi = im;
		}}
	log("bqs: np=%d, nz=%d", np, nz);
	PZ_END(bqs, 3);
}

static void cheb_calc_ri(double *re, double *im, double ripl, int n) {
	if (ripl<1e-10) { *re = *im = 1.0; return; }
	if (ripl>0.29) ripl = 0.29;
	double eps = 1.0 / (1.0-ripl);
	eps = sqrt(eps*eps-1.0);
	double nu = 0.5 * asinh(1.0/eps) / (double)n;
	double kk = 1.0 / cosh( 0.5*acosh(1.0/eps)/(double)n );
	*re = sinh(nu) * kk;
	*im = cosh(nu) * kk;
}

//? {{{!._chF}}}
//? Butterworth/Chebyshev filter
//? in - filter input
//? +-fq - cutoff freq. (>0: low-pass, <0:high-pass)
//? #pol - number of poles
//? ripl - passband ripple 0.0 ... 0.29 (0: Butterworth f.)
//? zr* - extra (low-p.) or missing (high-p.) zeros, 0...(#pol-1)
//? dcy - decay (out multiplied by exp(-t*dcy))
//? CPU usage depends on number of poles.
//? non-zero "zr*" turns the lp/hp filter to asymmetric bandpass

QUICK_PZFILT(ChebF) {
	int np = (int)lround(inb[1][0]), wf = 0;
	if (np<1) return 0; else if (np>30) np=30;
	int np_skip = (int)lround(inb[3][0]);
	double ripl = inb[2][0], step=0.5*M_PI/(double)np, fq = **inb, fqw = fq_warp(fq), cm = fqw*fqw;
	if (ripl < -1e-11) ripl = -ripl, wf = 2; else fq = fqw;
	set_n_d(np, inb[4][0]);
	RECF_PZItem * pzn = new RECF_PZItem[np];
	double ang = ((double)np-0.5)*step, re_fac, im_fac;
	cheb_calc_ri(&re_fac, &im_fac, ripl, np);
	if (fq<0.0) {
		fq = -fq;
		for (int i=0; i<np; i++) {
			double pr = -sin(ang) * fq * re_fac;
			double pi = cos(ang) * fq * im_fac;
			double c = fq*fq/(pr*pr+pi*pi);
			pzn[i].pr = c*pr; 
			pzn[i].pi = c*pi;
			pzn[i].zr = pzn[i].zi = (i>=np_skip) ? 0.0 : NAN;
			pzn[i].c = (i>=np_skip) ? 1.0 : cm;
			ang -= step;
		}
	} else {
		for (int i=0; i<np; i++) {
			pzn[i].pr = -sin(ang) * fq * re_fac;
			pzn[i].pi = cos(ang) * fq * im_fac;
			pzn[i].zr = pzn[i].zi = (i>=np_skip) ? NAN : 0.0;
			pzn[i].c = (i<np_skip) ? 1.0 : cm;
			ang -= step;
		}
	}
	if (np_skip<0) paraconj(pzn, np);
	rfpz_transform(m_ab, pzn, np, wf+1); return 0;
}

struct ChebPoleLs : BoxInst_B1 { static scf_t sc_ini, sc_f1; double * m_p; };

BX_SCALC(ChebPoleLs::sc_ini) {
	SCALC_BXI(ChebPoleLs); int np = bxi->m_arg, no = 2*np;
	double cr, ci, fq = inb[0][0], ripl = inb[1][0], step = 0.5*M_PI/(double)np,
	       ang = ((double)np-0.5)*step;
	cheb_calc_ri(&cr, &ci, ripl, np);
	double *p = (double*)bxi->alloc0(no*sizeof(double));
	if (fq<0.0) {
		fq = -fq; double fq2 = fq*fq;
		for (int i=0; i<np; i++, ang-=step) {
			double pr = -sin(ang)*fq*cr, pi = cos(ang)*fq*ci, pa2 = pr*pr + pi*pi,
			       c = fq2/pa2, pa = c * sqrt(pa2), aa = (-2.0/M_PI)*atan(pr/pi);
			p[2*i] = pa; p[2*i+1] = aa;
	}} else {
		for (int i=0; i<np; i++, ang-=step) {
			double pr = -sin(ang) * fq * cr, pi = cos(ang) * fq * ci,
			       pa2 = pr*pr + pi*pi, pa = sqrt(pa2), aa = (-2.0/M_PI)*atan(pr/pi);
			p[2*i] = pa; p[2*i+1] = aa;
		}}
	log_n("chls ("); for (int i=0; i<no; i+=2) log_n(" %.15g,%.15g", p[i],p[i+1]); log(" )");
	CALC_FW(sc_f1);
}

BX_SCALC(ChebPoleLs::sc_f1) { SCALC_BXI(ChebPoleLs); double *p = (double*)bxi->m_p0;
			      for (int i=0,no=2*bxi->m_arg; i<no; i++) outb[i][0]=p[i]; return 0; }

//? {{{!._arcF}}}
//? Harmonic arc(s of poles&zeroes) filter
//? in - filter input
//? fq1 - base frequency
//? [#p] - list of pole count for each frequency (fq1,...)
//? [#z] - list of zero count for each frequency (0,...)
//? pl0, pl1 - pole angles for low frq (from-to)
//? ph1, ph1 - pole angles for high frq
//? zl0, zl1, zh0, zh1 - same as above, for zeroes
//? [sc] - 6 scale types : [p0 p1 p z0 z1 z] (see map01)
//? stp1 - gap between pole & next zero fq. (*fq1)
//? stp2 - gap between zero & next pole fq. (*fq1)
//? dcy - decay (attenuation) mul. exp(-t*dcy)
//? (CPU usage is approx. linear with total number of poles)
//? ==> .!b.misc.map01
// in fq1 [np] [nz] pl0 pl1 ph0 ph1 zl0 zl1 zh0 zh1 [p0 p1 p z0 z1 z] stp1 stp2 dcy

QUICK_PZFILT(ArcRF) {
	NAN_UNPK_8(npole, inb[1], 0);
	NAN_UNPK_8(nzero, inb[2], 0);
	NAN_UNPK_8(sfun,  inb[11], 1);
	int nls = max_i(npole_n, nzero_n), ptot = 0, ztot = 0, pflg = 0, zflg = 0;
	double fq1 = inb[0][0], st1 = fq1 * inb[12][0], st2 = fq1 * inb[13][0], nfq = fq_warp(fq1),
	       sst = 1.0/((double)(nls-1)*(st1+st2)), ss1 = sst * st1, ss2 = sst * st2, sv = 0.0,
	       pfq[16], zfq[16], pss[16], psv[16], zss[16], zsv[16];
	Scale01 pzsc[4], psc[16], zsc[16];
	for (int j,i=0; i<4; i++) j=i+2*(i>>1), pzsc[i].set_all(inb[j+3][0], inb[j+5][0], sfun_x[i+(i>>1)]);
	int pcf = (nzero_x[0]<0) && (nzero_x[0]=0, 1);
	for (int n1, n2, i=0; i<nls; i++) {
		ptot += (n1=npole_x[i]); pflg |= !!n1 << (15-i);
		ztot += (n2=nzero_x[i]); zflg |= !!n2 << (15-i);
		if (i) { if (n2) zsc[i].set_all(pzsc[2].f(sv), pzsc[3].f(sv), sfun_x[5]),
				 zsv[i] = .5 / (double)n2, zss[i] = zsv[i]+zsv[i];
		 	 sv += ss2; zfq[i] = fq_warp(fq1); fq1 += st2; }
		if (n1) psc[i].set_all(pzsc[0].f(sv), pzsc[1].f(sv), sfun_x[2]), 
			psv[i] = .5 / (double)n1, pss[i] = psv[i]+psv[i];
		sv += ss1; pfq[i] = fq_warp(fq1); fq1 += st1;
	}
	int pix = ptot - 1, pzdiff = ptot - ztot; if (pzdiff<0) return RTE_FZZZ;
	RECF_PZItem * pzi = new RECF_PZItem[ptot];
	while (pflg) { BVFOR_JM(pflg) {
		int zi, i = 15 - j;
		pzi[pix].spp(pfq[i], psc[i].f(psv[i])); psv[i] += pss[i];
		if (!--npole_x[i]) pflg &= ~(1<<j);
		if (nzero_x[i]) zi = i;
		else if (nzero_x[i+1]) zi = i+1;
		else if (pzdiff) --pzdiff, zi = -1;
		else zi = __builtin_ffs(zflg), zi = zi ? 16-zi : -1;
		if (zi<0) { pzi[pix].zi = pzi[pix].zr = NAN; goto norm; }
		if (!--nzero_x[zi]) zflg &= ~(32768>>zi);
		if (zi) pzi[pix].spz(zfq[zi], zsc[zi].f(zsv[zi])), zsv[zi] += zss[zi]; 
		else pzi[pix].zi = pzi[pix].zr = 0.0;
norm:		pzi[pix].c = 1.0 / pzi[pix].eval1(0.0, nfq); --pix;
	}}
	if (pix != -1) log("ArcRF BUG!!: pix = %d (exp. -1)", pix);
	set_n_d(ptot, inb[14][0]);
	if (pcf) paraconj(pzi, ptot); rfpz_transform(m_ab, pzi, ptot); return 0;
}

static void bq_init(ANode * rn) {
	char nm[8]; memcpy(nm,"=bq*1\0 ", 8);
	qmk_box(rn, nm, QMB_A_BX(BiQuadSeqF), 1, 10, 1, "bqF", "i-o*R*1", 30,
			"in$fq1$[sc]$[md]$nfq$n1$f1<$>f1$q1<$>q1$n2$f2<$>f2$q2<$>q2$"
			  "n3$f3<$>f3$q3<$>q3$n4$f4<$>f4$q4<$>q4$n5$f5<$>f5$q5<$>q5",
			  "out", "uuu3%F");
	for (int i=2; i<6; i++) ++nm[4], qmk_box(rn, nm, QMB_A_BX(BiQuadSeqF), i, 5*i+5, 1, "bqsF", "1");
	nm[3] = '0', nm[4] = '1';
	qmk_box(rn, nm, QMB_A_BX(BiQuad1F), 1, 6, 1, "bqF", "i-i:o*R*1",
			4, "in$fq1$[#]$nfq", 2, 3, "f$q", "out", "uuu3%F");
	for (int i=2; i<14; i++) i==10 ? (++nm[3], nm[4]=48) : ++nm[4],
		qmk_box(rn, nm, QMB_A_BX(BiQuad1F), i, 2*i+4, 1, "bqF", "1");

}

static void pls_init(ANode * rn) {
	char nm[8]; memcpy(nm,"chb-ls2\0 ", 8);
	qmk_box(rn, nm, QMB_A_BX(ChebPoleLs), 2, 2, 36, "chL", "i*o-R*1", "fq$ripl", 18,
			"f1$q1$f2$q2$f3$q3$f4$q4$f5$q5$f6$q6$f7$q7$f8$q8$f9$q9", "uuu%F%");
	for (int i=3; i<10; i++) ++nm[6], qmk_box(rn, nm, QMB_A_BX(ChebPoleLs), i, 2, 32+2*i, "chL", "1");
}

static void eq_init(ANode * rn) {
	qmk_box(rn, "=eql-ls", QMB_A_BX(EqlzLsF),0,6,33,"eqLF", "i*R*", "in$fq1$step$wid$2div$[a]", "uuu3%F");
	char nm[8]; memcpy(nm,"=eql1\0", 6);
	qmk_box(rn, nm, QMB_A_BX(EqlzF), 1, 5, 33, "eqF", "i-R*1", 29,
			"in$fq1$f*1$wid1$amp1$f*2$wid2$amp2$f*3$wid3$amp3$f*4$wid4$amp4$f*5$wid5$amp5$"
			"f*6$wid6$amp6$f*7$wid7$amp7$f*8$wid8$amp8$f*9$wid9$amp9", "uuu3%F");
	for (int i=2; i<10; i++) 
		++nm[4], qmk_box(rn, nm, QMB_A_BX(EqlzF), i, 3*i+2, 33, "eqF", "1");
}

void b_filt_pz_init(ANode * rn) {
	qmk_box(rn, "=band", QMB_A_BX(SBandF), 0, 3, 1, "bF", "i*o*R*1", "in$fq$Q", "out", "uuu3%F");
	qmk_box(rn, "=cheb", QMB_A_BX(ChebF), 0, 6, 1, "chF", "1i-", 0x105, "+-fq$#pol$ripl$zr*$dcy");
	qmk_box(rn, "=arc", QMB_A_BX(ArcRF), 0, 16, 1, "arcF", "1+i*", "zza", "in$fq1$[#p]$[#z]$"
			"pl0$pl1$ph0$ph1$zl0$zl1$zh0$zh1$[sc]$stp1$stp2$dcy");
	qmk_box(rn, "=hpp", QMB_A_BX(HPPF), 0, 8,  1, "hppF", "1i*", "in$fq$Q$np$z/2p$stp$cfac$dcy");
	qmk_box(rn, "=hppsc", QMB_A_BX(HPPscF), 0, 12, 1, "hpsF", "1i*",
			"in$fq$np$z/2p$stp$Q0$Q1$Qs$c0$c1$cs$dcy");
	bq_init(qmk_dir(rn, "bq"));
	eq_init(qmk_dir(rn, "eql"));
	pls_init(qmk_dir(rn, "pls"));
}
