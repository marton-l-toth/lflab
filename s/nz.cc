#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glob.h"
#include "cfgtab.inc"
#include "util.h"
#include "box0.h"
#include "distr.h"

struct vzent_t { int iacc,rsrv; double x0, xu, yu; };
struct vzdat_t { dfun1_t f; vzent_t e[256]; };

rnd_dsc_t rnd_dsc[] = { // data gen. by ../t/nzar.c
	{0.001774861363360951, 'g', 0}, // x1:0.9999999999999704 iacc:0.958538
	{0.0005035812084291209, 'e', 8}, // x1:0.9999999999999536 iacc:0.948282
	{0.0002482277951613336, 'E', 9}, // x1:0.9999999999999356 iacc:0.971488
	{0.00618528540070983, 'c', 10}, // x1:0.9999999999999991 iacc:0.981215
	{0.004236443289877119, 'h', 11}, // x1:0.9999999999999996 iacc:0.969821
};

static uint64_t xoro_glob[2];
static int seed_glob;
static int iacc_ny[2];

// contains code from Sebastiano Vigna (Univ. of Milan)
// see his original code (public domain) at http://xoroshiro.di.unimi.it/

static uint64_t splitmix64(uint64_t *p) {
	uint64_t z = (*p += UINT64_C(0x9E3779B97F4A7C15));
	z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
	return z ^ (z >> 31);
}

static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

#define XORO_STEP(X0) (s1^=s0, (X0), s0=rotl(s0,55)^s1^(s1<<14), s1=rotl(s1,36))
void xoro_jump(uint64_t *s) {
        uint64_t x0 = 0, x1 = 0, s0 = s[0], s1 = s[1];
        for (uint64_t m, j = 0xbeac0467eba5facb; j; j>>=1) m=-(j&1), x1^=(s1&m), XORO_STEP(x0^=s0&m);
        for (uint64_t m, j = 0xd86b048b86aa9922; j; j>>=1) m=-(j&1), x1^=(s1&m), XORO_STEP(x0^=s0&m);
        s[0] = x0; s[1] = x1;
}

void xoro_s64(uint64_t *s, int k) {
	uint64_t s1 = (uint64_t)k;
	splitmix64(&s1); s[0] = splitmix64(&s1); s[1] = splitmix64(&s1); 
}

void gwn_debug(int step) {
	int n = iacc_ny[0], y = iacc_ny[1], tot = n+y; if (!tot) return;
	double prc = 100.0 / (double)tot;
	log("perfstat: iacc: tot:%d y:%d(%g%%) n:%d(%g%%)", tot, y, prc*(double)y,
								 n, prc*(double)n); }

//? {{{!._nz}}}
//? Simple noise generator
//? ty - noise type
//? sd - random seed
//? ---
//? type: 0:normal(Gaussian) 1:linear 2:2*lin 3:3*lin 4:lin^2
//? 5:lin^3 6:lin^4 7:lin^6 8:exponential (signed) 
//? 9: exp. (unsigned) 10:half-circle 11:(x^2-1)^2
//? ---
//? seed: rounded to integer, <0 and NaN values are reserved
//? 1...2147483647 gives a reproducible (*) random output
//? 0 gives the default random seed, configurable only from
//? command line with -rseed and -rstep (see lflab -h)
//? ---
//? (*) for "sample-discarding" generators (0,8,9,10,11) very
//? small rounding differences may lead to altered output data.
//? (insertion/deletion of a few samples)

class NZ0Box : public BoxInst {
	public: static void mk_dsc(int ix, int ty, double ar);
		NZ0Box() : m_psc(sc_ini) {}
		void ini(const double *pty, const double *psd);
		CALC_TODO;
	protected:
		static scf_t sc_ini, sc_vz, sc_x1, sc_x2a, sc_x3a,
			     sc_x2m, sc_x3m, sc_x4m, sc_x6m;
		static void* m0_dsctab[16];
		static const sc_t m0_sftab[16];
		uint64_t m_xoro[2];
		void * m_dsc;
};

class NZ0Buf {
	public:
		NZ0Buf() : m_pv(m_v), m_ix(15) {}
		inline void ini(const double *pty, const double *psd) { m_b.ini(pty,psd); }
		inline double v() { return (m_ix<15) ? m_v[++m_ix]
						     : (m_b.calc(0,0,&m_pv,16), m_v[m_ix=0]); }
	protected:
		NZ0Box m_b;
		double m_v[16], *m_pv;
		int m_ix;
};

//? {{{!._nz2}}}
//? Two simple noise generators, multiplied
//? ty1, ty2 - noise types
//? sd1, sd2 - seed values
//? ==> .!b.noise.nz0 -- details for simple noise gen.

class NZ0x2Box : public BoxInst {
	public: NZ0x2Box() : m_psc(sc0) {}
		CALC_TODO;
	protected:
		static scf_t sc0, sc1;
		NZ0Box m_b0, m_b1;
};

//? {{{!._rpt}}}
//? Randomized pulse train
//? fq1 - base fq
//? fnz - fq noise -- fq is fq1 * (1.0 + fnz*nz1())
//? vnz - vol noise -- imp vol is (1.0 - vnz) + vnz*nz2()
//? rsrv - reserved for future use
//? ty1,ty2,sd1,sd2 - parameters to nz1() and nz2()
//? ==> .!b.noise.nz0 -- details for simple noise gen.

class RndPt : public BoxInst {
	public: RndPt() : m_psc(sc0), m_acc(1.00000001) {}
		CALC_TODO;
	protected:
		static scf_t sc0, sc1;
		NZ0Buf m_b0, m_b1;
		double m_acc, m_fmul;
};

void NZ0Box::ini(const double *pty, const double *psd) {
	int ty = (int)lround(*pty) & 15, seed = (int)lround(*psd);
	if      (seed)		 xoro_s64(m_xoro, seed);
	else if (CFG_RND_STEP.i) xoro_s64(m_xoro, seed_glob++);
	else  memcpy(m_xoro, xoro_glob, 16), xoro_jump(xoro_glob);
	m_dsc = m0_dsctab[ty];  m_psc = m0_sftab [ty];
}

BX_SCALC(NZ0Box::sc_ini) { SCALC_BXI(NZ0Box); bxi->ini(inb[0],inb[1]); return (bxi->m_psc)(bxi,0,0,outb,n); }

#define NZ0_START SCALC_BXI(NZ0Box); uint64_t r, s0 = bxi->m_xoro[0], s1 = bxi->m_xoro[1]; double *to = *outb
#define NZ0_END   bxi->m_xoro[0] = s0; bxi->m_xoro[1] = s1; return 1

BX_SCALC(NZ0Box::sc_vz) {
	NZ0_START, *tol = to+n;
	const vzdat_t *dsc = (const vzdat_t*)bxi->m_dsc;
	const vzent_t *z0 = dsc->e;  dfun1_t fp = dsc->f;
	for(;;){r = s0 + s1;
		int k = (r>>56); const vzent_t * z = z0 + k;
		int xi, yi; // xi = (int)(r>>25) & 0x7fffffff, yi = r&0x1ffffff;
		XORO_STEP((xi = (int)(r>>25) & 0x7fffffff, yi = r&0x1ffffff));
		double x = z->x0 + (double)xi*z->xu;
		if ((
#ifdef NZ_DEBUG
		(k=(yi<z->iacc), ++iacc_ny[k], k)
#else
                (yi<z->iacc)
#endif
				|| (double)yi * z->yu < (*fp)(x)) && (*to=x,++to>=tol)) break;
        }
	NZ0_END;
}

#define I52_TO_D(X,E) (r>>=12, r|=(0x3ff0000000000001L+((int64_t)(E)<<52)), memcpy(&(X), &r, 8))
#define XORO_STEP_2D(X,E,X0) (r = s0+s1, XORO_STEP(I52_TO_D(X,E)), (X)+=(X0))

BX_SCALC(NZ0Box::sc_x1) {
	NZ0_START;
	for (int i=0; i<n; i++) XORO_STEP_2D(to[i], 1, -3.0);
	NZ0_END;  }

BX_SCALC(NZ0Box::sc_x2a) {
	NZ0_START, z;
	for (int i=0; i<n; i++) XORO_STEP_2D(z, 0, -3.0), XORO_STEP_2D(to[i], 0, z);
	NZ0_END;  }

BX_SCALC(NZ0Box::sc_x3a) {
	NZ0_START, y, z;
	for (int i=0; i<n; i++) XORO_STEP_2D(z, 1, -9.0), XORO_STEP_2D(y, 0, z),
				XORO_STEP_2D(z, 1,    y), to[i] = .3333333333333333 * z;
	NZ0_END;  }

BX_SCALC(NZ0Box::sc_x2m) {
	NZ0_START, y, z;
	for (int i=0; i<n; i++) XORO_STEP_2D(z, 1, -3.0), XORO_STEP_2D(y, 1, -3.0), to[i] = y*z;
	NZ0_END;  }

BX_SCALC(NZ0Box::sc_x3m) {
	NZ0_START, y, z;
	for (int i=0; i<n; i++) XORO_STEP_2D(z, 1, -3.0), XORO_STEP_2D(y, 1, -3.0), z*=y,
				XORO_STEP_2D(y, 1, -3.0), to[i] = y*z;
	NZ0_END;  }

BX_SCALC(NZ0Box::sc_x4m) {
	NZ0_START, y, z;
	for (int i=0; i<n; i++) XORO_STEP_2D(z, 1, -3.0),       XORO_STEP_2D(y, 1, -3.0), z*=y,
				XORO_STEP_2D(y, 1, -3.0), z*=y, XORO_STEP_2D(y, 1, -3.0), to[i] = y*z;
	NZ0_END;  }

BX_SCALC(NZ0Box::sc_x6m) {
	NZ0_START, y, z;
	for (int i=0; i<n; i++) XORO_STEP_2D(z, 1, -3.0),       XORO_STEP_2D(y, 1, -3.0), z*=y,
				XORO_STEP_2D(y, 1, -3.0), z*=y, XORO_STEP_2D(y, 1, -3.0), z*=y,
				XORO_STEP_2D(y, 1, -3.0), z*=y, XORO_STEP_2D(y, 1, -3.0), to[i] = y*z;
	NZ0_END;  }

const BoxInst::sc_t NZ0Box::m0_sftab[16] = { sc_vz, sc_x1, sc_x2a, sc_x3a, sc_x2m, sc_x3m, sc_x4m, sc_x6m,
	sc_vz, sc_vz, sc_vz, sc_vz, sc_x1, sc_x1, sc_x1, sc_x1 };
void * NZ0Box::m0_dsctab[16];

void NZ0Box::mk_dsc(int ix, int ty, double ar) {
	static const double xmul = 1.0/(double)((1u<<31u)), ymul = 1.0/(double)(1<<25);
	DISTR_NF(ty); double x0, x1 = 0.0;
	vzdat_t * zd = (vzdat_t*)(m0_dsctab[ix] = nf_alloc(sizeof(vzdat_t)));  zd->f = f; 
	for (int i=0; i<n; i++) {
		x0 = x1; x1 += ar / (*f)(x0);
		vzent_t * p = zd->e + i;
		p->xu = xmul*(x1-x0); p->x0 = x0 + .5*p->xu;
		p->yu = ymul * (*f)(p->x0);
		p->iacc = (int)lround( (*f)(x1) / p->yu );
	}
	if (ty&32) { for (int i=0; i<n; i++) {
		vzent_t *p = zd->e+i, *q = p+n;
		memcpy(q, p, n); q->x0 = -q->x0; q->xu = -q->xu;
	}}}

BX_SCALC(NZ0x2Box::sc1) { 
	SCALC_BXI(NZ0x2Box); double tmp[n], *q0 = *outb, *q1 = tmp;
	bxi->m_b0.calc(0, 0, &q0, n);
	bxi->m_b1.calc(0, 0, &q1, n);  for (int i=0; i<n; i++) q0[i] *= q1[i];  return 1;  }

BX_SCALC(NZ0x2Box::sc0) { 
	SCALC_BXI(NZ0x2Box); bxi->m_b0.ini(inb[0], inb[2]); bxi->m_b1.ini(inb[1], inb[3]);
	return (*(bxi->m_psc = &sc1)) (bxi, inflg, inb, outb, n);  }

#define RNDPT1 (v = vnzp[i&vnzmsk], v=1.0-v+v*bxi->m_b0.v(), \
		fmul = sample_length*(1.0+bxi->m_b1.v()*fnzp[i&fnzmsk]))
BX_SCALC(RndPt::sc1) {
	SCALC_BXI(RndPt);
        int i; double d, fq, v, *fqp = *inb, *q = *outb, fmul = bxi->m_fmul, acc = bxi->m_acc;
        PREP_INPUT(fnz, 1); PREP_INPUT(vnz, 2);
        if (inflg&1) for (i=0; i<n; i++) q[i] = acc<1.0 ? 0.0 : (RNDPT1, acc-=1.0, v), acc += fqp[i]*fmul;
	else for (i=0,fq=*fqp,d=fq*fmul; i<n; i++) q[i] = acc<1.0 ? 0.0 : (d=fq*RNDPT1, acc-=1.0, v), acc+=d;
        bxi->m_fmul = fmul; bxi->m_acc = acc;
        return 1;
}

BX_SCALC(RndPt::sc0) {
	SCALC_BXI(RndPt); bxi->m_b0.ini(inb[4], inb[6]); bxi->m_b1.ini(inb[5], inb[7]);
	return (*(bxi->m_psc = &sc1)) (bxi, inflg, inb, outb, n);  }

void b_noise_init(ANode *rn) {
	int seed = CFG_RND_SEED.i; if (!seed) seed = zero_sec;
	xoro_s64(xoro_glob, seed);
	for (int i=0; i<ASIZ(rnd_dsc); i++)
		NZ0Box::mk_dsc(rnd_dsc[i].ix, rnd_dsc[i].ty, rnd_dsc[i].ar);
	qmk_box(rn, "nz0", QMB_ARG0(NZ0Box), 0, 2, 33, "nz", "i*o*R*1", "ty$sd", "out", "zzzOOO");
	qmk_box(rn, "nz0*2", QMB_ARG0(NZ0x2Box), 0, 4, 33, "nz2", "1i*", "ty1$ty2$sd1$sd2");
	qmk_box(rn, "rnd-pt", QMB_ARG0(RndPt), 0, 8, 33, "rpt", "i*o*R*", "fq1$fnz$vnz$rsrv$ty1$ty2$sd1$sd2",
			"out", "z%%OOO");
}
