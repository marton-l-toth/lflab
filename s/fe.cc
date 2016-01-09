#include "util.h"
#include "box0.h"
#include "util2.h"
#include "glob.h"

struct FE1 { double y, v, cs[4], w, wgt; };
struct FEblk {
	static FEblk * mk0(int m, int i0, int iw);
	int m, i0, iw;
	FE1 *p; 
	double *cs1, *iv, lastcs[2];
	void cs_ini_1(int i);
	void cs_ini_x(int i, double x);
	void iv_ini(double z0, double z1);
	void geom_ini(double *x012, double *y012, double pkx, double pky);
};

//? {{{!._fe}}}
//? Linear finite-element based filter (string/membrane/etc.)
//? in - filter input
//? #s - number of segments (CPU usage approx. linear with this)
//? opt - reserved for future use
//? x2, x1, x0 - coeff. for surface (quadratic) function
//? y2, y1, y0 - coeff. for thickness (q.) function
//? su0, su1 - suspension coeff. (1.0 - normal, 0.0 - floating)
//? ipos, iwid - input position & width (0..1)
//? pkx, pky - pickup position (0..1) and distance
//? F/y - elastic force coeff.
//? F/v - nonelastic force coeff.
//? dmp - dampening
//? F/y, F/v and dmp can be non-constant

//  0    0    1   2  g1/3 1  2  3  4  5  6  7  8   9   10   11   12  13  p1[57] 1  2 
// (in | hhv hhp hhf) #s opt x2 x1 x0 y2 y1 y0 su0 su1 ipos iwid pkx pky F/y   ie dmp
class FEBox : public BoxInst {
        public:
                FEBox(int flg) : m_flg(flg), m_blk(0) {}
                virtual ~FEBox() { if (m_blk) free(m_blk); }
                virtual int calc(int inflg, double** inb, double** outb, int n);
		inline int geom_pos() const { return 1 + (m_flg & 2); }
		inline int par_pos() const { return 15 + (m_flg & 2); }
        protected:
		FEblk * ini_blk(double **inb);
		int m_flg;
		FEblk * m_blk;
};

static inline void ply_qi(double *q, const double *p, int n) {
	for (int i=0; i<n; i++) q[i] = p[i] * smalldiv63[i+1]; }
static inline double ply_ev(const double *p, int n, double x) {
	double y = p[n-1]; for (int i=n-2; i>=0; i--) y*=x,y+=p[i]; return y; }

void FEblk::geom_ini(double *x012, double *y012, double pkx, double pky) {
	double w0_4[5]; memset(w0_4, 0, 40); pky *= pky;
	for (int i=0; i<3; i++) for (int j=0; j<3; j++) w0_4[i+j] += x012[i]*y012[j];
	double c0=0.0, ts=0.0, hstp = .5/(double)(m-2),   sprf[3]; ply_qi(sprf, x012, 3);
	double x =0.0, tw=0.0, c1 = ply_ev(w0_4, 5, 0.0), wprf[5]; ply_qi(wprf, w0_4, 5);
	double stp = hstp+hstp, csmax = 0.0, swtot = 0.0, *q = cs1;
	for (int i=0; i<m-2; i++, q+=2) {
		x += stp; c0 = c1; c1 = ply_ev(w0_4, 5, x);
		double ts2 = x * ply_ev(sprf, 3, x), surf = ts2 - ts, d2 = x - hstp - pkx, u2,
		       tw2 = x * ply_ev(wprf, 5, x), wgt  = tw2 - tw, u0 = c0/wgt, u1 = c1/wgt;
		q[0] = u0; q[1] = u1; p[i+1].wgt = wgt;
		if (     (u2=u0   +u1)>csmax) csmax = u2;
		if (i && (u2=q[-1]+u0)>csmax) csmax = u2;
		swtot += (p[i+1].w = surf / (d2*d2+pky)); p[i+1].y = p[i+1].v = 0.0;
	}
	swtot = 1.0 / swtot; for (int i=1; i<m-1; i++) p[i].w *= swtot;
	double csmul[2]; csmul[0] = -(csmul[1] = 1.0 / csmax);
	for (int i=0; i<2*m; i++) cs1[i] *= csmul[i&1];
}

FEblk * FEblk::mk0(int m, int i0, int iw) {
	static const int bsiz = (sizeof(FEblk)+7) & 0xfff8;
	int dtot = 8 * (iw + 2*m - 4);
	char *raw = (char*)malloc(bsiz + dtot + m*sizeof(FE1)), *pdbl = raw + bsiz;
	FEblk * r = (FEblk*)raw; r->m = m; r->i0 = i0; r->iw = iw;
	FE1   * p = r->p  = (FE1*)(pdbl + dtot);
	r->lastcs[0] = r->lastcs[1] = NAN; p[0].y = p[0].v = p[m-1].y = p[m-1].v = 0.0;
	r->iv = (double*)pdbl; r->cs1 = r->iv + iw; return r;
}

void FEblk::cs_ini_1(int i) { if (!approx_cmp(lastcs[i], 1.0)) return; else lastcs[i] = 1.0;
	for (int ii=2*i,j=0; j<m-2; j++) p[j+1].cs[ii]=cs1[2*j], p[j+1].cs[ii+1]=cs1[2*j+1]; }
void FEblk::cs_ini_x(int i, double x) { if (!approx_cmp(lastcs[i], x*=x)) return; else lastcs[i] = x;
	for (int ii=2*i,j=0; j<m-2; j++) p[j+1].cs[ii]=x*cs1[2*j], p[j+1].cs[ii+1]=x*cs1[2*j+1]; }

void FEblk::iv_ini(double z0, double z1) { double j = (double)i0, y, a = 0.0;
	for (int i=0; i<iw; i++, j+=1.0) y=(j-z0)*(j-z1), a += p[i+i0].wgt * (iv[i]=y*y);
	a = 1.0/a; for (int i=0; i<iw; i++) iv[i] *= a; }

FEblk * FEBox::ini_blk(double **inb) {
	double **pg = inb + 1 + (m_flg&2);
	double x012[3]; for (int i=0; i<3; i++) x012[i] = pg[4-i][0];
	double y012[3]; for (int i=0; i<3; i++) y012[i] = pg[7-i][0];
	int m0 = (int)lround(**pg), m = ((unsigned int)(m0-1)<250u) ? m0+2 : (m0>1?252:3),
	    ipm0 = (int)ceil(pg[10][0]*(double)(m-2)), ipm = ipm0<1 ? 1 : (min_i(ipm0, m-2));
	if (m_flg&2) {
		(m_blk = FEblk::mk0(m, ipm, 0))->geom_ini(x012, y012, pg[12][0], pg[13][0]);
	} else {
		double ipmd = (double)ipm, ihw = 0.5 * pg[11][0], z0 = ipmd - ihw, z1 = ipmd + ihw;
		int z0i = (int)ceil (z0); if (z0i<1  ) z0i = 1;   else if (z0i>ipm) z0i = ipm;
		int z1i = (int)floor(z1); if (z1i>m-2) z1i = m-2; else if (z1i>ipm) z1i = ipm;
		if (z1i<z0i) z0i = z1i = ipm;
		m_blk = FEblk::mk0(m, z0i, z1i-z0i+1);
		m_blk->geom_ini(x012, y012, pg[12][0], pg[13][0]);
		m_blk->iv_ini(z0, z1);
	}
	m_blk->cs1[0] *= pg[8][0]; m_blk->cs1[2*m-5] *= pg[9][0];
	return m_blk;
}

#define FC_00(IX) case IX: { 
#define FC_0a(IX) FC_00(IX) PREP_INPUT(in,0);
#define FC_0h(IX) FC_00(IX) PREP_INPUT(hhv,0); PREP_INPUT(hhp,1); PREP_INPUT(hhf,2); FE1 *qhh=p+blk->i0;
#define FC_2 for (int j, i=0; i<n; i++) { double a=0.0;  if (parf&4) att = att2mul(*(attp++));
#define FC_3 for (j=0,k=blk->i0; j<blk->iw; j++) p[k+j].v += inp[i&inmsk] * blk->iv[j];
#define FC_Tc(T) for (j=0; j<m-1; j++) t1[j] = p[j+1].T-p[j].T;
#define FC_Tv(T,I) double k##T = par[I][i]; k##T *= k##T; for (j=0; j<m-1; j++) t1[j] = k##T*(p[j+1].T-p[j].T);
#define FC_su1x (q->v = att * (q->v + t1[j-1] * q->cs[0] + t1[j] * q->cs[1]) )
#define FC_5vy for (j=1, q=p+1; j<m-1; j++,q++) a += q->w * (q->y += FC_su1x);
#define FC_5v  for (j=1, q=p+1; j<m-1; j++,q++) FC_su1x;
#define FC_7   for (j=1, q=p+1; j<m-1; j++,q++) a+= q->w*(q->y+= (q->v+= t1[j-1]*q->cs[2]+t1[j]*q->cs[3]));
#define FC_8 double vd=hhvp[i&hhvmsk]-qhh->v, hhp0=hhpp[i&hhpmsk], hhp=(m_flg&1)?hhp0*hhfp[i&hhfmsk] : hhp0;\
	if (vd>0.0) vd<hhp ? (m_flg|=1) : (vd= hhp, m_flg&=~1); \
	else       -vd<hhp ? (m_flg|=1) : (vd=-hhp, m_flg&=~1);   qhh->v+=vd; qhh->y+=vd;
#define FC_Z to[i] = a; }} return 1;

int FEBox::calc(int inflg, double** inb, double** outb, int n) {
	FEblk * blk = m_blk ? m_blk : ini_blk(inb);
	int m = blk->m, k = par_pos(), parf = inflg>>k, ty = parf&1;
	double s0, t1[m], *to = *outb, **par = inb+k, *attp = par[2], att = (parf&4) ? 0.0 : att2mul(*attp);
	FE1 *q,*p = blk->p;
	if (ty) blk->cs_ini_1(0); else blk->cs_ini_x(0, par[0][0]);
	if (parf&2) { if (inb[1]!=zeroblkD) ty+=4, blk->cs_ini_1(1); }
	else	    { if (fabs(s0=par[1][0])>1e-280) ty+=2, blk->cs_ini_x(1,s0); }
	if ((inflg&1)?*inb!=zeroblkD:**inb>1e-280) ty += 6 + 3*(m_flg&2);
	switch(ty) {
		FC_00( 0) FC_2      FC_Tc(y)   FC_5vy 			   FC_Z;
		FC_00( 1) FC_2      FC_Tv(y,0) FC_5vy 			   FC_Z;
		FC_00( 2) FC_2      FC_Tc(y)   FC_5v  FC_Tc(v)   FC_7	   FC_Z;
		FC_00( 3) FC_2	    FC_Tv(y,0) FC_5v  FC_Tc(v)   FC_7	   FC_Z;
		FC_00( 4) FC_2      FC_Tc(y)   FC_5v  FC_Tv(v,1) FC_7	   FC_Z;
		FC_00( 5) FC_2      FC_Tv(y,0) FC_5v  FC_Tv(v,1) FC_7	   FC_Z;
		FC_0a( 6) FC_2 FC_3 FC_Tc(y)   FC_5vy 			   FC_Z;
		FC_0a( 7) FC_2 FC_3 FC_Tv(y,0) FC_5vy 			   FC_Z;
		FC_0a( 8) FC_2 FC_3 FC_Tc(y)   FC_5v  FC_Tc(v)   FC_7      FC_Z;
		FC_0a( 9) FC_2 FC_3 FC_Tv(y,0) FC_5v  FC_Tc(v)   FC_7	   FC_Z;
		FC_0a(10) FC_2 FC_3 FC_Tc(y)   FC_5v  FC_Tv(v,1) FC_7	   FC_Z;
		FC_0a(11) FC_2 FC_3 FC_Tv(y,0) FC_5v  FC_Tv(v,1) FC_7	   FC_Z;
		FC_0h(12) FC_2      FC_Tc(y)   FC_5vy 		      FC_8 FC_Z;
		FC_0h(13) FC_2      FC_Tv(y,0) FC_5vy 		      FC_8 FC_Z;
		FC_0h(14) FC_2      FC_Tc(y)   FC_5v  FC_Tc(v)   FC_7 FC_8 FC_Z;
		FC_0h(15) FC_2	    FC_Tv(y,0) FC_5v  FC_Tc(v)   FC_7 FC_8 FC_Z;
		FC_0h(16) FC_2      FC_Tc(y)   FC_5v  FC_Tv(v,1) FC_7 FC_8 FC_Z;
		FC_0h(17) FC_2      FC_Tv(y,0) FC_5v  FC_Tv(v,1) FC_7 FC_8 FC_Z;
		default: log("BUG: fe/calc: invalid ty %d", ty); return RTE_BUG;
	}}

static const char * fe_inm2 = "#s$opt$x2$x1$x0$y2$y1$y0$su0$su1$ipos$iwid$pkx$pky$F/y$F/v$dmp";
void b_filt_fe_init(ANode * rn) {
	qmk_box(rn, "=fe",  QMB_ARG1(FEBox), 0, 18, 1, "_fe", "*i-i-", "0zz%II", 
			1, "in", 0x111, fe_inm2);
	//qmk_box(rn, "fe/hh", QMB_ARG1(FEBox), 2, 20, 1, "_fe", "*i-i-", "0zz%II",   TODO
	//		3, "hhv$hhp$hhf", 0x311, fe_inm2);
}
//  0    0    1   2  g1/3 1  2  3  4  5  6  7  8   9   10   11   12  13  p1[57] 1  2 
// (in | hhv hhp hhf) #s opt x2 x1 x0 y2 y1 y0 su0 su1 ipos iwid pkx pky ffac  st dmp
