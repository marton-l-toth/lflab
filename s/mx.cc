#include "mx.h"
#include "util.h"
#include "tab7.h"
#include "box0.h"
#include "glob.h"
#include "wrap.h"

class BoxGen; 
int box_mxc_notify(BoxGen *p, int ky, int flg), trk_cut_time(BoxGen *p, int t);
typedef void (*tfunini_t) (double *, int, int);
typedef void (*tfunmul1_t) (double *, double *, int);
typedef void (*tfunmul2_t) (double *, double *, double *, int);

#define TFUN_ENT(c, nm) {c,{ wrtfun_1_##nm, wrtfun_2_##nm, wrtfun_i_##nm }}
#define TFUN_DEF(nm, iex, uex) \
void wrtfun_i_##nm (double *q, int n, int ud) { (void)(iex); } \
void wrtfun_1_##nm (double *q, double *x, int n) { for(int i=0;i<n;x[i++]*=*q,(uex)); } \
void wrtfun_2_##nm (double *q, double *x, double *y, int n) { \
	for(int i=0; i<n; x[i]*=*q, y[i]*=*q, (uex), i++); }
#define TFUN_S (1-2*(ud&1))
#define MXLOG if (debug_flags & DFLG_MX) log
#define MXLOGN if (debug_flags & DFLG_MX) log_n
#define QLOGV(S) ((debug_flags & DFLG_MX) ? log(S) : (void)0)

#define MXIARG(V,I,C) MxItem * V = mx_ptr(I); \
	if (V->m_ty!=*#C) return V->m_ty ? MXE_EXP##C : MXE_NOSUCH

struct TFunDesc { tfunmul1_t f1; tfunmul2_t f2; tfunini_t fi; };

struct TFunTab {
	TAB7_DECL(f, TFunDesc);
	static int ini(double *to, int c, int n, int ud) {
		if (c>96) c-=32, ud+=2;  if ((c=f_ent_i(c))) m0_f_tab[c].f.fi(to, n, ud); return c; }
	static inline void tf12(int f, double *a, double *x, double *y, int n) {
		if(y) (m0_f_tab[f].f.f2) (a, x, y, n);
		else  (m0_f_tab[f].f.f1) (a, x,    n); }
};

static const double 
	tf_exp_q0[4] = {1.45519152283669e-11, 1.0, 3.814697265625e-06, 1.0},
	tf_exp_xf[2] = {24.953298500158, 12.476649250079},
	tf_sgm_q1[4] = {68719476735.0, 1.45519152285786e-11, 262143.0, 3.81471181759574e-06},
	tf_sgm_xf[2] = {49.906597000287, 24.9532908707489},
	tf_gau_xf[2] = {4.99532766694618, 3.53223006754642};
	
TFUN_DEF(lin, (ud&=1, q[0]=(double)ud, q[1]=(double)(1-2*ud)/(double)n), (q[0]+=q[1]))
TFUN_DEF(sq , (ud&=1, q[0]=(double)ud, q[1]=1.0/(double)n, q[1]*=q[1], q[2]=q[1]+q[1],
	       ud && (q[1]*=(double)(1-2*n))),
	      (q[0]+=q[1], q[1]+=q[2]))
TFUN_DEF(exp, (q[0]=tf_exp_q0[ud], q[1]=exp(tf_exp_xf[ud>>1]/(double)(n*TFUN_S))),
	      (q[0]*=q[1]))
TFUN_DEF(sgm, (q[1] = tf_sgm_q1[ud], q[2] = exp(tf_sgm_xf[ud>>1]/(double)((-n)*TFUN_S)),
	       q[0] = 1.0/(1.0+q[1])),
	      (q[0] = 1.0 / (1.0 + (q[1] *= q[2]))))
TFUN_DEF(gau, (q[0] = tf_gau_xf[ud>>1]/(double)(n+1), q[0] *= -q[0], q[2] = exp(2.0*q[0]),
	       q[1] = exp((ud&1) ? q[0] : (double)(-1-2*n)*q[0]), q[0] = tf_exp_q0[ud]),
	      (q[0]*=q[1], q[1]*=q[2]))

TAB7_DEF(TFunTab, f) { {'-',{0, 0, 0}},
	TFUN_ENT('L',lin), TFUN_ENT('Q',sq ), TFUN_ENT('E',exp), TFUN_ENT('S',sgm),
	TFUN_ENT('G',gau),
	{0, {0,0,0}}} ;

struct MxItem;

struct MxInP {
	double ** pp;
	union {
		struct { double *p[8], v[8]; } p8v8;
		struct { MxItem *p2; double v[15]; } iv15;
	} u;
	void ini(const double *p, int n);
	void del2();
	int fi_cmp(const double *p, int n);
};

struct MxTrg {
	double *pt[2], *pp5[6];
	int siz, lrc01, fdly; // lrc01: 32-chkdelay
	MxTrg(int n, int f, int flt, double *tmp, double *o0, double *o1);
	inline void prep(double *q, int dly, int ae) {  for (int i=0 ; i<dly; i++) q[i] = 0.0;
							for (int i=ae; i<siz; i++) q[i] = 0.0; }
	int addb(int flg, int dly, int n, const double * v2);
	int ret();
};

struct MxItBox {
	BoxInst * bx;
	int t, t1, t2, t3, delay;
	unsigned short ocfg, ctk, ctp, ctn; 
	char up, dn, abflg; unsigned char trec; // abf: (1,2:ou) 4:dir 8:ste 16:flt 32:vol=1
	double upx[3], dnx[3], v[2];
	MxInP in;
};

struct MxItFilt {
	BoxInst *bx[2];
	BoxModel *mdl;
	double vlim;
	int tlim, t; int abflg; // abf+: 3:mode 16384:running
	double v[2], v1, lr;
	unsigned short ocfg, up;
	MxInP in;
};

struct MxItCtrl {
	BoxGen * bref;
	int nb, uk;
	unsigned short bx[100];
};

struct MxItem {
	unsigned short m_id, m_pv, m_nx;
	char m_ty, m_x8; 
	union {
		struct MxItem * nxfr;
		double v[31];
		double *p[31];
		struct { double *p[23], v[8]; } pv;
		MxItBox b;
		MxItFilt f;
		MxItCtrl c;
		struct { char flg, rsrv[7]; short mi[120]; } r;
		struct { BoxModel * m; short fi[120]; } m;
		struct { BoxInst * bxi; int flg; unsigned char dat[224]; } l;
	} u;
	void b_ini(const char * updnnino, const double * arg, int delay, int tlim);
	int b_calc(MxTrg * mtg, int n);
	int b_evp(double *p, double *q, int n);
	void b_ccut();
	void b_fbv(MxItem * up, int ocf, double lr0, double v0);
	MxItem * r_mfind(BoxModel * mdl, int force);
	MxItem * m_ffind(int ni, double *arg, int force, int up);
	void rf_bfor(MxTrg * mtg, int n);
	int f_calc(MxTrg * mtg, int n);
	int r_calc(double *to1, double *to2, int n, int f);
	void rf_clear();
	void r2_clear();
	void f_ini_1(double lr);
	void f_ini_2();
	void f_1to2() { u.f.abflg |= 10; u.f.abflg &= ~36; u.f.bx[1] = u.f.mdl->mk_box(); }
	void f_reset() { rf_clear(); u.f.abflg=0; delete(u.f.bx[0]); delete(u.f.bx[1]); u.f.bx[0]=u.f.bx[1]=0; }
	int r_isemp() const { return !m_x8 && m_nx==m_id; }
	int c_find(int k, int f = 0);
	int c_ins(MxItem * that, int k);
	int c_cut(MxItem * that);
	int c_stop_j(int j, int f);
	void c_cut1j(int j) { int n = --u.c.nb - j; if (n) memmove(u.c.bx+j,u.c.bx+j+1,2*n); }
	void c_sane();
	void debug();
};

static MxItem *mx_ptab[512], *mx_fptr = 0;
static int mx_ptix = 0;
static int mx_live0 = 0;

void MxItem::debug() {
	log_n("(%c/%x", m_ty ? m_ty : '0', m_id); switch(m_ty) {
		case 'l': log_n("%c%x", 43+2*!u.l.bxi, u.l.flg); break;
		default: break;
	} log_n(")");
}
/////// trk rec /////////////////////////////////////////////////////

static BoxGen * trec_bx[256];
static unsigned short trec_bi[256];
static int trec_alloc() { int k = *trec_bi; return k ? (*trec_bi=trec_bi[k], k) : MXE_REC255; }
static void trec_free(int j) { trec_bx[j] = 0; trec_bi[j] = *trec_bi; *trec_bi = j; }
static void trec_fin_1(int j) { wrap_set_trec(trec_bx[j], 0); trec_free(j); }
static void trec_fin_2(MxItem * b, int j) { b->u.b.trec = 0, trec_fin_1(j); }
static void trec_fin_t(MxItem * b, int j, int t) {
	trk_cut_time(trec_bx[j],t); b->u.b.trec = 0; trec_free(j); }

/////// alloc / aux //////////////////////////////////////////////////////

static MxItem * mx_newblk() {
	if (mx_ptix == 512) return 0;
	MxItem * p = mx_ptab[mx_ptix] = (MxItem*)alloc_32k();
	for (int i=0, k=mx_ptix<<7; i<128; i++) p[i].m_id=i+k, p[i].m_ty=0, p[i].u.nxfr=p+i+1;
	++mx_ptix;
	p[127].u.nxfr = 0; return mx_fptr = p;
}

static inline MxItem * mx_alloc(char ty, unsigned short pv=0, unsigned short nx=0) {
	MxItem * p = mx_fptr ? mx_fptr : mx_newblk();
	MXLOG("mx_alloc: %p (0x%x, %c)", p, p->m_id, ty);
	return mx_fptr = p->u.nxfr, p->m_ty = ty, p->m_pv = pv, p->m_nx = nx, p; }

static inline void mx_free(MxItem *p) { MXLOG("mx_free: %p (0x%x, %c)", p, p->m_id, p->m_ty); p->m_ty = 0; p->u.nxfr = mx_fptr; mx_fptr = p; }
static inline MxItem * mx_ptr(int i) { return mx_ptab[(i>>7)&511] + (i&127); }

void MxInP::ini(const double *p, int n) {
	if (n<9) {
		memcpy(u.p8v8.v, p, 8*n); pp = u.p8v8.p;
		for(int i=0; i<n; i++) pp[i] = u.p8v8.v+i;
	} else if (n<24) {
		MxItem *p2 = u.iv15.p2 = mx_alloc('I');
		int n2 = min_i(n, 15); pp = p2->u.pv.p;
		memcpy(u.iv15.v, p, 8*n2);
		for (int i=0; i<n2; i++) pp[i] = u.iv15.v+i;
		if (n2==n) return; 
		memcpy(p2->u.pv.v, p+15, 8*n-120);
		for (int i=15; i<n; i++) pp[i] = p2->u.pv.v + (i-15);
	} else {
		MxItem *p3 = mx_alloc('j'), *p2 = u.iv15.p2 = mx_alloc('i', p3->m_id);
		pp = p2->u.p; memcpy(p3->u.v, p, 8*n);
		for (int i=0; i<n; i++) pp[i] = p3->u.v + i;
	}}

void MxInP::del2() {
	if (pp==u.p8v8.p) return;
	MxItem *p2 = u.iv15.p2; if (!p2) return;
	int i3 = p2->m_pv; if (i3) mx_free(mx_ptr(i3));
	mx_free(p2); u.iv15.p2 = 0;
}

int MxInP::fi_cmp(const double *p, int n) {
	if (n<15 || n>22) return approx_cmp_v(p, pp[1], n);
	int r = approx_cmp_v(p, u.iv15.v+1, 14);
	return r ? r : approx_cmp_v(p+14, pp[15], n-14);
}

static int mx_chkv(double *p, double *q, int * pt, int tlim, double lim, int n) {
	int pos = 0, t3 = *pt;
	while(1) {
		int k0 = min_i(pos + tlim - t3, n), k = k0;
		if (q) while (k>pos && fabs(p[k-1])<lim && fabs(q[k-1])<lim) --k;
		else   while (k>pos && fabs(p[k-1])<lim) --k;
		if (k==pos) return (((*pt=t3+k0-k)>=tlim)<<15) + k0;
		if (k0==n) return *pt = k0-k, k0;
		t3 = k0-k; pos = k0;
	}}

static void mx_bdel_z (MxItem * p) { if (p->u.b.trec) trec_fin_1(p->u.b.trec);
	p->u.b.in.del2(); delete(p->u.b.bx); mx_free(p); }
static void mx_bdel_rf(MxItem * p) { int pvi = p->m_pv, nxi = p->m_nx;
	mx_ptr(pvi)->m_nx = nxi;  mx_ptr(nxi)->m_pv = pvi; mx_bdel_z(p); }
static void mx_bdel_c (MxItem * p) { if (p->u.b.ctp) p->b_ccut();  mx_bdel_z (p); }
static void mx_bdel   (MxItem * p) { if (p->u.b.ctp) p->b_ccut();  mx_bdel_rf(p); }

static void mx_fdel(MxItem * p) { p->f_reset(); p->u.f.in.del2(); mx_free(p); }

MxTrg::MxTrg(int n, int f, int flt, double *tmp, double *o0, double *o1) : siz(n), fdly(0) {
	double *q = tmp + 2*n;  pt[0] = tmp;  pt[1] = tmp+n;
	if (f>1)     pp5[0] = o0, pp5[1] =  o1, pp5[2] =  q, q +=   n, lrc01 = 1280;
	else if (o1) pp5[0] =  q, pp5[1] =  o1, pp5[2] = o0, q +=   n, lrc01 = 512*f+256;
	else         pp5[0] =  q, pp5[1] = q+n, pp5[2] = o0, q += 2*n, lrc01 = 512*f;
	if (flt) pp5[3] = q, pp5[4] = q+n, pp5[5] = 0;
	else     pp5[3]   =  pp5[4]    =   pp5[5] = 0;
}

int MxTrg::ret() {
	int md = (lrc01>>8)&7, n = siz, lrc;
	double *p0, *p1, *q;
	switch(lrc = lrc01&7) {
		case 0: return md<2 ? 0 : (md>4 ? (memset(pp5[0],0,8*n), memset(pp5[1],0,8*n), 2)
						: (memset(pp5[2],0,8*n), 1));
		case 1: if (!(md&1)) goto halfz;
			memset(pp5[1],0,8*n); if (md<4) memcpy(pp5[2], pp5[0], 8*n);  return 2;
		case 2: if (md&1) return memset(pp5[md&4?0:2],0,8*n), 2; else goto halfz;
		case 3: if (!(md&1))  { q = pp5[2], p0 = pp5[0], p1=pp5[1];
					for (int i=0;i<n;i++) q[i]=.5*(p0[i]+p1[i]); return 1; }
			if (md<4) memcpy(pp5[2],pp5[0],8*n);   return 2;
		case 4: return md>4 ? (memcpy(pp5[0], pp5[2], 8*n), memcpy(pp5[1], pp5[2], 8*n), 2) : 1;
		case 5: if (!(md&1)) goto halfp; if (md>4) goto m5_2;
			q = pp5[2]; p0 = pp5[0]; memcpy(pp5[1], q, 8*n); goto qp2;
		case 6: if (!(md&1)) goto halfp; if (md>4) goto m5_2;
			q = pp5[1]; p0 = pp5[2]; goto qp2;
		case 7: q = pp5[2]; p0 = pp5[0]; p1=pp5[1];
			if (!(md&1)) { for (int i=0;i<n;i++) q[i] += .5*(p0[i]+p1[i]); return 1; }
			if (md<4) for (int i=0;i<n;i++) p1[i] += q[i], q[i] += p0[i];
			else	  for (int i=0;i<n;i++) p1[i] += q[i], p0[i] += q[i];  return 2;
	}
qp2: 	for (int i=0;i<n;i++) q[i] += p0[i]; return 2;
halfz:  q = pp5[2]; p0 = pp5[lrc-1]; for (int i=0;i<n;i++) q[i]  = .5*p0[i]; return 1;
halfp:  q = pp5[2]; p0 = pp5[lrc-5]; for (int i=0;i<n;i++) q[i] += .5*p0[i]; return 1;
m5_2:	p0 = pp5[2]; q = pp5[lrc-5]; for (int i=0;i<n;i++) q[i] +=    p0[i]; return 2;
}

#define MXSW(X0,X1) k  = ((lrc01&m0) ? 0 : (prep(q0, dly, ae), lrc01|=m0, 1)); q0 += dly;        \
		    k += ((lrc01&m1) ? 0 : (prep(q1, dly, ae), lrc01|=m1, 2)); q1 += dly;         \
		    switch(k) { case 0: for (int i=0;i<n;i++) q0[i] += (X0), q1[i] += (X1); break; \
				case 1: for (int i=0;i<n;i++) q0[i]  = (X0), q1[i] += (X1); break; \
				case 2: for (int i=0;i<n;i++) q0[i] += (X0), q1[i]  = (X1); break; \
				case 3: for (int i=0;i<n;i++) q0[i]  = (X0), q1[i]  = (X1); break; }
int MxTrg::addb(int flg, int dly, int n, const double * v2) {
	static const char tab[32]={0,040,0,0, 0,012,0,0, 0,040,030,070, 0,010,011,020,
	                           0,043,0,0, 0,013,0,0, 0,043,033,073, 0,013,014,023};
	int k, m1, ae = dly+n, tf = tab[flg&31]+(flg&32), j = tf&7, m0 = 1<<j, op = tf>>3;
	double v0, v1, *p0, *p1, *q1, *q0 = pp5[j]; 
	if (dly<fdly) fdly = dly;
	MXLOG("addb: flg=0x%x, op=%d, j=%d, dly=%d, n=%d", flg, op, j, dly, n);
	switch (op&7) { case 4: q1=pp5[j+1]; m1=2*m0; goto c12;    case 1:  goto c11;
			case 3: q1=pp5[j+1]; m1=2*m0; goto c12x;   case 5:  goto c11u;
			case 2: q1=pp5[j+1]; m1=2*m0; goto c22;    case 6:  goto c22u;
			case 7: q1=pp5[j+1]; m1=2*m0; goto xmix;   default: goto oops; }
c11u:   p0=pt[0]; k = m0 & lrc01;
	if (!k) { prep(q0, dly, ae); lrc01 |= m0; memcpy(q0+dly, p0, 8*n); }
	else    { q0 += dly; for(int i=0; i<n; i++) q0[i] += p0[i]; }
	return flg&32768;
c22u:   q1=pp5[j+1]; m1=2*m0; p0 = pt[0]; p1 = pt[1];
	if (lrc01&m0) { q0+=dly; for(int i=0;i<n;i++) q0[i] += p0[i]; } 
	else { lrc01 |= m0; prep(q0,dly,ae); memcpy(q0+dly, p0, 8*n); }
	if (lrc01&m1) { q1+=dly; for(int i=0;i<n;i++) q1[i] += p1[i]; } 
	else { lrc01 |= m1; prep(q1,dly,ae); memcpy(q1+dly, p1, 8*n); }
	return flg&32768;
c11:	k = (~flg&1); p0 = pt[k]; v0 = v2[k]; k = m0 & lrc01; 
	if (!k) { prep(q0, dly, ae); lrc01 |= m0; q0 += dly; for (int i=0; i<n; i++) q0[i]  = v0*p0[i]; }
	else	{                    	          q0 += dly; for (int i=0; i<n; i++) q0[i] += v0*p0[i]; }
	return flg&32768;
c12x:	v0 = -v2[1]; v1 = v2[0]; p0 = pt[1]; goto c12y;
c12: 	v0 =  v2[0]; v1 = v2[1]; p0 = pt[0];
c12y:	MXSW(v0*p0[i], v1*p0[i]) return flg&32768;
c22:	v0 =  v2[0]; v1 = v2[1]; p0 = pt[0]; p1 = pt[1];
	MXSW(v0*p0[i], v1*p0[i]) return flg&32768;
xmix:   v0 =  v2[0]; v1 = v2[1]; p0 = pt[0]; p1 = pt[1];
	MXSW(v0*p0[i]-v1*p1[i], v1*p0[i]+v0*p1[i]) return flg&32768;
oops:   return log("mx/addb: BUG: flg=0x%x, n=%d", flg, n), 1;
}

/////// box /////////////////////////////////////////////////////////

void MxItem::b_ini(const char * updnnino, const double * arg, int delay, int tlim) {
	/*u.b.flg = (updnnino[3])>1;*/ int t3d = sec2samp(arg[4]);
	u.b.t1 = sec2samp(arg[2]); u.b.t2 = min_i(sec2samp(arg[3]), tlim);
	     u.b.up = TFunTab::ini(u.b.upx, updnnino[0], u.b.t1, 0);
	if ((u.b.dn = TFunTab::ini(u.b.dnx, updnnino[1], t3d   , 1))) u.b.t3 = u.b.t2 + t3d;
	else u.b.t3 = 0, u.b.dnx[0] = arg[4];
	u.b.in.ini(arg+6, updnnino[2]);
	u.b.delay = delay; u.b.ctp = u.b.t = 0; u.b.trec = 0;
}

int MxItem::b_evp(double *p, double *q, int n) {
	int t2, t = u.b.t - n;
	if (u.b.up && (t2 = u.b.t1-t) > 0) TFunTab::tf12(u.b.up, u.b.upx, p, q, min_i(t2, n));
	if (!u.b.dn) return mx_chkv(p,q, &u.b.t3, u.b.t2, u.b.dnx[0], n);
	if (t >= u.b.t2) TFunTab::tf12(u.b.dn, u.b.dnx, p, q, n);
	else if ((t2 = n+t-u.b.t2) > 0) TFunTab::tf12(u.b.dn, u.b.dnx, p+n-t2, q?q+n-t2:0, t2);
	return n;
}

int MxItem::b_calc(MxTrg * mtg, int n) {
	if (n<=0) return 0;
	int dly = u.b.delay; if (dly) { if(dly>=n) return u.b.delay-=n, 0; else u.b.delay=0, n -= dly; }
	int rv = (u.b.dn && u.b.t+n>=u.b.t3) ? (n = u.b.t3-u.b.t, 32768) : 0;
	if (n<=0) return QLOGV("mxstop1"), rv; else u.b.t += n;
	int of = u.b.bx -> calc_nzo2(u.b.ocfg, mtg->pt[0], mtg->pt[1], 0, u.b.in.pp, n);
	if (!m_id) bug("mx/b: wtf3");
	if (of<=0) return of ? (log("m/f#%x: removing b#%x: \"%s\"(%d)", u.b.up, m_id, err_str(of), of),
				gui_errq_add(of), 1)
			     : (!u.b.dn && (u.b.t3+=n) >= u.b.t2 && (QLOGV("mxstop3"),1) );
	int k = b_evp(mtg->pt[of==2], (of==3)?mtg->pt[1]:0, n); n = k&8191; rv |= (k&32768);
	if (rv && (debug_flags & DFLG_MX)) log("mxstop2 %d", rv);
	return mtg->addb(u.b.abflg|rv|of, dly, n, u.b.v);
}

void MxItem::b_ccut() {
	int pvi = u.b.ctp, nxi = u.b.ctn; u.b.ctp = 0;
	MxItem * q = mx_ptr(pvi); 
	if (q->m_ty == 'b') return (void)((q->u.b.ctn=nxi) && (mx_ptr(nxi)->u.b.ctp=pvi));
	if (q->m_ty != 'c') return log("BUG: b_ccut: invalid pv id:%d ty:0x%x", pvi, q->m_ty);
	if (nxi) return (void)(q->u.c.bx[q->c_find(u.b.ctk, 4)] = nxi, mx_ptr(nxi)->u.b.ctp = pvi);
	int k = u.b.ctk, f = 4 + 64*!q->c_find(k, 6);
	if (q->u.c.uk==k) q->u.c.uk = 0;
	BoxGen *bx = q->u.c.bref; 
	if (f&64) mx_free(q); if (bx) box_mxc_notify(bx, k, f);
}

void MxItem::b_fbv(MxItem * up, int ocf, double lr0, double v0) {
	u.b.ocfg = ocf;  int abf, fab, md, ty = up->m_ty;
	if (ty=='r') { if (~ocf&0x7c00) { abf = 8; u.b.v[1] = v0; goto one; } 
		       if (abf=0, fabs(lr0)>1e-14) { goto cs; } else { u.b.v[1] = 0.0; goto one; }}
	if (ty!='f') bug("mx/b_fbv: this=%p id=0x%x up(%p):id=0x%x ty=0x%x", this,m_id, up,up->m_id,ty);
	fab = up->u.f.abflg; md = !!(~ocf&0x7c00); abf = 16+8*md; md += 2*fab;
	MXLOG("fbv: md=%o, abf=0x%x  fab=0x%x", md&7, abf, fab);
	switch(md&7) {
		case 0: up->f_ini_1(lr0); u.b.v[1] = 0.0; goto one;
		case 1: up->f_ini_2();    u.b.v[1] =  v0; goto one;
		case 2: if (fabs(lr0-=up->u.f.lr)<1e-14) { u.b.v[1] = 0.0; goto one; }
			up->f_1to2(); goto cs2;
		case 3: up->f_1to2();
		case 7: lr0 = -1.0-up->u.f.lr; v0 *= (.5*M_SQRT2); goto cs2;
		case 4: goto cs;
		case 5: u.b.v[1] = v0; goto one;
		case 6: if (fabs(lr0-=up->u.f.lr)<1e-14) { u.b.v[1] = 0.0; goto one; }  goto cs2;
	}
cs:	v0 *= M_SQRT2; lr0 += 1.0; 
cs2:    lr0 *= (.25*M_PI); u.b.v[0] = v0*cos(lr0); u.b.v[1] = v0*sin(lr0); u.b.abflg = abf; return;
one:	u.b.abflg = abf | ((fabs((u.b.v[0]=v0)-1.0) < 1e-14) ? 36 : 4);   return;
}

/////// filter //////////////////////////////////////////////////////
// flg: 1-force arg(f): tlim vlim i_1 ... i_ni-1
//arg(b): v0 v1 u h d lr i_0 ... i_ni-1
MxItem * MxItem::r_mfind(BoxModel * mdl, int force) {
	int lo = 0, hi = m_x8 - 1, md, k;
	MxItem *r;
	long ref = (long)mdl, cur;
	while (lo<=hi) {
		cur = (long)(r = mx_ptr(u.r.mi[md = (lo+hi)>>1]))->u.m.m;
		if (cur==ref) return r;
		if (cur<ref) hi=md-1; else lo=md+1;
	}
	if (!force || m_x8>=120) return 0;
	r = mx_alloc('m', m_id); r->m_x8 = 0; r->u.m.m = mdl;
	if ((k=m_x8-lo)) memmove(u.r.mi+lo+1, u.r.mi+lo, 2*k);
	u.r.mi[lo] = r->m_id; ++m_x8;
	return r;
}

MxItem * MxItem::m_ffind(int ni, double *arg, int force, int up) {
	MxItem *r; int lo, hi, md, k;
	for (lo=0, hi=m_x8-1; lo<=hi; (k<0) ? (hi=md-1) : (lo=md+1))
		if (!(k = (r=mx_ptr(u.m.fi[md = (lo+hi)>>1]))->u.f.in.fi_cmp(arg+2, ni-1)))
			return force && ((arg[1]<r->u.f.vlim) && (r->u.f.vlim = arg[1]),
					 k = sec2samp(arg[0]),
					 (k > r->u.f.tlim) && (r->u.f.tlim = k)), r;
	if (!force || m_x8>=120) return 0;
	r = mx_alloc('f', m_id); r->u.f.abflg = 0; r->u.f.bx[0] = r->u.f.bx[1] = 0;
	if (debug_flags & DFLG_MX) { 
		log_n("m_ffind/new:"); for (int i=0; i<ni+1; i++) log_n(" %.15g", arg[i]); log(""); }
	r->u.f.tlim = sec2samp(arg[0]); r->u.f.vlim = arg[1]; r->u.f.up = up;
	r->u.f.in.ini(arg+1, ni); r->m_pv = r->m_nx = r->m_id; r->u.f.mdl = u.m.m;
	if ((k=m_x8-lo)) memmove(u.m.fi+lo+1, u.m.fi+lo, 2*k);
	u.m.fi[lo] = r->m_id; ++m_x8;
	return r;
}

void MxItem::rf_bfor(MxTrg * mtg, int n) {
	int i0 = m_id, icur = m_nx;
	while (icur != i0) {
		MxItem *p = mx_ptr(icur); icur = p->m_nx;
		if (p->b_calc(mtg, n)) mx_bdel(p); 
	}}

void MxItem::f_ini_1(double lr) {
	u.f.lr = lr; u.f.bx[0] = u.f.mdl->mk_box();
	if (fabs(lr)<1e-14) u.f.v[0] = u.f.v[1] = 1.0, u.f.abflg = 37;
	else lr+=1.0, lr*=.25*M_PI, u.f.v[0]=M_SQRT2*cos(lr), u.f.v[1]=M_SQRT2*sin(lr), u.f.abflg=1; }

void MxItem::f_ini_2() { u.f.bx[0] = u.f.mdl->mk_box();
	u.f.abflg = 46;  u.f.bx[1] = u.f.mdl->mk_box();  u.f.v[0] = u.f.v[1] = 1.0; }

int MxItem::f_calc(MxTrg * mtg, int n) {
	if (!m_id) bug("mx/f: wtf???");
	int ifg0, ifg1, dly, r2, chkf = (m_nx==m_id), oc = u.f.ocfg, abf = u.f.abflg&32764, f2 = abf&8;
	double *pi0, *pi1, **ppi=u.f.in.pp, *pt0 = mtg->pt[0], *pt1 = mtg->pt[1], zr = 0.0;
	if (chkf){ if (u.f.t >= u.f.tlim) { MXLOG("mx/f_calc: t(%d)>=tlim(%d)", u.f.t, u.f.tlim); return 1; }
		   if (!(abf&16384)) return u.f.t+=n, 0; else pi0 = pi1 = &zr, ifg0 = ifg1 = dly = 0; }
	else	 { u.f.t = 0; mtg->lrc01 &= ~24; mtg->fdly = ~abf & 16384; rf_bfor(mtg, n);
		   if ((n -= (dly=mtg->fdly)) <= 0) return 0; else u.f.abflg |= 16384;
		   pi0 = (ifg0 = (mtg->lrc01>>3)&1) ? mtg->pp5[3]+dly : &zr;
		   pi1 = (ifg1 = (mtg->lrc01>>4)&1) ? mtg->pp5[4]+dly : &zr; }
	int r0 =      (*ppi = pi0, u.f.bx[0]->calc_nzo2(oc, pt0, 0, ifg0, ppi, n)),
	    r1 = f2 ? (*ppi = pi1, u.f.bx[1]->calc_nzo2(oc, pt1, 0, ifg1, ppi, n)) : 0;
	if ((r2 = r0|(2*r1)) <= 0) return r2 ? (gui_errq_add(min_i(r0 , r1), "mx/f_calc"), 1)
					     : ((u.f.t+=n) >= u.f.tlim);
	if (chkf) chkf = r0 ? mx_chkv(pt0, r1?pt1:0, &u.f.t, u.f.tlim, u.f.vlim, n)
			    : mx_chkv(pt1,    0,     &u.f.t, u.f.tlim, u.f.vlim, n), n = chkf&8191;
	if (chkf &= 32768) MXLOG("mx/f_calc: done, n=%d, tlim=%d, vlim=%.15g", n, u.f.tlim, u.f.vlim);
	return mtg->addb(abf|chkf|r2, dly, n, u.f.v);
}

void MxItem::r2_clear() {
	for (int i=0, n=m_x8; i<n; i++) {
		MxItem * pm = mx_ptr(u.r.mi[i]);
		for (int j=0, m=pm->m_x8; j<m; j++) mx_fdel(mx_ptr(pm->u.m.fi[j])); 
		mx_free(pm); // TODO: model unref
	} m_x8 = 0; }

#define MXR_ZERO (f ? (memset(to0,0,8*n), (f>1) ? (memset(to1,0,8*n),2) : 1) : 0)
int MxItem::r_calc(double *to0, double *to1, int n, int f) {
	int fu = m_nx!=m_id, nfm, nfm0 = m_x8;  if (!(fu|nfm0)) return MXR_ZERO;
	double tmp[(3 + !to1 + 2*!!nfm0)*n]; 
	MxTrg mtg(n, f, nfm0, tmp, to0, to1);
	if (fu) { rf_bfor(&mtg, n); if (!nfm0) return mtg.ret(); }
	for (int i=nfm=0; i<nfm0; i++) {
		unsigned short fi, mi = u.r.mi[i];
		MxItem * pm = mx_ptr(mi);
		int nf = 0, nf0 = pm->m_x8;
		for (int j=0; j<nf0; j++) {
			MxItem * pf = mx_ptr(fi = pm->u.m.fi[j]);
			if (pf->m_ty!='f') bug("mx: r_calc: f expected");
			if (pf->f_calc(&mtg, n)) mx_fdel(pf);
			else if (nf!=j) pm->u.m.fi[nf++] = fi;
			else ++nf;
		}
		if (!(pm->m_x8 = nf)) mx_free(pm); // TODO: model unref
		else if (nfm!=i) u.r.mi[nfm++] = mi;
		else nfm++;
	}
	return m_x8 = nfm, mtg.ret();
}

void MxItem::rf_clear() { MxItem * p;
	for (int i0=m_id, i=m_nx; i!=i0; p=mx_ptr(i), i=p->m_nx, mx_bdel_c(p));
	m_nx = m_pv = m_id; }

/////// control /////////////////////////////////////////////////////

void MxItem::c_sane() {
	int i, j, k0 = -1, b0;
	for (i=0; i<u.c.nb; i++) {
		int b = u.c.bx[i];
		MxItem * q = mx_ptr(b); if (q->m_ty != 'b') bug("mx/c_sane: box exp 1");
		int k = q->u.b.ctk; if (k<=k0) bug("mx/c_sane: ooordah!"); else k0 = k;
		if (q->u.b.ctp != m_id) bug("mx/c_sane: prev0");
		j = 0; while (1) {
			if (q->u.b.ctn==0) break;
			b0 = q->m_id; q = mx_ptr(q->u.b.ctn); 
			if (q->m_ty != 'b') bug("mx/c_sane: box exp 2");
			if (q->u.b.ctp != b0) bug("mx/c_sane: prev1");
			if (q->u.b.ctk != k) bug("mx/c_sane: key");
			if (++j > 999) bug("mx/c_sane: loop");
		}}}

int MxItem::c_find(int k, int f) {
	int d, md, lo = 0, hi = u.c.nb - 1;
	while (lo <= hi) if (md = (lo+hi)>>1, !(d = k - mx_ptr(u.c.bx[md])->u.b.ctk)) goto found;
			 else (d<0) ? (hi=md-1) : (lo=md+1);
	if (f!=1) { if (f&4) bug("mx/c_find"); return MXE_CTLU; }
	if (u.c.nb>99) return MXE_CFULL;
	if ((d = u.c.nb - lo)) memmove(u.c.bx+lo+1, u.c.bx+lo, 2*d);
	++u.c.nb; return lo+256;
found:  return (f&2) ? (c_cut1j(md), u.c.nb) : md;
}

int MxItem::c_ins(MxItem * that, int k) {
	c_sane();
	int i, j, uf = k & 65536;
	uf && (i=u.c.uk, u.c.uk = (k&=65535), i) && (j = c_find(i, 4), c_stop_j(j, 1)) &&
		(c_cut1j(j), u.c.bref) && box_mxc_notify(u.c.bref, i, 4);
	if ((j = c_find(that->u.b.ctk = k, 1)) < 0) return j;
	that->u.b.ctp = m_id;
	if (j&256) that->u.b.ctn = 0, u.c.bx[j&255] = that->m_id;
	else	   i = that->u.b.ctn = u.c.bx[j], mx_ptr(i)->u.b.ctp = u.c.bx[j] = that->m_id;
	return c_sane(), u.c.bref ? box_mxc_notify(u.c.bref, k, 2 + (uf>>13)) : 0;
}

int MxItem::c_stop_j(int j, int f) {
	MxItem * bn;
	if (f==2) {
		j = u.c.bx[j]; do bn=mx_ptr(j), j=bn->u.b.ctn, mx_bdel_rf(bn); while(j);
		return 1;
	}
	unsigned short * pl = u.c.bx + j; 
	int td, jj, keepf = 1, lastid = m_id, i = *pl;
	while (1) {
		bn = mx_ptr(i);
		if (!bn->u.b.dn || !bn->u.b.t2) { if (f) goto kill; else goto keep; }
		if ((td = bn->u.b.t - bn->u.b.t2)<0) {
			if ((jj=bn->u.b.trec)) trec_fin_t(bn, jj, bn->u.b.t);
			if (bn->u.b.t2==bn->u.b.t3) goto kill;
			bn->u.b.t2+=td, bn->u.b.t3+=td; }
keep:		if (!keepf) keepf = 1, bn->u.b.ctp = lastid, *pl=i;
		lastid = i; pl = &bn->u.b.ctn; if ((i = *pl)) continue; else return *pl = 0;
kill:		keepf = 0; i = bn->u.b.ctn; mx_bdel_rf(bn);
		if (i) continue; else return (lastid==m_id) ? 1 : *pl = 0;
	}}

/////// 16-bit output ///////////////////////////////////////////////

static int mx_v8_min[8], mx_v8_diff[8], mx_v8_mul[8] = {
	556, 606, 661, 721, 786, 857, 935, 1020 };

static void mx_v8_init() { for (int i=0; i<8; i++) {
	double x = (double)(mx_v8_mul[i]), ymin = 128.0/x, ymax = (double)(131069*64)/x;
	memcpy(mx_v8_min+i,  (char*)&ymin + 4, 4);
	memcpy(mx_v8_diff+i, (char*)&ymax + 4, 4); mx_v8_diff[i] -= mx_v8_min[i];
}}

static void mx_d2s(short *q, au16w_t * cfg, int vd, const double *p, int n) {
        const char *s = (const char*)p + 4; n *= 8;
        int vol = cfg->vol + vd, vh = vol >> 3, vl = vol&7, shbase = 1050 - vh, mul = mx_v8_mul[vl],
            stp = cfg->nch, kmin = mx_v8_min[vl] - (vh<<20);
        unsigned int kdiff = (unsigned int) mx_v8_diff[vl];
        int i, k, sg, kd, v; for (i=0; i<n; i+=8, q+=stp)
                memcpy(&k, s+i, 4), sg = k>>31, kd = (k&=0x7fffffff)-kmin,
                *q = ((unsigned int)kd > kdiff) ? ((kd<0) ? 0 : (++cfg->ovf, (short)((32767^sg)-sg)))
                        : (v = (mul*(1048576|(k&1048575))) >> (shbase - (k>>20)),
                           (short)( (((v+1)>>1)^sg) - sg));
}

/////// live obj ////////////////////////////////////////////////////

static void mx_ldel(MxItem *p) {
	int pvi = p->m_pv, nxi = p->m_nx;
	if      (!pvi) { if (!nxi) mx_live0 = 0; else mx_ptr(mx_live0=nxi) -> m_pv = 0; }
	else           { mx_ptr(pvi)->m_nx = nxi; if (nxi) mx_ptr(nxi) -> m_pv = pvi; }
	mx_free(p);
}

static void mx_l_boxrm(MxItem *p) { p->u.l.bxi = 0;
	if (!(p->u.l.flg & MXLF_WIN)) mx_ldel(p); else gui_closewin(MX_L_WIN(p->m_id)); }

static int mx_l_closewin(MxItem *p) {
	int rv = p->m_nx; return p->u.l.bxi ? (void)(p->u.l.flg &= ~MXLF_WIN) : mx_ldel(p), rv; }

static void mx_l_debug() {
	MxItem *p; log_n("mixer live objs:");
	for (int i = mx_live0; i; i=p->m_nx) (p=mx_ptr(i))->debug();
	log(""); }

/////// export //////////////////////////////////////////////////////

int mx_mkroot() {
	MxItem * p = mx_alloc('r');
	return p ? (p->m_x8=0, p->m_pv=p->m_nx=p->m_id, p->u.r.flg=0, p->m_id) : MXE_GFULL; }

int mx_mkctl(BoxGen * bx) {
	MxItem * p = mx_alloc('c');
	return p ? (p->u.c.bref=bx, p->u.c.nb=p->u.c.uk=0, p->m_id) : MXE_GFULL; }

int mx_mklive(BoxInst * bxi) {
	MxItem * p = mx_alloc('l'); if (!p) return MXE_GFULL;
	if (mx_live0) mx_ptr(mx_live0)->m_pv = p->m_id;
	return p->m_pv=0, p->m_nx = mx_live0, p->u.l.flg = MXLF_WIN, p->u.l.bxi = bxi, mx_live0 = p->m_id; }

int mx_clear(int i) { 
	MxItem * p = mx_ptr(i);
	switch(p->m_ty) {
		case 'r': p->r2_clear(); p->rf_clear(); return 0;
		case 'f': p->f_reset(); return 0;
		default: return p->m_ty ? MXE_EXPRF : MXE_NOSUCH; }}

int mx_del(int i) {
	MxItem * p = mx_ptr(i);
	switch(p->m_ty) {
		case 'r': p->r2_clear(); p->rf_clear(); mx_free(p); return 0;
		case 'f': p->f_reset(); mx_free(p); return 0;
		case 'b': mx_bdel(p); return 0;
		default: return p->m_ty ? MXE_EXPBRF : MXE_NOSUCH; }}

int mx_r_isemp(int ix) { return mx_ptr(ix)->r_isemp(); }

// arg: tlim vlim x[1] ... x[ni-1]
int mx_add_filter(int trgi, BoxModel * mdl, int ni, double * arg, int osel) {
	if (debug_flags & DFLG_MX) {
		log_n("mx_add_filter(%d):",trgi); for (int i=0; i<ni+1; i++) log_n(" %.15g", arg[i]); }
	MXIARG(rnod, trgi, r);
	MxItem *mnod = rnod->r_mfind(mdl, 1); if (!mnod) return MXE_RFULL;
	MxItem *fnod = mnod->m_ffind(ni, arg, 1, trgi);
	MXLOG(" ...res: %d", fnod ? fnod->m_id : -123456789);
	return fnod ? (fnod->u.f.ocfg = osel, fnod->m_id) : MXE_MFULL;
}

// arg: v1 v2 u h d lr x[0] ... x[ni-1]
int mx_add_box(int trgi, BoxInst * bxi, const char * updnnino, const double * arg, int ocf,int delay,int tlim){
	if (trgi & ~65535) return MXE_WRONGID;
	MxItem *up = mx_ptr(trgi); 
	if (debug_flags & DFLG_MX) { log_n("mx_add_box(%d):",trgi); for (int i=0; i<updnnino[2]+6; i++)
		log_n(" %.15g", arg[i]); log(""); }
	MxItem *old0 = mx_ptr(up->m_nx), *rn = mx_alloc('b', trgi, old0->m_id);
	rn->u.b.bx = bxi; rn->u.b.up = trgi;
	up->m_nx = old0->m_pv = rn->m_id; rn->u.b.ocfg = ocf;
	rn->b_fbv(up, ocf, arg[5], arg[0]*arg[1]);
	rn->b_ini(updnnino, arg, delay, tlim);
	MXLOG("addbx(%d) res: %d, ocfg:0x%x, abf:0x%x", trgi, rn->m_id, rn->u.b.ocfg, rn->u.b.abflg);
	return rn->m_id;
}

int mx_c_stop(int ci, int ky, int f) {
	MXIARG(cn,ci,c); BoxGen * br = cn->u.c.bref; 
	if (ky>0) {
		int j = cn->c_find(ky&=65535,0); if (j<0) return j;
		if (!cn->c_stop_j(j, f)) return 0;
		cn->c_cut1j(j); int nf = cn->u.c.nb ? 4 : (mx_free(cn), 68);
		return br ? box_mxc_notify(br, ky, nf) : 0;
	}
	for (int n=cn->u.c.nb, i=0, j=0; i<n; i++) {
		int k = mx_ptr(cn->u.c.bx[i])->u.b.ctk;
		if (!cn->c_stop_j(i, f)) { cn->u.c.bx[j++] = cn->u.c.bx[i]; continue; }
		--cn->u.c.nb; if (br) box_mxc_notify(br, k, 4+64*!cn->u.c.nb);
		if (!cn->u.c.nb) return mx_free(cn), 0;
	}
	return 0;
}

int mx_c_dump_keys(char * to, int ci) {
	MxItem *c = mx_ptr(ci); int j = 0;
	if (c->m_ty!='c') return log("BUG: mx_c_dump_keys: id=%d, ty=0x%x", c->m_id, c->m_ty), 0;
	for (int i=0, n=c->u.c.nb; i<n; i++) {
		int k = mx_ptr(c->u.c.bx[i])->u.b.ctk, kh = k>>6, kl = k&63;
		if ((unsigned int)(--kh) < 51 && kl < 51) to[j]=kh+48, to[j+1]=kl+48, j+=2;
	}
	return j;
}

int mx_calc(int ix, double *to1, double *to2, int n, int f) {
	MXIARG(rn, ix, r); return rn->r_calc(to1, to2, n, f); }

int mx_au16_cfg(au16w_t * to, int nch, const char * s) {
	int e = (nch<1 && (nch=1)) || (nch>16 && (nch=16)), of = 0, ncp = 0;
	to->ovf = 0; to->nch = nch;
	for (int j, m, c, i=0; i<nch; i++) {
		if (!(c=s[i])) { of |= 128; break; }
		for (j=0; j<5; j++) if (c=="lracs"[j]) goto found;
		e |= (c!='z'); of |= 128; continue;
found:		if (of & (m=1<<j)) to->cp[ncp++] = (unsigned char)(16*i + to->trg[j]);
		else of |= m, to->trg[j] = i;
	}
	to->cp[ncp] = 0; to->oflg = of; return e ? MXE_CFGPAR : 0;
}

int mx_calc_int(int ix, short * to, au16w_t * cfg, fa_writer * fa, int n) {
	MXIARG(rn, ix, r);
	double tmp0[n]; unsigned char *s;
	int r, i, j, k, n2, vd = 0, ec = 0, oflg = 0, s_nch = 0, f_nch = fa ? fa->nch : 0;
	if (cfg && (s_nch = cfg->nch, n2 = s_nch*n, (oflg = cfg->oflg) & 128)) memset(to, 0, 2*n2);
	if ((oflg&3)||f_nch>1) {
		double tmp1[n], *p1;
		if ((r = rn->r_calc(tmp0, tmp1, n, 0)) <= 0) goto ze; else p1 = (r>1) ? tmp1 : tmp0;
		if (oflg&1) mx_d2s(to+(int)cfg->trg[0], cfg, vd, tmp0, n);
		if (oflg&2) mx_d2s(to+(int)cfg->trg[1], cfg, vd,   p1, n);
		if (f_nch!=1) { if (f_nch) ec = fa_add12(fa, tmp0, p1, 2*n); 
			        if (!(oflg&28)) goto cp; }
		if (r==1 && f_nch!=1) vd += 8; else for (i=0; i<n; i++) tmp0[i] += p1[i];
	} else  if ((r = rn->r_calc(tmp0,  0  , n, 0)) <= 0) { goto ze; }
	for (i=0; i<3; i++, vd-=4) if (oflg&(4<<i)) mx_d2s(to+(int)cfg->trg[i+2], cfg, vd, tmp0, n);
	if (f_nch==1) ec = fa_add12(fa, tmp0, 0, n);
cp:	if (cfg) for (s=cfg->cp; (k=*s); s++) for (j=k&15, k=(k>>4)-j; j<n2; j+=s_nch) to[j+k]=to[j];
	goto done;
ze:	if (r) return r;
	glob_flg |= GLF_SILENCE;
	if (fa) fa_add12(fa, zeroblkD, f_nch?zeroblkD:0, f_nch*n);
	if (oflg && !(oflg&128)) memset(to, 0, 2*n2);
done:	return ec<0 ? MXE_HCPFAIL : 0;
}

int mx_c_add(int ci, int bi, int ky) { MXIARG(cn,ci,c); MXIARG(bn,bi,b); return cn->c_ins(bn, ky); }
int mx_c_unlink(int ci) { MXIARG(cn,ci,c); cn->u.c.bref = 0; return 0; }

int mx_c_bpm_ugly_hack(int ci, int bp10m) {
 	MXIARG(cn, ci, c); int r = cn->c_find(1, 0);
	if (r<0) return log("mx_c_bpm_ugly_hack: box not found"), MXE_WTF;
	mx_ptr(cn->u.c.bx[r])->u.b.in.pp[0][0] = 0.1 * (double)bp10m; return 0;
}

unsigned char * mx_l_dat(int li) { MxItem * p = mx_ptr(li); return (p->m_ty == 'l') ? p->u.l.dat : 0; }

int mx_l_op(int li, int ix, int val) {
	if (li<0) { if (ix!=255 || val!=1) return MXE_L_ALL;
		    int i = mx_live0; while (i>0) i = mx_l_closewin(mx_ptr(i));
		    return i; }
	MXIARG(ln, li, l);
	if (ix<224) return (ix<0) ? MXE_L_IX : (ln->u.l.dat[ix] = val, 0);
	if (ix!=255) return MXE_L_IX;
	switch(val) { case 1:  return mx_l_closewin(ln), 0;
		      case 2:  return mx_l_boxrm(ln), 0;
		      case 3:  return ln->u.l.flg;
		      default: return MXE_L_FF; 
	}}

int mx_tr_add(int bi, BoxGen * bx) { MXIARG(bn, bi, b); int r = trec_alloc();
	return r<0 ? r : (trec_bx[r] = bx, trec_bi[r] = bi, bn->u.b.trec = r); }

void mx_tr_rm(int tri) { if (!trec_bx[tri]) log("BUG: mx_tr_rm(%d): nothing to remove", tri);
			 else mx_ptr(trec_bi[tri])->u.b.trec = 0, trec_free(tri); }

int mx_debug(const char *s) { switch(*s) {
	case 'L': return mx_l_debug(), 0;
	default:  return MXE_INVDBG; }}

void mx_init() {
	if (sizeof(MxItem)!=256) bug("mx/size");
	TFunTab::f_init();  mx_v8_init();
	for (int i=0; i<512; i++) mx_ptab[i] = (MxItem*)zeroblkC;
	for (int i=0; i<255; i++) trec_bi[i] = i+1;
	int r; if ((r=mx_mkroot()) != 0) bug("ini/mx_mkroot(): exp:0, got:%d", r);
	       if ((r=mx_mkroot()) != 1) bug("ini/mx_mkroot(): exp:1, got:%d", r);
}
