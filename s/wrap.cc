#include <limits.h>
#include <fcntl.h>
#include <sys/types.h>

#include "cmd.h"
#include "glob.h"
#include "gp.h"
#include "guistub.h"
#include "mx.h"
#include "util.h"
#include "util2.h"
#include "cfgtab.inc"

/////// decl ////////////////////////////////////////////////////////////////

#define WRF_NOCON  0x40000000 // bx 2m
#define WRF_PASS   0x40000000 // bx 
#define WRF_W_CH2  0x20000000 // bx
#define WRF_W_PLOT 0x10000000 // bx
#define WRF_S_NI   0x0000007c // bx
#define WRF_S_NO   0x00000f80 // bx
#define WRF_F_NI   0x0001f000 // bx
#define WRF_F_NO   0x003e0000 // bx
#define WRF_W_ALL  0x30000000
#define WRF_W_SH   28
#define WRF_SKIPV0 0x02000000 // 2m
#define WRF_MONO   0x01000000 // 2m
#define WRF_AVJOB  0x0300000f // 2m
#define WRF_SETOC  0x00800000 // wi
#define WRF_OC     0x00400000 // wi
#define WRF_CFILT  0x00010000 // wi
#define WRF_XFILT  0x00020000 // wi
#define WRF_XFILT_SH 17
#define WRF_V89    256        // v11
#define WRF_V10    512        // v11
#define WRF_V11A   768        // v11
#define WRF_GRMODE 3	      // bx
#define DBGC (debug_flags&DFLG_WRAP)

#define WR_AVCONF ((unsigned char*)(bxw_rawptr[1].c))
#define WR_TAB   (bxw_rawptr[2].c[6])
#define WR_WLG   (bxw_rawptr[2].c[7])
#define WR_SLFLG (bxw_rawptr[2].c[5])

int wrap_dump_keytab(char * to, unsigned int * bv, short ** pk);

// 1:lbl 2:src 4:fun 8:xf  32:c 64:v0 128:v1   256:(xf=='=') 512: slref 1024: unflg
class DWrapGen;
class WrapScale {
        public:
                static void st_init();
		int xop() const { return m_xop; }
                inline double v1(double x) { return (*m0_sfun[(int)m_fi])(this, x); }
                inline double v11(const double *p) { double x = (*m0_sfun[(int)m_fi])(this, p[(int)m_si]);
                        return (m_xop=='.') ? x : x * p[8+(m_xop>64)+(m_xop&1)]; }
                inline bool noncon() const { return !!m_fch; }
		inline int uses_xf() const { return m_xop != '.'; }
		inline int si() const { return m_si; }
                int cmd(const char *s);
		void draw(int oid, int ix, int flg, const double * pc, BoxGen * bx);
		void ini_default(int k);
		int save(SvArg * sv, int i);
		void wl_1(int flg), wl_2(int flg);
		void debug() { log_n("(%.15g:%.15g,%.15g,s%dx%cf%d,%c)", m_v0, m_v1, m_vx,
					m_si, m_xop?m_xop:'Z', m_fi, m_fch); }
        protected:
                typedef void (*ifun_p) (WrapScale*), (ifun_t) (WrapScale*);
                typedef double (*sfun_p) (WrapScale*,double), (sfun_t) (WrapScale*,double);
		typedef struct { ifun_p f; const char * s; } ifun_s;
                static ifun_t icon, ic1, ilin, ihrm, ilg, isq, irsq, icb, ircb;
                static sfun_t sc1, slin, slg, ssq0, ssq1, srsq0, srsq1, shrm0, shrm1,
                              ssq0z, ssq1z, scb0z, scb1z, scb0, scb1, srcb0, srcb1;
                TAB7_DECL(sfi, ifun_s);
                inline void updf() { (*sfi_ent(m_fch)->f.f) (this); }
                inline int sxop(int c) { if (c==63) return 0; if (c==m_xop) return 8; 
				         return  c=='='?(m_xop='.',296):(m_xop=c,40); }
                inline int sfun(int c) {return c==m_fch ? 0 : ((*sfi_ent(m_fch=c)->f.f)(this), 36);}
		inline const char * fnm() const { return sfi_ent(m_fch)->f.s; }
		inline const char * snm() const { return "x\0 y\0 s1\0s2\0s3\0s4\0s5\0s6" + 3*m_si; }
		inline int sch() const { return m_si + (m_si>1 ? 47 : 120); }
		inline int fch() const { return m_fch ? m_fch : 'c'; }
                int ssrc(int c);
                static sfun_p m0_sfun[17];
                double m_v0, m_v1, m_vx;
                char m_fch, m_fi, m_si, m_xop, rsrv[4];
};

class WrapScVec : public SOB {
	public: 
		SOBDEF_64(WrapScVec);
		WrapScVec(int arg) : SOB(arg) { memcpy(updn, "-A", 2);
				memset(qs, 0, 20);  pps[0] = pps[1] = pps[2] = 0; }
		WrapScVec(const WrapScVec * that, int arg);
		~WrapScVec() { clear(); }
		int save2(SvArg * sv);
		void clear();
		WrapScale * item(int ix, int f);
		void w_ud(int oid, int flg);
		void ini_default(int k);
		void debug2();
		int cmdi(int ix, const char * cmd);

		char updn[2], rsrv[2];
		unsigned char qs[20]; // rsss rsss
		WrapScale ** pps[3];
};

class WrapCore : public SOB {
	public:
		SOBDEF_64(WrapCore);
		typedef struct { short * pk[4]; unsigned int n5, bv[4]; } key_t;
		WrapCore(int uarg) : SOB(uarg) {}
		WrapCore(const WrapCore * that, int uarg) : SOB(uarg), xfd(that->xfd) { 
			if (DBGC) log("wrcore copy: %p -> %p", that, this);
			memcpy(tf01, that->tf01, 16); memcpy(grdim, that->grdim, 10); }
		~WrapCore() { del_keytab(); }
		void ini_default(int _) { tf01[0]=tf01[2]=0.0; tf01[1]=1.0; tf01[3]=22050.0;
			xfd=63; memcpy(grdim, "\031\031""333333\012\012", 10); } 
                int save2(SvArg * p);
		void debug2() { log("WrapCore: sorry"); } // TODO
		void grcmd(const char * i8) { gui2.c4('#', 'g', grdim[0]+48, grdim[1]+48);
			gui2.c4(grdim[8]+48, grdim[9]+48, i8[0]+48, i8[1]+48); }
		void upd_grid(const char * i8, int gui9) {
			if (gui9&1) gui2.setwin(gui9|11,'w'), gui2.w0(), grcmd(i8);
			if (gui9&8) gui2.setwin(gui9| 9,'#'), gui2.w0(), grcmd(i8); }
		void slcmd(const char * i8, int flg) { gui2.c1('!');
			for (int k=1,i=2; i<8; i++,k*=4)
				gui2.c2(flg&k ? i8[i]+48 : 44, flg&(2*k) ? grdim[i]+47 : 44); }
		int has_key(int i) { return ktab && (ktab->bv[(i>>5)&3] & (1u<<(i&31))); }
		int get_key(int i) { return ktab->pk[(i>>5)&3][i&31]; }
		void set_key_nz(int i, int v);
		void rm_key(int i);
		void set_key(int i, int v) { if (v) set_key_nz(i,v); else rm_key(i); }
		void del_keytab();
		int dump_keytab(char * to);
		float tf01[4];
		char grdim[10], xfd;
		key_t * ktab;
};
SOB_INIFUN(WrapCore, 1)

class WrapAutoVol : public SOB {
        public: 
		friend class WrapAVReader;
                static void mk0();
                static WrapAutoVol* mk(int uarg) { return new WrapAutoVol(uarg); }
                static void save_start(int sv_id);
                static bool save_end(FILE * f);
                static int cmd(char *s);
		static void del(WrapAutoVol * p) { delete(p); }
		static int save2a(SvArg *sv, ANode * nd, unsigned char * c8);
                WrapAutoVol(int uarg) : SOB(uarg), m_dat(0) { m_xy12rt[0] = 0;  }
                ~WrapAutoVol() { delete(m_dat); }
                int save2(SvArg * p);
		void debug2() { log("WrapAV: sorry"); } // TODO
                int upd_xy12rt(const unsigned char * p);
                void get_xy12rt(unsigned char * to) { memcpy(to, m_xy12rt, 8); }
                double vol(double x, double y, double z, double w);
                int vol_ix(int i, int j, int k, int l) { return l + m_xy12rt[3] * (k + m_xy12rt[2] * (j + m_xy12rt[1]*i)); }
                double vol_i(int i, int j, int k, int l) { return (double)m_dat[vol_ix(i,j,k,l)] / 46.0; }
                short * dat() { return m_dat; }
                int fill_data(char * s, int ver); // version: 0 1 2
                int pred(int i, int j, int k, int l);
                int pred2(int i, int j, int k, int l);
		unsigned char * xy12rt() { return m_xy12rt; }
		int parse_dim(const char * s);
		void set_c8(unsigned char * _) {}
		WrapAutoVol * copy(int uarg) { return new WrapAutoVol(uarg); } // no copy needed
        protected:
                static void unref2(WrapAutoVol* p, int ix);
                static int freeblk();
                static const int eol = -1234567890;

                short * m_dat;
                unsigned char m_xy12rt[8];
};

class WrapAVReader : public AReader {
	public:
		WrapAVReader(WrapAutoVol * av) : m_av(av) {}
		virtual int line(char * s);
	protected:
                LWArr<char> m_dat_str;
		WrapAutoVol * m_av;
};

class WrapAVJob : public Job {
        public: 
                WrapAVJob(JobQ::ent_t * ent, DWrapGen * wbx, WrapAutoVol * vtab, int mxid);
                virtual ~WrapAVJob() {}
                virtual int run1();
		virtual void abort();
        protected:
                DWrapGen * m_bx;
                WrapAutoVol * m_vt;
               	int m_mxid; 
                double m_fst, m_fstep;
                double m_max;
                unsigned char m_xy12ij[8], m_xy12ij_lim[8];
                double m_ixmul[4];
                short * m_trg;
};

class WrapSOB : public SOB {
	public:
		SOBDEF_64(WrapSOB);
		WrapSOB(int uarg) : SOB(uarg) {}
		WrapSOB(const WrapSOB * that, int uarg);
		void ini_default(int k);
                int save2(SvArg * p) { return p->st2=-1, 1; }
                void debug2();
		int icmd(const char * s, int flg, int gui9, const char * i8, BoxGen ** sf);
		void v(double *to, double * v11, int ix, int nf);
		void prep11(double *v11, int flg, const char * i8);
		WrapScale * sit_ro(int tix) { return m_scl[(tix>>6)&1].ro()->item(tix, 0); }
		WrapScale * sit_rw(int tix) { return SOB_RW(scl[(tix>>6)&1])->item(tix, 1); }
		int sit_cmd(int tix, const char * s) { if (DBGC) log("sit_cmd: 0x%x, \"%s\"", tix, s); 
						       return SOB_RW(scl[(tix>>6)&1]) -> cmdi(tix, s); }
		void unflg();
		void wl(int oid, int ix0, int n, int flg, const char *i8, BoxGen * bx);
		SOB_RW_F0(core, WrapCore) SOB_RW_F0(avol, WrapAutoVol)
		SOB_RW_F1(con, DblVec)    SOB_RW_F1(scl, WrapScVec)

		SOB_p<WrapCore> m_core;
		SOB_p<WrapAutoVol> m_avol;
		SOB_p<DblVec> m_con[2];
		SOB_p<WrapScVec> m_scl[2];
};

class AWrapGen : public BoxGen {
	public:
		friend void wrap_set_trec(BoxGen * bx, int v);
		AWrapGen(ABoxNode * nd);
		AWrapGen(ABoxNode * nd, const AWrapGen * that);
                virtual ~AWrapGen();
		virtual int add2mx_txdlv(int trg, int xflg, int dly, int lim, const double *vs) = 0;
		virtual const char * eff_rgb() = 0;
                virtual int n_in() const { return 0; }
                virtual int n_out() const { return 0; }
                virtual const char * cl_name() { return "wrap"; }
		virtual int  aux_window();
		virtual void box_window();
		virtual bool io_alias_perm() const { bug("wr: io_alias_perm() called"); return 0; }
		virtual int df_ui_ix() const { return 2; }
		virtual int mini(char * to);
		virtual void set_model() { bug("wr: set_model() called"); }
		virtual int mxc_notify(int k,int f){ if(k>63) w_gr_xyk(k&63,(k>>6)-1,f); if(f&64) m_mxctl=0;
						     return 0; }
		inline int pflg() { return m_bflg & WRF_PASS; }
		void delayed_clip_upd();
		int key_op(int k, int op, const char * xys, int nof);
	protected:
		typedef struct { int f; char xy[8]; } qsav_t;
		typedef int (acmd_t) (AWrapGen *, const char *, CmdBuf *);
		static int slf_conv(int k);
		static acmd_t c_c2k, c_cut, c_gcl, c_gmd, c_gr, c_ky, c_pl, c_stp, c_tf, c_wav, c_win;
		virtual WrapCore * core_ro() = 0; 
		virtual WrapCore * core_rw() = 0;
		virtual int get_nm2(char * to) = 0;
		virtual int show_tab_2(sthg * bxw_rawptr, int i) = 0;
		virtual int wlg(sthg * bxw_rawptr, int ix, int flg) = 0;
		virtual void w_col0(int flg) = 0;
		virtual int sl_bv() = 0;
		inline void qsav(qsav_t * q) { q->f = m_bflg; memcpy(q->xy, m_xys6, 8); }
		inline void qrst(qsav_t * q) { m_bflg = q->f; memcpy(m_xys6, q->xy, 8); }
		inline int qset(const char *s, int l) { for (int k,i=0; i<l; m_xys6[i++] = k) 
			                if ((unsigned int)(k=s[i]-48)>50u) return BXE_PARSE; return 0;}
		inline void wdat_cons_aw(sthg * bxw_rawptr) { WR_TAB=0; WR_WLG=0; WR_SLFLG=64; }
		inline int add2ctl(int mxbi, int mxky) {
			return mx_c_add(m_mxctl?m_mxctl:(m_mxctl=mx_mkctl(this)), mxbi, mxky); }
		int save2_aw(SvArg * sv);
		void gr_rgbc() { char buf[4]; gui2.c1('M'); get_nm2(buf); gui2.sn(buf, 2); 
			         gui2.sn(eff_rgb(), 6); }
		int plot_t(double t0, double t1, int n);
		int plot_f(double t0, double t1, double f0, double f1, int n, bool zpad);
		int mx1(int f=0) { int k, r = mx_mkroot(); if (r<0) return r;
			return (k=add2mx_txdlv(r,f,0,0,0))<0 ? (mx_del(r),k) : r; }
		int qcopy(int tf, int nof) { return tf ? trk_glob_paste(this, nof)
				: Node::mk(0, ClipNode::kcp(1), 0, 'W', nof|NOF_NOIDX, 0, this); }
		int batch_mono(double * to, int skip, int n);
		int write_a20();
		void w_gr_xyk(int x0, int y0, int f);
		void w_mini();
		void w_gmd() { gui2.setwin(w_oid(), 'w'); gui2.wupd('m'); gui2.c2('+', 48+(m_bflg&3)); }
		void w_a20(int flg);
		void w_tlim(int flg);
		void w_dim(int f20);
		void w_tab0(int f20);
		void w_slbv(int flg); // 1: force 2: setwin
		int show_tab(int i);
		char m_xys6[8];
		int m_bflg; 
		unsigned short m_mxctl; unsigned char m_trec, m_rsrv;
};

class DWrapGen : public AWrapGen {
	public:
		static void st_init() { cmd_init(); }
		DWrapGen(ABoxNode * nd);
		DWrapGen(ABoxNode * nd, const DWrapGen * that);
                virtual ~DWrapGen() {}
		virtual int add2mx_txdlv(int trg, int xflg, int dly, int lim, const double *vs);
		virtual int save2(SvArg * sv); 
		virtual void spec_debug();
		virtual void notify_nio(BoxGen * bx);
		virtual int save_sob(SvArg *p);
		virtual void wdat_cons(sthg * p);
		virtual int start_job_3(JobQ::ent_t * ent, char * arg);
		virtual const char * eff_rgb();
		WrapAutoVol * avol() { return m_sob.ro()->m_avol.ro(); }
		int avj_state() { return jobq.jst(m_node, 1); }
		int qdiff(DWrapGen * that); 
		int sob_from(int ix, BoxGen * bx0, int bxf);
		AReader * avreader() { return new WrapAVReader(SOB_RW(sob)->avol_rw()); }
        protected:
		virtual WrapCore * core_ro() { return m_sob.ro()->m_core.ro(); } 
		virtual WrapCore * core_rw() { return SOB_RW(sob)->core_rw(); }
		virtual int get_nm2(char * to);
		virtual int sl_bv();
		virtual int show_tab_2(sthg * bxw_rawptr, int i);
		virtual int wlg(sthg * bxw_rawptr, int ix, int flg);
		virtual void w_col0(int flg);
		void upd_conn() {Node::set_conn_2(m_node, BoxGen::node0(m_sfbx[0]),BoxGen::node0(m_sfbx[1]));}
		BXCMD_DECL(DWrapGen) c_in, c_bw, c_sf, c_ud, c_vt, c_so;
		int set_sf(int ff, BoxGen * bx);
		int set_sf_2(int ff, BoxGen * bx);
		void upd_updn(int flg);
		int grid_cmd(const char * s);
		void gr_rgbc(){ char buf[4]; get_nm2(buf); gui2.c3('M',buf[0],buf[1]); gui2.sn(eff_rgb(),6); }
		void w_sob(int ix, int sl);
		void w_avol(int f, unsigned char * s);
		int wupd1(int col, int ix1, int ix2 = 333);
		int av_guiconf(int c, const char * s);
		int avj_cmd(int stp);
		SOB_p<WrapSOB> m_sob;
		BoxGen * m_sfbx[2];
};

static int slbv_2(int r,int y) { return r>>=2, r | (64&(((0x0e13173f>>(4*(bitcnt_8(r)&6)))&63)-y)); }
/////// core SOB (grid config) //////////////////////////////////////////////

void WrapCore::del_keytab() { if (ktab) {
	for (int i=0; i<4; i++) if (ktab->bv[i]) ANode::f64(ktab->pk[i]);
	ANode::f64(ktab); ktab = 0; }}

void WrapCore::rm_key(int i) {
	if (!ktab) return; int ih = (i>>5) & 3, il = i & 31;
	if (!(ktab->bv[ih]&=~(1u<<il)) && (ANode::f64(ktab->pk[ih]), !--ktab->n5)) ANode::f64(ktab), ktab=0; }

void WrapCore::set_key_nz(int i, int v) {
	int ih = (i>>5) & 3, il = i & 31;
	if (!ktab) (ktab=(key_t*)ANode::z64())->pk[ih] = (short*)ANode::z64(), ktab->n5 = 1;
	else if (!ktab->bv[ih]) ktab->pk[ih] = (short*)ANode::z64(), ++ktab->n5;
	ktab->bv[ih] |= (1u << il);  ktab->pk[ih][il] = v;
}

int WrapCore::save2(SvArg * sv) {
	char buf[600]; memcpy(buf, "X$#*", 4);
	for (int i=0; i<10; i++) buf[i+4] = grdim[i]+48;   buf[14] = 10;
	int l = (xfd == 63) ? 15 : 15+sprintf(buf+15,"X$i%02xX\n", wr_ixtr_r(xfd));
	if (ktab) memcpy(buf+l, "X$k", 3), l += 4+wrap_dump_keytab(buf+l+3, ktab->bv, ktab->pk), buf[l-1] = 10;
	return sv->st2=-1, sv->out->sn(buf, l);
}

/////// scale ///////////////////////////////////////////////////////////////

static const double div1tab[52] = { 0.0, 0.0, 1.0, 0.5, 0.333333333333333, 0.25, 0.2,
0.166666666666667, 0.142857142857143, 0.125, 0.111111111111111, 0.1, 0.0909090909090909,
0.0833333333333333, 0.0769230769230769, 0.0714285714285714, 0.0666666666666667, 0.0625,
0.0588235294117647, 0.0555555555555556, 0.0526315789473684, 0.05, 0.0476190476190476,
0.0454545454545455, 0.0434782608695652, 0.0416666666666667, 0.04, 0.0384615384615385,
0.037037037037037, 0.0357142857142857, 0.0344827586206897, 0.0333333333333333,
0.032258064516129, 0.03125, 0.0303030303030303, 0.0294117647058824, 0.0285714285714286,
0.0277777777777778, 0.027027027027027, 0.0263157894736842, 0.0256410256410256, 0.025, 
0.024390243902439, 0.0238095238095238, 0.0232558139534884, 0.0227272727272727, 
0.0222222222222222, 0.0217391304347826, 0.0212765957446809, 0.0208333333333333,
0.0204081632653061, 0.02 };

static double dummy_v[8];
static double * dummy_pv[8];
static WrapScale dummy_sc[2];
static WrapScale * dummy_ps[8] = { dummy_sc, dummy_sc, dummy_sc, dummy_sc,
				   dummy_sc, dummy_sc, dummy_sc, dummy_sc };
static WrapScale ** dummy_pps[3] = { dummy_ps, dummy_ps, dummy_ps };

void WrapScale::st_init() { sfi_init(); dummy_sc[0].m_xop = dummy_sc[1].m_xop = '.'; }
static void wrap_fill(double * to, WrapScale *** ppps, double ** ppc, const unsigned char * pf,
                const double * p11, int ix, int nf, unsigned char * vfl = 0) {
        int i, j, k, f, f2, cf = -!(nf&WRF_NOCON), ix1 = nf & 63, f1v = !(nf & WRF_XFILT);
        if (!ix1) return; else ix1 += ix;
	unsigned char junk[8], *q; if (!vfl) vfl = junk; else memset(vfl, 0, 8);
        WrapScale **pps, *ps;
        double *pc;
i16:    if (!(pps = ppps[ix>>4])) pps = dummy_ps;
i8:     pc = (f = pf[k=ix>>3]) ? ppc[k] : dummy_v;
	q = vfl + k; f2 = f & cf;
i1:     i = ix & 7, j = ix & 15;
	*(to++) = (!(f2&(1<<i)) && (ps=pps[j>>1]) && (ps+=(j&1))->noncon()) ? 
		(*q|=((f1v||ps->uses_xf())<<i), ps->v11(p11)) : pc[i];
        if (++ix >= ix1) return; if(ix&7)goto i1; if(ix&8)goto i8; goto i16;
}

void WrapScale::icon(WrapScale *p) { p->m_fi = 99; p->m_fch = 0; }
void WrapScale::ic1 (WrapScale *p) { p->m_fi = 0; p->m_vx = p->m_v0; }
void WrapScale::ilin(WrapScale *p) { p->m_fi = 1; p->m_vx = p->m_v1 - p->m_v0; }
void WrapScale::ilg (WrapScale *p) { p->m_fi = 2; p->m_vx = log(p->m_v1 / p->m_v0); }
#define IFUN_ORD double a0, a1; \
        if ((a1=fabs(p->m_v1)) < (a0=fabs(p->m_v0))) { \
                if (a1<1e-280) { if (a0<1e-280) goto z2; else goto z1; } goto x1; \
        } else { \
                if (a0<1e-280) { if (a1<1e-280) goto z2; else goto z0; } goto x0; } \
z2:     p->m_fi = 0; p->m_vx = 0.0; return

void WrapScale::isq (WrapScale *p) {
        IFUN_ORD;
x0:     p->m_fi = 3; p->m_vx = sqrt(p->m_v1/p->m_v0)-1.0; return;
x1:     p->m_fi = 4; p->m_vx = sqrt(p->m_v0/p->m_v1)-1.0; return;
z0:     p->m_fi = 5; return;
z1:     p->m_fi = 6; return;
}

void WrapScale::icb (WrapScale *p) {
        IFUN_ORD;
x0:     p->m_fi = 7; p->m_vx = cbrt(p->m_v1/p->m_v0)-1.0; return;
x1:     p->m_fi = 8; p->m_vx = cbrt(p->m_v0/p->m_v1)-1.0; return;
z0:     p->m_fi = 9; return;
z1:     p->m_fi = 10; return;
}

void WrapScale::ihrm (WrapScale *p) {
        IFUN_ORD;
x0:     p->m_fi = 12; p->m_vx = p->m_v1/p->m_v0 - 1.0; return;
x1:     p->m_fi = 11; p->m_vx = p->m_v0/p->m_v1 - 1.0; return;
z0:;z1: goto z2;
}

void WrapScale::irsq (WrapScale *p) {
        IFUN_ORD;
x0:     p->m_fi = 14; p->m_vx = sqrt(p->m_v1/p->m_v0) - 1.0; return;
x1:     p->m_fi = 13; p->m_vx = sqrt(p->m_v0/p->m_v1) - 1.0; return;
z0:;z1: goto z2;
}

void WrapScale::ircb (WrapScale *p) {
        IFUN_ORD;
x0:     p->m_fi = 16; p->m_vx = cbrt(p->m_v1/p->m_v0) - 1.0; return;
x1:     p->m_fi = 15; p->m_vx = cbrt(p->m_v0/p->m_v1) - 1.0; return;
z0:;z1: goto z2;
}

double WrapScale::sc1 (WrapScale *p, double x) { return p->m_vx; }
double WrapScale::slin(WrapScale *p, double x) { return p->m_v0 + p->m_vx * x; }
double WrapScale::slg (WrapScale *p, double x) { return p->m_v0 * exp(p->m_vx * x); }
#define SCLFUN_2(NM, OP) \
double WrapScale::s##NM##0(WrapScale *p, double x) { return x = 1.0+p->m_vx*x, p->m_v0 OP; } \
double WrapScale::s##NM##1(WrapScale *p, double x) { return x = 1.0+p->m_vx*(1.0-x), p->m_v1 OP; }
#define SCLFUN_Z(NM, EX) \
double WrapScale::s##NM##0z(WrapScale *p, double x) { return p->m_v1 * (EX); } \
double WrapScale::s##NM##1z(WrapScale *p, double x) { return x=1.0-x, p->m_v0 * (EX); } 
SCLFUN_Z(sq,   x*x)   SCLFUN_Z(cb,   x*x*x)
SCLFUN_2(sq,  *x*x)   SCLFUN_2(cb,  *x*x*x)
SCLFUN_2(hrm, /x) SCLFUN_2(rsq, /(x*x)) SCLFUN_2(rcb, /(x*x*x))

void WrapScale::ini_default(int k) { switch(k) {
	case 0:  m_v0 = 0.00707596612706177; m_v1 = 1.0; m_si = 7; m_xop = 'A'; sfun('l'); return;
	case 5:  m_v0 = -1.0;                m_v1 = 1.0; m_si = 6; m_xop = '.'; sfun('-'); return;
	default: m_v0 =                      m_v1 = 0.0; m_si = 0; m_xop = '.'; sfun('\'');return; 
}}

int WrapScale::ssrc(int c) {
        if ((c = ((c+(c<64))&7)) == m_si) return 0;
	m_si = c; return 0x222; }

int WrapScale::cmd(const char *s) { switch(*(s++)) {
        case '<': m_v0 = at0f(s); updf(); return 32;
        case '>': m_v1 = at0f(s); updf(); return 32;
        case 's': return ssrc(*s);
        case 'F': return sxop(*s);
        case 'f': return sfun(*s);
	case '*': return s += parse_num(&m_v0, s), s += parse_num(&m_v1, s), 
                         s += parse_sep(s, ':'), ssrc(s[0]) | sfun(s[1]) | sxop(s[2]);
	case '?': return log("v0:%.15g v0:%.15g v0:%.15g fi:%d fch:'%c' sch:'%c'",
				  m_v0, m_v1, m_vx, m_fi, fch(), sch()), 0;
	default: log("wrscale: invalid command: \"%s\"", s-1); return BXE_CENUM;
}}

int WrapScale::save(SvArg * sv, int i) {
	return m_fch ? sv->out->pf("X$i%02x*%s:%s:%c%c%c\n", i, dbl2str_s(1,m_v0), 
				    dbl2str_s(2,m_v1), sch(), fch(), m_xop=='.'?'?':m_xop) : 1; }

void WrapScale::wl_1(int flg) { if (flg &  2) gui2.sz(snm()), gui2.c1(36);
				if (flg &  4) gui2.sz(fnm()), gui2.c1(36); }
void WrapScale::wl_2(int flg) { if (flg & 64) gui2.hdbl(m_v0);
				if (flg &128) gui2.hdbl(m_v1); }

/////// scl-vec /////////////////////////////////////////////////////////////

WrapScVec::WrapScVec(const WrapScVec * that, int arg) : SOB(arg) {
	if (DBGC) log("scvec copy: %p -> %p", that, this);
	memcpy(updn, that->updn, 2); memcpy(qs, that->qs, 20);
	for (int i=0; i<3; i++) { 
		if (!that->pps[i]) { pps[i] = 0; continue; }
		pps[i] = (WrapScale**) ANode::a64(); 
		for (int j=0; j<8; j++) pps[i][j] = (WrapScale*)ANode::cp64(that->pps[i][j]);
	}}

void WrapScVec::ini_default(int k) {
	if (!k) item(0,1)->ini_default(0), item(5,1)->ini_default(5), qs[0]=7, qs[2]=96; }

SOB_INIFUN(WrapScVec, 2)

void WrapScVec::clear() { for (int i=0; i<3; i++) { if (pps[i]) {
	for (int j=0; j<8; j++) if (pps[i][j]) ANode::f64(pps[i][j]), pps[i][j] = 0;
	ANode::f64(pps[i]), pps[i] = 0; }}}

int WrapScVec::save2(SvArg * sv) {
	int ec, r = 0, k = sv->st&1;
	if (!k) { if (updn[0] != '-') { CHKERR(sv->out->pf("X$U%c\n", updn[0])); }
		  if (updn[1] != 'A') { CHKERR(sv->out->pf("X$D%c\n", updn[1])); }}
	for (int i0=64*k, i=0; i<3; i++, i0+=16) {
		WrapScale ** pp = pps[i]; if (!pp) continue;
		for (int j=0; j<8; j++) { if (pp[j]) { CHKERR(pp[j][0].save(sv,wr_ixtr_r(i0+2*j)));
						       CHKERR(pp[j][1].save(sv,wr_ixtr_r(i0+2*j+1))); }}}
	return sv->st2=-1, r;
}

WrapScale * WrapScVec::item(int ix, int f) { 
	int i = (ix>>4) & 3, j = (ix&14)>>1;
	if (!pps[i]) { if (!f) return dummy_sc; pps[i] = (WrapScale**)ANode::z64(); goto mk; }
	if (pps[i][j]) return pps[i][j] + (ix&1); else if (!f) return dummy_sc;
mk:     return (pps[i][j] = (WrapScale*)ANode::cp64(dummy_sc)) + (ix&1);
}

int WrapScVec::cmdi(int ix, const char * s) {
	WrapScale * it = item(ix, 1);
	int i, j, r = it->cmd(s); 
	if (r&2) i=(ix>>1)&31, j=4*(ix&1), qs[i] &= ~(7<<j),
		 qs[i] |= ((it->si()&7)<<j);
	return r;
}

WrapScale::sfun_p WrapScale::m0_sfun[17] = {
        /* 0 */ &sc1, &slin, &slg, /* 3 */ &ssq0, &ssq1, &ssq0z, &ssq1z,
        /* 7 */ &scb0, &scb1, &scb0z, &scb1z, /* 11 */ &shrm0, &shrm1,
        /* 13 */ &srsq0, &srsq1, &srcb0, &srcb1 };

TAB7_DEF(WrapScale, sfi) { {'\'',{&icon, "con"}}, {'"', {&ic1,  "cn1"}}, {'-', {&ilin, "lin"}},
                           {'q', {&isq,  "sq" }}, {'c', {&icb,  "cub"}}, {'l', {&ilg,  "log"}},
			   {'Q', {&irsq, "/sq"}}, {'C', {&ircb, "/cu"}}, {'h', {&ihrm, "1/x"}},
			   {0, {0,0}} };

void WrapScVec::w_ud(int oid, int flg) { gui2.setwin(oid,'w');
	if (flg&1) gui2.wupd_c0('E', 't', 2), gui2.c2('/',  updn[0]);
	if (flg&2) gui2.wupd_c0('E', 't', 3), gui2.c2('\\', updn[1]); }

void WrapScVec::debug2() {
	log_n("WrapScVec: u:%d d:%d", updn[0],updn[1]);
	for (int i=0; i<3; i++) { if (pps[i]) {
		for (int j=0; j<8; j++) { WrapScale * q = pps[i][j]; if (q) {
			if (q[0].noncon()) log_n(" %d", 16*i+2*j  ), q[0].debug();
			if (q[1].noncon()) log_n(" %d", 16*i+2*j+1), q[1].debug(); }}}}
	log("");
}

/////// autovol /////////////////////////////////////////////////////////////

int WrapAutoVol::parse_dim(const char * s) {
	unsigned char t[6];
	for (int i=0; i<5; i++) t[i] = s[i] - 48;
	if (s[5]!=':') return VTE_PARSE;
	t[5] = atoi(s+6); return upd_xy12rt(t);
}

int WrapAutoVol::fill_data(char *s, int ver) {
        B91Reader tr; tr.init(s);
        int cnt = 0;
        int nx=m_xy12rt[0], ny=m_xy12rt[1], n1=m_xy12rt[2], n2=m_xy12rt[3];
        for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                        for (int k=0; k<n1; k++) {
                                for (int l=0; l<n2; l++) {
                                        if (vol_ix(i,j,k,l)!=cnt) return VTE_IXDIFF;
                                        int pre = pred2(i,j,k,l);
                                        m_dat[cnt++] = pre + tr.get_short_k(ver);
                                }}}}
        return 0;
}

int WrapAutoVol::pred(int i, int j, int k, int l) {
        if      (l) return m_dat[vol_ix(i,j,k,l-1)];
        else if (k) return m_dat[vol_ix(i,j,k-1,l)];
        else if (j) return m_dat[vol_ix(i,j-1,k,l)];
        else if (i) return m_dat[vol_ix(i-1,j,k,l)];
        else return 0;
}

int WrapAutoVol::pred2(int i, int j, int k, int l) {
        if (!(i|j|k|l)) return 0;
        int ix=vol_ix(i,j,k,l), ixn=0, ixd[4];
        if (i) ixd[ixn++] = vol_ix(1,0,0,0);
        if (j) ixd[ixn++] = vol_ix(0,1,0,0);
        if (k) ixd[ixn++] = vol_ix(0,0,1,0);
        if (l) ixd[ixn++] = vol_ix(0,0,0,1);
        if (ixn==1) return m_dat[ix - ixd[0]];
        int prlg_sum = 0, prlg_cnt = ixn * (ixn-1) / 2;
        for (int i=0; i<ixn; i++) 
                for (int j=i+1; j<ixn; j++) 
                        prlg_sum += m_dat[ix-ixd[i]] + m_dat[ix-ixd[j]] - m_dat[ix-ixd[i]-ixd[j]];
        return (prlg_sum + (prlg_cnt/2)) / prlg_cnt;
}
 
int WrapAutoVol::save2(SvArg * p) {
	Clock clk; int c1,c2,c3,c4; clk.reset();
	AOBuf * f = p->out;
        int total = 1; for (int i=0; i<4; i++) total *= m_xy12rt[i];
        short difftab[total];
        int ec, r = 1, cnt = 0;
        int nx=m_xy12rt[0], ny=m_xy12rt[1], n1=m_xy12rt[2], n2=m_xy12rt[3];
        for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                        for (int k=0; k<n1; k++) {
                                for (int l=0; l<n2; l++) {
                                        if (vol_ix(i,j,k,l)!=cnt) bug("voltab/save2: ix4(%d)!=cnt(%d)",vol_ix(i,j,k,l),cnt);
                                        int pre = pred2(i,j,k,l);
                                        int dif = (m_dat[cnt] - pre);
                                        difftab[cnt++] = dif;
                                }}}}
	c1 = clk.reset();
        if (debug_flags & DFLG_VOLTAB) log_sortedblk(difftab, total, 1, "difftab2:");
        int cost[3]; for (int i=0; i<3; i++) {
                cost[i] = 0;
                for (int j=0; j<total; j++) cost[i] += b91_cost(i,difftab[j]);
        }
	c2 = clk.reset();
        int ty = cost[1]<cost[2] ? cost[1]<cost[0] : 2*(cost[2]<cost[0]);
        if (debug_flags & DFLG_VOLTAB) log("voltab_cost(%d): %d %d %d --> %d", p->cn->id(), cost[0], cost[1], cost[2], ty);
        B91Writer wr;
        for (int i=0; i<total; i++) wr.put_short_k(ty, difftab[i]);
        int l = wr.n_bytes();
        char * s = wr.get_str();
	char buf[20]; memcpy(buf, "X$V:", 4);
	for (int i=0; i<5; i++) buf[4+i] = m_xy12rt[i]+48; buf[9] = ':';
	sprintf(buf+10, "%03d", m_xy12rt[5]); memcpy(buf+13, "\nN$V", 4);
	CHKERR(f->sn(buf, 17));
	c3 = clk.reset();
        while (l) {
		CHKERR(f->sn("\n<", 2));
                int k = (l<75) ? l : 75;
		CHKERR(f->sn(s, k)); s+=k; l-=k;
	}
	CHKERR(f->sn("!\n\"\n#\n"+2*ty, 2)); 
	c4 = clk.reset(); log("avol/save clk: %d %d %d %d", c1, c2, c3, c4);
	p->st2 = -1; return r;
}

int WrapAutoVol::upd_xy12rt(const unsigned char * p) {
        int k = 1; 
	for (int i=0; i<5; i++) {
		if (p[i]>9) return VTE_9;
		k *= (m_xy12rt[i] = p[i]);
	}
	m_xy12rt[5] = p[5];
        delete[](m_dat); m_dat = new short[k];
        for (int i=0; i<k; i++) m_dat[i] = -30000;
	return 0;
}

double WrapAutoVol::vol(double x, double y, double z, double w) {
        double xyzw[4];
        double tt[4];
        xyzw[0] = x; xyzw[1] = y; xyzw[2] = z; xyzw[3] = w;
        int ix[8]; // x0 y0 z0 w0 x1 y1 z1 w1
        for (int i=0; i<4; i++) {
                double ixd = xyzw[i]*((double)m_xy12rt[i]-1.0);
                if (fabs(ixd-round(ixd)) < 0.00001) {
                        ix[i] = ix[i+4] = (int)lround(ixd);
                        tt[i] = 0.5;
                } else {
                        ix[i] = (int)floor(ixd); if(ix[i]<0) ix[i]=0;
                        ix[i+4] = (int)ceil(ixd); if(ix[i+4]>=m_xy12rt[i]) ix[i+4]=m_xy12rt[i]-1;
                        tt[i] = ixd - (double)ix[i];
                }
        }
        double val[8];
        for (int i=0; i<8; i++)
                val[i] = (1.0-tt[3]) * vol_i(ix[0+(i&4)], ix[1+2*(i&2)], ix[2+4*(i&1)], ix[3])
                        + tt[3]      * vol_i(ix[0+(i&4)], ix[1+2*(i&2)], ix[2+4*(i&1)], ix[7]);
        for (int i=0; i<4; i++)
                val[i] = (1.0-tt[0]) * val[i] + tt[0] * val[i+4];
        for (int i=0; i<2; i++)
                val[i] = (1.0-tt[1]) * val[i] + tt[1] * val[i+2];
        return exp ((1.0-tt[2])*val[0] + tt[2]*val[1]);
}

WrapAVJob::WrapAVJob(JobQ::ent_t * ent, DWrapGen * wbx, WrapAutoVol * vtab, int mxid) :
	Job(ent), m_bx(wbx), m_vt(vtab), m_mxid(mxid), m_fst(0.0), m_max(0.0), m_trg(vtab->dat()) {
	memset(m_xy12ij, 0, 8); vtab->get_xy12rt(m_xy12ij_lim); 
	int t = 1; for (int i=0; i<6; i++) t *= m_xy12ij_lim[i];
        m_fstep = 1000.0 / (double) t;
	for (int k,i=0; i<4; i++) k = m_xy12ij_lim[i]-1, m_ixmul[i] = k ? 1.0 / (double)(k) : 0.0;
	m_trg = vtab->dat();
}

int WrapAVJob::run1() {
        int i, r, k = 5; while (k>3 && !m_xy12ij[k]) --k;
        double mv, buf[882];
        switch(k) {
                case 3:
			m_max = 0.0;
                case 4: 
			for (i=0; i<4; i++) buf[i] = m_ixmul[i] * (double)m_xy12ij[i];
			mx_clear(m_mxid);
			m_bx -> add2mx_txdlv(m_mxid, WRF_AVJOB|WRF_NOCON, 0, 0, buf);
                default:
			r = mx_calc(m_mxid, buf, buf+441, 441, 0);
                        mv = m_max; for (int i=0,n=r*441; i<n; i++) {
                                double v = fabs(buf[i]); if (v>mv) mv = v; }
                        m_max = mv;
        }
        int rv = (int)floor(m_fst += m_fstep); if ((unsigned int)rv>999u) rv = 999;
        if (++m_xy12ij[5]<m_xy12ij_lim[5]) return rv; else m_xy12ij[5] = 0;
	mx_clear(m_mxid);
        if (++m_xy12ij[4]<m_xy12ij_lim[4]) return rv; else m_xy12ij[4] = 0;
        if (m_max==0.0) {
                log("ERROR: av: max=0 at %d,%d,%d,%d", m_xy12ij[0], m_xy12ij[1], m_xy12ij[2], m_xy12ij[3]);
		gui2.errq_add(BXE_AVZERO); return JQE_FAIL; }
        int vi = (int)lround(-46.0 * log(m_max));
        if (vi>32767) vi = 32767; else if (vi<-32767) vi = -32767;
        *(m_trg++) = vi;
        for (int i=3; i>=0; i--)
                if (++m_xy12ij[i]<m_xy12ij_lim[i]) return rv; else m_xy12ij[i] = 0;
	log("av: m_fst = %.13g", m_fst);
        return JQE_DONE;
}

void WrapAVJob::abort() { if (m_mxid) mx_del(m_mxid); }

int WrapAVReader::line(char *s) {
	int l = strlen(s); if (!l) return 0;
	int k = 0; if (s[l-1]<36) k=1, --l;
	int n0 = m_dat_str.n();
	if (l) m_dat_str.resize(n0+l), memcpy(m_dat_str.p(n0), s, l);
	return k ? m_av->fill_data(m_dat_str.p(), s[l]-33) : 1;
}

/////// main shared obj /////////////////////////////////////////////////////

WrapSOB::WrapSOB(const WrapSOB * that, int uarg) : SOB(uarg) {
	if (DBGC) log("wrsob copy: %p -> %p", that, this);
	m_core.from(that->m_core); m_avol.from(that->m_avol); m_u8_refcnt |= (that->m_u8_refcnt&0x7f000000);
	for (int i=0; i<2; i++) m_con[i].from(that->m_con[i]), m_scl[i].from(that->m_scl[i]);
}

void WrapSOB::unflg() {
	for (int i=0; i<2; i++) {
		WrapScVec * scl = m_scl[i].ro();
		DblVec    * con = SOB_RW(con[i]);
		for (int j=0; j<40; j++) {
			WrapScale ** psi = scl->pps[j>>4], *sit;
			if (psi && (sit=psi[(j>>1)&7]) && sit[j&1].noncon() && 
				(con->bv[j>>3] & (1<<(j&7))) ) con->rm(j);
		}}}

void WrapSOB::v(double *to, double *v11, int ix, int nf) { int i = (ix>>6)&1;
	wrap_fill(to, m_scl[i].ro()->pps, m_con[i].ro()->p, m_con[i].ro()->bv, v11, ix&63, nf); }

void WrapSOB::prep11(double *v11, int flg, const char * i8) {
	WrapCore * core = m_core.ro();
	for (int i=0; i<8; i++) if (!(flg & (1<<i)))
		v11[i] = div1tab[(int)core->grdim[i]] * (double)i8[i];
	v11[10] = ((flg&512) && m_avol.ro()) ? m_avol.ro()->vol(v11[0], v11[1], v11[2], v11[3]) : 1.0;
	if (!(flg&256) || core->xfd==63) v11[8] = v11[9] = 1.0;
	else v(v11+8, v11, core->xfd, (flg&WRF_PASS)+1), v11[9] = 1.0 / v11[8];
}

void WrapSOB::wl(int oid, int ix0, int n, int flg, const char *i8, BoxGen *bx) {
	if (DBGC) log("wr/wl: ix0=%d n=%d flg=0x%x", ix0, n, flg);
	if (n<=0) return;
	int ix = wr_ixtr(ix0), i = (ix>>6) & 1, j = ix & 63;
	int skf = flg & WRF_CFILT;
	double v11[11], in[36];
	WrapScVec * scl = m_scl[i].ro();
	double ** pcon = 0;
	WrapCore * core = m_core.ro();
	unsigned char vfl[8];
	if (flg & 32) {
		int aflg = 0;
		for (int c,k=0; k<n; k++) if ((c=sit_ro(ix+k)->xop())!='.') aflg |= 1 + (c=='A');
		prep11(v11, (aflg<<8)|(flg&WRF_PASS), i8);
		pcon = m_con[i].ro()->p; if (!pcon) pcon = dummy_pv;
		wrap_fill(in, scl ? scl->pps : dummy_pps, pcon, m_con[i].ro()->bv, 
				v11, j, n|(flg&(WRF_PASS|WRF_XFILT)), vfl);
	}
	gui2.setwin(oid, 'w');
	for (int i2=ix,k=0; k<n; k++,i2++) {
		if (skf && !(vfl[(i2>>3)&7]&(1u<<(i2&7)))) continue;
		gui2.t0(); gui2.hexn(ix0=wr_ixtr_r(i2), 2); gui2.hexn(flg&255, 2);
		WrapScale * sit = sit_ro(i2);
		if (flg & 1 & ~(ix0>>5)) gui2.ionm(bx->node(), 0, ix0), gui2.c1(36);
		sit -> wl_1(flg);
		if (flg & 8) gui2.c2(i2==core->xfd ? '=' : sit->xop(), 36);
		if (flg & 32)  gui2.hdbl(in[k]);
		sit -> wl_2(flg);
	}
}

int WrapSOB::icmd(const char * s, int flg, int gui9, const char * i8, BoxGen ** sf) {
	WrapCore * core = m_core.ro(); 
	int ix0 = 16*(s[0]&7) + hxd2i(s[1]), ix = wr_ixtr(ix0), k,
	    i = (ix>>6)&1, j = ix&63, cf, r=2*(ix==core->xfd);
	if (j > 39) return BXE_IDX; else s += 2;
	if (*s=='X') return SOB_RW(core)->xfd = ix, 0; // save file only
	else if (*s=='=') r |= (flg & WRF_NOCON) ? (unflg(), flg &= ~WRF_NOCON, 1) : 0,
		     cf = (s[1]=='x' && !s[2]) ? (SOB_RW(con[i])->rm(j), s+=2, 4) 
				    	       : (parse_num(SOB_RW(con[i])->addp(j), s+1), 0);
	else         cf = sit_cmd(ix, s);
	if (cf<0) return cf;
	if (cf&8) {
		if (!(cf & 256)) (r&2) && (SOB_RW(core)->xfd = 63);
		else if ( ((r^=2)&2) && (k = core->xfd, SOB_RW(core)->xfd = ix, gui9&1)) 
			wl(gui9|11, wr_ixtr_r(k), 1, 8, 0, sf[i]);
	}
	if (gui9 & 1) wl(gui9|11, ix0, 1, cf+(flg&WRF_PASS), i8, sf[i]);
	return r | ((cf>>7)&4);
}

void WrapSOB::debug2() {
	log("WrapSOB: flg:0x%x core s_scl f_scl s_con f_con av", m_u8_refcnt>>24);
	m_core.debug(); m_scl[0].debug(); m_scl[1].debug();
			m_con[0].debug(); m_con[1].debug(); m_avol.debug(); }

void WrapSOB::ini_default(int k) {
	m_core.set(WrapCore_default(k));
	for (int i=0; i<2; i++) m_con[i].set(DblVec_default(i+1)),
			        m_scl[i].set(WrapScVec_default(i)); }
SOB_INIFUN(WrapSOB, 1)

/////// box (abs) ////////////////////////////////////////////////////////////

AWrapGen::AWrapGen(ABoxNode * nd) : BoxGen(nd), m_bflg(0), m_mxctl(0) {
	memcpy(m_xys6, "\0\0\0\0\0\0\031\062", 8);
	nd->set_ui_d(2); nd->setbx_0(this); }

AWrapGen::AWrapGen(ABoxNode * nd, const AWrapGen * p) : 
	BoxGen(nd), m_bflg(p->m_bflg), m_mxctl(0) {
	memcpy(m_xys6, p->m_xys6, 8); 
	nd->set_ui_f(p->node()); nd->setbx_0(this);
}

AWrapGen::~AWrapGen() { if (m_mxctl) mx_c_unlink(m_mxctl);
		        if (m_trec) mx_tr_rm(m_trec);  }

int AWrapGen::slf_conv(int k) {
	int v = 0, d = (k&64) ? (03221111 >> (3*bitcnt_8(k&63))) & 7 : 9;
	for (int i=0; i<6; i++) if (k&(1<<i)) v |= 1 << (2*i+(--d<0));    return v; }

int AWrapGen::aux_window() {
	gui2.cre(w_oid(9), '#'); gui2.c1('T'); gui2.nname(m_node);
	gui2.c1(36); gr_rgbc();
	WrapCore * cr = core_ro();
	cr->slcmd(m_xys6, 4095); w_slbv(1); cr->grcmd(m_xys6);
	if (cr->ktab) gui2.c1('B'), gui2.wr_keylist(cr->ktab->bv, cr->ktab->pk); 
	return 16*'#'+9;
}

int AWrapGen::save2_aw(SvArg * sv) {
	BXSV2_HEAD; int k, l = 12;
	char buf[24]; memcpy(buf, "X$G", 3); buf[3] = 51 + 3*!(m_bflg&WRF_NOCON);
	for (int i=0; i<8; i++) buf[i+4] = m_xys6[i]+48;
	if ((k=m_bflg&3)) memcpy(buf+l, "\nX$m", 4), buf[16] = 48 + (m_bflg&3), l = 17;
	buf[l] = 10; CHKERR(f->sn(buf, l+1));
	if ((k=m_node->etc()->i[0])) { CHKERR(sv->out->pf("X$T%x\n", k)); }
	return r;
}

void AWrapGen::w_tlim(int f) {
	int k = m_node->etc()->i[0], t = k&INT_MAX, mf = k<0;
	gui2.setwin(w_oid(), 'w'); 
	if (f&1) gui2.wupd_i1('Y', mf,  20);
	if (f&2) gui2.wupd_d('Y', (double)t/40320.0, 22);
}

int AWrapGen::show_tab(int i) {
	BXW_GET; if (i<0) i=WR_TAB; else if (i==WR_TAB) return 0; else WR_TAB=i;
	gui2.setwin(w_oid(), 'w'); gui2.wupd_0('Y', ".W"); gui2.c1(48+i);
	switch(i) {
		case 0: w_tab0(0xaaaaa); return 0;
		case 1: w_a20(-1); return 0;
		case 2: w_tlim(-1); return 0;
		default: return show_tab_2(bxw_rawptr, i);
	}}

void AWrapGen::box_window() {
	gui2.cre(w_oid(), 'w'); gui2.own_title(); show_tab(-1); w_mini(); w_gmd();
	WrapCore * cr = core_ro();
	const float * p = cr->tf01;
	for (int i=0; i<4; i++) gui2.wupd_d((0x46665474>>(8*i))&127, (double)p[i]);
	cr->upd_grid(m_xys6, w_oid(1));
	sthg * q = m_node->wdat_raw();
	int ec; for (int i=0; (ec=wlg(q, i, pflg()))>=0; i++) ;
	if (ec!=BXE_WLGDONE) gui_errq_add(ec, "wrap/wlg");
}

// rel b1 cp set   sb1 cp2 set* cb1   kill stop tggl ply  uniq rsrv rsrv nop
int AWrapGen::key_op(int k, int op, const char * xys, int nof) {
	static const unsigned int optr[4] = { 0xc65a32bf, 0xc65a32b9, 0xc65b32af, 0xb65a32cf };
	if (op<8) { if (op<0) return BXE_CENUM; else op = (optr[m_bflg&3]>>(4*op)) & 15; }
	if (op==15) return 0;
	char buf[8], xysav[8]; int ec, nb, bfsav = m_bflg;
	if (k<0) { 
		nb = 2+6*(k&1); for (int j,i=0;i<nb;i++) if ((j=xys[i]-48)<0) return BXE_PARSE; else buf[i]=j;
		k = 64*buf[1]+buf[0]+64;
	} else {
		if (k&65536) { WrapCore *cr = core_ro(); k &= 65535;
			       if (!cr->has_key(k)) return BXE_UNDEFKEY; else k = cr->get_key(k); }
		if (k<64) nb = 0; else nb = 2, buf[0] = k&63, buf[1] = (k>>6)-1;
	}
	if ((unsigned int)(op-8)<3u){ if (m_mxctl && (ec=mx_c_stop(m_mxctl,k,1+(op==8)))!=MXE_CTLU) return ec;
				      if (op==10) ++op; else return EEE_NOEFF; }
	if (op!=3 && op!=6) memcpy(xysav, m_xys6, 8);
	else if (!(nof&NOF_FORCE) && !m_node->perm(DF_EDBOX)) return NDE_PERM;
	if (nb) { if (op!=6) m_bflg |= WRF_NOCON; for (int i=0; i<nb; i++) m_xys6[i] = buf[i]; }
	switch(op) {
		case  2: ec = qcopy(0, nof); break;
		case  3: case 6: 
			 if (!nb) return EEE_NOEFF;
			 w_gr_xyk(m_xys6[0], m_xys6[1], 1); w_mini(); 
			 if (wnfl()) w_col0(pflg());
			 return 1;
		case  5: ec = qcopy(1, nof & ~NOF_FGUI); break;
		case 11: ec = add2mx_txdlv(0,0,0,0,0); if (ec>=0) ec = add2ctl(ec, k);       break;
		case 12: ec = add2mx_txdlv(0,0,0,0,0); if (ec>=0) ec = add2ctl(ec, k|65536); break;
		default: ec = BXE_CENUM; break;
	}
	if (nb) m_bflg = bfsav, memcpy(m_xys6, xysav, 8);
	return ec;
}

void AWrapGen::delayed_clip_upd() {
	if (!wnfl(8)) return; else m_node->winflg_and(~8);
	ANode * cl = m_node->up(); switch(cl->cl_id()) {
		case 'C': static_cast<ClipNode*>(cl)->show_newbox(m_node); return;
		case 't': trk_cond_pm(cl->box0(), m_node, '+'); return;
		default:  log("BUG: wr/delayed_clip_upd: up() is no clipb or trk"); return;
	}}

int AWrapGen::mini(char *to) {
	get_nm2(to); for (int i=0; i<4; i++) to[i+2] = i2aA(m_xys6[i]);
	to[6] = i2aA(m_xys6[7]); to[7] = ':';
	memcpy(to+8, eff_rgb(), 6); return 14;
}

void AWrapGen::w_a20(int flg) { int k; gui2.setwin(w_oid(),'w'); 
	if (flg&(k=WRF_W_CH2 )) gui2.wupd_i1('Y', !!(m_bflg&k),  10);
	if (flg&(k=WRF_W_PLOT)) gui2.wupd_i1('Y', !!(m_bflg&k),  11); }

void AWrapGen::w_mini() {
	if (wnfl()) gui2.setwin(w_oid(), 'w'), gui2.wupd_c0('w', 's'), gui2.bxmini(this);
	ANode * up = m_node->up(); 
	if (up->cl_id()=='C' && up->winflg(8)) static_cast<ClipNode*>(up)->draw_1(m_node); // TODO: trk
}

// f: 1-sel 2-on 4-off 8-uniq
void AWrapGen::w_gr_xyk(int x0, int y0, int f) {
	for (int g9 = m_node->gui9(), i=1, x=x0+48, y=y0+48; i<9; i+=7) { if (g9 & i) {
		int c = 0x91623ee >> i; gui2.setwin(g9|((c>>16)&15), c&127); gui2.wupd('#');
		if (f&1) gui2.c3('b', x, y);
		if (f&6) gui2.c3('K'+16*(f&2), x, y);
		if (f&8) gui2.c3('u', (f&4)?33:x, y);
	}}}

void AWrapGen::w_tab0(int f20) { gui2.setwin(w_oid(), 'w'); const char *s = core_ro()->grdim;
	for (int i=0,j=2; i<10; i++,j*=4) if (f20&j) gui2.wupd_c48('Y', s[i], i); }

// flg: 16:slupd 17:xry 18:slf 19:yry
void AWrapGen::w_dim(int f20) {
	int wf = wnfl(2560); if (!wf) return;
	WrapCore * cr = core_ro();
	int grf = f20 & 0xa000f;
	if (f20 & 0x14055) w_mini();
	if (wf & 512) {
		gui2.setwin(w_oid(9), '#'); gui2.t0();
		int slf = (f20>>4) & (f20&65536 ? 0xfff : 0xaaa);
		if (slf) cr->slcmd(m_xys6, slf);
		w_slbv(0);
		if (grf) cr->grcmd(m_xys6);
	}
	if (wf & 2048) {
		BXW_GETV("wr/dim"); gui2.setwin(w_oid(), 'w');
		if (!WR_TAB && (f20&0xaaaaa)) w_tab0(f20);
		if (grf) gui2.w0(), cr->grcmd(m_xys6);
		if (f20 & 0x5555) w_col0(pflg());
	}}

void AWrapGen::w_slbv(int flg) {
	if (!m_node->winflg(512)) return; BXW_GETV("slbv"); 
	int sf = sl_bv(); 
	if (sf!=WR_SLFLG) WR_SLFLG = sf; else if (!(flg&1)) return;
	if (flg&2) gui2.setwin(w_oid(9),'#'), gui2.t0();
	gui2.c1('+'); gui2.hexn(slf_conv(sf), 3);
}

int AWrapGen::write_a20() {
	int n, ec, ch2 = !!(m_bflg&WRF_W_CH2), mxr = mx1(ch2?0:WRF_MONO), nf;
	if (m_bflg&WRF_W_PLOT) {
		float * q = core_ro()->tf01;
		int skip = (int)lround((double)q[0] * (double)sample_rate);
		if  ((nf = (int)lround((double)q[1] * (double)sample_rate) - skip)<1) return BXE_ZEROWAV;
		while (skip>1023) mx_calc(mxr, junkbuf, ch2?junkbuf:0, 1024, 0), skip -= 1024;
		if    (skip)      mx_calc(mxr, junkbuf, ch2?junkbuf:0, skip, 0);
	} else { nf = sample_rate * CFG_AO_TLIM.i; }
	if (mx_r_isemp(mxr)) return BXE_ZEROWAV;
	fa_writer out; if ((ec=fa_start(&out, 1+ch2))<0) return ec;
	do n = min_i(nf, 1024), mx_calc_int(mxr, 0, 0, &out, n), nf -= n; while (nf && !mx_r_isemp(mxr));
	mx_del(mxr); gui_acv_op(out.id); return fa_end(&out)<0 ? EEE_ERRNO : 0;
}

int AWrapGen::batch_mono(double * to, int skip, int n) {
	Clock clk; clk.reset();
	int mxid = mx1(WRF_MONO); if (mxid<0) return mxid;
	while (skip>=1024) mx_calc(mxid, junkbuf, 0, 1024, 1), skip -= 1024;
	if (skip) mx_calc(mxid, junkbuf, 0, skip, 1);
	while (n>=1024) mx_calc(mxid, to, 0, 1024, 1), n -= 1024, to += 1024;
	if (n) mx_calc(mxid, to, 0, n, 1);
	mx_del(mxid); return log("batch_mono: %d", clk.get()), 0;
}

int AWrapGen::plot_t(double t0, double t1, int n) {
	log("plot_t: %.15g %.15g %d", t0, t1, n);
        int i0 = sec2samp(t0), i1 = sec2samp(t1), len = i1 - i0;
        if (len<1) return log("plot_t: length = %d samples, sorry.", len), BXE_RANGE;
	double * res = new double[len];
	int r = batch_mono(res, i0, len); if (r<0) return r;
        if (len <= n) {
                PlotPar_arr par(res, len, t0, t1);
                Gnuplot::sg()->setfun1(0, arrfun1, &par);
                Gnuplot::sg()->plot1(1, t0, t1, len);
        } else {
                double * stat = new double [ 3 * n ];
                samp_stat(res, len, n, false, 0.0,
                                stat, stat + n, stat + 2 * n);
                PlotPar_arr p_min(stat, n, t0, t1);
                PlotPar_arr p_max(stat+n, n, t0, t1);
                PlotPar_arr p_avg(stat+2*n, n, t0, t1);
                Gnuplot::sg()->setfun1(0, arrfun1, &p_min);
                Gnuplot::sg()->setfun1(1, arrfun1, &p_max);
                Gnuplot::sg()->setfun1(2, arrfun1, &p_avg);
                Gnuplot::sg()->plot1(7, t0, t1, n);
                delete[] (stat);
        }
        delete[] (res); return 0;
}

int AWrapGen::plot_f(double t0, double t1, double f0, double f1, int n, bool zpad) {
	log("plot_f: %.15g %.15g %.15g %.15g %d %c", t0, t1, f0, f1, n, "ny"[zpad]);
        int i0 = sec2samp(t0), i1 = sec2samp(t1), len = i1 - i0;
        if (len<1) return log("plot_f: length = %d samples, sorry.", len), BXE_RANGE;
        int siz = 64, bits = 6; 
        if (len > (1<<24)) {
                log("plot(F): len cut to %g", sample_length*(double)(1<<24));
                siz = len = (1<<24);
                bits = 24;
        } else {
                while (len > siz) siz<<=1, ++bits;
                if (!zpad) len = siz;
        }
        double * re = new double[ 2 * siz ];
        double * im = re + siz; 
	int r = batch_mono(re, i0, len); if (r<0) return r;
        for (int i=len; i<2*siz; i++) re[i] = 0.0;
        double * res = fft(re, im, bits, false);
        double * res2 = res + siz;
	int fmx = -1; double fmv = -1.0;
        for (int i=0; i<=siz/2; i++)
                if ((res[i] = sqrt(res[i]*res[i] + res2[i]*res2[i]))>fmv) fmv = res[i], fmx = i;
	log("plot_f: max @ %g Hz", (double)fmx / (double)siz * (double)sample_rate);
        int ix0 = (int)round(f0 * sample_length * (double)(siz));
        int ix1 = (int)round(f1 * sample_length * (double)(siz));
        if (ix1<siz-1) ++ix1;
        int n0 = ix1 - ix0;
        if (n > n0) n = n0;
	if (n<2 || n0<2) return log("plot_f: n=%d, n0=%d", n, n0), BXE_RANGE;
        bool statflg = (n < n0);
        double * stat = new double [ (statflg ? 3 : 2) * n ];
        samp_stat(res+ix0, n0, n, false, 0.0, 0, 0, stat);
        samp_stat(res+ix0, n0, n, true, 55.0, 0,
                        statflg ? stat+2*n : 0, stat+n);
        PlotPar_arr p_avg(stat, n, f0, f1);
        PlotPar_arr p_lavg(stat+n, n, f0, f1);
        PlotPar_arr p_lmax(stat+2*n, n, f0, f1);
        Gnuplot::sg()->setfun1(0, arrfun1, &p_avg);
        Gnuplot::sg()->setfun1(1, arrfun1, &p_lavg, true);
        if (statflg) Gnuplot::sg()->setfun1(2, arrfun1, &p_lmax, true);
        Gnuplot::sg()->plot1(statflg ? 7 : 3, f0, f1, n);

        delete[] (res); delete[] (stat); return 0;
}

#define CH(X) BXCMD_H(AWrapGen, X)

CH(gcl){return (s[1]=='.') ? (s[2]==51 ? p->m_node->draw_window(25) 
			               : p->key_op(-2, s[2]-48, s+3, cb->cnof()))
			   : p->key_op(-1, s[1]-48, s+2, cb->cnof()); }

CH(gmd){p->m_bflg&=~3, p->m_bflg|=(s[1]&3); if (p->wnfl()) p->w_gmd(); return 0; }
CH(c2k){return Node::mk(0, ClipNode::kcp(1), 0, 'W', cb->cnof()|NOF_NOIDX, 0, p); } // TODO
CH(stp){return p->m_mxctl ? mx_c_stop(p->m_mxctl, 0, s[1]&3) : EEE_NOEFF; }

CH(tf){	float *q = p->core_rw()->tf01; int flg = s[1]; double x; s += 2;
	for (int i=0, j=1; i<4; i++, j+=j) if (flg&j) s+=parse_num(&x,s), q[i] = (float)x;
	return 0; }

CH(pl){	if (s[1]<58) return p->key_op(11, s[1]-48, 0, cb->cnof());
	if (s[1]=='P') return p->add2mx_txdlv(0,0,0,0,0);
	const float * q = p->core_ro()->tf01; switch(s[1]) {
		case 'T': return p->plot_t(q[0], q[1], 512);
		case 'F': return p->plot_f(q[0], q[1], q[2], q[3], 512, false);
		default: return BXE_CENUM;
}}

CH(win){return (s[1]=='t') ? (p->show_tab(s[2]&7), 0) 
			   : p->wlg(p->m_node->wdat_raw(), s[1]&3, (8+(s[1]&4))<<20); }

CH(ky){ WrapCore * cr = p->core_rw();
	for (++s; s[0]>47 && s[1]>47 && s[2]>47 && s[3]>47; s+=4)
		cr->set_key(hex2(s)&127, 64*s[2] + s[3] - 3120);   return 0; }

CH(wav){int k = 0, wf = p->m_node->winflg(2048);
	switch(s[1]) {
		case 'F': p->m_bflg &= ~(k=WRF_W_ALL); p->m_bflg |= (((s[1]-48)<<WRF_W_SH)&k); goto flgw;
		case '2': k = WRF_W_CH2; goto flgs;
		case '-': k = WRF_W_PLOT; goto flgs;
		case 'W': case 0: return p->write_a20();
		default: return BXE_CENUM;
	}
flgs:	if (s[2]&1) p->m_bflg |= k; else p->m_bflg &= ~k;
flgw:	if (wf) p->w_a20(k);
	return 0;
}

CH(cut){int *q = p->m_node->etc()->i; double x = 0.0;
	switch(s[1]) {
		case '_': (s[2]&1) ? (*q|=INT_MIN) : (*q&=INT_MAX); break;
		case 'T': *q &= INT_MIN; parse_num(&x, s+2);
			  *q |= (INT_MAX & (int)lround(x*40320.0)); return 0;
		default : *q = atoi_h(s+1); break;
	}
	if (p->m_node->winflg(2048)) p->w_tlim(1);    return 0;
}

CH(gr){	int k, i = s[2] & 7, wdf = 65536;
	WrapCore * core = p->core_rw();
	switch (s[1]) {
		case 'd': wdf |= intv_cmd_c(core->grdim+i, s+3, 2, 51, 0x33330801) << (1+2*i); break;
		case 'R': i&=1; core->grdim[i+8] = (s[3]-48) & 127; wdf |= 131072 << (2*i); break;
		case '*':
			  for (i=0; i<8 && (unsigned int)(k=s[i+2]-50)<50u; i++) core->grdim[i] = k+2;
			  if (i==8 && s[10]) core->grdim[8] = (s[10]-48) & 127,
				  	     core->grdim[9] = (s[11]-48) & 63, i+=2;
			  wdf |= ((1<<(2*i))-1) & 0xaaaaa; break;
		default: return BXE_CENUM;
	}
	p->w_dim(wdf); return 0;
}

#define AW_CTAB {'+'+256, (cmd_t)c_c2k}, {'P'+256, (cmd_t)c_pl }, {'t', (cmd_t)c_tf}, \
	        {'W'+256, (cmd_t)c_win}, {'A'+256, (cmd_t)c_wav}, {'#', (cmd_t)c_gr}, {'T', (cmd_t)c_cut}, \
	        {'.'+256, (cmd_t)c_stp}, {'G'+256, (cmd_t)c_gcl}, {'k', (cmd_t)c_ky}, {'m', (cmd_t)c_gmd}

/////// box (wr) ////////////////////////////////////////////////////////////

DWrapGen::DWrapGen(ABoxNode * nd) : AWrapGen(nd) {
	m_sob.set(WrapSOB_default(0)); m_sfbx[0] = m_sfbx[1] = 0;
	delayed_clip_upd(); }

DWrapGen::DWrapGen(ABoxNode * nd, const DWrapGen * p) : AWrapGen(nd, p) {
	m_sob.from(p->m_sob); memcpy(m_sfbx, p->m_sfbx, 2*sizeof(BoxGen*));
	upd_conn(); delayed_clip_upd(); }

int DWrapGen::get_nm2(char * to) {
	BoxGen *bxs = m_sfbx[0], *bxf = m_sfbx[1];
	if (!bxs) return to[0]=to[1]=63, 2;;
	if (!bxf) return bxs->node()->get_nm2(to);
	if (bxs==box_bookmark[1] || bxs==box_bookmark[2]) return bxf->node()->get_nm2(to);
	char buf[4]; bxs->node()->get_nm2(to); bxf->node()->get_nm2(buf); to[1] = buf[0]; return 2;
}

int DWrapGen::sl_bv() {
	WrapSOB * sob = m_sob.ro();
	unsigned char *p = sob->m_scl[0].ro()->qs,
	     	      *q = sob->m_scl[1].ro()->qs;
	int sni=((m_bflg>>2)&31)+6, fni=((m_bflg>>12)&31)+2, r = 0, y = sob->m_core.ro()->grdim[1];
	for (int i=0; 2*i  <sni; i++) r |= 1<<( p[i]    &7);
	for (int i=0; 2*i+1<sni; i++) r |= 1<<((p[i]>>4)&7);
	for (int i=0; 2*i  <fni; i++) r |= 1<<( q[i]    &7);
	for (int i=0; 2*i+1<fni; i++) r |= 1<<((q[i]>>4)&7);
	return slbv_2(r, y);
}

int DWrapGen::add2mx_txdlv(int trg, int flg, int dly, int lim, const double *vs) {
	if (!m_sfbx[0]) return EEE_NOEFF;
	double in[36], v11[11];
	int ni, no, nv = 0, ncf = (flg|=pflg()) & WRF_NOCON, v0f = !!(flg & WRF_SKIPV0),
	    trf = !trg && trk_rec_trg;
	BVFOR_JM(flg&255) v11[j] = vs[nv++];
	WrapSOB * sob = m_sob.ro();
	sob->prep11(v11, flg | (768^(v0f<<9)), m_xys6);
	if (m_sfbx[1]) {
		sob->v(in, v11, 64, ncf + (ni=(m_bflg>>12)&31) + 1);
		int ocf = 0x7c01; // TODO (?)
		if ((trg = mx_add_filter(trg, m_sfbx[1]->model(), ni, in, ocf)) < 0) return trg;
	}
	char udio[4], *uds = sob->m_scl[0].ro()->updn; udio[0] = uds[0]; udio[1] = uds[1];
	udio[2] = ni = (m_bflg>>2)&31; 
	udio[3] = no = (m_bflg>>7)&31;
	sob->v(in + v0f, v11, v0f, ncf + ni + 6 - v0f);
	if (v0f) in[0] = in[1] = 1.0;
	if (flg & WRF_MONO) in[5] = 0.0, in[0] *= M_SQRT2;
	int ocb = no>1 ? 0x402 : 0x7c01; // TODO
	int r = mx_add_box(trg, m_sfbx[0]->mk_box(), udio, in, ocb, dly, lim); if (r<0) return r;
	if (trf) trk_rec(m_node, r);
	return r;
}

int DWrapGen::qdiff(DWrapGen * that) {
	return (this!=that) && (m_sob.ro() != that->m_sob.ro() || memcmp(m_xys6, that->m_xys6, 8)); }

int DWrapGen::set_sf_2(int ff, BoxGen * bx) {
	int ec, sh = 2 + (10 &- ff), ni, no, msk = ~(1023<<sh);
	if (bx) ni = bx->n_in(), no = bx->n_out();
	if (ff && bx && (!ni||no!=1)) return (bx==m_sfbx[1]) ? (set_boxp(m_sfbx+1,0), m_bflg&=msk, BXE_FILTRM) 
						             : BXE_FILTNIO;
	if ((ec = set_boxp(m_sfbx + ff, bx)) < 0) return ec;
	m_bflg &= msk; if (bx) m_bflg |= (32*no+ni)<<sh;
	return 0;
}

int DWrapGen::set_sf(int ff, BoxGen * bx) {
	BoxGen ** ppbx = m_sfbx + (ff&=1);
	if (bx==this || bx==*ppbx) return 0;
	int ec, wf = wnfl(2048);
	if (bx && bx->node()->cl_id()=='w') {
		DWrapGen * that = dynamic_cast<DWrapGen*> (bx); if (!that) return BXE_SORRY;
		if ((ec=set_sf_2(ff, that->m_sfbx[ff]))<0) return ec;
		sob_from(4+ff, that, 1);
		sob_from(2+ff, that, 1);
		w_mini(); w_slbv(2); if (wf) { BXW_GET; wlg(bxw_rawptr, 0, 0); wlg(bxw_rawptr, ff+1, 0); } }
	else {  if ((ec=set_sf_2(ff, bx))<0) return ec;
		w_mini(); w_slbv(2); if (wf) { BXW_GET; wlg(bxw_rawptr, ff+1, 0); } }
	return 0;
}

void DWrapGen::notify_nio(BoxGen * bx) {
	if (DBGC) log("wr%d/notify_nio: %d", id(), bx->id());
	int ec = 0, ni, no, flg = (bx==m_sfbx[0]) + 2*(bx==m_sfbx[1]);
	if (flg&1) m_bflg &= ~0xffc, m_bflg |= 4*(32*bx->n_out()+bx->n_in());
	if (flg&2) ((no=bx->n_out())==1 && (ni=bx->n_in())) ? 
		(m_bflg&=~0x3ff000, m_bflg |= (32*no+ni)<<12) : (set_sf(1, 0), ec = BXE_FILTRM);
	if (ec) gui_errq_add(ec, "wrap/nio"); if (!flg || !wnfl()) return;
	sthg * q = m_node->wdat_raw(); if (!q) return gui_errq_add(BXE_WTF, "wrap/nio/w");
	if (flg&1) wlg(q,1,0); if (flg&2) wlg(q,2,0); 
}

int DWrapGen::sob_from(int ix, BoxGen * bx0, int bxf) {
	DWrapGen * that = dynamic_cast<DWrapGen*> (bx0); if (!that) return NDE_EXPWRAP;
	int ec, ff = ix&1;
	if (!that) return bx0 ? BXE_TYDIFF : BXE_ARGNBX;
	if (!ix) { if (bxf && ((ec=set_sf_2(0, that->m_sfbx[0]))<0 || (ec=set_sf_2(1, that->m_sfbx[1]))<0))
			return ec;
		   m_sob.from(that->m_sob); return 0; }
	WrapSOB *sob = SOB_RW(sob), *sob2 = that->m_sob.ro();
	switch (ix) {
		case 1: 	sob->m_core   .from(sob2->m_core);    return 0;
		case 4: case 5: sob->m_con[ff].from(sob2->m_con[ff]); return 0;
		case 6: 	sob->m_avol   .from(sob2->m_avol);    return 0;
		case 2: case 3: if (bxf && (ec=set_sf_2(ff, that->m_sfbx[ff]))<0) return ec;
				sob->m_scl[ff].from(sob2->m_scl[ff]); return 0;
		default: return BXE_IDX;
	}}

int DWrapGen::start_job_3(JobQ::ent_t * ent, char * arg) {
	int ec; switch(ent->i3f) {
		case 1:
			ent->plttwwii = 0x6b775924;
			WrapAutoVol * av; av = SOB_RW(sob)->avol_rw();
			if (arg && *arg && *arg!='w') {
				if ((ec = av->parse_dim(arg))<0) return ec;
			} else {
				BXW_GET; if ((ec = av->upd_xy12rt(WR_AVCONF))<0) return ec;
			}
			ec = mx_mkroot();
			return ec<0 ? ec : (ent->p = new WrapAVJob(ent, this, av, ec), 0);
		case 9: return NDE_PERM;
		default: return JQE_UNDEF;
	}}

int DWrapGen::save_sob(SvArg *p) { switch(p->st) {
	case 16:	  return m_sob.save(p);
	case 17:	  return m_sob.ro()->m_core.save(p);
	case 18: case 19: return m_sob.ro()->m_scl[p->st&1].save(p);
	case 20: case 21: return m_sob.ro()->m_con[p->st&1].save(p);
	case 22:	  return m_sob.ro()->m_avol.save(p);
	case 23: case 24: p->st = 3; return 1;
	default: log("wr: invalid SOB state %d, skipping", p->st); return p->st2=-1, 1;
}}

int DWrapGen::save2(SvArg * sv) {
	BXSV2_HEAD; CHKERR(save2_aw(sv));
	if (m_sfbx[0]) { CHKERR(f->sn("X$b", 3)); CHKERR(m_sfbx[0]->node()->sv_path(10)); }
	if (m_sfbx[1]) { CHKERR(f->sn("X$>", 3)); CHKERR(m_sfbx[1]->node()->sv_path(10)); }
	return r;
}

void DWrapGen::spec_debug() {
	log("wrap: flg=0x%x", m_bflg);
	m_sob.debug();
}

int DWrapGen::wlg(sthg * bxw_rawptr, int ix, int flg) {
	if (!bxw_rawptr) return NDE_NOWIN; if ((unsigned int)ix > 2u) return ix<0 ? BXE_RANGE : BXE_WLGDONE;
	int oid = w_oid(), gc = (0x605a5345>>8*ix) & 127, pf = pflg(),
	    of = (flg & WRF_SETOC) ? ((flg&WRF_OC) ? (WR_WLG|=(1<<ix),1) : (WR_WLG&=~(1<<ix),0))
		    		   : (WR_WLG >> ix) & 1;
	WrapSOB * sob = m_sob.ro();
	gui2.setwin(oid, 'w'); 
	if (of)                      gui2.wupd_0(gc, "S><$XW", 4), gui2.c1(48 + ix);
	else gui2.wupd_0(gc, ".+1"), gui2.wupd_0(gc, "S<>$XW", 4), gui2.c1(52 + ix);
	if (ix) {
		BoxGen * bx = m_sfbx[--ix];
		int ni = (m_bflg>>(2+10*ix))&31; ni -= (ni && ix);
		if (of) gui2.wupd_0(gc, ".+"), gui2.c1(49+ni), 
			sob->wl(oid, 65*ix, ni, pf|255, m_xys6, bx);
		gui2.ref_title(gc, bx ? bx->node() : 0, 3, "source\0 filter"+8*ix);
	} else {
		sob->m_scl[0].ro()->w_ud(oid, 3);
		if (of) gui2.wupd_0(gc, ".+9"), sob->wl(oid, 32, 6, pf|254, m_xys6, 0),
						sob->wl(oid, 96, 2, pf|254, m_xys6, 0);
	}
	return 0;
}

const char * DWrapGen::eff_rgb() {
	return (m_sfbx[1] ? m_sfbx[1]->node() : (m_sfbx[0] ? m_sfbx[0]->node() : m_node)) -> own_rgb(); }

void DWrapGen::wdat_cons(sthg * bxw_rawptr) {
	WrapAutoVol * av = m_sob.ro()->m_avol.ro();
	memcpy(WR_AVCONF, av ? av->xy12rt() : (const unsigned char*)"\x9\x9\x9\x9\x1\x32", 6);
	wdat_cons_aw(bxw_rawptr); }

void DWrapGen::w_sob(int ix, int sl) {
	int wf = wnfl(2560); if (!wf) return;
	if (ix<2) { if (wf&512) aux_window(); if (wf&2048) box_window(); return; }
	if (wf&512) {
		if (ix==1) aux_window();
		if (sl) w_slbv(2); 
		if (wf==512) return;
	}
	BXW_GETV("wr/w_sob");
	switch (ix) {
		case 1: if (!WR_TAB) w_tab0(-1); return;
		case 2: case 3: wlg(bxw_rawptr, 0, 0); wlg(bxw_rawptr, ix-1, 0); return;
		case 4: case 5: w_col0(pflg()); return;
		case 6:	memcpy(WR_AVCONF, m_sob.ro()->m_avol.ro()->xy12rt(), 6); w_avol(63, WR_AVCONF); return;
		default:	log("WTF: invalid index for w_sob (%d)", ix); return;
	}}

void DWrapGen::w_col0(int flg) {
	BXW_GETV("w_col0"); WrapSOB * sob = m_sob.ro();
	int gf = WR_WLG; if (!gf) return;
	int wlf = 32 + WRF_CFILT + ((flg&1) << WRF_XFILT_SH); flg &= WRF_PASS; flg |= wlf;
	int oid = w_oid(), gf1 = gf & 1, ns = (gf&2) ? ((m_bflg>> 2) & 31)     : 0,
	    				 nf = (gf&4) ? ((m_bflg>>12) & 31) - 1 : 0;
	if (DBGC) log("w_col: gf=%d, ns=%d, nf=%d", gf, ns, nf);
	if (gf&3) sob->wl(oid,    32*gf1, 6*gf1+ns, flg, m_xys6, m_sfbx[0]);
	if (gf&5) sob->wl(oid, 65+31*gf1, 2*gf1+nf, flg, m_xys6, m_sfbx[1]);
}
					     
int DWrapGen::av_guiconf(int c, const char * s) {
	BXW_GET; int i, k = 0x8040109;
	switch(c) {
		case 'x': case 'y': case 'z': i = c-'x'; break;
		case 'w': i = 3; break;  
		case 'r': i = 4; break;
		case 't': i = 5; k = 0x320a01ff; break;
		default: return BXE_UCMD;
	}
	intv_cmd_uc(WR_AVCONF+i, s, 1, k&255, k>>8);
        if (wnfl() && WR_TAB==3) w_avol(1<<i, WR_AVCONF);
	return 0;
}

void DWrapGen::w_avol(int f, unsigned char * s) { gui2.setwin(w_oid(), 'w');
	BVFOR_JM(f) gui2.wupd_i('Y', s[j], j+30); }

int DWrapGen::show_tab_2(sthg * bxw_rawptr, int i) { switch(i) {
	case 3: w_avol(63, WR_AVCONF); return 0;
	default:return BXE_CENUM;
}}

#undef CH
#define CH(X) BXCMD_H(DWrapGen, X)

CH(bw){	BoxGen * bx = p->m_sfbx[(s[1]>>2)&1]; return bx ? bx->node()->draw_window(16) : EEE_NOEFF; }

CH(in){ WrapSOB * sob = SOB_RWP(p,sob); 
	int g9 = p->m_node->gui9(), r = sob->icmd(s+1, p->m_bflg, g9, p->m_xys6, p->m_sfbx);
	if (r<0) return r;
	if  (r&4) p->w_slbv(2); 
	if  (r&1) p->m_bflg &= ~WRF_NOCON;
	if ((r&2) && (g9&1)) p->w_col0(p->pflg()|1);    return 0; }

CH(sf){ ANode * nd = 0; BoxGen * bx = 0;
	if (s[1]!=48||s[2]) {   if (!(nd = cb->lookup(s+1))) return BXE_ARGLU;
				if (!(bx = nd->box0())) return BXE_ARGNBX; }
	return p->set_sf((*s&16)>>4, bx); }

CH(ud){	int ix = (*s & 1) ^ 1, k = s[1];
	if ( (k<66) ? (k!=45+20*ix) : (k>122 || !((1<<(k&31))&0xa10a0)) ) return BXE_CENUM;
	SOB_RWP(p,sob)->scl_rw(0)->updn[ix] = k;
	if (p->wnfl()) p->m_sob.ro()->m_scl[0].ro()->w_ud(p->w_oid(), 1+ix);
	return 0; }

CH(so){	if (!s[1]||!s[2]) return BXE_NOARG;
	ANode * nd = cb->lookup(s+2);
	int i = s[1]&7,  ec = nd ? p->sob_from(i, nd->box0(), cb->before(0,5)) : BXE_ARGLU;
	if (ec>=0) p->w_sob(i, 1), p->w_slbv(2); return ec; }

CH(vt){	if (s[1] > 'q') return p->av_guiconf(s[1], s+2);
	WrapSOB * sob = SOB_RWP(p,sob); int ec = 0;
	switch(s[1]) {
		case 'D': sob->m_avol.set(0); break;
		case ':': ec = sob->avol_rw()->parse_dim(s+2); break;
		default:  return BXE_UCMD;
	}
	if (ec<0 || !p->wnfl()) return ec;
	BXW_GETP(p); if (WR_TAB==2) p->w_sob(6, 0);
	return ec;
}

BXCMD_DEF(DWrapGen) {    {8192+'\\', 0}, AW_CTAB, 
	{'U', c_ud}, {'b', c_sf}, {'V',c_vt},  {'i',c_in}, {'B'+256, c_bw}, 
	{'D', c_ud}, {'>', c_sf}, {'O', c_so}, {0, 0}    };

///////////// export /////////////////////////////////////////////////////////////

#define WARG  (static_cast<AWrapGen*>(bx))
#define WARGD (static_cast<DWrapGen*>(bx))
int setbox_wrap(ABoxNode * nd, BoxGen * _) { return new (ANode::a64()) DWrapGen(nd), 3; }
int setbox_wrap_qcp (ABoxNode * nd, BoxGen * bx) { return new (ANode::a64()) DWrapGen(nd, WARGD), 3; }
const char * wrap_rgb(BoxGen * bx) { return WARG -> eff_rgb(); }
int wrap_2mx_txdlv(BoxGen * bx, int trg, int xflg, int dly, int lim, double *v) { 
	return WARG -> add2mx_txdlv(trg, xflg, dly, lim, v); }
int wrap_nd_2mx(ABoxNode * bnd, int trg, double bpm, int dly) {
	BoxGen * bx = bnd->box(); if (!bx) return 0;
	int tf = bnd->etc()->i[0]; if(tf){ if (DBGC) log("nd2mx: tf0=%d, bpm=%g", tf, bpm); 
				           if(tf<0)return 0; else tf=(int)lround((double)tf*natural_bpm/bpm); }
	if (DBGC) log("nd2mx: tf=%d", tf); return  WARG -> add2mx_txdlv(trg, 0, dly, tf, 0); }
int wrap_qdiff(BoxGen * bx, BoxGen * b2) { return WARGD -> qdiff(static_cast<DWrapGen*>(b2)); } // TODO
AReader * wrap_avreader(BoxGen * bx) { return WARGD -> avreader(); } // TODO: typechk(?)
int wrap_dump_keytab(char * to, unsigned int * bv, short ** pk) {
	int v, n = 0; for (int i=0; i<4; i++) BVFOR_JM(bv[i])
		to[n] = hexc1(2*i+(j>>4)), to[n+1] = hexc1(j&15), v = pk[i][j],
		to[n+2] = 48+(v>>6), to[n+3] = 48+(v&63), n += 4;
	return n; }
int wrap_key_op(BoxGen * bx, int ky, int op, const char *s, int nof) { return WARG -> key_op(ky,op,s,nof); }
void wrap_set_trec(BoxGen * bx, int j) { WARG->m_trec = j; }
void wrap_init() { WrapScale::st_init(); DWrapGen::st_init(); }
