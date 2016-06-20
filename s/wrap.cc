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
#define WRF_SHADOW 0x08000000 // bx
#define WRF_MKSH   0x04000000 // bx
#define WRF_SLUPD  0x02000000 // bx
#define WRF_B_RSRV 0x01c00000 // bx (rsrv)
#define WRF_MVC (WRF_W_CH2|WRF_W_PLOT|WRF_MKSH|WRF_SLUPD|WRF_B_RSRV)
#define WRF_CLID(X) ('w' - (((X)>>25)&4))
#define WRF_S_NI   0x0000007c // bxD
#define WRF_S_NO   0x00000f80 // bxD
#define WRF_F_NI   0x0001f000 // bxD
#define WRF_F_NO   0x003e0000 // bxD
#define WRF_VFL    0x0000037c // bxS
#define WRF_IADJS  0x00010000 // bxS, 2mx
#define WRF_IADJT  0x00020000 // bxS, 2mx
#define WRF_COMP   0x00000400 // v8tr(S)
#define WRF_W_ALL  0x30000000
#define WRF_W_SH   28
#define WRF_NOREC  0x00040000 // 2m
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
#define WRF_GRMODE 3	      // bx*
#define DBGC  (debug_flags&DFLG_WRAP)
#define DBGCV (debug_flags&DFLG_VOLTAB)

#define WR_AVCONF ((unsigned char*)(bxw_rawptr[1].c))
#define WR_TAB   (bxw_rawptr[2].c[6])
#define WR_WLG   (bxw_rawptr[2].c[7])
#define WR_SLFLG (bxw_rawptr[2].c[5])

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
static const unsigned char i8tab[8] = {2,3,4,5,6,7,8,9};
static const double d8tab[8] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
static const int wr_flg_mv[4] = { WRF_MVC|3, 0, WRF_IADJS|WRF_IADJT|WRF_MVC|1023, 12 };
#define WR_FLG_MV(P) (wr_flg_mv+2*!!((P)->m_bflg&WRF_SHADOW))

static int slbv_2(int r,int y) { return r>>=2, r | (64&(((0x0e13173f>>(4*(bitcnt_8(r)&6)))&63)-y)); }
int wrap_dump_keytab(char * to, unsigned int * bv, short ** pk);

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
                inline int sxop(int c) { return (c==63 || c==m_xop) ? 0 : (m_xop=c, 0x28); } 
                inline int sfun(int c) { return c==m_fch ? 0 : ((*sfi_ent(m_fch=c)->f.f)(this), 0x24); }
		inline int ssrc(int c) { return ((c = ((c+(c<64))&7)) == m_si) ? 0 : (m_si=c, 0x22); }
		inline const char * fnm() const { return sfi_ent(m_fch)->f.s; }
		inline const char * snm() const { return "x\0 y\0 s1\0s2\0s3\0s4\0s5\0s6" + 3*m_si; }
		inline int sch() const { return m_si + (m_si>1 ? 47 : 120); }
		inline int fch() const { return m_fch ? m_fch : 'c'; }
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
		WrapCore(int uarg) : SOB(uarg), ktab(0) {}
		WrapCore(const WrapCore * that, int uarg);
		~WrapCore() { del_keytab(); }
		void ini_default(int _) { tf01[0]=tf01[2]=0.0; tf01[1]=1.0; tf01[3]=22050.0;
			xfd=63; memcpy(grdim, "\031\031""333333\012\012", 10); } 
                int save2(SvArg * p);
		void debug2() { log("WrapCore: xfd=0x%x(0x%x)", xfd, wr_ixtr_r(xfd)); } 
		int sv_tf01(char * to);
		void grcmd(const char * i8) { gui2.c4('#', 'g', grdim[0]+48, grdim[1]+48);
			gui2.c4(grdim[8]+48, grdim[9]+48, i8[0]+48, i8[1]+48); }
		void upd_grid(const char * i8, int gui9) {
			if (gui9&1) gui2.setwin(gui9|11,'w'), gui2.w0(), grcmd(i8);
			if (gui9&8) gui2.setwin(gui9| 9,'#'), gui2.w0(), grcmd(i8); }
		void slcmd(const char * i8, int flg) { gui2.c1('!');
			for (int i=2; i<8; i++, flg>>=1)
				(flg&1) ? gui2.c2(i8[i]+48, grdim[i]+47) : gui2.c2(44,44); }
		void i2v(double *v11, const char *i8, int flg) {
			BVFOR_JM(~flg&255) v11[j] = div1tab[(int)grdim[j]] * (double)i8[j]; }
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
                int fill_data(char * s, int ver, int cflg); // version: 0 1 2
                int pred2(int i, int j, int k, int l);
                int pred2o(int i, int j, int k, int l);
                int pred2t(int i, int j, int k, int l);
		unsigned char * xy12rt() { return m_xy12rt; }
		int parse_dim(const char * s);
		void set_c8(unsigned char * _) {}
		WrapAutoVol * copy(int uarg) { return new WrapAutoVol(uarg); } // no copy needed
        protected:
                static void unref2(WrapAutoVol* p, int ix);
                static int freeblk();
                static const int eol = -1234567890;

                short * m_dat;
		int m_hood[16];
                unsigned char m_xy12rt[8];
};

class WrapAVReader : public AReader {
	public:
		WrapAVReader(WrapAutoVol * av, int cflg) : m_av(av), m_cflg(cflg) {}
		virtual int line(char * s);
	protected:
                LWArr<char> m_dat_str;
		WrapAutoVol * m_av;
		int m_cflg;
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
		static const int sflg = 1;
		SOBDEF_64(WrapSOB);
		WrapSOB(int uarg) : SOB(uarg) {}
		WrapSOB(const WrapSOB * that, int uarg);
		void ini_default(int k);
                int save2(SvArg * p) { return p->st2=-1, 1; }
                void debug2();
		void v(double *to, const double * v11, int ix, int nf);
		void prep11(double *v11, int flg, const char * i8);
		WrapScale * sit_ro(int tix) { return m_scl[(tix>>6)&1].ro()->item(tix, 0); }
		WrapScale * sit_rw(int tix) { return SOB_RW(scl[(tix>>6)&1])->item(tix, 1); }
		int sit_cmd(int tix, const char * s) { if (DBGC) log("sit_cmd: 0x%x, \"%s\"", tix, s); 
						       return SOB_RW(scl[(tix>>6)&1]) -> cmdi(tix, s); }
		void unflg();
		int icmd(const char *s, int *pbf);
		void wl(int oid, int ix0, int n, int flg, const char *i8, BoxGen * bx);
		SOB_RW_F0(core, WrapCore) SOB_RW_F0(avol, WrapAutoVol)
		SOB_RW_F1(con, DblVec)    SOB_RW_F1(scl, WrapScVec)

		SOB_p<WrapCore> m_core;
		SOB_p<WrapAutoVol> m_avol;
		SOB_p<DblVec> m_con[2];
		SOB_p<WrapScVec> m_scl[2];
};

class SWrapTab : public SOB {
	public: 
		SOBDEF_64(SWrapTab);
		SWrapTab(int arg) : SOB(arg) { memset(cc,0,8); memcpy(ix, i8tab, 8); ce[0]=ce[1]=0; }
		SWrapTab(const SWrapTab * that, int arg) : SOB(arg) {          		memcpy(cc,that->cc,8);
			for (int i=0;i<2;i++) ce[i]= (double*)ANode::cp64(that->ce[i]); memcpy(ix,that->ix,8);}
		void ini_default(int k) { }
		~SWrapTab() {  }
		int save2(SvArg * sv);
		void debug2() {}
		double v1(int j, const double **src2);
		void w_jline(int j, double con);
		void w1(int j, int k);
		void set_xy(int j, int f, int v) { ix[j] = f ? (ix[j]&15)+16*v : (ix[j]&240)+v; } 
		void set_ab(int j, int f, double v) { 
			(ce[f] ? ce[f] : (ce[f] = (double*)(f ? ANode::z64() : ANode::cp64(d8tab))))[j] = v; }
		void set_iif(int j, int v) { unsigned char m=1<<j; if (v) iif|=m; else iif&=~m; }
		int set_all(int j, const char * a0, const char * a1);
		int ibv(int f8) { int r=0; BVFOR_JM(f8) r |= (1<<(ix[j]&15)) | (1<<(ix[j]>>4)); return r>>2; }

		unsigned char ix[8]; 
		  signed char cc[8];
		unsigned char iif, rsrv7[7], *rsrv8;
		double *ce[2];
};

class SWrapSOB : public SOB {
	public:
		static const int sflg = 1;
		SOBDEF_64(SWrapSOB);
		SWrapSOB(int uarg) : SOB(uarg) {}
		SWrapSOB(const SWrapSOB * that, int uarg);
		void ini_default(int k);
                int save2(SvArg * p) { return log("swrsob/save %d", p->st),  p->st2=-1, 1; }
                void debug2();
		SOB_RW_F0(con, DblVec) 
		SOB_RW_F0(tab, SWrapTab)

		SOB_p<DblVec> m_con;
		SOB_p<SWrapTab> m_tab;
};

class AWrapGen : public BoxGen {
	public:
		friend void wrap_set_trec(BoxGen * bx, int v); friend class SWrapGen;
		AWrapGen(ABoxNode * nd);
		AWrapGen(ABoxNode * nd, const AWrapGen * that);
                virtual ~AWrapGen();
		virtual int add2mx_txdlv(int trg, int xflg, int dly, int lim, const double *vs) = 0;
		virtual int get_nm2(char * to) = 0;
		virtual WrapCore * core_ro() = 0; 
		virtual WrapCore * core_rw() = 0;
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
		virtual int ifflg() const { return BIF_QCP; }
		inline int pflg() { return m_bflg & WRF_PASS; }
		void delayed_clip_upd();
		void set_tl(int t = INT_MAX) { m_node->etc()->i[0] = t; }
		int key_op(int k, int op, const char * xys, int nof);
	protected:
		typedef struct { int f; char xy[8]; } qsav_t;
		typedef int (acmd_t) (AWrapGen *, const char *, CmdBuf *);
		static int slf_conv(int k);
		static acmd_t c_c2k, c_cut, c_gcl, c_gmd, c_gr, c_ky, c_pl, c_stp, c_tf, c_wav, c_win, c_rvt,
			      c_xfd, c_flg;
		virtual int show_tab_2(sthg * bxw_rawptr, int i) = 0;
		virtual int wlg(sthg * bxw_rawptr, int ix, int flg) = 0;
		virtual void w_col0() = 0;
		virtual int sl_bv() = 0;
		inline void qsav(qsav_t * q) { q->f = m_bflg; memcpy(q->xy, m_xys6, 8); }
		inline void qrst(qsav_t * q) { m_bflg = q->f; memcpy(m_xys6, q->xy, 8); }
		inline int qset(const char *s, int l) { for (int k,i=0; i<l; m_xys6[i++] = k) 
			                if ((unsigned int)(k=s[i]-48)>50u) return BXE_PARSE; return 0;}
		inline void wdat_cons_aw(sthg * bxw_rawptr) { WR_TAB=0; WR_WLG=0; WR_SLFLG=64; }
		inline int add2ctl(int mxbi, int mxky) {
			return mx_c_add(m_mxctl?m_mxctl:(m_mxctl=mx_mkctl(this)), mxbi, mxky); }
		int save2_aw(SvArg * sv);
		void gr_rgbc(){ char buf[4]; get_nm2(buf); gui2.c3('M',buf[0],buf[1]); gui2.sn(v_rgb(),6); }
		int plot_t(double t0, double t1, int n, int flg);
		int plot_f(double t0, double t1, double f0, double f1, int n, int flg);
		int mx1(int f=0) { int k, r = mx_mkroot(); if (r<0) return r;
			return (k=add2mx_txdlv(r,f,0,INT_MAX,0))<0 ? (mx_del(r),k) : r; }
		int qcopy(int tf, int nof) { return tf ? trk_glob_paste(this, nof)
			: Node::mk(0, ClipNode::kcp(1), 0, 'W'^((m_bflg>>24)&4), nof|NOF_NOIDX, 0, this); }
		int batch_calc(double * to0, double * to1, int skip, int n, int nch);
		int write_a20();
		int lim8(const char * q);
		void w_gr_xyk(int x0, int y0, int f);
		void w_mini();
		void w_gmd() { gui2.setwin(w_oid(), 'w'); gui2.wupd('m'); gui2.c2('+', 48+(m_bflg&3)); }
		void w_a20(int flg);
		void w_tlim(int flg);
		void w_dim(int gis6, int wf);
		void w_tab0(int f20);
		void w_slbv(int flg); // 1: force 2: setwin
		int show_tab(int i);
		int wlg_vis(int m) { if (!wnfl()) return 0; BXW_GET; return (WR_WLG &  m); }
		int tab_vis(int j) { if (!wnfl()) return 0; BXW_GET; return (WR_TAB == j); }
		
		char m_xys6[8];
		int m_bflg; 
		unsigned short m_mxctl;
		unsigned char m_trec, m_rsrv;
};

class DWrapGen : public AWrapGen {
	public:
		static void st_init() { cmd_init(); }
		DWrapGen(ABoxNode * nd);
		DWrapGen(ABoxNode * nd, const DWrapGen * that);
                virtual ~DWrapGen() {}
		virtual int add2mx_txdlv(int trg, int xflg, int dly, int lim, const double *vs);
		virtual WrapCore * core_ro() { return m_sob.ro()->m_core.ro(); } 
		virtual WrapCore * core_rw() { return SOB_RW(sob)->core_rw(); }
		virtual int save2(SvArg * sv); 
		virtual void spec_debug();
		virtual void notify_nio(BoxGen * bx);
		virtual BoxGen * qcp3(ABoxNode * nd) { return new (ANode::a64()) DWrapGen(nd, this); }
		virtual int save_sob(SvArg *p);
		virtual void wdat_cons(sthg * p);
		virtual int start_job_3(JobQ::ent_t * ent, char * arg);
		virtual const char * v_rgb();
		virtual int get_nm2(char * to);
		WrapAutoVol * avol() { return m_sob.ro()->m_avol.ro(); }
		int avj_state() { return jobq.jst(m_node, 1); }
		int qdiff(DWrapGen * that); 
		int sob_from(int ix, BoxGen * bx0, int bxf);
		AReader * avreader(int cflg) { return new WrapAVReader(SOB_RW(sob)->avol_rw(), cflg); }
        protected:
		virtual int sl_bv();
		virtual int show_tab_2(sthg * bxw_rawptr, int i);
		virtual int wlg(sthg * bxw_rawptr, int ix, int flg);
		virtual void w_col0() { w_col0_d(0); }
		void upd_conn() {Node::set_conn_2(m_node, BoxGen::node0(m_sfbx[0]),BoxGen::node0(m_sfbx[1]));}
		BXCMD_DECL(DWrapGen) c_in, c_bw, c_sf, c_ud, c_vt, c_so;
		int set_sf(int ff, BoxGen * bx);
		int set_sf_2(int ff, BoxGen * bx);
		int xfd_chk(int ff, int tf);
		void upd_updn(int flg);
		int grid_cmd(const char * s);
		void w_sob(int ix, int sl);
		void w_col0_d(int flg); // 1:xfdep 0:noncon
		void w_avol(int f, unsigned char * s);
		int wupd1(int col, int ix1, int ix2 = 333);
		int av_guiconf(int c, const char * s);
		int avj_cmd(int stp);
		SOB_p<WrapSOB> m_sob;
		BoxGen * m_sfbx[2];
};

class SWrapGen : public AWrapGen {
	public:
		static void st_init() { cmd_init(); }
		SWrapGen(ABoxNode * nd);
		SWrapGen(ABoxNode * nd, const SWrapGen * that);
		SWrapGen(ABoxNode * nd, AWrapGen * fr, int flg);
                virtual ~SWrapGen() {}
		virtual BoxGen * qcp3(ABoxNode * nd) { return new (ANode::a64()) SWrapGen(nd, this); }
		virtual int save2(SvArg * sv); 
		virtual int save_sob(SvArg *p);
		virtual int get_nm2(char * to) { return (m_trg) ? m_trg->get_nm2(to) : (to[0]=to[1]=63, 2); }
		virtual int add2mx_txdlv(int trg, int xflg, int dly, int lim, const double *vs);
		virtual WrapCore * core_ro() { return m_core.ro(); } 
		virtual WrapCore * core_rw() { return SOB_RW(core); }
		virtual const char * v_rgb() { return m_trg ? m_trg->v_rgb() : "EEE%%%"; }
		virtual void spec_debug();
		virtual void wdat_cons(sthg * p) { return wdat_cons_aw(p); }
		int sob_from(int ix, BoxGen * bx0);
		int set_trg(AWrapGen * q);
	protected:
		virtual int show_tab_2(sthg * bxw_rawptr, int i);
		virtual int wlg(sthg * bxw_rawptr, int ix, int flg);
		virtual void w_col0() { w_col0_s(m_bflg&1020); }
		virtual int sl_bv() {
		        return slbv_2(m_ssob.ro()->m_tab.ro()->ibv((m_bflg>>2)&255), m_core.ro()->grdim[1]); }
		void v8tr(double *to, const double *src, int flg);
		void upd_conn() {Node::set_conn_2(m_node, BoxGen::node0(m_trg), 0); }
		void w_sob(int ix, int sl) {  }
		void w_trg2(){ gui2.setwin(w_oid(),'w'); gui2.ref_title('E', BoxGen::node0(m_trg),3,"target");}
		void w_trg() { gui2.setwin(w_oid(),'w'); gui2.own_title(); w_mini(); w_trg2(); }
		void w_vfl() { BXW_GETV("vfl"); if (WR_WLG&1) gui2.setwin(w_oid(),'w'), gui2.t_sn("&",1),
							      gui2.hexn(m_bflg>>2,2); }
		void w_jtab(int flg);
		void w_col0_s(int filt);
		void jt_show() { gui2.wupd_0('E', ".*"), gui2.hex4((m_bflg&1020)+3); }
		void cur_c0(double *to8, int flg) { 
			double v8[8]; m_core.ro()->i2v(v8, m_xys6, 0); v8tr(to8, v8, pflg()|flg); }
		BXCMD_DECL(SWrapGen) c_trg, c_so, c_vfl, c_tab, c_scf;
		SOB_p<WrapCore> m_core;
		SOB_p<SWrapSOB> m_ssob;
		AWrapGen * m_trg;
};

SOB_INIFUN(WrapScVec, 2)
SOB_INIFUN(WrapCore,  1)
SOB_INIFUN(WrapSOB,   1)
SOB_INIFUN(SWrapTab,  1)
SOB_INIFUN(SWrapSOB,  1)
/////// core SOB (grid config) //////////////////////////////////////////////

WrapCore::WrapCore(const WrapCore * that, int uarg) : SOB(uarg), xfd(that->xfd) {
	if (DBGC) log("wrcore copy: %p -> %p", that, this);
	memcpy(tf01, that->tf01, 16); memcpy(grdim, that->grdim, 10); 
	key_t * oktab = that->ktab; if (!oktab) { ktab = 0; return; }
	(ktab = (key_t*)ANode::a64())->n5 = oktab->n5;
	for (int i=0;i<4;i++) if ((ktab->bv[i]=oktab->bv[i])) ktab->pk[i] = (short*)ANode::cp64(oktab->pk[i]);
}

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

int WrapCore::sv_tf01(char * to) {
	float *q = tf01, *q0 = WrapCore_default(0)->tf01;
	int r = 0, flg = 0; 
	for (int i=0; i<4; i++) if (fabs(q[i]-q0[i])>1e-5) flg |= (1<<i);
	if (!flg) return 0; else memcpy(to, "X$t0", 4), to[3]+=flg, r = 4;
	for (int i=0; i<4; i++) if (flg & (1<<i)) to[r] = '#', doub2hx(to+r+1, q[i]), r+=17;
	return to[r]=10, r+1;
}

int WrapCore::save2(SvArg * sv) {
	char buf[600]; memcpy(buf, "X$#*", 4);
	for (int i=0; i<10; i++) buf[i+4] = grdim[i]+48;   buf[14] = 10;
	int l = (xfd == 63) ? 15 : 15+sprintf(buf+15,"X$x%02x\n", wr_ixtr_r(xfd));
	if (ktab) memcpy(buf+l, "X$k", 3), l += 4+wrap_dump_keytab(buf+l+3, ktab->bv, ktab->pk), buf[l-1] = 10;
	return sv->st2=-1, l+=sv_tf01(buf+l), sv->out->sn(buf, l);
}

/////// scale ///////////////////////////////////////////////////////////////

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
	default: log("wrscale: invalid command: \"%s\"", s-1); return EEE_PARSE;
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
	if ((r&0x80000002)==2)  i=(ix>>1)&31, j=4*(ix&1), qs[i] &= ~(7<<j),
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

#define AVLOOP4(F) for (int i=0; i<nx; i++) { for (int j=0; j<ny; j++) {     \
		       for (int k=0; k<n1; k++) { for (int l=0; l<n2; l++) {  \
			   if (vol_ix(i,j,k,l)!=cnt) return VTE_IXDIFF;        \
			   m_dat[cnt++] = F (i,j,k,l) + tr.get_short_k(ver); }}}}

int WrapAutoVol::fill_data(char *s, int ver, int cflg) {
        B91Reader tr; tr.init(s);
        int cnt=0, nx=m_xy12rt[0], ny=m_xy12rt[1], n1=m_xy12rt[2], n2=m_xy12rt[3];
	if (DBGCV) log("fill_data: compat flg = %d", cflg);
	if (cflg) { AVLOOP4(pred2o); }
	else 	  { AVLOOP4(pred2);  } 
	return 0;
}

#define PREDN(X,Y) m_dat[ix+m_hood[4*X+Y]]

#define PRED1(I) return PREDN(I,I)
#define PRED2(I,J) return PREDN(I,I) + PREDN(J,J) - PREDN(I,J) 
#define PRED3(I,J,K) return ((2*(PREDN(I,I)+PREDN(J,J)+PREDN(K,K))  \
			       -(PREDN(I,J)+PREDN(I,K)+PREDN(J,K)))*5461+8192) >> 14
		  
int WrapAutoVol::pred2(int i, int j, int k, int l) {
	int ix = vol_ix(i,j,k,l), flg = !!i+2*!!j+4*!!k+8*!!l;
	switch(flg) {
		case  0: return 0;
		case  1: PRED1(0);
		case  2: PRED1(1);
		case  4: PRED1(2);
		case  8: PRED1(3);
		case  3: PRED2(0,1);
		case  5: PRED2(0,2);
		case  9: PRED2(0,3);
		case  6: PRED2(1,2);
		case 10: PRED2(1,3);
		case 12: PRED2(2,3);
		case  7: PRED3(0,1,2);
		case 11: PRED3(0,1,3);
		case 13: PRED3(0,2,3);
		case 14: PRED3(1,2,3);
		case 15: return ((PREDN(0,0)+PREDN(1,1)+PREDN(2,2)+PREDN(3,3)+1) -
		    ((4096+5461*(PREDN(0,1)+PREDN(0,2)+PREDN(0,3)+PREDN(1,2)+PREDN(1,3)+PREDN(2,3)))>>14)) >>1;
		default: return bug("pred2b"), 0;
	}}

int WrapAutoVol::pred2o(int i, int j, int k, int l) {
        if (!(i|j|k|l)) return 0;
        int ix=vol_ix(i,j,k,l), ixn=0, ixd[4];
        if (i) ixd[ixn++] = -m_hood[0];
        if (j) ixd[ixn++] = -m_hood[5];
        if (k) ixd[ixn++] = -m_hood[10];
        if (l) ixd[ixn++] = -m_hood[15];
        if (ixn==1) return m_dat[ix - ixd[0]];
        int prlg_sum = 0, prlg_cnt = ixn * (ixn-1) / 2;
        for (int i=0; i<ixn; i++) 
                for (int j=i+1; j<ixn; j++) 
                        prlg_sum += m_dat[ix-ixd[i]] + m_dat[ix-ixd[j]] - m_dat[ix-ixd[i]-ixd[j]];
        return (prlg_sum + (prlg_cnt/2)) / prlg_cnt;
}
 
int WrapAutoVol::pred2t(int i, int j, int k, int l) {
	static int mdif = 0, c0=0, cm=0, cp=0;
	int x = pred2o(i,j,k,l), y= pred2(i,j,k,l), d = x-y; 
	if (!d) c0++; else if (d<0) cm++; else cp++;
	if (d) d=abs(d), log("pred2(%d,%d,%d,%d): %d != %d (max:%d, <:%d =:%d >:%d)", i,j,k,l, x,y, d>mdif ? (mdif=d):mdif, cm, c0, cp);
	return x;
}

int WrapAutoVol::save2(SvArg * p) {
	Clock clk; int c1=0,c2=0,c3=0,c4=0; if (DBGCV) clk.reset();
	AOBuf * f = p->out;
        int total = 1; for (int i=0; i<4; i++) total *= m_xy12rt[i];
        short difftab[total];
        int ec, r = 1, cnt = 0;
        int nx=m_xy12rt[0], ny=m_xy12rt[1], n1=m_xy12rt[2], n2=m_xy12rt[3];
        for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                        for (int k=0; k<n1; k++) {
                                for (int l=0; l<n2; l++) {
                                        if (vol_ix(i,j,k,l)!=cnt) bug("voltab/save2: ix4(%d)!=cnt(%d)",
									     vol_ix(i,j,k,l), cnt);
                                        int pre = pred2(i,j,k,l);
                                        int dif = (m_dat[cnt] - pre);
                                        difftab[cnt++] = dif;
                                }}}}
	
        if (DBGCV) c1 = clk.reset(), log_sortedblk(difftab, total, 1, "difftab2:");
        int cost[3];    cost[0] = b91_cost0(difftab, total);
			cost[1] = b91_cost1(difftab, total);	
			cost[2] = b91_cost2(difftab, total);	
	if (DBGCV) c2 = clk.reset();
        int ty = cost[1]<cost[2] ? cost[1]<cost[0] : 2*(cost[2]<cost[0]);
        if (debug_flags & DFLG_VOLTAB) log("voltab_cost(%d): %d %d %d --> %d", p->cn->id(), cost[0], cost[1], cost[2], ty);
        B91Writer wr; wr.put_short_tpn(ty, difftab, total);
        int l = wr.n_bytes();
        char * s = wr.get_str();
	char buf[20]; memcpy(buf, "X$V:", 4);
	for (int i=0; i<5; i++) buf[4+i] = m_xy12rt[i]+48; buf[9] = ':';
	sprintf(buf+10, "%03d", m_xy12rt[5]); memcpy(buf+13, "\nN$V", 4);
	CHKERR(f->sn(buf, 17));
	if (DBGCV) c3 = clk.reset();
        while (l) {
		CHKERR(f->sn("\n<", 2));
                int k = (l<75) ? l : 75;
		CHKERR(f->sn(s, k)); s+=k; l-=k;
	}
	CHKERR(f->sn("!\n\"\n#\n"+2*ty, 2)); 
	if (DBGCV) c4 = clk.reset(), log("avol/save clk: %d %d %d %d", c1, c2, c3, c4);
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
	m_hood[15] = -1; for (int i=3; i>0; i--) m_hood[5*i-5] = m_hood[5*i]*m_xy12rt[i];
	for (int i=0;i<4;i++) for (int j=0;j<i;j++) m_hood[4*i+j] = m_hood[4*j+i] = m_hood[5*i]+m_hood[5*j];
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
	m_trg = vtab->dat(); if (!m_trg) bug("avjob / zero trg");
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
			m_bx -> add2mx_txdlv(m_mxid, WRF_AVJOB|WRF_NOCON, 0, INT_MAX, buf);
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
	return k ? m_av->fill_data(m_dat_str.p(), s[l]-33, m_cflg) : 1;
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

void WrapSOB::v(double *to, const double *v11, int ix, int nf) { int i = (ix>>6)&1;
	wrap_fill(to, m_scl[i].ro()->pps, m_con[i].ro()->p, m_con[i].ro()->bv, v11, ix&63, nf); }

void WrapSOB::prep11(double *v11, int flg, const char * i8) {
	WrapCore * core = m_core.ro();  core->i2v(v11, i8, flg);
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

// 1:lbl 2:src 4:fun 8:xf  32:c 64:v0 128:v1   256:(xf=='=') 512: slref 1024: unflg
// ret: ff:col ff00:row ff0000:row2 1<<24:dep
int WrapSOB::icmd(const char *s, int *pbf) { 
	WrapCore * core = m_core.ro(); 
	int ix0 = 16*(s[0]&7) + hxd2i(s[1]), ix = wr_ixtr(ix0),
	    i = (ix>>6)&1, j = ix&63, xfd0 = core->xfd, r = ((xfd0==ix)<<24) | (ix0<<8);
	if (j > 39) return BXE_IDX; else s += 2;
	switch(*s) {
		case 'X': return SOB_RW(core)->xfd = ix, 0; // (old) save file only
		case '=': if (*pbf & WRF_NOCON) unflg(), *pbf &= ~WRF_NOCON;
			  return (s[1]=='x' && !s[2]) ? (SOB_RW(con[i])->rm(j), r|4)
						      : (parse_num(SOB_RW(con[i])->addp(j), s+1), 0);
		case 'F': if (s[1]!='=') return r | sit_cmd(ix,s) | (xfd0==ix ? (SOB_RW(core)->xfd=63,8) : 0);
			  return (xfd0==ix) ? 0 : (SOB_RW(core)->xfd = ix, r | 0x1000008 | ((xfd0^63)<<16)
					  				     | sit_cmd(ix,"F.") );
		default:  return r | sit_cmd(ix, s);
	}}

void WrapSOB::debug2() {
	log("WrapSOB: flg:0x%x core s_scl f_scl s_con f_con av", m_u8_refcnt>>24);
	m_core.debug(); m_scl[0].debug(); m_scl[1].debug();
			m_con[0].debug(); m_con[1].debug(); m_avol.debug(); }

void WrapSOB::ini_default(int k) {
	m_core.set(WrapCore_default(k));
	for (int i=0; i<2; i++) m_con[i].set(DblVec_default(i+1)),
			        m_scl[i].set(WrapScVec_default(i)); }

/////// swrap SOBs //////////////////////////////////////////////////////////

int SWrapTab::save2(SvArg * sv) {
	int ec, r = 1;
	const double * q0 = ce[0] ? ce[0] : d8tab, *q1 = ce[1];
	if(q1) { for (int j=0; j<8; j++) { CHKERR(sv->out->pf("X$j%c%c%02x%02x%s$%s\n", j+48, 42+((iif>>j)&1),
						    cc[j],ix[j], dbl2str_s(0,q0[j]), dbl2str_s(1,q1[j]))); }}
	else   { for (int j=0; j<8; j++) { CHKERR(sv->out->pf("X$j%c%c%02x%02x%s\n",    j+48, 42+((iif>>j)&1),
						    cc[j],ix[j], dbl2str_s(0,q0[j]))); }}
	return sv->st2 = -1, r;
}

double SWrapTab::v1(int j, const double **src2) {
	int k = ix[j], k0 = k&15, k1 = k>>4, ii = iif & (1<<j);
	const double * src = src2[!!ii];
	double x, v = 0.0;
	if (k0) x = !ce[0] ? 1.0 : ce[0][j], v += ((k0-=2)<0) ? x : x*src[k0];
	if (k1 &&    ce[1])   x =  ce[1][j], v += ((k1-=2)<0) ? x : x*src[k1];
	return ii ? (v + (double)cc[j])  :   (x = cc[j] ? v+.01*cc[j] : v, x<0.0?x:(x>1.0?1.0:x));
}

int SWrapTab::set_all(int j, const char * a0, const char * a1) {
	set_iif(j, *(a0++) & 1);
	if (!a0[0]||!a0[1]||!a0[2]||!a0[3]||!a0[4]) return BXE_PARSE;
	int k = qh4rs(a0); cc[j] = k>>8; ix[j] = k&255;
	set_ab(j, 0, at0f(a0+4)); set_ab(j, 1, a1 ? at0f(a1) : 0.0); return 0;
}

void SWrapTab::w_jline(int j, double con) {
	gui2.t0(); gui2.c2('-'-2*((iif>>j)&1), j|48); gui2.hex4(256*ix[j]+((unsigned char) cc[j]));
	gui2.hdbl(con); gui2.hdbl(ce[0] ? ce[0][j] : 1.0); gui2.hdbl(ce[1] ? ce[1][j] : 0.0);
}

void SWrapTab::w1(int j, int k) {
	gui2.wupd('E', 8*j+k+16);
	if (k==4) (k=cc[j])<0 ? gui2.c4('x', 48+((-k)>>4), hexc1((-k)&15), '-')
			      : gui2.c3('x', 48+(  k >>4), hexc1( k  &15));
	else gui2.c2('c', 48+(k ? ix[j]>>(4*k-20) : (iif>>j)&1) );
}

SWrapSOB::SWrapSOB(const SWrapSOB * that, int uarg) : SOB(uarg) {
	if (DBGC) log("s/wrsob copy: %p -> %p", that, this);
	m_con.from(that->m_con); m_tab.from(that->m_tab); m_u8_refcnt |= (that->m_u8_refcnt&0x7f000000); }

void SWrapSOB::debug2() {
	log("SWrapSOB: flg:0x%x core", m_u8_refcnt>>24);
}	

void SWrapSOB::ini_default(int k) {
	log("SWrapSOB::ini_default %d", k);
	m_con.set(DblVec_default(0));
	m_tab.set(SWrapTab_default(0));
}

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
	gui2.c1(36); gr_rgbc(); if (m_bflg&WRF_SLUPD) gui2.c1('u');
	WrapCore * cr = core_ro();
	cr->slcmd(m_xys6, 63); w_slbv(1); cr->grcmd(m_xys6);
	if (cr->ktab) gui2.c1('B'), gui2.wr_keylist(cr->ktab->bv, cr->ktab->pk); 
	return 16*'#'+9;
}

int AWrapGen::save2_aw(SvArg * sv) {
	BXSV2_HEAD; int k, m, l = 12, bf = m_bflg;
	char buf[24]; memcpy(buf, "X$G", 3); buf[3] = 51 + 3*!(bf&WRF_NOCON);
	for (int i=0; i<8; i++) buf[i+4] = m_xys6[i]+48;
	if ((k=packflg(bf, WR_FLG_MV(this)))>=0)
		m = qh4(k), memcpy(buf+l, "\nX$F", 4), buf[l+4] = 48+(k>>16), memcpy(buf+l+5, &m, 4), l += 9;
	buf[l] = 10; CHKERR(f->sn(buf, l+1));
	if ((k=m_node->etc()->i[0])!=INT_MAX) { CHKERR(sv->out->pf("X$T%x\n", k)); }
	CHKERR(m_node->sv_wr_backref()); return r;
}

void AWrapGen::w_tlim(int f) {
	int k = m_node->etc()->i[0], t = k&INT_MAX, mf = k<0;
	gui2.setwin(w_oid(), 'w'); 
	if (f&1) gui2.wupd_i1('Y', mf,  20);
	if (f&2) gui2.wupd_d('Y', t==INT_MAX ? -1.0 : (double)t/40320.0, 22);
	if (f&4) gui2.wupd_i1('Y', !!(m_bflg&WRF_MKSH),  23);
	if (f&8) gui2.wupd_i1('Y', !!(m_bflg&WRF_SLUPD), 24);
}

int AWrapGen::show_tab(int i) {
	BXW_GET; if (i<0) i=WR_TAB; else if (i==WR_TAB) return 0; else WR_TAB=i;
	gui2.setwin(w_oid(), 'w'); gui2.wupd_0('Y', ".W"); gui2.c1(48+i);
	switch(i) {
		case 0: w_tab0(-1); return 0;
		case 1: w_a20(-1); return 0;
		case 2: w_tlim(-1); return 0;
		default: return show_tab_2(bxw_rawptr, i);
	}}

void AWrapGen::box_window() {
	gui2.cre(w_oid(), 'w'); if (m_node->cl_id()=='s') gui2.c1('s');
	gui2.own_title(); show_tab(-1); w_mini(); w_gmd();
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
	char buf[8], xysav[8]; int ec, nb, bfsav = m_bflg, tlim = INT_MAX;
	if (k<0) { 
		nb = 2+6*(k&1); for (int j,i=0;i<nb;i++) if ((j=xys[i]-48)<0) return BXE_PARSE; else buf[i]=j;
		k = 64*buf[1]+buf[0]+64;
	} else {
		if (k&65536) { WrapCore *cr = core_ro(); k &= 65535;
			       if (!cr->has_key(k)) return BXE_UNDEFKEY; else k = cr->get_key(k); }
		if (k<64) { nb = 0; if ((m_bflg&3)==1) tlim = min_i(2*sample_rate, m_node->etc()->i[0]); }
		else { nb = 2, buf[0] = k&63, buf[1] = (k>>6)-1; }
	}
	if ((unsigned int)(op-8)<3u){ if (m_mxctl && (ec=mx_c_stop(m_mxctl,k,1+(op==8)))!=MXE_CTLU) return ec;
				      if (op==10) ++op; else return EEE_NOEFF; }
	if (op!=3 && op!=6) memcpy(xysav, m_xys6, 8);
	else if (!(nof&NOF_FORCE) && !m_node->perm(DF_EDBOX)) return NDE_PERM;
	if (nb) { if (op!=6) m_bflg |= WRF_NOCON; for (int i=0; i<nb; i++) m_xys6[i] = buf[i]; }
	switch(op) {
		case  2: ec = qcopy(0, nof); break;
		case  3:
		case  6: return nb ? (w_gr_xyk(m_xys6[0], m_xys6[1], 1), w_mini(), w_col0(), 0) : EEE_NOEFF;
		case  5: ec = qcopy(1, nof & ~NOF_FGUI); break;
		case 11: ec = add2mx_txdlv(0,0,0,tlim,0); if (ec>=0) ec = add2ctl(ec, k);       break;
		case 12: ec = add2mx_txdlv(0,0,0,tlim,0); if (ec>=0) ec = add2ctl(ec, k|65536); break;
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

int AWrapGen::lim8(const char * q) { 
	int i,r=0; for (i=0; i<8; i++) if (m_xys6[i]>=q[i]) m_xys6[i]=q[i]-1, ++r;  return r; }

int AWrapGen::mini(char *to) {
	get_nm2(to); for (int i=0; i<4; i++) to[i+2] = i2aA(m_xys6[i]);
	to[6] = i2aA(m_xys6[7]); to[7] = ':' + !!(m_bflg&WRF_SHADOW);
	memcpy(to+8, v_rgb(), 6); return 14;
}

void AWrapGen::w_a20(int flg) { int k; gui2.setwin(w_oid(),'w'); 
	if (flg&(k=WRF_W_CH2 )) gui2.wupd_i1('Y', !!(m_bflg&k),  10);
	if (flg&(k=WRF_W_PLOT)) gui2.wupd_i1('Y', !!(m_bflg&k),  11); }

void AWrapGen::w_mini() {
	if (wnfl()) gui2.setwin(w_oid(), 'w'), gui2.wupd_c0('w', 's'), gui2.bxmini(this);
	ABoxNode *nd = m_node; ANode *up = m_node->up();
	switch(up->cl_id()) {
		case 'C': if (up->winflg(8)) static_cast<ClipNode*>(up)->draw_1(nd);   return;
		case 't': if (up->winflg(2048)) trk_cond_pm(static_cast<ABoxNode*>(up)->box(), nd, '+');return;
		default:  return;
	}}

// f: 1-sel 2-on 4-off 8-uniq
void AWrapGen::w_gr_xyk(int x0, int y0, int f) {
	for (int g9 = m_node->gui9(), i=1, x=x0+48, y=y0+48; i<9; i+=7) { if (g9 & i) {
		int c = 0x91623ee >> i; gui2.setwin(g9|((c>>16)&15), c&127); gui2.wupd('#');
		if (f&1) gui2.c3('b', x, y);
		if (f&6) gui2.c3('K'+16*(f&2), x, y);
		if (f&8) gui2.c3('u', (f&4)?33:x, y);
	}}}

void AWrapGen::w_tab0(int f10) { gui2.setwin(w_oid(), 'w'); const char *s = core_ro()->grdim;
	BVFOR_JM(f10&1023) gui2.wupd_c48('Y', s[j], j); }

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

int AWrapGen::batch_calc(double *to0, double *to1, int skip, int n, int nch) {
	Clock clk; if (DBGC) clk.reset();
	int nzch = 0, mxid = mx1(0); if (mxid<0) return mxid;
	while (skip>0) mx_calc(mxid, junkbuf, junkbuf, min_i(skip, 1024), 0), skip -= 1024;
	for (int n1, j=0; (n1=min_i(n-j,1024)) > 0; j+=n1) { 
		int k, r = mx_calc(mxid, to0+j, to1?to1+j:0, n1, nzch);
		if ((unsigned int)r > 2u) return (r<0) ? r : (log("BUG: wr/batch: r=%d",r), BXE_WTF);
		if (r>nzch) { switch(k=2*r+nzch, nzch=r, k) { 
			case 4: memset(to1,   0, 8*j);
			case 2: memset(to0,   0, 8*j); break;
			case 5: memcpy(to1, to0, 8*j); break;
			default: return log("BUG: wr/batch: k=%d",k), BXE_WTF; }}}
	mx_del(mxid); if (DBGC) log("batch_mono: %d", clk.get()); return nzch;
}

int AWrapGen::plot_t(double t0, double t1, int n, int flg) {
	if (DBGC) log("plot_t: %.15g %.15g %d f:%d", t0, t1, n, flg);
        int i0 = sec2samp(t0), i1 = sec2samp(t1), len = i1 - i0, f7 = flg&7;
        if (len<1) return log("plot_t: length = %d samples, sorry.", len), BXE_RANGE;
	double res[f7?2*len:len];
	int k, r = batch_calc(res, f7?res+len:0, i0, len, 0);
	if (r<=0) return r ? r : (qstat.store(zeroblkD, 1), BXE_ZPLOT);
	if (f7 && r==1) log("plot_t: sound is centered, drawing mono..."), flg = f7 = 0;
        if (len <= n) {
		qstat.store(res, f7?2*len:len);
                PlotPar_arr par0(res, len, t0, t1), par1(res+len, len, t0, t1);
		if (f7<3) k=1, Gnuplot::sg()->setfun1(0, arrfun1, f7==2 ? &par1:&par0, 0, "out\0L\0  R"+4*f7);
		else      k=3, Gnuplot::sg()->setfun1(0, arrfun1, &par0, 0, "L"),
			       Gnuplot::sg()->setfun1(1, arrfun1, &par1, 0, "R");
		return Gnuplot::sg()->plot1(k, t0, t1, len), 0;
        } else {
		int nst = (0x62333>>(4*f7))&7; if (!nst) return BXE_PARSE;
		PlotPar_arr par[nst]; double st[n*nst]; const char * spt;
		if (f7<3) samp_stat(f7==2?res+len:res, len, n, 0, 0.0, st+n,st+2*n,st),
			  spt="avg\0min\0max\0agL\0mnL\0mxL\0agR\0mnR\0mxR\0"+12*f7;
		else if (f7==3) samp_stat(res,     len, n, 0, 0.0, 0,0, st), spt="agL\0agR",
				samp_stat(res+len, len, n, 0, 0.0, 0,0, st+n);
		else samp_stat(res,     len, n, 0, 0.0, st+2*n, st+4*n, st),
		     samp_stat(res+len, len, n, 0, 0.0, st+3*n, st+5*n, st+n),
		     spt="agL\0agR\0mnL\0mnR\0mxL\0mxR";
		if (!f7) qstat.store(st+n, n); else if (f7==3) qstat.store(st, 2*n);
		for (int i=0; i<nst; i++) par[i].s(st+i*n, n, t0, t1),
				          Gnuplot::sg()->setfun1(i, arrfun1, par+i, 0, spt+4*i);
		return Gnuplot::sg()->plot1((1<<nst)-1, t0, t1, n), 0;
	}}

int AWrapGen::plot_f(double t0, double t1, double f0, double f1, int n, int flg) {
	if (DBGC) log("plot_f: %.15g %.15g %.15g %.15g %d %d", t0, t1, f0, f1, n, flg);
        int i0 = sec2samp(t0), i1 = sec2samp(t1), len = i1 - i0, f3 = flg&3;
        if (len<1) return log("plot_f: length = %d samples, sorry.", len), BXE_RANGE;
        int siz = 64, bits = 6;
        if (len>(1<<24)) { log("plot(F): len cut to %g", sample_length*(double)(1<<24));
			   siz = len = (1<<24), bits = 24; }
	else 		 { while (len>siz) siz+=siz,++bits;   if (!(flg&4)) len=siz; }
        double *buf = (double*)malloc(16*siz), *re, *im;
	int r = batch_calc(buf, f3?buf+siz:0, i0, len, 0);
	if (r<=0) return r ? r : (qstat.store(zeroblkD, 1), BXE_ZPLOT);
	if (f3&r&2) re=buf+siz, im=buf; else re=buf, im=buf+siz;
	memset(im, 0, 8*siz); if (len<siz) memset(re+len, 0, 8*(siz-len));
	double *res = fft(re, im, bits, false), *res2 = res + siz;  free(buf);
	int fmx = -1; double fmv = -1.0;
        for (int i=0; i<=siz/2; i++)
                if ((res[i] = sqrt(res[i]*res[i] + res2[i]*res2[i]))>fmv) fmv = res[i], fmx = i;
	log("plot_f: (max@Hz) %s%.9g", DEBUG_UTXT(15), (double)fmx / (double)siz * (double)sample_rate);
        int ix0 = (int)round(f0 * sample_length * (double)(siz));
        int ix1 = (int)round(f1 * sample_length * (double)(siz));
        if (ix1<siz-1) ++ix1;
        int n0 = ix1 - ix0;
        if (n > n0) n = n0;
	if (n<2 || n0<2) return log("plot_f: n=%d, n0=%d", n, n0), BXE_RANGE;
        bool statflg = (n < n0);
        double stat[ (statflg ? 3 : 2) * n ];
        samp_stat(res+ix0, n0, n, false, 0.0, 0, 0, stat);
	qstat.store(stat, n);
        samp_stat(res+ix0, n0, n, true, 55.0, 0, statflg ? stat+2*n : 0, stat+n);
        PlotPar_arr p_avg(stat, n, f0, f1);
        PlotPar_arr p_lavg(stat+n, n, f0, f1);
        PlotPar_arr p_lmax(stat+2*n, n, f0, f1);
        Gnuplot::sg()->setfun1(0, arrfun1, &p_avg, 0, "avg\0lft\0rgt\0BUG\0A/z\0L/z\0R/z\0BUG"+4*flg);
        Gnuplot::sg()->setfun1(1, arrfun1, &p_lavg, 1, "lavg");
        if (statflg) Gnuplot::sg()->setfun1(2, arrfun1, &p_lmax, 1, "lmax");
        Gnuplot::sg()->plot1(statflg ? 7 : 3, f0, f1, n);

        delete[] (res); return 0;
}

#define CH(X) BXCMD_H(AWrapGen, X)

CH(gcl){return (s[1]=='.') ? (s[2]==51 ? p->m_node->draw_window(25) 
			               : p->key_op(-2, s[2]-48, s+3, cb->cnof()))
			   : p->key_op(-1, s[1]-48, s+2, cb->cnof()); }

CH(gmd){p->m_bflg&=~3, p->m_bflg|=(s[1]&3); if (p->wnfl()) p->w_gmd(); return 0; }
CH(c2k){return min_i(0, p->qcopy(0, cb->cnof())); } 
CH(stp){return p->m_mxctl ? mx_c_stop(p->m_mxctl, 0, s[1]&3) : EEE_NOEFF; }
CH(xfd){return p->core_rw()->xfd = wr_ixtr(atoi_h(s+1)), 0; }
CH(flg){return unpkflg(&p->m_bflg, atoi_h(s+1), WR_FLG_MV(p)), p->m_bflg&=~WRF_B_RSRV, 0; }

CH(tf){	float *q = p->core_rw()->tf01; int flg = s[1]; double x; s += 2;
	for (int i=0, j=1; i<4; i++, j+=j) if (flg&j) s+=parse_num(&x,s), q[i] = (float)x;
	return 0; }

CH(pl){	if (s[1]<58) return p->key_op(11, s[1]-48, 0, cb->cnof());
	if (s[1]=='P') return p->add2mx_txdlv(0,0,0,INT_MAX,0);
	const float * q = p->core_ro()->tf01; switch(s[1]) {
		case 'T': return p->plot_t(q[0], q[1], 		   512, s[2]&7);
		case 'F': return p->plot_f(q[0], q[1], q[2], q[3], 512, s[2]&7);
		default: return BXE_CENUM;
}}

CH(win){return (s[1]=='t') ? (p->show_tab(s[2]&7), 0) 
			   : p->wlg(p->m_node->wdat_raw(), s[1]&3, (8+(s[1]&4))<<20); }

CH(ky){ WrapCore * cr = p->core_rw();
	for (++s; s[0]>47 && s[1]>47 && s[2]>47 && s[3]>47; s+=4)
		cr->set_key(hex2(s)&127, 64*s[2] + s[3] - 3120);   return 0; }

CH(wav){int j, k = 0, wf = p->wnfl(), *q = &p->m_bflg;
	switch(s[1]) {
		case 'F': p->m_bflg &= ~(k=WRF_W_ALL); p->m_bflg |= (((s[1]-48)<<WRF_W_SH)&k); goto flgw;
		case '2': k = WRF_W_CH2; goto flgs;
		case '-': k = WRF_W_PLOT; goto flgs;
		case 'W': case 0: return p->write_a20();
		case 'S': return k=WRF_MKSH,  (  s[2]&1) ? (*q|=k) : (*q&=~k), wf ? (p->w_tlim(4), 0) : 0;
		case 'U': return k=WRF_SLUPD, (j=s[2]&1) ? (*q|=k) : (*q&=~k), p->tab_vis(2)&&(p->w_tlim(8),0),
			         p->wnfl(512)  ? (gui2.setwin(p->w_oid(9),'#'), gui2.t_sn("Uu"+j,1), 0) : 0;
		default: return BXE_CENUM;
	}
flgs:	if (s[2]&1) *q |= k; else *q &= ~k;
flgw:	if (wf) p->w_a20(k);
	return 0;
}

CH(cut){int *q = p->m_node->etc()->i; double x = 0.0;
	switch(s[1]) {
		case '_': (s[2]&1) ? (*q|=INT_MIN) : (*q&=INT_MAX); break;
		case 'T': *q &= INT_MIN; parse_num(&x, s+2);
			  *q |= (x<-1e-5 ? INT_MAX : (INT_MAX & (int)lround(x*40320.0))); return 0;
		default : *q = atoi_h(s+1); break;
	}
	if (p->m_node->winflg(2048)) p->w_tlim(1);    return 0;
}

CH(gr){	int k, i = s[2] & 7, g9 = p->m_node->gui9(), f10 = 0;
	WrapCore * cr = p->core_rw();
	switch (s[1]) {
		case 'd': intv_cmd_c(cr->grdim+i, s+3, 2, 51, 0x33330801); f10 |= 1<<i; break;
		case 'R': i&=1; cr->grdim[i+8] = (s[3]-48) & 127;          f10 |= 256<<i;  break;
		case '*':
			  for (i=0; i<8 && (unsigned int)(k=s[i+2]-50)<50u; i++) cr->grdim[i] = k+2;
			  if (i==8 && s[10]) cr->grdim[8] = (s[10]-48) & 127,
				  	     cr->grdim[9] = (s[11]-48) & 63, i+=2;
			  f10 = 1023; break;
		default: return BXE_CENUM;
	}
	if (p->lim8(cr->grdim) || (f10&128)) p->w_mini();
	if (!(g9&9)) return 0;
	if (f10&255) p->w_col0();
	if (f10&771) cr->upd_grid(p->m_xys6, g9);
	if (p->tab_vis(0)) p->w_tab0(f10);
	if ((g9&8) && (f10>>=2, f10&=63)) gui2.setwin(g9|9, '#'), gui2.t0(), cr->slcmd(p->m_xys6, f10);
	return 0;
}

CH(rvt){ANode * nd = cb->lookup(s+1); if (!nd) return BXE_ARGLU;
	if (nd->cl_id()!='s') return NDE_EXPWRAPS; 
	return cb->perm(nd, DF_EDBOX) ? STC_BOX(nd, SWrap)->set_trg(p) : NDE_PERM; }

#define AW_CTAB {'+'+256, (cmd_t)c_c2k}, {'P'+256, (cmd_t)c_pl }, {'t', (cmd_t)c_tf}, {'x', (cmd_t)c_xfd}, \
	        {'W'+256, (cmd_t)c_win}, {'A'+256, (cmd_t)c_wav}, {'#', (cmd_t)c_gr}, {'T', (cmd_t)c_cut}, \
	        {'.'+256, (cmd_t)c_stp}, {'G'+256, (cmd_t)c_gcl}, {'k', (cmd_t)c_ky}, {'m', (cmd_t)c_gmd}, \
		{'<'+256, (cmd_t)c_rvt},                                              {'F', (cmd_t)c_flg}

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
	if (bxs==box_bookmark[1] || bxs==box_bookmark[3]) return bxf->node()->get_nm2(to);
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

#define TXD_GET(T,C) BVFOR_JM(flg&255) x = vs[nv++], T[j] = (flg&(256<<j)) ? \
	(x=((flg&WRF_IADJT)?x+(double)m_xys6[j]:x)*div1tab[(int)C->grdim[j]], x<0.0?0.0:(x>1.0?1.0:x)) : x
int DWrapGen::add2mx_txdlv(int trg, int flg, int dly, int lim, const double *vs) {
	if (!m_sfbx[0]) return EEE_NOEFF;
	double x, in[36], v11[11];
	int ni, no, nv = 0, ncf = (flg|=pflg()) & WRF_NOCON, v0f = !!(flg & WRF_SKIPV0),
	    trf = !(flg&WRF_NOREC) && !trg && trk_rec_trg;
	WrapSOB * sob = m_sob.ro(); 
	TXD_GET(v11, sob->m_core.ro());
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
	if (flg & WRF_MONO) in[5] = 0.0; // in[0] *= M_SQRT2;
	int ocb = no>1 ? 0x402 : 0x7c01; // TODO
	int r = mx_add_box(trg, m_sfbx[0]->mk_box(), udio, in, ocb, dly, lim); 
	return (r<0 || !trf) ? r : trk_rec(m_node, r);
}

int DWrapGen::qdiff(DWrapGen * that) {
	return (this!=that) && (m_sob.ro() != that->m_sob.ro() || memcmp(m_xys6, that->m_xys6, 8)); }

int DWrapGen::xfd_chk(int ff, int tf) {
	int xfd, xfd1 = m_sob.ro()->m_core.ro()->xfd;
	if (xfd1==63 || (xfd1>>6)!=ff) return 0; else xfd = wr_ixtr_r(xfd1);
	return  (!(xfd&32) && ((m_bflg>>(2+10*ff))&31)<=(xfd&31)) ||
		(tf&&m_sob.ro()->m_scl[ff].ro()->item(xfd1&63,0)->uses_xf()) ?
			(SOB_RW(sob)->core_rw()->xfd=63) : 0;
}
	
int DWrapGen::set_sf_2(int ff, BoxGen * bx) {
	int ec, sh = 2 + (10 &- ff), ni, no, msk = ~(1023<<sh);
	if (bx) ni = bx->n_in(), no = bx->n_out();
	if (ff && bx && (!ni||no!=1)) return (bx==m_sfbx[1]) ? (set_boxp(m_sfbx+1,0), m_bflg&=msk, BXE_FILTRM) 
						             : BXE_FILTNIO;
	if ((ec = set_boxp(m_sfbx + ff, bx)) < 0) return ec;
	m_bflg &= msk; if (bx) m_bflg |= (32*no+ni)<<sh;
	return 0;
}

#define SFUPD(X) w_mini(); w_slbv(2); if (wf) { BXW_GET; X ; wlg(bxw_rawptr, ff+1, 0); if (cf) w_col0_d(1); }
int DWrapGen::set_sf(int ff, BoxGen * bx) {
	BoxGen ** ppbx = m_sfbx + (ff&=1);
	if (bx==this || bx==*ppbx) return EEE_NOEFF;
	int ec, cf, wf = wnfl(2048);
	if (bx && bx->node()->is_wrap()) {
		DWrapGen * that = dynamic_cast<DWrapGen*> (bx); if (!that) return BXE_SORRY; //TODO (?)
		if ((ec=set_sf_2(ff, that->m_sfbx[ff]))<0) return ec;
		sob_from(4+ff, that, 1);
		sob_from(2+ff, that, 1);  cf = xfd_chk(ff,1);  SFUPD(wlg(bxw_rawptr, 0, 0)); }
	else {  if ((ec=set_sf_2(ff, bx))<0) return ec; else cf = xfd_chk(ff,0);  SFUPD( ); }
	return 0;
}

void DWrapGen::notify_nio(BoxGen * bx) {
	if (DBGC) log("wr%d/notify_nio: %d", id(), bx->id());
	int ec = 0, ni, no, flg = (bx==m_sfbx[0]) + 2*(bx==m_sfbx[1]);
	if (flg&1) m_bflg &= ~0xffc, m_bflg |= 4*(32*bx->n_out()+bx->n_in()), xfd_chk(0,0);
	if (flg&2) ((no=bx->n_out())==1 && (ni=bx->n_in())) ? 
		(m_bflg&=~0x3ff000, m_bflg |= (32*no+ni)<<12) : (set_sf(1, 0), ec = BXE_FILTRM), xfd_chk(1,0);
	if (ec) gui_errq_add(ec, "wrap/nio"); if (!flg || !wnfl()) return;
	sthg * q = m_node->wdat_raw(); if (!q) return gui_errq_add(BXE_WTF, "wrap/nio/w");
	if (flg&1) wlg(q,1,0); if (flg&2) wlg(q,2,0); 
}


int DWrapGen::sob_from(int ix, BoxGen * bx0, int bxf) {
	int cl = bx0->node()->cl_id();
	if(ix==1) return ((cl|4)=='w') ? (SOB_RW(sob)->m_core.set(static_cast<AWrapGen*>(bx0)->core_ro()), 0)
				       : NDE_EXPWRAP;
	if (cl!='w') return NDE_EXPWRAPD;
 	DWrapGen * that = static_cast<DWrapGen*> (bx0);
	int ec, ff = ix&1;
	if (!ix) { if (bxf && ((ec=set_sf_2(0, that->m_sfbx[0]))<0 || (ec=set_sf_2(1, that->m_sfbx[1]))<0))
			return ec;
		   m_sob.from(that->m_sob); return 0; }
	WrapSOB *sob = SOB_RW(sob), *sob2 = that->m_sob.ro();
	switch (ix) {
		case 4: case 5: sob->m_con[ff].from(sob2->m_con[ff]); return 0;
		case 6: 	sob->m_avol   .from(sob2->m_avol);    return 0;
		case 2: case 3: if (bxf && (ec=set_sf_2(ff, that->m_sfbx[ff]))<0) return ec;
				sob->m_scl[ff].from(sob2->m_scl[ff]); return 0;
		default: return BXE_IDX;
	}}

int DWrapGen::start_job_3(JobQ::ent_t * ent, char * arg) {
	int ec; switch(ent->i4f) {
		case 1:
			ent->plttwwii = 0x6b775924;
			WrapAutoVol * av; av = (glob_flg&GLF_AVOLSHR) ? m_sob.ro()->m_avol.ro() 
								      : SOB_RW(sob)->avol_rw();
			if (arg && *arg && *arg!='w') {
				if ((ec = av->parse_dim(arg))<0) return ec;
			} else if (wnfl()) {
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

const char * DWrapGen::v_rgb() {
	if (m_sfbx[1]) return m_sfbx[1]->v_rgb();
	if (m_sfbx[0]) return m_sfbx[0]->v_rgb(); return "KKK%%%"; }

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
		case 4: case 5: w_col0_d(-1); return;
		case 6:	memcpy(WR_AVCONF, m_sob.ro()->m_avol.ro()->xy12rt(), 6); w_avol(63, WR_AVCONF); return;
		default:	log("WTF: invalid index for w_sob (%d)", ix); return;
	}}

void DWrapGen::w_col0_d(int f) {
	WrapSOB * sob = m_sob.ro();
	int gf = wlg_vis(255); if (gf<=0) return;
	int flg = 32 | WRF_CFILT | ((f&1) << WRF_XFILT_SH) | pflg();
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
	int wf = p->m_node->winflg(2560), r = sob->icmd(s+1, &p->m_bflg);
	if (!wf || r<=0) return r; else if (DBGC) log("icmd ret: 0x%x", r);
	if ((wf>>8) & r & 2) p->w_slbv(2);
	if (!(wf & 2048)) return r;
	int oid = p->w_oid(), pf = p->pflg(), col = r&255, j1 = (r>>8)&255, j2 = (r>>16)&255;
		sob->wl(oid, j1, 	       1, col|pf, p->m_xys6, 0);
	if (j2) sob->wl(oid, wr_ixtr_r(j2^63), 1, col|pf, p->m_xys6, 0);
	if (r&(1<<24)) p->w_col0_d(1);
	return 1;
}

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
		case ':': if (debug_flags&DFLG_VOLTAB) log("vt-conf: #%x", p->m_node->id());
			  ec = sob->avol_rw()->parse_dim(s+2); break;
		default:  return BXE_UCMD;
	}
	if (ec<0 || !p->wnfl()) return ec;
	BXW_GETP(p); if (WR_TAB==2) p->w_sob(6, 0);
	return ec;
}

BXCMD_DEF(DWrapGen) {    {8192+'\\', 0}, AW_CTAB , 
	{'U', c_ud}, {'b', c_sf}, {'V',c_vt},  {'i',c_in}, {'B'+256, c_bw}, 
	{'D', c_ud}, {'>', c_sf}, {'O', c_so}, {0, 0}    };

/////// box (sh) ////////////////////////////////////////////////////////////

SWrapGen::SWrapGen(ABoxNode * nd) : AWrapGen(nd) {
	m_ssob.set(SWrapSOB_default(0)); m_trg = 0; 
	m_core.set(WrapCore_default(0)); m_bflg |= WRF_SHADOW|12; delayed_clip_upd(); }

SWrapGen::SWrapGen(ABoxNode * nd, const SWrapGen * p) : AWrapGen(nd, p) {
	m_ssob.from(p->m_ssob); m_core.from(p->m_core); m_trg = p->m_trg; 
	upd_conn(); delayed_clip_upd(); }

SWrapGen::SWrapGen(ABoxNode * nd, AWrapGen * p, int flg) : AWrapGen(nd) {
	m_bflg |= WRF_SHADOW|12; m_ssob.set(SWrapSOB_default(0)); m_trg = p; Node::conn(nd, p->node());
	m_core.set(p->core_ro()); memcpy(m_xys6, p->m_xys6, 8); delayed_clip_upd(); }

int SWrapGen::set_trg(AWrapGen* aw) {
	if (aw==m_trg) return 0;
	int ec = set_boxp((BoxGen**)&m_trg, aw);
	return (ec<0) ? ec : (wnfl() ? (w_trg(),0) : 0);
}

int SWrapGen::sob_from(int ix, BoxGen * bx0) {
	int cl = bx0->node()->cl_id();
	if (!ix) return ((cl|4)=='w') ? (m_core.set(static_cast<AWrapGen*>(bx0)->core_ro()), 0) : NDE_EXPWRAP;
	if (cl!='s') return NDE_EXPWRAPS;
	SWrapGen * that = static_cast<SWrapGen*>(bx0);
	if (ix==1) return m_ssob.from(that->m_ssob), 0;
	//SWrapSOB *sob = SOB_RW(ssob), *sob2 = that->m_ssob.ro();
	switch (ix) {
		default: return BXE_IDX;
	}}

void SWrapGen::v8tr(double *to, const double *src, int flg) {
	WrapCore * core = m_core.ro();      const char *grd = core->grdim;
	SWrapSOB * ssob = m_ssob.ro();      int cof, k = 0, cp = flg&WRF_COMP, f8 = (flg>>2)&255;
	SWrapTab * tab  = ssob->m_tab.ro(); double iarg[8], *pcon;
	if (flg&WRF_NOCON) {			cof =          0; pcon =         0; } 
	else { DblVec *con = ssob->m_con.ro();  cof = con->bv[0]; pcon = con->p[0]; }
	int ia8f=0; BVFOR_JM(f8&tab->iif) { int i = tab->ix[j]; ia8f |= (1<<(i&15)) | (1<<(i>>4)); } ia8f>>=2;
	if (flg&WRF_IADJS) { BVFOR_JM(ia8f) { iarg[j] = src[j] * (double)(grd[j]-1) + (double)m_xys6[j]; }}
	else		   { BVFOR_JM(ia8f) { iarg[j] = src[j] * (double)(grd[j]-1); }}
	const double *pp[2]; pp[0] = src; pp[1] = iarg;
	BVFOR_JM(f8) to[cp ? k++ : j] = (cof&(1<<j)) ? pcon[j] : tab->v1(j, pp);
}

int SWrapGen::add2mx_txdlv(int trg, int flg, int dly, int lim, const double *vs) {
	if (!m_trg) return EEE_NOEFF;
	double x, o8[8], v8[8];
	WrapCore * core = m_core.ro(); 
	SWrapSOB * ssob = m_ssob.ro();
	int bf = m_bflg, nv = 0, oflg = bf&1020,
	    pf = bf & (WRF_PASS|WRF_IADJS|WRF_IADJT), iif = ssob->m_tab.ro()->iif;
	TXD_GET(v8, core); 
	m_core.ro()->i2v (v8, m_xys6, flg|pf);
	v8tr(o8, v8, pf|oflg|WRF_COMP);
	int r = m_trg->add2mx_txdlv(trg, WRF_NOCON|WRF_NOREC|pf|(iif<<8)|(oflg>>2), dly, lim, o8);
	return (r<0 || (flg&WRF_NOREC) || trg || !trk_rec_trg) ? r : trk_rec(m_node, r);
}

int SWrapGen::save_sob(SvArg *p) { switch(p->st) {
	case 16: return m_core.save(p);
	case 17: return m_ssob.save(p);
	case 18: return m_ssob.ro()->m_tab.save(p);
	case 19: return m_ssob.ro()->m_con.save(p);
	case 20: case 21: case 22: case 23: case 24: p->st = 3; return 1;
	default: log("swr: invalid SOB state %d, skipping", p->st); return p->st2=-1, 1;
}}

int SWrapGen::save2(SvArg * sv) {
	BXSV2_HEAD; CHKERR(save2_aw(sv)); if (!m_trg) return r;
	ABoxNode * tnd = m_trg->node();
	if (tnd->is_save_trg()) { Node::eg_f_mv_to0(m_node, tnd); }
	else { CHKERR(f->sn("X$>", 3)); CHKERR(m_trg->node()->sv_path(10)); }
	return r;
}

int SWrapGen::wlg(sthg * bxw_rawptr, int ix, int flg) {
	if (!bxw_rawptr) return NDE_NOWIN; if ((unsigned int)ix > 0u) return ix<0 ? BXE_RANGE : BXE_WLGDONE;
	int oid = w_oid(), gc = (0x605a5345>>8*ix) & 127,
	    of = (flg & WRF_SETOC) ? ((flg&WRF_OC) ? (WR_WLG|=(1<<ix),1) : (WR_WLG&=~(1<<ix),0))
		    		   : (WR_WLG >> ix) & 1;
	// SWrapSOB * ssob = m_ssob.ro();
	gui2.setwin(oid, 'w'); if (!of) gui2.wupd_0(gc, ".+1");
	switch (ix) {
		case 0:  if (of) jt_show(), gui2.t_sn("&",1), gui2.hexn(m_bflg>>2,2), w_jtab(m_bflg&1020);
			 w_trg2(); break;
		default: return BXE_WTF;
	}
	gui2.wupd(gc,4); gui2.sn("S<>$XW  S><$XW"+8*of, 6); gui2.c1(52 + ix - 4*of); return 0;
}

void SWrapGen::w_jtab(int flg) {
	SWrapTab * tab = m_ssob.ro()->m_tab.ro();
	DblVec   *  dv = m_ssob.ro()->m_con.ro();
	double v[8], *pcon = dv->p[0];
	cur_c0(v, flg); flg>>=2; gui2.setwin(w_oid(), 'w');
	int i=0, bv = (m_bflg&WRF_NOCON) ? 0 : dv->bv[0];
	BVFOR_JM(flg) tab->w_jline(j, (bv&(1<<j)) ? pcon[j] : v[i]), i++;
}

void SWrapGen::w_col0_s(int f8) {
	if (wlg_vis(1) <= 0) return;
	double v[8]; cur_c0(v, f8); f8>>=2;
	gui2.setwin(w_oid(), 'w'); gui2.t0(); gui2.c3('!', hexc1(f8>>4), hexc1(f8&15));
	BVFOR_JM(f8) gui2.hdbl(v[j]);
}

int SWrapGen::show_tab_2(sthg * bxw_rawptr, int i) { switch(i) {
	case 3: return gui2.setwin(w_oid(),'w'), gui2.wupd_i1('Y', !!(m_bflg&WRF_IADJS), 30),
					   	 gui2.wupd_i1('Y', !!(m_bflg&WRF_IADJT), 31), 0;
	default:return BXE_CENUM;
}}

void SWrapGen::spec_debug() {
	log("swrap: flg=0x%x", m_bflg);
}

#undef CH
#define CH(X) BXCMD_H(SWrapGen, X)

CH(trg){if (s[1]==48&&!s[2]) return p->set_trg(0);
	ANode * nd = cb->lookup(s+1); if (!nd) return BXE_ARGLU;
	return nd->is_wrap() ? p->set_trg(STC_BOX(nd, AWrap)) : NDE_EXPWRAP; }

CH(so){	if (!s[1]||!s[2]) return BXE_NOARG;
	ANode * nd = cb->lookup(s+2);
	int i = s[1]&7,  ec = nd ? p->sob_from(i, nd->box0()) : BXE_ARGLU;
	if (ec>=0) p->w_sob(i, 1), p->w_slbv(2); return ec; }

CH(vfl){if (s[1]==42&&s[2]) { p->m_bflg &= ~255; p->m_bflg |= hex2(s+2); if (p->wnfl()) p->w_vfl(); return 0; }
	int c, m, j = s[1]-48; if (j&~7) return BXE_PARSE; else m = 4<<j;
	if ((c=s[2]&1)) p->m_bflg |= m; else p->m_bflg &= ~m;
	if (p->wlg_vis(1)) {
		gui2.setwin(p->w_oid(),'w'); gui2.wupd_i1('E', c, 8+j); p->jt_show(); if (c) p->w_jtab(m); }
	return p->w_slbv(2), 0; }

CH(scf){int fj, fm, k;
	switch(s[1]){case 'S': fm = WRF_IADJS; fj = 30; break;
		     case 'T': fm = WRF_IADJT; fj = 31; break;
		     default: return BXE_CENUM; }
	if ((k=s[2]&1)) p->m_bflg |= fm; else p->m_bflg &= ~fm;
	if (p->tab_vis(3)) gui2.setwin(p->w_oid(),'w'), gui2.wupd_i1('Y', k, fj);
	return 0; }

CH(tab){int x, j = s[1]-48, k = 0;
	if (j&~7) return EEE_RANGE; SWrapSOB * sob = SOB_RWP(p,ssob);
	switch(s[2]){case 'v':  if (p->m_bflg&WRF_NOCON) p->m_bflg&=~WRF_NOCON, 
							 sob->m_con.set(DblVec_default(0));
			        return *sob->con_rw()->addp(j) = at0f(s+3), 0;
		     case 'c':	intv_cmd_sc(sob->tab_rw()->cc+j, s+3, -125, 125, 0x190501); k = 4; break;
		     case 'a':case 'b': sob->tab_rw()->set_ab(j, s[2]-97, at0f(s+3)); break;
		     case 'x':case 'y': if ((unsigned int)(x=s[3]-48)>9u) return BXE_RANGE;
					sob->tab_rw()->set_xy(j, k=s[2]-120, s[3]-48); k += 21; break;
		     case 'i':	  	sob->tab_rw()->set_iif(j, s[3]&1); k = 16; break;
		     case '*':case '+': return   sob->tab_rw()->set_all(j, s+2, cb->tok());
		     default : return BXE_CENUM; }
	if (p->wlg_vis(1)) {    gui2.setwin(p->w_oid(),'w'); p->w_col0_s(4<<j);
				if (k) sob->m_tab.ro()->w1(j, k&15); }
	if (k&16) p->w_slbv(2);
	return 0;
}

BXCMD_DEF(SWrapGen) { {8192+'\\', 0}, AW_CTAB, 
	{'>', c_trg}, {'O', c_so}, {'v', c_vfl}, {'j', c_tab}, {'c', c_scf}, {0, 0} };

///////////// export /////////////////////////////////////////////////////////////

#define WARG  (static_cast<AWrapGen*>(bx))
#define WARGD (static_cast<DWrapGen*>(bx))
int setbox_wrap(ABoxNode *nd, BoxGen* _) { return (new (ANode::a64()) DWrapGen(nd))->set_tl(), 3; }
int setbox_shwr(ABoxNode *nd, BoxGen* _) { return (new (ANode::a64()) SWrapGen(nd))->set_tl(), 3; }
int setbox_shtg(ABoxNode *nd, BoxGen * bx) { return (bx->node()->is_wrap()) ?
	((new (ANode::a64()) SWrapGen(nd,WARG,0))->set_tl(), 3) : NDE_EXPWRAP; }
int wrap_2mx_txdlv(BoxGen * bx, int trg, int xflg, int dly, int lim, double *v) { 
	return WARG -> add2mx_txdlv(trg, xflg, dly, lim, v); }
int wrap_nd_2mx(ABoxNode * bnd, int trg, double bpm, int dly) {
	BoxGen * bx = bnd->box(); if (!bx) return 0;
	int tf = bnd->etc()->i[0]; if(tf != INT_MAX) {
		if (DBGC) log("nd2mx: tf0=%d, bpm=%g", tf, bpm); 
		if(tf<0)return 0; else tf=(int)lround((double)tf*natural_bpm/bpm); }
	if (DBGC) log("nd2mx: tf=%d", tf); return  WARG -> add2mx_txdlv(trg, 0, dly, tf, 0); }
int wrap_qdiff(BoxGen * bx, BoxGen * b2) { return WARGD -> qdiff(static_cast<DWrapGen*>(b2)); } // TODO
AReader * wrap_avreader(BoxGen * bx, int cflg) { return WARGD -> avreader(cflg); } 
int wrap_dump_keytab(char * to, unsigned int * bv, short ** pk) {
	int v, n = 0; for (int i=0; i<4; i++) BVFOR_JM(bv[i])
		to[n] = hexc1(2*i+(j>>4)), to[n+1] = hexc1(j&15), v = pk[i][j],
		to[n+2] = 48+(v>>6), to[n+3] = 48+(v&63), n += 4;
	return n; }
int wrap_key_op(BoxGen * bx, int ky, int op, const char *s, int nof) { return WARG -> key_op(ky,op,s,nof); }
void wrap_set_trec(BoxGen * bx, int j) { WARG->m_trec = j; }
void wrap_init() { WrapScale::st_init(); DWrapGen::st_init(); SWrapGen::st_init(); }
