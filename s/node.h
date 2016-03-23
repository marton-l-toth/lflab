#ifndef __qwe_node_h__
#define __qwe_node_h__

#include <new>
#include "util.h"
#include "tab7.h"
#include "glob.h"
#include "job.h"
#include "aobuf.h"
#include "wrap.h"
#include "trk.h"
#include "box0.h"

#define NODE_NAME_LEN 22

#define DF_EDBOX 1
#define DF_EDWRP 2
#define DF_EDDIR 4
#define DF_RMBD 5
#define DF_SAVE 010
#define DF_HROOT 020
#define DF_RESERVED 040
#define DF_EDWRP_T 3
#define DF_DELBOX 5
#define DF_ALL 077

#define WF_JOBF 1
#define WF_LT   2
#define WF_RT   4
#define WF_CLIP 8
#define WF_OID  0xfe1e
#define WF_OIDW 0xfe18
#define WF_WI1  0x10000
#define WF_WI6  0x3f0000
#define WF_WI8  0xff0000
#define WF_JOB1 0x1000000
#define WF_JOB5 0x1f000000
#define WF_JOB6 0x1f000001
#define WF_SHL  0x20
#define WF_LSEL 0x40
#define WF_RSEL 0x80
#define WF_2SEL 0xc0
#define WF_XSEL 0xe0
#define WF_ROOT 0x100
#define WF_SVNC 0x20000000
#define WF_SVFR 0x40000000

#define NC_DIR   0
#define NC_CLIP  32
#define NC_RSRV  64
#define NC_NBOX  96
#define NC_ANONW 255
#define NC_TRKW  254
#define NC_COMP  253
#define NC_CLIPW 252

#define NOF_FORCE  0x40000000
#define NOF_NOIDX  0x20000000
#define NOF_STRICT 0x10000000
#define NOF_OVRRD  0x08000000
#define NOF_FGUI   0x04000000
#define NOF_PARSE  0x02000000
#define NOF_YES    0x01000000
#define NOF_ERRMSG 0x00800000
#define NOF_THROW  0x00400000
#define NOF_FLAGS  0x7fe00000
#define NOF_IDX    0x000fffff
#define NOF_BATCH  0x11400000
#define NOF_BATCHV 0x11c00000
#define NOF_BATCHF 0x11a00000
#define NOF_CLI    0x00800000
#define NOF_GUI    0x04800000
#define NOF_EXPEOF 0x00200000
#define SVF_FORCE 1
#define SVF_COPY  2
#define SVF_WRAP  64
#define SOBF_LIB   0x100000
#define CHKERR(x) if ((ec = (x))<0) return ec; else r += ec
#define ND_NAME(s,n) char s[24]; s[n->get_name(s)] = 0
#define CLIP_DRAW_F(x) if (!((x)->winflg()|WF_CLIP)) (x)->draw(); else
#define VFUN_I(nm,r) virtual int nm() { return (r); }

class ANode; class ADirNode; class NDirNode; class ClipNode; class ABoxNode;
class BoxEdge; class BoxGen; class CmdBuf;

typedef int (sbfun_t)(ABoxNode *nd, BoxGen *bx),
	    (*sbfun )(ABoxNode *nd, BoxGen *bx);

struct SvArg {
        ANode *cn, *rn, *r2n, *wl;
        AOBuf * out;
        unsigned int vis; short st2; char st, flg;
	int nxup(int x);
};

struct trk_24 { char ty, ct; short i; int j; ANode * pv; int rsrv[2]; }; 
struct dir_24 { char ty, ct, n, s[21]; };
struct clp_24 { char ty, ct, i, rsrv[21]; };

class ANode {
        public: // static
                friend class Node; friend class NDirNode; friend class ClipNode;
                friend class ABoxNode; friend class LBoxNode; friend class TBoxNode;
                union u24_t { trk_24 t; clp_24 c; dir_24 d; char s[24]; unsigned char u[24]; };
                static void st_init();
                static char  * a64() { char *p = m0_f64 ? m0_f64 : blk64();
                                       m0_f64 = *(char**)p; return p; }
                static char  * aN () { ANode *p = m0_fN ? m0_fN : blkN();
                                       m0_fN = p->m_next; return (char*)p; }
                static void f64(void *p) { return;if ((unsigned long)p<1024) bug("f64-wtf??"); *(char**)p = m0_f64; m0_f64 = (char*)p; }
                static void f64c(void *p) { if (p) f64(p); }
                static void fN (ANode *p) { p->m_u24.d.ty = 0; p->m_next = m0_fN; m0_fN = p; }
                static char * z64() { return (char*)memset(a64(), 0, 64); }
                static char * cp64(const void *p) { return (char*)(p ? memcpy(a64(), (const char*)p, 64) : 0); }
                static inline ANode * root() { return (ANode*)m0_pnb[0]; }
                static ANode * lookup_n_q(int n) { return (ANode*)(m0_pnb[(n>>8)&4095]+128*(n&255)); }
                static ANode * lookup_n(int n) { ANode*p= (ANode*)(m0_pnb[(n>>8)&4095]+128*(n&255));
                                                 return p->m_u24.u[0] ? p : 0; }
                static void sv_dump1(const char *s);
		static void wi_clear();
		static void wi_debug();
        protected:
                typedef int (*nm_fun_t)(char*,const union u24_t*), (nm_fun_dt)(char*,const union u24_t*);
                static char * blk64();
                static ANode * blkN();
		static nm_fun_dt nm_0, nm_A, nm_C, nm_W, nm_b, nm_T,
		                 n2_0, n2_A, n2_C,             n2_T;
                static void del(ANode* p) { p->del2(); fN(p); }

                static void sv_start(AOBuf* out, ANode* r, int flg);
                static int  sv_write(int lim);
                static void sv_end();
		static void wi_init();

                static char *m0_f64;
                static ANode *m0_fN;
                static char *m0_pnb[4096];
                static int m0_nb_c;
                static int m0_lock;

                static nm_fun_t m0_nmfun[16];
                static ANode * m0_lr[2];
                static ANode * m0_glob_awlist;

                static int m0_wi_df, m0_wi_bf;
                static int m0_wi_d[256];
                static sthg m0_wi_b[256];

                static SvArg m0_sv;
        public: // this->
                virtual BoxGen * box0() { return 0; }
                virtual int set_perm(int mask, int perm) { return NDE_SETPERM; }
                virtual int perm_ed() = 0;
                virtual int perm_del() = 0;
                virtual int cond_del() = 0;
                virtual int save1() = 0; // ret: measure / error
                virtual ANode * sn(const char **pp) = 0;
                virtual ANode * sn_list(ANode ** pwl) = 0;
                virtual int perm(int flg = DF_ALL) = 0;
                virtual void debug(int flg) = 0;
                virtual int ccmd(CmdBuf * cb) = 0;
                virtual int wdat_alloc() = 0;
                virtual int wdat_free() = 0;
		virtual const char * rgb() { return "zz%z%%"; }
                virtual int start_job_2(JobQ::ent_t * ent, char * arg) { return JQE_UNDEF; }
                inline int sob_rw_arg() const { return (glob_flg&GLF_LIBMODE) + m_id; }
                int id() const { return m_id; }
                int cl_id() const { return m_u24.u[0]; }
                int get_name(char * to) const { return (*m0_nmfun[   m_u24.u[1]&7 ])(to, &m_u24); }
                int get_nm2 (char * to) const { return (*m0_nmfun[8+(m_u24.u[1]&7)])(to, &m_u24); }
                const char * s_name();
                ANode* up()   const { return m_up; }
                ANode* next() const { return m_next; }
                bool is_dir() const { return (m_u24.u[0] & 120)== 64; }
                bool is_box() const { return (m_u24.u[0] > 94); }
                int winflg(int m = -1) const { return m_winflg & m; }
                int gui9() const { return 16*m_id + ((m_winflg>>11)&1) + ((m_winflg>>6)&8); }
                void winflg_or(int x) { m_winflg |= x; }
                void winflg_and(int x) { m_winflg &= x; }
                int get_path(char* to, int max);
                int sv_path(int cz = -1);
                int get_path_uf(char* to, int max);
                int draw_window(int x);
                void close_window(int x);
                int title_arg(char * to, int wid = 512); // (id)rgbRGBpath
                int is_save_trg();
                int sv_cre(int flg = 0); // 1:no NL   2:force
                int start_job(int ix, char * arg, int flg = 0); // 1:force
                void set_job_result(int ix3, int res);
                int kill_job(int ix, int ec = JQE_KILL) { return jobq.kill(this, ix, ec); }
		const trk_24* cth() const { return &m_u24.t; }
        protected:
                ANode() : m_winflg(0), m_visitor(0), m_up(0), m_next(0) { }
                virtual int add(ANode * that, const char * nm, int i = NOF_PARSE, int j = -1) = 0;
                virtual int rm(ANode * that) = 0;
                virtual const char * pseudonym() { return "??"; }
                virtual void del2() = 0;
                virtual int draw_window_2(int x) { return NDE_NOWIN; }
                int get_path_2(char* to, int max);
                void set_name(const char * s) { int i=0; while (s[i]) m_u24.d.s[i]=s[i],i++; m_u24.d.n=i; }
                void a_debug();

                int m_id, m_winflg;
                unsigned int m_visitor, m_fill;
                ANode *m_up, *m_next;
		union u24_t m_u24;
};

// DC must def: int save2(SvArg * sv); 
// 		void debug2();
// static DC* mk(int);   ||  SOBDEF_[STD|64](DC)
// static void del(DC*); ||  DC(int); ~DC();
// DC* copy(int) const;  ||  DC(const DC*, int);

#define SOBDEF_STD(NM)  static NM * mk(int uarg) { return new NM(uarg); } \
                        static void del(NM * p) { delete(p); } \
                        NM * copy(int x) const { return new NM(this, x); }
#define SOBDEF_64(NM)   static NM * mk(int uarg) { return new (ANode::a64()) NM(uarg); } \
                        static void del(NM * p) { p->~NM(); ANode::f64(p); }  \
                        NM * copy(int x) const { return new (ANode::a64()) NM(this, x); }
#define SOB_RW(NM) (m_##NM.rw_2(sob_rw_arg()))
#define SOB_RWP(P,NM) ((P)->m_##NM.rw_2((P)->sob_rw_arg()))
#define SOB_RW_F0(NM,TY) inline TY * NM##_rw() { return SOB_RW(NM); }
#define SOB_RW_F1(NM,TY) inline TY * NM##_rw(int i) { return SOB_RW(NM[i]); }
#define SOB_INIFUN(TY,N) TY * TY##_default(int i) { \
        static SOB_p<TY> store[N]; \
        if ((unsigned int)i >= N##u) log("BUG: %s_default(%d) called (N:%d, ->0)", #TY, i, N), i=0; \
        TY * p = store[i].ro(); if (!p) (p=store[i].rw_2(1048576))->ini_default(i);   return p; }

class SOB {
        public:
                inline int sob_rw_arg() const { return (glob_flg&GLF_LIBMODE) + (m_u8_flg4_own&0xfffff); }
                int is_save_trg(SvArg * sv) {
                        return !(m_u8_flg4_own&SOBF_LIB) && m_visitor!=sv->vis; }
                void mark(int arg, unsigned int vis) { // ndid|SOBF_LIB
                        m_u8_flg4_own &= 0xffe00000; m_u8_flg4_own |= (unsigned int)arg; m_visitor = vis; }
                int   ref() { return ++m_u8_refcnt & 0xffffff; }
                int unref() { return --m_u8_refcnt & 0xffffff; }
                int shared() const { return (int)(m_u8_refcnt & 0xfffffe); }
                ANode * owner() { int i = m_u8_flg4_own&0xfffff; return i ? ANode::lookup_n(i) : 0; }
		void debug1() { log("SOB: p=0x%x, flg_own=0x%x, refcnt=%d", this, m_u8_flg4_own & 0xffffff,
									          m_u8_refcnt   & 0xffffff); }
        protected:
                SOB(int flg_own) : m_visitor(0u),
                        m_u8_flg4_own(flg_own), m_u8_refcnt(1u), m_u32(0u) {}
                SOB(const SOB* p, int flg_own) : m_visitor(0u),
                        m_u8_flg4_own((p->m_u8_flg4_own&0xffe00000)|flg_own),
                        m_u8_refcnt((p->m_u8_refcnt & 0xff000000)+1), m_u32(p->m_u32) {}
                ~SOB() {}
                unsigned int m_visitor, m_u8_flg4_own, m_u8_refcnt, m_u32;
};

inline static int sobref_h(SvArg * sv) { int x = 19*(sv->st>>4) + 0x304f2445 + ((sv->st&7)<<24);
					 return sv->out->sn((const char*)&x, 4); }

template <class T> class SOB_p {
        public: 
                SOB_p() : m_p(0) {}
                SOB_p(const SOB_p<T>& that) { (m_p = that.m_p)->ref(); }
                ~SOB_p() { unref_p(m_p); }
                inline T * ro() const { return m_p; }
                T * rw_2(int uarg) {
                        if (!m_p) return (m_p = T::mk(uarg));
                        if (!m_p->shared()) return m_p;
                        m_p->unref(); return (m_p = m_p->copy(uarg));
                }
		inline T * rw_o0() { return rw_2(1048576); }
                void set(T * p) { p!=m_p && (unref_p(m_p), m_p = p) && p->ref(); }
                void from(const SOB_p<T>& p) { set(p.ro()); }
                int save(SvArg * sv) {
                        if (!m_p) return sv->st2=-1, 1;
                        int ec, r = 1; ANode * xnd;
			if (sv->flg & SVF_COPY) {
				CHKERR(sobref_h(sv)); CHKERR(sv->out->pf("#%x\n",sv->cn->id())); goto qdone; }
                        if (sv->st2 || m_p->is_save_trg(sv)) {
                                CHKERR(m_p->save2(sv));
                                m_p->mark(sv->cn->id(), sv->vis); return r;
                        }
                        if (!(xnd = m_p->owner())) return sv->st2=-1, 1;
                        CHKERR(sobref_h(sv));   CHKERR(xnd->sv_path('\n'));
                 qdone: sv->st2 = -1; if (!(sv->st&7)) sv->st += 7;
                        return r;
                }
		void debug() { if (m_p) m_p->debug1(), m_p->debug2(); else log("SOB: NULL"); }
        protected:
                static void unref_p(T* p) { if (p && !p->unref()) T::del(p); }
                T * m_p;
};

class AReader { public: virtual int line(char * s) = 0; virtual ~AReader() {} }; // ret: 0:done <0:err

char * bigblk(int n);
char * alloc_32k();

class MiniDir {
        public: 
                MiniDir() { m_ent[0] = m_ent[1] = 0u; }
                int ty(int i) { return (int) (m_ent[i]>>25); }
                ANode * nd(int i) { return ANode::lookup_n((int)(m_ent[i]&0xfffff)); }

                int add1(int i, ANode * nd);
                int add2(int i, ANode * nd);
                int rm(int i);
        protected:
                unsigned int m_ent[16];
};

class BoxDesc : public SOB {
        public:
                SOBDEF_64(BoxDesc);
                BoxDesc(int arg) : SOB(arg) {}
                BoxDesc(const BoxDesc * that, int arg) : SOB(arg) { memcpy(m_u.ky, that->m_u.ky, 48);}
                ~BoxDesc() {}
		int dsc(char * to);
		union { char ** pp[4]; char ky[48]; } m_u;
};

class DblVec : public SOB {
        public:
                SOBDEF_64(DblVec);
                DblVec(int arg) : SOB(arg) { memset(bv, 0, 5); }
                DblVec(const DblVec * that, int arg);
                ~DblVec() { clear(); }
                int save2(SvArg * sv);
                void ini_default(int k);
                void clear();
                double * addp(int i) { int k = i>>3, j = i&7; if (!bv[k]) p[k] = (double*)ANode::z64();
                         bv[k] |= (1<<j); return p[k] + j; }
		void add_bv(int bv, const double * p) { BVFOR_JM(bv) *addp(j) = *(p++); }
                void rm(int i) { int k = i>>3; if (bv[k] && !(bv[k]&=255^(1<<(i&7)))) ANode::f64(p[k]); }
		void debug2(); 

                unsigned char bv[5];
                double * p[5];
};
DblVec * DblVec_default(int k);

struct NameVecAux {
        char * p(int ix) { int k = ix-14; return k<0 ? m_s14+4*ix : m_p+4*k; }
        char *pf(int ix) { int k = ix-14; return k<0 ? m_s14+4*ix : (m_p?m_p:(m_p=ANode::a64())) + 4*k; }
        char *m_p, m_s14[56];
};

class NameVec : public SOB {
        public:
                SOBDEF_64(NameVec);
                NameVec(int uarg) : SOB(uarg) { m_u32 = 0u; }
                NameVec(const NameVec * that, int uarg) : SOB(uarg) {
                        m_u32 = that->m_u32; memcpy(m_patt, that->m_patt, 8);
                        memcpy(m_s8, that->m_s8, 32); m_aux = (NameVecAux*)(ANode::cp64(that->m_aux)); }
                ~NameVec() { if (m_aux) { if (m_aux->m_p) ANode::f64(m_aux->m_p); ANode::f64(m_aux); } }
                int save2(SvArg * sv);
                void ini_default(int k);
                int get_nm(char *to, int ix);
                void set_nm(int bv, const char *s); // $-sep
                void set_nm_1(int ix, const char *s, int l = 4);
                int set_pt(int i, const char *s, int d);
		void qpt(int i, int j, const char *s) { s+=set_pt(0,s,-i); set_pt(1,j<99?s+(*s==36):"",-j); }
		int ls(char * to, int i0);
		void debug2(); 
        protected:
                char * p0(int i) { int j = i-8; return j<0 ? m_s8+4*i : m_aux->p(j); }
                NameVecAux * m_aux;
                char m_s8[32];
                signed char m_patt[8]; // EeeeOooo  E/O: ix(0)/ix(1)
};

class BoxUI : public SOB {
        public:
                SOBDEF_64(BoxUI);
                BoxUI(int uarg) : SOB(uarg) {}
                BoxUI(const BoxUI * that, int uarg);
		~BoxUI() {}
                int save2(SvArg * sv) { return sv->st2=-1, sv->out->pf("E$G%.6s\n", m_rgb); }
                void ini_default(int k);
		int cmd(const char *s);
                SOB_RW_F0(dv, DblVec) SOB_RW_F1(nm, NameVec)
		int draw_window_2(ANode * nd);
		void w_rgb(int oid);
		void debug2(); 

                SOB_p<DblVec> m_dv;
                SOB_p<NameVec> m_nm[2];
		SOB_p<BoxDesc> m_dsc;
                char m_rgb[6];
};
BoxUI * BoxUI_default(int k);

class ADirNode : public ANode {
	public:
		friend class Node;
		virtual int size() const = 0;
		virtual int gui_list(char *to, int flg) = 0; // 1: type (no 'x')
		virtual int perm(int flg = DF_ALL) { return perm_d(flg); }
		virtual int perm_ed() { return perm_d(DF_EDDIR); }
		virtual int perm_del() { return m_up && ((m_pm_val | ~m_pm_msk) & DF_EDDIR) 
						     && m_up->perm(DF_EDDIR); }
		virtual int set_perm(int mask, int perm) { m_pm_msk = mask&DF_ALL; m_pm_val = perm&m_pm_msk; return 0; }
		virtual int ccmd(CmdBuf * cb);
		virtual int save1(); // ret: measure / error
		virtual int wdat_alloc();
		virtual int wdat_free();
		virtual const char * rgb() { return "%%%ppp  %%pppp"+8*(m_u24.s[0]&1); }
		ADirNode * get_hroot() {
			for (ADirNode * p = this; 1; p=static_cast<ADirNode*>(p->up()))
				if (p->m_pm_msk&p->m_pm_val&DF_HROOT) return p; }
	protected:
		ADirNode() : m_pm_msk(0), m_pm_val(0) {}
		void ad_debug();
		int perm_d(int flg);
		unsigned char m_pm_msk, m_pm_val, m_siz;
};

class ABoxNode : public ANode {
	public:
		friend class Node;
		friend sbfun_t setbox_wrap, setbox_wrap_qcp, setbox_graph, setbox_calc, setbox_it;
		friend ANode * qmk_box(ANode * up, const char * nm, qmb_arg_t qa, int k, int ni, int no,
				                                    const char * cl, const char * fmt, ...);
		virtual int perm(int msk = DF_ALL) { return perm_b(msk); }
		virtual int perm_ed() { return perm_b(DF_EDBOX); }
		virtual int perm_del() { return m_up->is_dir() ? (m_up->perm(DF_RMBD)==DF_RMBD)
							       : perm_b(DF_EDBOX); }
		virtual int cond_del();
		virtual BoxGen * box0() { return m_box; }
		virtual const char * rgb();
		virtual int ccmd(CmdBuf * cb);
		virtual int save1(); // ret: measure / error
		virtual int wdat_alloc();
		virtual int wdat_free();
		virtual int start_job_2(JobQ::ent_t * ent, char * arg);

		BoxGen * box() const { return m_box; }
		void setbx_0(BoxGen * bx) { m_box = bx; }
		void unset_model_rec();
		const unsigned char * rgb_ro();
		void mwin_head(int wid);
		char * ionm(int io, int j, int f = 0);
		int ui_cmd(CmdBuf * cb);
		int find_fw(ABoxNode * to, LWArr<int>* rpath) {
			return this==to ? NDE_LOOPB1 : (to->m_et0 ? find_fw_2(to, rpath) : 0); }
		int find_fw_2(ABoxNode * to, LWArr<int>* rpath);
		void nio_change();
		sthg * wdat_raw() { int i = (m_winflg>>16)&255; return i ? m0_wi_b + 4*i : 0; }
		void set_ui_d(int k) { m_ui.set(BoxUI_default(k)); }
		void set_ui_f(ABoxNode * nd) { m_ui.from(nd->m_ui); }
		const char * own_rgb() { return m_ui.ro()->m_rgb; }
		int get_ionm(char *to, int io, int j);
		int dsc(char * to);
		void ab_debug(int flg);
		sthg * etc() { return &m_etc; }
	protected:
		ABoxNode(int ty) : m_box(0), m_ef0(0), m_efz(0), m_et0(0), m_etz(0) { m_u24.s[0] = ty; }
		virtual void del2();
		virtual int draw_window_2(int x);
		int save_rgb() { return 1; }
		int save_sob();
		int save_cfg();
		int perm_b(int msk) { ANode *p = m_up; while (!p->is_dir()) p = p->m_up; return p->perm(msk);}
		void qc_st8(int trg, int j);
		void qc_rgb(int ty, const char * s);
		void qc_dfv(int bv, double * dv);
		void qc_ptn(int trg, int i, int j, const char *s);
		void qc_iot(int trg, int bv, const char *s);

		unsigned char * rgb_rw();
		BoxGen * m_box;
		BoxEdge *m_ef0, *m_efz, *m_et0, *m_etz;
		sthg m_etc;
		SOB_p<BoxUI> m_ui;
}; 

struct BoxEdge {
        BoxEdge(ABoxNode *_fr, ABoxNode *_to) : cnt(1), fr_bx(_fr->box()), fr(_fr), to(_to),
                                              fr_p(0), to_p(0) {}
        int rsrv, cnt;
        BoxGen *fr_bx;
        ABoxNode *fr, *to;
        BoxEdge *fr_p, *fr_n, *to_p, *to_n;
};

class ClipNode : public ADirNode { // name: i_nnxy12
        public: 
		friend class Node;
                static ClipNode * kcp(int i) { return m0_kcp[i]; }  
		virtual int size() const { return bitcnt(m_map); }
		virtual int cond_del() { return !m_map ? 0 : NDE_NONEMP; }
		virtual int gui_list(char *to, int flg);
		virtual int ccmd(CmdBuf * cb) { return cmd(cb); }
		virtual ANode * sn(const char **pp) { int k = **pp=='*' ? m_sel : b32_to_i(**pp);
			return k<0 || !(m_map&(1u<<k)) ? 0 : (++*pp, ent_j(k)); }
		virtual ANode * sn_list(ANode ** pwl);
		virtual void debug(int flg);
		ABoxNode * ent_j(int j) { return (ABoxNode*)(m0_pnb[(int)m_eh[j]] + 128*m_el[j]); }
		ABoxNode * ent_jv(int j){ return (m_map&(1u<<j)) ? ent_j(j) : 0; }
		ABoxNode * ent_sel() { return ent_jv(m_sel); }
		BoxGen * bx_j(int j) { return (m_map&(1u<<j)) ? ent_j(j)->box() : 0; }
		BoxGen * bx_sel() { return bx_j(m_sel); }
		int keyop_j(int j, int ky, int op, const char *s, int nof);
		int keyop_f(int ky, int op, const char *s, int nof);
		int cmd(CmdBuf * cb);
		int add_wb(BoxGen * wb, int flg); //flg: 1:gui 2:sel 4:forcedup r: -1:full -2:!wb -3:nodup
		void del_wb(int i, bool guiflg);
		//void mov(int i, int j);
		int xchg(int i, int j);
		void draw();
		int sel(int k);
		void show_newbox(ABoxNode * nd);
		void draw_1(ABoxNode * nd);
        protected:
		static ClipNode * m0_kcp[3];
		static const char * wbname(BoxGen * wb, int ix); // p2static
                ClipNode() : m_map(0), m_sel(0), m_flg(3) {
			m_eh = (short*)a64(); m_u24.s[0] = 'C'; }
		virtual void del2() { f64((char*)m_eh); }
		virtual int add(ANode * that, const char * nm, int i = NOF_PARSE, int j = -1);
		virtual int rm(ANode * that);
		virtual int draw_window_2(int x) { return (x&&(x-3)) ? NDE_IW4 : (draw(), 16*'K'+3); }
                int wfind(BoxGen * bx);
		int find_free(int i0);

		unsigned int m_map;
		char m_sel,  m_flg; // 1:dup 2:au
		short * m_eh;
		unsigned char m_el[32];
};

class Node {
	public:
		static void init();
		static ADirNode * root() { return static_cast<ADirNode*> (ANode::root()); }
		static ANode * lr(int i) { return ANode::m0_lr[i]; }
		static ADirNode * lr_dir(int i) { ANode * p = ANode::m0_lr[i]; 
			return dynamic_cast<ADirNode*> (p->is_dir() ? p : p->m_up); }
		static void set_lr(int i, ANode * p);
		static ANode* lookup_cb (CmdBuf * cb, const char * path);
		static ANode* lookup_n(int id) { return ANode::lookup_n(id); }
		static int parse_target(CmdBuf * cb, char ** ppname, ANode** ppto);

		static int mk(ANode ** rr, ANode * up, const char * name, int ty, int i, int j = -1, BoxGen* from = 0);
		static int move(ANode * p, ANode * to, const char * name, int i, int j = -1);
		static int copy(ANode * p, ANode * to, const char * name, int i, int j = -1);
		static int del(ANode * p, int nof);

		static int conn(ABoxNode * fr, ABoxNode * to);
		static int disconn(ABoxNode * fr, ABoxNode * to);
		static int set_conn_2(ABoxNode * fr, ABoxNode * to1, ABoxNode * to2);
		static void eg_splice(BoxEdge * p);

		static void shl_add(ANode * nd);
		static void shl_rm(ANode * nd);
		static void slr_invd(int flg) { m0_slr_flg |= flg; }
		static int slr_flg() { return m0_slr_flg; }
		static int slr_upd_2(char * to);
		static int chkwin(int oid);
		static bool being_crawled() { return !!ANode::m0_sv.rn; }
		static void trk_chk_ord(ANode *tn, ANode *p, ANode *q, const char * msg);
		static ANode *trk_ij(ABoxNode * tn0, int i, int j);
		static ANode *trk_fwf(ANode *q), *trk_fwb(ANode *q);

		static int obj_help(int cl);
		static void lib_start(), lib_end();
		static int lib_cfg(ANode * nd);
		static int save_batch(ADirNode * dir, const char* fn, int flg);
	protected:
		static sbfun_t sb_btin, sb_trk;
		static BoxEdge * find_edge(ABoxNode * fr, ABoxNode * to);
		static int slr_gui_list(char * to, int ix);
		static void ini_echo(ANode * e0);
		static int conn2(ABoxNode * fr, ABoxNode * to);
		static int disconn2(ABoxNode * fr, ABoxNode * to, BoxEdge * p);
		static void contrib_init(ANode * r);

		static unsigned int m0_visitor;
		static int m0_slr_flg;
		static int m0_shl_n, m0_shl[32];
		static ADirNode * m0_clibroot;
	private:
		virtual void dummy() = 0;
};

#endif // __qwe_node_h__
