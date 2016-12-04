#include <ctype.h>
#include <limits.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "glob.h"
#include "pt.h"
#include "cmd.h"
#include "trk.h"
#include "guistub.h"
#include "util.h"
#include "util2.h"
#include "contrib.h"
#include "cfgtab.inc"

#ifdef LF_C_MEMDEBUG
int nd_count64;
int nd_mem_debug() {
	log("a64 total: %d", nd_count64);
	return 0;
}
#else
int nd_mem_debug() { log("compiled w/o mem_debug"); return NDE_SORRY; }
#endif

SOB_INIFUN(DblVec, 3)
SOB_INIFUN(NameVec, 6)
SOB_INIFUN(BoxUI, 5)

///////// decl ///////////////////////////////////////////////////////////////

class DscReader : public AReader {
	public:
		DscReader(BoxDesc * dsc, int n) : m_dsc(dsc), m_i(0), m_n(n) {}
		virtual int line(char * s);
	protected:
		BoxDesc * m_dsc;
		int m_i, m_n;
};

class NDirNode : public ADirNode {
        public: 
                friend class Node;
                virtual ANode * sn(const char **pp);
		virtual ANode * sn_list(ANode ** pwl);
                virtual int cond_del() { return ( (-m_siz)>>6 ) & NDE_NONEMP ; }
                virtual int gui_list(char *to, int flg);
                virtual int size() const { return m_siz; }
                virtual int start_job_2(JobQ::ent_t * ent, char * arg);
                virtual void debug(int flg);
        protected:
                static int hash2(const char * s, int sep);
                static inline int hash(const          char * s) { return hash2(             s, 0); }
                static inline int hash(const unsigned char * s) { return hash2((const char*)s, 0); }
                static inline int hash_dot(const char * s) { return hash2(s, '.'); }
                inline static int uicmp(const char *s, int h12, unsigned int trg) {
                        int r = h12 - (int)(trg>>20u); return r ? r : cmp_2(s, (int)trg); }
                static int cmp_2(const char *s, int trg);

                NDirNode() { m_siz = 0; m_u24.s[0] = 'D'; m_e[0] = m_small; m_e[1] = 0; }
                ANode * ent_i(int k) { return ANode::lookup_n_q(m_e[k>>4][k&15]); }
                virtual int add(ANode * that, const char * nm, int i = NOF_PARSE, int j = -1);
                virtual int rm(ANode * that);
                virtual void del2() {}
                inline int find(const char* name, int hash) { return m_siz ? find2(name, hash, 0, m_siz-1) : 64; }
                int find2(const char* name, int hash, int lo, int hi);
                int incr();
                void decr();
                void up1(int ix, int len);
                void dn1(int ix, int len);
                int find_sn(int id);
                int locmv(ANode * that, const char * nm);
                int add2(ANode * that, int h, const char * nm);

                unsigned int *m_e[2];
                unsigned int m_small[8];
};

class LBoxNode : public ABoxNode {
	public:
		friend class Node;
		virtual ANode * sn(const char **pp) { return 0; }
		virtual ANode * sn_list(ANode ** pwl) { return 0; }
		virtual void debug(int flg) { ab_debug(flg); }
	protected:
		LBoxNode(int ty) : ABoxNode(ty) {}
		virtual int add(ANode * that, const char * nm, int i = -1, int j = -1) { return NDE_LEAF; }
		virtual int rm(ANode * that) { return (that->m_up==this) ? NDE_WTF : NDE_RMWHAT; }
};

class TBoxNode : public ABoxNode {
	public:
		virtual ANode * sn(const char **pp);
		virtual ANode * sn_list(ANode ** pwl) {
			if (!pwl) return 0;
			ANode *g1 = trk_bkm_find(m_box, INT_MAX),
			      *w0 = trk_bkm_find(m_box, 0); 
			if (debug_flags&DFLG_TRK) log("tb/snl: w0=%p g1=%p pwl=%p", w0, g1, pwl);
			g1->m_next = *pwl; *pwl = w0; return 0;
		}
		virtual void debug(int flg);
		ANode * find_ijf(int i, int j, int flg, ANode * nd = 0);
		TBoxNode(int ty) : ABoxNode(ty) {}
		virtual int add(ANode * that, const char * nm, int i = -1, int j = -1);
		virtual int rm(ANode * that);
		void chk_ord(ANode * p, ANode *q, const char * s);
		void do_ins(ANode * q, int i, int j, ANode *nx);
		void do_cut(ANode * q);
		void ins_guard(ANode * that, int j, ANode * nx);
};

class AGuardNode : public ANode {
        public: 
                virtual int perm(int flg = DF_ALL) { return 0; }
                VFUN_I(perm_ed,0) VFUN_I(perm_del,0) VFUN_I(wdat_free,0)
                VFUN_I(wdat_alloc,NDE_GUARD)  VFUN_I(cond_del,0)
                virtual ANode * sn(const char **pp) { return 0; }
                virtual ANode * sn_list(ANode ** pwl) { return 0; }
                virtual int save1() { return bug("cannot save guard node"), 1; }
                virtual int ccmd(CmdBuf* cb) { return NDE_GUARD; }
                virtual int add(ANode*, const char*, int, int) { return NDE_GUARD; }
                virtual int rm(ANode*) { return NDE_GUARD; }
                virtual void del2() {}
};

class TGuardNode : public AGuardNode {
        public: 
                TGuardNode() { m_u24.t.ty = '!'; }
                TGuardNode(ANode * up, int ty) { m_up = up; m_u24.t.ty = '!'; m_u24.t.i = ty; m_next = 0; }
                virtual void debug(int i) { a_debug(); log("t.guard, i:0x%x, j:0x%x", m_u24.t.i, m_u24.t.j); }
};

///////// TODO: minidir //////////////////////////////////////////////////////

///////// dbl-vec/nm-vec/ui //////////////////////////////////////////////////

DblVec::DblVec(const DblVec * that, int arg) : SOB(arg) { memcpy(bv, that->bv, 5);
        for (int i=0; i<5; i++) if (bv[i]) p[i] = (double*)ANode::cp64(that->p[i]); }

void DblVec::clear() { for (int i=0; i<5; i++) if (bv[i]) bv[i]=0, ANode::f64(p[i]); }

void DblVec::debug2() {
	log_n("DblVec: [");
	for (int i=0; i<5; i++) {
		double * q = p[i]; int k = bv[i];
		for (int j=0; j<8; j++) if (k&(1<<j)) log_n(" %d:%.15g", 8*i+j, q[j]); }
	log(" ]");
}

int DblVec::save2(SvArg * sv) { switch(sv->st) {
	case 19: return save_sw(sv);
	case 20: return save_dw(sv, 0);
	case 21: return save_dw(sv, 1);
	default: return log("BUG: dblvec/sv.st=%d, skipped",sv->st), sv->st2=-1, 1;
}}

int DblVec::save_dw(SvArg * sv, int k) {
        int ec, r = 1, i0 = 64*k;
        for (int i=0; i<5; i++, i0+=8) {
                double * q = p[i];
		BVFOR_JM(bv[i]) { CHKERR(sv->out->pf("X$i%02x=%s\n", wr_ixtr_r(i0+j), dbl2str_s(0,q[j]))); }}
        return sv->st2=-1, 5*r;
}

int DblVec::save_sw(SvArg * sv) {
	int ec, r=1; double * q = p[0];
	BVFOR_JM(bv[0]) { CHKERR(sv->out->pf("X$j%cv%s\n", 48+j, dbl2str_s(0,q[j]))); }
	return sv->st2=-1, 5*r;
}

void DblVec::ini_default(int k) { switch(k) {
        case 1: *addp(1) = 1.0, *addp(3) = .5, *addp(4) = 1e-5; return;
        case 2: *addp(0) = .5, *addp(1) = 1e-5; return;
        default: return;
}}

void NameVec::debug2() {
	log_n("NameVec p0:%d,\"%3s\" p1:%d,\"%3s\" (", m_patt[0], m_patt+1, m_patt[4], m_patt+5);
	for (int i=0; i<32; i++) {
		unsigned int msk = 1u<<i;
		if (m_u32&msk) {
			char buf[8]; buf[get_nm(buf,i)] = 0; log_n("%s", buf);
		}
		if (!(m_u32 & ~(msk+msk-1u))) break;
		log_n("|");
	}
	log(")");
}

int NameVec::save2(SvArg * sv) { char buf[160]; buf[0] = 'E';
	int l = ls(buf+1,32*(sv->st&1)); buf[l+1] = 10; return sv->st2=-1, sv->out->sn(buf, l+2); }

int NameVec::ls(char * to, int i0) {
	int r = 0; BVFOR_JM(m_u32) {
		int k = i0 + (int)j;
		to[r] = 36; to[r+1] = 48+(k>>4); to[r+2] = hexc1(k&15); r += 3;
		const char * s = p0(j); for (int i=0; i<4 && s[i]; i++) to[r++] = s[i];
	}
	return r;
}

int NameVec::get_nm(char *to, int ix) {
        if (m_u32 & (1u << (ix&=31))) {
                const char * s = p0(ix);
                for (int i=0; 1; to[i] = s[i], i++) if (i==4 || !s[i]) return i;
        }
        int j; char * n3;
        if (!m_patt[5])    n3 = (char*)(m_patt+1)  , j = ix  + m_patt[0];
        else j = 4*(ix&1), n3 = (char*)(m_patt+j+1), j = (ix + m_patt[j])>>1;
        int tf = ((unsigned int)j>9u), i = 0;
        while (i+tf < 3 && n3[i]) to[i] = n3[i], ++i;
        if (tf) j<20 ? (j<0 ? (to[i++]='-', j = -j) : (to[i++]=49, j-=10))
                     : (tf = (j>29), to[i++] = 50+tf, j -= 10*(tf+2));
        to[i++] = 48+j; return i;
}

void NameVec::set_nm_1(int ix, const char *s, int l) {
	unsigned int msk = 1u << (ix &= 31);
        if (!l || !s || !*s) return (void) (m_u32 &= ~msk); else m_u32 |= msk;
        if ((unsigned int)l>4u) l = 4;
        int i = 0, k = (ix&=31) - 8;
        char * q = k<0 ? m_s8+4*ix : (m_aux?m_aux:(m_aux=(NameVecAux*)ANode::z64()))->pf(k);
        do { q[i] = s[i]; i++; } while (i<l && s[i] && s[i]!=36);
        if (i<4) q[i] = 0;
}

void NameVec::set_nm(int bv, const char *s) { while (1) {
        while (*s==36) ++s; if (!*s) return;
	int j = __builtin_ffs(bv) - 1; if (j<0) return;
        const char * s0 = s; while (++s, *s && *s!=36);
        set_nm_1(j, s0, s-s0); bv &= ~(1<<j);
}}

int NameVec::set_pt(int i, const char *s, int d) {
        i = 4*(i&1); m_patt[i] = d;
        int j=0; for (;j<3 && (m_patt[i+j+1] = (s[j]==36?0:s[j])); j++); return j; }

void NameVec::ini_default(int k) { switch (k) {
        case 0: set_pt(0, "in", 0); return;
        case 1: set_pt(0, "out", 0); return;
        case 2: set_pt(0, "x", 0); return;
        case 3: set_pt(0, "y", 0); return;
        case 4: set_pt(0, "s", -1); set_nm(3, "x$y"); return;
	case 5: set_nm(63, "BPM$from$to$rpt#$r.fr$r.to"); return;
}}

void BoxDesc::set_n(int k) {
	if (k==n) return;
	if (k>n) { do { int n2=n>>3, n3=n&7; if (!n3) pp[n2] = (char**)ANode::z64();
			pp[n2][n3] = (char*)ANode::a64(); ++n; } while(k>n); }
	else	 { do { --n; int n2=n>>3, n3=n&7; ANode::f64(pp[n2][n3]); pp[n2][n3]=0;
			if (!n3) ANode::f64(pp[n2]),pp[n2]=0;  } while(k<n); }}

int BoxDesc::save2(SvArg * sv) {
	if (!n) return log("BUG: saving empty boxdesc"), 1;
	char buf[3072]; int len = 5; memcpy(buf, "E$C0\n", 5); buf[3] += n;
	for (int i=0; i<n; i++) {
		const char * s = ln(i); int l1 = strlen(s);
		buf[len++] = '<'; memcpy(buf+len, s, l1); buf[len+l1] = 10; len+=l1+1;
	}
	return sv->st2=-1, sv->out->sn(buf, len);
}

int DscReader::line(char*s) {
	if (!m_i) m_dsc->set_n(m_n);
	int l = strlen(s); char *q = m_dsc->ln(m_i++);
	if (l<64) memcpy(q,s,l+1); else memcpy(q,s,63), q[63] = 0;
	return m_i < m_n;
}
AReader * ABoxNode::dsc_reader(int n) {
	BoxUI * ui = SOB_RW(ui); return new DscReader(SOB_RWP(ui,dsc), n); }

BoxUI::BoxUI(const BoxUI * that, int uarg) : SOB(uarg) {
        memcpy(m_rgb, that->m_rgb, 6);
        m_dv. from(that->m_dv ); m_nm[0].from(that->m_nm[0]);
        m_dsc.from(that->m_dsc); m_nm[1].from(that->m_nm[1]); }

int BoxUI::dump_dsc(char *to, int flg) {
	BoxDesc * dsc = m_dsc.ro(); 
	if (!dsc) return flg ? (memcpy(to, "1no description",15), 15) : 0;
	int n=dsc->n, r = flg?(*to=48+n,1):0;
	for (int i=0; i<n; i++) {
		if (i) to[r++] = '$';
		const char * s = dsc->ln(i); int l = strlen(s); memcpy(to+r,s,l); r+=l;
	}
	return r;
}

void BoxUI::ini_default(int k) {
        static const char * rgb = "uuu%%e  %%%uu0  Pz%%0%  zz%%z%  %%%fYu";
        memcpy(m_rgb, rgb+8*k, 6);
        m_dv.set(DblVec_default(0));
        m_nm[0].set(NameVec_default( (0x54200>>(4*k)) & 15 ));
        m_nm[1].set(NameVec_default( (0x11311>>(4*k)) & 15 ));
}

int BoxUI::cmd(const char *s) { int i; switch(s[0]) {
	case 'G': return (s[1]==33) || memcpy(m_rgb, s+1, 6);
	case '0': case '1': case '2': case '3':
		  i = *s - 48; if (!s[1]) return NDE_PARSE;
		  return SOB_RW(nm[i>>1])->set_nm_1(16*(i&1)+hxd2i(s[1]), s+2), 0;
	default:  return NDE_PARSE;
}}

void BoxUI::w_rgb(int oid) { if (oid) gui2.setwin(oid, 'C'), gui2.t0();
			    gui2.c1('R'); gui2.sn(m_rgb, 6); }

int BoxUI::draw_window_2(ANode * nd) {
	int oid = 16*nd->id() + 12;
	gui2.cre(oid, 'C'); w_rgb(0); gui2.c1('T'); gui2.nname(nd);
	gui2.namevec(m_nm[0].ro(), 0); 
	gui2.namevec(m_nm[1].ro(), 32);     return oid;
}

void BoxUI::debug2() {
	log("BoxUI: inm, onm, defv");
	m_nm[0].debug(); m_nm[1].debug(); m_dv.debug();
}

void ABoxNode::set_ui_rgb(const char * rgb) { memcpy(SOB_RW(ui)->m_rgb, rgb, 6); }
void ABoxNode::set_ui_lbl(int t, int i, const char * lbl) {
	SOB_RWP(SOB_RW(ui), nm[t])->set_nm_1(i, lbl, strlen(lbl)); }

///////// abs. node //////////////////////////////////////////////////////////

char * ANode::m0_f64;
ANode * ANode::m0_fN;
char * ANode::m0_pnb[4096];
int ANode::m0_nb_c = 0;
int ANode::m0_lock = 0;

int ANode::m0_wi_df, ANode::m0_wi_bf;
int ANode::m0_wi_d[256];
sthg ANode::m0_wi_b[256];

SvArg ANode::m0_sv;
ANode * ANode::m0_lr[2];
ANode::nm_fun_t ANode::m0_nmfun[16]={&nm_0, &nm_A, &nm_b, &nm_C, &nm_T, &nm_0, &nm_V, &nm_W,
				     &n2_0, &n2_A, &nm_b, &n2_C, &n2_T, &n2_0, &n2_V, &nm_W };
ANode * ANode::m0_glob_awlist = 0;

#define NMFUN(NM) int ANode::NM(char *to, const union ANode::u24_t * u)
NMFUN(nm_0) { *(int*)to = 0x21475542; return 4; }
NMFUN(n2_0) { to[0] = to[1] = '?'; return 2; }
NMFUN(nm_A) { int l = u->d.n; memcpy(to, u->d.s, l); return l; }
NMFUN(n2_A) { int l = u->d.n; const char * s = u->d.s; if(l>2 && (*s=='^'||*s=='='||*s=='~')) ++s,--l;
	to[0] = s[0], to[1]=l>1?s[l-1]:' '; return 2; }
NMFUN(nm_V) { *to = '?'; return 1; }
NMFUN(n2_V) { *to = '+'; to[1] = '?'; return 2; }
NMFUN(nm_C) { *to = i_to_b32(u->c.i); return 1; }
NMFUN(n2_C) { *to = '+'; to[1] = i_to_b32(u->c.i); return 2; }
NMFUN(nm_W) { *to = 'w'; to[1] = i_to_b32(u->c.i); return 2; }
NMFUN(nm_b) { int x=u->c.i; *to=48+(x>>4); to[1]=hexc1(x&15); return 2; }
NMFUN(nm_T) { *(int*)to = qh4(16*u->t.i)+0xa000000; *(int*)(to+4) = qh4(u->t.j>>16); 
					            *(int*)(to+8) = qh4(u->t.j&65535); return 12; }
NMFUN(n2_T) { to[0] = to[1] = '+'; return 2; }

const char * ANode::s_name() { static char s[24]; return m_u24.s[1]=='A' ? m_u24.d.s : (s[get_name(s)]=0,s); }
const char * ANode::path255(int l) { static char s[256]; s[get_path_uf(s,l&255)] = 0; return s; }
void ANode::close_window(int x) { int m=1<<(x&15); if ((m_winflg&m) && !((m_winflg&=~m)&WF_OID)) wdat_free(); }

int SvArg::nxup(int x) {
	ANode * q;      if (cn->winflg(WF_ROOT)) goto root;
	q = cn->next(); if (flg&SVF_WRAP) 	 goto skip_g;
	return q ? (cn = q, st = 0, x) : (cn = cn->up(), st = 1, x);
root: 	flg |= SVF_WRAP; q = wl;
skip_g: while (q && q->cl_id()<48) q = q->next(); 
	return (cn=q) ? (st = 0, x) : 0;
}

char * ANode::blk64() {
        char * p = alloc_32k(), *q = p, *q2, *qlim = p + 32704;
        do q2 = q + 64, *(char**)q = q2; while ((q=q2)<qlim);
        *(char**)q = 0; return p;
}

ANode * ANode::blkN() {
        char * p = m0_pnb[m0_nb_c] = alloc_32k(), *q = p;
        int cid = 256*(m0_nb_c++);
        ANode * nd;
        do nd = (ANode*)q, q += 128, nd->m_id = cid, nd->m_u24.s[0] = 0,
                nd->m_next = (ANode*)q; while(++cid & 255);
        nd -> m_next = 0; return (ANode*)p;
}

void ANode::wi_init() {
        for (int i=1; i<255; i++) m0_wi_d[i] = 1048576 + (i+1);
        m0_wi_d[255] = 1048576; m0_wi_df = 1; m0_wi_d[0] = 0;
        for (int i=1; i<63; i++) m0_wi_b[4*i].i[0] = 1048576 + (i+1);
        m0_wi_b[252].i[0] = 1048576; m0_wi_bf = 1;
}

void ANode::wi_clear() { const int msk = ~(WF_OID|WF_WI8);
	for (int k,i=1; i<255; i++) if ((k=m0_wi_d[i])<1048576) lookup_n_q(k)->winflg_and(msk);
	for (int k,i=1; i<64; i++) if ((k=m0_wi_b[4*i].i[0])<1048576) lookup_n_q(k)->winflg_and(msk);
	wi_init(); }

void ANode::wi_debug() {
	int dfc = 0, bfc = 0;  ANode * nd;
	log_n("dir:");
	for (int k,i=1; i<255; i++) if ((k=m0_wi_d[i])>1048575)  ++dfc; else
		nd = lookup_n_q(k), log_n(" 0x%x(%s)", k, nd->cl_id() ? nd->s_name() : "...0");
	log_n(" (free: %d)\nbox:", dfc);
	for (int k,i=1; i< 64; i++) if ((k=m0_wi_b[4*i].i[0])>1048575)  ++bfc; else
		nd = lookup_n_q(k), log_n(" 0x%x(%s)", k, nd->cl_id() ? nd->s_name() : "...0");
	log(" (free: %d)", bfc);
}

void ANode::st_init() {
        static int cnt = 0; if (++cnt>1) bug("WTF: ANode::st_init() called %d times", cnt);
        for (int i=0; i<4096; i++) m0_pnb[i] = zeroblkC;
	wi_init(); memset(&m0_sv, 0, sizeof(m0_sv));
}

int ANode::get_path(char* to, int max) {
        int r = get_path_2(to, max-1) & ~(1<<30);
        if (!r) to[r++] = '.';
        to[r] = 0;
        return r;
}

int ANode::get_path_uf(char * to, int max) {
	ANode *rt = root(), *nd = this;
	if (nd==rt) return *to='.', 1;
	if (max<44) return max = min_i(max,6), memcpy(to, "BUG44!", max), max;
	char nm[24]; int k = --max;
	while (1) {
		if (!nd) return memcpy(to, "BUG!!!", 6), 6;
		int n1 = nd->get_name(nm);
		if (n1>=k) { if (nd->m_up==rt) break; else if (k>5) memcpy(to+5, nm+(n1-k+5), k-5);
			     do nd=nd->m_up; while (nd->m_up!=rt);   break; }
		memcpy(to+(k-=n1), nm, n1), to[--k] = '.';
		if ((nd=nd->m_up)==rt) return k ? (memmove(to, to+k, max-=k), max) : max;
	}
	to[0] = '.'; memcpy(to+1+nd->get_name(to+1),"...",3); return max;
}

int ANode::get_path_2(char* to, int max) {
        if (!m_up) return 0;
        int r = m_up->get_path_2(to, max);
        return (r+22>max) ? r|(1<<30) : (to[r++]='.', r + get_name(to+r));
}

int ANode::draw_window(int x) {
        int k, x4 = x&15, nf = !(m_winflg & WF_OID);
        if (nf) wdat_alloc();
        if ((k = draw_window_2(x4))<0) {
                if (nf) wdat_free();
                return k;
        }
        x4 = k & 15; k >>= 4;
        int m1 = 1<<x4, prs = -(x>>4) & m_winflg & m1;
        if (prs) /*gui2.setwin(16*m_id + x4, k),*/ gui2.lx0('P');
        return m_winflg |= m1;
}

int ANode::title_arg(char * to, int wid) {
        char * p = to;
        *(p++) = '('; p += hx5(p, m_id); *(p++) = ')';
        memcpy(p, rgb(), 6); p += 6;
        return (p-to) + get_path(p, wid);
}

int ANode::is_save_trg() {
        if (m_visitor==m0_sv.vis) return 0; 
        if (!(m0_sv.flg&SVF_COPY)) return !!perm(DF_SAVE); 
        ANode * q; for (q = this; !q->winflg(WF_ROOT); q = q->m_up);
        return (q==m0_sv.rn);
}

void ANode::set_job_result(int ix3, int res) {
        res = ((res&1020)==1016) ? res&3 : 3;
        m_winflg &= ~(1|WF_JOB5);
        m_winflg |= (4*(ix3&7) + res) << 24;
}

int ANode::start_job(int ix, char * arg, int flg) {
        if (m_winflg&1) return JQE_DUP;
        ix &= 7; if (!(flg&1) && !perm_ed()) ix += 8;
        int ec = jobq.launch(this, ix, arg);
        if (ec>=0) m_winflg &= ~WF_JOB5, m_winflg |= 1 + ((ec&31)<<24);
        return ec;
}

int ANode::sv_cre(int flg) {
        if (!m_up) return NDE_NOROOT;
        if ((m0_sv.flg & SVF_COPY) && m0_sv.rn==this) return 1;
        AOBuf * out = m0_sv.out;
        unsigned int vis = m0_sv.vis;
        ANode *p, *pp[48];
        int r = 1, ec, nn = 0;
	if (flg & 2) {
		for (p=this; p->m_id; p=p->m_up)
			if (nn==48) return NDE_48; else pp[nn++] = p;
	} else {
		for (p=this; (p->m_visitor!=vis||!(p->m_winflg&WF_SVNC)) && !(p->m_winflg&WF_ROOT); p=p->m_up) {
			if ((debug_flags & DFLG_SAVE) && p->m_visitor==vis && !(p->m_winflg&WF_SVNC)) log("sv_cre/wflg: 0x%x", p->m_id);
			if (nn==48) return NDE_48; else pp[nn++] = p;
		}
	}
        if (!nn) return NDE_NOSVTRG;
        out->sn("C", 1); p->sv_path();
        char buf[24]; buf[0] = 36;
        for (int i=nn-1; i>=0; i--) {
		int ty = pp[i]->cl_id(); if (ty=='_') return NDE_BTCOPY;
                if ( ((buf[1] = ty)&120) != 64 && i ) return NDE_USVCONT;
                pp[i]->m_visitor = vis; pp[i]->m_winflg |= WF_SVNC; 
		if (debug_flags & DFLG_SAVE) log("svnc_flg set: 0x%x", pp[i]->m_id); 
                int l = pp[i]->get_name(buf+2); out->sn(buf, l+2); r += l + 2;
        }
        if (!(flg&1)) ++r, out->sn("\n", 1);
        return (ec = out->err()) ? ec : r;
}

int ANode::sv_path(int cz) {
        ANode *p=this, *pp[48];
        AOBuf * out = m0_sv.out;
        char buf[24];
        if (!(p->m_up)) return out->sn(".", 1);
        int ec, l, r = 1, nn = 0, cp = m0_sv.flg & SVF_COPY;
        if (cp) {
		for (p = this; !(p->m_winflg & WF_ROOT); p = p->m_up)
			if (nn==48) return NDE_48; else pp[nn++] = p;
                if (!p->m_up) return *buf='#', l = hx5(buf+1, m_id), buf[l+1]=cz,  out->sn(buf, l+1+(cz>=0));
                r += out->sn("=1", 2);
        } else {
		for (p = this; p->m_id; p = p->m_up) 
			if (nn==48) return NDE_48; else pp[nn++] = p;
	}
        *buf = '.';
        for (int i=nn-1; i>=0; i--) l = pp[i]->get_name(buf+1), r += out->sn(buf, l+1);
        if (cz>=0) r += out->sn((char*)&cz, 1);
        return (ec = out->err()) ? ec : r;
}

void ANode::sv_start(AOBuf * out, ANode * nd, int flg) {
        if (m0_sv.rn) { log("BUG: sv_start called twice!"); return; }
	ClipNode::kcp(0)->m_visitor = ++m0_sv.vis; 
	ClipNode::kcp(0)->m_winflg |= WF_SVNC; 
	m0_sv.flg = flg; m0_sv.out = out; m0_sv.wl = 0;
        (m0_sv.rn = m0_sv.cn = nd) -> winflg_or(WF_ROOT); nd->winflg_and(~WF_SVFR);
        m0_sv.st = 0;
}

int ANode::sv_write(int lim) {
        while(1) {
                int r = m0_sv.cn->save1();
                if (r<=0) return r;
                if (lim>0 && (lim-=r)<=0) return 1;
        }}

void ANode::sv_end() {
        if (m0_sv.out) delete(m0_sv.out), m0_sv.out = 0;
        ANode * nd = m0_sv.rn; if (nd && nd->id()) nd->winflg_and(~WF_ROOT);
        m0_sv.cn = m0_sv.rn = 0;
}

void ANode::sv_dump1(const char *s) {
        char buf[256]; buf[m0_sv.cn->get_path(buf, 255)] = 0;
        log("%s/save1: #%x(%c:%s), st=%d, st2=%d", s, m0_sv.cn->id(), m0_sv.cn->cl_id(),
                                  buf, m0_sv.st, m0_sv.st2); }

void ANode::a_debug() {
        log("p=%p id=%x, ty=0x%x(%c), ty2=0x%x(%c), up=0x%x(%p) name=\"%s\" wf=0x%x",
                        this, m_id, m_u24.s[0], m_u24.s[0], m_u24.s[1], m_u24.s[1],
                        m_up?m_up->m_id:-1, m_up, s_name(), m_winflg);
}

void nd0_init() { ANode::st_init(); }

///////// abs. dir ///////////////////////////////////////////////////////////

int ADirNode::perm_d(int msk) {
	msk &= DF_ALL; glob_flg &= ~GLF_SAVED;
	int m2=0, prm=0;
	for (ADirNode * p = this; m2!=msk; p=dynamic_cast<ADirNode*>(p->m_up))
		prm |= (p->m_pm_val&~m2), m2 |= (msk&p->m_pm_msk);
	return prm & msk;
}

int ADirNode::ccmd(CmdBuf * cb) {
	char * s = cb->a1(); ANode * nd;
	int ec = Node::mk(&nd, cb->cnode(), s+1, *s, cb->cnof()|NOF_PARSE);
	cb->set_curnode(ec<0 ? 0 : nd);
	return ec;
}

int ADirNode::save1() {
	if (debug_flags & DFLG_SAVE) sv_dump1("dir");
	char * pst = &m0_sv.st;
	if (*pst || (!perm(DF_SAVE) && !(m0_sv.flg&SVF_COPY))) return m0_sv.nxup(1);
	int ec = 0; ANode * snl = sn_list(&m0_sv.wl);
	if (snl) return m0_sv.cn = snl, ec + size() + 1;
	if ((ec = sv_cre()) < 0) return (ec==NDE_NOSVTRG) ? m0_sv.nxup(1) : ec;
	return m0_sv.nxup(ec+1);
}

void ADirNode::ad_debug() {
	a_debug();
	log("perm_mask=%o, perm_val=%o", m_pm_msk, m_pm_val); 
	char gl[1024]; 
	gl[gui_list(gl,0)] = 0; log("gui_list(0)={%s}", gl);
	gl[gui_list(gl,1)] = 0; log("gui_list(1)={%s}", gl);
}

int ADirNode::wdat_alloc() {
	if (m_winflg & WF_WI8) return 0;
	if (!m0_wi_df) return NDE_DIR255;
	int i = m0_wi_df; m0_wi_df = m0_wi_d[i] & 255; m0_wi_d[i] = m_id;
	m_winflg |= 65536 * i;
	return i;
}

int ADirNode::wdat_free() {
	int i = m_winflg & WF_WI8; m_winflg &= ~WF_WI8;
	if (!i) return 0; else i >>= 16;
	m0_wi_d[i] = m0_wi_df + 1048576; m0_wi_df = i;
	return 1;
}

///////// dir ////////////////////////////////////////////////////////////////

#define DIRLOOP2 for (int b=0; b<2; b++) if (winflg(2*b+2))
#define ADDL(J) (q = lookup_n_q(J), q->is_wrap() ? (q->m_next=wl, wl=q) : (q->m_next=r, r=q))
#define RMTHAT do { int ec; if (that->m_up && (ec=that->m_up->rm(that))<0) return ec; } while(0)

ANode * NDirNode::sn_list(ANode ** pwl) {
	int n = m_siz; if (!n) return 0;
	unsigned int * pe = m_e[0];
	ANode *q, *r = 0, *wl = pwl ? *pwl : 0;
	if (n>16) { for (int i=0; i<16; i++) ADDL(pe[i]);   pe = m_e[1]; n -= 16;
		    for (int i=0; i<n ; i++) ADDL(pe[i]);  }
	else {      for (int i=0; i<n ; i++) ADDL(pe[i]);  }
	if (pwl) *pwl = wl; return r;
}

int NDirNode::hash2(const char * s, int sep) {
        static const unsigned int cmap[8] = {0, ~0x4010u, ~0u, 0x7fffffff, 0,0,0,0};
        if (!s || !*s || *s==sep) return NDE_ZNAME;
        int c, i, h; 
	for (i=h=0, c=*s; i<21; h = (h<<5) + ((h>>7)&31) + c, c=s[++i])
		if (!(cmap[(unsigned char)c>>5] & (1U<<(c&31)))) 
			return sep ? ((c&&c!=sep) ? NDE_XNAME : (h&4095)+(i<<12))
				   : (c ? NDE_XNAME : (h&4095));
	return NDE_LNAME;
}

int NDirNode::cmp_2(const char *s, int to) {
	char * q = lookup_n_q(to)->m_u24.s + 2;
	int i, r, l = *(q++);
	for (i=0; i<=l; i++) if ((r = s[i]-q[i])) goto diff;
	return 0;
diff:	int sep = (s[i]==46);
	return (i==l) ? sep^1 : ((-sep)|r);
}

int NDirNode::incr() { 
	if (!(m_siz&8)) return (m_siz==32) ? NDE_FULL : 0;
	if (m_siz==8) { if (!m_e[1]) m_e[1]=m_e[0], memcpy(m_e[0]=(unsigned int*)a64(), m_small, 32); }
	else { 		if (m_e[1]==m_small) 	    memcpy(m_e[1]=(unsigned int*)a64(), m_small, 32); }
	return 0;
}

void NDirNode::decr() { switch(m_siz) {
	case 24: if (m_e[1]!=m_small) memcpy(m_small,m_e[1], 32), f64((char*)m_e[1]), m_e[1] = m_small; return;
	case 8:  if (m_e[1]) m_e[1]=0,memcpy(m_small,m_e[0], 32), f64((char*)m_e[0]), m_e[0] = m_small; return;
	default: return;
}}

int NDirNode::find2(const char* name, int h12, int lo, int hi) {
	int ofs=0, r;
	if ((lo^hi)&16) {
		if (!(r = uicmp(name, h12, m_e[0][15]))) return 15;
		if (r<0) hi=14, r=0; else lo=0, hi-=16, r=1, ofs=16;
	} else {
		ofs = (lo&16); r = lo>>4; lo&=15; hi&=15;
	}
	unsigned int * p = m_e[r];
	while(lo<=hi) {
		int md = (lo+hi)>>1;
		if (!(r = uicmp(name, h12, p[md]))) return md + ofs;
		if (r<0) hi = md-1; else lo = md+1;
	}
	return 64 + lo + ofs;
}

void NDirNode::up1(int ix, int len) {
	if (!len) return;
	if (! (((ix-1)^(ix+len-1))&16) ) {
		unsigned int * p = m_e[(ix>>4)&1] + (ix&15) - 1;
		for (int i=0; i<len; i++) p[i] = p[i+1]; return;
	}
	unsigned int *p = m_e[0], *q = m_e[1];
	for (int i=ix-1; i<15; i++) p[i] = p[i+1];
	p[15] = q[0]; len += ix-17;
	for (int i=0; i<len; i++) q[i] = q[i+1];
}

void NDirNode::dn1(int ix, int len) {
	if (!len) return;
	if (! ((ix^(ix+len))&16) ) {
		unsigned int * p = m_e[(ix>>4)&1] + (ix&15);
		for (int i=len; i>0; i--) p[i] = p[i-1]; return;
	}
	unsigned int *p = m_e[0], *q = m_e[1];
	for (int i=ix+len-16; i>0; i--) q[i] = q[i-1];
	q[0] = p[15];
	for (int i=15; i>ix; i--) p[i] = p[i-1];
}

ANode * NDirNode::sn(const char **pp) {
	const char * s = *pp;
	int h = hash_dot(s); if (h<0) return 0;
	int i = find(s, h&4095); if (i&64) return 0;
	*pp += (h>>12); return ent_i(i);
}

int NDirNode::find_sn(int id) {
	unsigned int uid = id, *p = m_e[0];
	if (m_siz<17) {
		for (int i=0; i<m_siz; i++) if ((p[i]&0xfffff) == uid) return i;
		return -1;
	}
	for (int i=0; i<16; i++) if ((p[i]&0xfffff) == uid) return i;
	p = m_e[1];
	for (int i=16; i<m_siz; i++) if ((p[i-16]&0xfffff) == uid) return i;
	return -1;
}

int NDirNode::add(ANode * that, const char * nm, int i, int j) {
	if (that->m_up == this) return locmv(that, nm);
	int r, h;
	char nmbuf[24];
	const char * nm2;
	if (!nm || !*nm) nmbuf[that->get_name(nmbuf)] = 0, h = hash(nm2 = nmbuf);
	else if ((h = hash(nm2 = nm)) < 0) return h;
	if ( !(m_siz&7) && (r=incr())<0) return r; 
	if ( (r=add2(that,h,nm2))!=NDE_NDUP || (i&NOF_STRICT) ) return r;
	if (nm2!=nmbuf) strcpy(nmbuf, nm2);
	int l = min_i(strlen(nmbuf), 18), c = '1';
	nmbuf[l] = '_'; nmbuf[l+2] = 0;
	while (1) {
		if (c=='9') c='a'; else if (c=='z') return NDE_WTF; else ++c;
		nmbuf[l+1] = c; h = hash(nmbuf);
		if ((r = add2(that, h, nmbuf)) != NDE_NDUP) return r;
	}
}

int NDirNode::locmv(ANode * that, const char * nm) {
	if (!nm || !*nm) return EEE_NOEFF;
	int /*h = hash(that->m_u24.d.s),*/ id = that->m_id, i0 = find_sn(id),
	    h2 = hash(nm);
//	log("locmv: \"%s\":%x -> \"%s\":%x", that->m_u24.d.s, h, nm, h2);
	if (h2<0) return h2;
	if (i0<0) return NDE_BASTARD;
	int r = uicmp(nm, h2, m_e[i0>>4][i0&15]);
	if (!r) return EEE_NOEFF;
	if (r<0) {
		r = i0 ? find2(nm, h2, 0, i0-1) : 64;
		if (r&64) r&=63; else return NDE_NDUP;
		dn1(r, i0-r); // log("locmv: dn %d %d", r, i0-r);
	} else {
		r = (i0<m_siz-1) ? find2(nm, h2, i0+1, m_siz-1) : 64+m_siz;
		if (r&64) r&=63; else return NDE_NDUP;
		up1(i0+1, r-i0-1); /* log("locmv: up %d %d", i0+1, r-i0-1); */ --r; --i0; 
	}
	m_e[r>>4][r&15] = ((unsigned int)h2<<20u) + (unsigned int)id;
	int l = strlen(nm); memcpy(that->m_u24.d.s, nm, l+1); that->m_u24.d.n = (char)l;
	DIRLOOP2 gui2.node_name(b, that);
	return r;
}
	
int NDirNode::add2(ANode * that, int h, const char * nm) {
	int l,ix = find(nm, h);
	if (!(ix&64)) return NDE_NDUP; else ix &= 63;
	RMTHAT; that->m_up = this; that->m_u24.d.ct = 'A';
	that->m_u24.d.n = l = strlen(nm); memcpy(that->m_u24.d.s, nm, l+1);
	dn1(ix, m_siz-ix); ++m_siz;
	m_e[ix>>4][ix&15] = ((unsigned int)h<<20u) + (unsigned int)(that->m_id);
	DIRLOOP2 gui2.node_name(b, that);
	return ix;
}

int NDirNode::rm(ANode * that) {
	if (that->m_up != this) return NDE_RMWHAT;
	int ix = find_sn(that->m_id); if (ix<0) return NDE_BASTARD;
	up1(ix+1, m_siz - ix - 1); 
	if (!(m_siz&7)) decr();
	--m_siz; 
	DIRLOOP2 gui2.node_rm(b, that);
	return 0;
}

int NDirNode::gui_list(char * to, int flg) {
	if (!m_siz) return 0;
	char * nm[m_siz]; int di = 0, bi = m_siz-1, ofs = (m_u24.d.s-(char*)this);
	for (int i=0; i<m_siz; i++) {
		char * s = (char*) ent_i(i);
		nm[ s[ofs-3]>94 ? bi-- : di++ ] = s;
	}
	ssort(nm, di, ofs); ssort(nm+di, m_siz-di, ofs);
	char * p = to;
	for (int i=0; i<m_siz; i++) {
		char * s = nm[i];
		if (i) *(p++) = 36;
		if (flg&1) *(p++) = s[ofs-3];
		else if (i==di) *(p++) = 'x';
		p += hx5(p, ((ANode*)s)->m_id); *(p++) = 36;
		memcpy(p, s+ofs, s[ofs-1]); p += s[ofs-1];
	}
	return p - to;
}

int NDirNode::start_job_2(JobQ::ent_t * ent, char * arg) {
	return NDE_SORRY;

}

void NDirNode::debug(int flg) {
	ad_debug();
	log(" nd_id hsh t node_name_size_20chr (siz=%d)", m_siz);
	for (int i=0; i<m_siz; i++) {
		unsigned int j = m_e[i>>4][i&15]; ANode * nd = ANode::lookup_n_q(j);
		log(" %05x %03x %c %s", j&0xfffff, j>>20u, nd->cl_id(), nd->s_name()); }
	log("");
}

///////// clipboard //////////////////////////////////////////////////////////

ClipNode * ClipNode::m0_kcp[3];

void ClipNode::del2() { 
	ClipNode **qq = m0_kcp, *dc = qq[0]; 
	int dci = dc->winflg(8) ? dc->id() : 0;
	for (int i=1; i<3; i++) if (qq[i]==this && (qq[i]=dc, dci)) gui2.clip_flg(dci,54+13*i,1);
	f64((char*)m_eh);
}

ANode * ClipNode::sn(const char **pp) { 
	int c = **pp, k = b32_to_i(c);
	if (k<0) { switch(c) {  case '*': k = m_sel; break;
				case '?': return m_extra ? (++*pp, lookup_n_q(m_extra)) : 0;
				default : return 0; }}
	return (m_map&(1u<<k)) ? (++*pp, ent_j(k)) : 0;
}

ANode * ClipNode::sn_list(ANode ** pwl) {
	ANode * r = m_extra ? lookup_n_q(m_extra) : 0;
	if (!pwl) return r;
	ANode *q, *wl = *pwl;
	BVFOR_JM(m_map) q = ent_j(j), q->m_next = wl, wl = q;
	*pwl = wl; return r;
}

int ClipNode::xchg(int i, int j) {
	if (i==j) return EEE_NOEFF;
	unsigned int t, m1 = 1u<<i, m2 = 1u<<j, m12 = m1|m2, m12a = m12 & m_map;
	if (!m12a) return EEE_NOEFF;
	ANode *p, *q;
	if (m12a==m12) { 
		t = m_eh[i], m_eh[i] = m_eh[j], m_eh[j] = t;
		t = m_el[i], m_el[i] = m_el[j], m_el[j] = t;
		(p = ent_j(i))->m_u24.c.i = i; 
		(q = ent_j(j))->m_u24.c.i = j;
		DIRLOOP2 gui2.node_name(b, p), gui2.node_name(b, q);
	} else {
		if ((m_map^=m12)&m1) i^=j, j^=i, i^=j;
		m_eh[j] = m_eh[i], m_el[j] = m_el[i];
		(p = ent_j(j))->m_u24.c.i = j;
		DIRLOOP2 gui2.node_name(b, p);
	}
	if (winflg(8)) gui2.clip_box(this, i, j, -1);
	return j;
}

int ClipNode::add_hlp(ANode * that) {
	if (m_extra) 	      return NDE_NDUP;     if (that->cl_id()!='h') return NDE_EXPHLP;
	if (that->m_up==this) return EEE_NOEFF;    if (this==m0_kcp[0])    return NDE_PERM;
	RMTHAT; that->m_up = this; that->m_u24.c.ct='V'; m_extra = that->id(); return 0;
}

int ClipNode::add(ANode * that, const char * nm, int i, int j) {
	int ix;
	if(i&NOF_PARSE) { if(nm?((ix=b32_to_i(*nm))<0):(ix=m_sel,0)) return *nm==63 ? add_hlp(that):NDE_PARSE;}
	else		{ if (~i&63) ix = (i&NOF_NOIDX) ? m_sel : i&31; else return add_hlp(that); }
	if (that->m_up==this) return xchg(that->m_u24.c.i, ix);
	unsigned int m1 = 1u << ix;
	if (!that->is_wrap()) return NDE_EXPWRAP;
	if (m_map & m1) {
		if (!~m_map) return NDE_FULL;
		if (i&NOF_STRICT) return NDE_NDUP;
		ix = __builtin_ffs(~(m_map | (m1-1))) - 1;
		if (ix < 0 && (ix = __builtin_ffs(~m_map)-1) < 0) return NDE_WTF;
		m1 = 1u << ix;
	}
	RMTHAT; m_map |= m1;
	if (ix<0 || ix>31) bug("clip/add: ix=%d", ix);
	m_eh[ix] = (unsigned short) (that->m_id >> 8);
	m_el[ix] = (unsigned char) (that->m_id & 255);
	that->m_u24.c.ct = 'c'; that->m_u24.c.i = ix; that->m_up = this;
	DIRLOOP2 gui2.node_name(b, that);
	if (winflg(8) || (i&NOF_FGUI)) show_newbox(ent_j(ix));
	return sel(ix);
}

int ClipNode::show_dsc(int f) { return !m_extra ? NDE_NODSC : f&1 ? lookup_n_q(m_extra)->draw_window(16) : 0; }
void ClipNode::draw_1(ABoxNode * nd) { gui2.clip_box(this, nd->m_u24.c.i); }
void ClipNode::show_newbox(ABoxNode * nd) { if (!nd->box()) return nd->winflg_or(8);
	if (winflg(8)) gui2.clip_box(this, nd->m_u24.c.i); else draw();  }

int ClipNode::rm(ANode * that) {
	int i = that->m_u24.c.i;
	unsigned int m1 = 1u << i; 
	if (!(m_map&m1) || ent_j(i)!=that) return (m_extra==that->id()) ? (m_extra=0) : NDE_RMWHAT;
	m_map &= ~m1; if (winflg(8)) gui2.clip_box(this, i);
	DIRLOOP2 gui2.node_rm(b, that);
	return 0;
}

int ClipNode::gui_list(char * to, int flg) {
	if (!m_map) return 0;
	char * p = (flg&1) ? to : (*to='x', to+1);
	for (unsigned int j,m = m_map; j=__builtin_ffs(m), j-- > 0; m &= ~(1u<<j)) {
		int k = 256*m_eh[j] + m_el[j];
		if (p>to+1) *(p++) = 36;
		if (flg&1) *(p++) = max_i(32, lookup_n_q(k)->cl_id());
		p += hx5(p, k); *(p++) = 36; *(p++) = i_to_b32((int)j);
	}
	return p - to;
}

int ClipNode::wfind(BoxGen * bx) {
	for (unsigned int j,m = m_map; j=__builtin_ffs(m), j-->0; m &= ~(1u<<j)) 
		if (!wrap_qdiff(bx, bx_j((int)j))) return (int)j;
	return -1;
}

int ClipNode::sel(int i) {
	if (m_sel==(i&=31)) return i; else m_sel = i;
	if (winflg(8)) gui2.clip_box(this, -1, -1, i);
	return i;
}

BoxGen * ClipNode::find_bx(BoxGen * cb, int stp, int msk) {
	int j0 = (cb->node()->m_u24.c.i), j0h = j0 & ~msk, j = j0;
	while ((j = ((j+stp)&msk)|j0h) != j0) if (m_map&(1u<<j)) return ent_j(j)->box();
	return cb;
}

int ClipNode::keyop_j(int j, int ky, int op, const char *s, int nof) {
	return (m_map&(1u<<j)) ? wrap_key_op(ent_j(j)->box(), ky, op, s, nof) : EEE_NOEFF; }

#define KO_COND(J) ((r = wrap_key_op(ent_j(J)->box(),ky,op,s,nof))!=BXE_UNDEFKEY)
int ClipNode::keyop_f(int ky, int op, const char *s, int nof) {
	unsigned int m0 = m_map, m1 = (1u<<m_sel);    ky |= 65536;
	int r, rot = m_sel&24;   if (m0&m1) { if (KO_COND(m_sel)) return r; else m0 &= ~m1; }
	BVFOR_JM((m0>>rot)|(m0<<(32-rot))) { if (KO_COND((j+rot)&31)) return r; } return EEE_NOEFF; }

int ClipNode::cmd(CmdBuf * cb) {
	char * s = cb->tok();
	int k, f = cb->cnof();
	BoxGen * p; ABoxNode * nd;
	ClipNode ** q;
	switch(*s) {
		case 'W': draw(); return 0;
		case 'D': k = 1; goto flg;
		case 'A': k = 2; goto flg;
		case 'C': q = m0_kcp+1; k=4; goto io;
		case 'P': q = m0_kcp+2; k=8; goto io;
		case 'X': k = 16; goto flg;
		case 'd': return (m_map & (1u<<m_sel)) ? Node::del(ent_j(m_sel), f) : EEE_NOEFF;
		case '1': case '4':
			  if (m_sel==(k=s[1]-48)) return (m_flg&2) ? keyop_j(k, 3, 10, 0, 0) : EEE_NOEFF;
			  keyop_j(m_sel, 3, 9, 0, 0); sel(k);
			  return (m_flg&2) ? keyop_j(m_sel, 3, 11, 0, 0) : 0;
		case '2':
			  k = (s[1]-48)&31; p = (m_map&(1u<<k)) ? bx_j(k) : 0;
			  if (p)  return Node::mk(0, m0_kcp[1], 0, 'W',
					          f|NOF_NOIDX|(NOF_OVRRD&-(this==m0_kcp[1])), 0, p);
			  else if ((p = m0_kcp[2]->bx_sel()))
				  return Node::mk(0, this, 0, 'W',
						  k|f|(NOF_OVRRD&-(this==m0_kcp[2])), 0, p);
			  else return EEE_NOEFF;
		case '3':
			  k = (s[1]-48)&31;
			  if (!(m_map&(1u<<k)) || (k!=m_sel && (m_flg&16))) goto xcg;
			  return ent_j(k) -> draw_window(0x19); 
		case '9': return (m_map&(1u<<(k=(s[1]-48)&31))) ? ent_j(k)->draw_window(0x1b) : EEE_NOEFF;
		case 'k': return keyop_f(hex2(s+1), 1, 0, 0);
		case 'K': return keyop_f(hex2(s+1), 0, 0, 0);
		case 'Z': if (!m_map) return EEE_NOEFF;
			  BVFOR_JM(m_map) if ((k=Node::del(ent_j(j), f))<0) return k;
			  if (winflg(8)) draw();    return 0;
		case 'G': k=0;  BVFOR_JM(m_map) if ((nd=ent_j(j))->cl_id()=='s') 
			  				k=min_i(k,swrap_grab_c(nd->box(), 0));  return k;
		default:  return GCE_UCLIP;
	}
xcg:
	return cb->cperm(DF_EDDIR) ? xchg(m_sel, k) : NDE_PERM;
flg:
	if (s[1]&1) m_flg |= k; else m_flg &= ~k;
	if (winflg(8)) gui2.clip_flg(m_id, s[0], s[1]&1);
	return 0;
io:
	if (s[1]&1) {
		if (*q==this) return 0;
		if (k==4 && !cb->cperm(DF_EDDIR)) return NDE_PERM;
		ClipNode *old = *q; *q = this;
		if (old->winflg(8)) gui2.clip_flg(old->id(),*s,0);
		if (winflg(8)) gui2.clip_flg(m_id,*s,1); 
	} else if (*q==this && *q!=m0_kcp[0]) {
		*q = m0_kcp[0];
		if (winflg(8)) gui2.clip_flg(m_id, *s, 0);
		if ((*q)->winflg(8)) gui2.clip_flg((*q)->id(), *s, 1);
		else if (f & NOF_FGUI) (*q)->draw();
	}
	return 0;
}

void ClipNode::draw() {
	gui2.cre(16*m_id+3, 'K', "!,");
	gui2.npath(this, 333); 
	gui2.c3(36, i_to_b32(m_sel), 48 + m_flg + 4*(this==m0_kcp[1]) + 8*(this==m0_kcp[2]));
	for (unsigned int j,m = m_map; j=__builtin_ffs(m), j-- > 0; m &= ~(1u<<j))
		gui2.clip_box_1(this, (int)j);
}

void ClipNode::debug(int flg) {
	ad_debug();
	log("map: %08x flg: %x", m_map, m_flg);
}

///////// box ////////////////////////////////////////////////////////////////

int ABoxNode::ccmd(CmdBuf * cb) { return m_box->cmd(cb); }
int ABoxNode::start_job_2(JobQ::ent_t * ent, char * arg) { return m_box->start_job_3(ent, arg); }

int ABoxNode::get_ionm(char *to, int io, int j) {
	if (m_box) { int k = m_box->v_get_ionm(to, io, j); if (k) return k; }
	return m_ui.ro()->m_nm[io&1].ro()->get_nm(to, j); }

const char * ABoxNode::rgb() { return m_box ? m_box->v_rgb() : "zz%z%%"; }

void ABoxNode::del2() { 
	BoxEdge * p; while ((p = m_ef0)) Node::disconn(this, p->to);
	if (is_wrap()) m_box->~BoxGen(), f64(m_box); else delete(m_box); }

int ABoxNode::save1() {
	SvArg * p = &m0_sv;
	if (!p->st && m_visitor==p->vis) return p->cn = (p->st = !m_next) ? m_up : m_next, 1;
	if (debug_flags & DFLG_SAVE) sv_dump1("box");
	int ec, r = 1, wf = p->flg & SVF_WRAP;
	BoxEdge * eg = 0; ANode * snl;

	if (p->st > 7) return save_sob();
	if ((p->st & 1)) { switch(p->st) { 
		case 1:goto s1;case 3:goto s3;case 5:goto s5;default:return NDE_WTF; }}
	m_visitor = p->vis; m_winflg &= 0x1fffffff; m_winflg |= p->st << 29;
	if (!wf) { for (eg = m_ef0; eg; eg = eg->fr_n) if (eg->to->is_save_trg()) goto edg; }
s5:	CHKERR(sv_cre()); m_winflg |= WF_SVNC; if (debug_flags & DFLG_SAVE) log("svnc_flg set1: 0x%x", m_id); 
	p->st = 8; p->st2 = 0; CHKERR(save_sob()); return r;
s3:	if ((snl = sn_list(&m0_sv.wl))) log("BUG: unexp. non-wrap subbox 0x%x, skipped", snl->id()); // TODO
	p->st = 1;
s1:	if ((ec=m_box->save2(&m0_sv)) < 0) return ec; else r += ec;
	if (debug_flags & DFLG_SAVE) log("bxsv/s1: 0x%x: from: %d, nx:%p", m_id, m_winflg & WF_SVFR, m_next);
	if (!(m_winflg & WF_SVFR)) return m0_sv.nxup(r);
	for (eg = m_et0->fr_n; eg; eg = eg->fr_n) if (eg->to->is_save_trg()) goto edg;
	if (debug_flags & DFLG_SAVE) log("bxsv/s1: 0x%x: ret 0x%p", m_id, m_et0->fr);
	p->cn = m_et0->fr; p->st = 5; return r;
edg:	Node::eg_mv_to0(eg); p->st = 2; p->cn = eg->to; return r + 5;
}

int ABoxNode::sv_wr_backref() {
	AOBuf * f = m0_sv.out;
	int ec, r=1;  ABoxNode * fnd;
	for (BoxEdge * eg = m_et0; eg && (fnd=eg->fr)->m_visitor==m0_sv.vis; eg = eg->to_n) {
		CHKERR(f->sn("X$<", 3)); CHKERR(fnd->sv_path(10)); }
	return r;
}

int ABoxNode::save_sob() {
	char * ps = &m0_sv.st;
	if (debug_flags & DFLG_SAVE) log("save_sob: id=0x%x, st=%d", id(), *ps);
	if (m0_sv.st2<0) ++*ps, m0_sv.st2=0;
	if (*ps>15) return m_box->save_sob(&m0_sv);
	switch(*ps) {
		case 8: 	  return m_ui.save(&m0_sv);
		case 9: 	  return m_ui.ro()->m_dv.save(&m0_sv);
		case 10: case 11: return m_ui.ro()->m_nm[*ps&1].save(&m0_sv);
		case 12:	  return m_ui.ro()->m_dsc.save(&m0_sv);
		default:	  return *ps=16, 1;
	}}

int ABoxNode::cond_del() {
	if (m_et0) return NDE_USEDBY;
	if (m_winflg&1) return NDE_BUSY;
	return m_box -> cond_del();
}

void ABoxNode::unset_model_rec() {
	for (BoxEdge * p = m_et0; p; p=p->to_n)
		if (p->fr_bx->unset_model_1()) p->fr->unset_model_rec();
}

void ABoxNode::nio_change() {
	for (BoxEdge * p = m_et0; p; p=p->to_n) p->fr_bx->notify_nio(m_box); }

int ABoxNode::ui_cmd(CmdBuf * cb) {
	char *s = cb->a1(); if (*s=='W') return draw_window(0x1c);
	if (!cb->cperm(DF_EDBOX)) return NDE_PERM;
	switch(*s) {
		case 'D': case 'E': {
			char buf[4096]; buf[0]='h'; buf[1]=*s; h5f(buf+2, m_id); buf[7]='.';
			int l = 9+get_path_uf(buf+8,70); buf[l-1] = '$';
			l += m_ui.ro()->dump_dsc(buf+l,0);
			buf[l++] = 10; return pt_iocmd_sn(buf, l); }
		case 'Q': {
			ANode * from = m0_lr[s[1]&1];  const char * snnm="?";
			if (from->cl_id()=='C' && !(from=from->sn(&snnm))) return NDE_NOSOBCP;
			if (!from->is_box()) return BXE_ARGNBX;
			BoxDesc * dsc = static_cast<ABoxNode*>(from)->m_ui.ro()->m_dsc.ro();
			return dsc ? (SOB_RW(ui)->m_dsc.set(dsc), cond_docw(), 0) : NDE_NOSOBCP; }
		case 'O': {
			ANode * nd2a = cb->lookup(s+2); if (!nd2a) return BXE_ARGLU;
			ABoxNode * nd2 = dynamic_cast<ABoxNode*>(nd2a); if (!nd2) return BXE_ARGNBX;
			int k = s[1] & 7; 
			switch(k) {
				case 0: m_ui.from(nd2->m_ui); return 0;
				case 1: SOB_RW(ui)->m_dv.from(nd2->m_ui.ro()->m_dv); return 0;
				case 2:case 3: SOB_RW(ui)->m_nm[k-2].from(nd2->m_ui.ro()->m_nm[k-2]); return 0;
				case 4: SOB_RW(ui)->m_dsc.from(nd2->m_ui.ro()->m_dsc); return 0;
				default: return NDE_WTF;
			}}
		default:{
			int ec = SOB_RW(ui)->cmd(s); 
			if ((ec&0x80000001)==1 && winflg(4096)) m_ui.ro()->w_rgb(16*m_id+12);
			return ec; 
			}}}

int ABoxNode::draw_window_2(int x) {
	if (!m_box) return NDE_NULLBOX; else Node::shl_add(this);
	switch(x) {
		case 0: case 0xb: return m_box->box_window(), 16*id()+0xb;
		case 0x9: return m_box->aux_window();
		case 0xe: return m_box->extra_window();
		case 0xc: return (m_box->ifflg()&BIF_GC) ? m_ui.ro()->draw_window_2(this) : NDE_NOGUI;
		case 0xd: return (m_box->ifflg()&BIF_GC) ? (m_box->doc_window(13),16*id()+13) : NDE_NOGUI;
		default: return BXE_CENUM;
	}}

int ABoxNode::find_fw_2(ABoxNode * to, LWArr<int>* rpath) {
	for (BoxEdge * p = m_ef0; p; p=p->fr_n) {
		if (p->to==to) {
			if (rpath) rpath->resize(1), *rpath->p() = m_id; return NDE_LOOPBX; }
		if (p->to->find_fw_2(to, rpath)) {
			if (rpath) rpath->add(m_id); return NDE_LOOPBX; }}
	return 0;
}

void ABoxNode::mwin_head(int wid) {
	int oid = 16*m_id + 11;
	if (!winflg(8)) gui2.cre(oid, m_u24.s[0]); else gui2.setwin(oid, m_u24.s[0]);
	gui2.wupd('_',0); gui2.c1('!'); gui2.sn(rgb(), 6);
	gui2.npath(this, wid);
}

int ABoxNode::wdat_alloc() {
	if (m_winflg & WF_WI6) return 0;
	if (!m0_wi_bf) return NDE_BOX63;
	int i = m0_wi_bf;
	sthg * p = m0_wi_b + 4*i;
	m0_wi_bf = p->i[0] & 63; p->i[0] = m_id;
	m_box -> wdat_cons(p);
	m_winflg |= 65536 * i;
	return i;
}

int ABoxNode::wdat_free() {
	int i = m_winflg & WF_WI6;
	if (!i) return 0; else i >>= 16, m_winflg &= ~WF_WI8;
	sthg * p = m0_wi_b + 4*i;
	m_box -> wdat_del(p);
	p->i[0] = m0_wi_bf + 1048576; m0_wi_bf = i;
	return 1;
}

void ABoxNode::qc_rgb(int c, const char * s) { memcpy(m_ui.rw_o0()->m_rgb+3*(c=='-'), s, 3+3*(c=='*')); }
void ABoxNode::qc_dfv(int bv, double * dv) { m_ui.rw_o0()->m_dv.rw_o0()->add_bv(bv, dv); }
void ABoxNode::qc_ptn(int trg, int i, int j, const char *s) {m_ui.rw_o0()->m_nm[trg=='o'].rw_o0()->qpt(i,j,s);}
void ABoxNode::qc_iot(int trg, int bv, const char *s) {  BoxUI * ui = m_ui.rw_o0();
	if ((trg|=32)=='t') bug("strncpy(ui->m_dsc.rw_o0()->m_u.ky, s, 47)"); 
	else ui->m_nm[trg=='o'].rw_o0()->set_nm(bv, s); }
void ABoxNode::qc_st8(int trg, int j) {
	static SOB_p<BoxUI> store[8];
	if (!(j&=7))       { if (trg=='R') trg='r';
			     if (!store[0].ro()) store[0].set(BoxUI_default(0)); }
	if ((trg|32)=='r') { if (trg&32) m_ui.from(store[j]); else store[j].from(m_ui); return; }
	BoxUI * ui = m_ui.rw_o0();
	switch(trg|32) {
		case 'i': case 'o': trg = (trg=='o'); ui->m_nm[trg].from(store[j].ro()->m_nm[trg]); return;
		case 't': ui->m_dsc.from(store[j].ro()->m_dsc); return;
		case 'v': ui->m_dv. from(store[j].ro()->m_dv ); return;
		default:  return log("BUG: qmk/st8: unexp. trg 0x%x (%c) (skipped)", trg, trg); 
	}}

int ABoxNode::dsc(char * to) {
	if (cl_id() != '_') return m_ui.ro()->dump_dsc(to, 1);
	ANode * p; int j,n = 1; *to = 33;
	for (p = m_up; !p->is_dir(); p = p->up()) ;
	p = static_cast<ADirNode*>(p)->get_hroot();
	n += ((j=p->id())<3) ? (j && (to[n] = 3+30*j,1)) : p->get_name(to+n);
	to[n++] = '.'; return n + s__cat(to+n, m_box->cl_name());
}

int ABoxNode::cmd_H(CmdBuf * p) {
	BoxUI* ui = SOB_RW(ui);
	char *av[33]; int n, e, rv = 0;
	for (n=0; n<33; n++) if (!(av[n]=p->tok())) break;
	if (!n) { ui->m_dsc.set(0); }
	else {  n -= (e = (n==33));
		BoxDesc * dsc = SOB_RWP(ui, dsc);  dsc->set_n(n);
		for (int i=0; i<n; i++) { char *s = av[i], *q = dsc->ln(i); int l = strlen(s); 
			          	  if (l>63) memcpy(q,s,63),q[63]=0,e|=2; else memcpy(q,s,l+1); }
		if (e) rv = NDE_DSC_Y-1+e;
	}
	return cond_docw(), rv;
}

int ABoxNode::show_dsc(int flg) {
	return !(m_box->ifflg()&BIF_GC) ? NDE_NOGUI :
	       !(m_ui.ro()->m_dsc.ro()) ? NDE_NODSC :
	       (flg&1) ? draw_window(0x1d) : 0;     }

void ABoxNode::ab_debug(int flg) {
	char nm[24]; nm[get_name(nm)] = 0;
	log_n(" ("); m_box->model0()->debug0(); log_n(")"); if (!flg) return; 
	a_debug(); if (!(flg&1)) goto d2;
	log("box=%p, box_cl=\"%s\"", m_box, m_box->cl_name()); 
	// TODO: minidir
	log_n("uses: ");
	for (BoxEdge * p = m_ef0; p; p=p->fr_n)
		if (p->fr != this) log_n(" (bug)"); else log_n(" %s", p->to->path255());
	log(""); log_n("used by:");
	for (BoxEdge * p = m_et0; p; p=p->to_n)
		if (p->to != this) log_n(" (bug)"); else log_n(" %s", p->fr->path255());
	log(""); 
d2: 	if (flg&4) m_ui.debug(); 
	if (m_box && (flg&2)) m_box -> spec_debug();
}

///////// trk ////////////////////////////////////////////////////////////////

static int trkn_parse(int * pi, int * pj, const char **pp) {
	int ti = 0, tj = 0;
	for (int i=0; i<3 && is_hx(**pp); i++, ++*pp) ti*=16, ti+=hxd2i(**pp);
	if (**pp!=':') return NDE_PARSE; else ++*pp;
	for (int i=0; i<8 && is_hx(**pp); i++, ++*pp) tj*=16, tj+=hxd2i(**pp);
	*pi = ti; *pj = tj; return 0;
}

ANode * TBoxNode::find_ijf(int i, int j, int nxf, ANode * nd) {
	ANode *r = nd ? nd : trk_bkm_find(m_box, j);
	int di, dj = r->cth()->j-j;
	if (!dj) { if (!(di = r->cth()->i-i)) return r;
		   if (di<0) { r=r->next(); goto up; } else { r=r->cth()->pv; goto dn; }}
	if (dj<0) { do r=r->next();    while ((dj=r->cth()->j-j)<0); if (dj) goto up2; goto up; }
	else	  { do r=r->cth()->pv; while ((dj=r->cth()->j-j)>0); if (dj) goto dn2; goto dn; }
up:	while (!(dj=r->cth()->j-j) && (di=r->cth()->i-i)<0) r=r->next();    if (di|dj) goto up2; return r;
dn:	while (!(dj=r->cth()->j-j) && (di=r->cth()->i-i)>0) r=r->cth()->pv; if (di|dj) goto dn2; return r;
up2:	return nxf ? r 	       : 0;
dn2:	return nxf ? r->next() : 0;
}

ANode * TBoxNode::sn(const char **pp) {
	int i, j, ec = trkn_parse(&i, &j, pp);
	return ec<0 ? 0 : find_ijf(i, j, 0); }

void TBoxNode::do_ins(ANode * q, int i, int j, ANode *nx) {
	ANode * pv=0; trk_24 *ti, *tn;
	if (!nx || !(pv=(tn=&nx->m_u24.t)->pv)) return bug("trk/do_ins: nx=%p, pv=%p", nx, pv);
	ti = &q->m_u24.t; pv->m_next = tn->pv = q; 
	ti->pv = pv; ti->ct = 't'; ti->i = (short)i; ti->j = j; q->m_up = this; q->m_next = nx; 
	TRK_CO2(pv, q, "ins1"); TRK_CO2(q, nx, "ins2");
}

void TBoxNode::do_cut(ANode * q) {
	ANode *nx = q->next(), *pv = q->cth()->pv; if (!nx || !pv) bug("trk/do_cut: nx=%p, pv=%p", nx, pv);
	pv->m_next = nx; nx->m_u24.t.pv = pv; TRK_CO2(pv, nx, "cut"); }

int TBoxNode::rm(ANode * that) {
	if (that->m_up != this) return NDE_RMWHAT;
	if (!that->is_wrap()) return log("BUG: trk/rm: %s 0x%x/0x%x", that->cl_id()==33?"guard":"unexp",
			that->id(), that->cl_id()), NDE_WTF;
	if (m_winflg & 2048) trk_cond_pm(m_box, that, '-');
	do_cut(that); trk_bkm_rm(m_box, that);
	return 0;
}

int TBoxNode::add(ANode * that, const char * nm, int i, int j) {
	int ty = that->m_u24.c.ty; if ((ty|4)!='w') return NDE_EXPWRAP;
	int k, ec, sf = !!(i&NOF_STRICT), wf = m_winflg & 2048;
	if (!(i&NOF_PARSE)) i&=NOF_IDX; else if ((ec = trkn_parse(&i, &j, &nm)) < 0) return ec;
	if ((unsigned int)(i-16)>4079u) return TKE_RANGEI;
	if ((unsigned int)(j) > 2147483520u) return TKE_RANGEJ;
	if (!sf) i&=4080;
	ANode *nx = find_ijf(i, j, 1);
	if (nx->cth()->j == j && nx->cth()->i == i) {
	 	if (sf) return NDE_TDUP; 
		for (++i,k=1; k<16; k++,i++) if (nx=nx->m_next, nx->cth()->j!=j || nx->cth()->i!=i) goto iok;
		return TKE_NOROOM;
iok:    ;}
	ANode *pv = nx->cth()->pv;
	if (nx==that || pv==that) { if (debug_flags&DFLG_TRK) log("trkadd/lmv");
				    trk_bkm_rm (m_box, that);
			  	    if (wf)trk_cond_pm(m_box, that, '-');
				    that->m_u24.t.i = i; that->m_u24.t.j = j; TRK_CO1(that, "lmv"); }
	else { if  (debug_flags&DFLG_TRK) log("trkadd/rm?"); RMTHAT; do_ins(that, i, j, nx); }
	trk_bkm_add(m_box, that); if (wf)trk_cond_pm(m_box, that, '+');
	TRK_CO1(that, "add"); return 0;
}

void TBoxNode::ins_guard(ANode * that, int j, ANode * nx) {
	if (j<0) { if (!nx) bug("ins_guard: 00"); else j = nx->cth()->j; }
	const trk_24 * th = that->cth(); 
	if (th->ty!=33 || (th->i&~15)) bug("ins_guard: c=0x%x, i=0x%x", th->ty, th->i);
	nx = find_ijf(th->i, j, 1, nx);
	if (nx==that || nx->cth()->pv==that) { that->m_u24.t.j = j; TRK_CO1(that,"glm"); return; }
	if (that->m_next) do_cut(that); 
	do_ins(that, th->i, j, nx);
}

void TBoxNode::chk_ord(ANode * p, ANode *q, const char * dsc) {
	const trk_24 *a = p->cth(), *b = q->cth();
	if (!a->ty || !b->ty) bug("chk_ord/%s: %p,%p ty:%02x/%02x", dsc, p,q,a->ty,b->ty);
	if (p->m_up!=this || q->m_up!=this) 
		bug("chk_ord/%s: upnodes %p,%p / exp. %p", dsc, p->m_up, q->m_up, this);
	if (p->m_next!=q || b->pv != p) 
		bug("chk_ord/%s: link err (p=%p,nx=%p,q=%p,q->pv=%p)", dsc, p, p->m_next, q, q->cth()->pv);
	int di = b->i - a->i, dj = b->j - a->j;
	if (dj<0) goto obug; if (dj>0) return;
	if (di<0) goto obug; if (di>0 || b->i<16) return;
	bug("chk_ord/%s (duplicate %x:%x)", dsc, a->i, a->j);
obug:	bug("chk_ord/%s (%x:%x before %x:%x)", dsc, a->i, a->j, b->i, b->j);
}

void TBoxNode::debug(int flg) {
	ab_debug(flg); log_n("trk-nodes: [");
	ANode * nd = trk_bkm_find(m_box, -1); 
	while (1) { log_n(" #%x", nd->id()); if (nd->m_u24.t.j==0x7fffffff) break; nd = nd->next(); }
	log("]"); }

void   tgn_del (ANode *tn, ANode *gn) { if (!gn) return;
	if (gn->up()!=tn || gn->cl_id()!=33) bug("tgn_del"); 
	if(gn->next())static_cast<TBoxNode*>(tn)->do_cut(gn); else bug("tgn_del/!n"); ANode::fN(gn);}

void   tgn_move(ANode *tn, ANode *gn, int j, ANode*nx) { if (gn->up()!=tn || gn->cl_id()!=33) bug("tgn_move");
	static_cast<TBoxNode*>(tn)->ins_guard(gn, j, nx); }

ANode* tgn_new (ANode *tn, int i, int j, ANode *nx) { TRK_SANE2(tn->box0(), 6);
	ANode * r = new (ANode::aN()) TGuardNode(tn, i); 
	static_cast<TBoxNode*>(tn)->ins_guard(r, j, nx); return r; }

void Node::trk_chk_ord(ANode *tn, ANode *p, ANode *q, const char * msg) {
	static_cast<TBoxNode*>(tn)->chk_ord(p,q,msg); }

ANode * Node::trk_ij(ABoxNode * tn0, int i, int j) { 
	return static_cast<TBoxNode*>(tn0) -> find_ijf(i,j,1); }

ANode * Node::trk_fwf(ANode *q) { while (q->cth()->i<15) q = q->next();
				  return q->cth()->i>15 ? q : 0; }
ANode * Node::trk_fwb(ANode *q) { while ((unsigned int)(q->cth()->i-1) < 15u) q = q->cth()->pv;
				  return q->cth()->i>15 ? q : 0; }

///////// static /////////////////////////////////////////////////////////////

unsigned int Node::m0_visitor = 0;
int Node::m0_slr_flg;
int Node::m0_shl_n = 0, Node::m0_shl[32];
ADirNode * Node::m0_clibroot = 0;

int Node::hier(ANode * up, ANode *dn) { while(1) { if (!dn) return 0; if (dn==up) return 1; dn=dn->m_up; }}

ANode * Node::lookup_cb(CmdBuf * cb, const char * s) {
	ANode * r; int i, j, c = *(s++);
	if (c&1) {
		if (c=='#') { r = ANode::lookup_n_q(atoi_h_ppc(&s)); if (!r->m_u24.s[0]) return 0; }
		else if (c!='=' || !(r = cb->var(*(s++)&7))) { return 0; } 
	} else if (c=='.') { r = ANode::root(); } else {
		if (c!='@') return 0; else c = *(s++);
		switch(c) {
			case 'L': r = ANode::m0_lr[0]; if (!r->is_dir()) r = r->m_up; break;
			case 'R': r = ANode::m0_lr[1]; if (!r->is_dir()) r = r->m_up; break;
			case 'l': r = ANode::m0_lr[0]; break;
			case 'r': r = ANode::m0_lr[1]; break;
			case 'k': r = ClipNode::m0_kcp[0]; break;
			case 'c': r = ClipNode::m0_kcp[1]; break;
			case 'p': r = ClipNode::m0_kcp[2]; break;
			default:  if ((unsigned int)(i=c-48) > 9u) return 0;
				  for (; (unsigned int)(j=*s-48)<10u; s++) i = 10*i + j;
				  r = ANode::lookup_n_q(i); if (!r->m_u24.s[0]) return 0;
				  break;
		}
	}
	do { while (*s=='.') ++s; } while (*s && (r = r->sn(&s)));
	if (cb->cnof() & NOF_FGUI) shl_add(r);
	return r;
}

void Node::set_lr(int i, ANode * p) { 
	int wf = WF_LSEL << (i&=1);
	ANode * p0 = ANode::m0_lr[i]; if (p==p0) return; else slr_invd(wf);
	if (p0) { p0->winflg_and(~wf); if (!p0->is_dir()) p0->m_up->winflg_and(~wf); }
	p -> winflg_or(wf); if (!p->is_dir()) p->m_up->winflg_or(wf);
	gui2.t2_sel(i, ANode::m0_lr[i] = p);
}

int Node::parse_target(CmdBuf * cb, char ** ppname, ANode** ppto) {
	char *s = *ppname, *pd = 0;
	for (;*s;s++) if (*s=='.') pd = s;
	if (!pd) return *ppto = 0, NDirNode::hash(*ppname);
	*pd = 0; *ppto = (pd>*ppname) ? lookup_cb(cb, *ppname) : ANode::root(); *pd = '.';
	return *ppto ? ( pd[1] ? NDirNode::hash(*ppname=pd+1) : NDE_ZNAME ) 
		     : NDE_TRGLU;
}

int Node::obj_help(int cl) { 
	ANode * nd = ANode::lookup_n_q(4);
	const char * s = 0;
	if (debug_flags & DFLG_GUICMD) log("obj_help: 0x%x", cl);
	switch(cl) {
		case 'C': s = "clipboard"; break;
		case 'D': nd = ANode::lookup_n_q(5), s = "object tree"; break;
		case 'w': case 's': s = "wrap/config"; break;
		case 'g': s = "graph-box"; break;
		case 'c': s = "calc-box"; break;
		case 't': s = "track"; break;
		case 'i': s = "iterated box"; break;
		case 'h': s = "text"; break;
		case 'w'+256: case 's'+256: s = "wrap/play"; break;
		case '_': s = "primitive box"; break;
		default: return NDE_WTF;
	}
	nd = nd->sn(&s); if (!nd) return NDE_WTF;
	return nd->draw_window(16);
}

sbfun_t setbox_wrap, setbox_shwr, setbox_shtg, setbox_graph, setbox_calc, setbox_it, setbox_hlp;

int Node::sb_btin(ABoxNode * nd, BoxGen * bx) { (nd->m_box = bx) -> set_node(nd); return 0; }
int Node::sb_trk (ABoxNode * nd, BoxGen * bx) {
	TGuardNode *g0, *g1; trk_24 * t24;
	t24 = &(g0 = new (ANode::aN())TGuardNode(nd, 0))->m_u24.t; t24->pv=0;  t24->ct='t'; t24->j=0;
	t24 = &(g1 = new (ANode::aN())TGuardNode(nd,15))->m_u24.t; t24->pv=g0; t24->ct='t'; t24->j=0x7fffffff;
	g1->m_next=0; g0->m_next=g1;
	return nd->m_box = trk_mk(nd, g0, g1), 4;
}

ABoxNode * Node::qcp2(BoxGen *bx) {
	if (!bx || !(bx->ifflg()&BIF_QCP)) return 0;
	ABoxNode *ndf = bx->node(), *ndt = new (ANode::aN()) LBoxNode(ndf->cl_id());
	return bx->qcp3(ndt) ? (memcpy(ndt->etc(),ndf->etc(),8), ndt->set_ui_f(ndf), ndt) 
			     : (ANode::fN(ndt), (ABoxNode*)0);
}
		
int Node::mk(ANode ** rr, ANode * up, const char * name, int ty, int i, int j, BoxGen* from) {
	if (!up) return NDE_NOUP;
	if (!(i&NOF_FORCE) && !up->perm_ed()) return NDE_PERM;
	if ( (glob_flg & (GLF_EMPTY|GLF_LIBMODE|GLF_INI0)) == GLF_EMPTY) glob_flg &= ~GLF_EMPTY;
	ANode * nd; ABoxNode * bnd;
	int ec; sbfun_t * sbf = 0;
	switch(ty) {
		case 'C':           nd = new (ANode::aN()) ClipNode(); goto ndok;
		case 'D': case 'd': nd = new (ANode::aN()) NDirNode(); goto ndok;
		case 'w': sbf = setbox_wrap;  goto lb;
		case 's': sbf = setbox_shwr;  goto lb;
		case 'g': sbf = setbox_graph; goto lb;
		case 'c': sbf = setbox_calc;  goto lb;
		case 'i': sbf = setbox_it;    goto lb;
		case 'h': sbf = setbox_hlp;   goto lb;
		case '_': sbf = sb_btin;      goto lb;
		case 'S': sbf = setbox_shtg; ty='s';  goto lb;
		case 't': nd = bnd = new (ANode::aN()) TBoxNode('t'); sbf = sb_trk;  goto ndok;
		case 'W': if (!(nd=bnd=qcp2(from))) return from ? NDE_NOQCP:NDE_WTF; goto ndok;
		case '!': return log("BUG: Node::mk(tguard)"), NDE_WTF;
		default: return NDE_UTYPE;
	}
lb:   	nd = bnd = new (ANode::aN()) LBoxNode(ty);
ndok:	if ((ec=up->add(nd, name, i, j)) < 0) return ANode::fN(nd), ec;
	if (sbf) { if ((ec=(*sbf)(bnd, from))<0) return ec; else bnd->m_ui.set(BoxUI_default(ec)); }
	m0_slr_flg |= up->winflg(WF_2SEL);
	(i&NOF_FGUI) && (shl_add(nd), ty!='W') && nd->draw_window(16);
	if (rr) *rr = nd; return 0;
}

int Node::move(ANode * p, ANode * to, const char * name, int i, int j) {
	ANode * up = p->up(); if (!up) return NDE_NOROOT; if (!to) to = up;
	if (!(i&NOF_FORCE)) { if (!up->perm_ed() || (up!=to && !to->perm_ed())) return NDE_PERM;
			      if (p->is_dir()) { ADirNode * q = static_cast<ADirNode*>(p);
				      		 if (q->m_pm_msk & ~q->m_pm_val & DF_EDDIR) return NDE_PERM; }}
	if (hier(p, to)) return NDE_HIERMV;
	int ec = to->add(p, name, i, j); if (ec<0) return ec;
	int flg = (p->m_winflg|up->m_winflg|to->m_winflg) & WF_XSEL; if (!flg) return ec;
	m0_slr_flg |= flg; if (!(p->winflg(WF_2SEL))) return ec;
	if (ANode::m0_lr[0]==p) set_lr(0, up);
	if (ANode::m0_lr[1]==p) set_lr(1, up);
	return ec;
}

int Node::copy(ANode * p, ANode * to, const char * name, int i, int j) {
	ANode *up = p->up(), *to2 = to ? to : up; 
	if (!up) return NDE_NOROOT;
	if (!(i&NOF_FORCE) && !(to2->perm_ed())) return NDE_PERM;
	if (hier(p, to)) return NDE_HIERCP;
	ANode * trg = 0; int tid, ec, cl = p->cl_id();
	if (BoxGen * bx = p->box0()) {
		if (cl=='_') return NDE_BTCOPY;
		if (bx->ifflg()&BIF_QCP) return (ec=mk(&trg, to2, name, 'W', i, j, bx))<0 ? ec : trg->id();
	}
	SvArg * sv = &ANode::m0_sv; if (sv->rn) return NDE_LOCK;
	if ((ec = Node::mk(&trg, to2, name, p->cl_id(), i&~NOF_FGUI, j, 0))<0) return ec; else tid = trg->id();
	CmdBuf * cb = new CmdBuf(); 
	cb->init(-1, NOF_FORCE|NOF_STRICT);
	cb->setvar(0, tid); cb->setvar(1, tid); 
	ANode::sv_start(cb, p, SVF_COPY);
	ec = ANode::sv_write(-1); ANode::sv_end(); if (ec<0) return ec;
	if (i&NOF_FGUI) trg->draw_window(16);
	return tid;
}

int Node::del(ANode * p, int flg) {
	int wf, ec; if (!p) return EEE_NOEFF;
	ANode * up = p->up(); if (!up) return NDE_NOROOT;
	if (!(flg&NOF_FORCE) && !p->perm_del()) return NDE_PERM;
	if ((ec = p->cond_del())) return ec;
	if ((ec = up->rm(p))<0) return ec;
	if ((wf = p->winflg())) {
		m0_slr_flg |= (up->m_winflg|p->m_winflg) & WF_XSEL;
		if ((wf&WF_LSEL) && ANode::m0_lr[0]==p) set_lr(0, up);
		if ((wf&WF_RSEL) && ANode::m0_lr[1]==p) set_lr(1, up);
		if (wf&WF_SHL) shl_rm(p);
		if (p->cl_id()!='!') { int oid = 16*p->id(); BVFOR_JM(wf&WF_OIDW) gui2.closewin(oid+(int)j); }
	}
	ANode::del(p); return 0;
}

BoxEdge * Node::find_edge(ABoxNode * fr, ABoxNode * to) {
	for (BoxEdge * p = fr->m_ef0; p; p = p->fr_n) if (p->to==to) return p;
	return 0;
}

int Node::conn(ABoxNode * fr, ABoxNode * to) {
        BoxEdge *p = find_edge(fr, to);
        if (p) return ++p->cnt;
        int ec = to -> find_fw(fr, 0); return (ec<0) ? ec : conn2(fr, to);
}	

#define E_CFR0 (p0 = p->fr_n = fr->m_ef0, fr->m_ef0 = p, (p0) ? (p0->fr_p = p) : (fr->m_efz = p), p->fr_p = 0)
#define E_CTO0 (p0 = p->to_n = to->m_et0, to->m_et0 = p, (p0) ? (p0->to_p = p) : (to->m_etz = p), p->to_p = 0)
#define E_DCFR (pp0=(p0=p->fr_p)?&p0->fr_n:&fr->m_ef0, pp1=(p1=p->fr_n)?&p1->fr_p:&fr->m_efz, *pp0=p1, *pp1=p0)
#define E_DCTO (pp0=(p0=p->to_p)?&p0->to_n:&to->m_et0, pp1=(p1=p->to_n)?&p1->to_p:&to->m_etz, *pp0=p1, *pp1=p0)

int Node::conn2(ABoxNode * fr, ABoxNode * to) {
        BoxEdge *p0, *p = new (ANode::a64()) BoxEdge(fr, to); return E_CFR0, E_CTO0, 1; }

int Node::disconn(ABoxNode * fr, ABoxNode * to) {
	BoxEdge * p = find_edge(fr, to); return p ? (--p->cnt ? p->cnt : disconn2(fr, to, p))
					          : (log("BUG: node/disconnect failed"), 0); }

int Node::disconn2(ABoxNode * fr, ABoxNode * to, BoxEdge * p) {
        BoxEdge *p0, *p1, **pp0, **pp1; return E_DCFR, E_DCTO, ANode::f64((char*)p), 0; }

void Node::eg_mv_to0(BoxEdge * p) { if (!(p->to_p)) return;
	ABoxNode *to = p->to; BoxEdge *p0, *p1, **pp0, **pp1; E_DCTO, E_CTO0; }

void Node::eg_f_mv_to0(ABoxNode * fr, ABoxNode * to) { 
	BoxEdge * p = find_edge(fr, to); if (!p) return log("BUG: node/eg_f_mv_to0: find failed");
	if (!p->to_p) return; BoxEdge *p0, *p1, **pp0, **pp1; E_DCTO, E_CTO0; }

int Node::set_conn_2(ABoxNode * fr, ABoxNode * to1, ABoxNode * to2) {
	BoxEdge *e1, *e2; 
	if (!to1) to1 = to2, to2 = 0;
	if (!to1) return 0;
	if (!(e1=fr->m_ef0)) goto co;
	if ((e2 = fr->m_efz) == e1) {
		if (e1->to==to1) { if (to1==to2) return e1->cnt = 2, 0;
				   e1->cnt = 1; to1 = to2; to2 = 0; goto co; }
		if (e1->to==to2) { e1->cnt = 1; to2 = 0; goto co; }
		disconn2(fr, e1->to, e1); goto co;
	}
	if (e1->to==to1) to1=to2, to2=0; else if (e1->to==to2) to2=0; else disconn2(fr, e1->to, e1);
	if (e2->to==to1) to1=to2, to2=0; else if (e2->to==to2) to2=0; else disconn2(fr, e2->to, e2);
co:     return (to1 && conn(fr, to1)) + (to2 && conn(fr, to2));
}

void Node::shl_add(ANode * nd) {
	//log("m0_shl_n: %d", m0_shl_n);
	int id;  if (!nd || !nd->is_box()) return; else id = nd->m_id;
	if (m0_shl_n && m0_shl[0]==id) return; else slr_invd(WF_SHL);
	for (int i=1; i<m0_shl_n; i++) if (m0_shl[i]==id)
		return (void) (memmove(m0_shl+1, m0_shl, 4*i), m0_shl[0] = id);
	if (m0_shl_n==32) ANode::lookup_n_q(m0_shl[--m0_shl_n])->winflg_and(~WF_SHL);
	memmove(m0_shl+1, m0_shl, 4*m0_shl_n++); m0_shl[0] = id; nd->winflg_or(WF_SHL);
}

void Node::shl_rm(ANode * nd) {
	int i, id = nd->m_id;
	nd -> winflg_and(~WF_SHL);
	for (i=0; i<m0_shl_n; i++) if (m0_shl[i]==id) goto found;
	return;
found:	int n = --m0_shl_n - i;
	if (n) memmove(m0_shl+i, m0_shl+i+1, 4*n);
	slr_invd(WF_SHL);
}

int Node::slr_gui_list(char * to, int ix) { 
	if (ix) return lr_dir(ix-1)->gui_list(to,1);
	char * p = to; int ty, sep = 0;
	for (int i=0; i<m0_shl_n; i++) {
		int id = m0_shl[i];
		ANode * nd = ANode::lookup_n_q(id);
		if (!(ty=nd->m_u24.s[0])) { log("OOPS: deleted obj 0x%x in shortlist", id); continue; }
		if (!sep) sep=36; else *(p++) = sep;
		*(p++) = ty; p += hx5(p, id);
		*(p++) = 36; p += id ? nd->get_name(p) : (memcpy(p,". (root)",8), 8);
	}
	return p - to;
}

int Node::slr_upd_2(char * to) {
	char * p = to;
	for (int i=0, f=WF_SHL; i<3; i++, f+=f) {
		if (!(m0_slr_flg&f)) continue;
		p[0] = 9; p[1] = 61; p[2] = "=<>"[i];
		p += 3 + slr_gui_list(p + 3, i);
	}
	m0_slr_flg = 0; return p - to;
}

int Node::chkwin(int oid) {
	int j = oid&15;
	if (j==7) {
		if (oid==0x17) return -1;
		if (oid==0x27) return -2;
		return 0;
	}
	ANode * nd = ANode::lookup_n_q(oid>>4);
	int ty = nd->cl_id();
	return (ty && nd->winflg(1<<j)) ? ty : 0;
}

int Node::save_batch(ADirNode * dir, const char * fn, int flg) {
	int ec = 0, asvf = 0, snf = 0, n_bk = CFG_SV_BACKUP.i;
	if (!fn || !*fn) {
		if ((!*(fn=save_file_name) && (ec=EEE_NONAME)) || (coward(fn) && (ec=EEE_COWARD)))
			return (flg|NOF_FGUI) ? (gui2.sn("\tf>W", 4), ec) : ec; }
	else if (!memcmp(fn, "/" , 2)) { asvf = 1; fn = QENV('a');   n_bk = CFG_ASV_BACKUP.i; }
	else if (!memcmp(fn, "//", 3)) { asvf = 1; fn = QENV('x'); n_bk = 0; }
	else if (coward(fn)) { return EEE_COWARD; }
	else if (is_asv_name(fn)) { return EEE_OVWASV; }
	else if (!dir->id()) { snf = 1; }
	if (backup(fn, n_bk)<0) gui2.errq_add2(EEE_ERRNO, EEE_BACKUP, fn);
	Clock clk; clk.reset();
	log("saving \"%s\" to \"%s\"...", dir->path255(), fn);
	if (ANode::m0_lock) { if (flg&NOF_FORCE) ANode::sv_end(); else return NDE_LOCK; }
	int exef = !!CFG_SV_EXEC.i, fd = creat(fn, 0644+(0111&-exef)); if (fd<0) return EEE_ERRNO;
	FILE * f = fdopen(fd, "w"); if (!f) return EEE_ERRNO;
	ANode::sv_start(new FOBuf(f), dir, flg);
	if (fprintf(f,"#%s\n", exef ? "!/usr/bin/lflab" : " lflab save file")<=0) return EEE_ERRNO;
	if (fprintf(f,"_V%d.%d\n",v_major,v_minor)<=0) return EEE_ERRNO;
	if (dir->id()) {
		if (fprintf(f,":F:R")<=0) return EEE_ERRNO;
		if ((ec=dir->sv_cre(2))<0) return ec;
		if (fprintf(f,"N$F\n")<=0) return EEE_ERRNO;
	} else {
		if (fprintf(f,":TN.$F\n")<=0) return EEE_ERRNO;
		ClipNode::kcp(0)->sn_list(&ANode::m0_sv.wl);
	}
	ec = ANode::sv_write(-1);
	ANode::sv_end(); ANode::m0_sv.out = 0; 
	log("...save (%s): %s (%d)", fn, ec<0 ? err_str(ec) : "done", ec<0?ec:clk.get());
	if (ec<0) { if (asvf) gui2.errq_add2(ec, EEE_ASVFAIL, "autosave"); return ec; }
	if (!dir->id()) glob_flg |= GLF_SAVED;
	if (snf) strncpy(save_file_name, fn, 1023), glob_flg &= ~GLF_RECOVER, gui2.savename();
	return ec;
}

int Node::lib_cfg(ANode * nd) {
	if (!(glob_flg & GLF_LIBMODE)) return 0; else if (!(nd->id())) return NDE_NOLIB;
	ADirNode * dir = dynamic_cast<ADirNode*>(nd);  if (!dir) return GCE_EXADir;
	(m0_clibroot = dir)->set_perm(DF_ALL, DF_ALL);  return 0;
}

void Node::lib_start() {
	if ((glob_flg&GLF_LIBMODE) || m0_clibroot) log("BUG: lib_start");
	glob_flg |= GLF_LIBMODE; root()->set_perm(DF_ALL, 0); ClipNode::kcp(0)->set_perm(DF_ALL, 0);
}

void Node::lib_end() {
	if (!(glob_flg&GLF_LIBMODE) || !m0_clibroot) log("BUG: lib_end");
	if (m0_clibroot) m0_clibroot->set_perm(DF_ALL, 0), m0_clibroot = 0;
	glob_flg &= ~GLF_LIBMODE; root()->set_perm(DF_ALL,DF_ALL); ClipNode::kcp(0)->set_perm(DF_ALL,DF_ALL);
}

///////// qmk / init /////////////////////////////////////////////////////////

static void qmk_fail(const char * ty, ANode * up, const char * nm, int ec) {
	bug("qmk_%s fail: up:%p(0x%x) nm:\"%s\": %s", ty, up, up->id(), nm, err_str(ec)); }

static ClipNode* qmk_clip(ANode * up, const char * nm) {
	ANode * ret0 = 0;
	int ec = Node::mk(&ret0, up, nm, 'C', 0, 0, 0);
	if (ec<0) qmk_fail("clip", up, nm, ec);
	ClipNode *ret  = dynamic_cast<ClipNode*> (ret0);
	if (!ret) qmk_fail("clip", up, nm, NDE_QMKCAST);
	return ret;
}

ANode * qmk_dir(ANode * up, const char * nm) {
	ANode * ret = 0;
	int ec = Node::mk(&ret, up, nm, 'D', NOF_STRICT|NOF_FORCE, 0, 0);
	if (ec<0 || !ret) qmk_fail("dir", up, nm, ec);
	return ret;
}

ANode * qmk_box(ANode * up, const char * nm, qmb_arg_t qa, int k, int ni, int no, 
		                                           const char * cl, const char * fmt, ...) {
	BoxGen * bx = new (ANode::a64()) PrimBoxGen(qa, k, ni, no, cl+=(*cl=='_'));
	ANode * ret = 0; ABoxNode * bnd = 0;
	if ((*cl&120)==48 && cl[1]<48) { BoxGen** qq = box_bookmark + (*cl&7);
				         if (*qq) log("BUG: duplicate box_bookmark: %s", nm); *qq = bx; }
	int ec = Node::mk(&ret, up, nm, '_', NOF_STRICT|NOF_FORCE, 0, bx);
	if (ec<0 || !ret) qmk_fail("box", up, nm, ec); else bnd = static_cast<ABoxNode*>(ret);
	if (!fmt || !*fmt) return ret;
	va_list vl; va_start(vl, fmt);
	double dv[30];
	int trg = 'r';
	for (const char *s = fmt; *s; s++) {
		if (*s>64) { trg = *s; continue; }
		if ((*s|7)==55) { bnd->qc_st8(trg, *s); continue; }
		if ((trg|32)=='r') { bnd->qc_rgb(*s, va_arg(vl, const char*)); continue; }
		int x, y; switch(*s) {
			case '.': x = va_arg(vl, int); y = 999; goto pt;
			case ':': x = va_arg(vl, int); y = va_arg(vl, int); goto pt;
			case '+': x = va_arg(vl, int); goto ls;
			case '-': x = va_arg(vl, int); y = x>>8; x = ((1<<(x&31))-1)<<y; goto ls;
			case '*': x = (1<<(trg=='o'?no:ni)) - 1; goto ls;
			default: log("qmk_box: ignoring invalid gui conf \"%s\"", s); goto q;
		}
pt:		bnd->qc_ptn(trg, x, y, va_arg(vl, const char*)); continue;
ls:		if (trg=='v')      { for (int n=bitcnt(x), i=0; i<n; i++) dv[i] = va_arg(vl, double);
				     bnd->qc_dfv(x, dv); }
		else if (trg=='V') { bnd->qc_dfv(x, va_arg(vl, double*)); }
		else               { bnd->qc_iot(trg, x, va_arg(vl, const char*)); }}
q:	va_end(vl); return ret;
}

typedef void (ndini_fun)(ANode*);
ndini_fun b_ar_init, b_b0_init, b_help_init, b_in_init, b_map_init,
	  b_filt_pz_init, b_filt_echo_init, b_filt_fe_init, b_filt_misc_init, b_filt_v_init;
#define B_INI(nm) b_##nm##_init(qmk_dir(_b, #nm))
#define B_INI1(up,nm) b_##up##_##nm##_init(_b_##up)
#define B_INI2(up,nm) b_##up##_##nm##_init(qmk_dir(_b_##up, #nm))

void Node::contrib_init(ANode * r) {
	ANode * cdir[26]; 
	int cd_bv = 0;
	for (contrib_ent * p = contrib_list; p->s; p++) {
		int c = p->s[0]|32, j = (c-1) & 31, m = 1<<j; if (j>25) bug("contrib name \"%s\"", p->s);
		ANode * d = qmk_dir((cd_bv&m) ? cdir[j] : (cd_bv|=m, cdir[j]=qmk_dir(r, (char*)&c)), p->s);
		d->set_perm(DF_HROOT, DF_HROOT); (*p->f) (d); }}

void Node::init() {
	ADirNode * r00t = new (ADirNode::aN()) NDirNode();
	if ((sizeof(ClipNode)|sizeof(NDirNode)|sizeof(LBoxNode)|sizeof(TBoxNode))&0xffffff80) bug("node size");
	if (r00t->m_id) bug("nonzero id (0x%x) for root node", r00t->m_id);
	ANode::m0_lr[0] = ANode::m0_lr[1] = r00t;
	r00t->m_winflg |= WF_ROOT; r00t->set_perm(DF_ALL, DF_ALL);
	ANode *_b = qmk_dir(r00t, "!b"), *_c,  *_b_help = qmk_dir(_b, "?"), 
	      *_b_filt = qmk_dir(_b, "filt");
	b_help_init(_b_help);
	B_INI(ar); B_INI(in); B_INI(map); b_b0_init(_b); b_filt_misc_init(_b_filt);
	B_INI1(filt, pz); B_INI1(filt, fe); B_INI2(filt, echo); B_INI2(filt, v);
	contrib_init(_c = qmk_dir(r00t, "!c"));
	ClipNode **kcp = ClipNode::m0_kcp; kcp[0] = kcp[1] = kcp[2] = qmk_clip(_b, "[]");
	_b->set_perm(DF_ALL, DF_HROOT); _b_help->set_perm(DF_HROOT, DF_HROOT);
	_c->set_perm(DF_ALL, 0);         (*kcp)->set_perm(DF_ALL, DF_ALL);
}

void nd_init() { Node::init(); }
