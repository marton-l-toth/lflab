#include "box.h"
#include "mx.h"
#include "cmd.h"
#include "guistub.h"
#include "asnd.h"
#include "util2.h"

#define TR_SEL_ID (bxw_rawptr[0].i[1])
#define TR_SEL_XY (bxw_rawptr[1].i)
#define TR_V_XY01 (bxw_rawptr[2].s)
#define TR_UPP    (bxw_rawptr[3].s[0])

BoxGen * trk_rec_trg = 0;
static long long rec_t0;
static int rec_line, rec_x0, rec_mode;
static double rec_rate;
static ANode * trgp_node = 0; static int trgp_i, trgp_j;

struct BKMK128 { int n, rsrv; unsigned char * lo; unsigned short * hi[4]; unsigned int bv[4];
                 int find(int i), set(int i, int v); void del(); };
struct BKMK32k { unsigned long long bv[4]; BKMK128 *** p[4];
                 int set(int i, int v), find(int i); void del8k(int i), del(); };
struct BKMK1M { unsigned int bv; BKMK32k * p[32];  BKMK1M() : bv(0) {}
                 int find(int i); void set(int i, int v), del(); };


class TrackModel;
class TrackInst : public BoxInst {
	public:
		TrackInst (TrackModel * m);
		virtual ~TrackInst();
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int ini(double **inb);
		void mxprep(int bflg, double* bpm, int n);
		void sane(int n);

		TrackModel * m_m;
		ANode *m_pn, *m_fn, *m_tn;
		int m_mxid, m_to, m_to2, m_fr2, m_rpc, m_md;
		double m_vt;
		
};

class TrackModel : public BoxModel {
	public:
		TrackModel(ANode * t_n, ANode * g_0, ANode * g_1) : BoxModel(2), tn(t_n), g0(g_0), g1(g_1) {}
		virtual BoxInst * mk_box() { return new TrackInst(this); }
		ANode *tn, *g0, *g1;
		BKMK1M bkm;
};

class TrackGen : public BoxGen {
	public:
		friend void track_init();
		TrackGen(ABoxNode * nd, ANode * g0, ANode * g1);
		virtual int n_in () const { return 6; }
		virtual int n_out() const { return 2; }
		virtual bool io_alias_perm() const { return true; }
		virtual const char* cl_name() { return "trk"; }
		virtual void set_model() { bug("track/setmodel called"); }
		virtual void box_window();
		virtual void wdat_cons(sthg * bxw_rawptr);
		virtual int save2(SvArg * sv);
		virtual void mxc_notify(int k,int f) { if(f&64) m_mxctl=0; if (wnfl()) w_ply(); }
		ANode* bkm_find(int j);
		void   bkm_add(ANode * nd);
		void   bkm_rm (ANode * nd);
		int cond_pm(ANode* nd, int pm);
	protected:
		void grab_gsel(sthg *bxw_rawptr){ trgp_node=m_node; trgp_j=TR_SEL_XY[0]; trgp_i=TR_SEL_XY[1]; }
		void sel0w(sthg *bxw_rawptr, int i, int j) {
			TR_SEL_XY[0]=j; TR_SEL_XY[1]=i; grab_gsel(bxw_rawptr); }
		void sel0wn(sthg * bxw_rawptr, ANode * nd) {
			TR_SEL_ID=nd->id(); sel0w(bxw_rawptr, nd->cth()->i, nd->cth()->j); }
 		BXCMD_DECL(TrackGen) c_rq, c_cx0, c_cx1, c_view, c_mv, c_stp, c_rec, c_upp, c_gcf, c_bpm, c_pl,
				     c_qk, c_cky;
		int draw_bx_ini(int x, int ybv);
		int draw_bx_1();
		void w1b_pm(ANode * wb, int pm);
		int cx_cmd(sthg * bxw_rawptr, CmdBuf * cb, int c, int id, int x, int y);
		int w_sel(sthg * bxw_rawptr);
		int start_rec();
		int stop_rec() { trk_rec_trg = 0; if (wnfl()) w_rec(); return 0; }
		int gui_h4(int * to);
		int play(int op, int pos);
		void w_bpm() { gui2.setwin(w_oid(), 't'); gui2.wupd_i('b', m_bp10m); }
		void w_ply() { gui2.setwin(w_oid(), 't'); gui2.wupd_i1('p', !!m_mxctl); }
		void w_rec() { gui2.setwin(w_oid(), 't'); gui2.wupd_i1('R', trk_rec_trg==this); }
		TrackModel m_m;
		ANode * m_drq_g;
		int m_drq_x;
		unsigned int m_drq_msk;
		int m_bp10m, m_mxctl;
		unsigned char m_gwfr[4], m_div[256];
};

int BKMK128::find(int i) {
        int i2 = i>>5, i5 = i&31;
        if (bv[i2]) BVFOR_JM(bv[i2]>>i5) return i+=j, 16*hi[i2][i5+j]+((lo[i>>1]>>4*(i&1))&15);
        for (i2++; i2<4; i2++) {
                BVFOR_JM(bv[i2]) return i=32*i2+j, 16*hi[i2][j]+((lo[i>>1]>>4*(i&1))&15); }
        return 0;
}

int BKMK128::set(int i, int v) {
        int i2 = i>>5, i5 = i&31;
        unsigned int msk = 1u << i5, *q = bv+i2;
        if (!v) return (*q & msk) ? (*q&=~msk, --n) : n;
        if (!*q) { ++n; *q = msk;  if (!hi[i2]) hi[i2] = (unsigned short*) ANode::a64();
                                   if (!lo) lo = (unsigned char*)ANode::a64(); }
        else if (!(*q & msk)) { ++n; *q |= msk; }
        hi[i2][i5] = (unsigned short)(v>>4);
        unsigned char *p = lo+(i>>1);
        return i&=1, i*=4, *p&=(240>>i), *p|=(v&15)<<i, n;
}
void BKMK128::del() { ANode::f64c(lo); for (int i=0; i<4; i++) ANode::f64c(hi[i]); ANode::f64(this); }

int BKMK32k::find(int i) {
        int k, i7 = i&127, i6 = (i>>=7)&63, i2 = (i>>6)&3;
        unsigned long long x = bv[i2] >> i6;
        if ((x&1) && (k = p[i2][i6>>3][i6&7]->find(i7))) return k;
        ++i6; x>>=1;
        while (1) {
                BVFOR_JM((unsigned int)(x&0xffffffffu)) return i6+=j   , p[i2][i6>>3][i6&7]->find(0);
                BVFOR_JM((unsigned int)(  x >> 32    )) return i6+=j+32, p[i2][i6>>3][i6&7]->find(0);
                if (i2==3) return 0; else i6 = 0, x = bv[++i2];
        }}

void BKMK32k::del8k(int i2) { BKMK128 ***qqq = p[i2]; for (int i=0; i<8; i++) if (qqq[i]) {
                for (int j=0; j<8; j++) ANode::f64c(qqq[i][j]); ANode::f64(qqq[i]); }
                ANode::f64(qqq); bv[i2] = 0; p[i2] = 0; }
int BKMK32k::set(int i, int v) {
        int i7 = i&127, i6 = (i>>=7)&63, i2 = (i>>6)&3;
        unsigned long long msk = 1ull << i6;
        BKMK128 ***qqq, **qq, *q = (bv[i2]&msk) ? p[i2][i6>>3][i6&7] : 0;
        if (!v) return !q || q->set(i7,0) || (q->del(), bv[i2]&=~msk) ||
                                        (del8k(i2), bv[i2^1]) || bv[i2^2] || bv[i2^3];
        if (q) return q->set(i7, v);
        qqq = (p[i2] ? p[i2] : (p[i2] = (BKMK128***) ANode::z64())) + (i6>>3);
        qq  = (*qqq  ? *qqq  : (*qqq  = (BKMK128**)  ANode::z64())) + (i6& 7);
        return bv[i2]|=msk,    (*qq   = (BKMK128*)   ANode::z64()) -> set(i7, v);
}
void BKMK32k::del() { for (int i=0; i<4; i++) if (p[i]) BKMK32k::del8k(i); ANode::f64(this); }

int BKMK1M::find(int i) {
        int k, i5 = i>>15, i15 = i&32767;
        unsigned int x = bv >> i5;
        if ((x&1) && (k = p[i5]->find(i15))) return k;
        x >>= 1; ++i5;
        BVFOR_JM(x) return p[i5+j]->find(0);
        return 0;
}

void BKMK1M::del() { BVFOR_JM(bv) p[j]->del(); }
void BKMK1M::set(int i, int v) {
        int i5 = i>>15, i15 = i&32767;
        unsigned int msk = 1u<<i5;
        if (!v) { if ((bv&msk) && !p[i5]->set(i15,0)) bv&=~msk, p[i5]->del();   return; }
        ((bv&msk) ? p[i5] : (bv |= msk, p[i5] = (BKMK32k*)ANode::z64())) -> set(i15, v);
}

TrackInst::TrackInst (TrackModel * m) : m_m(m), m_pn(0), m_fn(0), m_tn(0), m_md(0) {
	log("trki: hello");
	BoxModel::ref(m); m_mxid=mx_mkroot(); }

TrackInst::~TrackInst() { if (m_mxid>0) mx_del(m_mxid); Node::del(m_pn,NOF_FORCE);
	log("trki: bye");
	Node::del(m_fn,NOF_FORCE); Node::del(m_tn,NOF_FORCE); BoxModel::unref(m_m); }

#define VTARG(x) ((int)lround(40320.0*inb[x][0]))
int TrackInst::ini(double **inb) {
	m_md = 1; m_rpc = (int)lround(inb[3][0]);
	m_vt = 40320.0*inb[1][0]; m_to = VTARG(2); 
	int vti = (int)lround(m_vt); if ((vti|m_to)<0) return BXE_RANGE;
	int ec = Node::mk(&m_pn, m_m->tn, 0, '!', 2|NOF_FORCE, vti);
	if (ec<0 || (ec = Node::mk(&m_tn, m_m->tn, 0, '!', 3|NOF_FORCE, m_to))<0) return ec;
	if (!m_rpc) return 0;
	m_fr2 = VTARG(4); m_to2 = VTARG(5);
	return Node::mk(&m_fn, m_m->tn, 0, '!', 1|NOF_FORCE, m_fr2);
}

void TrackInst::sane(int n) {
	int k = 9999; ANode *pv = m_m->g0, *nd = pv->next();
	while (1) {
		if (!--k) bug("trk/sane/loop(%d)", n);
		if (nd->cth()->pv!=pv) bug("trk/sane/link(%d)", n);
		if (nd->cth()->j < pv->cth()->j) bug("trk/sane/ord(%d)", n);
		if ((pv=nd, nd=nd->next())==m_m->g1) return; 
	}}

void TrackInst::mxprep(int bflg, double* bpm, int n) {
	ANode * nd = m_pn -> next();
	double vt = m_vt, vtl = (double)nd->cth()->j - 0.001;
	sane(2);
	for (int i=0; i<n; i++, vt+=*bpm*0.015238095238095238, bpm+=bflg) {
		sane(99);
		for (; vt>=vtl; nd=nd->next(), vtl = (double)nd->cth()->j - 0.001) {
			while (nd->cth()->i == 2) nd = nd->next();
			if (nd->cl_id()=='w') { wrap_2mx(static_cast<ABoxNode*>(nd)->box(), m_mxid, 0, i); 
				                if (debug_flags & DFLG_TRK) log("trk2mx: %d", m_mxid); }
			else if (nd==m_tn) {
				if (!m_rpc) { m_md = 3; return; }
				vt += (double) (m_fr2 - m_to);
				sane(0);
				if (m_md==1) m_md = 2, Node::move(m_tn, m_m->tn, 0, 7|NOF_FORCE, m_to=m_to2);
				sane(1);
				nd = m_fn;  }
			else if (nd->cth()->i == 4095) { m_md = 3; return; }
			else if (nd->cl_id()!='!') {
				return m_md = 4, log("BUG: trk: 0x%x, cl 0x%x", nd->id(), nd->cl_id()); }
		}}
	while (nd->cth()->i == 2) nd = nd->next();
	if (nd->cl_id()=='!' && nd->cth()->i == 4095) m_md = 3;
	else m_vt = vt, sane(299), m_pn->trk_insbf(nd, (int)lround(vt)+1), sane(3);
}

int TrackInst::calc(int inflg, double** inb, double** outb, int n) {
	//log("trk/calc: md=%d, n=%d", m_md, n);
	int ec;switch (m_md) {
		case 0: if ((ec = ini(inb))<0) return m_md = 4, ec;  sane(111);
		case 1: case 2: mxprep(inflg&1, *inb, n);
		case 3: break;
		case 4: outb[0][0] = outb[1][0] = 0.0; return 0;
		default: bug("trk: invalid m_md (%d)", m_md); return 0;
	}
	switch(ec = mx_calc(m_mxid, outb[0], outb[1], n, 0)) {
		case 0:	if (m_md==3 && mx_r_isemp(m_mxid)) m_md=4, mx_del(m_mxid), m_mxid=-1;
			outb[0][0] = outb[1][0] = 0.0; return 0;
		case 1: if (outb[1]!=junkbuf) memcpy(outb[1], outb[0], 8*n);
		case 2: return 3;
		default: return ec<0 ? ec : TKE_WTF;
	}
}

/////////////////////////////////////

TrackGen::TrackGen(ABoxNode *nd, ANode *g0, ANode *g1) : BoxGen(nd), m_m(nd,g0,g1), m_drq_g(0), 
	m_bp10m(600), m_mxctl(0) {
	m_model = &m_m; memcpy(m_gwfr, TRK_DEF_GWFR, 4); m_div[0] = 3; memset(m_div+1, 0, 255); }

int TrackGen::gui_h4(int * to) {
	int k, n = 0;
	for (int i=0; i<4; i++) if ((k=m_gwfr[i])!=TRK_DEF_GWFR[i]) to[n++] = qh4((i<<8)+k) + 47;
	if ((k=m_div[0]) != 3) to[n++] = qh4(k);
	for (int i=1; i<256; i++) if ((k=m_div[i])) to[n++] = qh4((i<<8)+k);
	return n;
}

void TrackGen::box_window() {
	int buf[260], n = gui_h4(buf);
        gui2.cre(w_oid(), 't'); gui2.c1('g'); gui2.sn((char*)buf, 4*n); gui2.c1('.');
	gui2.own_title(); w_sel(m_node->wdat_raw()); w_bpm(); w_ply(); w_rec();
}

void TrackGen::wdat_cons(sthg * bxw_rawptr) {
	TR_SEL_ID = TR_SEL_XY[0] = 0; TR_SEL_XY[1] = 16; TR_V_XY01[1] = TR_V_XY01[3] = -1; TR_UPP = 420; }

int TrackGen::w_sel(sthg * bxw_rawptr) {
	char s[20]; s[3] = '^';
	int nid = TR_SEL_ID, *sxy = TR_SEL_XY;
	*(int*)(s+ 4) = qh4(16*sxy[1]+(nid>>16));  *(int*)(s+ 8) = qh4(nid   &65535);
	*(int*)(s+12) = qh4(sxy[0]>>16);           *(int*)(s+16) = qh4(sxy[0]&65535);
	gui2.setwin(w_oid(), 't'); gui2.t_sn(s+3, 17); if (!nid) return 0;
	ANode * nd = ANode::lookup_n_q(nid); if (nd->cl_id()!='w') return BXE_WTF;
	gui2.bxmini(nd->box0()); return 0;
}

int TrackGen::cx_cmd(sthg * bxw_rawptr, CmdBuf * cb, int c, int id, int x, int y) {
	ANode * nd; if (!id || (nd = ANode::lookup_n_q(id))->cl_id()!='w') id = 0, nd = 0;
	int ec; switch(c) {
		case '1': return (TR_SEL_ID=id) ? sel0w(bxw_rawptr, nd->cth()->i, nd->cth()->j)
			  			: sel0w(bxw_rawptr,y, x),  w_sel(bxw_rawptr);
		case '4': return nd ? wrap_2mx(nd->box0(), 0) : EEE_NOEFF;
		case '6': return nd ? nd->draw_window(0x1b) : BXE_NOARG;
		case '9': return nd ? nd->draw_window(0x19) : BXE_NOARG;
		case 'C': return TR_SEL_ID ? Node::copy(ANode::lookup_n_q(TR_SEL_ID), m_node, 0, y|(cb->cnof()&~NOF_FGUI), x) : EEE_NOEFF;
		case 'M': ec = (id=TR_SEL_ID) ? Node::move(ANode::lookup_n_q(TR_SEL_ID), m_node, 0, y|(cb->cnof()&~NOF_FGUI), x) : EEE_NOEFF;
			  if (ec>=0) sel0w(bxw_rawptr, y, x), TR_SEL_ID=id, w_sel(bxw_rawptr);    return ec;
		case 'N': return Node::mk(0, m_node, 0, 'w', y|cb->cnof(), x);
		case 'v': return (nd = ClipNode::kcp(2)->ent_sel()) ? Node::copy(nd, m_node, 0, y|(cb->cnof()&~NOF_FGUI), x) : EEE_NOEFF;
		case 'c': return nd ? Node::copy(nd, ClipNode::kcp(1), 0, NOF_NOIDX|cb->cnof()) : EEE_NOEFF;
		case 'x': return nd ? Node::move(nd, ClipNode::kcp(1), 0, NOF_NOIDX|cb->cnof()) : EEE_NOEFF;
		default: return BXE_CENUM;
}}

int TrackGen::play(int op, int pos) {
	if (!m_mxctl) { if (!op) return EEE_NOEFF; }
	else { mx_c_stop(m_mxctl,1,2); if (m_mxctl) return TKE_WTF; if (!op) return 0; }
	m_mxctl = mx_mkctl(this);
	double v[12]; v[0] = v[1] = 1.0; v[2] = v[5] = 0.0; v[3] = 10.0; v[4] = 1e-5;
	v[6] = 0.1 * (double)m_bp10m; v[9] = (op&1) ? 0.0 : 9999.0;
	v[7] = v[10] = (double)pos; v[8] = v[11] = (op==1) ? 53000.0 : (double)(pos+(op>>1));
	int bi = mx_add_box(0, m_m.mk_box(), "-A\x06\x02", v, 0x402); 
	return (bi<0) ? bi : mx_c_add(m_mxctl, bi, 1);
}

#define CH(X) BXCMD_H(TrackGen, X)

CH(cx0){BXW_GETP(p); return p->cx_cmd(bxw_rawptr, cb, s[1], TR_SEL_ID, TR_SEL_XY[0], TR_SEL_XY[1]); }
CH(upp){BXW_GETP(p); TR_UPP = atoi_h(s+1); return 0; }
CH(gcf){return trk_g_parse(s+1, p->m_div, p->m_gwfr), 0; }
CH(pl) {return s[1]>47 ? p->play(s[1]-48, atoi_h(s+2)) : BXE_PARSE; }
CH(qk) {BXW_GETP(p); return p->grab_gsel(bxw_rawptr), ClipNode::kcp(2)->keyop_f(atoi_h(s+1),5,0,cb->cnof()); }
CH(cky){int x = atoi_h(s+1); log("trk_qk: atoi_h(s+1) = %d", x); if (x!=53) return BXE_CENUM;
	BXW_GETP(p); if (!TR_SEL_ID) return EEE_NOEFF;
	ANode *nd0 = ANode::lookup_n_q(TR_SEL_ID), *nd = nd0->next();
	int ec = Node::move(nd0, ClipNode::kcp(1), 0, NOF_NOIDX|cb->cnof()); if (ec<0) return ec;
	while (nd->cth()->j<0x7fffffff && nd->cl_id()!='w') nd=nd->next();
	if (nd->cl_id()=='w') TR_SEL_ID=nd->id(), p->sel0wn(bxw_rawptr, nd), p->w_sel(bxw_rawptr);
	return ec;
}

CH(bpm){if (!intv_cmd(&p->m_bp10m,s+1,1,9999,0x64640a01)) return 0;
	if (p->wnfl()) p->w_bpm(); if (p->m_mxctl) mx_c_bpm_ugly_hack(p->m_mxctl, p->m_bp10m); return 0; }

CH(cx1){int id, x, y, c = s[1]; if (!c) return BXE_NOARG; else id = atoi_h(s+2);
	if (!(s=cb->tok())) return BXE_NOARG; else y = atoi_h(s);
	if (!(s=cb->tok())) return BXE_NOARG; else x = atoi_h(s);
	BXW_GETP(p); return p->cx_cmd(bxw_rawptr, cb, c, id, x, y);
}

CH(view){int i, a[4]; 
	if ((unsigned int)(a[2] = s[1]-48) > 31 ||
	    (unsigned int)(a[3] = s[2]-48) > 31   ) return BXE_PARSE;
	const char *s1, *s2; if (!(s1=cb->tok()) || !(s2=cb->tok())) return BXE_NOARG;
	if ( ((a[0]=atoi_h(s1)) | (a[1]=atoi_h(s2))) & 0xffffe000 ) return BXE_RANGE;
	BXW_GETP(p); for (i=0; i<4; i++) TR_V_XY01[i] = (short)a[i]; return 0;
}

CH(rq){	if (s[1]==',') return p->draw_bx_1();
	const char * s2 = cb->tok(); if (!s2) return BXE_NOARG;
	int ec = p->draw_bx_ini(atoi_h(s+1), atoi_h(s2)); return (ec<0) ? ec : p->draw_bx_1();
}

CH(stp){BXW_GETP(p); if (!TR_SEL_ID) return BXE_NOARG;
	ANode * nd = ANode::lookup_n_q(TR_SEL_ID);
	switch(s[1]) {
		case '>':
			for (nd=nd->next(); nd->cth()->j<0x7fffffff && nd->cl_id()!='w'; nd=nd->next());
			break;
		case '<':
			for (nd=nd->cth()->pv; nd->cth()->j>=0 && nd->cl_id()!='w'; nd=nd->cth()->pv);
			break;
		default: return BXE_CENUM;
	}
	if (nd->cl_id()=='w') p->sel0wn(bxw_rawptr, nd), p->w_sel(bxw_rawptr);
	return 0;
}

CH(mv){	if (!s[1] || !s[2]) return BXE_NOARG;
	BXW_GETP(p);
	int ec, x = TR_SEL_XY[0], y = TR_SEL_XY[1], id = TR_SEL_ID, v = ((0x64640a01>>(8*(s[3]&3)))&255) * ((s[2]&2)-1);
	ANode * nd = id ? ANode::lookup_n_q(id) : 0;
	switch (s[1]) {
		case 'l': y += 16*v; break;
		case 'b': x += 40320*v; break;
		case 'd': x += TR_UPP*hex2(s+4)*v; break;
		case 'p': x += TR_UPP*v; break;
		case 'u': x += v; break;
		default: return BXE_CENUM;
	}
	if (y<16) y = 16; else if (y>4080) y = 4080;
	if (x<0) x = 0; else if (x>2096640000) x = 2096640000;
	if(id && (ec = Node::move(ANode::lookup_n_q(id), p->m_node, 0, y|(cb->cnof()&~NOF_FGUI), x))<0) return ec;
	if (id) p->sel0wn(bxw_rawptr, nd); else p->sel0w(bxw_rawptr, y, x);
	p->w_sel(bxw_rawptr);
	return 0;
}

CH(rec){if (trk_rec_trg == (BoxGen*)p) return p->stop_rec();
	if (trk_rec_trg) static_cast<TrackGen*>(trk_rec_trg)->stop_rec();
	return p->start_rec();
}

int TrackGen::start_rec() {
	BXW_GET; trk_rec_trg = this; rec_mode = 0; rec_x0 = TR_SEL_XY[0]; rec_line = TR_SEL_XY[1];
	rec_rate = 67.2 * (double)m_bp10m * sample_length;
	if (wnfl()) w_rec(); return 0; }

void TrackGen::w1b_pm(ANode * wb, int pm) {
	const trk_24 * th = wb->cth(); int id = wb->id(); gui2.c1(pm); gui2.hex4(16*th->i+(id>>16));
	gui2.hex4(id&65535); gui2.hex8(th->j); if (pm&2) gui2.bxmini(wb->box0()); }

int TrackGen::cond_pm(ANode * wb, int pm) {
	if (pm=='+' && !wb->box0()) return wb->winflg_or(8), 0;
	BXW_GET; const trk_24 * th = wb->cth();
	if (pm=='-' && wb->id()==TR_SEL_ID) TR_SEL_ID = 0, w_sel(bxw_rawptr);
	log("trk/cond_pm: i=%d, rng: %d...%d j=%d, rng: %d...%d", th->i, 128*TR_V_XY01[2], 128*TR_V_XY01[3]+127,
			th->j, (TR_V_XY01[0]<<18)-20000, (TR_V_XY01[1]<<18)+262143);
	if (th->i<128*TR_V_XY01[2] || th->i>128*TR_V_XY01[3]+127) return 0;
	if (th->j<(TR_V_XY01[0]<<18)-20000 || th->j>(TR_V_XY01[1]<<18)+262143) return 0;
	gui2.setwin(w_oid(), 't'); gui2.wupd('t'); w1b_pm(wb, pm); return 1;
}

int TrackGen::draw_bx_ini(int x, int ybv) {
	if (m_drq_g) return EEE_STATE;
	if (!ybv) return EEE_NOEFF;
	int ec = Node::mk(&m_drq_g, m_node, 0, '!', 9|NOF_FORCE, (x+1)<<18);  if (ec<0) return ec;
	m_drq_x = x<<18; m_drq_msk = ybv; return 0;
}

int TrackGen::draw_bx_1() {
	if (!m_drq_g) return EEE_STATE;
	ANode * q = m_drq_g->cth()->pv;  int k = 0;
	gui2.setwin(w_oid(), 't'); gui2.wupd('t');
	while (1) {
		if (q->cth()->j < m_drq_x) goto done;
		if (q->cl_id()!='w') { ++k; q = q->cth()->pv; continue; }
		if (k>10000) goto cont;
		w1b_pm(q,'+'); k += 32; q = q->cth()->pv;
	}
done:	gui2.c1('*'); gui2.hex4(m_drq_x>>18); gui2.hex8((int)(m_drq_msk)); 
	Node::del(m_drq_g, NOF_FORCE); m_drq_g = 0; return 0;
cont:   gui2.c1(','); m_drq_g -> trk_insbf(q->next()); return 0;
}

ANode * TrackGen::bkm_find(int j) {
	if ((j+1)<2) return j ? (j<0 ? m_m.g0 : m_m.g1) : m_m.g0->next();
	int j20 = j>>12, ni = m_m.bkm.find(j20);
	if (!ni) return m_m.g1;
	ANode *q, *r = ANode::lookup_n_q(ni);
	while ((q = r->cth()->pv)->cth()->j >= j) r = q;
	return r;
}

void TrackGen::bkm_add(ANode * nd) {  const trk_24 * p = nd->cth();  int j = p->j;
	if (p->ty=='w' && ((p->pv->cth()->j ^ j) & 0xfffff000)) m_m.bkm.set(j>>12, nd->id()); }

void TrackGen::bkm_rm(ANode * nd) {   const trk_24 * p = nd->cth();  int j = p->j;
	if (p->ty=='w' && ((p->pv->cth()->j ^ j) & 0xfffff000)) m_m.bkm.set(j>>12,
			((j^(nd=nd->next())->cth()->j)&0xfffff000) ? 0 : nd->id());  }

int TrackGen::save2(SvArg * sv) {
	BXSV2_HEAD; int buf[300], n;
	buf[0] = 0x67245800; n = gui_h4(buf+1); buf[n+1] = 10; CHKERR(f->sn((char*)buf+1, 4*n+4));
	return r;
}

BXCMD_DEF(TrackGen) { {8192+'\\',0}, {'R'|256, c_rq}, {'C'|256, c_cx1}, {'c'|256, c_cx0},
	{'V'|256, c_view}, {'m'|256, c_mv}, {'s'|256, c_stp}, {'r', c_rec}, {'u'|256, c_upp},
	{'g'|256, c_gcf}, {'b', c_bpm}, {'p'|256, c_pl}, {'Q', c_qk}, {'K', c_cky}, {0, 0} };

//export
BoxGen * trk_mk(ABoxNode * nd, ANode * g0, ANode * g1) { return new TrackGen(nd,g0,g1); }
ANode * trk_bkm_find(BoxGen * abx, int j) { return static_cast<TrackGen*>(abx)->bkm_find(j); }
void trk_bkm_add(BoxGen * abx, ANode * nd) { static_cast<TrackGen*>(abx)->bkm_add(nd); }
void trk_bkm_rm(BoxGen * abx, ANode * nd) { static_cast<TrackGen*>(abx)->bkm_rm(nd); }
int trk_cond_pm(BoxGen * abx, ANode * nd, int pm) { return static_cast<TrackGen*>(abx)->cond_pm(nd, pm); }
int trk_glob_paste(BoxGen *bx, int nof) { return !trgp_node ? EEE_NOEFF : 
	Node::mk(0, trgp_node, 0, 'W', trgp_i|nof, trgp_j, bx); }

int trk_rec(ANode * nd, int mxbi) {
	ANode *cpnd; BoxGen *cpbx, *bx0 = static_cast<ABoxNode*>(nd)->box();
	int t = (rec_mode) ? (int)lround(rec_rate * (double)(snd0.total_played() - rec_t0))
			   : (rec_t0 = snd0.total_played(), rec_mode=1, 0),
	    r = Node::mk(&cpnd, trk_rec_trg->node(), 0,'W', rec_line, rec_x0+t, bx0);
	if (r<0 || (r = mx_tr_add(mxbi, cpbx=static_cast<ABoxNode*>(cpnd)->box()))<0) return r;
	wrap_set_trec(cpbx, r); return 0;
}

void track_init() { TrackGen::cmd_init(); }
