#include "util.h"
#include "combo.h"
#include "guistub.h"
#include "cmd.h"
#include "glob.h"
#include "util2.h"

class CFIBoxModel : public BoxModel {
        public:
                CFIBoxModel(ConStore *pc, char *q, int ff) :
			pcs(pc), foif(ff), fdsc((unsigned char*)q), cdsc((unsigned char*)(q+(ff&31))) {}
                virtual BoxInst * mk_box();
                ModelPtr fm, cm;
                ConStore * pcs;
                int ivck_bv[4], foif, rsrv;
                unsigned char xdsc[2], *fdsc, *cdsc;
};

class CFIBoxInst : public BoxInst {
        public:
		static scf_t sc_f0, sc_f1, sc_fa0, sc_fa1;
                CFIBoxInst(CFIBoxModel * m) : BoxInst(sc_f0), m_nr(-1), m_m(m), m_blk(0) { BoxModel::ref(m); }
                ~CFIBoxInst();
        protected:
                int ini(double** inb);
		int prep_in(int inflg, double** inb);
                void fill_con(double**to, int bv, unsigned char * ix) {
                        ConStore * pcs = m_m->pcs; BVFOR_JM(bv) to[j] = pcs->p(ix[j]&63); }
                int m_nr, m_inc_bv;
                CFIBoxModel * m_m;
                double ** m_in1;
                BoxInst ** m_ppbx;
		char * m_blk;
};

class CFIBoxGen : public BoxGen {
        public: 
                friend void itb_init();
                CFIBoxGen(ABoxNode * nd);
                virtual int n_in() const { return m_ni; }
                virtual int n_out() const { return 1; }
                virtual const char * cl_name() { return "iter"; }
                virtual bool io_alias_perm() const { return false; }
		virtual int ifflg() const { return BIF_GC; }
                virtual int save2(SvArg * sv);
                virtual void notify_nio(BoxGen * bx);
                virtual void box_window() { w(-1); }
                virtual void spec_debug();
        protected:
		static void in_chk(unsigned char *p, int n, int ty, int lim);
                BXCMD_DECL(CFIBoxGen) c_b09, c_s09, c_r09, c_ni, c_ck, c_gc, c_gf, c_gx, c_sc, c_sf,
				      c_cnn, c_val, c_ui;
                virtual void set_mdl();
		inline int cf() const { return 2*!!m_fbx + !!m_cbx; }
		inline int w2l(int i) { return x2x(i,  0); }
		inline int l2w(int i) { return x2x(i, 16); }
		int save_ioc(AOBuf *f, int i, int j, int v);
		int conn(int ii, int ij, int oi, int oj);
		int scon(int ii, int ij, double v, int f);
		int orng(int i, int j); // 0:ok
		int set_ni(int n);
		void set_nif(int n), set_nic(int n), set_noc(int n);
                void w(int flg);
		void w_gr(), w_i1(int x), w_gn(int ix, BoxGen*bx, int ni, int no, unsigned char *dsc);
		int x2x(int i, int k);
		int gr_el(unsigned char *to);
                int set_fb(BoxGen * bx), set_cb(BoxGen * bx);
		int set_fb_2(BoxGen * bx, int ni);
		void rm_fb();
		int upd09();
		int iosel(int i, int j);
		unsigned char * ijp(int i, int j);
		void au_rgb(), au_lbl();
		int cf_win(int j) { return ((m_cbx&&!j) ? m_cbx : m_fbx?m_fbx:this)->node()->draw_window(16); }
                int foif() const { log("foif: m_nif=%d m_nic=%d m_noc=%d", m_nif, m_nic, m_noc);
			return m_fbx ? (m_fbx->io_alias_perm()<<30) + (m_flg<<15) + 1024*m_noc + 32*m_nic+m_nif
				     : (m_flg<<15) + 1; }
		void dsc_debug(const char * head, const unsigned char *p, int n);
                BoxGen *m_fbx, *m_cbx;
                unsigned char m_ni, m_flg, m_nif, m_nic, m_noc, m_f09, m_dseq,
			      m_sel[2], m_xin[2], m_fin[30], m_cin[30];
		double m_wnum;
                ConStore m_cs;
};

BoxInst * CFIBoxModel::mk_box() { return new CFIBoxInst(this); }
CFIBoxInst::~CFIBoxInst() { BoxModel::unref(m_m); for (int i=0;i<m_nr;i++) delete (m_ppbx[i]); free(m_blk); }

int CFIBoxInst::ini(double** inb) {
        CFIBoxModel *mo=m_m; if (!mo->fm.nz()) return m_nr = 0;
        int nri=mo->xdsc[0], nr=m_nr = ivlim((int)lround((nri&64) ? mo->pcs->v(nri&63) : inb[nri][0]), 0,1024);
        if (!nr) return 0;
	int rv = (nr==1) ? 1 : (2+(mo->foif>>30));
        int ff = mo->foif, nfi = ff&31, nco = (ff>>10)&31, sz_pco = nr*nco*sizeof(double);
        char *q = m_blk = (char*)malloc((nr+nfi)*sizeof(void*) + sz_pco);
        double *pco   = (double*)  q; q += sz_pco;
        m_ppbx =        (BoxInst**)q; q += nr *sizeof(void*);
        double **in1 = m_in1 = (double**) q;
        fill_con(in1, mo->ivck_bv[1], mo->fdsc);
	mo->fm.mk_boxv(m_ppbx, nr);
        if (!mo->cm.nz()) return m_inc_bv = 0, rv;
        int nci = (ff>>5)&31; unsigned char *cdsc = mo->cdsc, *fdsc = mo->fdsc;
        double *tci[nci]; fill_con(tci, mo->ivck_bv[3], cdsc);
        BVFOR_JM(((1<<nci)-1)^mo->ivck_bv[3]) tci[j] = inb[(int)cdsc[j]];
        double *ppco[nco]; for (int j=0; j<nco; j++) ppco[j] = junkbuf;
	IFDBGX(EXPR) log("ff=0x%x nco=%d bv[2]=0x%x", ff, nco, mo->ivck_bv[2]);
        BVFOR_JM(mo->ivck_bv[2]) { int k = fdsc[j]&31; ppco[k] = in1[j] = pco+k*nr; 
				   if (k>=nco) bug("k(%d)>=nco(%d)",k,nco);
				   IFDBGX(EXPR) log("bv[2]: j=%d k=%d",j,k); }
        BoxInst * cb = mo->cm.mk_box();
        int ibv = 0, cret = cb -> calc(0, tci, ppco, nr); delete(cb); if (cret<0) return cret;
        BVFOR_JM(mo->ivck_bv[2]) { int k = fdsc[j]&31; ibv |= ((cret>>k)&1)<<j; }
        return m_inc_bv = ibv, rv;
}

int CFIBoxInst::prep_in(int inflg, double** inb) {
	CFIBoxModel *mo = m_m; int rv = 0; double **in1 = m_in1; unsigned char *dsc = mo->fdsc;
	in1[0]=inb[0]; BVFOR_JM(mo->ivck_bv[0]) { int k=dsc[j]; in1[j]=inb[k], rv|=(((inflg>>k)&1)<<j); }
	return rv;
}

BX_SCALC(CFIBoxInst::sc_f0) {
	static sc_t sctab[4] = { sc_cp0, sc_f1, sc_fa0, sc_fa1 };
	SCALC_BXI(CFIBoxInst); int r = bxi->ini(inb); if (r<0) return r; CALC_FW(sctab[r]); }

#define CFI_C0  SCALC_BXI(CFIBoxInst); \
		double **in1 = bxi->m_in1; int of = inflg&1, f1_bv = bxi->prep_in(inflg, inb)

BX_SCALC(CFIBoxInst::sc_f1) { CFI_C0; return bxi->m_ppbx[0] -> calc(f1_bv|of, in1, outb, n); }

#define CCALC1(J,Q) if ((of = bxi->m_ppbx[J]->calc(f1_bv|of, in1, Q, n)) < 0) return of; else of &= 1

BX_SCALC(CFIBoxInst::sc_fa1) {
	CFI_C0; int incbv = bxi->m_inc_bv, nr = bxi->m_nr;
	CCALC1(0, outb); in1[0] = outb[0];
	for (int k=1; k<nr; k++) { BVFOR_JM(incbv) ++in1[j]; CCALC1(k, outb); }
	--nr; BVFOR_JM(incbv) in1[j] -= nr; return of;
}

BX_SCALC(CFIBoxInst::sc_fa0) {
	CFI_C0; int incbv = bxi->m_inc_bv, nr = bxi->m_nr, j0 = nr&1; double tmp[n], *pt = tmp;
        if(j0) {           CCALC1(0,outb);                   in1[0]=*outb; BVFOR_JM(incbv) ++in1[j]; }
        for (int k=j0;;) { CCALC1(k, &pt); ++k;              in1[0]=tmp;   BVFOR_JM(incbv) ++in1[j];
                           CCALC1(k,outb); if(++k>=nr)break; in1[0]=*outb; BVFOR_JM(incbv) ++in1[j]; }
	--nr; BVFOR_JM(incbv) in1[j] -= nr; return of;
}

CFIBoxGen::CFIBoxGen(ABoxNode *nd) : BoxGen(nd), m_fbx(0),m_cbx(0), m_ni(1),m_nif(1),m_nic(0),m_noc(0),
	m_f09(0), m_dseq(0), m_wnum(0.0) { m_xin[0]=70; }

void CFIBoxGen::set_mdl() {
	int nic = m_nic, nif = m_nif;
	CFIBoxModel * mo = m_mdlp.mkc1<CFIBoxModel>(&m_cs, nif+nic, foif());
	int *ivck = mo->ivck_bv;  mo->xdsc[0] = m_xin[0];unsigned char *q = mo->fdsc;
	memset(ivck,0,4*sizeof(int)); for (int j=1; j<nif; j++) ivck[(q[j]=m_fin[j])>>6] |= 1<<j;
	if (ivck[3]) bug("cfi/setm: nid=0x%x, bug_bv=0x%x", id(), ivck[3]);
	q = mo->cdsc; for (int j=0; j<nic; j++) ivck[3] |= ((q[j]=m_cin[j])>>6) << j;
	if (m_fbx) m_fbx->mdl_cpto(&mo->fm);
	if (m_cbx) m_cbx->mdl_cpto(&mo->cm);
}

int CFIBoxGen::save_ioc(AOBuf *f, int i, int j, int v) {
	return ((v&192)==64) ? f->pf("X$=%c%c%.15g\n", i+48, i_to_b32(j), m_cs.v(v&63))
			     : f->pf("X$+%c%c%c%c\n",  i+48, i_to_b32(j), 48+(v>>7), i_to_b32(v&31)); }

int CFIBoxGen::save2(SvArg * sv) {
	BXSV2_HEAD;
	CHKERR(f->pf("X$i%d\n", m_ni));
	if (m_cbx) { CHKERR(f->sn("X$c", 3)); CHKERR(m_cbx->node()->sv_path(10)); }
	if (m_fbx) { CHKERR(f->sn("X$f", 3)); CHKERR(m_fbx->node()->sv_path(10)); }
	CHKERR(save_ioc(f, 3, 1, m_xin[0]));
	for (int i=0; i<m_nic; i++) { CHKERR(save_ioc(f, 1, i, m_cin[i])); }
	for (int i=1; i<m_nif; i++) { CHKERR(save_ioc(f, 2, i, m_fin[i])); }
	return r;
}

void CFIBoxGen::in_chk(unsigned char *p, int n, int ty, int lim) {
	for (int x,i=0; i<n; i++) if (((x=p[i])>>6)==ty && (x&63)>=lim) p[i]=69; }

int CFIBoxGen::set_ni(int n) {
	return (n>=m_ni) ? ((n>m_ni) ? (unset_model(), m_ni=n, BCR_NIO|BCR_WIN) : 0)
			 : (unset_model(), in_chk(m_fin, m_nif, 0, m_ni=n),
					   in_chk(m_cin, m_nic, 2, n), BCR_NIO|BCR_WIN); }

void CFIBoxGen::set_noc(int n) { if (n<m_noc) in_chk(m_fin, m_nif, 0, n);   m_noc=n; }
void CFIBoxGen::set_nif(int n) { for (int i=m_nif; i<n; i++) m_fin[i]=69;   m_nif=n; }
void CFIBoxGen::set_nic(int n) { for (int i=m_nic; i<n; i++) m_cin[i]=69;   m_nic=n; }
void CFIBoxGen::rm_fb() { m_nif=1; set_boxp(&m_fbx,0); unset_model(); w(511); gui_errq_add(BXE_FILTRM,"cfi"); }

int CFIBoxGen::set_fb_2(BoxGen *bx, int ni) {
	if (bx!=m_fbx) { int ec = set_boxp(&m_fbx,bx); if (ec<0) return ec; }
	else if (ni==m_nif) { return 0; }
	unset_model(); set_nif(ni); w(511); return 0; }

int CFIBoxGen::set_fb(BoxGen *bx) {
	if (!bx) return set_fb_2(0, 1);   
	int ni = bx->n_in(), e = !ni | (bx->n_out()-1); 
	return e ? (bx==m_fbx ? (rm_fb(),0) : BXE_FILTNIO) : set_fb_2(bx, ni); }

int CFIBoxGen::set_cb(BoxGen *bx) {
	int ec, ni, no;
	if (bx==m_cbx) { if (bx && (((ni=bx->n_in())-m_nic)|(no=bx->n_out()-m_noc))) goto chg; else return 0; }
	if ((ec = set_boxp(&m_cbx,bx)) < 0) return ec;
	if (bx) ni=bx->n_in(), no=bx->n_out(); else ni=no=0;
chg:	return unset_model(), set_nic(ni), set_noc(no), BCR_WIN;
}

int CFIBoxGen::upd09() {
	int nxi, ec, sc = m_f09 & 15, nif = m_nif;
	IFDBGX(EXPR) log("upd09: %s -- f9=0x%x, nif=%d", m_node->path255(), m_f09, nif);
	if (!sc) { m_ni=nif+1; m_xin[0]=1; for (int i=1; i<nif; i++) m_fin[i] = i+1; goto done; }
	if ((ec=set_cb(box_bookmark[8+(m_f09>>4)])) < 0) return ec;
	nxi = 3+(sc==1);
	m_xin[0] = 1; for (int i=1; i<nif-1; i++) m_fin[i] = i+nxi; m_fin[nif-1] = 128;
	for (int i=0; i<3; i++) m_cin[i] = i+1; m_cin[3] = sc + (sc==1 ? 3 : 64);
	m_ni = nif + nxi - 1; 
done:   unset_model(); au_rgb(); au_lbl(); return BCR_NIO|BCR_WIN;
}

void CFIBoxGen::notify_nio(BoxGen * bx) {
	int ec; if (bx==m_fbx && (ec=set_fb(bx))<0) gui_errq_add(ec, "cfi/nio");
	if (bx==m_cbx) set_cb(bx); }

void CFIBoxGen::w_i1(int x) {
	static const int lbl[4] = {0xc7fc0, 0, 0xc33a0, 0xc3e40};
	int x6=x&63, x2=x>>6; if (x2==1) gui2.hdbl(m_cs.v(x6)); else gui2.b32n(lbl[x2]+x6,4); }

int CFIBoxGen::gr_el(unsigned char *to) {
	unsigned char *q = to; int x, x2, k = 64*(1+!!m_cbx), k2 = k+64*!!m_fbx;
	if ((x2=(x=m_xin[0])>>6)!=1) *q = x&31, q[1]=k2+1, q+=2;
	if (m_cbx) { for (int i=0;i<m_nic;i++) if ((x2=(x=m_cin[i])>>6)!=1) *q=      (x&31), q[1]=64+i, q+=2;}
	if (m_fbx) { for (int i=1;i<m_nif;i++) if ((x2=(x=m_fin[i])>>6)!=1) *q=32*x2+(x&31), q[1]=k +i, q+=2;}
	return ((q-to)>>1) + ( m_fbx ? (q[0]=0, q[1]=q[2]=k, q[3] = k2, 2) : (q[0]=0, q[1]=k2, 1) );
}

void CFIBoxGen::w_gn(int ix, BoxGen*bx, int ni, int no, unsigned char *dsc) {
	ABoxNode * bnd = bx->node(); gui2.gn_start(ix, bx, ni, no); 
	for (int j=0; j<ni; j++) w_i1(dsc[j]), gui2.ionm(bnd, 0, j), gui2.c1(36);
	for (int j=0; j<no; j++) gui2.ionm(bnd, 1, j), gui2.c1(36);   }

void CFIBoxGen::w_gr() {
	int sq = (++m_dseq) & 255, nb0 = !!m_fbx + !!m_cbx, ni = m_ni;
	unsigned char el[128]; int ecnt = gr_el(el);
	gui2.wupd_0('g', "i"); gui2.hexn(sq, 2); gui2.c4('_',48,48,50+nb0);
	gui2.c3('_',48,48); gui2.hexn(ecnt, 2);
	gui2.gn_start(-1, 0, 0, ni); for (int j=0; j<ni; j++) gui2.ionm(m_node, 0, j), gui2.c1(36);
	if (m_cbx) w_gn(0,       m_cbx, m_nic, m_noc, m_cin);
	if (m_fbx) w_gn(!!m_cbx, m_fbx, m_nif, 1,     m_fin);
	gui2.gn_start(nb0, 0, 2, 0); w_i1(m_fbx ? 192 : 0); gui2.sn("out$",  4);
				     w_i1(m_xin[0]); 	    gui2.sn("(#b)$", 5);
	gui2.wupd_0('g', "z");
	Dot * to = Dot::sg();
	if (!to->active()) return log("gr/to_dot: error: dot not running!");
	to->ghead(16*m_node->id()+11, 'i', 'g', -1, sq, nb0+2, ecnt);
	to->gnode(0, 0, m_ni, -2);
	if (m_cbx) to->gnode(1,		m_nic, m_noc, text_wid16(m_cbx->node()->s_name()) >> 5);
	if (m_fbx) to->gnode(1+!!m_cbx, m_nif, 1,     text_wid16(m_fbx->node()->s_name()) >> 5);
	to->gnode(nb0+1, 2, 0, -1);
	for (int i=0; i<ecnt; i++) to->gedge(el[2*i]>>6, el[2*i]&31, el[2*i+1]>>6, el[2*i+1]&31);
	to->gend();
}

int CFIBoxGen::iosel(int i0, int j0) {
	int i = (i0-47) &- (i0!='v'), j = (j0-=48)&31, k = (j0>>5)&1;
	gui2.setwin(w_oid(), 'i'); m_sel[k] = 64*i+j;
	if (!k) { unsigned char*p = ijp(i,j);  if (!p) return EEE_NOEFF;
		  int x=*p; if ((x&192)==64) gui2.wupd_d('x', m_wnum = m_cs.v(x&63)); }
	gui2.wupd_c0('g', '+'); gui2.c4(48+2*k, 48, 48+i, i_to_b32(j)); return 0; }

int CFIBoxGen::x2x(int i, int k) { static const signed char t[32] = {0,3,-1,-1, 0,1,3,-1, 0,2,3,-1, 0,1,2,3,
								     0,-1,-1,1, 0,1,-1,2, 0,-1,1,2, 0,1,2,3};
			    	   return t[4*cf()+(i&3)+k] |- !!(i&~3); }

int CFIBoxGen::orng(int i, int j) { switch(i) { case 0: return !j || j>=m_ni;
						case 1: return j>=m_noc;
						default:return 1; }}		

unsigned char * CFIBoxGen::ijp(int i, int j) { switch(i) { case 1: return j<m_nic    ? m_cin+j : 0;
							   case 2: return j&&j<m_nif ? m_fin+j : 0;
							   case 3: return j==1	     ? m_xin   : 0;
							   default : return 0; }}
int CFIBoxGen::conn(int ii, int ij, int oi, int oj) {
	if (ii<=oi || (ii==3 && oi)) return BXE_ITINV;
	unsigned char * q = ijp(ii, ij); if (!q || orng(oi, oj)) return BXE_ITNONE;
	int o8 = 64*(oi+!!oi) + oj; if (*q==o8) return 0;
	*q = o8; if (wnfl()) w(2);
	return 0;
}

int CFIBoxGen::scon(int ii, int ij, double v, int f) {
	unsigned char * q = ijp(ii, ij); if (!q) return BXE_ITNONE;
	int wf = wnfl(), cf = ((*q&192)==64);  if (!(f|cf)) return 0;
	if (cf) { m_cs.rm(*q&63); *q=64+m_cs.add(v); if (wf) gui2.setwin(w_oid(),'i'), gui2.grc(ii, ij, v); }
	else	{ 		  *q=64+m_cs.add(v); if (wf) w(2); }
	return 0;
}

void CFIBoxGen::w(int flg) {
	int wf = m_node->winflg(); if (!((wf|flg)&2048)) return;
	if (wf&2048) gui2.setwin(w_oid(), 'i'); else gui2.cre(w_oid(), 'i');
	if (flg&1) gui2.own_title();
	if (flg&2) w_gr(), m_sel[0]=m_sel[1]=255;
	if (flg&4) gui2.wupd_d('x', m_wnum);
	if (flg&8) gui2.wupd_i('i', m_ni);
}

void CFIBoxGen::au_rgb() {
	if (m_fbx && (m_node->set_ui_rgb(m_fbx->v_rgb()), wnfl())) gui2.setwin(w_oid(),'i'), gui2.own_title();}

#define COND_LBL(K, S) if ((k=(K))<30) m_node->set_ui_lbl(0, k, (S))
void CFIBoxGen::au_lbl() {
	char nm[8]; int k;
	for (int j=1; j<m_nif; j++) COND_LBL(m_fin[j], (nm[m_fbx->get_ionm(nm, 0, j)]=0, nm)); 
	for (int j=0; j<m_nic; j++) COND_LBL(m_cin[j], (nm[m_cbx->get_ionm(nm, 0, j)]=0, nm)); 
	COND_LBL(m_xin[0], "(#b)");
	if (wnfl()) w(2);
}

void CFIBoxGen::dsc_debug(const char * head, const unsigned char *p, int n) {
	char buf[1024],*q=buf; for (int x,i=0; i<n; i++) { switch ((x=p[i])>>6) {
		case 0:  q += sprintf(q, " i%d", x&63); break;
		case 1:  q += sprintf(q, " %g", m_cs.v(x&63)); break;
		case 2:  q += sprintf(q, " c%d", x&63); break;
		default: q += sprintf(q, " BUG%d", x); break; }}
	*q=0; log("%s%s", head, buf); }

void CFIBoxGen::spec_debug() {
	log("filt: %s", m_fbx ? m_fbx->node()->path255() : "(none)");
	log("ctrl: %s", m_cbx ? m_cbx->node()->path255() : "(none)");
	dsc_debug("#,f-inp :", m_fin, m_nif);
	dsc_debug("c-inp   :", m_cin, m_nic);
}

#define CH(X) BXCMD_H(CFIBoxGen, X)
CH(b09){BoxGen * bx = cb->lookup_b(s+1); if (!bx) return BXE_ARGLU;
	int ec = p->set_fb(bx); return ec<0 ? ec : p->upd09(); }

CH(s09){p->m_f09 &= ~15; p->m_f09 |= ivlim(s[1]-48, 0, 8); return p->upd09(); }
CH(r09){return (((p->m_f09>>4)^s[1])&1) ? (p->m_f09 ^= 16, p->upd09()) : 0;     }
CH(gc){ SEQCHK; BoxGen * bx = cb->lookup_b(s+3); return bx ? p->set_cb(bx) : BXE_ARGLU; } 
CH(sc){         BoxGen * bx = cb->lookup_b(s+1); return bx ? p->set_cb(bx) : BXE_ARGLU; } 
CH(gf){ SEQCHK; BoxGen * bx = cb->lookup_b(s+3); return bx ? p->set_fb(bx) : BXE_ARGLU; } 
CH(sf){         BoxGen * bx = cb->lookup_b(s+1); return bx ? p->set_fb(bx) : BXE_ARGLU; } 
CH(ni){ int k = p->m_ni; intv_cmd(&k, s+1, 1, 30); return p->set_ni(k); }
CH(cnn){return p->conn(s[1]&3, b32_to_i(s[2]), s[3]&3, b32_to_i(s[4])); }
CH(val){return p->scon(s[1]&3, b32_to_i(s[2]), atof(s+3), 1); }
CH(ui){ if (s[1]&1) p->au_rgb(); if (s[1]&2) p->au_lbl();  return 0; }
	
CH(ck){ SEQCHK; int k = 2*(*s-'1') + (s[5]==42);
	return !k ? p->iosel(s[4], s[5]) : (k==5) ? p->cf_win(s[4]-48) : EEE_NOEFF; }

CH(gx){ SEQCHK; int i8=p->m_sel[0], o8=p->m_sel[1], ii=i8>>6, ij=i8&31; switch(s[3]) {
	case 'v': return (i8==255) 	? BXE_GIOSEL : p->scon(ii, ij, p->m_wnum=atof(s+4), 0);
	case '=': return (i8==255)      ? BXE_GIOSEL : p->scon(ii, ij, p->m_wnum, 1);
	case '+': return ((i8|o8)==255) ? BXE_GIOSEL : p->conn(ii, ij, o8>>6, o8&31);
	default:  return BXE_CENUM; }}

BXCMD_DEF(CFIBoxGen) { {8192+'\\',0}, {'b',c_b09}, {'s',c_s09}, {'r',c_r09}, {'1'|256,c_ck}, {'3'|256,c_ck},
	{'G',c_gx}, {'C',c_gc}, {'F',c_gf}, {'i',c_ni}, {'+',c_cnn}, {'=',c_val}, {'c',c_sc}, {'f',c_sf},
	{'U',c_ui}, {0,0} };

int setbox_cfi(ABoxNode * nd, BoxGen * _) { nd->m_box = new CFIBoxGen(nd); return 1; }
void itb_init() { CFIBoxGen::cmd_init(); }
