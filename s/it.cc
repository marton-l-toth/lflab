#include "util.h"
#include "combo.h"
#include "guistub.h"
#include "cmd.h"
#include "glob.h"
#include "util2.h"

#define IBF_ALI 0x40000000
#define IBF_REL 0x10000

class ItBoxInst : public BoxInst {
	public:
		ItBoxInst(ModelPtr b1m, int flg) : m_b1m(b1m), m_ppbx(0), m_zvec(0), m_flg(flg), m_n(-1) {}
		virtual ~ItBoxInst() { for (int i=0; i<m_n; i++) delete(m_ppbx[i]); 
				       free(m_ppbx); free(m_zvec); }
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int ini(double **par);
		ModelPtr m_b1m;
		BoxInst ** m_ppbx;
		double * m_zvec;
		int m_flg, m_n;
};

class ItBoxModel : public BoxModel {
        public: ItBoxModel(int _) {}
		virtual BoxInst * mk_box() { log("b1m/mk=%p", box1m.rawmp()); return new ItBoxInst(box1m, flg); }
                ModelPtr box1m;
		int flg;
};

class ItBoxGen : public BoxGen {
        public: 
                friend void itb_init();
                ItBoxGen(ABoxNode * nd) : BoxGen(nd), m_bx1(0), m_ni(1), m_scl(0), m_flg(0) {}
                virtual int n_in() const { return m_ni + (m_scl ? 2+(m_scl==1) : 1); }
                virtual int n_out() const { return 1; }
                virtual const char * cl_name() { return "iter"; }
		virtual bool io_alias_perm() const { return false; }
		virtual int save2(SvArg * sv);
		virtual void notify_nio(BoxGen * bx) { set_bx1(bx); }
		virtual int v_get_ionm(char *to, int io, int j);
		virtual void box_window() { upd_window(-1); }
		virtual void spec_debug() { log("itbox: sb=#%x, ni = %d", 
						m_bx1 ? m_bx1->node()->id() : 0, m_ni); }
		virtual const char * v_rgb() { return m_bx1 ? m_bx1->v_rgb() : "KKK%%%"; }
	protected:
		BXCMD_DECL(ItBoxGen) c_bx1, c_scl, c_rel;
		virtual void set_mdl();
		void upd_window(int flg);
		void upd2(int flg) { unset_model(), m_node->nio_change(), upd_window(flg); }
		int set_bx1_2(BoxGen * bx, int ni);
		int set_bx1(BoxGen * bx);
		int mflg() const { return (m_bx1->io_alias_perm()<<30) + 32*m_scl + (m_flg<<16) + m_bx1->n_in()-1; }
		BoxGen * m_bx1;
		char m_ni, m_scl, m_flg;
};

void ItBoxGen::set_mdl() {
	ItBoxModel *mo = m_mdlp.mk0<ItBoxModel>(0);
	if (m_bx1) m_bx1->mdl_cpto(&mo->box1m), mo->flg=mflg(); else mo->flg=0; }

// in i(1)...i(skip) x(0)...x(nx) i(skip+1)...i(ni1-1)

int ItBoxInst::ini(double **par) { log("b1m/ii=%p", m_b1m.rawmp());
	int n = m_n = m_b1m.nz() ? ivlim((int)lround(**par), 0, 1024) : 0; if (!n) return 0;
	m_ppbx = (BoxInst**) malloc(n*sizeof(BoxInst*));
	m_b1m.mk_boxv(m_ppbx, n);
	int sc = m_flg & 0x1e0; if (!sc) return n; else sc >>= 5;
	double v0 = par[1][0], v1 = (m_flg&IBF_REL) ? v0*par[2][0] : par[2][0];
	m_zvec = (double*)malloc(8*n);
	Scale01 scl; scl.set_all(v0, v1, sc==1 ? (int)lround(par[3][0]) : sc-5);
	if (n==1) return m_zvec[0] = scl.f(0.5), 1;
	double stp = 1.0/((double)(n-1)), x = stp;
	m_zvec[0] = v0; m_zvec[n-1] = v1;
	for (int i=1; i<n-1; i++, x+=stp) m_zvec[i] = scl.f(x);
	return n;
}

#define CALC1(J,Q) if ((r = m_ppbx[J]->calc(ifg1|of, in1, Q, n)) < 0) return r; else of = r&1
int ItBoxInst::calc(int inflg, double** inb, double** outb, int n) {
	int nr = (m_n>=0) ? m_n : ini(inb+1); if (!nr) return BOX_CP0;
	int r, flg = m_flg, of = inflg&1, sc = (m_flg>>5)&15, ni = (flg&31) - !!(sc),
	    xin = sc ? 3+(sc==1) : 1, ifg1 = (inflg>>xin) & ~1, ali = m_flg&IBF_ALI;
	double *in1[ni+2], **qq = 0;
	if (ni) memcpy(in1+1, inb+1+xin, sizeof(double*)*ni);
	if (sc) *(qq=in1+ni+1) = m_zvec;
	if (nr==1) return in1[0]=inb[0], m_ppbx[0] -> calc(ifg1|of, in1, outb, n);
	if (ali) {
		in1[0] = inb[0]; CALC1(0, outb); in1[0] = outb[0];
		for (int j=1; j<nr; j++) { if (qq) ++*qq; CALC1(j, outb); }
	} else {
		double tmp[n], *pt = tmp; 
		int j0 = nr & 1; 
		if(j0){ if (of) memcpy(tmp, *inb, 8*n); else tmp[0] = **inb;
			in1[0] = tmp; CALC1(0, outb); in1[0] = outb[0]; if(qq)++*qq; }
		else  { in1[0] = inb[0]; }
		for (int j=j0; j<nr; j+=2) { CALC1(j  , &pt) ; in1[0] = tmp;     if(qq)++*qq;
					     CALC1(j+1, outb); in1[0] = outb[0]; if(qq)++*qq; }
	}
	return of;
}

int ItBoxGen::save2(SvArg * sv) {
	BXSV2_HEAD; 
	if (m_bx1)   { CHKERR(f->sn("X$b", 3)); CHKERR(m_bx1->node()->sv_path(10)); }
	if (m_scl)   { CHKERR(f->pf("X$s%c\n", 48+m_scl)); }
	if (m_flg&1) { CHKERR(f->sn("X$r1\n", 5)); }
	return r;
}

int ItBoxGen::set_bx1_2(BoxGen * bx, int ni) { 
	int f = 6, r = set_boxp(&m_bx1,bx); if (r<0) return r;
	if (ni<2 && m_scl) m_scl = 0, f += 8; 
	return m_ni=ni, upd2(f), 0;
}

int ItBoxGen::set_bx1(BoxGen * bx) {
	if (!bx) return set_bx1_2(0, 1);
	int ni = bx->n_in();
	if (!ni || ni>29 || bx->n_out() != 1) return bx==m_bx1 ? (set_bx1_2(0, 1), BXE_FILTRM) : BXE_FILTNIO;
	return set_bx1_2(bx, ni);
}

int ItBoxGen::v_get_ionm(char *to, int io, int j) {
	if (io) return memcpy(to, "out", 3), 3; else j &= 31;
	int nx = m_scl ? (3+(m_scl==1)) : 1, k = j-nx;
	if (k>0) return m_bx1->node()->get_ionm(to, 0, k);
	if (!j) return to[0]='i', to[1]='n', 2;
	if (j==3 && (m_flg&1)) ++j; else --j;
	return memcpy(to, "(#b)(z0)(z1)(z%)(z*)" + 4*j, 4), 4;
}

void ItBoxGen::upd_window(int flg) {
	if (flg& 1) gui2.cre(w_oid(), 'i'); else if (wnfl()) gui2.setwin(w_oid(), 'i'); else return;
	if (flg& 2) gui2.own_title();
	if (flg& 4) gui2.ref_title('3', m_bx1 ? m_bx1->node() : 0, -1, "filter");
	if (flg& 8) gui2.wupd_c48('s', m_scl);
	if (flg&16) gui2.wupd_i1('r', m_flg&1);
}

#undef CH
#define CH(X) BXCMD_H(ItBoxGen, X)

CH(bx1){ ANode * nd = 0; BoxGen * bx = 0;
	 if (!(nd = cb->lookup(s+1))) return BXE_ARGLU;
	 return (bx = nd->box0()) ? p->set_bx1(bx) : BXE_ARGNBX; }

CH(scl){ p->m_scl = ivlim(s[1]-48, 0, 8&-(p->m_ni>1)); p->upd2(8); return 0;  }
CH(rel){ if ((p->m_flg^s[1])&1) p->m_flg ^= 1, p->upd2(16); return 0; }

BXCMD_DEF(ItBoxGen) { {8192+'\\',0}, {'b',c_bx1}, {'s',c_scl}, {'r',c_rel}, {0,0} };

// export

void itb_init() { ItBoxGen::cmd_init(); }
int setbox_it(ABoxNode * nd, BoxGen * _) { nd->m_box = new ItBoxGen(nd); return 1; }
