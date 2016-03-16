#include "util.h"
#include "combo.h"
#include "guistub.h"
#include "cmd.h"
#include "glob.h"
#include "util2.h"

class ItBoxInst : public BoxInst {
	public:
		ItBoxInst(BoxModel * b1m, int flg) : m_b1m(b1m), m_ppbx(0), m_flg(flg), m_n(-1) { 
								BoxModel::ref(b1m); }
		virtual ~ItBoxInst() { BoxModel::unref(m_b1m); if (m_n<1) return; 
				       for (int i=0; i<m_n; i++) delete(m_ppbx[i]);  free(m_ppbx); }
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int ini(double in1);
		BoxModel * m_b1m;
		BoxInst ** m_ppbx;
		int m_flg, m_n;
};

class ItBoxModel : public BoxModel {
        public:
                ItBoxModel(BoxModel * b1m, int f) : box1m(b1m), flg(f) { BoxModel::ref(b1m); }
                virtual ~ItBoxModel() { BoxModel::unref(box1m); }
                virtual BoxInst * mk_box() { return new ItBoxInst(box1m, flg); }
	protected:
                BoxModel * box1m;
		int flg;
};

class ItBoxGen : public BoxGen {
        public: 
                friend void itb_init();
                ItBoxGen(ABoxNode * nd) : BoxGen(nd), m_bx1(0), m_ni(2) {}
                virtual ~ItBoxGen() { set_boxp(&m_bx1, 0); }
                virtual void set_model() { m_model = m_bx1 ? new ItBoxModel(m_bx1->model(), 
								 32*m_bx1->io_alias_perm()+m_bx1->n_in()-1)
						           : new ItBoxModel(0, 0); }
                virtual int n_in() const { return m_ni; }
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
		BXCMD_DECL(ItBoxGen) c_bx1;
		void upd_window(int flg);
		int set_bx1_2(BoxGen * bx, int ni);
		int set_bx1(BoxGen * bx);

		BoxGen * m_bx1;
		int m_ni;
};

int ItBoxInst::ini(double in1) {
	int n = m_n = m_b1m ? ivlim((int)lround(in1), 0, 1024) : 0;  if (!n) return 0;
	m_ppbx = (BoxInst**) malloc(n*sizeof(BoxInst*));
	for (int i=0; i<n; i++) m_ppbx[i] = m_b1m->mk_box();
	return n;
}

#define CALC1(J,Q) if ((r = m_ppbx[J]->calc(ifg1|of, in1, Q, n)) < 0) return r; else of = r&1
int ItBoxInst::calc(int inflg, double** inb, double** outb, int n) {
	int nr = (m_n>=0) ? m_n : ini(inb[1][0]);
	if (!nr) return BOX_CP0;
	int r, ni = m_flg & 31, ali = m_flg & 32, ifg1 = (inflg>>1) & ~1, of = inflg & 1;
	double *in1[ni+1]; if (ni) memcpy(in1+1, inb+2, ni*sizeof(double*));
	if (nr==1) return in1[0]=inb[0], m_ppbx[0] -> calc(ifg1|of, in1, outb, n);
	if (ali) {
		in1[0] = inb[0]; 	 CALC1(0, outb); in1[0] = outb[0];
		for (int j=1; j<nr; j++) CALC1(j, outb);
	} else {
		double tmp[n], *pt = tmp; 
		int j0 = nr & 1; 
		if(j0){ if (of) memcpy(tmp, *inb, 8*n); else tmp[0] = **inb;
			in1[0] = tmp; CALC1(0, outb); in1[0] = outb[0]; }
		else  { in1[0] = inb[0]; }
		for (int j=j0; j<nr; j+=2) { CALC1(j  , &pt) ; in1[0] = tmp;
					     CALC1(j+1, outb); in1[0] = outb[0]; }
	}
	return of;
}

int ItBoxGen::save2(SvArg * sv) {
	BXSV2_HEAD; 
	if (m_bx1) { CHKERR(f->sn("X$b", 3)); CHKERR(m_bx1->node()->sv_path(10)); }
	return r;
}

int ItBoxGen::set_bx1_2(BoxGen * bx, int ni) { 
	int r = set_boxp(&m_bx1,bx);
	if (r>=0) unset_model(), m_ni=ni, m_node->nio_change(), upd_window(6);
	return r; }

int ItBoxGen::set_bx1(BoxGen * bx) {
	if (!bx) return set_bx1_2(0, 2);
	int ni = bx->n_in();
	if (!ni || ni>29 || bx->n_out() != 1) return bx==m_bx1 ? (set_bx1_2(0, 2), BXE_FILTRM) : BXE_FILTNIO;
	return set_bx1_2(bx, ni + 1);
}

int ItBoxGen::v_get_ionm(char *to, int io, int j) {
	if (io) return memcpy(to, "out", 3), 3; else j &= 31;
	switch(j) { case 0: return to[0]='i', to[1]='n', 2;
		    case 1: return to[0]='#', to[1]='b', 2;
		    default:return m_bx1 ? m_bx1->node()->get_ionm(to,io,j-1) : 0; }}

void ItBoxGen::upd_window(int flg) {
	if (flg&1) gui2.cre(w_oid(), 'i'); else if (wnfl()) gui2.setwin(w_oid(), 'i'); else return;
	if (flg&2) gui2.own_title();
	if (flg&4) gui2.ref_title('3', m_bx1 ? m_bx1->node() : 0, -1, "filter");
}

#define CH(X) BXCMD_H(ItBoxGen, X)

CH(bx1){ ANode * nd = 0; BoxGen * bx = 0;
	 if (!(nd = cb->lookup(s+1))) return BXE_ARGLU;
	 return (bx = nd->box0()) ? p->set_bx1(bx) : BXE_ARGNBX; }

BXCMD_DEF(ItBoxGen) { {8192+'\\', 0}, {'b', c_bx1}, {0, 0} };

// export

void itb_init() { ItBoxGen::cmd_init(); }
int setbox_it(ABoxNode * nd, BoxGen * _) { nd->m_box = new ItBoxGen(nd); return 1; }
