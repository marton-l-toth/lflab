#include "util.h"
#include "glob.h"
#include "cmd.h"
#include "guistub.h"
#include "pt.h"

class HelpBoxGen : public BoxGen {
	public:
		HelpBoxGen(ABoxNode * nd) : BoxGen(nd) {}
		virtual ~HelpBoxGen() {}
		virtual int n_in() const { return 0; }
		virtual int n_out() const { return 0; }
                virtual bool io_alias_perm() const { bug("hlp: io_alias_perm() called"); return 0; }
		virtual const char * cl_name() { return "help"; }
		virtual void box_window() { return doc_window(11); }
	protected:
                virtual void set_mdl() { bug("hlp: set_model() called"); }
};

BoxGen * box_bookmark[10];
int box_mxc_notify(BoxGen *p, int ky, int flg) { return p->mxc_notify(ky, flg); }

void BoxInst::rmcon(int flg, double **pp, int n) { BVFOR_JM(flg) {
	double *p=pp[j], x=*p; for (int t=1; t<n; t++) p[t] = x; }}

BX_SCALC(BoxInst::sc_bug)  { bug("sc_bug(%p,%d,%p,%p,%d)", abxi, inflg, inb, outb, n); }
BX_SCALC(BoxInst::sc_zero) { return **outb=0.0, 0; }
BX_SCALC(BoxInst::sc_cp0)  { return BOX_CP0; }

int BoxInst::calc_nzo2(int ocfg, double *o0, double *o1, int inflg, double **inb, int n) {
	if (ocfg==0x7c01) {
		int r = calc(inflg, inb, &o0, n); if (r<0 || (r&=1)) return r;
		double x = *o0; if (fabs(x)<1e-270) return 0;
		for (int i=1; i<n; i++) o0[i]=x;    return 1;
	}
	int no = ocfg&31, j0 = (ocfg>>5)&31, j1 = (ocfg>>10)&31, om0 = 1<<j0, om1 = (1<<j1) & 0x3fffffff;
	double x, *oa[no]; 
	for (int i=0; i<no; i++) oa[i] = junkbuf;   oa[j0] = o0; if (j1<no) oa[j1] = o1;
	int r = calc(inflg, inb, oa, n); if (r<0) return r; else r = 2*!!(r&om1) + !!(r&om0);
	if (!(r&1) && fabs(x=*o0)>=1e-270) { r|=1; for (int i=1; i<n; i++) o0[i]=x; }
	if (!om1) return r; if (j0==j1) return r ? (memcpy(o1, o0, 8*n), 3) : 0;
	if (!(r&2) && fabs(x=*o1)>=1e-270) { r|=2; for (int i=1; i<n; i++) o1[i]=x; }
	return r;
}

void BoxModel::debug0() { log_n("%p,%d", this, this?refcount:0); }
void BoxModel::del(BoxModel *p) { IFDBGX(MODELDEL) log("model_del: %p",p); p->~BoxModel(); free(p); }

void BoxGen::spec_debug() { log("no class-specific debug info"); }
void BoxGen::box_window() { log("BUG: undefined box_window() p:%p, #%x", m_node, m_node?m_node->id():-1); }
void BoxGen::doc_window(int id4) { 
	m_node->winflg_or(1<<id4); gui2.cre(w_oid(id4), 'D', ""); gui2.bn_dsc(m_node); gui2.own_title(); 
	gui2.wupd_c0('+','M'); gui2.hex4(id4==13?0x3c0:(m_mdlp.nz()?0x33f:0x21)); }

int BoxGen::set_boxp(BoxGen ** pp, BoxGen * to) {
	int ec = 0; if (to && (ec=Node::conn(m_node, to->node())) < 0) return ec;
	if (*pp) ec = Node::disconn(m_node, (*pp)->node());
	*pp = to; return ec;
}

PrimBoxGen::PrimBoxGen(BoxModel * mdl, int ni, int no, const char * cl) : BoxGen(mdl), m_ni(ni), m_no(no) { 
	m_cl[0] = '_'; cl += (*cl=='_');
	int i; for (i=0; i<4 && cl[i] && cl[i]!='.'; i++) m_cl[i+1] = cl[i]; m_cl[i+1] = 0;  }

PrimBoxGen::PrimBoxGen(qmb_arg_t qa, int k, int ni, int no, const char * cl) 
	: BoxGen((*qa)(m_mdl_spc, k)), m_ni(ni), m_no(no) { 
		m_cl[0] = '_'; cl += (*cl=='_'); int i;
		for (i=0; i<4 && cl[i] && cl[i]!='.'; i++) m_cl[i+1] = cl[i]; m_cl[i+1] = 0; }

void PrimBoxGen::box_window() { doc_window(11); if (this==box_bookmark[2]) pt_show_lic(); }
int  PrimBoxGen::n_in() const { return m_ni; }
int  PrimBoxGen::n_out() const { return m_no&31; }
bool PrimBoxGen::io_alias_perm() const { return !!(m_no&32); }

double stat_con[32] = {-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
                       .1,.2,.3,.4,.5,.6,.7,.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6};
double *pstat_con[32];

ConStore * ConStore::cp(char * to) const {
        ConStore * r = (ConStore*)to; memcpy(r, this, sizeof(ConStore));
        if (m_n) memcpy(r->m_p = (double*)(to+((sizeof(ConStore)+7)&~7)), m_p, 8*m_n);
        return r; }

void ConStore::rm(int j) { if ((j-=32) >=0 ) m_p[j] = (double)m_fh, m_fh = j; }

int ConStore::add(double x) {
        int j = ((int)x + 5)&15;       if (stat_con[j]==x) return j;
        j = (((int)(10.0*x)-1)&15)+16; if (stat_con[j]==x) return j;
        if (m_fh>=0) return j=m_fh, m_fh=(int)m_p[j], m_p[j]=x, j+32;
        if (m_a==m_n) m_p = (double*)(m_p ? realloc(m_p, 8*(m_a*=2)) : malloc(8*(m_a=8)));
        m_p[j=m_n++] = x; return j+32;
}

int setbox_hlp(ABoxNode * nd, BoxGen * _) { nd->m_box = new HelpBoxGen(nd); return 0; }
void reg_bn(ANode * nd, int i) { BoxGen **qq = box_bookmark+i; 
				 if (*qq) log("BUG: duplicate box_bkm: %s", nd->s_name()); *qq = nd->box0(); }
void dat_init() { for (int i=0; i<32; i++) pstat_con[i] = stat_con+i; }
