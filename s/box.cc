#include "util.h"
#include "glob.h"
#include "cmd.h"
#include "guistub.h"
#include "pt.h"

BoxGen * box_bookmark[8];
void box_mxc_notify(BoxGen *p, int ky, int flg) { p->mxc_notify(ky, flg); }

void BoxInst::rmcon(int flg, double **pp, int n) { BVFOR_JM(flg) {
	double *p=pp[j], x=*p; for (int t=1; t<n; t++) p[t] = x; }}

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
void BoxModel::unref(BoxModel *p) {
	if (p && !--p->refcount) {
		if (debug_flags&DFLG_MODELDEL) log("model_del: %p",p); delete(p); }}

     BoxGen::~BoxGen()    { unset_model(); }
void BoxGen::spec_debug() { log("no class-specific debug info"); }
void BoxGen::box_window() { log("BUG: undefined box_window() p:%p, #%x", m_node, m_node?m_node->id():-1); }
void BoxGen::doc_window(int id4) { m_node->winflg_or(1<<id4); gui2.cre(w_oid(id4), 'D', "");
	gui2.bn_dsc(m_node); gui2.own_title(); }

int BoxGen::set_boxp(BoxGen ** pp, BoxGen * to) {
	int ec = 0; if (to && (ec=Node::conn(m_node, to->node())) < 0) return ec;
	if (*pp) ec = Node::disconn(m_node, (*pp)->node());
	*pp = to; return ec;
}

PrimBoxGen::PrimBoxGen(BoxModel * mdl, int ni, int no, const char * cl) : m_ni(ni), m_no(no) { 
	m_model = mdl; m_cl[0] = '_'; cl += (*cl=='_');
	int i; for (i=0; i<4 && cl[i] && cl[i]!='.'; i++) m_cl[i+1] = cl[i]; m_cl[i+1] = 0;  }

PrimBoxGen::PrimBoxGen(qmb_arg_t qa, int k, int ni, int no, const char * cl) : m_ni(ni), m_no(no) { 
	m_model = (*qa)(m_mdl_spc, k); m_cl[0] = '_'; cl += (*cl=='_');
	int i; for (i=0; i<4 && cl[i] && cl[i]!='.'; i++) m_cl[i+1] = cl[i]; m_cl[i+1] = 0;  }

void PrimBoxGen::box_window() { doc_window(11); if (!memcmp(m_cl+1, "abou", 4)) pt_show_lic(); }
int  PrimBoxGen::n_in() const { return m_ni; }
int  PrimBoxGen::n_out() const { return m_no&31; }
bool PrimBoxGen::io_alias_perm() const { return !!(m_no&32); }
