#include "util.h"
#include "glob.h"
#include "cmd.h"
#include "guistub.h"
#include "pt.h"

class BoxModel_SL : public BoxModel {
//        public: BoxModel_SL(BoxInst::sc_t f, int k) : BoxModel(0,2){(m_bx=(BoxInst*)nf_alloc(16))->fdk(f,0,k);}
        public: BoxModel_SL(BoxInst::sc_t f, int k) : BoxModel(0,2){
			(m_bx=(BoxInst*)nf_alloc(16))->fdk(f,0,k); log("stateless-box %p %x",m_bx,k);}
                virtual BoxInst * place_box(void *p) { return m_bx; }
        protected:
		BoxInst * m_bx;
};

class BoxModel_B0 : public BoxModel {
	public: BoxModel_B0(BoxInst::sc_t f, int k, int siz) : BoxModel(siz, 2) { m_proto.fdk(f, 1, k); }
		virtual BoxInst * place_box(void *p) { BoxInst *q = (BoxInst*)p;
						       memcpy(q, &m_proto, sizeof(BoxInst)); return q; }
	protected:
		BoxInst m_proto;
};

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

void BoxInst::del(BoxInst *p) { if(p&&p->m_dt) { if (p->m_dt>1u) dtor2(p); free(p); }} 

void BoxInst::dtor2(BoxInst *p) { return (p->m_dt==2) ? free(((BoxInst_B1*)p)->m_p0)
						      :    (*((BoxInst_BU*)p)->m_dtf)(p); }

double *zero30[30], *junk30[30];
BoxGen * box_bookmark[10];
int box_mxc_notify(BoxGen *p, int ky, int flg) { return p->mxc_notify(ky, flg); }

void BoxInst::rmcon(int flg, double **pp, int n) { BVFOR_JM(flg) {
	double *p=pp[j], x=*p; for (int t=1; t<n; t++) p[t] = x; }}

BX_SCALC(BoxInst::sc_bug)  { bug("sc_bug(%p,%d,%p,%p,%d)", abxi, inflg, inb, outb, n); }

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

void BoxModel::debug0() const { log("%p,%d,%d", this, refcount, m_size); }
void BoxModel::del(BoxModel *p) { IFDBGX(MODELDEL) log("model_del: %p",p); p->~BoxModel(); free(p); }
void ModelPtr::debug() const { if (m_p) m_p->debug0(); else log("no model"); }


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

int BoxInst::sc_ctor_2(BoxInst *p, BoxInst::sc_t f0, int siz) { p->m_psc=f0; p->m_arg=siz; return BXE_CTOR; }

static BoxModel * place_b_mdl(char *q, qmb_arg_t qa, int k) {
	BoxInst b; b.m_arg=k; int r = qa(&b, 0, zero30, junk30, 1);
	return (r==BXE_CTOR) ? (BoxModel*) (new (q) BoxModel_B0(b.m_psc, k, b.m_arg))
			     : (BoxModel*) (new (q) BoxModel_SL(qa,	    k)	    ); }

PrimBoxGen::PrimBoxGen(qmb_arg_t qa, int k, int ni, int no, const char * cl) 
	: BoxGen(place_b_mdl(m_mdl_spc, qa, k)), m_ni(ni), m_no(no) { 
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
void dat_init() { for (int i=0; i<32; i++) pstat_con[i] = stat_con+i;
		  for (int i=0; i<30; i++) zero30[i] = zeroblkD;
		  for (int i=0; i<30; i++) junk30[i] = junkbuf;  }
