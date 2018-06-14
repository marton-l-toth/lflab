#include "util.h"
#include "combo.h"
#include "guistub.h"
#include "cmd.h"
#include "glob.h"
#include "util.h"
#include "util2.h"
#include "pt.h"

class GRBX_Node {
        public: 
                GRBX_Node(BoxGen* bx, int k=0);
                ~GRBX_Node() { delete[]iarg; }
                BoxGen* box;
                int ni, no;
                int *iarg, *oarg;
};

class GraphBoxGen : public BoxGen {
        public:
                friend void graph_init();
                GraphBoxGen(ABoxNode * nd);
                ~GraphBoxGen();
                virtual void set_mdl();
                virtual int n_in() const { return m_n_in - m_n_fb; }
                virtual int n_out() const { return m_n_out - 2*m_n_fb; }
                virtual const char * cl_name() { return "graph"; }
                //virtual int cmd(CmdBuf* cb);
                virtual bool io_alias_perm() const { return false; }
                virtual int save2(SvArg * sv);
                virtual void notify_nio(BoxGen * bx);
                virtual void box_window();
		virtual int ifflg() const { return BIF_GC; }
                int n_subbox() { return m_nodes.n() - 1; }
                BoxGen * subbox(int i) { return m_nodes[i]->box; }
                int set_n_in(int n);
                int set_n_out(int n);
                int n_box() { return m_nodes.n() - 1; }
                int n_edg();
		void decr_no(int i, int no);
                int add_box(int pos, BoxGen* bx);
                int del_box(int pos);
		int rplc_box(int pos, BoxGen * bx);
                int cp_box(int pos);
                bool xchg_box(int pos /* , pos+1 */ );
                // tgtbox: 0...n_box()-1 -- box n_box():output 
                // srcbox: 0...n_box()-1 -- box -1:input -2: const
                int conn(int tgtbox, int tgtix, int srcbox, int srcix, double v=0.0);
                void write_dot(FILE * fp);
                GRBX_Node** get_nodes() { return m_nodes.p(); }
                double get_con(int i) { return m_cs.v(i); }
		void w_rpflg() { gui2.setwin(w_oid(), 'g'); gui2.wupd_i1('r', m_op); }
                void shuffle();
                void to_dot(int seq, int ne);
        protected:
                BXCMD_DECL(GraphBoxGen) c_ni, c_no, c_nfb, c_conn, c_ibx, c_dbx, c_shfl, c_dbg,
                                        c_gnd, c_g2, c_ck, c_uic, c_rpl;
                typedef void (rcf_ft)(int*,int,int), (*rcf_t)(int*,int,int);
		static rcf_ft rcf_ins, rcf_cut, rcf_xcg, rcf_nin, rcf_fb;

                void dump_model();
                void reconn(int pos, int i0, rcf_t fun);
                void reorder(int * ix);
                bool top_ord(bool rnd = false);
                void debug();
                void draw_node(int i);
                void click(const char * s);
                void sel3(int t, int i, int j);
                int sel_ix(int t) { return t=m_sel[t]>>5, t>1021?-1:t; }
                int gui2_cmd(const char * s);
                int gnd_cmd(ANode * nd);
		void uinm1(int j);
		int uicfg(int ch);
                void upd_con(int i, int j);
                int inm_2(char*q,int i){ int k=i-n_in(); return k<0 ? get_ionm(q,0,i) : lb3(q,2,k); }
                int onm_2(char*q,int i){ int k=i-n_out();return k<0 ? get_ionm(q,1,i) : lb3(q,k&1,k>>1); }
		int lb3(char*q,int t,int j){ return *(q++)="!><"[t], j<10 ? (*q=48+j,2):(*q=49,q[1]=38+j,3);}
                LWArr<GRBX_Node*> m_nodes;
		ConStore m_cs;
                int m_n_in, m_n_out, m_n_fb;
                short m_dseq;
                short m_sel[3];
                double m_guitmp;
		char m_op; // 0:ins 1:repl
};

/////////////////////////////////////////////////////////////////////////////////////////

GRBX_Node::GRBX_Node(BoxGen* bx, int k) : box(bx), iarg(0) {
	if (bx) ni = bx->n_in(), no = bx->n_out();
	else ni = k, no = 0;
	int n = ni + no; if (!n) return;
	iarg = new int[ni+no];
	for (int i=0; i<ni; i++) iarg[i]=0x5fffc;
	oarg = iarg + ni;
}

GraphBoxGen::GraphBoxGen(ABoxNode * nd) : BoxGen(nd), m_n_in(0), m_n_out(1), m_n_fb(0),
	m_dseq(0), m_guitmp(0.0), m_op(0) { m_nodes.add(new GRBX_Node(0,1)); }

GraphBoxGen::~GraphBoxGen() {
	for (int i=0; i<n_box(); i++) 
		delete m_nodes[i];
	m_nodes.resize(0); // default destructor would not work
}

void GraphBoxGen::reconn(int pos, int i0, rcf_t fun) {
	for (int i=i0; i<m_nodes.n(); i++) {
		int ni = m_nodes[i]->ni;
		for (int j=0; j<ni; j++)
			(*fun)(m_nodes[i]->iarg+j,pos,i);
	}}

void GraphBoxGen::rcf_ins(int *p, int pos, int ix) {
	int k = *p&0xffff;
	if (k>pos && k<0xfff0) ++(*p);
}
void GraphBoxGen::rcf_cut(int *p, int pos, int ix) { 
	int k = *p&0xffff;
	if (k==pos) *p = 0x5fffc;
	else if (k>pos && k<0xfff0) --(*p);
}
void GraphBoxGen::rcf_xcg(int *p, int pos, int ix) {
	int k = *p&0xffff;
	if (k==pos) ix==pos ? (*p=0x5fffc) : ++(*p);
	else if (k==pos+1) --(*p);
}
void GraphBoxGen::rcf_nin(int *p, int pos, int ix) {
	if ((*p&0xffff)==0xfffe && (*p>>16) >= pos) *p=0x5fffc;
}
void GraphBoxGen::rcf_fb(int *p, int pos, int ix) { // pos: old<<16 + new
	if ((*p&0xffff)==0xfffe && *p >= (pos&0x7fff0000)) *p += (pos<<16) - (pos&0x7fff0000); }

int GraphBoxGen::add_box(int pos, BoxGen* bx) {
	if (n_box()==0xfff0) return BXE_GFULL; if (!bx->n_out()) return BXE_ZOUT;
	BoxGen * p2 = 0;
	int ec = set_boxp(&p2, bx); if (ec<0) return ec;
	unset_model();
	if (pos<0) pos=0;
	else if (pos>n_box()) pos = n_box();
	m_nodes.ins(pos, new GRBX_Node(bx));
	reconn(pos, pos, rcf_ins);
	return BCR_WIN;
}

int GraphBoxGen::del_box(int pos) {
	if (pos<0 || pos>=n_box()) return BXE_IDX;
	unset_model();
	BoxGen * bx = m_nodes[pos] -> box;
	delete m_nodes[pos];
	m_nodes.cut(pos);
	reconn(pos, pos, rcf_cut);
	int ec = set_boxp(&bx, 0); return ec<0 ? ec : BCR_WIN;
}

int GraphBoxGen::rplc_box(int pos, BoxGen * bx) {
	if (pos<0 || pos>=n_box()) return BXE_IDX; else unset_model();
	GRBX_Node * gnd = m_nodes[pos];  if (gnd->box==bx) return EEE_NOEFF;
	int ec = set_boxp(&gnd->box, bx); if (ec<0) return ec;
	int ni2 = bx->n_in(), no2 = bx->n_out(), *iarg2 = new int[ni2+no2];
	memcpy(iarg2, gnd->iarg, 4*min_i(gnd->ni, ni2));
	for (int i=gnd->ni; i<ni2; i++) iarg2[i] = 0x5fffc;
	if (no2 < gnd->no) decr_no(pos, no2);
	gnd->ni = ni2; gnd->no = no2; 
	delete[](gnd->iarg); gnd->oarg = (gnd->iarg = iarg2) + ni2;
	return BCR_WIN;
}

int GraphBoxGen::cp_box(int pos) {
	if (pos<0 || pos>=n_box()) return BXE_IDX;
	if (n_box()==0xfff0) return BXE_GFULL;
	unset_model();
	BoxGen * bx = 0;
	int ec = set_boxp(&bx, m_nodes[pos]->box); if (ec<0) return ec;
	m_nodes.ins(pos, new GRBX_Node(bx));
	GRBX_Node * p = m_nodes[pos], *p0 = m_nodes[pos+1];
	for (int i=0; i<p->ni; i++) {
		int k = p0->iarg[i];
		p->iarg[i] = ((k&0xffff)!=0xfffc || k<0x20fffc) ? k
			: 0xfffc + (m_cs.add(m_cs.v(k>>16))<<16);
	}
	reconn(pos, pos, rcf_ins);
	return 0;
}

int GraphBoxGen::conn(int tgtbox, int tgtix, int srcbox, int srcix, double v) {
	int nb = n_box();
	int r;
	if (tgtbox<0 || tgtbox>nb) return BXE_IDX;
	GRBX_Node * gnd = m_nodes[tgtbox];
	if (tgtix<0 || tgtix>=gnd->ni) return BXE_IDX;
	int *ap = gnd->iarg + tgtix;
	if (srcbox<-2||srcbox>=nb) return BXE_IDX;
	if (srcix<0) return BXE_IDX;
	switch(srcbox) {
		case -1:
			if (srcix>=m_n_in) return BXE_IDX;
			r = (srcix<<16) + 0xfffe;
			break;
		case -2:
			r = 0xfffc + (m_cs.add(v)<<16); break;
			break;
		default:
			if (srcix>=m_nodes[srcbox]->no) return BXE_IDX;
			r = (srcix<<16) + srcbox;
			break;
	}
	int sav = *ap;
	*ap = r;
	if (srcbox < tgtbox || top_ord()) {
		if ((sav & 0xffff)==0xfffc) m_cs.rm(sav>>16);
		unset_model(); return BCR_WIN;
	} else { *ap = sav; return BXE_ACYC; }}

int GraphBoxGen::n_edg() {
	int r = 0, nb = n_box();
	for (int i=0; i<=nb; i++) {
		GRBX_Node * nd = m_nodes[i];
		for (int j=0; j<nd->ni; j++) 
			if ((nd->iarg[j]&0xffff)!=0xfffc) ++r;
	}
	return r;
}


void GraphBoxGen::to_dot(int seq, int ne) {
	int nb = n_box();
	Dot * to = Dot::sg();
	if (!to->active()) return log("gr/to_dot: error: dot not running!");
	to->ghead(16*m_node->id()+11, 'g', 'g', -1, seq, nb+2, ne);
	if (m_n_in) to->gnode(0, 0, m_n_in, -2);
	for (int i=0; i<nb; i++) {
		GRBX_Node * nd = m_nodes[i];
		to->gnode(i+1, nd->ni, nd->no, 
			       (text_wid16(nd->box->node()->s_name())) >> 5);
	}
	to->gnode(nb+1, m_n_out, 0, -1);
	for (int i=0; i<=nb; i++) {
		GRBX_Node * nd = m_nodes[i];
		for (int j=0; j<nd->ni; j++) {
			int ag = nd->iarg[j], n0 = ag&0xffff, o0 = ag>>16;
			if (n0==0xfffc) continue;
			if (n0==0xfffe) n0=0; else ++n0;
			to->gedge(n0, o0, i+1, j);
		}
	}
	to->gend();
}

void GraphBoxGen::draw_node(int i) {
	GRBX_Node * nd; int ni, no; char buf[8];
	ABoxNode * bnd; BoxGen * bx;
	if (i<0) nd = 0, ni = 0, no = m_n_in, bx = 0, bnd = 0;
	else nd = m_nodes[i], ni = nd->ni, no = nd->no, bx = nd->box, bnd = bx ? bx->node() : 0;
	gui2.setwin(w_oid(), 'g'); gui2.gn_start(i, bx, ni, no);
	for (int j=0; j<ni; j++) {
		int ag = nd->iarg[j], ag1 = ag&0xffff, ag2=ag>>16;
		if (ag1==0xfffc) gui2.hdbl(m_cs.v(ag2));
		else gui2.b32n( 0xc0000 + 32*(ag1&0x3ff) + ag2, 4 );
		if (bnd) gui2.ionm(bnd, 0, j); else gui2.sn(buf, onm_2(buf, j));
		gui2.c1(36);
	}
	if (bnd) for (int j=0; j<no; j++) gui2.ionm(bnd, 1, j), gui2.c1(36);
	else	 for (int j=0; j<no; j++) gui2.sn(buf, inm_2(buf, j)), gui2.c1(36);
}

void GraphBoxGen::box_window() {
	int nb = n_box(), ne = n_edg();
	int sq = (++m_dseq) & 255;
	gui2.cre(w_oid(), 'g'); gui2.own_title();
	gui2.wupd_0('g', "i"); gui2.hexn(sq, 2); gui2.c1('_');
	gui2.hexn(nb+2, 3); gui2.c1('_'); gui2.hexn(ne, 4);
	for (int i=!m_n_in-1; i<=nb; i++) draw_node(i);
	gui2.wupd_0('g', "z");
	gui2.wupd_i2('i', m_n_in); gui2.wupd_i2('o', m_n_out); gui2.wupd_i2('f', m_n_fb);
	w_rpflg();
	to_dot(sq, ne);
	m_sel[0] = m_sel[1] = m_sel[2] = -1;
}

int GraphBoxGen::set_n_in(int n) {
	if (n==m_n_in) return 0;
	unset_model();
	if (m_n_fb) reconn( ((m_n_in-m_n_fb)<<16) + n - m_n_fb,0,rcf_fb);
	if (n<m_n_in) reconn(n,0,rcf_nin);
	m_n_in = n; return BCR_NIO|BCR_WIN;
}

int GraphBoxGen::set_n_out(int n) {
	if (n==m_n_out) return 0;
	unset_model();
	GRBX_Node * obn = m_nodes[n_box()];
	int* p = n ? new int[n] : 0;
	int n_cp = n<m_n_out ? n : m_n_out;
	if (n_cp) memcpy(p, obn->iarg, n_cp*sizeof(int));
	for (int i=m_n_out; i<n; i++) p[i] = 0x5fffc;
	if (obn->iarg) delete[](obn->iarg);
	obn->iarg = p;
	m_n_out = obn->ni = n;
	return BCR_NIO|BCR_WIN;
}

void GraphBoxGen::reorder(int * ix) { // int[n_box()+1]  
	int nb = n_box();
	int ix_1[nb+1];
	for (int i=0; i<=nb; i++) ix_1[ix[i]] = i;
	GRBX_Node** gnp = m_nodes.forget();
	m_nodes.resize(nb+1);
	for (int i=0; i<=nb; i++) {
		GRBX_Node * gn = m_nodes[i] = gnp[ix[i]];
		for (int j=0; j<gn->ni; j++) {
			int k = gn->iarg[j];
			if ((k&0xffff) < 0xfffc) 
				gn->iarg[j] = (k&~0xffff) | ix_1[k&0xffff];
		}
	}
	delete[] (gnp);
}

void GraphBoxGen::shuffle() {
	int nb = n_box(), ix[nb+1];
	for (int i=0; i<=nb; i++) ix[i] = i;
	::shuffle(ix, nb); reorder(ix);
	if (!top_ord(true)) log("BUG: gr/shuffle: top_ord() failed");
}

bool GraphBoxGen::top_ord(bool rnd) {
	int nb = n_box();
	int out_sum[nb];
	int zo_0=0, zo_n=0, zo_list[nb];
	LWArr<int> ixix;

	for (int i=0; i<nb; i++) out_sum[i] = 0;
	for (int i=0; i<nb; i++) {
		GRBX_Node * gn = m_nodes[i];
		for (int j=0; j<gn->ni; j++) {
			int k = gn->iarg[j] & 0xffff;
			if (k<0xfffc) ++out_sum[k];
		}
	}
	for (int i=nb-1; i>=0; i--) 
		if (!out_sum[i]) zo_list[zo_n++] = i;
	while (zo_0 < zo_n) {
		GRBX_Node * gn = m_nodes[zo_list[zo_0++]];
		if (rnd) {
			ixix.resize(gn->ni);
			for (int j = 0; j < gn->ni; j++) ixix[j] = j;
			::shuffle(ixix.p(), gn->ni);
		}
		for (int j0=0; j0<gn->ni; j0++) {
			int j = rnd ? ixix[j0] : j0;
			int k = gn->iarg[j] & 0xffff;
			if (k<0xfffc && !--out_sum[k])
				zo_list[zo_n++] = k;
		}
	}
	if (zo_n < nb) return false;
	int ix[nb+1];
	ix[nb] = nb;
	for (int i=0; i<nb; i++) ix[i] = zo_list[nb-i-1];
	reorder(ix);
	return true;
}

void GraphBoxGen::upd_con(int i, int j) { int k = m_nodes[i]->iarg[j]; if ((k&0xffff)==0xfffc) 
					  	gui2.setwin(w_oid(), 'g'), gui2.grc(i+1, j, m_cs.v(k>>16)); }
void GraphBoxGen::uinm1(int j) {
	BoxGen * bx = m_nodes[j]->box;  if (!bx) return;
	int ni = m_nodes[j]->ni, *p = m_nodes[j]->iarg;
	char buf[8];
	for (int k,i=0; i<ni; i++) if (((k=p[i])&0xffff)==0xfffe)
		k>>=16, buf[bx->get_ionm(buf, 0, i)]=0, m_node->set_ui_lbl(0, k, buf);
}

int GraphBoxGen::uicfg(int ch) {
	int j; switch(ch) {
		case 'i': return (j=sel_ix(1))>=0 ? (uinm1(j), BCR_WIN) : BXE_GBXSEL;
		case 'r': if ((j=sel_ix(1)) < 0) return BXE_GBXSEL;
			  m_node->set_ui_rgb(m_nodes[j]->box->v_rgb());
			  if (m_node->winflg(2048)) gui2.setwin(w_oid(), 'g'), gui2.own_title(); 
			  return 0;
		case 'I': if (!(j=n_box())) return EEE_NOEFF; for (; j>=0; j--) uinm1(j);  return BCR_WIN;
		default: return BXE_CENUM;
	}}

void GraphBoxGen::sel3(int t, int i, int j) {
	j &= ((t-1)|(1-t));
	short v = 32*(i&1023) + j;
	if (m_sel[t]==v) return;
	m_sel[t] = v;
	gui2.setwin(w_oid(), 'g'); 
	gui2.wupd_0('g', "+0\0+1\0+2\0"+3*t); gui2.b32n(32*((i+1)&1023)+j, 3);
	if (t) return;
	int k = m_nodes[i]->iarg[j];
	if ((k&0xffff)==0xfffc) gui2.wupd_d('x', m_guitmp = m_cs.v(k>>16));
	else gui2.wupd_s('x', ""), m_guitmp = 0.0;
}

void GraphBoxGen::set_mdl() {
	LWArr<int> free_tbuf, cparg;
	int nb = n_box(), n_tmp = 1;
	for (int i=0; i<=nb; i++) for (int j=0,n=m_nodes[i]->no; j<n; j++) m_nodes[i]->oarg[j] = 0xffff;
	GRBX_Node * outnd = m_nodes[nb];
	for (int i=0; i<outnd->ni; i++) { // prepare copy boxes 
		int *q, ag = outnd->iarg[i], cur = (i<<16) + 0xfffd;
		if ((ag&0xffff)>0xfffb) cparg.add(ag, cur), free_tbuf.add(cur);
		else if (*(q = m_nodes[ag&0xffff]->oarg+(ag>>16))==0xffff) *q = cur;
		else cparg.add(*q, cur), free_tbuf.add(cur);
	}
	int n_iot = cparg.n(), n_cp = n_iot>>1;
	for (int i=nb-1; i>=0; i--) { // set oargs / alloc tmps
		GRBX_Node * nd = m_nodes[i]; 
		int ni = nd->ni, no = nd->no, perm = nd->box->io_alias_perm();
		n_iot += ni + no;
		if (perm) for(int c,j=0;j<nd->no;j++) if(((c=nd->oarg[j])&0xfffd)==0xfffd) free_tbuf.add(c);
		for (int k, t, ag, ty, *oap, j=0; j<nd->ni; j++) 
			if((ty=(ag=nd->iarg[j])&0xffff)<0xfffc && *(oap=m_nodes[ty]->oarg+(ag>>16))==0xffff)
				*oap = (k = free_tbuf.n()) ? (t = free_tbuf[k-1], free_tbuf.resize(k-1), t)
							   : ((n_tmp++)<<16) + 0xffff;
		if(!perm) for(int c,j=0;j<nd->no;j++) if(((c=nd->oarg[j])&0xfffd)==0xfffd) free_tbuf.add(c); 
	}
	int nb2 = nb+n_cp, nc = m_cs.n(), n_cb = (nc+63)>>5, s1 = (nc+n_cb+nb2)*sizeof(void*) + 2*(nb2+n_iot);
	ComboBoxModel * mdl = m_mdlp.mkc1<ComboBoxModel>(&m_cs, s1, nb2);
	mdl->n_t = n_tmp; mdl->niof = 65536*n_in() + 256*n_out() + m_n_fb;
	for (int i=0; i<nb; i++) {
		GRBX_Node * nd = m_nodes[i];
		BoxModel::ref(mdl->boxm[i] = nd->box->rawmp());
		mdl->dsc_nio(nd->ni, nd->no);
		for (int ag,ty,j=0; j<nd->ni; j++)
			mdl->dsc_arg( (ty=((ag=nd->iarg[j])&0xffff))>0xfffb ? ag :
					m_nodes[ty]->oarg[ag>>16] );
		for (int j=0; j<nd->no; j++) mdl->dsc_arg(nd->oarg[j]);
	}
	BoxModel *cpm = box_bookmark[3]->rawmp(), **mpp = mdl->boxm + nb;
	for (int i=0; i<n_cp; i++) BoxModel::ref(mpp[i]=cpm), mdl->dsc_nio(1,1),
				   mdl->dsc_arg(cparg[2*i]), mdl->dsc_arg(cparg[2*i+1]);
	int d0=2*(nb+n_cp+n_iot), d1=mdl->eob-mdl->iolist; if (d0!=d1) bug("gr/model: %d!=%d", d0, d1);
}

void GraphBoxGen::debug() {
	int nb = n_box();
	char buf[1025];
	log("-----GR_debug: %d boxes:",nb);
	ComboBoxModel * mp = static_cast<ComboBoxModel*> (rawmp());
	for (int i=0; i<=nb; i++) {
		char * p = buf;
		GRBX_Node * gn = m_nodes[i];
		const char * name = gn->box ? gn->box->cl_name() : "out";
		int oid = gn->box ? gn->box->node()->id() : 0;
		for (int j=0; j<gn->ni; j++) {
			int x = gn->iarg[j];
			switch(x&0xffff) {
				case 0xfffc: p += sprintf(p, " c%d", x>>16); break;
				case 0xfffe: p += sprintf(p, " i%d", x>>16); break;
				default: p += sprintf(p, " %d/%d", x&0xffff, x>>16); break;
			}

		}
		log("%02d.%s(%02d) > in:%s out:%s",i,name,oid,buf,p);
	}
	mp -> dump();
}

int GraphBoxGen::save2(SvArg * sv) {
	BXSV2_HEAD;
	CHKERR(xprintf(f,"X$<%d\n", n_in()));
	CHKERR(xprintf(f,"X$>%d\n", n_out()));
	if (m_n_fb) { CHKERR(xprintf(f,"X$@%d\n", m_n_fb)); }
	for (int i=0; i<n_box(); i++) {
		CHKERR(xprintf(f,"X$iz$"));
		CHKERR(m_nodes[i]->box->node()->sv_path());
		CHKERR(f->sn("\n", 1));
	}
	for (int i=0; i<=n_box(); i++) {
		GRBX_Node * nd = m_nodes[i];
		for (int j=0; j<nd->ni; j++) {
			int r, ag = nd->iarg[j];
			switch(ag&0xffff) {
				case 0xfffc:
					CHKERR(xprintf(f, "X$+%d:%d=%s\n", i, j, dbl2str_s(0,m_cs.v(ag>>16))));
					break;
				case 0xfffe:
					CHKERR(xprintf(f, "X$+%d:%d<%d\n", i, j, ag>>16));
					break;
				default:
					CHKERR(xprintf(f, "X$+%d:%d<%d:%d\n", i, j, ag&0xffff, ag>>16));
					break;
			}}}
	return r;
}

void GraphBoxGen::decr_no(int i, int no) {
	for (int j=i+1, nb = n_box(); j<=nb; j++) {
		GRBX_Node * q = m_nodes[j];
		for (int k=0,*ap=q->iarg; k<q->ni; k++,ap++) 
			(*ap&0xffff)==i && (*ap>>16)>=no && (*ap=0x5fffc); }}

void GraphBoxGen::notify_nio(BoxGen * bx) {
	IFDBGX(GRTMP) log("gr%d/notify_nio: %d", id(), bx->id());
	int nb = n_box();
	unset_model();
	for (int i=0; i<nb; i++) {
		GRBX_Node * p = m_nodes[i];
		if (p->box != bx) continue;
		int ni = bx->n_in(), no = bx->n_out();
		if (ni==p->ni && no==p->no) continue;
		if (no<p->no) decr_no(i, no);
		int cp=min_i(ni,p->ni), *av=new int[ni+no];
		if (cp) memcpy(av, p->iarg, cp*sizeof(int));
		for (int j=cp; j<ni; j++) av[j] = 0x5fffc;
		delete[] (p->iarg); p->iarg = av;
		p->oarg = av + ni; p->ni = ni; p->no = no;
	}
	if (wnfl()) box_window(); 
}

#define CH(X) BXCMD_H(GraphBoxGen, X)
CH(uic) { return p->uicfg(s[1]); }
CH(dbx) { return p->del_box(atoi(s+1)); }
CH(shfl){ return p->shuffle(), BCR_WIN; }
CH(dbg) { return p->debug(), 0; }
CH(ni)  { int k = p->n_in();  intv_cmd(&k, s+1, 0, 30-  p->m_n_fb); return p->set_n_in (k+  p->m_n_fb); }
CH(no)  { int k = p->n_out(); intv_cmd(&k, s+1, 1, 30-2*p->m_n_fb); return p->set_n_out(k+2*p->m_n_fb); }

CH(nfb) { int k = p->m_n_fb;
	  if (!intv_cmd(&k, s+1, 0, k + min_i(30-p->m_n_in, (30-p->m_n_out)>>1))) return 0;
	  int kd = k - p->m_n_fb; if (!kd) return 0;
	  return p->set_n_in(p->m_n_in+kd), p->set_n_out(p->m_n_out+2*kd), p->m_n_fb=k, BCR_WIN; }

CH(conn){ int t_b, t_i, s_b, s_i, scix = 0;
	  double v = 0.0;
	  int k = sscanf(s+1," %d : %d < %d : %d", &t_b,&t_i,&s_b,&s_i);
	  if (k==4) return p->conn(t_b,t_i,s_b,s_i);
	  if (k==3) return p->conn(t_b,t_i,-1,s_b);
	  if (sscanf(s+1, " %d : %d = %n", &t_b, &t_i, &scix)<2 || !scix) return BXE_PARSE;
	  parse_num(&v, s+1+scix); return p->conn(t_b,t_i,-2,0,v); }

CH(ibx) { int k = s[1]=='z' ? p->n_box() : atoi(s+1);
	  char * path = cb->tok(); if (!path) return BXE_NOARG;
	  BoxGen* bx; if (!(bx = cb->lookup_b(path))) return BXE_ARGLU;
	  return p->add_box(k, bx); }

CH(gnd) { SEQCHK; BoxGen * bx = cb->lookup_b(s+3); if (!bx) return BXE_ARGLU;
	  if (!cb->cperm(DF_EDBOX)) return NDE_PERM;
	  return (p->m_op && p->m_sel[1]>=0) ? p->rplc_box(p->sel_ix(1), bx) : p->add_box(p->n_box(), bx); }

CH(ck) {  SEQCHK; int i = 32*b32_to_i(s[3]) + b32_to_i(s[4]),
	  	      t = (s[5]==42) ? 1 : 2*(s[5]>79), j = (s[5]-16)&31;
	  if (*s=='1') return p->sel3(t, i, j), 0;
	  if (t==1 && i<p->n_box()) return p->m_nodes[i]->box->node()->draw_window(16), 0;
	  return EEE_NOEFF;  }

CH(rpl){  p->m_op = s[1]&1;  if (p->wnfl()) p->w_rpflg(); return 0; }

CH(g2) {SEQCHK; if (!cb->cperm(DF_EDBOX)) return NDE_PERM;
	short * sel = p->m_sel;
	int i, j, ec;
	switch (*(s+=3)) {
		case 'X': return (sel[1]<0) ? BXE_GBXSEL : p->del_box(p->sel_ix(1));
		case '+': return ((sel[0]|sel[2])<0) ? BXE_GIOSEL :
			         p->conn(p->sel_ix(0), sel[0]&31, p->sel_ix(2), sel[2]&31);
		case '=': 
			  if (sel[0]<0) return BXE_GIOSEL;
			  i = p->sel_ix(0), j = p->m_sel[0] & 31;
			  if ((p->m_nodes[i]->iarg[j] & 0xffff) == 0xfffc) return EEE_NOEFF;
			  return p->conn(i, j, -2, 0, p->m_guitmp);
		case 'v': 
			  p->m_guitmp = at0f(s+1);
			  if (sel[0] < 0) return 0;
			  i = p->sel_ix(0), j = sel[0] & 31;
			  if ((p->m_nodes[i]->iarg[j] & 0xffff) != 0xfffc) return 0;
			  if ((ec = p->conn(i, j, -2, 0, p->m_guitmp)) > 0) p->upd_con(i,j);
			  return 0;
		default:  return BXE_CENUM;
	}}

BXCMD_DEF(GraphBoxGen) { {8192+'\\',0}, {'<',c_ni}, {'>',c_no}, {'@',c_nfb}, {'+',c_conn},
	{'i',c_ibx}, {'d',c_dbx}, {'s',c_shfl}, {'d'|256,c_dbg}, {'B'|256,c_gnd}, {'G',c_g2},
	{'1'|256,c_ck}, {'3'|256,c_ck}, {'U',c_uic}, {'r',c_rpl}, {0,0} };

void graph_init() { GraphBoxGen::cmd_init(); }
int setbox_graph(ABoxNode * nd, BoxGen * _) { nd->m_box = new GraphBoxGen(nd); return 1; }
int dot_dead() { return Dot::dot_dead(); }
