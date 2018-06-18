#include <fcntl.h>

#include "pt.h"
#include "util.h"
#include "combo.h"
#include "glob.h"

class FeedbackBoxInst : public BoxInst {
        public:
		static scf_t sc_first, sc_rest;
                FeedbackBoxInst(BoxInst * sub, int niof);
                virtual ~FeedbackBoxInst() { delete(m_sub); }
        protected:
                int calc1(int inflg, double** inb, double** outb);
                int calc2(int inflg, double** inb, double** outb, int i0, int n);
                BoxInst * m_sub;
                int m_ni, m_no, m_nf, m_chs, m_omsk;
                double * m_buf;
                int * m_tjb;
};

class ComboBoxInst : public BoxInst {
        public:
		static scf_t sc_f;
                ComboBoxInst(ComboBoxModel * model);
                virtual ~ComboBoxInst();
        protected:
                BoxInst ** m_bxpp;
                ComboBoxModel * m_m;
};

FeedbackBoxInst::FeedbackBoxInst(BoxInst * sub, int niof) : 
	BoxInst(sc_first), m_sub(sub), m_ni(niof>>16), m_no((niof>>8)&31),
	m_nf(niof&15), m_chs(-1), m_omsk(((1<<m_no)-1) | ((0xaaaaaaa&((1<<(2*m_nf))-1))<<m_no)),
	m_buf(0), m_tjb(0) {}

BX_SCALC(FeedbackBoxInst::sc_first) {
	SCALC_BXI(FeedbackBoxInst); bxi->m_psc = sc_rest; 
	int i = bxi->calc1(inflg, inb, outb); if (i<0) return i;
	while (i<n) i += bxi->calc2(inflg, inb, outb, i, n-i);   return (1<<bxi->m_no)-1;  }

BX_SCALC(FeedbackBoxInst::sc_rest) {
	SCALC_BXI(FeedbackBoxInst); int i = 0;
	while (i<n) i += bxi->calc2(inflg, inb, outb, i, n-i);   return (1<<bxi->m_no)-1;  }

int FeedbackBoxInst::calc1(int inflg, double** inb, double** outb) {
	int nf = m_nf, ni1 = m_ni+nf, no1 = m_no+2*nf, tlim = 10*sample_rate,
	    bix[nf], idly[nf];
	double dly[nf], y0[nf], *pin[ni1], *pou[no1];
	memcpy(pin, inb,  m_ni*sizeof(double*));
	memcpy(pou, outb, m_no*sizeof(double*));
	for (int i=0; i<nf; i++) pin[m_ni+i] = zeroblkD;
	for (int i=0; i<nf; i++) pou[m_no+2*i] = dly + i, pou[m_no+2*i+1] = y0 + i;
	int r = m_sub->calc(inflg | ((1<<ni1) - (1<<m_ni)), pin, pou, 1); if (r<1) return r;
	m_chs = 512;
	for (int i=0; i<nf; i++) {
		int t = (int)lround((double)sample_rate * dly[i]);
		if (t<1) t = 1; else if (t>tlim) t = tlim;
		if (t<m_chs) m_chs = t;
		idly[i] = t;
	}
	int bufs = 0, bs0 = 2*m_chs - 1;
	for (int i=0; i<nf; i++) bix[i] = bufs, bufs += idly[i] + bs0;
	char * buf = (char*) malloc(bufs*sizeof(double) + 3*nf*sizeof(int));
	m_buf = (double*)buf; m_tjb = (int*)(buf+bufs*sizeof(double));
	for (int i=0, *p=m_tjb; i<nf; i++, p+=3) {
		int t = p[0] = idly[i];
		double * q = m_buf + (p[2]=bix[i]);
		memset(q, 0, 8*t); q[t] = y0[i]; p[1] = 1;
	}
	return 1;
}

int FeedbackBoxInst::calc2(int inflg, double** inb, double** outb, int i0, int n) {
	int k, ni1 = m_ni+m_nf, no1 = m_no+2*m_nf, c = m_chs, cpn[m_nf];
	double *pin[ni1], *pou[no1], *cpt[m_nf], *cps[m_nf];
	memset(cpn, 0, m_nf*sizeof(int));
	if (n>c) n = c;
	for (int i=0, f=inflg; i<m_ni; i++,f>>=1) pin[i] = inb[i] + (i0 &- (f&1));
	for (int i=0; i<m_no; i++) pou[i] = outb[i] + i0;
	for (int i=0; i<m_nf; i++) {
		int *p = m_tjb + 3*i, t = p[0], j = p[1], jj = (j>=c ? j-c : j+t);
		double * q = m_buf + p[2];
		pin[m_ni+i] = q + j; pou[m_no+2*i] = junkbuf; pou[m_no+2*i+1] = q + jj;
		if ((k=j+n-c-t) > 0) memcpy(q+c+t, q, 8*k);
		else if ((k=jj+n-c-t) > 0) cpn[i]=8*k, cpt[i]=q, cps[i]=q+c+t;
		j += n; if ((k=j-c-t)>=0) j = k; p[1] = j;
	}
	int r = m_sub->calc(inflg | ((1<<ni1) - (1<<m_ni)), pin, pou, n);
	if (r<0) return r; else BoxInst::rmcon(~r&m_omsk, pou, n);
	for (int i=0; i<m_nf; i++) if (cpn[i]) memcpy(cpt[i], cps[i], cpn[i]);
	return n;
}

ComboBoxModel::ComboBoxModel(ConStore * cs, char * q, int nb) : n_bx(nb), pcs(cs) {
	int nc = cs->n(), ncb = (nc+63)>>5, o1=nc*sizeof(void*), o2 = o1+ncb*sizeof(void*);
	double **ppcon = (double**)q;   pppcon = (double***)(q+o1);  boxm = (BoxModel**)(q+o2);
	iolist = eob = (unsigned char*)(q+o2+nb*sizeof(void*)); n_cb = ncb;
	for (int i=0; i<nc; i++) ppcon[i] = pcs->p(i+32);
	*pppcon = pstat_con; for (int i=1; i<ncb; i++) pppcon[i] = ppcon + 32*(i-1);
	IFDBGX(EXPR) log("nb=%d nc=%d n_cb=%d o1=%d o2=%d", nb, nc, n_cb, o1, o2);
}

ComboBoxModel::~ComboBoxModel() { for (int i=0; i<n_bx; i++) BoxModel::unref(boxm[i]); }

BoxInst * ComboBoxModel::mk_box() {
	ComboBoxInst * cb = new ComboBoxInst(this); // log("combo: niof = 0x%x", niof);
	return (niof&31) ? (BoxInst*)new FeedbackBoxInst(cb, niof) : (BoxInst*)cb;
}

ComboBoxInst::ComboBoxInst(ComboBoxModel * model) : BoxInst(sc_f), m_m(model) {
	int nb = m_m->n_bx; if (!nb) { m_bxpp = 0; return; }
	m_bxpp = (BoxInst**) malloc(nb*sizeof(BoxInst*));
	for (int i=0; i<nb; i++) m_bxpp[i] = model->boxm[i]->mk_box();
	BoxModel::ref(model);
}

ComboBoxInst::~ComboBoxInst() {
        for (int i=0,n=m_m->n_bx; i<n; i++) delete m_bxpp[i];
	BoxModel::unref(m_m); }

BX_SCALC(ComboBoxInst::sc_f) {
	SCALC_BXI(ComboBoxInst); ComboBoxModel * mdl = bxi->m_m; BoxInst ** bxpp = bxi->m_bxpp;
	int n_b = mdl->n_bx, n_cb = mdl->n_cb, n_t = mdl->n_t, n_tb = (n_t+31)>>5, nblk = 2+n_cb+n_tb;
	unsigned int flg[n_tb];
	double **pblk[nblk], *ptmp[n_t], tmp[n_t*n], *p = tmp, **pp = ptmp, ***ppp = pblk + 2 + n_cb,
	       *iarg[32], *oarg[32];
	pblk[0] = inb; pblk[1] = outb; flg[0] = (unsigned int)inflg; flg[1] = 0u;
	for (int i=0; i<n_cb; i++) pblk[2+i] = mdl->pppcon[i], flg[2+i] = 0u;
	for (int i=0; i<n_t;  i++) ptmp[i] = p, p += n;
	for (int i=0; i<n_tb; i++) ppp[i] = pp, pp += 32;
	const unsigned char *s = mdl->iolist;
	for (int i=0; i<n_b; i++) {
		int ifg = 0, ni = *(s++), no = *(s++), omsk = (1<<no)-1, ih, il;
		for (int j=0; j<ni; j++) iarg[j] = pblk[ih=s[0]][il=s[1]],
					 ifg |= ((flg[ih]>>il)&1) << j, s += 2;
		for (int j=0; j<no; j++) oarg[j] = pblk[ih=s[2*j]][il=s[2*j+1]];
		int rv = bxpp[i]->calc(ifg, iarg, oarg, n); if (rv<0) return rv;
		BVFOR_JM( rv & omsk) flg[s[2*j]] |=  (1u << s[2*j+1]);
		BVFOR_JM(~rv & omsk) flg[s[2*j]] &= ~(1u << s[2*j+1]);
		s += 2*no;
	}
	return (int) flg[1];
}

Dot * Dot::m0_sg = 0;
Dot * Dot::sg() {
        if (!m0_sg) (m0_sg = new Dot()) -> start();
        return m0_sg; }

int Dot::dot_dead() {
        errtemp_cond("dot_dead");
        if (!m0_sg) log("BUG: unexp. dot_dead() (!sg) ");
        else close(m0_sg->m_pipe), m0_sg->m_pid = 0, delete(m0_sg), m0_sg = 0;
        return 0;
}

int Dot::start() {
        if (m_pid>0) { log("dot already started: pid %d", m_pid); return -3; }
        const char * path = "dot"; // TODO: config(?)
        m_pid = launch(path, "!8>Tt", &m_pipe, "-Tplain-ext", (char*)0);
        return (m_pid<0) ? (log("failed to start %s (ret: %d, ue: %s)\n", path,m_pid,strerror(errno)), -1) : 0;
}

void Dot::flush() {
        static char * t_dot = 0;
        IFDBGX(GRTMP) {
                if (!t_dot) t_dot = (char*)malloc(QENVL('w')+8), memcpy(t_dot, QENV('w'), QENVL('w')),
                                                                  memcpy(t_dot+QENVL('w'), "/gr.dot", 8);
                int fd = creat(t_dot, 0644);
                if (fd<0) perror(t_dot);
                else write(fd, m_buf, m_n), close(fd);
        }
        write(m_pipe, m_buf, m_n); m_n = 0;
}

void Dot::ghead(int tid, int tt, int wid, int wix, int seq, int nn, int ne) {
        if (m_n) { log("dot/ghead: n=%d, dropped", m_n); }
        char * p = m_buf;
        memcpy(p, "graph qw { node [shape=rectangle]; C", 36); p += 36;
        p += hx5(p, tid); *(p++) = '_';
        HEX2(p, tt); HEX2(p, wid);
        if (wix<0) p[0] = p[1] = 'z', p += 2; else HEX2(p, wix);
        *(p++)='_'; HEX2(p,seq); *(p++)='_'; p+=hx5(p, nn); *(p++)='_'; p+=hx5(p, ne);
        memcpy(p, " }\ndigraph lofasz { node [shape=record]; edge[arrowhead=none,arrowtail=none]\n", 77);
        m_n = (p - m_buf) + 77;
        m_nids = new char[4*nn];
}

void Dot::gnode(int ix, int ni, int no, int w) {
        char * p = m_buf + m_n; *(p++) = 'b';
        p[0] = i_to_b32((ix>>5)&31); p[1] = i_to_b32(ix&31);
        p[2] = i_to_b32(ni); p[3] = i_to_b32(no);
        memcpy(m_nids + 4*ix, p, 4); p+=4;
        memcpy(p, " [label=\"", 9); p+=9;
        int i;
        if (w==-1) {
                for (i=0; i<ni; i++) {
                        if (i) *(p++) = '|';
                        p[0]='<'; p[1]='i'; p[2] = i_to_b32(i);
                        memcpy(p+3, ">mml\\nmml", 9); p += 12;
                }
        } else if (w==-2) {
                for (i=0; i<no; i++) {
                        if (i) *(p++) = '|';
                        p[0]='<'; p[1]='o'; p[2] = i_to_b32(i);
                        memcpy(p+3, ">mml", 4); p += 7;
                }
        } else {
                *(p++) = '{';
                for (i=0; i<ni; i++) {
                        p[0] = '|'-!i; p[1] = '<'; p[2] = 'i';
                        p[3] = i_to_b32(i);
                        memcpy(p+4, ">mml", 4); p += 8;
                }
                if (ni) p[0]='}', p[1]='|', p+=2;
                p += lorem(p, w); *(p++) = 'l'; *(p++) = '|';
                for (i=0; i<no; i++) {
                        p[0] = '|' - !i;
                        p[1]='<'; p[2]='o'; p[3] = i_to_b32(i);
                        memcpy(p+4, ">mml", 4); p += 8;
                }
                p[0] = p[1] = '}'; p += 2;
        }
        memcpy(p, "\"]\n", 3); m_n = (p+3) - m_buf;
}

void Dot::gedge(int n0, int o0, int n1, int i1) {
        char * p = m_buf + m_n;
        p[0] = 'b'; memcpy(p+1, m_nids+4*n0, 4); p[5]=':';
        p[6]='o'; p[7] = i_to_b32(o0); memcpy(p+8, "->b", 3);
        memcpy(p+11, m_nids+4*n1, 4); p[15]=':'; p[16]='i'; p[17] = i_to_b32(i1);
        p[18] = '\n'; m_n += 19;
}

void Dot::gend() {
        m_buf[m_n++]='}'; m_buf[m_n++]='\n'; flush();
        delete[](m_nids); m_nids = 0; }
