#include "util.h"
#include "combo.h"
#include "glob.h"

class ComboBoxInst : public BoxInst {
        public: 
                ComboBoxInst(ComboBoxModel * model);
                virtual ~ComboBoxInst();
                virtual int calc(int inflg, double** inb, double** outb, int n);
        protected:
                BoxInst ** m_bxpp;
                ComboBoxModel * m_m;
};

FeedbackBoxInst::FeedbackBoxInst(BoxInst * sub, int niof) : m_sub(sub), m_ni(niof>>16), m_no((niof>>8)&31),
	m_nf(niof&15), m_chs(-1), m_omsk(((1<<m_no)-1) | ((0xaaaaaaa&((1<<(2*m_nf))-1))<<m_no)),
	m_buf(0), m_tjb(0) {}

int FeedbackBoxInst::calc(int ifg, double** ib, double** ob, int n) {
	int i = 0; if (m_chs<0 && (!n || (i=calc1(ifg, ib, ob, n))<0)) return i;
	while (i<n) i += calc2(ifg, ib, ob, i, n-i);   return (1<<m_no)-1;  }

int FeedbackBoxInst::calc1(int inflg, double** inb, double** outb, int n) {
	int ni1 = m_ni+m_nf, no1 = m_no+2*m_nf, tlim = 10*sample_rate,
	    bix[m_nf], idly[m_nf];
	double dly[m_nf], y0[m_nf], *pin[ni1], *pou[no1];
	memcpy(pin, inb,  m_ni*sizeof(double*));
	memcpy(pou, outb, m_no*sizeof(double*));
	for (int i=0; i<m_nf; i++) pin[m_ni+i] = zeroblkD;
	for (int i=0; i<m_nf; i++) pou[m_no+2*i] = dly + i, pou[m_no+2*i+1] = y0 + i;
	int r = m_sub->calc(inflg | ((1<<ni1) - (1<<m_ni)), pin, pou, 1); if (r<1) return r;
	m_chs = 512;
	for (int i=0; i<m_nf; i++) {
		int t = (int)lround((double)sample_rate * dly[i]);
		if (t<1) t = 1; else if (t>tlim) t = tlim;
		if (t<m_chs) m_chs = t;
		idly[i] = t;
	}
	int bufs = 0, bs0 = 2*m_chs - 1;
	for (int i=0; i<m_nf; i++) bix[i] = bufs, bufs += idly[i] + bs0;
	char * buf = (char*) malloc(bufs*sizeof(double) + 3*m_nf*sizeof(int));
	m_buf = (double*)buf; m_tjb = (int*)(buf+bufs*sizeof(double));
	for (int i=0, *p=m_tjb; i<m_nf; i++, p+=3) {
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

ComboBoxModel::ComboBoxModel(int nb, int nc, int nt, int nio) : n_bx(nb), n_cb((nc+31)>>5), n_t(nt) {
	int opc = 8*nc, oppc = opc + nc*sizeof(double*), obm = oppc + n_cb*sizeof(double**),
	    odsc = obm + nb*sizeof(BoxModel*), blksiz = odsc + 2*(nb+nio);
	char * blk = (char*) malloc(blksiz);
	pcon = (double*)blk; pppcon = (double***)(blk+oppc); boxm = (BoxModel**)(blk+obm); 
	iolist = eob = (unsigned char*)blk+odsc;
	double **ppcon = (double**)(blk+opc);
	for (int i=0; i<nc; i++) ppcon[i] = pcon + i;
	for (int i=0; i<n_cb; i++) pppcon[i] = ppcon + 32*i;
}

ComboBoxModel::~ComboBoxModel() { free(pcon); }

BoxInst * ComboBoxModel::mk_box() {
	ComboBoxInst * cb = new ComboBoxInst(this); // log("combo: niof = 0x%x", niof);
	return (niof&31) ? (BoxInst*)new FeedbackBoxInst(cb, niof) : (BoxInst*)cb;
}

ComboBoxInst::ComboBoxInst(ComboBoxModel * model) : m_m(model) {
	int nb = m_m->n_bx; if (!nb) { m_bxpp = 0; return; }
	m_bxpp = (BoxInst**) malloc(nb*sizeof(BoxInst*));
	for (int i=0; i<nb; i++) m_bxpp[i] = model->boxm[i]->mk_box();
	BoxModel::ref(model);
}

ComboBoxInst::~ComboBoxInst() {
        for (int i=0,n=m_m->n_bx; i<n; i++) delete m_bxpp[i];
	BoxModel::unref(m_m); }

int ComboBoxInst::calc(int inflg, double** inb, double** outb, int n) {
	int n_b = m_m->n_bx, n_cb = m_m->n_cb, n_t = m_m->n_t, n_tb = (n_t+31)>>5, nblk = 2+n_cb+n_tb;
	unsigned int flg[n_tb];
	double **pblk[nblk], *ptmp[n_t], tmp[n_t*n], *p = tmp, **pp = ptmp, ***ppp = pblk + 2 + n_cb,
	       *iarg[32], *oarg[32];
	pblk[0] = inb; pblk[1] = outb; flg[0] = (unsigned int)inflg; flg[1] = 0u;
	for (int i=0; i<n_cb; i++) pblk[2+i] = m_m->pppcon[i], flg[2+i] = 0u;
	for (int i=0; i<n_t;  i++) ptmp[i] = p, p += n;
	for (int i=0; i<n_tb; i++) ppp[i] = pp, pp += 32;
	BoxInst ** bxpp = m_bxpp; const unsigned char *s = m_m->iolist;
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
