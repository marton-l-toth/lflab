#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "glob.h"
#include "box0.h"
#include "mx.h"
#include "errtab.inc"

#define IBF_STARTED 1
#define IBF_DETACH  2

class AInputBox : public BoxInst {
	public:
		AInputBox() : m_flg(0), m_mxid(0), m_dat(0) {}
		virtual ~AInputBox() { if (m_mxid>0 && mx_l_op(m_mxid,255,2)<0) gui_errq_add(RTE_MXLDFAIL); }
		virtual int ini(double ** inb) = 0;
		virtual int calc2(double **pto, int n) = 0;
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		void w_slider(double **inb, int n, int bv);
		int m_flg, m_mxid;
		unsigned char * m_dat;
};

static int vlim_calc(double *to, double *pos, double trg, double vlim, double vlim_1, int n) {
	double y = *pos;
	int sg = approx_cmp(trg, y); if (!sg) return *to=*pos, 0;
	double v = sg<0 ? -vlim : vlim, n2 = (sg<0 ? y-trg : trg-y) * vlim_1;
	if (n2 > (double)n) { for (int i=0; i<n; i++) to[i] = (y += v);  *pos = y; return 1; }
	int n2i = (int)floor(n2);
	for (int i=0;   i<n2i; i++) to[n] = (y += v);
	for (int i=n2i; i<n;   i++) to[i] = trg;
	*pos = trg; return 1;
}

int AInputBox::calc(int inflg, double** inb, double** outb, int n) {
	if (!(m_flg & IBF_STARTED)) {
		if ((m_mxid = mx_mklive(this)) < 0) return m_mxid;
		if (!(m_dat = mx_l_dat(m_mxid))) return RTE_BUG;
		int ec = ini(inb); if (ec<0) return ec; else m_flg |= IBF_STARTED;
	}
	if (!(m_flg & IBF_DETACH) && !(mx_l_op(m_mxid,255,3)&MXLF_WIN)) return RTE_IWCLOSE;
	return calc2(outb, n);
}

void AInputBox::w_slider(double **inb, int n, int bv) {
	double v[n]; int ix = 0;
	BVFOR_JM(bv) if (v[ix] = inb[j][0], ++ix>n) break;
	gui_sliderwin(MX_L_WIN(m_mxid), n, v, m_dat);
}

//? {{{!._slg}}}
//? input box: speed-limited slider group
//? titl - window title
//? [s0] - start pos (0...100) for sliders 0-6 
//? [s7] - start pos (0...100) for sliders 7-12
//? lb<i> - label for slider #i
//? vl<i> - speed limit (1/s) for slider #i
//? a speed limit of 0 (or above the samp. rate) means no limit
//? if [s0] or [s7] are shorter on non-lists, the default is 50
//? all inputs are expected to be constant

//? {{{!._slgA}}}
//? input box: speed-limited & smoothed slider group
//? titl - window title
//? [s0] - start pos (0...100) for sliders 0-6
//? [s7] - start pos (0...100) for sliders 7-8
//? lb<i> - label for slider #i
//? vl<i> - speed limit (1/s) for slider #i
//? dmp<i> - dampening: 0=unfiltered, otherwise as for =accN(md=0)
//? a speed limit <= 0 (or above the samp. rate) means no limit
//? if [s0] or [s7] are shorter on non-lists, the default is 50
//? all inputs are expected to be constant

class VLSBox : public AInputBox {
	public:
		VLSBox(int x) : m_arg(x), m_no(x&31) {
			m_store = (double*)malloc( (24+((x&64)>>2)) * m_no ); }
		~VLSBox() { free(m_store); }
		virtual int ini(double ** inb);
		virtual int calc2(double **pto, int n);
	protected:
		int m_arg, m_no, m_lim_bv, m_acc_bv;
		double *m_store;
};

int VLSBox::ini(double **inb) {
	int no = m_no, ff = !!(m_arg & 64), ac1 = (no>7-2*ff);
	double *pos = m_store, *vlim = pos+no, *vlim_1 = vlim+no;
	long long nan;
	memset(m_dat, 50, no);
	memcpy(&nan, inb[1], 8), nan_unpk((char*)m_dat,   0, nan, 0);
	if (ac1) memcpy(&nan, inb[2], 8), nan_unpk((char*)m_dat+7, 0, nan, 0);
	w_slider(inb, no, ac1 ? 0x2aaaaaa9 : 0x5555);
	m_lim_bv = m_acc_bv = 0;
	for (int i=0; i<no; i++) {
		pos[i] = .01 * (double)m_dat[i];
		double x = inb[2+ac1+2*i][0];
		if (x>1e-8 && (x*=sample_length)<1.0) vlim[i]=x, vlim_1[i]=1.0/x, m_lim_bv|=(1<<i);
	}
	return 0;
}
// static int vlim_calc(double *to, double *pos, double trg, double vlim, double vlim_1, int n) {

int VLSBox::calc2(double **pto, int n) {
	int no = m_no;
	double *pos = m_store, *vlim = pos+no, *vlim_1 = vlim+no;
	int rv = 0, cbv = ((1<<no)-1) ^ m_lim_bv;
	BVFOR_JM(cbv) pto[j][0] = .01 * m_dat[j];
	BVFOR_JM(m_lim_bv) rv |= vlim_calc(pto[j], pos+j, .01*m_dat[j], vlim[j], vlim_1[j], n) << j;
	if (!(m_arg&64)) return rv;
	return RTE_BUG; // TODO: acc
}

#define LVD(J) "$lb" #J "$vl" #J "$dmp" #J
void b_in_init(ANode *r) {
	qmb_arg_t qa = QMB_ARG1(VLSBox);
	ANode *dv = qmk_dir(r, "vls"), *dva = qmk_dir(r, "vlsA");
	char nm[8]; memcpy(nm, "vls01", 6);
	qmk_box(dv, nm, qa, 1, 4, 1, "slg", "i-i:R*1", 2, "titl$[s0]", 2, 2, "lb$vl", "z%%%%z"); ++nm[4];
	for (int i=2; i<=7; i++) qmk_box(dv, nm, qa, i, 2+2*i, i, "slg", "1"), ++nm[4];
	qmk_box(dv, nm, qa, 8, 13, 8, "slg", "i-i:R*1", 3, "titl$[s0]$[s1]", 3, 3, "vl$lb", "z%%%%z");
	for (int i=9; i<=13; i++) (i==10?(nm[3]=49,nm[4]=48):++nm[4]), 
				  qmk_box(dv, nm, qa, i, 3+2*i, i, "slg", "1");
}

