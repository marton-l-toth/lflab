#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "glob.h"
#include "box0.h"
#include "mx.h"
#include "errtab.inc"

#define IBF_STARTED 1
#define IBF_DETACH  2

struct GuiIn {
	GuiIn() : flg(0) {}
	int flg, mxid; unsigned char * dat;
	int ini(BoxInst* bxi) { if ((mxid = mx_mklive(bxi)) < 0) return mxid;
				return (dat = mx_l_dat(mxid)) ? (flg|=IBF_STARTED,0) : RTE_BUG; }
	void bye() { int e; if (mxid>0 && (e=mx_l_op(mxid,255,2))<0) gui_errq_add2(e,RTE_MXLDFAIL); }
	inline int cls() { return (!(flg&IBF_DETACH)) && !(mx_l_op(mxid,255,3)&MXLF_WIN); }
	void w_slider(double **inb, int n, int bv);
};

class AVlimBox : public BoxInst {
	public:
		AVlimBox() : m_run(0) {}
		virtual int ini(double ** inb) = 0;
		virtual int set_trg(double * to) = 0;
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		char m_no, m_vlpos, m_vlstep, m_run;
		double * m_st;
};

//? {{{!._slg}}}
//? input box: speed-limited slider group (GUI)
//? titl - window title
//? [s0] - start pos (0...100) for sliders 0-6 
//? [s7] - start pos (0...100) for sliders 7-12
//? lb<i> - label for slider #i
//? vl<i> - speed limit (1/s) for slider #i
//? a speed limit of 0 (or above the samp. rate) means no limit
//? if [s0] or [s7] are shorter on non-lists, the default is 50
//? ---
//? all inputs are expected to be constant (vl<i> only consantt per
//? time-slice so they can be connected to an output of a different
//? (speed-unlimited) control box)
//? ---
//? limitation: input boxes are not compatible with recording
//? to track (in the recorded track, the sound will be as if
//? all sliders were left on starting position)

class VLSBox : public AVlimBox {
	public:
		VLSBox(int x);
		~VLSBox() { m_gi.bye(), free(m_st); }
		virtual int ini(double ** inb);
		virtual int set_trg(double *to);
	protected:
		GuiIn m_gi;
		int m_arg;
};

void GuiIn::w_slider(double **inb, int n, int bv) {
	double v[n]; int ix = 0;
	BVFOR_JM(bv) if (v[ix] = inb[j][0], ++ix>n) break;
	gui_sliderwin(MX_L_WIN(mxid), n, v, dat);
}


int AVlimBox::calc(int inflg, double** inb, double** outb, int n) {
	if (!m_run) m_run=1, ini(inb);
	int ec, r=0, no = m_no, vstp = m_vlstep;
	double **vlim=inb+(int)m_vlpos, trg[no];  if ((ec=set_trg(trg))) return ec;
	for (int k=0,vi=0,m=1; k<no; k++, m*=2, vi+=vstp) {
		double vl = vlim[vi][0] * sample_length, *q = outb[k], 
		       y = m_st[k], t = trg[k], dif = t-y, adif = fabs(dif);
		if (vl<1e-9) vl = 2.0;
		if (adif<vl) { *q = m_st[k] = trg[k]; continue; }
		int n2 = (int)floor(adif/vl);  r |= m;  if (dif<0.0) vl = -vl;
		if (n2<0) bug("n2=%d, adif=%.15g, vl=%.15g no=%d k=%d y=%.15g t=%.15g", n2, adif, vl, no, k, y, t);
		if (n2>=n) { for (int i=0;  i<n;  i++) q[i]=(y+=vl);   m_st[k]=y; }
		else { 	     for (int i=0;  i<n2; i++) q[i]=(y+=vl);
			     for (int i=n2; i<n;  i++) q[i]=t;          m_st[k]=t; }}
	for (int k=0; k<no; k++) log_n(" %d:%g->%g", k, m_st[k], trg[k]); log("");
	return r;
}

VLSBox::VLSBox(int x) : m_arg(x) { m_no=(x&=31); m_vlpos=3+(x>7); m_vlstep=2; m_st=(double*)malloc(8*x); }

int VLSBox::ini(double **inb) {
	int ec = m_gi.ini(this); if (ec<0) return ec;
	int no = m_no, ac1 = (no>7);
	long long nan; memset(m_gi.dat, 50, no);
		 memcpy(&nan, inb[1], 8), nan_unpk((char*)m_gi.dat,   0, nan, 0);
	if (ac1) memcpy(&nan, inb[2], 8), nan_unpk((char*)m_gi.dat+7, 0, nan, 0);
	m_gi.w_slider(inb, no, ac1 ? 0x2aaaaaa9 : 0x5555);
	for (int i=0; i<no; i++) m_st[i] = .01*(double)m_gi.dat[i];
	return 0;
}

int VLSBox::set_trg(double * to) { 
	if (m_gi.cls()) return RTE_IWCLOSE; unsigned char * dat = m_gi.dat; int no = m_no;
	for (int i=0; i<no; i++) to[i] = .01 * (double)dat[i];    return 0; }

void b_in_init(ANode *r) {
	qmb_arg_t qa = QMB_ARG1(VLSBox);
	ANode *dv = qmk_dir(r, "vls");
	char nm[8]; memcpy(nm, "vls01", 6);
	qmk_box(dv, nm, qa, 1, 4, 1, "slg", "i-i:R*1", 2, "titl$[s0]", 2, 2, "lb$vl", "%zzz%%"); ++nm[4];
	for (int i=2; i<=7; i++) qmk_box(dv, nm, qa, i, 2+2*i, i, "slg", "1"), ++nm[4];
	qmk_box(dv, nm, qa, 8, 19, 8, "slg", "1i-i:R1", 3, "titl$[s0]$[s1]", 3, 3, "vl$lb");
	for (int i=9; i<=13; i++) (i==10?(nm[3]=49,nm[4]=48):++nm[4]), 
				  qmk_box(dv, nm, qa, i, 3+2*i, i, "slg", "1");
}
