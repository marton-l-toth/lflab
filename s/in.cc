#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "glob.h"
#include "box0.h"
#include "mx.h"
#include "errtab.inc"
#include "midi.h"

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
		static scf_t sc_f0, sc_f1;
		AVlimBox() : BoxInst(sc_f0) {}
		virtual int ini(double ** inb) = 0;
		virtual int set_trg(double * to) = 0;
	protected:
		char m_no, m_vlpos, m_vlstep, m_rsrv;
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
//? all inputs are expected to be constant (vl<i> only constant per
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

//? {{{!._vlm}}}
//? input box: speed-limited MIDI input box
//? dev# - (logical) device ID (0..31)
//? chan - MIDI channel (0..15)
//? base - base key/control ID
//? [i0] - list of key/ctrl id (relative to base) for sliders 0-5
//? [i6] - list of key/ctrl id (relative to base) for sliders 6-11
//? vl<i> - speed limit (1/s) for slider #i
//? a speed limit of 0 (or above the samp. rate) means no limit
//? ---
//? inputs should be constant, record to track does not work, see
//? ==> .!b.in.vls.vls01 -- (last 2 paragraphs)

class MidiVLBox : public AVlimBox {
	public: 
		MidiVLBox(int x);
		~MidiVLBox() { free(m_ix); }
		virtual int ini(double **inb);
 		virtual int set_trg(double *q) { const unsigned int*p = blk(); for (int i=0,no=m_no; i<no; i++)
						 q[i]=.007874015748031496*(double)(p[m_ix[i]]&127); return 0;}
	protected:
		inline const unsigned int * blk() { return mi_getblk(m_ldev, m_ch); }
		signed char m_ldev, m_ch;
		unsigned char * m_ix;
};

void GuiIn::w_slider(double **inb, int n, int bv) {
	double v[n]; int ix = 0;
	BVFOR_JM(bv) if (v[ix] = inb[j][0], ++ix>n) break;
	gui_sliderwin(MX_L_WIN(mxid), n, v, dat);
}

BX_SCALC(AVlimBox::sc_f0) { SCALC_BXI(AVlimBox); bxi->ini(inb); CALC_FW(sc_f1); }

BX_SCALC(AVlimBox::sc_f1) {
	SCALC_BXI(AVlimBox);
	int r=0, no = bxi->m_no, vstp = bxi->m_vlstep;
	double **vlim=inb+(int)bxi->m_vlpos, *pst = bxi->m_st, trg[no];
	int ec=bxi->set_trg(trg); if(ec) return ec;
	for (int k=0,vi=0,m=1; k<no; k++, m*=2, vi+=vstp) {
		double vl = vlim[vi][0] * sample_length, *q = outb[k], 
		       y = pst[k], t = trg[k], dif = t-y, adif = fabs(dif);
		if (vl<1e-9) vl = 2.0;
		if (adif<vl) { *q = pst[k] = trg[k]; continue; }
		int n2 = (int)floor(adif/vl);  r |= m;  if (dif<0.0) vl = -vl;
		if (n2<0) bug("n2=%d, adif=%.15g, vl=%.15g no=%d k=%d y=%.15g t=%.15g", n2, adif, vl, no, k, y, t);
		if (n2>=n) { for (int i=0;  i<n;  i++) q[i]=(y+=vl);   pst[k]=y; }
		else { 	     for (int i=0;  i<n2; i++) q[i]=(y+=vl);
			     for (int i=n2; i<n;  i++) q[i]=t;         pst[k]=t; }}
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

MidiVLBox::MidiVLBox(int x) : m_ldev(99) { m_no=(x&=31); int y=(x+7)&~7; m_vlpos=4+(x>6); m_vlstep=1;
	 				   m_ix = (unsigned char*)malloc(y+8*x); m_st = (double*)(m_ix+y); }

int MidiVLBox::ini(double **inb) {
	m_ldev = (int)lround(**inb) & 31; m_ch = (int)lround(*inb[1]) & 15;
	int bs = (int)lround(*inb[2]);
		    NAN_UNPK_32(i0, inb[3],0); for (int i=0;i<6;i++) m_ix[i  ] = (unsigned char)(bs+i0_x[i]);
	if(m_no>6){ NAN_UNPK_32(i6, inb[4],0); for (int i=0;i<6;i++) m_ix[i+6] = (unsigned char)(bs+i6_x[i]); }
	return set_trg(m_st); }

void b_in_init(ANode *r) {
	qmb_arg_t qa = QMB_ARG1(VLSBox);
	ANode *dv = qmk_dir(r, "vls");
	char nm[8]; memcpy(nm, "vls01", 6);
	qmk_box(dv, nm, qa, 1, 4, 1, "slg", "i-i:R*1", 2, "titl$[s0]", 2, 2, "lb$vl", "%zzz%%"); ++nm[4];
	for (int i=2; i<=7; i++) qmk_box(dv, nm, qa, i, 2+2*i, i, "slg", "1"), ++nm[4];
	qmk_box(dv, nm, qa, 8, 19, 8, "slg", "1i-i:R1", 3, "titl$[s0]$[s1]", 3, 3, "vl$lb");
	for (int i=9; i<=13; i++) bump_dec2(nm+3), qmk_box(dv, nm, qa, i, 3+2*i, i, "slg", "1");
	dv = qmk_dir(r, "vlm"); memcpy(nm, "vlm01", 6); qa = QMB_ARG1(MidiVLBox);
 	qmk_box(dv, nm, qa, 1, 5, 1, "vlm", "i-i.R*1", 4, "dev#$chan$base$[i0]", 4, "vl", "%zzz%%"); ++nm[4];
	for (int i=2; i<=6; i++) qmk_box(dv, nm, qa, i, 4+i, i, "vlm", "1"), ++nm[4];
 	qmk_box(dv, nm, qa, 7, 11, 7, "vlm", "1i-i.R1", 0x401, "[i6]", 5, "vl");
	for (int i=8; i<=12; i++) bump_dec2(nm+3), qmk_box(dv, nm, qa, i, 5+i, i, "vlm", "1");
}
