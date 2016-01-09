#ifndef __qwe_combo_h__
#define __qwe_combo_h__

#include "box0.h"

#define AX_TY_IN   0x1000
#define AX_TY_OUT  0x2000
#define AX_TY_CON  0x3000
#define AX_TY_TMP  0x4000
#define AX_TY_JUNK 0x5000
#define AX_TY_MASK 0x7000
#define AX_IX_MASK 0x0fff

class ComboBoxModel : public BoxModel {
        public: 
                ComboBoxModel(int nb, int nc, int nt, int nio);
                virtual ~ComboBoxModel();
		virtual BoxInst * mk_box();
		inline void dsc_nio(int i, int o) { *(eob++) = i, *(eob++) = o; }
		inline void dsc_arg(int a) { int t = (a+1)&3, j = a>>16;
			(t&2) ? dsc_nio(3-t, j) : dsc_nio(2+((t-1)&n_cb)+(j>>5), j&31); }
		void dump() {} // TODO 
		int n_bx, n_cb, n_t, niof;
		double *pcon, ***pppcon;
		BoxModel ** boxm;
		unsigned char *iolist, *eob;
};

class FeedbackBoxInst : public BoxInst {
	public:
		FeedbackBoxInst(BoxInst * sub, int niof);
		virtual ~FeedbackBoxInst() { delete(m_sub); }
		virtual int calc(int inflg, double** inb, double** outb, int n);
	protected:
		int calc1(int inflg, double** inb, double** outb, int n);
		int calc2(int inflg, double** inb, double** outb, int i0, int n);
		BoxInst * m_sub;
		int m_ni, m_no, m_nf, m_chs, m_omsk;
		double * m_buf;
		int * m_tjb;
};

#endif // __qwe_combo_h__
