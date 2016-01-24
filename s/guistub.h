#ifndef __qwe_guistub_h__
#define __qwe_guistub_h__

#include "box.h"

#define D_GUIFUN0(GF,CL,CF) inline void GF(CL *obj) { m_bufp += obj->CF(m_bufp); }
#define D_GUIFUN1(GF,CL,CF) inline void GF(CL *obj, int x) { m_bufp += obj->CF(m_bufp,x); }
int mx_c_dump_keys(char*, int), wrap_dump_keytab(char*, unsigned int *, short**);

class GuiStub {
	public:
		typedef int (*fun_t) (void*, char*, int);
		static int gui_dead(int pid, int stat, int td);
		GuiStub() : m_pid(0), m_inpipe(-1), m_outpipe(-1), m_tpipe(-1), m_lwi(-1), 
			    m_errq_n(0), m_bufp(0) { m_errq_t0[0] = m_errq_t0[1] = 0; }
		// ~GuiStub () { if (m_pid) stop(); }
		int start();
		void stop();
		int outpipe() const { return m_outpipe; }
		int tpipe() const { return m_tpipe; }
		int flush();
		void flush_all();
		int pending() const {return m_bufp!=m_buf0 || glob_flg!=m_gf0 || (m_errq_n|Node::slr_flg());}
		inline int cfl() { return ((m_bufp-m_buf0)&0xffff8000) && flush(); }
		void invd() { m_lwi = -1; }
		void brk() { *(m_bufp++)='\n'; m_lwi = -1; }
		void clear() { m_bufp = m_buf0; m_lwi = -1; }

		inline void c1(int c) { *(m_bufp++)= c; }
		inline void c2(int c, int d) { *(m_bufp++) = c; *(m_bufp++) = d; }
		inline void c3(int c, int d, int e) { *(m_bufp++) = c; *(m_bufp++) = d; *(m_bufp++) = e; }
		inline void c4(int c, int d, int e, int f) { *(m_bufp++) = c; 
			   *(m_bufp++) = d; *(m_bufp++) = e; *(m_bufp++) = f; }
		inline void hex1(int x) { *(m_bufp++) = hexc1(x); }
		void sn(const char * s, int n) { memcpy(m_bufp, s, n); m_bufp += n; }
		void sz(const char * s) { m_bufp += s__cat(m_bufp, s); }
		void hex5(int x) { m_bufp += hx5(m_bufp, x); }
		void hex4(int x) { *(int*)m_bufp = qh4(x); m_bufp += 4; }
		void hex8(int x) { hex4((int)((unsigned int)x>>16)); hex4(x&65535); }
		void hexn(int x, int w) { while (--w>=0) *(m_bufp++) = hexc1   ((x>>(4*w))&15); }
		void hexs(unsigned const char * s, int k);
		void b32n(int x, int w) { while (--w>=0) *(m_bufp++) = i_to_b32((x>>(5*w))&31); }
		void hdbl(double x) { doub2hx(m_bufp, x); m_bufp += 16; }
		void pf(const char * fmt, ...);
		void fun(fun_t pf, void *obj, int arg) { m_bufp += (*pf)(obj, m_bufp, arg); }

		D_GUIFUN0(nname,   ANode,  get_name)
		D_GUIFUN0(bxmini,  BoxGen, mini)
		D_GUIFUN0(bn_dsc,  ABoxNode,dsc) 
		D_GUIFUN1(npath,   ANode,  get_path_uf)
		D_GUIFUN1(namevec, NameVec,ls)
		void ionm(ABoxNode * bn, int io, int j) { m_bufp += bn->get_ionm(m_bufp, io, j); }
		void slr_upd() { m_bufp += Node::slr_upd_2(m_bufp); }
		void mx_keylist(int ci) { m_bufp += mx_c_dump_keys(m_bufp, ci); }
		void wr_keylist(unsigned int *bv, short ** pk) { m_bufp += wrap_dump_keytab(m_bufp, bv, pk); }
		void setwin(int oid, int ty) { m_cwi=oid; m_cwt=ty; ty<<=24; m_th=ty+0x245509; m_wh=ty+0x245709; }
		void cre(int oid, int ty, const char *s = 0);
		void t_sz(const char *s) { t0(); sz(s); }
		void t_sn(const char *s, int n) { t0(); sn(s,n); }
		void t_pf(const char * fmt, ...);

		void t0() { cfl(); m_bufp += (m_cwi==m_lwi) ? (*(int*)m_bufp = m_th, 4) : lx0('U'); }
		void w0() { cfl(); m_bufp += (m_cwi==m_lwi) ? (*(int*)m_bufp = m_wh, 4) : lx0('W'); }
		int lx0(int c) { c2(9, c); hexn(m_lwi=m_cwi, 6); c2(36, m_cwt); return 0; }
		void closewin(int id) { c2(9,'Z'); hexn(id, 6); c2(36, 63); }
		void errq_cfl();
		void wupd(int wwt, int wwix = -1) { w0(); c1(wwt); if (wwix>=0) c2(hexc1(wwix>>4),hexc1(wwix&15)); }
		void wupd_s(int wwt, const char *s, int wwix = -1) { wupd(wwt,wwix); c1('t'); sz(s); }
		void wupd_cs(int wwt, int c, const char *s, int wwix = -1) { wupd(wwt,wwix); c1(c); sz(s); }
		void wupd_ct(int wwt, const char *s, int wwix = -1) { wupd(wwt,wwix); c2('t',','); sz(s); }
		void wupd_0(int wwt, const char *s, int wwix = -1) { wupd(wwt,wwix); sz(s); }
		void wupd_c0(int wwt, int c, int wwix = -1) { wupd(wwt,wwix); c1(c); }
		void wupd_i(int wwt, int x, int wwix = -1) { wupd(wwt,wwix); c1('x'); if (x&~0xfffff) hexn(x,8); else hex5(x); }
		void wupd_i1(int wwt, int x, int wwix = -1) { wupd(wwt,wwix); c2('x', hexc1(x)); }
		void wupd_c48(int wwt, int x, int wwix = -1) { wupd(wwt,wwix); c2('c', x+48); }
		void wupd_ls(int wwt, int x, int wwix = -1) { wupd(wwt,wwix); c2('+', x+48); }
		void wupd_i2(int wwt, int x, int wwix = -1) { wupd(wwt,wwix); c1('x'); c2(hexc1(x>>4), hexc1(x&15)); }
		void wupd_d(int wwt, double x, int wwix = -1) { wupd(wwt,wwix); c1('@'); hdbl(x); }

		int tree_expand(int lr, ADirNode * dir);
		void root_expand() { tree_expand(0, Node::root()); tree_expand(1, Node::root()); }
		void tree_force(int lr, ADirNode * dir) {
			if (dir && !dir->winflg(2+2*lr)) {
				tree_force(lr, dynamic_cast<ADirNode*>(dir->up())); tree_expand(lr, dir); }}
		void tree_sel(int lr, ANode * nd);
		void t2_sel(int lr, ANode * nd);

		void clip_box(ClipNode * cl, int i, int j = -1, int sel = -1);
		void clip_box_1(ClipNode * cl, int i);
		void clip_flg(int id, int fid, int v01);
		void node_name(int i, ANode * nd);
		void node_rm(int i, ANode * nd);
		void own_title(int flg = 3);
		void ref_title(int wwt, ANode * nd, int wwix = -1, const char * defstr = 0);
		void j_upd(int wwt, int st, int wwix = -1);
		void errq_add(int ec);
		void mcfg_ud(int wch, cfg_ent * pc, const char *sdef, int ldef);
		void mcfg_win(int flg);
		void savename();
		void vol();
		void set_tlog();
	protected:
		int m_pid, m_inpipe, m_outpipe, m_tpipe, m_gf0;
		int m_cwt, m_lwi, m_cwi;
		int m_th, m_wh;
		int m_errq_v[32], m_errq_n, m_errq_t0[2];
		char m_buf0[65536], *m_bufp;
};

extern GuiStub gui2;

#endif // __qwe_guistub_h__
