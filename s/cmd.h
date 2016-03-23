#ifndef __qwe_cmd_h__
#define __qwe_cmd_h__

#include "box.h"

class CmdTok {
	public: CmdTok() : m_s0(0), m_p(0) {}
		char * t(char *s = 0);
		char * un();
		void rs() { m_s0 = 0; }
	protected: char *m_s0, *m_p;
};

class CmdBuf : public AOBuf {
	public:
		friend class CmdTab;
		static void st_init();
		static int read_batch(const char * name, int nof);
		typedef int (*conv_fun) (CmdBuf*, char*, int);

		CmdBuf() : m_iname(0), m_try(0), m_fd(-1), m_nof0(0), m_sv_M(v_major), m_sv_m(v_minor),
			   m_errcnt(0), m_cont(0), m_buf(0) {}
		~CmdBuf() { if (m_buf) free(m_buf); }
		virtual int vpf(const char * fmt, va_list ap);
		virtual int sn(const char * s, int n);
		void init(int fd, int nof, int px = 0, const char * inm = 0, int bs = 4096, int rs = 1024);
		void cfg(int nof, int rs = 0) { if (nof>=0) m_nof0 = nof; if (rs) m_rsiz = rs; }
		bool is_gui() { return m_prefix=='~'; }
		int read_f();

		char * tok(char * s = 0) { return m_c_tok.t(s); }
		void set_curnode(ANode * nd) { *m_nd_var = (m_c_node=nd) ? nd->id() : -1; }
		ANode * dlookup(int *p) { if (*p<0) return 0; 
			ANode * q = ANode::lookup_n_q(*p); return q->cl_id() ? q : (*p = -1, (ANode*)0); }
		ANode * curnode() { return dlookup(m_nd_var); }
		char * a0() const { return m_c_a0; }
		char * a1() const { return m_c_a1; }
		char * a2() const { return m_c_a2; }
		int  cnof() const { return m_c_nof; }
		ANode * cnode() const { return m_c_node; }
		int fd() const { return m_fd; }
		void show_error(int ec);
		int perm(ANode * nd, int f) { return (m_c_nof&NOF_FORCE) ? f : nd->perm(f); }
		int cperm(int f) { return perm(m_c_node, f); }
		ANode * lookup(const char * s) { return Node::lookup_cb(this, s); }
  		ANode * var(int i) { return i ? dlookup(m_nd_var+i) : m_c_node; }
		void setvar(int i, int id) { m_nd_var[i&7] = id; }
		int curnode_ccmd() { ANode * p = curnode(); return p ? p->ccmd(this) : GCE_WTF; }
		bool before(int i, int j) const { return m_sv_M<i || (m_sv_M==i && m_sv_m<j); }
	protected:
		static int cf_i2p(CmdBuf *p, char *s, int l),
			   cf_p2i(CmdBuf *p, char *s, int l);
		int rpl_cp(char *to, const char *s, conv_fun f);
		char * untok() { return m_c_tok.un(); }
		int chunk(int len);
		int bprep(int siz);
		int line(char * s);
		int fdok(int ec, int ty);

		CmdTok m_c_tok;
		ANode * m_c_node;
		char *m_c_a0, *m_c_a1, *m_c_a2;
		int m_c_nof;
		int m_nd_var[8];

		char * m_iname;
		int m_try, m_fd, m_nof0, m_prefix,
		    m_bsiz, m_rsiz, m_sv_M, m_sv_m,
		    m_rpos, m_cpos, m_lineno, m_errcnt;
		AReader * m_cont;

		char * m_buf;
};

#endif // __qwe_cmd_h__
