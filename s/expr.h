#ifndef __qwe_expr_h__
#define __qwe_expr_h__

void calc_init();

class NANpun {
	public:
		NANpun(int x19, int x32) { set(x19, x32); }
		NANpun(double x) { set(x); }
		NANpun(const double *p) { set(p); }
		static double mk(int x19, int x32 = 0) { m0_st.set(x19,x32); return m0_st.x(); }
		void set(double x) { memcpy(m_i,&x,8); }
		void set(const double *p) { memcpy(m_i,p,8); }
		void set(int x19, int x32) { set19(x19); set32(x32); }
		void set32(int x32) { m_i[0] = x32; }
		void set19(int x19) { m_i[1] = 0x7ff80000 | (0x7ffff & x19); }
		int x19() const { return m_i[1] & 0x7ffff; }
		int x32() const { return m_i[0]; }
		double x() const { double r; memcpy(&r, m_i, 8); return r; }
		void dump(FILE *f) { fprintf(f,"(%g %d %d)",x(),x19(),x32()); }
	protected:
		static NANpun m0_st;
		int m_i[2];
};

struct CalcStkOp { int f, i; union {double x; int j[2];} xj; };

class CalcXF {
	public:
		virtual double xf0(int ac, double * av) { return debug("xf0", ac, av); }
		virtual double xf1(int ac, double * av) { return debug("xf1", ac, av); }
		virtual double xf2(int ac, double * av) { return debug("xf2", ac, av); }
		virtual double xf3(int ac, double * av) { return debug("xf3", ac, av); }
		virtual double xf4(int ac, double * av) { return debug("xf4", ac, av); }
		virtual double xf5(int ac, double * av) { return debug("xf5", ac, av); }
		virtual double xf6(int ac, double * av) { return debug("xf6", ac, av); }
		virtual double xf7(int ac, double * av) { return debug("xf7", ac, av); }
		virtual ~CalcXF() {}
	protected:
		double debug(const char* txt, int ac, double * av);
};

void calc_run1(CalcStkOp *p, double *px, double *py, CalcXF * xf = 0);
int calc_run_b(CalcStkOp *p, int iflg, double **ppx, double **py, int n, int nx);

class CalcExpr {
	public:
		enum err_t { _ok=0, _syn=-1, _chk=-2 };
		CalcExpr() :  m_xn(0), m_yn(0), m_zn(0), m_xf_flg(0), m_ec('?'), m_xf(0), m_pop(0) {}
		void vconf(int xn, int yn, int zn, bool rel = false);
		void set_xf(CalcXF * xf) { m_xf = xf; }
		void set(const char * s);
		int n_main() { return m_ec ? !m_pop : m_i_main.n() + m_pop + !!m_i_aux.n(); }
		int n_aux()  { return m_ec ? 0 : m_i_aux.n(); }
		int add_main(CalcStkOp* q, int ixdif);
		int add_aux (CalcStkOp* q);
		int xmax() { return m_xmax; }
		int ymax() { return m_ymax; }
		int zmax() { return m_zmax; }
		int xf_flg() { return m_xf_flg; }
		void set_pop(bool x) { m_pop = x; }
		int ec() { return m_ec; }
		void upd_xyzmax();
		void chk_xyz();
		const char * str() { return m_exp.p(); }
		void dump(FILE *f);
	protected:
		int m_xmax, m_ymax, m_zmax;
		int m_xn, m_yn, m_zn;
		int m_xf_flg;
		LWArr<char> m_exp;
		LWArr<double> m_tree;
		LWArr<CalcStkOp> m_i_main;
		LWArr<CalcStkOp> m_i_aux;
		int m_ec;
		CalcXF * m_xf;
		bool m_pop;
};

class CalcEL {
	public:
		CalcEL() : m_nx(0), m_nz(0), m_np(0), m_ny(0) {}
		CalcExpr * ins(int ty, int i); // ty: z|p|y
		bool cut(int ty, int i);
		int nx(int ty) const {switch(ty){case 'x':return m_nx; case'z':return m_nz;case'p':return m_np;
			                   case'y':return m_ny;}return-1;}
		CalcExpr * ex(int i) { return (0<=i && i<m_el.n()) ? m_el[i] : 0; }
		CalcExpr * ex2(int ty, int i) { return ex(eix(ty,i)); }
		void set_nx(int n) { m_nx = n; for (int i=0; i<m_el.n(); i++) m_el[i]->vconf(n,-1,-1); }
		int code_len() { return sizeof(CalcStkOp) * (len_m() + len_a()); }
		void code(CalcStkOp * to);
		int xf_flg() { int r = 0; for (int i=0; i<m_el.n(); i++) r|=m_el[i]->xf_flg(); return r; }
		void dump(FILE *f);
		bool set_l(int ty, int n);
	protected:
		int len_m() { int r=2; for (int i=0; i<m_el.n(); i++) r+=m_el[i]->n_main(); return r; }
		int len_a() { int r=0; for (int i=0; i<m_el.n(); i++) r+=m_el[i]->n_aux(); return r; }
		int eix(int ty, int i) const { switch(ty) {case 'z':return i; case 'p':return m_nz+i;
		 	                             case 'y':return m_nz+m_np+i; default:return -1; } }
		int * pnx(int ty) {switch(ty){case'z':return &m_nz;case'p':return &m_np;case'y':return &m_ny;}return 0;}
		LWArr<CalcExpr*> m_el;
		int m_nx, m_nz, m_np, m_ny;
};

#define cx_fun1(fun) \
static void c1_##fun(CalcStkOp * op, CalcEnv1 * env) { *calc_sp1 = (fun)(*calc_sp1)

#define cx_op2(name,o,yv,o2) \
static void c1_##name(CalcStkOp * op, CalcEnv1 * env){double*p=--calc_sp1;*p o##= p[1];}\
static void cb_##name(CalcStkOp * op, CalcEnvB * env, int i0, int i1) { \
        CalcStkElem *sy = calc_sp_b, *sx = calc_sp_b-1; \
        double *px = sx->p, *py = sy->p, *to = sx->a; \
        -- calc_sp_b; \
        if (!px) { \
                if (!py) { \
                        sx->p = 0; sx->c = sx->c o sy->c; \
                } else { \
                        double x = sx->c; \
                        for (int i=i0; i<i1; i++) to[i] = x o py[i]; \
                        sx -> p = to; \
                } \
        } else { \
                if (!py) { \
                        double y = yv ; \
                        for (int i=i0; i<i1; i++) to[i] = px[i] o2 y; \
                } else { \
                        for (int i=i0; i<i1; i++) to[i] = px[i] o py[i]; \
                } \
                sx->p = to; \
        } \
}

class CalcWrap
{
	public:
		int n_exp(int ty) const { return m_l.nx(ty); }
		int n_tmp() { return n_exp('z'); }
		int n_pop() { return n_exp('p'); }
		const char * exp_str(int ty, int i) { return m_l.ex2(ty,i) -> str(); }
		void set_exp(int ty, int i, const char * s) { m_l.ex2(ty,i) -> set(s); }
		int ec(int ty, int i) { return m_l.ex2(ty,i) -> ec(); }
		int ecg(int ty, int i) { int k = ec(ty,i); return k ? k : ' '; }
	protected:
		CalcEL m_l;
};

#endif // __qwe_expr_h__
