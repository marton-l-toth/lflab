#include <limits.h>

#include "util.h"
#include "expr.h"
#include "glob.h"
#include "box0.h"


#define CALC_BLKSIZ 256

struct CalcStkElem { double *p; double c; double a[CALC_BLKSIZ]; };

struct CalcEnv1 { 
	double *px, *py0, *sp0; 
	int nt, ny; 
	CalcXF * xf; 
	CalcStkOp * aux;
	double ixv[16];
};

class CalcEnvB {
	public:
		CalcEnvB(int ifg, double ** _ppx, double **_ppy, CalcStkElem * sp, int nx) :
			ppx(_ppx), ppy0(_ppy), sp0(sp), iflg(ifg) {
				BVFOR_JM(((1<<nx)-1) & ~ifg) con[j] = ppx[j][0]; }
		double **ppx, **ppy0, con[30];
		CalcStkElem * sp0;
		int nt, ny;
		int xy_offs;
		int iflg;
		CalcStkOp * aux;
};

typedef void (*c1fun_t) (CalcStkOp * op, CalcEnv1 * env);
typedef void (*cbfun_t) (CalcStkOp * op, CalcEnvB * env, int i0, int i1);
typedef int (*ccgfun_t) (LWArr<CalcStkOp>* to, LWArr<CalcStkOp>* aux, int f, double *p, int len);

struct SymTabEnt { int nm, fid; };
class  FunTabEnt {
	public:
		enum flg_t { _i=1, _j0=2, _j1=4, _x=8, _ac=16 };
		const char * sym_s;
		int sym_i;
		c1fun_t fun1;
		cbfun_t fun_b;
		ccgfun_t fun_cg;
		int flg;
		int acmin, acmax; // min!=max --> stkop.i==ac
		void dump(FILE *f);
}; 

extern FunTabEnt funtab[];
static c1fun_t * fun1_tab;
static cbfun_t * fun_b_tab;
static int ft_end = 0;

void dump_op(FILE *f, CalcStkOp * op, const char * pref = 0) {
	if (!pref) pref="";
	FunTabEnt * q = funtab + op->f;
	fprintf(f,"%s%s", pref, q->sym_s);
	if (q->flg & 1) fprintf(f,"_%d",op->i);
	if (q->flg & 2) fprintf(f,"_%d",op->xj.j[0]);
	if (q->flg & 4) fprintf(f,"_%d",op->xj.j[1]);
	if (q->flg & 8) fprintf(f,"_%g",op->xj.x);
}

void dump_ops(FILE *f, const char *pref, CalcStkOp *p, int n) {
	fprintf(f,"%s",pref);
	for (int i=0; i<n; i++) dump_op(stderr, p+i, " ");
	fflush(f);
}

static LWArr<SymTabEnt> symtab;

static double calc_stk1[4096];
static double * calc_sp1 = calc_stk1;
static CalcStkElem calc_stk_b[1024];
static CalcStkElem * calc_sp_b = calc_stk_b;

static CalcXF default_xf;

NANpun NANpun::m0_st(0,0);
char sym6_tab[64];
char sym6_r_tab[256];

static void run1_e(CalcStkOp *op, CalcEnv1* env) {
	c1fun_t fun; while ((fun = fun1_tab[op->f])) (*fun) (op++, env);
	if (op->f != ft_end) bug("run/1: null func for %s(%d)", funtab[op->f].sym_s, op->f);
}

void calc_run1(CalcStkOp *p, double *px, double *py, CalcXF * xf) {
	CalcEnv1 env; if (!xf) xf = &default_xf;
	env.px = px; env.py0 = py; env.sp0 = calc_sp1 + 1; env.xf = xf;
	run1_e(p, &env);
	if (env.ny) memcpy(py, env.sp0+env.nt, env.ny*sizeof(double));
	IFDBGX(EXPR) {
		fprintf(stderr, " (run1 nt:%d ny:%d stk:%d) ", env.nt, env.ny, (int)(calc_sp1 - env.sp0 + 1));
		fflush(stderr);
	}
	calc_sp1 = env.sp0 - 1;
}

static void run_b_e(CalcStkOp *op, CalcEnvB* env, int i0, int i1) {
	cbfun_t fun; while ((fun = fun_b_tab[op->f])) (*fun) (op++, env, i0, i1);
	if (op->f != ft_end) bug("run/b: null func for %s(%d)", funtab[op->f].sym_s, op->f);
}

int calc_run_b(CalcStkOp *op, int iflg, double **ppx, double **ppy, int n, int nx) {
	CalcEnvB env(iflg, ppx, ppy, calc_sp_b+1, nx);
	int i0 = 0, i1 = 0, r = 0;
	while (i0<n) {
		int n2 = CALC_BLKSIZ;
		if ( (i1 = i0 + n2) > n) i1 = n, n2 = n - i0;
		env.xy_offs = i0;
		run_b_e(op, &env, 0, n2);
		CalcStkElem * sy0 = env.sp0 + env.nt;
		for (int j=0,m=1; j<env.ny; j++,m+=m) {
			double x, *q, *to = ppy[j] + i0;
			if (sy0[j].p) { 
				if (!(r&m)) for(q=ppy[j]+1,x=q[-1],r|=m; q<to; q++) *q=x;
				memcpy(to, sy0[j].p, n2*sizeof(double));
			} else {
				double x = sy0[j].c;
				for (int k=0; k<n2; k++) to[k] = x;
			}
		}
		i0 = i1; calc_sp_b = env.sp0 - 1;
	}
	return r;
}

static void sym6_init()
{
	for (int i=0; i<10; i++) sym6_tab[i] = '0'+i;
	for (int i=0; i<26; i++) sym6_tab[10+i] = 'a'+i;
	memcpy(sym6_tab+36,"!#%&*+-./:;<=>?@^_{|}~ABCDEF",28);
	for (int i=0; i<256; i++) sym6_r_tab[i] = 99;
	for (int i=0; i<60; i++) {
		int k = sym6_tab[i]; sym6_r_tab[k] = i;
		if (k>='a' && k<='z') sym6_r_tab[k-32] = i;
	}
}

const char * i_to_sym6(int k) {
	static char ret[6]; ret[5]=0;
	int i = 5;
	if (k<=0) return k?"!NEG!":"!NUL!";
	if (k>0x39e79e79) return "!OVF!";
	while (k) ret[--i] = sym6_tab[k&63], k>>=6;
	return ret+i;
}

int sym6_to_i(const char * s) {
	int k, r = 0;
	while ( (k=sym6_r_tab[(unsigned char)*(s++)]) < 99 ) r = 64*r+k;
	return r;
}

static double mk_struc(LWArr<double>* to, int len, double * arg)
{
	int n = to->n();
	to -> resize(n + len);
	double * p = to->p() + n;
	memcpy(p, arg, len*sizeof(double));
	return NANpun::mk(512+len, n);
}

static double mk_f0(LWArr<double>* to, int f) {
	double x = NANpun::mk('$', f);
	return mk_struc(to, 1, &x);
}

static double mk_f1(LWArr<double>* to, int f, double x1) {
	double x[2]; x[0] = NANpun::mk('$', f);
	x[1] = x1;
	return mk_struc(to, 2, x);
}

static double mk_f2(LWArr<double>* to, int f, double x1, double x2) {
	double x[3]; x[0] = NANpun::mk('$', f);
	x[1] = x1; x[2] = x2;
	return mk_struc(to, 3, x);
}

static double mk_f0(LWArr<double>* to, const char * f) 
	{ return mk_f0(to, sym6_to_i(f)); } 
static double mk_f1(LWArr<double>* to, const char * f, double x1) 
	{ return mk_f1(to, sym6_to_i(f), x1); } 
static double mk_f2(LWArr<double>* to, const char * f, double x1, double x2) 
	{ return mk_f2(to, sym6_to_i(f), x1, x2); } 
int kindof_int(double x) {
	if (fabs(x)<1e-300) return 0;
	if (fabs(x-round(x))<(x*1e-13)) return (int)lround(x);
	return INT_MAX;
}

static void c1_xf0(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf0(k,calc_sp1); }
static void c1_xf1(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf1(k,calc_sp1); }
static void c1_xf2(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf2(k,calc_sp1); }
static void c1_xf3(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf3(k,calc_sp1); }
static void c1_xf4(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf4(k,calc_sp1); }
static void c1_xf5(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf5(k,calc_sp1); }
static void c1_xf6(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf6(k,calc_sp1); }
static void c1_xf7(CalcStkOp * op, CalcEnv1 * env) { int k=op->i; calc_sp1+=1-k; *calc_sp1=env->xf->xf7(k,calc_sp1); }

static void c1_pop1(CalcStkOp * op, CalcEnv1 * env) { --calc_sp1; }
static void cb_pop1(CalcStkOp * op, CalcEnvB * env, int i0, int i1) { --calc_sp_b; }

static void c1_env(CalcStkOp * op, CalcEnv1 * env) { env->ny = op->i; env->nt = op->xj.j[0]; }
static void cb_env(CalcStkOp * op, CalcEnvB * env, int i0, int i1) { env->ny = op->i; env->nt = op->xj.j[0]; }

static void c1_aux(CalcStkOp * op, CalcEnv1 * env) { env->aux = op + op->i; }
static void cb_aux(CalcStkOp * op, CalcEnvB * env, int i0, int i1) { env->aux = op + op->i; }

static void c1_con(CalcStkOp * op, CalcEnv1 * env) { *(++calc_sp1) = op->xj.x; }
static void cb_con(CalcStkOp * op, CalcEnvB * env, int i0, int i1) {
	CalcStkElem *p = ++calc_sp_b;
	p->p = 0; p->c = op->xj.x;
}

static void c1_ldx(CalcStkOp * op, CalcEnv1 * env) { *(++calc_sp1) = env->px[op->i]; }
static void cb_ldx(CalcStkOp * op, CalcEnvB * env, int i0, int i1) {
	CalcStkElem *p = ++calc_sp_b;
	if (env->iflg & (1<<op->i)) p->p = env->ppx[op->i] + env->xy_offs;
	else p->p=0, p->c = env->con[op->i];
}

static void c1_ldz(CalcStkOp * op, CalcEnv1 * env) { *(++calc_sp1) = env->sp0[op->i]; }
static void c1_ldy(CalcStkOp * op, CalcEnv1 * env) { *(++calc_sp1) = env->sp0[op->i + env->nt]; }
static void cb_ldz(CalcStkOp * op, CalcEnvB * env, int i0, int i1) {
	CalcStkElem *p=++calc_sp_b, *q=env->sp0+op->i; p->p=q->p, p->c=q->c; }
static void cb_ldy(CalcStkOp * op, CalcEnvB * env, int i0, int i1) {
	CalcStkElem *p=++calc_sp_b, *q=env->sp0+op->i+env->nt; p->p=q->p, p->c=q->c; }

static void c1_ldi(CalcStkOp * op, CalcEnv1 * env) { *(++calc_sp1) = env->ixv[op->i]; }

static void c1_neg(CalcStkOp * op, CalcEnv1 * env) { *calc_sp1 = - *calc_sp1; }
static void cb_neg(CalcStkOp * op, CalcEnvB * env, int i0, int i1) {
	double *p = calc_sp_b->p;
	if (!p) { calc_sp_b->c = -calc_sp_b->c; return; }
	double *q = calc_sp_b->a;
	for (int i=i0; i<i1; i++) q[i] = -p[i];
}

cx_op2(add,+,sy->c,+)
cx_op2(sub,-,sy->c,-)
cx_op2(mul,*,sy->c,*)
cx_op2(div,/,1.0/sy->c,*)

static void c1_pow(CalcStkOp * op, CalcEnv1 * env)
{
	double * p = --calc_sp1;
	*p = pow(*p, p[1]);
}

static void cb_pow(CalcStkOp * op, CalcEnvB * env, int i0, int i1)
{
	CalcStkElem *sy = calc_sp_b, *sx = calc_sp_b-1;
	double *px = sx->p, *py = sy->p, *to = sx->a;
	-- calc_sp_b;
	if (!px) {
		if (!py) {
			sx->p=0; sx->c = pow(sx->c,sy->c);
		} else {
			double lx = log(sx->c);
			for (int i=i0; i<i1; i++) to[i] = exp(lx * py[i]);
			sx->p = to;
		}
	} else {
		if (!sy->p) {
			double y = sy->c;
			int k = kindof_int(y);
			if (!k) {
				for (int i=i0; i<i1; i++) to[i] = 1.0;
			} else if (k==INT_MAX) {
				for (int i=i0; i<i1; i++) to[i] = pow(px[i], y);
			} else if (k>0) {
				for (int i=i0; i<i1; i++) to[i] = ipow(px[i], k);
			} else {
				for (int i=i0; i<i1; i++) to[i] = 1.0 / ipow(px[i], -k);
			}
		} else {
			for (int i=i0; i<i1; i++) to[i] = pow(px[i], py[i]);
		}
		sx->p = to;
	}
}

static void c1_if(CalcStkOp * op, CalcEnv1 * env) {
	int sub = calc_sp1[-1] < calc_sp1[0] ? op->i : op->xj.j[0];
	calc_sp1 -= 2;
	run1_e(env->aux + sub, env);
}

static void cb_if(CalcStkOp * op, CalcEnvB * env, int i0, int i1) {
	int n = i1 - i0; if (!n) return;
	CalcStkElem * trg = --calc_sp_b;
	int iflg = !!trg[0].p + 2*!!trg[1].p;
	int rl_n = 0, rl_last = i0, rlen[n];
	int cc, c0;
	double x;
	switch(iflg) {
		case 0: c0 = (trg[0].c < trg[1].c); break;
		case 1: x = trg[1].c; cc = c0 = (trg[0].p[i0] < x);
			for (int i=i0+1; i<i1; i++) if (cc != (trg[0].p[i] < x)) 
				cc ^= 1, rlen[rl_n++] = i-rl_last, rl_last=i;
			break;
		case 2: x = trg[0].c; cc = c0 = (x < trg[1].p[i0]);
			for (int i=i0+1; i<i1; i++) if (cc != (x < trg[1].p[i])) 
				cc ^= 1, rlen[rl_n++] = i-rl_last, rl_last=i;
			break;
		case 3: cc = c0 = (trg[0].p[i0] < trg[1].p[i0]);
			for (int i=i0+1; i<i1; i++) if (cc != (trg[0].p[i] < trg[1].p[i])) 
				cc ^= 1, rlen[rl_n++] = i-rl_last, rl_last=i;
			break;
	}
	rlen[rl_n++] = i1-rl_last;
	CalcStkOp * sub[2]; sub[1] = env->aux + op->i; sub[0] = env->aux + op->xj.j[0];
	for (int i=0; i<rl_n; c0^=1, i0+=rlen[i++]) {
		--calc_sp_b; run_b_e(sub[c0], env, i0, i0+rlen[i]);
		if (trg->p==trg->a) continue;
		if (!trg->p) {
			x=trg->c; trg->p = trg->a;
			for (int j=0; j<rlen[i]; j++) trg->a[i0+j] = x;
		} else {
			memcpy(trg->a+i0, trg->p+i0, rlen[i]*sizeof(double));
			trg->p = trg->a;
		}
	}
	if (i0 != i1) bug("cb_if: i0(%d) != i1(%d)", i0, i1);
	if (calc_sp_b!=trg) bug("cb_if: sp-trg=%d",calc_sp_b-trg);
}

static void c1_for01(CalcStkOp * op, CalcEnv1 * env) {
	CalcStkOp * sub = env->aux + op->i;
	int ixix = op->xj.j[0];
	int n = calc_sp1[0];
	if (!n) return;
	double x0 = n>1 ? 0.0 : 0.5;
	double xd = 1.0 / (double)(n-1);
	for (int i=0; i<n; i++, x0+=xd) {
		env->ixv[ixix] = x0;
		--calc_sp1; run1_e(sub, env);
	}
}

double get_sym(const char **pp, int len=5) {
	int k, r = 0;
	static int s_pi = -1; if (s_pi<0) s_pi = sym6_to_i("pi");
	while ( (k=sym6_r_tab[(unsigned char)*((*pp)++)])<99 && (len--) ) r = 64*r+k;
	-- (*pp);
	if (r==s_pi) return M_PI;
	return r ? NANpun::mk('$', r) : NANpun::mk('E');
}

int dump_token(FILE *f, double x) {
	if (x==x) return fprintf(f, "%g", x);
	NANpun pun(x);
	int ty = pun.x19();
	int v = pun.x32();
	if (ty=='$') return fprintf(f,"$%s",i_to_sym6(v));
	if (!v && ty!='i' && (ty<'x'||ty>'z')) return fputc(ty, f) != EOF;
	if (ty>511) return fprintf(f,"[%d:%d]",v,ty-512);
	return fprintf(f,"%c%d",ty,pun.x32());
}

static void expr_relativize(double *p, int n) {
	for (int i=1; i<n; i++) {
		if (p[i]==p[i]) continue;
		NANpun pun(p+i);
		if (pun.x19()<512) continue;
		pun.set32(pun.x32()-i);
		p[i] = pun.x();
	}
}
static void expr_relativize(LWArr<double> *a) { expr_relativize(a->p(), a->n()); }

static void expr_dump(FILE *f, double *p) {
	if (!p) { fprintf(f,"[NULL]"); return; }
	if (*p==*p) { dump_token(f,*p); return; }
	NANpun pun(p); int k = pun.x19() - 512;
	if (k<0) { dump_token(f,*p); return; }
	p += pun.x32();
	for (int i=0; i<k; i++) {
		fputc(i?' ':'(', f);
		expr_dump(f,p+i);
	}
	fputc(')',f);
}

static double get_num(const char** pp, bool frac) {
        char buf[21];
        bool e_flg = frac;
        bool dot_flg = frac;
        bool sg_flg = false;
        int i = 0;
        for ( ;; buf[i] = **pp, i++, (*pp)++ ) {
                if (i==20) return NANpun::mk('E');
                int c = **pp;
                if (sg_flg) {
                        sg_flg = false;
                        if (c=='+' || c=='-') continue;
                }
                if (c>='0' && c<='9') continue;
                if (c=='E'||c=='e') {
                        if (!e_flg) break;
                        e_flg = dot_flg = false;
                        sg_flg = true;
                        continue;
                }
                if (c=='.' && dot_flg) {
                        dot_flg = false;
                        continue;
                }
                break;
        }
        buf[i] = 0;
        return frac ? atof(buf) : atoi(buf);
}

static int get_int(const char ** pp) { return (int)lround(get_num(pp, false)); }

static double get_token(const char** pp)
{
        while (**pp==' ' || **pp=='\t' || **pp=='\n')
                ++(*pp);
        int c0 = **pp;
        if (c0=='.' || (c0>='0'&&c0<='9'))
                return get_num(pp, true);
        ++ (*pp);
        switch(c0) {
                case 0:
                        return NANpun::mk('e');
                case 'x':
                case 'y':
                case 'z':
		case 'i':
			if ('0'<=**pp && **pp<='9') return NANpun::mk(c0, get_int(pp));
			--(*pp); return get_sym(pp); 
                case '+':
                case '-':
                case '*':
                case '/':
                case '^':
                case '(':
                case ',':
                case ')':
                        return NANpun::mk(c0);
                default:
                        --(*pp); return get_sym(pp);
        }
}

static bool is_nantok(double x, int i) {
	if (x==x) return false;
	NANpun pun(x); return pun.x19()==i;
}
static bool is_nantok_2(double x, int i, int j) {
	if (x==x) return false;
	NANpun pun(x); return pun.x19()==i||pun.x19()==j;
}

#define PREC_MINUS 0
static const char * calc_2ops[] = { "+-", "*/", 0 };
static const char * calc_op2fun[] = { "add", "sub", "mul", "div" };
static double get_op2exp(LWArr<double> * tree, double** pptok, int prec);

static double get_func(LWArr<double> * tree, double** pptok, double head)
{
	LWArr<double> av; av.add(head);
        if (!is_nantok(**pptok, '('))
                return log("calc/func:'(' expected"), NANpun::mk('E');
        ++(*pptok);

        while(1) {
                double x = **pptok;
		if (is_nantok(x,')')) { ++(*pptok); break; }
		if (av.n()>1) {
			if (!is_nantok(x,',')) 
				return log("calc/func: ',' or ')' expexted"), NANpun::mk('E');
			++(*pptok);
		}
                x = get_op2exp(tree, pptok, 0);
                if (is_nantok(x,'E')) return x; else av.add(x);
        }
	return mk_struc(tree, av.n(), av.p());
}

static double get_atom(LWArr<double> * tree, double** pptok)
{
        double x = **pptok;
        if (x==x) return ++(*pptok), x;
	NANpun pun(x); int ty = pun.x19();
        switch (ty) {
                case '(':
                        ++(*pptok);
                        x = get_op2exp(tree, pptok, 0);
                        if (is_nantok(x,'E')) return x;
                        if (!is_nantok(**pptok,')')) {
                                log("calc/atom: closing ')' not found");
                                return NANpun::mk('E');
                        }
                        return ++(*pptok), x;
		case 'x': case 'y': case 'z': case 'i':
                        return ++(*pptok), x;
                default:
                        break;
        }
        if (ty!='$') return log("calc/atom: expected num|var|func|'('"), NANpun::mk('E');
        return ++(*pptok), get_func(tree, pptok, x);
}

static double get_powexp(LWArr<double> * tree, double ** pptok)
{
        double base = get_atom(tree, pptok);
        if (is_nantok(base,'E') || !is_nantok(**pptok,'^')) return base;
        ++(*pptok);
        double up = get_powexp(tree, pptok);
        if (base==base && up==up) {
                double x = pow(base, up);
                return (x==x) ? x : NANpun::mk('E');
        }
        if (is_nantok(up,'E')) return up;
	return mk_f2(tree, "pow", base, up);
}

static double get_op2exp(LWArr<double> * tree, double** pptok, int prec)
{
        const char * op = calc_2ops[prec];
        double cur;
        if (prec==PREC_MINUS && is_nantok(**pptok, '-')) {
                ++(*pptok);
                cur = get_op2exp(tree, pptok, prec + 1);
                if (cur==cur) cur = -cur;
                else if (is_nantok(cur, 'E')) return cur;
                else cur = mk_f1(tree, "neg", cur);
        } else {
                if (!op) return get_powexp(tree, pptok);
                cur = get_op2exp(tree, pptok, prec + 1);
                if (is_nantok(cur,'E')) return cur;
        }
        while(1) {
                double x = **pptok;
                if (x==x) return cur;
		NANpun pun(x); int k = pun.x19(); if (k>127) return cur;
                int op_ix = strchr_i(op, k); if (op_ix<0) return cur;
                ++ (*pptok);
                double arg2 = get_op2exp(tree, pptok, prec + 1);
                if (is_nantok(arg2,'E')) return arg2;
		cur = mk_f2(tree, calc_op2fun[2*prec+op_ix], cur, arg2);
        }
}

static bool get_exp(LWArr<double> * tree, double** pptok)
{
        tree -> resize(0);
        tree -> add(0.0);
        double x = get_op2exp(tree, pptok, 0);
        if (is_nantok(x,'E')) {
                tree->clear();
                return false;
        }
        tree->p()[0] = x;
	expr_relativize(tree);
        return true;
}

static int symtab_lookup(int ni) {
	if (ni<=0) return -1;
	SymTabEnt * p = symtab.p();
	int lo = 0, hi = symtab.n()-1;
	while(lo<hi) {
		int mid = (lo+hi)>>1;
		int k = p[mid].nm;
		if (ni==k) return p[mid].fid;
		if (ni<k) hi = mid - 1;
		else lo = mid + 1;
	}
	return (lo==hi && p[lo].nm==ni) ? p[lo].fid : -1;
}

#define FT_CON  0
#define FT_LDX  1 // LDY=LDX+1, LDZ=LDX+2
#define FT_POP1 4
#define FT_ENV  5
#define FT_IF   6
#define FT_AUX  7
#define FT_LDI  8

static int exp_codegen_2(LWArr<CalcStkOp>* i_to, LWArr<CalcStkOp>* i_aux, double * ex);
static int exp_codegen_aux(LWArr<CalcStkOp>* i_aux, double * ex);
static int cgf_def(LWArr<CalcStkOp>* i_to, LWArr<CalcStkOp>* i_aux, int f, double *ps, int k) {
	int r, ac0 = funtab[f].acmin, ac1 = funtab[f].acmax;
	if (k-1<ac0) return '<'; if (k-1>ac1) return '>';
        for (int i=1; i<k; i++) 
		if ((r = exp_codegen_2(i_to, i_aux, ps+i))) return r;
	CalcStkOp * p = i_to->add();
	p->f = f;
	if (ac0!=ac1 || (funtab[f].flg & FunTabEnt::_ac)) p->i = k-1;
	return 0;
}

static int ccg_if(LWArr<CalcStkOp>* i_to, LWArr<CalcStkOp>* i_aux, int f, double *ps, int k)
{
	if (k!=5) return "<>"[k>5];
	int r, a1, a2;
	if ((r = exp_codegen_2(i_to, i_aux, ps+1))) return r;
	if ((r = exp_codegen_2(i_to, i_aux, ps+2))) return r;
	if ((a1 = exp_codegen_aux(i_aux, ps+3)) < 0) return -a1;
	if ((a2 = exp_codegen_aux(i_aux, ps+4)) < 0) return -a2;
	CalcStkOp * p = i_to->add();
	p->f = FT_IF; p->i = a1; p->xj.j[0] = a2;
	return 0;
}

static int ccg_for(LWArr<CalcStkOp>* i_to, LWArr<CalcStkOp>* i_aux, int f, double *ps, int k)
{
	if (k!=4) return "<>"[k>4];
	int r, aux;
	if (ps[1]==ps[1]) return 'a';
	NANpun ixn(ps+1); if (ixn.x19()!='i') return 'a';
	int ixix = ixn.x32();
	if ((r = exp_codegen_2(i_to, i_aux, ps+2))) return r;
	if ((aux = exp_codegen_aux(i_aux, ps+3)) < 0) return -aux;
	CalcStkOp * p = i_to->add();
	p->f = f; p->i = aux; p->xj.j[0] = ixix;
	return 0;
}

static int exp_codegen_2(LWArr<CalcStkOp>* i_to, LWArr<CalcStkOp>* i_aux, double * ex)
{
	CalcStkOp * p;
	if (*ex==*ex) {
		p = i_to->add();
		p->f = FT_CON; p->xj.x = *ex;
		return 0;
	}
	NANpun pun(ex); int i, k = pun.x19();
	switch(k) {
		case 'x':
		case 'y':
		case 'z':
			p = i_to->add();
			p->f = FT_LDX + k - 'x';
			p->i = pun.x32();
			return 0;
		case 'i':
			i = pun.x32();
			if (i<0||i>15) return 'i';
			p = i_to->add();
			p->f = FT_LDI; p->i = i;
			return 0;
		default:
			break;
	}
	k -= 512; if (k<=0) return k?'e':'E';
	double * ps = ex + pun.x32();
	pun.set(ps);
	if (pun.x19() != '$') return 's';
	int si = pun.x32();
	int f = symtab_lookup(si);
	if (f<0) return 'f';
	ccgfun_t cg = funtab[f].fun_cg;
	return cg ? (*cg)   (i_to, i_aux, f, ps, k) 
		  : cgf_def (i_to, i_aux, f, ps, k);
}

static int exp_codegen_aux(LWArr<CalcStkOp>* i_aux, double * ex) {
	LWArr<CalcStkOp> i_tmp;
	int k = exp_codegen_2(&i_tmp, i_aux, ex);
	if (k) return -k;
	int i0 = i_aux->n(), n = i_tmp.n();
	i_aux->resize(i0+n); memcpy(i_aux->p(i0), i_tmp.p(), n*sizeof(CalcStkOp));
	i_aux -> add() -> f = ft_end;
	return i0;
}

static int exp_codegen(LWArr<CalcStkOp>* i_main, LWArr<CalcStkOp>* i_aux, double *ex) {
	return exp_codegen_2(i_main, i_aux, ex);
}

//                   _i=1, _j0=2, _j1=4, _x=8, _ac=16 \/
FunTabEnt funtab[] = {
	{"con",     0, c1_con,    cb_con,    0,        8,  0, -1 }, //0
	{"ldx",     0, c1_ldx,    cb_ldx,    0,        1,  0, -1 }, 
	{"ldy",     0, c1_ldy,    cb_ldy,    0,        1,  0, -1 },
	{"ldz",     0, c1_ldz,    cb_ldz,    0,        1,  0, -1 },
	{"pop1",    0, c1_pop1,   cb_pop1,   0,        0,  0, -1 }, //4
	{"env",     0, c1_env,    cb_env,    0,        3,  0, -1 },    
	{"if<",	    0, c1_if,     cb_if,     ccg_if,   3,  4, 4  },
	{"aux",	    0, c1_aux,    cb_aux,    0,        1,  1, 1  },
	{"ldi",     0, c1_ldi,    0,         0,        1,  0, -1 }, //8
	{"add",	    0, c1_add,    cb_add,    0,        0,  2, 2  }, 
	{"sub",	    0, c1_sub,    cb_sub,    0,        0,  2, 2  },    
	{"mul",	    0, c1_mul,    cb_mul,    0,        0,  2, 2  },
	{"div",	    0, c1_div,    cb_div,    0,        0,  2, 2  }, //12
	{"pow",	    0, c1_pow,    cb_pow,    0,        0,  2, 2  }, 
	{"neg",	    0, c1_neg,    cb_neg,    0,        0,  1, 1  },
	{"xf0",	    0, c1_xf0,    0,         0,        0,  0, 15 },     
	{"xf1",	    0, c1_xf1,    0,         0,        0,  0, 15 }, //16    
	{"xf2",	    0, c1_xf2,    0,         0,        0,  0, 15 },     
	{"xf3",	    0, c1_xf3,    0,         0,        0,  0, 15 },     
	{"xf4",	    0, c1_xf4,    0,         0,        0,  0, 15 },     
	{"xf5",	    0, c1_xf5,    0,         0,        0,  0, 15 },     
	{"xf6",	    0, c1_xf6,    0,         0,        0,  0, 15 },     
	{"xf7",	    0, c1_xf7,    0,         0,        0,  0, 15 },     
	{"end",     0, 0,         0,         0,        0,  0, -1 } 
};

static int ste_cmp(const void * p, const void * q) {
	return ((const SymTabEnt*)p)->nm - ((const SymTabEnt*)q)->nm; }

static void funtab_init() {
	FunTabEnt *p;
	for (p = funtab; 1; p++) {
		p->sym_i = sym6_to_i(p->sym_s);
		if (!p->fun1 && !p->fun_b && !p->fun_cg) break;
		if (p->acmax<0) continue;
		SymTabEnt * q = symtab.add();
		q->nm = p->sym_i; q->fid = p-funtab;
	}
	ft_end = p - funtab;
	fun1_tab  = new c1fun_t[ft_end+1];
	fun_b_tab = new cbfun_t[ft_end+1];
	for (int i=0; i<=ft_end; i++) fun1_tab[i]=funtab[i].fun1, fun_b_tab[i]=funtab[i].fun_b;
	qsort(symtab.p(), symtab.n(), sizeof(SymTabEnt), ste_cmp);
}

void calc_init() { sym6_init(); funtab_init(); }

double CalcXF::debug(const char* txt, int ac, double * av) {
	fprintf(stderr," [cXF:%s", txt?txt:"(null)");
	for (int i=0; i<ac; i++) fprintf(stderr," %g",av[i]);
	fprintf(stderr,"] "); fflush(stderr);
	return 0.0;
}

void CalcExpr::vconf(int xn, int yn, int zn, bool rel) {
	if (rel) { m_xn+=xn, m_yn+=yn, m_zn+=zn; } else {
		xn>=0 && (m_xn=xn);  yn>=0 && (m_yn=yn); zn>=0 && (m_zn=zn);
	}
	chk_xyz();
}

void CalcExpr::chk_xyz() {
	if (m_ec && (m_ec<'x'||m_ec>'z')) return;
	if (m_xmax >= m_xn) m_ec='x';
	else if (m_ymax >= m_yn) m_ec='y';
	else if (m_zmax >= m_zn) m_ec='z';
	else m_ec = 0;
}

void CalcExpr::upd_xyzmax() {
	int mx[3]; mx[0]=mx[1]=mx[2]=-1;
	int n = m_tree.n();
	for (int i=0; i<n; i++) {
		double v = m_tree[i];
		if (v==v) continue;
		NANpun pun(v); int k=pun.x19();
		if (k<'x' || k>'z') continue;
		int j = pun.x32(), *p = mx + k-'x';
		if (j > *p) *p = j;
	}
	m_xmax = mx[0]; m_ymax = mx[1]; m_zmax = mx[2];
}

void CalcExpr::set(const char * s)
{
	int l = strlen(s);
	if (l && s[l-1] == '\n') --l;
	m_exp.resize(l+1); if (l) memcpy(m_exp.p(), s, l); m_exp[l] = 0;
	LWArr<double> tok; double x;
	do tok.add(x=get_token(&s)); while(!is_nantok_2(x,'e','E'));
	if (is_nantok(x,'E')) { m_ec='?'; return; }
	m_tree.clear();
	double * tp = tok.p(); 
	if (!get_exp(&m_tree, &tp)) { m_ec='('; return; }
	upd_xyzmax();
	m_i_main.resize(0); m_i_aux.resize(0);
	int r = exp_codegen(&m_i_main, &m_i_aux, m_tree.p());
	if (r) { m_ec = r; return; }
	m_ec = 0; chk_xyz();
	return;
}

int CalcExpr::add_main(CalcStkOp* q, int ixdif) {
	if (m_ec) return m_pop ? 0 : (q->f = FT_CON, q->xj.x = 0.0, 1);
	int r = m_i_aux.n() ? (q->f=FT_AUX, q->i=ixdif, 1) : 0;
	int nf = m_i_main.n(); memcpy(q+r, m_i_main.p(), nf*sizeof(CalcStkOp)); r+=nf;
	return m_pop ? (q[r].f=FT_POP1, r+1) : r;
}

int CalcExpr::add_aux(CalcStkOp* q) {
	int n = m_i_aux.n(); return (!n|m_ec) ? 0 : (memcpy(q, m_i_aux.p(), n*sizeof(CalcStkOp)), n); }

void CalcEL::code(CalcStkOp * to) {
	int nm = len_m(), na = len_a();
	to->f = FT_ENV; to->i = m_ny; to->xj.j[0] = m_nz;
	int n = m_el.n(), ai = nm, mi = 1;
	for (int i=0; i<n; i++) mi += m_el[i]->add_main(to+mi,ai-mi),  ai += m_el[i]->add_aux(to+ai);
	to[mi++].f = ft_end;
	if (mi!=nm) bug("cEL/code: main(%d) != exp_main(%d)", mi, nm);
	if (ai-nm!=na) bug("cEL/code: aux (%d) != exp_aux (%d)", ai-nm, na);
}

void CalcExpr::dump(FILE *f)
{
	fprintf(f,"ec:%c xyz_max:%d %d %d", m_ec?m_ec:'0', m_xmax, m_ymax, m_zmax);
	fprintf(f," xyz_n:%d %d %d",m_xn, m_yn, m_zn);
	fprintf(f," xf_f:%x pop:%c",m_xf_flg, "ny"[m_pop]);
	fprintf(f," str:\"%s\" exp:",m_exp.n()?m_exp.p():"");
	expr_dump(f,m_tree.p());
	dump_ops(f, " main:",m_i_main.p(), m_i_main.n());
	if (m_i_aux.n()) dump_ops(f," aux:",m_i_aux.p(), m_i_aux.n());
	fputc('\n',f); fflush(f);
}

CalcExpr * CalcEL::ins(int ty, int i)
{
	int * pn = pnx(ty); if (!pn) return 0; 
	if (i<0) i=0; else if (i>*pn) i=*pn;
	CalcExpr * p = new CalcExpr();
	int ix = eix(ty, i); 
	m_el.ins(ix, p);
	++ *pn;
	p -> vconf(m_nx, ty=='y'?i:0, ty=='z'?i:m_nz);
	p -> set("0");
	if (ty=='p') { p->set_pop(1); return p; }
	if (ty=='z') for (int i=ix+1; i<m_el.n(); i++) m_el[i]->vconf(0,0,1,1);
	else for (int i=ix+1; i<m_el.n(); i++) m_el[i]->vconf(0,1,0,1);
	return p;
}

bool CalcEL::cut(int ty, int i)
{
	int * pn = pnx(ty); if (!pn) return 0; 
	if (i<0 || i>=*pn) return 0;
	int ix = eix(ty, i); 
	m_el.cut(ix);
	-- *pn;
	if (ty=='p') return 1;
	if (ty=='z') for (int i=ix; i<m_el.n(); i++) m_el[i]->vconf(0,0,-1,1);
	else for (int i=ix; i<m_el.n(); i++) m_el[i]->vconf(0,-1,0,1);
	return 1;
}

bool CalcEL::set_l(int ty, int n) {
	int *pn = pnx(ty); if (!pn || n<0) return 0;
	if (n==*pn) return 0;
	if (n<*pn) while(cut(ty, *pn - 1) && *pn>n);
	else while(ins(ty, *pn) && *pn<n);
	return 1;
}

void CalcEL::dump(FILE *f)
{
	fprintf(f, "cEL:m_nz, m_np, m_ny: %d %d %d\n", m_nz, m_np, m_ny);
	int n = m_el.n();
	if (m_nz + m_np + m_ny != n) fprintf(f,"BUG:nz+np+ny!=%d\n",n);
	for (int i=0; i<m_el.n(); i++) {
		int j=i, c='z';
		if (j>=m_nz) (j-=m_nz)<m_np ? (c='p') : (j-=m_np, c='y');
		fprintf(f,"%c%d: ",c,j); m_el[i]->dump(f);
	}
	int nc = code_len() / sizeof(CalcStkOp);
	if (nc) { CalcStkOp cod[nc]; code(cod); dump_ops(f,"code:", cod, nc); fputc('\n',f); }
	else { fprintf(f,"no code\n"); }
}
