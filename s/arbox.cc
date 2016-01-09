#include "util.h"
#include "util2.h"
#include "expr.h"
#include "box.h"
#include "glob.h"
#include "cmd.h"
#include "guistub.h"

//? {{{!._a2}}}
//? stateless box performing a binary arithmetic operation

typedef double (*abxf_cc)(double,double);
typedef void (*abxf_cv)(double*,double,double*,int);
typedef void (*abxf_vv)(double*,double*,double*,int);

STATELESS_BOX_0(Ar2BoxAD) { double *o = outb[0]; switch(inflg&3) {
	case 0: return **outb = *inb[0] + *inb[1], 0;
	case 1: { double *p=inb[0],y=*inb[1]; for (int i=0;i<n;i++) o[i] = p[i]+y; return 1; }
	case 2: { double x=*inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = x+q[i]; return 1; }
	case 3: { double *p=inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = p[i]+q[i]; return 1; }}}
STATELESS_BOX_0(Ar2BoxSB) { double *o = outb[0]; switch(inflg&3) {
	case 0: return **outb = *inb[0] - *inb[1], 0;
	case 1: { double *p=inb[0],y=*inb[1]; for (int i=0;i<n;i++) o[i] = p[i]-y; return 1; }
	case 2: { double x=*inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = x-q[i]; return 1; }
	case 3: { double *p=inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = p[i]-q[i]; return 1; }}}
STATELESS_BOX_0(Ar2BoxML) { double *o = outb[0]; switch(inflg&3) {
	case 0: return **outb = *inb[0] * *inb[1], 0;
	case 1: { double *p=inb[0],y=*inb[1]; if (fabs(y)<1e-280) return *o=0.0, 0; 
					      for (int i=0;i<n;i++) o[i] = p[i]*y; return 1; }
	case 2: { double x=*inb[0],*q=inb[1]; if (fabs(x)<1e-280) return *o=0.0, 0; 
					      for (int i=0;i<n;i++) o[i] = x*q[i]; return 1; }
	case 3: { double *p=inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = p[i]*q[i]; return 1; }}}
STATELESS_BOX_0(Ar2BoxDV) { double *o = outb[0]; switch(inflg&3) {
	case 0: return **outb = *inb[0] / *inb[1], 0;
	case 1: { double *p=inb[0],y=*inb[1]; if (fabs(y)<1e-280) return *o=1.0/0.0, 0; else y = 1.0/y;
					      for (int i=0;i<n;i++) o[i] = p[i]*y; return 1; }
	case 2: { double x=*inb[0],*q=inb[1]; if (fabs(x)<1e-280) return *o=0.0, 0; 
					      for (int i=0;i<n;i++) o[i] = x/q[i]; return 1; }
	case 3: { double *p=inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = p[i]/q[i]; return 1; }}}

//? {{{!._xmx}}}
//? cross-mixer, output is (1-t)*x + t*y
#define XMX_L1(X,Y,J) x=(X); y=(Y); p = inb[J]; for (int i=0; i<n; i++) to[i] = x + y*p[i]; return 1
STATELESS_BOX_0(CrossMix) {
	double x, y, *p, *q, *r, *to = *outb;
	switch(inflg & 7) {
		case 0: return *to = *inb[1] * *inb[2] + *inb[0] * (1.0-*inb[2]), 0;
		case 1: XMX_L1(inb[1][0]*inb[2][0], 1.0-inb[2][0], 0);
		case 2: XMX_L1(inb[0][0]*(1.0-inb[2][0]), inb[2][0], 1);
		case 3: x = inb[2][0]; y = 1.0 - x; p=inb[0], q=inb[1];
			for (int i=0; i<n; i++) to[i] = y*p[i] + x*q[i];    return 1;
		case 4: XMX_L1(inb[0][0], inb[1][0] - x, 2);
		case 5: y = inb[1][0]; p = inb[0]; r = inb[2];
			for (int i=0; i<n; i++) to[i] = (1.0-r[i])*p[i] + r[i]*y;  return 1;
		case 6: x = inb[0][0]; q = inb[1]; r = inb[2];
			for (int i=0; i<n; i++) to[i] = (1.0-r[i])*x + r[i]*q[i];  return 1;
		case 7: p=inb[0], q=inb[1], r=inb[2]; 
			for (int i=0; i<n; i++) to[i] = p[i]*(1.0-r[i]) + q[i]*r[i]; return 1;
	}}

//? {{{!._smx}}}
//? mixer, output is a1*b1 + a2*b2 + ... + a<n>*b<n>
#define SMXC(J,X) case J: for (int i=0;i<n;i++) to[i] = (X);  return 1
#define SMXc(J) (x[J]*p[J][i])
#define SMXV(J) (q[J][i]*r[J][i])
STATELESS_BOX_1(SMix) {
	double c = 0.0, *to = outb[0], x[4], *p[4], *q[4], *r[4];
	int flg = 0, ni = m_arg*2, f = inflg;
	for (int j,i=0; i<ni; i+=2, f>>=2) { switch(f&3) {
		case 0: flg |= 64; c += inb[i][0] * inb[i+1][0]; break;
		case 1: j = flg&3; x[j] = inb[i+1][0]; p[j] = inb[ i ]; ++flg; break;
		case 2: j = flg&3; x[j] = inb[ i ][0]; p[j] = inb[i+1]; ++flg; break;
		case 3: j = (flg>>3)&3; q[j] = inb[i]; r[j] = inb[i+1]; flg+=8; break;
	}}
	switch(flg) {
		SMXC(0002, SMXc(0) + SMXc(1));
		SMXC(0003, SMXc(0) + SMXc(1) + SMXc(2));
		SMXC(0004, SMXc(0) + SMXc(1) + SMXc(2) + SMXc(3));
		SMXC(0011, SMXV(0) + SMXc(0));
		SMXC(0012, SMXV(0) + SMXc(0) + SMXc(1));
		SMXC(0013, SMXV(0) + SMXc(0) + SMXc(1) + SMXc(2));
		SMXC(0020, SMXV(0) + SMXV(1));
		SMXC(0021, SMXV(0) + SMXV(1) + SMXc(0));
		SMXC(0022, SMXV(0) + SMXV(1) + SMXc(0) + SMXc(1));
		SMXC(0030, SMXV(0) + SMXV(1) + SMXV(2));
		SMXC(0031, SMXV(0) + SMXV(1) + SMXV(2) + SMXc(1));
		SMXC(0040, SMXV(0) + SMXV(1) + SMXV(2) + SMXV(3));
		case 0100: return *to=c, 0;
		SMXC(0101, c + SMXc(0));
		SMXC(0102, c + SMXc(0) + SMXc(1));
		SMXC(0103, c + SMXc(0) + SMXc(1) + SMXc(2));
		SMXC(0110, c + SMXV(0));
		SMXC(0111, c + SMXV(0) + SMXc(0));
		SMXC(0112, c + SMXV(0) + SMXc(0) + SMXc(1));
		SMXC(0120, c + SMXV(0) + SMXV(1));
		SMXC(0121, c + SMXV(0) + SMXV(1) + SMXc(0));
		SMXC(0130, c + SMXV(0) + SMXV(1) + SMXV(2));
		default: return RTE_BUG;
	}}

extern double fq_warp(double fq);
//? {{{!._a1}}}
//? stateless box calculating a single-variable function
FUN1_BOX(Ar1Abs, fabs(x))
FUN1_BOX(Ar1Sqrt, sqrt(x))
FUN1_BOX(Ar1Cbrt, cbrt(x))
FUN1_BOX(Ar1FqWarp, fq_warp(x))

class CalcBoxGen : public BoxGen {
	public: 
		CalcBoxGen(ABoxNode * nd) : BoxGen(nd) { m_l.ins('y',0); }
		virtual int cmd(CmdBuf* cbf);
		virtual void box_window();
		virtual const char * cl_name() { return "calc"; }
		virtual int n_in() const { return m_l.nx('x'); }
		virtual int n_out() const { return m_l.nx('y'); }
		virtual int save2(SvArg * sv);
		virtual void spec_debug() { log("=======calcbox:%p=======",this); m_l.dump(stderr); }
		virtual int df_ui_ix() const { return 1; }
		virtual bool io_alias_perm() const { return false; }

		int n_tmp() const { return m_l.nx('z'); }
		void set_n_in(int n) { m_l.set_nx(n); }
		void set_n_out(int n) { m_l.set_l('y', n); }
		void set_n_tmp(int n) { m_l.set_l('z', n); }
		const char * exp_str(int ty, int i) { CalcExpr * e = m_l.ex2(ty,i); return e ? e->str(): 0; }
		void set_exp(int ty, int i, const char * s) { CalcExpr * e = m_l.ex2(ty,i);
			if (e) e -> set(s); else log("WARN: invalid exp i=%d ty=0x%x(%c)", i, ty, ty); }
		int ec(int ty, int i) { return m_l.ex2(ty,i) -> ec(); }
		int ecg(int ty, int i) { int k = ec(ty,i); return k ? k : ' '; }
		void upd_xyz(int c, int n);
		void upd_line(int t, int i, int flg); // 1:ent
	protected:
		void set_model();
		CalcEL m_l;
};

class CalcBoxModel : public BoxModel {
        public: 
		CalcBoxModel(CalcStkOp * oopp, int nnxx) : op(oopp), nx(nnxx) {}
                virtual ~CalcBoxModel() { if (op) delete[] (op); }
		virtual BoxInst * mk_box();
		CalcStkOp * op; int nx;
};

class CalcBoxInst : public BoxInst {
	public:
		CalcBoxInst(CalcBoxModel * mdl) : m_m(mdl), m_op(mdl->op) { BoxModel::ref(mdl); }
		virtual int calc(int inflg, double** inb, double** outb, int n) {
			return calc_run_b(m_op, inflg, inb, outb, n, m_m->nx); }
		virtual ~CalcBoxInst() { if (m_m) BoxModel::unref(m_m); }
	protected:
		CalcBoxModel * m_m;
		CalcStkOp * m_op;
};

int CalcBoxGen::cmd(CmdBuf* cb) {
	char * c0 = cb->a1(), *c1;
	int k, v, v0, cx;
	switch (*c0) {
		case 'X': case 'Y': case 'Z':
			v0 = v = m_l.nx(cx = *c0 + 32);
			if (!intv_cmd(&v, c0+1, (*c0=='Y'), 30))  return 0;
			unset_model();
			if (*c0=='X') m_l.set_nx(v); else m_l.set_l(cx, v);
			if (*c0!='Z') m_node->nio_change();
			if (wnfl()) { 
				upd_xyz(cx, v);
				if (cx!='x') for (int i=v0; i<v; i++) upd_line(cx, i, 1);
				if (cx!='y' && v0!=v) {
					for (int i=0; i<n_tmp(); i++) upd_line('z', i, 0);
					for (int i=0; i<n_out(); i++) upd_line('y', i, 0);
				}}
			return 1;
		case 'y':
		case 'z':
			unset_model();
			k = atoi(c0+1);
			c1 = cb->tok(); if (!c1) return BXE_NOARG;
			set_exp(*c0, k, c1 ? c1 : "");
			if (wnfl()) upd_line(*c0, k, 0);
			return 0;
		default:
			return BXE_UCMD;
	}}

void CalcBoxGen::set_model() { m_model = new CalcBoxModel(m_l.code(), m_l.nx('x')); }

int CalcBoxGen::save2(SvArg * sv) {
	BXSV2_HEAD;
	CHKERR(xprintf(f,"X$X%d\nX$Y%d\nX$Z%d\n",n_in(), n_out(), n_tmp()));
	for (int i=0; i<n_tmp(); i++) {
		CHKERR(xprintf(f,"X$z%d$%s\n", i, m_l.ex2('z',i)->str())); }
	for (int i=0; i<n_out(); i++) {
		CHKERR(xprintf(f,"X$y%d$%s\n", i, m_l.ex2('y',i)->str())); }
	return r;
}

BoxInst * CalcBoxModel::mk_box() { return new CalcBoxInst(this); }

void CalcBoxGen::upd_xyz(int c, int n) {
	gui2.setwin(w_oid(), 'c');
	gui2.wupd_c48(c, n); if (c=='x') return;
	gui2.wupd_c0(c-32, '.'); gui2.c2('+', 48+n);
}

void CalcBoxGen::upd_line(int t, int i, int flg) {
	gui2.setwin(w_oid(), 'c');
	int e = ec(t, i);
	if (flg) gui2.wupd_s(t-32, m_l.ex2(t,i)->str(), 3*i+1);
	gui2.wupd_0(t-32, "Czzz", 3*i+2); gui2.sn("k%%:%a%:"+4*(!e), 4); gui2.c1(e?e:43);
}

void CalcBoxGen::box_window() {
	m_node->winflg_or(2048);
	gui2.cre(w_oid(), 'c'); gui2.own_title();
	int no = n_out();
	int nt = n_tmp();
	upd_xyz('x', n_in()); upd_xyz('z', nt); upd_xyz('y', no);
	for (int i=0; i<nt; i++) upd_line('z', i, 1);
	for (int i=0; i<no; i++) upd_line('y', i, 1);
}

void b_ar_init(ANode * rn) {
	qmk_box(rn, "+", QMB_ARG0(Ar2BoxAD), 0, 2, 1, "a2", "i*o*R*1", "x$y", "x+y", "Pz%%0%"); 
	qmk_box(rn, "-", QMB_ARG0(Ar2BoxSB), 0, 2, 1, "a2", "1o*", "x-y");
	qmk_box(rn, "*", QMB_ARG0(Ar2BoxML), 0, 2, 1, "a2", "1o*", "x*y");
	qmk_box(rn, "/", QMB_ARG0(Ar2BoxDV), 0, 2, 1, "a2", "1o*", "x/y");
	ANode *f1 = qmk_dir(rn, "f1"), *mx = qmk_dir(rn, "mix");
	qmk_box(f1, "abs", QMB_ARG0(Ar1Abs), 0, 1, 33, "a1", "1i*o*R1", "x", "f(x)");
	qmk_box(f1, "sqrt", QMB_ARG0(Ar1Sqrt), 0, 1, 33, "a1", "1");
	qmk_box(f1, "cbrt", QMB_ARG0(Ar1Cbrt), 0, 1, 33, "a1", "1");
	qmk_box(f1, "fq_warp", QMB_ARG0(Ar1FqWarp), 0, 1, 33, "a1", "1");
	qmk_box(mx, "xmix", QMB_ARG0(CrossMix), 0, 3, 33, "xmx", "i*R*1", "x$y$t", "qqq%%q");
	qmk_box(mx, "mix2", QMB_ARG1(SMix), 2, 4, 33, "smx", "1i*R1", "a1$b1$a2$b2");
	qmk_box(mx, "mix3", QMB_ARG1(SMix), 3, 6, 33, "smx", "1i*R1", "a1$b1$a2$b2$a3$b3");
	qmk_box(mx, "mix4", QMB_ARG1(SMix), 4, 8, 33, "smx", "1i*R1", "a1$b1$a2$b2$a3$b3$a4$b4");
}

int setbox_calc(ABoxNode * nd, BoxGen * _) { nd->m_box = new CalcBoxGen(nd); return 2; }
