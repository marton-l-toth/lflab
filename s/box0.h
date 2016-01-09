#ifndef __qwe_box0_h__
#define __qwe_box0_h__

#include <new>

class ANode; class ABoxNode; class BoxModel;
extern int sample_rate; extern double sample_length;

typedef BoxModel * (*qmb_arg_t) (char *, int);
ANode * qmk_dir(ANode * up, const char * nm);
ANode * qmk_box(ANode * up, const char * nm, qmb_arg_t qa, int k, int ni, int no, 
					     const char * cl, const char * fmt, ...);
#define QMB_ARG0(CL) (&Pr0BoxModel<CL >::cons_wrap)
#define QMB_ARG1(CL) (&Pr1BoxModel<CL >::cons_wrap)

const char * iolbl_name(int ty, int ix);
void bug(const char *, ...);
void log(const char *, ...);
int nan_unpk(char * to8, int * to32, long long xl, int flg);
inline static int sec2samp(double x) { return (int)lround((double)sample_rate * x); }
inline static double att2mul(double x) { return exp(-sample_length*x); }

#define STATELESS_BOX_0(NM) class NM : public BoxInst { public: \
	virtual int calc(int inflg, double** inb, double** outb, int n); }; \
	int     NM::calc(int inflg, double** inb, double** outb, int n)

#define STATELESS_BOX_1(NM) class NM : public BoxInst { public: int m_arg; \
	NM(int arg) : m_arg(arg) {} \
	virtual int calc(int inflg, double** inb, double** outb, int n); }; \
	int     NM::calc(int inflg, double** inb, double** outb, int n)

#define FUN1_BOX(NM, X) STATELESS_BOX_0(NM) { double x; \
	if (!(inflg&1)) return x = **inb, **outb = (X), 0; \
	double *p=*inb, *q=*outb; for (int i=0; i<n; i++) x=p[i], q[i]=(X); return 1; }

#define NAN_UNPK_X(nm, p, d, ty, a1, a2) long long nm##_ll; memcpy(&nm##_ll, p, 8); \
	        ty nm##_x[16]; int nm##_n = nan_unpk(a1, a2, nm##_ll, 0) & 255; \
        if (nm##_n>63) nm##_n = 0; for (int i=nm##_n; i<16; i++) nm##_x[i] = d
#define NAN_UNPK_8(nm,  p, d) NAN_UNPK_X(nm, p, d, char, nm##_x, 0)
#define NAN_UNPK_32(nm, p, d) NAN_UNPK_X(nm, p, d, int,  0, nm##_x)

#define BOX_CP0 ( (*inb==*outb) ? (inflg&1) : ((inflg&1)?(memcpy(*outb,*inb,8*n),1):(**outb=**inb,0)))

#define PREP_INPUT(NM, I) double * NM##p, NM##v; int NM##msk = (inflg&(1<<(I))) ? \
	(NM##p = inb[I], -1) : (NM##p = &NM##v, NM##v = *inb[I], 0)

class BoxInst {
        public:
		static void rmcon(int flg, double **pp, int n); 
                virtual int calc(int inflg, double** inb, double** outb, int n) = 0;
                virtual bool done() const { return false; };
                virtual ~BoxInst() {};
		int calc_nzo2(int ocfg, double *o0, double *o1, int inflg, double **inb, int n);
};

class BoxModel {
        public:
                static void   ref(BoxModel *p) { p && ++p->refcount; }
                static void unref(BoxModel *p);
                BoxModel() : refcount(1) {}
                BoxModel(int rc) : refcount(rc) {}
                virtual ~BoxModel() {}
                virtual BoxInst * mk_box() = 0;
                void debug0();
        private:
                int refcount;
};

template <class BXI> class Pr0BoxModel : public BoxModel {
        public: static BoxModel * cons_wrap(char *to, int k) { return new (to) Pr0BoxModel<BXI>(); }
                virtual BoxInst * mk_box() { return new BXI(); }
	private:Pr0BoxModel() : BoxModel(2) {} };

template <class BXI> class Pr1BoxModel : public BoxModel {
        public: static BoxModel * cons_wrap(char *to, int k) { return new (to) Pr1BoxModel<BXI>(k); }
                virtual BoxInst * mk_box() { return new BXI(arg); }
	private:Pr1BoxModel(int k) : BoxModel(2), arg(k) {}
                int arg; };

#endif // __qwe_box0_h__
