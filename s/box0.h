#ifndef __qwe_box0_h__
#define __qwe_box0_h__

inline void* operator new(size_t, void* p) { return p; }

class ANode; class ABoxNode; class BoxModel; class BoxGen;
extern int sample_rate; extern double sample_length;
extern BoxGen * box_bookmark[10];
extern double *zero30[30], *junk30[30];

#define QMB_A_SL(CL) (&CL##_calc)
#define QMB_A_BX(CL) (&CL::sc_ctor<CL >)

const char * iolbl_name(int ty, int ix);
void bug(const char *, ...) __attribute__ ((noreturn));
void log(const char *, ...);
int nan_unpk(char * to8, int * to32, long long xl, int flg);
void reg_bn(ANode * nd, int i);
inline static int sec2samp(double x) { return (int)lround((double)sample_rate * x); }
inline static double att2mul(double x) { return exp(-sample_length*x); }

#define BX_SCALC(F) int F(BoxInst * abxi, int inflg, double** inb, double** outb, int n)

#define STATELESS_BOX(NM) static BX_SCALC(NM##_calc)

#define FUN1_BOX(NM, X) STATELESS_BOX(NM) { double x; \
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

#define SCALC_BXI(T) T* bxi = static_cast<T*>(abxi)
#define CALC_FW(F) return (*(bxi->m_psc=(F))) (bxi, inflg, inb, outb, n)
#define LOCAL_BOX(V,M) BoxModel *V##_mdl = (M); char V##_spc[V##_mdl->size()]; \
		       BoxInst  *V = V##_mdl->place_box(V##_spc);

struct BoxInst {
	typedef int(*sc_t)(BoxInst*, int, double**, double**, int);	typedef void (*dtor_t)(BoxInst*);
	typedef int(scf_t)(BoxInst*, int, double**, double**, int);	typedef void (dtorf_t)(BoxInst*);
	static int  sc_ctor_2(BoxInst *p, BoxInst::sc_t f0, int siz);
	static void rmcon(int flg, double **pp, int n); 
	static void dtor(BoxInst *p) { if0(p->m_dt>1u) dtor2(p); }
	static void dtor2(BoxInst *p);
	static void del(BoxInst* p);
	static scf_t sc_bug, sc_zero, sc_cp0;

	inline void fdk(sc_t f, int d, int k) { m_psc=f; m_dt=d; m_arg=k; }
	inline int calc(int inflg, double** inb, double** outb, int n) {
		return (*m_psc)(this, inflg, inb, outb, n); }
//	inline int calc(int inflg, double** inb, double** outb, int n) {
//		int r = (*m_psc)(this, inflg, inb, outb, n); if (r<ERRTAB_LAST) bug("ec:%p %d",this,r); return r; }
	int calc_nzo2(int ocfg, double *o0, double *o1, int inflg, double **inb, int n);

	sc_t m_psc;
	unsigned char m_dt, m_c1, m_c2, m_c3;
	int m_arg;
};

typedef int (*qmb_arg_t) (BoxInst*, int, double**, double**, int);
ANode * qmk_dir(ANode * up, const char * nm);
ANode * qmk_box(ANode * up, const char * nm, qmb_arg_t qa, int k, int ni, int no, 
					     const char * cl, const char * fmt, ...);
class BoxModel {
        public:
		typedef int(*sc_t)(BoxInst*, int, double**, double**, int);
		static void del (BoxModel *p);
                static inline void   ref   (BoxModel *p) { if (p)    ++p->refcount; }
                static inline void unref   (BoxModel *p) { if (p && !--p->refcount) del(p); }
                BoxModel(int siz, int rc) : refcount(rc), m_size(siz) {}
                virtual ~BoxModel() {}
                virtual BoxInst * place_box(void *to) = 0;
                inline  BoxInst * mk_box() { return place_box(m_size ? malloc(m_size) : 0); }
		inline int size() const { return m_size; }
                void debug0() const;
        private:
                int refcount;
	protected:
		int m_size;
};

struct BoxInst_B0 : BoxInst { template<class B> static int sc_ctor(BoxInst *q, int,double**,double**,int) { 
						      return sc_ctor_2(q, B::sc_ini, sizeof(B)); }};

struct BoxInst_B1 : BoxInst_B0 { void *m_p0;  inline void* alloc0(int k)   { m_dt=2; return m_p0=malloc(k); }};
struct BoxInst_BU : BoxInst_B0 { dtor_t m_dtf;inline void set_dtf(dtor_t f){ m_dt=3; m_dtf = f; }};

struct Impulse1 : BoxInst_B0 { static scf_t sc_ini; };
#define sc_imp1 Impulse1::sc_ini

extern double stat_con[32], *pstat_con[32];
class ConStore {
        public:
                ConStore() : m_a(0), m_n(0), m_fh(-1), m_p(0) {}
                ~ConStore() { free(m_p); }
                int siz() const { return ((sizeof(ConStore)+7)&~7) + 8*m_n; }
                int n() const { return m_n; }
                ConStore * cp(char *to) const;
                int add(double x);
                void rm(int j);
                double  v(int j) const { return j-=32, (j<0?stat_con+32:m_p)[j]; }
                double *p(int j) const { return j-=32, (j<0?stat_con+32:m_p)+j ; }
        protected:
                short m_a, m_n, m_fh, m_rsrv;
                double * m_p;
};

class ModelPtr {
	public:
		ModelPtr() : m_p(0) {}
		ModelPtr(const ModelPtr& rhs) : m_p(rhs.m_p) { BoxModel::ref(m_p); }
		ModelPtr(BoxModel *p)	      : m_p(p)	     { BoxModel::ref(m_p); }
		~ModelPtr() { BoxModel::unref(m_p); }
		inline void cp_to(ModelPtr *pp) { BoxModel::ref(pp->m_p = m_p); }
		inline BoxModel * rawmp() { return m_p; }
		template <class M> M* mk0(int x) { M * r = new (malloc(sizeof(M))) M(x); m_p=r; return r; }
		template <class M> M* mk1(int s1, int x) { 
			int s0 = (sizeof(M)+7)&~7; char *q = (char*)malloc(s0+s1);
			M *r = new (q) M(q+s0, x); m_p = r; return r; }
		template <class M> M* mk2(int s1, int s2, int x) {
			int s0 = (sizeof(M)+7)&~7, s01 = s0+((s1+7)&~7);
			char *q = (char*)malloc(s01+s2); M *r = new (q) M(q+s0,q+s01,x); m_p=r; return r; }
		template <class M> M* mkc1(const ConStore * cs, int s1, int x) { 
			int s0 = (sizeof(M)+7)&~7, s0c = s0 + cs->siz(); char *q = (char*)malloc(s0c + s1);
			M *r = new (q) M(cs->cp(q+s0), q+s0c, x); m_p=r; return r; }
		void z()  { if (m_p) BoxModel::unref(m_p), m_p = 0; }
		void z1() {          BoxModel::unref(m_p), m_p = 0; }
		bool nz() const { return !!m_p; }
		void debug() const;
		BoxInst * mk_box() { return m_p->mk_box(); }
		inline int size() const { return m_p->size(); }
		void mk_boxv(char* to, int n) { for (int i=0,sz=size();i<n;i++) m_p->place_box(to), to+=sz; }
	private:
		BoxModel * m_p;
};

#endif // __qwe_box0_h__
