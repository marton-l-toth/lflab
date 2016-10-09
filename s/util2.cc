#include <fcntl.h>
#include <unistd.h>
#include <time.h>

#include "node.h"
#include "util2.h"
#include "glob.h"
#include "pt.h"
#include "cfgtab.inc"

void bye(int x) {
	if (x&256 || (glob_flg&(GLF_FINAL_ASV|GLF_INI0|GLF_INI1))) fprintf(stderr, "skipping autosave\n");
	else glob_flg|=GLF_FINAL_ASV, Node::save_batch(Node::root(), "//"+!(glob_flg&GLF_FSTATE), NOF_FORCE);
	if (x&512) pt_iocmd_sn("r\n", 2);
	exit(x&255);
}

void bug(const char * fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	fprintf(stderr,"%c[1;33m",27);
	fprintf(stderr,"BUG!!!\n");
	vfprintf(stderr, fmt, ap);
	fprintf(stderr,"\n");
	perror("last unix error");
	fprintf(stderr,"%c[0m\n",27);
	// Node::save_batch("__autosave_bug", NOF_FORCE);
	va_end(ap);
	abort();
}

static unsigned int icv_stockmul[4] = { 0xf0f00501, 0xf0190501, 0xff400801, 0xf0640a01 };
template <class T>
int icmd_t(T * p, const char * arg, int min, int max, int mul4=0x01010101) {
        int x; if ((mul4|3)==-1) mul4 = (int)icv_stockmul[mul4&3];
        switch(arg[0]) {
                case '+': x = *p + ((mul4>>(8*(arg[1]&3)))&255); break;
                case '-': x = *p - ((mul4>>(8*(arg[1]&3)))&255); break;
                case ',': return 2 + 4 * (arg[1] & 3);
                case 'x': x = atoi_h(arg+1); break;
                case 'd': ++arg;
                default: x = atoi(arg); break;
        }
        if (x<min) x=min; else if (x>max) x=max;
        return (x!=*p) && (*p=x, 1);
}

int intv_cmd    (int           *p, const char * arg, int min, int max, int mul4)  {
   return icmd_t<int>(          p,		arg,	 min, 	  max, 	   mul4); }
int intv_cmd_c  (char 	       *p, const char * arg, int min, int max, int mul4)  {
   return icmd_t<char>(         p,		arg,	 min, 	  max, 	   mul4); }
int intv_cmd_sc (  signed char *p, const char * arg, int min, int max, int mul4)  {
   return icmd_t<  signed char>(p,		arg,	 min, 	  max, 	   mul4); }
int intv_cmd_uc (unsigned char *p, const char * arg, int min, int max, int mul4)  {
   return icmd_t<unsigned char>(p,		arg,	 min, 	  max, 	   mul4); }
int intv_cmd_cfg(cfg_ent       *q, const char * arg, 		       int mul4)  {
   return icmd_t<int>(        &q->i,	        arg,    q->i_m,  q->i_M,   mul4); }

int intv_cmd_b(unsigned int *bv, int b0, int nb, const char * arg, int mul4, int min, int max) {
	int msk = (1<<nb)-1, v = (*bv>>b0) & msk, v0 = v;
	return intv_cmd(&v, arg, min, max?max:msk, mul4) && (*bv ^= (unsigned int)(v0^v)<<b0, 1); 
}

/////////////////////////////// voltab compression /////////////////////////////////////////////////

// sh0: 0 10sx 110sxx 1110sxxxx 11110sx6 11111sx15
// sh1: 00 01sx 10sxxx 110sx6 111sx15
// sh2: 00 01sxx 10sxxxxx 110sx9 111sx15

#define B91_SEL0(X,A,B,C,D,E,F) ( (X<=2) ? (X?B:A) : ((X<=22)?(X<=6?C:D): (X<=86?E:F)) )
#define B91_SEL1(X,A,B,C,D,E) ( (X<=2) ? (X?B:A) : ((X<=10)?C:(X<= 74?D:E)) )
#define B91_SEL2(X,A,B,C,D,E) ( (X<=4) ? (X?B:A) : ((X<=36)?C:(X<=548?D:E)) )

int b91_cost0(const short * q, int n) {
	int v = 0; for (int k,i=0; i<n; i++) k = abs(q[i]), v += B91_SEL0(k, 1, 4, 6, 9, 12, 21);
	return v; }

int b91_cost1(const short * q, int n) {
	int v = 0; for (int k,i=0; i<n; i++) k = abs(q[i]), v += B91_SEL1(k, 2, 4, 6, 10, 19);
	return v; }

int b91_cost2(const short * q, int n) {
	int v = 0; for (int k,i=0; i<n; i++) k = abs(q[i]), v += B91_SEL2(k, 2, 5, 8, 13, 19);
	return v; }

void B91Writer::put_bebin(int x, int n) {
	while (1) {
		int j,k = min_i(n, 13-bits); 
		if (k<8) { j = (x>>(n-k)) & ((1<<k)-1); j = (int)bitrev8[j] >> (8-k);
			   cur |= j<<bits; if ((bits+=k)==13) flush13(); if (!(n-=k)) return; }
		else     { cur |= bitrev8[(x>>(n-=8))&255] << bits; bits += 8;
			   if (bits==13) flush13(); if (!n) return; }}}

// k: 0 (4*ssh) 0 (5*tot) (21*base)
#define B91_PSH(Z) for (int i=0; i<n; i++) { \
	int x = p[i], xa = abs(x), k = (Z); \
	put_bebin(xa+(k&07777777)+((x<0)<<(k>>27)), (k>>21)&31); }

void B91Writer::put_short_tpn(int ty, const short *p, int n) { switch(ty) {
	case 0: B91_PSH(B91_SEL0(xa, 001010000000, 001040000007, 002060000055,
				     004110000671, 006140007351, 017257577651)); return;
	case 1: B91_PSH(B91_SEL1(xa, 001020000000, 001040000003, 003060000035,
				     006120001365, 017231577665)); return;
	case 2: B91_PSH(B91_SEL2(xa, 001020000000, 002050000007, 005100000173,
				     011150013733, 017231576733)); return;
	default: log("BUG: b91/kpn: ty=%d", ty); return; }}

int B91Reader::get_short0() {
        int base, bits; fill20();
        if (!get_bit()) return 0;
        if (!get_bit()) base=bits=1;
        else if (!get_bit()) base=3, bits=2;
        else if (!get_bit()) base=7, bits = 4;
        else get_bit() ? (fill13(), base=0x57, bits=15) : (base=0x17, bits=6); 
	return get_bit() ? -base-get_bebin(bits) : base+get_bebin(bits);
}

int B91Reader::get_short1() {
        int base, bits; fill20();
        if (!get_bit()) {
                if (!get_bit()) return 0;
                base=1; bits=1;
        } else {
                if (!get_bit()) base=bits=3;
                else get_bit() ? (base=75,bits=15) : (base=11,bits=6);
        }
        return get_bit() ? -base-get_bebin(bits) : base+get_bebin(bits);
}

int B91Reader::get_short2() {
        int base, bits; fill20();
        if (!get_bit()) { if (!get_bit()) return 0; base=1, bits=2; }
        else if (!get_bit()) { base = bits = 5;  }
        else { if (get_bit()) base=0x225, bits=15; else base=0x25, bits=9; }
        return get_bit() ? -base-get_bebin(bits) : base+get_bebin(bits);
}

///// buffer clock
#define BLLN 1000000000
#define BCLK00 struct timespec ts; memcpy(&ts, &m_ts, sizeof(ts))
#define BCLK01 int dn = ((int)m_ts.tv_nsec&~127) - ((int)ts.tv_nsec&~127), ds = m_ts.tv_sec - ts.tv_sec, \
                   cf = (ds|(dn&0x80000000))
#define BCLK0 BCLK00; clock_gettime(m_ty, &m_ts); BCLK01
#define BCLK1(E) (dn = BLLN + ((ds==1) ? (dn<=0 ? dn : ((E)|=1,128)) : (ds=1+(ds<0), (E)|=ds, ds<<7)))
#define BCLKN(C) (m_ix+=2, m_buf[m_ix&m_ix_msk]  = ((m_ix<<m_seq_sh)&(3u<<30)) + (C))
#define BCLKJ(R) if (dn>m_j_max) m_j_max=dn; return m_j_ct-=dn, (R)

int BufClock::ini(int bits, int flg, int jst, int jsn) {
	static clockid_t cv[2] = { CLOCK_MONOTONIC, CLOCK_MONOTONIC_RAW };
        if (!(m_buf = (unsigned int *)map_wdir_shm('@', 8<<bits, 3))) return -1;
	m_bits = bits, m_ix_msk = (2<<bits)-1, m_seq_sh = 29-bits; m_cf_jst=jst, m_cf_jsn=jsn;
	m_buf[0] = 0x40000021; m_g_ix = m_ix = 2u<<bits; m_gcnt = 0;
	if (clock_gettime(m_ty=cv[1& flg], &m_ts)>=0) return 0;
	if (clock_gettime(m_ty=cv[1&~flg], &m_ts)>=0) return 1; else return -2;
}

int BufClock::f2play(int qf) {
	int nf, t0 = ev('P'+qf);
	if (t0>m_cf_half) { if (qf && t0>2*m_cf_full) set(m_cf_half-1); nf = 0; }
	else if ((nf=(m_cf_full-t0)/m_cf_nspf) > m_cf_fmax) { nf = m_cf_fmax; if (qf && t0<0) set(m_cf_half); }
	return *pa() = nf;
}

int BufClock::set(int t) {
	BCLK00; if (clock_gettime(m_ty, &m_ts)<0) return -1; BCLK01; int e;
	return m_t = t, cf ? (e=0,BCLK1(t),e|=m_err) : (e=m_err), m_err=0, *pt()+=dn, BCLKN('K'), *pa()=t, e; }

void BufClock::bcfg(int rate, int bs, int bs2, int fmax) {
        double npf = 1e9 / (double)rate; int e,f,fe;
        m_cf_nspf  = (int)lround(     npf);   m_cf_empty = e = (int)lround(npf*(double)(  bs+bs2));
        m_cf_ns16f = (int)lround(16.0*npf);   m_cf_full  = f = (int)lround(npf*(double)(2*bs+bs2));
        m_cf_fmax = fmax; fe = f-e; m_cf_jtmin = m_cf_half = e+(fe>>1); m_cf_stmin = e+(fe>>3); }

int BufClock::ev(int c) { BCLK0; if (cf) BCLK1(m_err); return *pt()+=dn, BCLKN(c), m_t-=dn; }

int BufClock::sel(int zf) { int t = ev('s'); return *pa() = (zf||t<m_cf_stmin) ? 0 : (t-m_cf_empty)>>10; }

int BufClock::j0() { BCLK0; if (cf) BCLK1(m_err);   int r = ((m_t-=dn) >= m_cf_jtmin);
                     return *pt() += dn, BCLKN('J'+32*r), jvi(), r;     }

int BufClock::jwr(int cf) { return *pt() += (m_cf_jst-m_j_ct),
                                   *pa() = m_j_max+(m_cf_jsn-1-m_j_cn), BCLKN('J'+32*cf), jvi(), cf; }

int BufClock::j1(int cont) {
        BCLK0; int qstp = ((cont-1)|(m_t-dn-m_cf_jtmin)), qwr = (--m_j_cn)|(m_j_ct-dn);
        if (!(cf | ((qstp|qwr)&0x80000000))) { m_t -= dn; BCLKJ(1); }
        if (!cf) { m_t -= dn; BCLKJ(jwr(qstp>=0)); }
        int ec = 0; if (BCLK1(ec), ec) { m_err|=ec; m_t-=dn; BCLKJ(jwr(0)); }
        cont &= ((m_t-=dn>m_cf_jtmin)); BCLKJ( (((cont-1)|m_j_ct|m_j_cn)<0) ? jwr(cont) : cont );
}

int packflg(int flg, const int * mv) {
	int fm = mv[0], f2 = flg & fm; if (f2==mv[1]) return -1;
	int k=0, r=0; BVFOR_JM(fm) r |= ((f2>>j)&1) << (k++);   return r; }

void unpkflg(int *to, int fpk, const int * mv) {
	int m1, fm = mv[0], k = 0; BVFOR_JM(fm) m1=1<<j, (fpk&(1<<(k++))) ? (*to|=m1) : (*to&=~m1); }


