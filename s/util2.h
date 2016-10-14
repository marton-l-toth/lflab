#ifndef __qwe_util2_h__
#define __qwe_util2_h__

class Clock {
        public: 
                int reset() {
                        int sec = tv.tv_sec;
                        int usec = tv.tv_usec;
                        gettimeofday(&tv, 0);
                        return 1000000*(tv.tv_sec - sec) + tv.tv_usec - usec;
                }
                int get() {
                        struct timeval t2;
                        gettimeofday(&t2, 0);
                        return 1000000*(t2.tv_sec - tv.tv_sec) + t2.tv_usec - tv.tv_usec;
                }
        protected:
                struct timeval tv;
};

extern void gui_tlog(int,int);
class BufClock {
        public: 
                int ini(int bits, int flg, int jst, int jsn); // flg: 1-monot.raw
                void bcfg(int rate, int bs, int bs2, int fmax);
                int set(int t);
                int ev(int c), ev2(int c, int a = 0), j0(), j1(int cont);
		int sel(int zf), f2play(int qf);
                inline int t() const { return m_t; }
                inline int add(int t) { return m_t += t; }
                inline unsigned int *pt(int j=0) { return m_buf+((m_ix+2*j  ) & m_ix_msk); }
                inline unsigned int *pa(int j=0) { return m_buf+((m_ix+2*j+1) & m_ix_msk); }
                inline int err(int z=1) { return z ? (z=m_err, m_err=0, z) : m_err; }
                inline int tcond(struct timespec *p, int min_ms = 0x80000000) {
                        return (1000*(m_ts.tv_sec-p->tv_sec)+(m_ts.tv_nsec-p->tv_nsec)/1000000>=min_ms)
                                && (memcpy(p, &m_ts, sizeof(m_ts)), 1); }
		inline int f2ns(int nf) { return (nf*m_cf_ns16f+8)>>4; }
		inline int set_f(int nf) { return set(f2ns(nf)); }
		inline int add_f(int nf) { return m_t += f2ns(nf); }
		inline void cls() { ev('x'); m_t = m_cf_half-1; }
		void gcond() {  extern void gui_tlog(int,int); if (++m_gcnt<4) return (void)ev(65+m_gcnt);;
				m_gcnt=0; gui_tlog(m_g_ix&m_ix_msk, m_ix+2-m_g_ix); m_g_ix=m_ix+2; ev('#'); }
		int wrk(int op, int n = 0);
        protected:
		inline int rnd5() { int r = m_rnd>31 ? m_rnd : random()|0x40000000; m_rnd=r>>5; return r&31; }
                inline void jvi() { m_j_ct = m_cf_jst, m_j_cn = m_cf_jsn-1, m_j_max = 0; }
                int jwr(int cont);
                struct timespec m_ts;
                unsigned int *m_buf;
                unsigned int m_ix, m_ix_msk, m_seq_sh;
                int m_j_cn, m_j_ct, m_j_max;
                int m_cf_full, m_cf_empty, m_cf_half, m_cf_stmin, m_cf_jtmin, m_cf_jst, m_cf_jsn;
                int m_cf_nspf, m_cf_ns16f, m_cf_fmax;
                clockid_t m_ty;
                int m_bits, m_t, m_gcnt, m_g_ix, m_err, m_rnd;
};

// sh0: 0 10sx 110sxx 1110sxxxx 11110sx6 11111sx15
// sh1: 00 01sx 10sxxx 110sx6 111sx15
// sh2: 00 01sxx 10sxxxxx 110sx9 111sx15

int b91_cost0(const short *q, int n), b91_cost1(const short *q, int n), b91_cost2(const short *q, int n);

class B91Reader {
        public: 
                void init(char *p) { s=p; bits = 0; cur=0;}
                bool eof() { return *s<36 || *s>126; }
                inline int get_bit()   { int r = (cur&1); cur>>=1; --bits; return r; }
		inline void get13() {
			if (*s<36 || *s>126) bits+=13;
			else if (s[1]<36 || s[1]>126) cur |= (*(s++)-36)<<bits, bits+=13;
			else cur |= ((*s-36)*91 + (s[1]-36))<<bits, s+=2, bits+=13; }
		inline void fill20() { if (bits<20) get13(); if (bits<20) get13(); }
		inline void fill13() { if (bits<20) get13(); }
                inline int get_bebin(int n) {
			int r=0; for (int i=0; i<n; i++) r = 2*r+(cur&1), cur>>=1; return bits-=n, r; }
                int get_short0();
                int get_short1();
                int get_short2();
        protected:
                char *s;
                int bits;
                unsigned int cur;
};

class B91Writer {
        public: 
                B91Writer() : bits(0), cur(0) {}
                void put_bit(int b) { cur |= (b<<(bits++)); if (bits==13) flush13(); }
                void put_bit(int b1, int b2) { put_bit(b1); put_bit(b2); }
                void put_bit(int b1, int b2, int b3) { put_bit(b1); put_bit(b2); put_bit(b3); }
                void put_bebin(int x, int n);
                void put_short_tpn(int ty, const short *p, int n);

                char * get_str() { flush13(); buf.add(0); return buf.forget(); }
                int n_bytes() { flush13(); return buf.n(); }
        protected:
                void flush13() { if (!bits) return; buf.add(36+cur/91, 36+cur%91); cur = bits = 0; }
                LWArr<char> buf;
                int bits, cur;
};

// arg: d<dec>|x<hex>|+i/-i|(i:0(1)|1(step1)|2(step2)|3(min/max))|<dec> 
// ret: 0:nothing happened 1:changed 2+4*flg:conf 

struct cfg_ent;
int intv_cmd    (int           *p, const char * arg, int min, int max, int mul4=0x01010101);
int intv_cmd_c  (         char *p, const char * arg, int min, int max, int mul4=0x01010101);
int intv_cmd_sc (  signed char *p, const char * arg, int min, int max, int mul4=0x01010101);
int intv_cmd_uc (unsigned char *p, const char * arg, int min, int max, int mul4=0x01010101);
int intv_cmd_cfg(cfg_ent       *p, const char * arg,                   int mul4=0x01010101);
int intv_cmd_b(unsigned int *bv, int b0, int nb, const char * arg, int mul4=0x01010101, int min=0, int max=0);

int  packflg(         int flg, const int * mv); // mv[0]:mask mv[1]:dflt-val
void unpkflg(int *to, int fpk, const int * mv);

extern BufClock clk0;

#endif // __qwe_util2_h__
