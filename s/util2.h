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

// sh0: 0 10sx 110sxx 1110sxxxx 11110sx6 11111sx15
// sh1: 00 01sx 10sxxx 110sx6 111sx15
// sh2: 00 01sxx 10sxxxxx 110sx9 111sx15

int b91_cost0(const short *q, int n), b91_cost1(const short *q, int n), b91_cost2(const short *q, int n);

class B91Reader {
        public: 
                void init(char *p) { s=p; bits = 0;}
                bool eof() { return *s<36 || *s>126; }
                int get_bit() { if (!bits) get_bit_2();
				int r = (cur&1); cur>>=1; --bits; return r; }
		void get_bit_2();
                int get_bebin(int n);
                int get_short0();
                int get_short1();
                int get_short2();
                int get_short_k(int k);
        protected:
                char *s;
                int bits;
                int cur;
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

#endif // __qwe_util2_h__
