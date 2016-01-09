#ifndef __qwe_uc1_h__
#define __qwe_uc1_h__

#define LOGPIPES "uxptk"
#define CMDPIPES "gi"
#define BVFOR_JMC(X) unsigned int j, m; for (m = (X); j=__builtin_ffs(m), j-- > 0; m &= ~(1u<<j))
#define TRK_DEF_GWFR ((const unsigned char*)"$%\0")

typedef union _sthg {
	double d; void *p; int i[2]; short s[4]; char c[8]; float f[2];
} sthg;

typedef struct _chtab {
	char tab[96];
	int n;
} chtab;

static inline int    min_i(int    x, int    y) { return x<y ? x : y; }
static inline double min_d(double x, double y) { return x<y ? x : y; }
static inline int    max_i(int    x, int    y) { return x>y ? x : y; }
static inline double max_d(double x, double y) { return x>y ? x : y; }
static inline int ivlim(int x, int m, int M) { return x<m ? m : (x>M ? M : x); }
static inline void clearbuf(double *p, int n) { double *p1=p+n; while(p<p1) *(p++) = 0.0; }

static inline const char * str0(const char *s) { return s?s:"(null)"; }
static inline double cut300(double x) { return fabs(x)>=1e-300 ? x : 0.0; }

static inline int is_hx(int c) { return (unsigned int)(c-48)<10u || (unsigned int)((c|32)-97)<6u; }
static inline int hxd2i(int c) { return 15&(c+9*(c>>6)); }
static inline int hex2(const char*s) { return 16*hxd2i(*s)+hxd2i(s[1]); }
static inline int hexc1(int x) { return x + 48 + (((9-x)>>4)&7); }

#define HEX2(p,v) ( *(p) = hexc1(((v)>>4)&15), (p)[1] = hexc1((v)&15), (p) += 2 )

static inline int i_to_b32(int x) { return x<10 ? '0'+x : 'a'-10+x ; }

static inline int atoi_h(const char *s) { int r = 0, x = 0;
	for (*s=='-'&&s++&&(--x); *s&80; s++) r=16*r+hxd2i(*s); return (r^x)-x; }

static inline int chtab_get(chtab * p, int c) { c-=32; if(c<0||c>95) c=95; return p->tab[c]; }
static inline int s__cat(char * to, const char * s) { int i = 0; while (s[i]) to[i] = s[i], i++; return i; }
static inline int bitcnt_8(int c) { return c=(c&0x55)+((c>>1)&0x55), c = (c&0x33)+((c>>2)&0x33), (c&15)+(c>>4);}

static inline int mul40320(int x) { return (x>87?315:"\x01\x05\a#\x03\x0f\x15i\t-?"[x>>3])<<(x&7); }
#define LOG_LIST "gpdx"

#ifndef QWE_UTILC_DEF

int backup(const char *fn, int k);
int bitcnt(unsigned int x);
void chtab_ini(chtab * p, int fill);
int chtab_force(chtab * p, int c); 
int b32_to_i(int x);

int i2aA(int i); 
int aA2i(int i);
int hx5(char * to, int v);
char * doub2hx(char *s, double x); 
double hx2doub(const char *s);

int nan_desc(long long xl);
int nan_unpk(char * to8, int * to32, long long xl, int flg);
int nan_2str(char * to, long long xl);
int dbl2str(char * to, double x);
int trk_g_parse(const char * s, unsigned char * div, unsigned char * gwfr);

#else

int backup(const char *fn, int k) {
        if (!k || !fn || !*fn) return 0; int l = strlen(fn);
        char nm1[l+4], nm2[l+4];
        memcpy(nm1, fn, l); memcpy(nm2, fn, l);
        nm1[l] = nm1[l+1] = nm2[l] = nm2[l+1] = '-'; nm1[l+3] = nm2[l+3] = 0;
        while (k--) {
                if (!k) nm1[l]=0; else nm1[l+2] = 48+k;
                nm2[l+2] = 49+k;
                if (rename(nm1, nm2)<0 && errno!=ENOENT) return -17; /*EEE_ERRNO*/
        } return 0; }

int bitcnt(unsigned int x) {
        x = (x&0x55555555) + ((x>>1)&0x55555555);
        x = (x&0x33333333) + ((x>>2)&0x33333333);
        x = (x&0x0f0f0f0f) + ((x>>4)&0x0f0f0f0f);
        x += x>>16; x += x>>8; return (int) (x & 255);
}

void chtab_ini(chtab * p, int fill) { int i; p->n = 0; for (i=0; i<96; i++) p->tab[i] = fill; }

int chtab_force(chtab * p, int c) {
	c -= 32; if (c<0 || c>94) return p->tab[95];
	return p->tab[c] - p->tab[95] ? p->tab[c] : (p->tab[c] = p->n++);
}

int b32_to_i(int x) { unsigned int y = x-48u; if (y<10u) return y;
	        return ((y = (x|32) - 97u) < 22u) ? (int)y + 10 : -1; }

int i2aA(int i) {
        if (i<0) return '-';
        if (i>51) return '+';
        if (i&1) return 'A'+(i>>1);
        else i = 'a'+(i>>1);
        return (i=='l') ? '1': i;
}

int aA2i(int i) {
        if (i=='1') i='l';
        if (i>='a' && i<='z') return 2*(i-'a');
        if (i>='A' && i<='Z') return 2*(i-'A')+1;
        return -1;
}

int hx5(char * to, int v) {
        int x, r = 0;
        x = (v>>16)&15; if (x) to[r++] = hexc1(x);
        x = (v>>12)&15; if (r|x) to[r++] = hexc1(x);
        x = (v>> 8)&15; if (r|x) to[r++] = hexc1(x);
        x = (v>> 4)&15; if (r|x) to[r++] = hexc1(x);
        x =   v  &  15;          to[r++] = hexc1(x);
        return r;
}

double hx2doub(const char *s) {
        char v[8]; double r; int i;
        for (i=0; i<8; i++) v[i] = 16*hxd2i(s[2*i]) + hxd2i(s[2*i+1]);
        memcpy(&r, v, 8); return r;
}

char * doub2hx(char *s, double x) {
        static char buf[17]; if (!s) s = buf;
        char v[8]; memcpy(v, &x, 8);
        int i; for (i=0; i<8; i++)
                s[2*i] = hexc1((v[i]>>4)&15), s[2*i+1] = hexc1(v[i]&15);
        return s;
}

/* 6 5 4 4 4 4 4 4 4 4
   3 0 9 8 7 6 5 4 3 2
   s 0 0 1 0 0 2x23
   s 0 0 1 0 1 0 0 0 lb42  6  7 3
   s 0 0 1 0 1 0 0 1 hex40 10 4 6
   s 0 0 1 0 1 0 1 0 oct42 14 3 7
   s 0 0 1 0 1 0 1 1 1x32
   s 0 0 1 0 1 1 1 hex44
   s 0 0 1 1 0 0 set45
   s 0 0 1 1 0 1 oct45
   s 0 0 1 1 1 0 5*9
   s 0 0 1 1 1 1 9*5
   s 0 1 7*7
   s 1 0 0 3*16
   s 1 0 1 4*12
   s 1 1 0 6*8
   s 1 1 1 8*6    */

static long long nan_pkdat[16] = {
        0x7ff947e000000000, 0x7ff9580000000020, 0x7ff9000000000017, 0x7ffc000000000010,
        0x7ffd00000000000c, 0x7ff9c00000000009, 0x7ffe000000000008, 0x7ffa000000000007,
        0x7fff000000000006, 0x7ff9e00000000005, 0x7ff9480000000004, 0x7ff9700000000004,
        0x7ff9560000000003, 0x7ff9540000000003, 0x7ff9500000000003, 0x7ff9a00000000003 };

static unsigned long long nan_dflt2[2] = {0x7ff8000000000000, 0xfff8000000000000};

#define NAN_2x23 0x7ff9000000000000
#define NAN_LBL  0x7ff9400000000000
#define NAN_1x32 0x7ff9580000000000
#define NAN_5x9  0x7ff9c00000000000
#define NAN_7x7  0x7ffa000000000000
#define NAN_SMSK 0x000007e000000000
#define NAN_MSK 0x7ff8000000000000

int nan_desc(long long xl) { /* -1:num / sg<<16 + bits<<8 + n, n: 254:dflt 253:invd */
        int ix, sg = (xl<0)<<16;
        if (~xl & NAN_MSK) return -1;
        xl &= 0x7fffffffffffffff;
        if (xl >= NAN_1x32) {
                if (xl < NAN_5x9) return ix = (int)((xl>>42)&24),
                        sg + 256*((0x4200301>>ix)&63) + ((0xb010f2d>>ix)&63);
                return sg + ( (xl<NAN_7x7)
                        ? (0x5090905 >> ((xl>>41)&16)) & 0xffff
                        : (ix = (int)((xl>>46)&28) - 8,
                          256*(((0x57bf66>>ix)&15)+1) + ((0x864377>>ix)&15) ));
        }
        if (xl<NAN_LBL) return sg + (xl<NAN_2x23 ? 253+!(xl^NAN_MSK) : 0x1702);
        ix = (int)((xl>>41)&12);
        sg += 256 * ((0x347>>ix)&7) + ((0xea6>>ix)&15);
        int b1 = (int)((xl>>36) & 127);
        while (b1&64) b1+=b1, sg--;
        return sg;
}

int nan_unpk(char * to8, int * to32, long long xl, int flg) {
        int dsc = nan_desc(xl), n = dsc&255; if (n>63) return dsc;
        int i, k, bits = (dsc>>8)&63, sg = dsc>>16;
        if (bits==1) return -4; // hmm.. set?
        if (!n) return 0;
        int sgm = (1<<(bits-1)) &- sg;
        long long msk = (1LL << bits) - 1LL;
        if (to8 && (!to32 || bits<8))
                for (i=0; i<n; i++, xl>>=bits) k = (int)(xl&msk), to8 [i] = k|-(k&sgm);
        else    for (i=0; i<n; i++, xl>>=bits) k = (int)(xl&msk), to32[i] = k|-(k&sgm);
        return dsc;
}

int nan_2str(char * to, long long xl) {
        union { int i[16]; char c[64]; } uu;
        int dsc = nan_unpk(uu.c, uu.i, xl, 0); if (dsc<0) return dsc;
        int i, j, n = dsc & 255, bits = (dsc>>8)&255, sg = dsc >> 16;
        if (n==254) return memcpy(to, "-nan"+1-sg, 3+sg), 3+sg;
        if (n==253) return memcpy(to, "-nan???"+1-sg, 6+sg), 6+sg;
        if (bits==1) return memcpy(to, "{tbd...}", 8), 8;
        if (bits<8 && n<7) return to[0]=to[n+1]=34, memcpy(to+1, uu.c, n), n+2;
        *to = '['; j = 1;
        if (bits<8) for (i=0; i<n; i++) j += sprintf(to+j, " %d"+!i, uu.c[i]);
        else        for (i=0; i<n; i++) j += sprintf(to+j, " %d"+!i, uu.i[i]);
        to[j] = ']'; return j+1;
}

int dbl2str(char * to, double x) {
        long long xl; memcpy(&xl, &x, 8);
        int r = nan_2str(to, xl); return r<0 ? sprintf(to, "%.15g", x) : r;
}

int trk_g_parse(const char * s, unsigned char * div, unsigned char * gwfr) {
	int i; for (i=0; s[i]>47; i+=4) *(s[i]=='_' ? gwfr+(s[i+1]&3) : div+hex2(s+i)) = hex2(s+i+2);
	return i + (s[i]=='.');   }

#endif // QWE_UTILC_DEF
#endif // __qwe_uc1_h__
