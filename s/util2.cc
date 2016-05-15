#include <fcntl.h>
#include <unistd.h>
#include "node.h"
#include "util2.h"
#include "glob.h"
#include "gp.h"
#include "pt.h"
#include "cfgtab.inc"

void bye(int x) {
	if (x&256 || (glob_flg&(GLF_FINAL_ASV|GLF_INI0|GLF_INI1))) fprintf(stderr, "skipping autosave\n");
	else glob_flg|=GLF_FINAL_ASV, Node::save_batch(Node::root(), "//"+!(glob_flg&GLF_FSTATE), NOF_FORCE);
	delete (Gnuplot::sg());
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

/////////////////////////////// voltab compression /////////////////////////////////////////////////

static int powtab6[] = {1, 6, 36, 216, 1296, 7776, 46656, 279936, 1679616, 10077696, 60466176, 362797056};

// sh0: 0 10sx 110sxx 1110sxxxx 11110sx6 11111sx15
// sh1: 00 01sx 10sxxx 110sx6 111sx15
// sh2: 00 01sxx 10sxxxxx 110sx9 111sx15
static int costtab0[] = { 0,1,  2,4,  6,6,   22,9,  86,12,  33333,21 };
static int costtab1[] = { 0,2,  2,4,  10,6,  74,10,  33333, 19 };
static int costtab2[] = { 0,2,  4,5,  36,8,  548,13, 33333, 19 };
static int * costtab_123[3] = { costtab0, costtab1, costtab2 };

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

void B91Reader::get_bit_2() {
	if (*s<36 || *s>126) cur=0, bits=-1;
	else if (s[1]<36 || s[1]>126) cur = *(s++)-36, bits=13;
	else cur = (*s-36)*91 + (s[1]-36), s+=2, bits=13; }

int B91Reader::get_bebin(int n) {
        int r = 0; for (int i=0; i<n; i++) r = r+r+get_bit();    return r; }

static unsigned char rev8[256] = {
          0, 128,  64, 192,  32, 160,  96, 224,  16, 144,  80, 208,  48, 176, 112, 240, 
          8, 136,  72, 200,  40, 168, 104, 232,  24, 152,  88, 216,  56, 184, 120, 248, 
          4, 132,  68, 196,  36, 164, 100, 228,  20, 148,  84, 212,  52, 180, 116, 244, 
         12, 140,  76, 204,  44, 172, 108, 236,  28, 156,  92, 220,  60, 188, 124, 252, 
          2, 130,  66, 194,  34, 162,  98, 226,  18, 146,  82, 210,  50, 178, 114, 242, 
         10, 138,  74, 202,  42, 170, 106, 234,  26, 154,  90, 218,  58, 186, 122, 250, 
          6, 134,  70, 198,  38, 166, 102, 230,  22, 150,  86, 214,  54, 182, 118, 246, 
         14, 142,  78, 206,  46, 174, 110, 238,  30, 158,  94, 222,  62, 190, 126, 254, 
          1, 129,  65, 193,  33, 161,  97, 225,  17, 145,  81, 209,  49, 177, 113, 241, 
          9, 137,  73, 201,  41, 169, 105, 233,  25, 153,  89, 217,  57, 185, 121, 249, 
          5, 133,  69, 197,  37, 165, 101, 229,  21, 149,  85, 213,  53, 181, 117, 245, 
         13, 141,  77, 205,  45, 173, 109, 237,  29, 157,  93, 221,  61, 189, 125, 253, 
          3, 131,  67, 195,  35, 163,  99, 227,  19, 147,  83, 211,  51, 179, 115, 243, 
         11, 139,  75, 203,  43, 171, 107, 235,  27, 155,  91, 219,  59, 187, 123, 251, 
          7, 135,  71, 199,  39, 167, 103, 231,  23, 151,  87, 215,  55, 183, 119, 247, 
         15, 143,  79, 207,  47, 175, 111, 239,  31, 159,  95, 223,  63, 191, 127, 255 };

void B91Writer::put_bebin(int x, int n) {
	while (1) {
		int j,k = min_i(n, 13-bits); 
		if (k<8) { j = (x>>(n-k)) & ((1<<k)-1); j = (int)rev8[j] >> (8-k);
			   cur |= j<<bits; if ((bits+=k)==13) flush13(); if (!(n-=k)) return; }
		else     { cur |= rev8[(x>>(n-=8))&255] << bits; bits += 8;
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
        int base, bits;
        if (!get_bit()) return 0;
        if (!get_bit()) base=bits=1;
        else if (!get_bit()) base=3, bits=2;
        else if (!get_bit()) base=7, bits = 4;
        else get_bit() ? (base=0x57, bits=15) : (base=0x17, bits=6); 
	return get_bit() ? -base-get_bebin(bits) : base+get_bebin(bits);
}

int B91Reader::get_short1() {
        int base, bits;
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
        int base, bits;
        if (!get_bit()) { if (!get_bit()) return 0; base=1, bits=2; }
        else if (!get_bit()) { base = bits = 5;  }
        else { if (get_bit()) base=0x225, bits=15; else base=0x25, bits=9; }
        return get_bit() ? -base-get_bebin(bits) : base+get_bebin(bits);
}

int B91Reader::get_short_k(int k) {
        if (!k) return get_short0();
        return k==1 ? get_short1() : get_short2();
}


