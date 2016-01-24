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

int b91_cost(int c, int k) { if ((k=abs(k)) > 32768) return -1;
        for (int *p = costtab_123[c]; 1; p+=2) if (k<=*p) return p[1]; }

int B91Reader::get_bit() {
        if (!bits) {
                if (*s<36 || *s>126) cur=0, bits=-1;
                else if (s[1]<36 || s[1]>126) cur = *(s++)-36, bits=13;
                else cur = (*s-36)*91 + (s[1]-36), s+=2, bits=13;
        }
        int r = (cur&1); cur>>=1; --bits;
        return r;
}

void B91Writer::put_int6(int x) {
        if (!x) { put_bit(0); put_bit(0); return; }
        int xa = abs(x), xs = (x<0);
        if (xa>=181398528) { printf("put_int6: value (%d) too big", xa); return; }
        int k = 0; while (2*xa > powtab6[k]) ++k;
        xa += xs*(powtab6[k]/2) - 1;
        int dig[12];
        for (int i=k-1; i>=0; i--) dig[i] = xa%6, xa=xa/6;
        for (int i=0; i<k; i++) {
                int d = dig[i]+2;
                put_bit(d>>2); put_bit((d>>1)&1), put_bit(d&1);
        }
        put_bit(0); put_bit(0);
}

int B91Reader::get_int6() {
        int nd = 0, acc = 0;
        while(1) {
                int b1 = get_bit();
                int b2 = get_bit();
                if (!(b1|b2)) break;
                int dig = 4*b1 + 2*b2 + get_bit() - 2;
                ++nd; acc = 6*acc + dig;
        }
        if (!nd) return 0;
        int half = powtab6[nd] / 2;
        return (acc<half) ? acc+1 : half - acc - 1;
}

int B91Reader::get_bebin(int n) {
        int r = 0; for (int i=0; i<n; i++) r = r+r+get_bit();    return r; }

void B91Writer::put_bebin(int x, int n) {
        for (int i=n-1; i>=0; i--) put_bit((x>>i)&1); }

void B91Writer::put_short0(int x) {
        if (!x) { put_bit(0); return; };
        int xa = abs(x), xs = (x<0); put_bit(1);
        if (xa>32768) { printf("put_smallint: value (%d) too big", xa); return; }
        if (xa <= 2) { put_bit(0, xs, xa==2); return; }
        put_bit(1); if (xa <= 6) { put_bit(0, xs); put_bebin(xa-3, 2); return; }
        put_bit(1); if (xa <= 22) { put_bit(0, xs); put_bebin(xa-7, 4); return; }
        put_bit(1); if (xa <= 86) { put_bit(0, xs); put_bebin(xa-23, 6); return; }
        put_bit(1, xs); put_bebin(xa-87, 15);
}

void B91Writer::put_short1(int x) {
        if (!x) { put_bit(0); put_bit(0); return; };
        int xa = abs(x), xs = (x<0);
        if (xa>32768) { printf("put_smallint: value (%d) too big", xa); return; }
        if (xa <= 2) { put_bit(0); put_bit(1); put_bit(xs); put_bit(xa==2); return; }
        put_bit(1);
        if (xa <= 10) { put_bit(0); put_bit(xs); put_bebin(xa-3, 3); return; }
        put_bit(1);
        if (xa <= 74) { put_bit(0); put_bit(xs); put_bebin(xa-11, 6); return; }
        put_bit(1); put_bit(xs); put_bebin(xa-75, 15);
}

void B91Writer::put_short2(int x) {
        if (!x) { put_bit(0, 0); return; };
        int xa = abs(x), xs = (x<0);
        if (xa>32768) { printf("put_smallint: value (%d) too big", xa); return; }
        if (xa <= 4) { put_bit(0, 1, xs); put_bebin(xa-1, 2); return; }
        put_bit(1); if (xa <= 36) { put_bit(0, xs); put_bebin(xa-5, 5); return; }
        put_bit(1); if (xa <= 548) { put_bit(0, xs); put_bebin(xa-37, 9); return; }
        put_bit(1, xs); put_bebin(xa-549, 15);
}

void B91Writer::put_short_k(int k, int x) {
        if (!k) { put_short0(x); return; }
        if (k==1) put_short1(x); else put_short2(x);
}

int B91Reader::get_short0() {
        int sg=1, bb;
        if (!get_bit()) return 0;
        if (!get_bit()) bb = 0x11;
        else if (!get_bit()) bb = 0x32;
        else if (!get_bit()) bb = 0x74;
        else bb = get_bit() ? 0x57f : 0x176;
        bb && get_bit() && (sg = -1);
        return sg * ((bb>>4) + get_bebin(bb&15));
}

int B91Reader::get_short1() {
        int sg, base, bits;
        if (!get_bit()) {
                if (!get_bit()) return 0;
                sg=get_bit(); base=1; bits=1;
        } else {
                if (!get_bit()) sg=get_bit(), base=bits=3;
                else get_bit() ? (base=75,bits=15) : (base=11,bits=6), sg=get_bit();
        }
        return (1-sg-sg) * (base + get_bebin(bits));
}

int B91Reader::get_short2() {
        int sg=1, bb;
        if (!get_bit()) bb = 0x12 &- get_bit();
        else if (!get_bit()) bb = 0x55;
        else bb = get_bit() ? 0x225f : 0x259;
        bb && get_bit() && (sg = -1);
        return sg * ((bb>>4) + get_bebin(bb&15));
}

int B91Reader::get_short_k(int k) {
        if (!k) return get_short0();
        return k==1 ? get_short1() : get_short2();
}


