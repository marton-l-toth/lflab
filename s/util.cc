#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <signal.h>
#include <cstdarg>
#include <errno.h>

#define QWE_UTILC_DEF
#include "uc0.h"
#include "uc1.h"
#include "util.h"
#include "util2.h"
#include "errtab.inc"
#include "cfgtab.inc"
#include "glob.h"
#include "pt.h"
//#include "lwarr.h"

const double smalldiv128[130] = {0.0,0.0,1.0,0.5,0.3333333333333333,0.25,0.2,0.1666666666666667,
0.1428571428571428,0.125,0.1111111111111111,0.1,0.09090909090909091,0.08333333333333333,0.07692307692307693,
0.07142857142857142,0.06666666666666667,0.0625,0.05882352941176471,0.05555555555555555,0.05263157894736842,
0.05,0.04761904761904762,0.04545454545454546,0.04347826086956522,0.04166666666666666,0.04,0.03846153846153846,
0.03703703703703703,0.03571428571428571,0.03448275862068965,0.03333333333333333,0.03225806451612903,0.03125,
0.0303030303030303,0.02941176470588235,0.02857142857142857,0.02777777777777778,0.02702702702702703,
0.02631578947368421,0.02564102564102564,0.025,0.02439024390243903,0.02380952380952381,0.02325581395348837,
0.02272727272727273,0.02222222222222222,0.02173913043478261,0.02127659574468085,0.02083333333333333,
0.02040816326530612,0.02,0.0196078431372549,0.01923076923076923,0.01886792452830189,0.01851851851851852,
0.01818181818181818,0.01785714285714286,0.01754385964912281,0.01724137931034483,0.01694915254237288,
0.01666666666666667,0.01639344262295082,0.01612903225806452,0.01587301587301587,0.015625,0.01538461538461539,
0.01515151515151515,0.01492537313432836,0.01470588235294118,0.01449275362318841,0.01428571428571429,
0.01408450704225352,0.01388888888888889,0.0136986301369863,0.01351351351351351,0.01333333333333333,
0.0131578947368421,0.01298701298701299,0.01282051282051282,0.01265822784810127,0.0125,0.01234567901234568,
0.01219512195121951,0.01204819277108434,0.0119047619047619,0.01176470588235294,0.01162790697674419,
0.01149425287356322,0.01136363636363636,0.01123595505617977,0.01111111111111111,0.01098901098901099,
0.0108695652173913,0.01075268817204301,0.01063829787234043,0.01052631578947368,0.01041666666666667,
0.01030927835051546,0.01020408163265306,0.0101010101010101,0.01,0.009900990099009901,0.009803921568627451,
0.009708737864077669,0.009615384615384616,0.009523809523809525,0.009433962264150943,0.009345794392523364,
0.009259259259259259,0.009174311926605505,0.00909090909090909,0.009009009009009009,0.008928571428571428,
0.008849557522123894,0.008771929824561403,0.008695652173913044,0.008620689655172414,0.008547008547008548,
0.008474576271186441,0.008403361344537815,0.008333333333333333,0.008264462809917356,0.00819672131147541,
0.008130081300813009,0.008064516129032258,0.008,0.007936507936507936,0.007874015748031496,0.0078125};

char * bigblk(int n) {
	return (char*)mmap(0, n, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE, -1, 0); }

char * alloc_32k() {
	static char *cur=0, *lim = 0;
	static int bs = 65536;
	if (cur==lim) {
		cur = bigblk(bs); lim = cur+bs;
		if (bs<(1<<24)) bs += bs;
	}
	char * p = cur; cur += 32768; return p;
}

int nfa_siz = 0; char *nfa_cur = 0;
char * nf_alloc_2(int x) { if (x>32768) return (char*)malloc(x);
			   char * p = alloc_32k();  nfa_cur = p+x; nfa_siz=32768-x; return p; }

int u_sleep(int x) {
	struct timespec ts; ts.tv_sec = 0; ts.tv_nsec = 1000*x; int ec;
	do ec = nanosleep(&ts, &ts); while (ec==EINTR); return ec; }

void shuffle(int * p, int n) { for (int j, t, i=0; i<n-1; i++) 
		j = i + (random() % (n-i)), t = p[i], p[i] = p[j], p[j] = t; }

static const char*fontwid_16= "566=:?<3668=5655::::::::::55===9@;;;<:9<<55:9=<=:=;:9<;A;9<656=889:9:96::3393?::::796:9=:99:5:=";

int text_wid16(const char * s) {
	int r=0; for (; *s; s++) r += (*s<32||*s>126) ? 10 : fontwid_16[*s-32];
	return r; }

int lorem(char * to, int l) {
	int r=0; while(l>3) to[r++]='m', l-=4;
	return r + (l && (to[r]=" .ta"[l]));
}

///////// arithmetic ////////////////////////////////////////////////////////

int approx_cmp(double x, double y) {
	if (x!=x) return (y==y) ? 1 : memcmp(&x,&y,8);
	if (y!=y) return -1;
	if (x==y) return 0;
	double d = x - y;
	int r = (d<0.0) ? (d=-d, -1) : 1;
	if (d<1e-280) return 0;
	if (d*1e14 < fabs(x)) return 0;
	return r;
}

int approx_cmp_v(const double *x, const double *y, int n) {
	for (int r,i=0; i<n; i++) { if ((r = memcmp(x+i,y+i,8))) {
		if (x[i]!=x[i]) return (y[i]==y[i]) ? -1 : r;
		if (y[i]!=y[i]) return 1;
		double d = x[i] - y[i];
		int r = (d<0.0) ? (d=-d, -1) : 1;
		if (d>=1e-280 && d*1e14 >= fabs(x[i])) return r;
	}}
	return 0;
}

double ipow(double x, int k) {
	double acc; if (!k) return 1.0;
	while (!(k&1)) x *= x, k>>=1;
	acc = x;
	while (k>>=1) { x *= x;  if (k&1) acc *= x; }
	return acc;
}

///////// debug tools ///////////////////////////////////////////////////////

class QSReader : public AReader {
	public:
		QSReader(int n);
		virtual ~QSReader() { free(m_buf); }
		virtual int line(char * s);
	protected:
		int m_siz, m_j;
		char * m_buf;
};

void log_sn(const char * pref, const char * str, int len) {
	if (pref) fprintf(stderr, "%s", pref), fflush(stderr); 
	write(2, str, len); if (pref) putc(10, stderr);
}

void log(const char * fmt, ...) {
        va_list ap; va_start(ap, fmt);
	vfprintf(stderr, fmt, ap); putc(10, stderr); va_end(ap);
}

void log_n(const char * fmt, ...) {
        va_list ap; va_start(ap, fmt);
        vfprintf(stderr, fmt, ap); va_end(ap);
}

void log_smallint(int i) { switch(i) {
		case 0: log_n("o"); return;
		case 1: log_n("/"); return;
		case 2: log_n("="); return;
		case 3: log_n("Y"); return;
		case 4: log_n("X"); return;
		case 5: log_n("H"); return;
		case 6: log_n("*"); return;
		default: log_n(" %d",i); return;
}}

int shortcmp(const void *p, const void *q) { return *(const short*)p - *(const short*)q; }
void log_sortedblk(short *p, int n, bool absf, const char * pref) {
	short q[n];
	if (pref) log_n("%s",pref);
	if (absf) for (int i=0; i<n; i++) q[i] = abs(p[i]);
	else for (int i=0; i<n; i++) q[i] = p[i];
	qsort(q, n, sizeof(short), shortcmp);
	for (int i=0; i<n; i++) log_smallint(q[i]);
	log("");
}

QSReader::QSReader(int n) : m_siz(12*n+20) { m_buf=(char*)malloc(m_siz); m_j=sprintf(m_buf,"Sc%d\n",n); }
int QSReader::line(char *s) {
	int j = m_j, l = strlen(s), l11 = l/11, r = l%11;   char * q = m_buf;
	if (l11<1 || j+l+8>m_siz) return BXE_PARSE;
	if (r && (r!=1 || s[--l]!='!')) return BXE_PARSE;
	q[j] = 'S'; q[j+1] = 'C'; q[j+2] = 48+l11; memcpy(q+j+3, s, l); q[j+3+l] = 10; j+=l+4; 
	return r ? (memcpy(m_buf+j,"SC0\n",4), pt_wrk_cmd(q,j+4)) : (m_j = j, 1);
}

int qstat_op(int op) { 
	static int n0=0, n1=0; int x, rv=0;
	if (op=='+') ++n0;
	else if (++n1, x=(op-48)) rv = QSE_BUG - ( ((unsigned int)x<5u) ? x : 0 );
	log("quickstat: %d/%d", n1, n0); return rv;
}

int qstat_cfg(int op, const char * arg) {
	int len; if (!arg || !*arg) arg="64", len=2; else len = strlen(arg);
	char buf[len+3]; buf[0] = 'S'; buf[1] = op; memcpy(buf+2, arg, len); buf[len+2] = 10;
	return pt_wrk_cmd(buf, len+3); }

int qstat_chk0(AReader ** ppr, const char *s) {
	if (*ppr) return GCE_CONT2;
	int k = atoi(s); if (k<1 || k>511) return EEE_RANGE;
	*ppr = new QSReader(k); return 0;
}

///////// sorting ASCII strings /////////////////////////////////////////////

typedef void (*ssort_t) (char**, char**, int);

#define sqless(i,j) (strcmp(fr[i]+k, fr[j]+k) < 0)
static void ssort_2x(char ** to, char ** fr, int k) {
	if (sqless(1,0)) to[0]=fr[1], to[1]=fr[0];
	else to[0]=fr[0], to[1]=fr[1]; }

#define perm3(x,y,z) (to[0]=fr[x], to[1]=fr[y], to[2]=fr[z])
static void ssort_3x(char ** to, char ** fr, int k) {
	if (sqless(0,1)) {
		if (sqless(1,2)) perm3(0,1,2);
		else if (sqless(0,2)) perm3(0,2,1);
		else perm3(2,0,1);
	} else {
		if (sqless(2,1)) perm3(2,1,0);
		else if (sqless(0,2)) perm3(1,0,2);
		else perm3(1,2,0);
	}}

static void ssort_4x(char ** to, char ** fr, int k) {
	int sf[2], md[2]; sf[0] = sqless(0,1); sf[1] = 2+sqless(2,3);
	int f1 = sqless(sf[0],     sf[1]), fM = sf[f1]; md[0] = sf[f1^1];
	int f2 = sqless(sf[1]^1, sf[0]^1), fm = sf[f2]^1; md[1] = sf[f2^1]^1;
	int j = ((md[0]^md[1])==1) ? (sf[md[0]>>1] ^ md[0]) & 1
				   : sqless(md[0], md[1]);
	to[0] = fr[fm]; to[1] = fr[md[j^1]]; to[2] = fr[md[j]]; to[3] = fr[fM];
}

#define perm5(x,y,z,u,v) (to[0]=fr[x], to[1]=fr[y], to[2]=fr[z], to[3]=fr[u], to[4]=fr[v])
static void ssort_5x(char ** to, char ** fr, int k) {
	int A = (9*sqless(0,1)) ^ 1;
	int B = (9*sqless(2,3)) ^ 19;
	if (sqless(B&7, A&7)) A ^= B, B ^= A, A ^= B;
	int a = A>>3, b = B>>3; A &= 7; B &= 7;
//	fprintf(stderr,"A:%d B:%d a:%d b:%d\n", A, B, a, b);
	if (sqless(4, B)) {
		if (sqless(B, a)) {
			if (sqless(a, b)) sqless(4, A) ? perm5(4,A,B,a,b) : perm5(A,4,B,a,b);
			else              sqless(4, A) ? perm5(4,A,B,b,a) : perm5(A,4,B,b,a);
		} else {
			if      (sqless(4, A)) perm5(4,A,a,B,b);
			else if (sqless(4, a)) perm5(A,4,a,B,b);
			else 		       perm5(A,a,4,B,b);
		}
	} else {
		if (sqless(4, b)) {
			if (sqless(a, 4)) sqless(a, B) ? perm5(A,a,B,4,b) : perm5(A,B,a,4,b);
			else		  sqless(a, b) ? perm5(A,B,4,a,b) : perm5(A,B,4,b,a);
		} else {
			if (sqless(a, b)) sqless(B, a) ? perm5(A,B,a,b,4) : perm5(A,a,B,b,4);
			else 	          sqless(4, a) ? perm5(A,B,b,4,a) : perm5(A,B,b,a,4);
		}
	}
}
static ssort_t ssort_tab[6] = { 0, 0, ssort_2x, ssort_3x, ssort_4x, ssort_5x };
#define N_SORTFUN 5

static void ssort_2(char ** to, char ** fr, int len, int k) {
	if (len<=N_SORTFUN) return (*ssort_tab[len])(to, fr, k);
	int next[len];
	char *tmpbuf[len], **tmp = tmpbuf;
	int head[96];
	unsigned int map[3]; 
	while(1) {
		int x, i, j;
		map[0] = map[1] = map[2] = 0u; j = 0;
		for (i=0; i<len; i++) {
			int c = (fr[i][k] - 31) & 255; c &= ((c-96) >> 9);
			int ch = c>>5, cl = c&31;
			unsigned int m1 = 1u << cl;
			int nflg = (int) (map[ch] & m1); nflg = ~( (nflg|-nflg) >> 31 );
			next[i] = head[c] | nflg; head[c] = i; map[ch] |= m1; j -= nflg;
		}
		if (j==1) {
			if (map[0] & 1u) { memcpy(to, fr, len*sizeof(char*)); return; }
			++k; continue;
		}
		int ti = 0;
		if (map[0] & 1u) {
			map[0] &= ~1u;
			for (x = head[0]; x>=0; x = next[x]) to[ti++] = fr[x];
		}
		int max_tj = 0, max_ti = 0;
		for (i=0; i<3; i++) {
			int j; unsigned int m; int * ph = head+32*i;
			for (m = map[i]; (j = __builtin_ffs(m)-1) >= 0; m &= ~(1u<<j)) {
				int x = ph[j], tj = 0;
				if (next[x]<0) { to[ti++] = fr[x]; continue; }
				for (; x>=0; x = next[x]) tmp[ti+(tj++)] = fr[x];
				if (tj <= N_SORTFUN) {
					(*ssort_tab[tj])(to+ti, tmp+ti, k+1);
				} else if (tj <= max_tj) {
					ssort_2(to+ti, tmp+ti, tj, k+1);
				} else {
					if (max_tj) ssort_2(to+max_ti, tmp+max_ti, max_tj, k+1);
					max_ti = ti; max_tj = tj;
				}
				ti += tj;
			}
		}
		if (!max_tj) return;
		char **pp = tmp+max_ti; tmp = fr+max_ti; fr = pp; to += max_ti;
		len = max_tj; ++k;
	} 
}

void ssort(char ** pp, int len, int k) {
	if (len<2) return;
	char * tmp[len]; memcpy(tmp, pp, len*sizeof(char*));
	ssort_2(pp, tmp, len, k);
}

///////// number parsing (incl. NaN) ////////////////////////////////////////

void wrnan_seq(void * to, const int * p, int n) {
        long long xl, tl;
        if ((unsigned int)n > 15u) n = 15;
        xl = nan_pkdat[n];
        int bits = xl & 63; xl ^= bits;
        long long msk = (1LL << bits) - 1LL;
        int i, j;
        for (i=j=0; i<n; i++, j+=bits)
                tl = (long long)p[i], xl |= (tl&(1LL<<63)) | ((tl&msk)<<j);
        memcpy(to, &xl, 8);
}

void wrnan_lbl(void * to, const char * s, int n) {
        int i,j;
        long long xl = NAN_LBL;
        if ((unsigned int)n > 6u) n = 6;
        for (i=j=0; i<n; i++, j+=7) xl |= (long long)(s[i]&127) << j;
        xl |= (NAN_SMSK<<n) & NAN_SMSK;
        memcpy(to, &xl, 8);
}

int parse_num(void * to, const char * s) {
        int i,j,k,n;
        char buf[24];
	const char *s0 = s;
        int dat[16], c, emd = 0, dfl = 0;
        while (*s==' ') ++s;
        if (*s==':') for (++s; *s==' '; s++);
        switch(*s) {
		case 0:   *(double*)to = 0.0; return 1;
                case 'n': case 'N':
                          memcpy(to, nan_dflt2, 8);
                          return 1+2*(((s[1]|32)=='a') && ((s[2]|32)=='n'));
                case '-':
                          if ((s[1]|32)=='n') return memcpy(to, nan_dflt2+1, 8),
                                  2+2*(((s[2]|32)=='a') && ((s[3]|32)=='n'));
		case '.': case '+': case 48: case 49: case 50: case 51:
                case 52: case 53: case 54: case 55: case 56: case 57:
                          dfl = ((buf[0] = *s) == '.');
                          for (i=1; i<24; buf[i++]=c, emd+=emd) {
                                  c = s[i];
                                  if (c>47 && c<58) continue;
                                  if (c=='.') { if (dfl++) break; else continue; }
                                  if ((c|33)=='e') { if (emd++) break; else continue; }
                                  if (c=='+' || c=='-') { if (emd!=2) break; else continue; }
				  break;
                          }
                          buf[i] = 0; *(double*)to = atof(buf);
                          return s - s0 + i;
                case '#':
                          for (i=1; i<17; i++) if (!(s[i]&80)) return *(double*)to=1.2345, i;
			  *(double*)to = hx2doub(s+1); return 17;
                case '"':
                          for (i=1; i<7; i++) if (s[i]<32 || s[i]==34 || s[i]>126) break;
                          wrnan_lbl(to, s+1, i-1); return i + (s[i]==34);
                case '[':
                          i=1; j=0;
                          while (1) {
                                  while (s[i]==' ') ++i;
                                  int c = s[i], sg = 0;
                                  if (!c || c==']') { ++i; break; }
                                  if (j==15) break;
                                  if (c=='-') ++i, sg = -1;
                                  else if (c=='+') ++i;
                                  for (n=0; (unsigned int)(k=s[i]-48)<10u; i++) n = 10*n + k;
                                  dat[j++] = (n^sg) - sg;
                          }
                          wrnan_seq(to, dat, j); return i;
                default: memcpy(to, nan_pkdat, 8); return 1;
        }}

const char * dbl2str_s(int ix, double x) {
	static char buf[192];
	char *p = buf + 48 * (ix&3);
	p[dbl2str(p, x)] = 0; return p;
}

///////// 20-bit audio out //////////////////////////////////////////////////

static inline int pml2(int x) { return x>>=7, (x*(0x56000-11*x))>>23; }

static void packhi_2a(char* to, const char *p, const char * q, int n) {   int i;
        if (!q) for (p+=4; n--; to+=4, p+=8) memcpy(to, p, 4);
        else for (i=0,p+=4,q+=4,n>>=1; i<n; i++) memcpy(to+8*i, p+8*i, 4), memcpy(to+8*i+4, q+8*i, 4); }

static int fa_w2k(fa_writer * fa, int fl) {
        int i,j,k, max = 0;
        unsigned short *q16 = fa->blk.lo;
        unsigned char  *q8  = fa->blk.hi;
        unsigned int *p = fa->buf, x, y;
        if (fl) for (i=fa->cur; i<818; i++) p[i] = 0u;
        for (i=0; i<818; i++) { if ((k=p[i]&0x7fffffff)>max) max = k; if (k>0x4fffffff) abort(); }
        int xp = max>>20; if (xp<0x380) xp = 0x380;
        for (i=0; i<818; i++) {
                x = p[i]; y = (x&0x80000000);
                k = xp - ((x>>20)&2047); if (k>20) { p[i] = 0x80000; continue; }
                if (k<0) return log("fa_w2k: inf/nan found"), -1;
                x = (262144u | ((x&0xfffff)>>2)) >> k;
                p[i] = y ? (0x80000-x) : (0x80000+x);
        }
        for (i=j=0; i<409; i++,j+=2) x = p[j], y = p[j+1], q16[j] = x&65535, q16[j+1] = y&65535,
                                     q8[i] = (x>>16)|((y>>12)&0xf0);
        unsigned short mxv = 256u*(xp>0x47f?255u:xp-0x380u) + pml2(max & 0xfffff);
        if (mxv>40000) fprintf(stderr, "wtf???\n"), abort();
        fa->blk.xp16 = fa->toc[k = fa->cnt++ & 1023] = mxv; fa->cur = 0;
        if (write(fa->fd, &fa->blk, 2048)<2048) return -1;
        if (fl) { for (i=k+1; i<1024; i++) fa->toc[i] = 0u; if (fa->nch==1) fa->cnt = ~fa->cnt;
                  if (write(fa->fd, fa->toc, 2048)<2048 || write(fa->fd, &fa->cnt,4)<4) return -1;
                  return close(fa->fd), 1; }
        if (k==1023 && write(fa->fd, fa->toc, 2048)<2048) return -1;
        return 0;
}

const char * au_file_name(const char *dir, int dlen, int id, const char *a1, const char *a2, const char *ext) {
	static char * bptr[2];  static int blen[2];
	int hx[2], j, al1, al2, ix = !memcmp(ext, "a20", 4), elen = ix ? 4 : strlen(ext),
	    siz = dlen + (al1=a1?strlen(a1)+1:0) + (al2=a2?strlen(a2)+1:0) + elen + 16;
	if (blen[ix]<siz) { for (blen[ix] += 8*!blen[ix]; blen[ix] < siz; blen[ix] <<= 1); 
		 	    free(bptr[ix]); bptr[ix] = (char*)malloc(blen[ix]); }
	char * q = bptr[ix]; memcpy(q, dir, dlen); q[dlen]='/';
	hx[0] = qh4((unsigned int)id>>16u) | 0x20202020; 
	hx[1] = qh4(	         id&65535) | 0x20202020; memcpy(q+dlen+1, hx, 8); j = dlen+9;
	if (al1) q[j] = '_', memcpy(q+j+1, a1, al1-1), j += al1;
	if (al2) q[j] = '_', memcpy(q+j+1, a2, al2-1), j += al2;
	q[j] = '.'; memcpy(q+j+1, ext, elen+1);
	return q;
}

const char * au_file_name(int id, int j) {
	static char *srcp = 0, *dstp = 0;
	const char *dir, *ext; int dlen, elen;
	switch(j&6) {
		case 0:
		case 6: if (*(dir=CFG_AO_DIR.s)) dlen = CFG_AO_DIR.i;
			else  dir=QENV('t'), 	 dlen = QENVL('t');     break;
		case 2: dir = ".", dlen = 1; break;
		case 4: dir = QENV('h'), dlen = QENVL('h'); break;
		default: return "WTF";
	}
	if (!j) ext=FA_SUFFIX, elen=FA_SUFFIX_ZLEN; else ext = ".wav\0.flac"+5*(j&1), elen = 5+(j&1);
	char **pp = j ? &dstp : &srcp, *p = *pp = (char*)realloc(*pp, dlen+9+elen);
	memcpy(p, dir, dlen); p[dlen] = '/'; sprintf(p+dlen+1, "%08x", id); memcpy(p+dlen+9, ext, elen);
	return p;
}

int fa_start(fa_writer * fa, int nch) {
	static int id = -1;
        fa->cnt = fa->fd = fa->cur = 0;
        if (id<0) id = time(NULL); else ++id;
	const char * nm = au_file_name(fa->id = id, 0);
	fa->nch = nch; 
	int r = fa->fd = creat(nm, 0600);
	return (r<0) ? (gui_errq_add(EEE_ERRNO), log("autmp: \"%s\": %s", nm, strerror(errno)), EEE_A20) : 0;
}

int fa_add12(fa_writer * fa, const double * x, const double * y, int n) { while (1) {
        int n2 = 818 - fa->cur, n22 = y ? n2>>1 : n2;
        if (n<n2) return n<1 ? 0 : (packhi_2a((char*)(fa->buf+fa->cur), (const char*)x, (const char*)y, n),
                        fa->cur+=n, 0);
        packhi_2a((char*)(fa->buf+fa->cur), (const char*)x, (const char*)y, n2); if (fa_w2k(fa,0)<0) return -1;
        n-=n2; x+=n22; if(y) y+=n22;
}}

int fa_end(fa_writer * fa) { return fa_w2k(fa, 1); }

//////////////////////// scale ///////////////////////////////////////////////////////////////

double Scale01::f0(double v0, double v1, int ty, double x) {
        switch(ty) {
                case  1:case '-': return (1.0-x)*v0 + x*v1;
                case  0:case 'l': return v0 * pow(v1/v0, x);
                case  3:case 'c': v0=cbrt(v0); v1=cbrt(v1); x=(1.0-x)*v0+x*v1; return x*x*x;
                case -3:case 'C': v0=1.0/cbrt(v0); v1=1.0/cbrt(v1); x=(1.0-x)*v0+x*v1; return 1.0/(x*x*x);
                case -1:case 'h': v0=1.0/v0; v1=1.0/v1; x=(1.0-x)*v0+x*v1; return 1.0/x;
                default: break;
        }
        double sg = 1.0;
        if (v0<0.0) v0 = -v0, v1 = -v1, sg = -sg;
        if (v1<0.0) return NAN;
        switch (ty) {
                case  2:case 'q': v0=sqrt(v0); v1=sqrt(v1); x=(1.0-x)*v0+x*v1; return sg*x*x;
                case -2:case 'Q': v0=1.0/sqrt(v0); v1=1.0/sqrt(v1); x=(1.0-x)*v0+x*v1; return sg/(x*x);
                default: return NAN;
        }}

void Scale01::set_ty_2(int ty) {
        if (ty<-3) return;
        if (ty<=3) { m_ty = ty; return; }
        for (int i=0; i<7; i++) {
                if (ty=="CQhl-qc"[i]) { m_ty = i-3; return; }}}

void Scale01::upd_t01() { switch(m_ty) {
        case -3: m_t0 = 1.0/cbrt(m_v0); m_t1 = 1.0/cbrt(m_v1) - m_t0; return;
        case -2: if (m_v0>=0.0&&m_v1>=0.0) m_sg=1.0, m_t0=1.0/sqrt(m_v0), m_t1=1.0/sqrt(m_v1)-m_t0;
                         else m_sg=-1.0, m_t0=1.0/sqrt(-m_v0), m_t1=1.0/sqrt(-m_v1)-m_t0;    return;
        case -1: m_t0 = 1.0/m_v0; m_t1 = 1.0/m_v1 - m_t0; return;
        case  0: m_t1 = log(m_v1/m_v0); return;
        case  1: m_t1 = m_v1 - m_v0; return;
        case  2: if (m_v0>=0.0&&m_v1>=0.0) m_sg=1.0, m_t0=sqrt(m_v0), m_t1=sqrt(m_v1)-m_t0;
                         else m_sg=-1.0, m_t0=sqrt(-m_v0), m_t1=sqrt(-m_v1)-m_t0;   return;
        case  3: m_t0 = cbrt(m_v0); m_t1 = cbrt(m_v1) - m_t0; return;
        default: return;
}}

double Scale01::f(double x) { double z; switch(m_ty) {
                case -3: z = m_t0 + x*m_t1; return 1.0/(z*z*z);
                case -2: z = m_t0 + x*m_t1; return m_sg/(z*z);
                case -1: return 1.0 / (m_t0 + x*m_t1);
                case  0: return m_v0 * exp(x*m_t1);
                case  1: return m_v0 + x*m_t1;
                case  2: z = m_t0 + x*m_t1; return m_sg*z*z;
                case  3: z = m_t0 + x*m_t1; return z*z*z;
                default: return 0.0;
}}

// file util
int is_asv_name(const char *s) { // /.../__asv.lf /.../__asv--x.lf
	int c, l = strlen(s);   const char * as = QENV('a');
	if (!memcmp(s, as, l+1)) return '0';
	return (!memcmp(s,as,l-6) && !memcmp(s+l-6,"--",2) && (c=s[l-4]) && !memcmp(s+l-3,".lf",4)) ? c : 0; }

int coward(const char * fn) {
	int r,fd; struct stat st;
	if ((r = stat(fn, &st))<0) return errno!=ENOENT && (gui_errq_add(EEE_ERRNO, fn), 1);
	if ((st.st_mode&S_IFMT)!=S_IFREG) return 1;
	if ((fd = open(fn, O_RDONLY)) < 0) return gui_errq_add(EEE_ERRNO, fn), 1;
	char buf[32]; memset(buf, 0, 32);
	r = (read(fd,buf,32)>0) && memcmp(buf, "# lflab save file\n_V", 20)
				&& memcmp(buf, "#!/usr/bin/lflab\n_V", 19)
				&& memcmp(buf, "_V0.4\n:TN.$F", 12)
				&& memcmp(buf, "_V0.4\n::F:RC", 12);
	close(fd); return r;
}
