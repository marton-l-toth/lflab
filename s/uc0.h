#ifndef __qwe_uc0_h__
#define __qwe_uc0_h__

#define FIFO_LIST "ptuxmskACT"
#define N_XAPP 3

static inline int    min_i(int    x, int    y) { return x<y ? x : y; }
static inline double min_d(double x, double y) { return x<y ? x : y; }
static inline int    max_i(int    x, int    y) { return x>y ? x : y; }
static inline double max_d(double x, double y) { return x>y ? x : y; }
static inline int ivlim(int x, int m, int M) { return x<m ? m : (x>M ? M : x); }
static inline void d59(char *p, int v) {int t = (v*13)>>7;   *(short*)p = (short)(t + ((v-10*t)<<8) + 0x3030);}
static inline void d99(char *p, int v) {int t = (v*205)>>11; *(short*)p = (short)(t + ((v-10*t)<<8) + 0x3030);}
static inline int hxd2i(int c) { return 15&(c+9*(c>>6)); }
static inline int hexc1(int x) { return x + 48 + (((9-x)>>4)&7); }
static inline int atoi_h(const char *s) { int r = 0, x = -(*s=='-');
	for (s-=x; *s&80; s++) r = 16*r+hxd2i(*s); return (r^x)-x; }
#define BVFOR_JMC(X) unsigned int j, m; for (m = (X); j=__builtin_ffs(m), j-- > 0; m &= ~(1u<<j))
#define qh4rs(P) qh4r(*(const int*)(P))
extern const unsigned char bitrev8[256];

#ifndef QWE_UTILC_DEF
extern int map_errno;
extern volatile char vstring[16];
extern char *xapp_env[N_XAPP], *xapp_dflt[N_XAPP];
extern const char **xapp_ls[N_XAPP];
void vstring_set(int i, int j);
void d999(char *p, int v);
int qh4(unsigned int x), qh4r(unsigned int x);
void h5f(char *to, int x);
const char * tpipe_name(int c);
int launch(const char * exe, const char * iocfg, ...);
void * map_wdir_shm(int c, size_t siz, int wrf); // 1:wr 2:cre
unsigned int * tlog_cp(const char *h16);
void set_fd(int *to, int from, int keep);
#else

volatile char vstring[16];
int map_errno = 0;
const char *xapp_env[N_XAPP] = { "LF_XTERM", "LF_X_ED_T", "LF_X_ED_X" };
const char *xapp_dflt[N_XAPP] = { "xterm", "nano", "xedit" };
static const char *xterm_ls[] = { "x-terminal-emulator", "xterm", "gnome-terminal", "konsole", "lxterminal",
	        		  "uxterm", "terminator", "rxvt", 0 },
		  *ed_t_ls [] = { "nano", "vim", "vi", "emacs", 0 },
		  *ed_x_ls [] = { "nedit", "xedit", "gvim", "gedit", "xemacs", 0 };
const char **xapp_ls[N_XAPP] = {xterm_ls, ed_t_ls, ed_x_ls};
const unsigned char bitrev8[256] = {
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

void vstring_set(int i, int j) { char buf[16]; j = snprintf(buf, 15, "vvrrssnn%d.%d", i, j); 
				for (i=0; i<j; i++) vstring[i] = buf[i]; }

int qh4(unsigned int x) {
        x = (x&0xff00) | ((x&255)<<24);
        x = (x&0xf000f00) | ((x>>12)&0xf000f);
        unsigned int ten = (x+0x6060606)&0x10101010;
        return (int)(x + 0x30303030 + (ten>>1) - (ten>>4));
}

int qh4r(unsigned int x) {
        unsigned int ten = x & 0x40404040;
        x = (x&0xf0f0f0f) + 9*(ten>>6);
        x = (x&0xf000f00) | ((x&0xf000f)<<12);
        return (int)((x&0xff00) | ((x>>24)&255));
}

void h5f(char *to, int x) {
	to[0] = hexc1(x>>16);
        x = (x&0xff00) | ((x&255)<<24);
        x = (x&0xf000f00) | ((x>>12)&0xf000f);
        unsigned int ten = (x+0x6060606)&0x10101010;
        *(int*)(to+1) = (int)(x + 0x30303030 + (ten>>1) - (ten>>4));
}

void d999(char *p, int v) { int ht = (v*205)>>11, h = (ht*205)>>11;
	*(short*)p = (short)(h + ((ht-10*h)<<8) + 0x3030), p[2] = v - 10*ht + 48; }

const char * tpipe_name(int c) {
        static char path[256];  static int ix = -1;
        if (ix<0) {
                const char * s =
#ifdef QENV
			QENV('t');
#else
			getenv("LF_TMPDIR");
#endif
                if (!s || !*s || (ix=strlen(s))>248) s="/tmp", ix = 4;
                memcpy(path, s, ix); memcpy(path+ix, "/_\0", 4); ix++;
        }
        path[ix] = c; return path;
}

void * map_wdir_shm(int c, size_t siz, int wrf) {
        static const int md[6] = { O_RDONLY,O_RDWR,O_RDONLY,O_RDWR|O_CREAT, PROT_READ,PROT_READ|PROT_WRITE };
        const char * fnm = tpipe_name(c); void *r = 0;
	int fd = open(fnm, md[wrf&3], 0600); if (fd<0) return  map_errno=errno, r;
	if ((wrf&2) && ftruncate(fd, siz)<0) return close(fd), map_errno=errno, r;
	if ((r = mmap(NULL, siz, md[4+(wrf&1)], MAP_SHARED, fd, 0)) == MAP_FAILED) map_errno=errno, r=0;
	close(fd); return r;
}

unsigned int * tlog_cp(const char *h16) {
	static unsigned int bits = 99, ec = 0, *buf = 0;
	if (bits>98)  { if (bits>99) return &bits;
			const char * s = getenv("LF_TLOG_BITS"); if (!s) return bits=0xffff0001, &bits;
			if (!(buf=(unsigned int *)map_wdir_shm('@', 8<<(bits=*s&31), 0))) 
				return bits=0xffff0002, &bits; }
	int i, h4[4]; for (i=0; i<4; i++) h4[i] = qh4rs(h16+4*i);
	int siz = 2<<bits, i0 = (h4[0]<<16) + h4[1], n = (h4[2]<<16) + h4[3];
	if ((i0|n) >= siz) return ec=0xffff0003, &ec;
	unsigned int *r = (unsigned int*)malloc(4*n+8); r[0] = n;
	return (i0+n<=siz) ?  memcpy(r+2, buf+i0, 4*n) 
			   : (memcpy(r+2, buf+i0, 4*(siz-i0)), memcpy(r+2+siz-i0, buf, 4*(n-siz+i0))), r;
}

void set_fd(int *to, int from, int keep) { if (*to!=from) *to<0 ? (*to=from)
	        : (close(*to), from<0 ? (*to=from) : (dup2(from, *to), (keep||close(from)))); }

// (iocfg-a*: (|)|. :none <|>:int* -|*|=|+->const char*) a1 a2 ... (char*)NULL 
// null(r), null(w), keep   ro,cre,app,rw,wr
// len(iocfg)>3 -> a1 is ("%04X",fd)*

#define LAUNCH_DBG //fprintf(stderr, "launch:%s:%d : vararg for '%c'\n", exe, i, k)
int launch(const char * exe, const char * iocfg, ...) {
        const char *av[256], *s;
        int k, ff, r, ac = 1, i, nf, fds[64], a0buf[64], a0cp = iocfg[0]=='!' && ++iocfg;
        for (i=0; i<64; i++) a0buf[i] = 0x45454545;
        va_list ap; va_start(ap, iocfg);
        for (i=0; i<60 && (k=iocfg[i]); i++) { switch(k) {
                case '<': if (pipe(fds+i)<0 || (ff=fcntl(fds[i],F_GETFD))<0
					    ||     fcntl(fds[i],F_SETFD,ff|FD_CLOEXEC)<0) return -1; LAUNCH_DBG;
                          *(va_arg(ap, int*)) = fds[i]; fds[i] = fds[i+1]; continue;
                case '>': if (pipe(fds+i)<0 || (ff=fcntl(fds[i+1],F_GETFD))<0
					    ||     fcntl(fds[i+1],F_SETFD,ff|FD_CLOEXEC)<0) return -1; LAUNCH_DBG;
                          *(va_arg(ap, int*)) = fds[i+1]; continue;
                case '-': case '+': case '*': case '=': LAUNCH_DBG;
                          av[i] = va_arg(ap, const char*); if (!av[i]) av[i] = "/dev/null"; continue;
                default: continue;
        }}
        av[nf=i] = a0cp ? exe : va_arg(ap, const char*);
        if (nf>3) av[nf+ac++] = (const char*)a0buf, a0buf[nf-3] = 0;
        while ((av[nf+ac] = va_arg(ap, const char*))) ++ac;
        va_end(ap); r = fork(); if (r) { // parent
		for (i=0; i<60 && (k=iocfg[i]); i++) if ((k|2)=='>') close(fds[i]);
		return r; }
        int md; //child
        for (i=0; i<nf; i++) {
                if ((unsigned int)(r=(k=iocfg[i])-48)<10u && r<i) {
                        if (i<3) close(i), dup2(r,i); else a0buf[i-3] = r<3 ? qh4(r) : a0buf[r-3];
                        continue; }
                switch(k) {
                        case '.': continue;
                        case '(': k = -1; md = O_RDONLY; goto fdchk;
                        case ')': k = -1; md = O_WRONLY|O_APPEND; goto fdchk;
                        case '<': k = fds[i]; md = O_RDONLY; goto fdchk;
                        case '>': k = fds[i]; md = O_WRONLY|O_APPEND; goto fdchk;
                        case '-': md = O_RDONLY; goto op;
                        case '+': md = O_WRONLY|O_APPEND|O_CREAT; goto op;
                        case '=': md = O_RDWR; goto op;
                        case '*': k = creat(av[i], 0600); md = O_WRONLY|O_APPEND; goto fdchk;
                        default:  md = O_WRONLY|O_APPEND|O_CREAT; av[i] = tpipe_name(k); goto op;
                }
op:             k = open(av[i], md, 0600);
fdchk:          if (k<0) k = open("/dev/null", md);
                if (k<0) continue; // really impossible...
                if (i<3 && k!=i) close(i), dup2(k, i), close(k);
                else a0buf[i-3] = qh4(k);
        }
        char abuf[16384], *av2[256], *q = abuf; av2[ac] = 0;
        for (i=0; i<ac; i++) {
                for (av2[i]=q, s=av[nf+i]; *s; s++,q++) *q = *s;    *(q++) = 0; }
        execvp(exe, av2);
        fprintf(stderr,"execvp \"%s\" failed: %s\n", exe,  strerror(errno)); fflush(stderr);
	raise(9); while(1) sleep(60); // call your <censored>, not the static destructors
}

#endif // QWE_UTILC_DEF
#endif // __qwe_uc0_h__
