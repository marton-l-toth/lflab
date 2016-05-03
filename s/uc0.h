#ifndef __qwe_uc0_h__
#define __qwe_uc0_h__

#define FIFO_LIST "ptuxmskACT"

static inline void d59(char *p, int v) {int t = (v*13)>>7;   *(short*)p = (short)(t + ((v-10*t)<<8) + 0x3030);}
static inline void d99(char *p, int v) {int t = (v*205)>>11; *(short*)p = (short)(t + ((v-10*t)<<8) + 0x3030);}
static inline int hxd2i(int c) { return 15&(c+9*(c>>6)); }

#define qh4rs(P) qh4r(*(const int*)(P))

#ifndef QWE_UTILC_DEF
extern volatile char vstring[16];
void vstring_set(int i, int j);
void d999(char *p, int v);
int qh4(unsigned int x), qh4r(unsigned int x);
const char * tpipe_name(int c);
int launch(const char * exe, const char * iocfg, ...);
#else

volatile char vstring[16];
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

void d999(char *p, int v) { int ht = (v*205)>>11, h = (ht*205)>>11;
	*(short*)p = (short)(h + ((ht-10*h)<<8) + 0x3030), p[2] = v - 10*ht + 48; }

const char * tpipe_name(int c) {
        static char path[256];  static int ix = -1;
        if (ix<0) {
                const char * s = getenv("LF_TMPDIR");
                if (!s || !*s || (ix=strlen(s))>248) s="/tmp", ix = 4;
                memcpy(path, s, ix); memcpy(path+ix, "/_\0", 4); ix++;
        }
        path[ix] = c; return path;
}

// (iocfg-a*: (|)|. :none <|>:int* -|*|=|+->const char*) a1 a2 ... (char*)NULL 
// null(r), null(w), keep   ro,cre,app,rw,wr
// len(iocfg)>3 -> a1 is ("%04X",fd)*

int launch(const char * exe, const char * iocfg, ...) {
        const char *av[256], *s;
        int k, ff, r, ac = 1, i, nf, fds[64], a0buf[64], a0cp = iocfg[0]=='!' && ++iocfg;
        for (i=0; i<64; i++) a0buf[i] = 0x45454545;
        va_list ap; va_start(ap, iocfg);
        for (i=0; i<60 && (k=iocfg[i]); i++) { switch(k) {
                case '<': if (pipe(fds+i)<0 || (ff=fcntl(fds[i],F_GETFD))<0
					    ||     fcntl(fds[i],F_SETFD,ff|FD_CLOEXEC)<0) return -1;
                          *(va_arg(ap, int*)) = fds[i]; fds[i] = fds[i+1]; continue;
                case '>': if (pipe(fds+i)<0 || (ff=fcntl(fds[i+1],F_GETFD))<0
					    ||     fcntl(fds[i+1],F_SETFD,ff|FD_CLOEXEC)<0) return -1; 
                          *(va_arg(ap, int*)) = fds[i+1]; continue;
                case '-': case '+': case '*': case '=':
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
                        case ')': k = -1; md = O_RDONLY; goto fdchk;
                        case '(': k = -1; md = O_WRONLY|O_APPEND; goto fdchk;
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
