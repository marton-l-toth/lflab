#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <unistd.h>
#include <signal.h>
#include <fcntl.h>
#include <time.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>

#define QWE_UTILC_DEF
#include "uc0.h"

//// auconv /////////////////////////////////////////////////////////////////

volatile int flacpid = 0;

void chld(int s) {
	int k, r = waitpid(flacpid, &k, WNOHANG), i = 1;
	if (r<0) perror("wait(flac)");
	else if (!r||r!=flacpid) fprintf(stderr, "BUG: flacpid=%d, wait ret: %d", flacpid, r);
	else if (WIFSIGNALED(k)) fprintf(stderr,"flac: signal %d\n", WTERMSIG(k));
	else if (WIFEXITED(k)) i=WEXITSTATUS(k), fprintf(stderr,"flac: exit %d (%s)\n",i,i?"error":"ok");
	else fprintf(stderr, "flac: %s\n", WIFSTOPPED(k) ? "stopped" : "what??");
	if (i) exit(i); else flacpid = 0;
}

static const short multab[256] = {
        4075,4063,4051,4040,4028,4016,4005,3993,3982,3971,3959,3948,3937,3925,3914,3903,
        3892,3881,3870,3859,3848,3837,3826,3815,3805,3794,3784,3773,3762,3752,3741,3731,
        3720,3710,3700,3690,3680,3669,3659,3649,3639,3629,3619,3609,3599,3589,3579,3569,
        3560,3550,3540,3531,3521,3511,3502,3492,3483,3473,3464,3454,3445,3436,3426,3417,
        3408,3399,3390,3381,3372,3362,3354,3345,3336,3327,3318,3309,3300,3291,3282,3274,
        3265,3256,3248,3239,3230,3222,3213,3205,3196,3188,3179,3171,3163,3154,3146,3138,
        3129,3121,3113,3105,3097,3089,3080,3072,3064,3056,3048,3040,3032,3024,3016,3008,
        3000,2993,2985,2977,2969,2961,2954,2946,2938,2931,2923,2915,2908,2900,2893,2885,
        2878,2870,2863,2855,2848,2840,2833,2825,2818,2811,2803,2796,2789,2782,2774,2767,
        2760,2753,2746,2739,2731,2725,2717,2710,2703,2696,2689,2682,2675,2668,2661,2654,
        2647,2640,2634,2627,2620,2613,2606,2599,2593,2586,2579,2572,2566,2559,2552,2545,
        2539,2532,2526,2519,2512,2506,2499,2493,2486,2479,2473,2466,2460,2453,2447,2440,
        2434,2428,2421,2415,2408,2402,2395,2389,2383,2376,2370,2364,2357,2351,2345,2338,
        2332,2326,2319,2313,2307,2301,2294,2288,2282,2276,2270,2263,2257,2251,2245,2239,
        2232,2226,2220,2214,2208,2202,2196,2189,2183,2177,2171,2165,2159,2153,2147,2140,
        2134,2128,2122,2116,2110,2104,2098,2092,2086,2080,2074,2067,2061,2055,2049,2043};

typedef struct { unsigned short xp16, lo[818]; unsigned char rsrv, hi[409]; } blk2k;

static int n_blk = -1, n_toc, infd = -1, n_ch = 2;
static const char *inname, *outname;

static unsigned short toc_ent(int i) {
	static unsigned short buf[1024];
	static int cur = -1;
	int j = i>>10; i&=1023;
	if (j!=cur) {
		cur = j; j = (j==n_toc-1) ? n_blk+n_toc-1 : 1025*j+1024;
		if (lseek(infd, j*2048, SEEK_SET)<0 || read(infd, buf, 2048)<2048) perror(inname);
	}
	return buf[i];
}

static void get_blk(blk2k * to, int i) {
	if (lseek(infd, (i+(i>>10))*2048, SEEK_SET)<0 || read(infd, to, 2048)<2048) perror(inname); }

static void wav_head(int fd, int len, int nch) {
        char buf[44];
        //                           scIsiz    aufmt#chan-samprate----byterate----blkal-by/s
        memcpy(buf, "RIFF----WAVEfmt \020\0\0\0\01\0\01\0\x44\xac\0\0\x88\x58\01\0\02\0\020\0data", 40);
        buf[22] = nch;
        int l = 2*len*nch, rate = 88200*nch;
        memcpy(buf+28, &rate, 4); buf[32] = 2*nch; memcpy(buf+40, &l, 4);
        l += 36; memcpy(buf+4, &l, 4);
        write(fd, buf, 44);
}

static int omin=0x7ffffff, omax = -0x7ffffff;
static void wr_blk(int fd, blk2k * b, int maxv) {
	int xp = 15 + (maxv>>8) - (b->xp16>>8), mul = multab[maxv&255],
	   	adj = 1<<(xp-1),  i, j;
	short o[818];
	if (xp>31) { for (i=0; i<818; i++) o[i] = 0; goto wr; }
	for (i=j=0; i<409; i++, j+=2) {
		int x  = (mul*((((int)b->hi[i]& 15)<<16) + (int)b->lo[j  ] - 524288)+adj) >> xp,
		    y  = (mul*((((int)b->hi[i]&240)<<12) + (int)b->lo[j+1] - 524288)+adj) >> xp;
		if (x<omin) omin=x; if (x>omax) omax=x;
		if (y<omin) omin=y; if (y>omax) omax=y;
		o[j] = (short)x; o[j+1] = (short)y;
	}
wr:     if (write(fd, o, 1636)<1636) perror(outname);
}

static blk2k blk;
static const char * usg = "[-r] infile outfile [skip [len]]\n  -r: remove infile on success\n"
"  if outfile ends with \".flac\", flac will be called to encode.\n";
static int a_main(int ac, char ** av) {
	int opt = 0; const char * nm = *av;
	for (; ac>1 && av[1][0]=='-'; ac--,av++) { switch(av[1][1]) {
		case 'r': opt|=1; continue;
		default: fprintf(stderr, "ignoring unknown option \"%s\"\n", av[1]); continue;
	}}
	if (ac<3 || ac>5) fprintf(stderr, "usage: %s %s\n", nm, usg), exit(1);
	infd = open(inname=av[1], O_RDONLY); if (infd<0) perror(av[1]), exit(1);
	int flen = lseek(infd, -4, SEEK_END); if (flen<0) perror(av[1]), exit(1);
	if (flen&2047) fprintf(stderr, "invalid file len %d\n", flen), exit(1);
	if (read(infd, &n_blk, 4)<4) perror("read"), exit(1);
	int fpblk = (n_blk<0) ? (n_blk^=-1, n_ch=1, 818) : 409;
	double secpblk = (double)(fpblk)/44100.0, blkpsec = 1.0 / secpblk;
	n_toc = (n_blk+1023)>>10;
	int skip = ac<4 ? 0 : (int)floor(blkpsec*atof(av[3]));
	int len  = ac<5 ? 0x7fffffff : (int)ceil(blkpsec*atof(av[4]));
	printf("nblk=%d, exp.len=%d, file len=%d\n", n_blk, 2048*(n_blk+n_toc),flen);
	if (skip>=n_blk) fprintf(stderr, "skipping all...\n"), exit(1);
	if (len>n_blk-skip) len = n_blk - skip;
	int i, k, max = 0;
	for (i=0; i<len; i++) if ((k=toc_ent(skip+i))>max) max = k;
	if (!max) fprintf(stderr, "it's all zeroes: %s\n", (!(opt&1) || ac>3) ? "nothing done" :
			   (unlink(inname)<0 ? strerror(errno) : "deleted") ), exit(90);
	int outfd, onl = strlen(outname=av[2]), flacmode = (onl>4 && !memcmp(outname+onl-5, ".flac", 5));
	if (flacmode) {
		char isz[24]; memcpy(isz, "--input-size=", 13); isz[13+sprintf(isz+13,"%d",len*1636)] = 0;
		if ((flacpid = launch("flac", "!>..", &outfd, "-o", outname, "--endian=l",
				      n_ch==1 ? "--channels=1" : "--channels=2",
				      "--bps=16", "--sign=s", "--sample-rate=44100", isz, "-", NULL)) < 0)
			return perror("flac"), 1;
		signal(SIGCHLD, chld);
	} else {
		outfd = creat(outname=av[2], 0600);  if (outfd<0) perror(av[2]), exit(1);
		wav_head(outfd, len*fpblk, n_ch);
	}
	for (i=0; i<len; i++) get_blk(&blk, i+skip), wr_blk(outfd, &blk, max);
	close(outfd); while (flacpid) sleep(60);
	fprintf(stderr, "written %d blocks (%f sec), min:%d, max:%d\n", len, 
			(double)len*409.0/44100.0, omin, omax);
	return (opt&1) ? (unlink(inname)<0 && (perror(inname),1)) : 0;
}

//// console ////////////////////////////////////////////////////////////////

int oops(const char *s) { if (*s==33) fprintf(stderr,"con/fatal: %s\n", s+1); else perror(s); 
			  sleep(5); return 1; }

int c_main() {
	const char * tdir = getenv("LF_TMPDIR"); if (!tdir) return oops("!LF_TMPDIR undefined");
	int i, infd, outfd, l = 0; while (tdir[l]) ++l; 
	char *p = malloc(l+3); for (i=0; i<l; i++) p[i] = tdir[i]; p[i++] = '/'; p[i+1] = 0;
	if ((p[i] = 'C', outfd = open(p, O_WRONLY)) < 0 ||
	    (p[i] = 'm',  infd = open(p, O_RDWR)  ) < 0    ) return oops(p);
	fd_set rset; struct timeval tv;
	char buf[4096];
	int reop_cnt = 0, c0buf[3];
	c0buf[0] = 0x635f0000; c0buf[1] = qh4(getpid()); c0buf[2] = 10; write(outfd, (char*)c0buf + 2, 7);
	while(1) {
		FD_ZERO(&rset); FD_SET(0, &rset); FD_SET(infd, &rset);
		tv.tv_sec=1; tv.tv_usec = 0;
		int r2, r = select(infd+1, &rset, 0, 0, &tv);
		if (r<0) return oops("select");
		if (!r) { if (reop_cnt) --reop_cnt; continue; }
		if (FD_ISSET(0, &rset)) {
			r = read(0, buf, 4096);
			if (r<=0) return oops(r?"stdin":"!EOF");
			if ((r2=write(outfd, buf, r)) !=r ) return oops(r2?"'C'":"!con: cmd write failed");
		}
		if (FD_ISSET(infd, &rset)) {
			if ((r = read(infd, buf, 4096)) > 0) { write(1, buf, r); continue; }
			if (++reop_cnt>5 || (close(infd), r=infd=open(p, O_RDWR))<0)
				return oops(r?p:"!reopen failed 5x");
			fprintf(stderr,"reopened %dx\n", reop_cnt);
		}
	}}

/**************** log i/o ****************************************************/

static int con_stat = 0, cxt_pid = 0, msg_fd = -1, errtemp = 0, killer_fd = 0, der_alte = -1;

void ts_9(char *to) {
        static time_t sec = 0; static char sbuf[8] = {48,48,48,48,48,46,0,0};
        struct timeval tv; gettimeofday(&tv, NULL);
        time_t sec2 = tv.tv_sec; if (sec2!=sec) {
                sec = sec2;
                struct tm loct; localtime_r(&sec, &loct);
                sbuf[0] = "0123456789ten123456789TE"[loct.tm_hour];
                d59(sbuf+1, loct.tm_min); d59(sbuf+3, loct.tm_sec);
        }
	memcpy(to, sbuf, 6); d999(to+6, tv.tv_usec/1000);
}

/*** output buf ***************/
static char logbuf[16384];
static int log_ix = 0, log_ox = 0, log_flg = 0;

static void LOG(const char * fmt, ...), kill_con();
static void log_sn(const char *s, int n, int conpfx) {
	int r = 0, ix0 = log_ix, n0 = n;
	if (n<1) { switch(n) {
		case -2: if (log_flg&1) write(msg_fd, logbuf+ix0, 16384-ix0);
			 if ((r = write(msg_fd, logbuf, ix0))<0) goto err; else return;
		case -1: if ((n=ix0 - log_ox)>0 &&
			      (r = write(2, logbuf+log_ox, n), log_ox=ix0, r<1)) goto lerr;
		case 0: return;
	}}
	if ((log_ix+=n)>=16384 && (log_flg|=1, log_ix&=16383)) memcpy(logbuf, s+(n=16384-ix0), log_ix);
	if (n) memcpy(logbuf + ix0, s, n);
	while ((log_ix^log_ox)&~4095) {
		int op2 = (log_ox&~4095) + 4096;
		    r   = write(2, logbuf + log_ox, op2 - log_ox);
		if (r<=0) { log_ox = log_ix; goto lerr; }
		log_ox = op2 & 16383;
	}
	if (msg_fd<0 || (n = n0 + conpfx)<1) return; else ix0 = log_ix - n;
	if (ix0<0 && (r=write(msg_fd, logbuf+(ix0&16383), -ix0), n+=ix0, ix0=0, r<1)) goto err;
	if ((r=write(msg_fd, logbuf+ix0, n))>0) return;
err:	kill_con(); return LOG("[console write err: %s]\n", r==-1 ? strerror(errno) : "zero ret");
lerr:   if (!con_stat) return;
	char errbuf[256]; n = snprintf(errbuf, 255, "log file write error: %d(%s)\n",
			r ? errno : 0, r ? strerror(errno) : "zero ret.");
	write(msg_fd, errbuf, n);
}

static void log_l(char *s, int n, int pfx) {
	char buf[12]; ts_9(buf); buf[9] = pfx; buf[10] = 32;
	log_sn(buf, 11, INT_MIN); log_sn(s, n+1, *s=='`'?INT_MIN:11); 
}

static void LOG(const char * fmt, ...) {
	char buf[1024]; ts_9(buf); buf[9] = 'i'; buf[10] = 32;
	va_list ap; va_start(ap, fmt);
	int r = vsnprintf(buf+11, 1000, fmt, ap); va_end(ap);
	if (r<1) memcpy(buf+11, "(??"")", r=4);
	buf[r+11] = 10; log_sn(buf, r+12, buf[11]=='`' ? INT_MIN : 0);
}

static void bye(int kf) {
	kill_con(); if (!kf) exit(1); char * exe = getenv("LF_CLEANUP");
	LOG("ioprc done, killer_fd = %d, cleanup script: \"%s\"", killer_fd, exe ? exe : "");
	log_sn(0,-1,0); if (killer_fd) close(killer_fd);
	if (exe) close(1), dup2(2,1), close(0), open("/dev/null", O_RDONLY),
	         execlp(exe, exe, (char*)0), perror("exec failed");
	exit(1);
}

/*** cmd **********************/

static void io_chld(int _) {
	int st, pid = waitpid(-1, &st, WNOHANG); if (pid<=0 || pid!=cxt_pid) goto done;
	msg_fd = -1; con_stat = cxt_pid = 0; write(1, "_c-2\n", 5);
done:	signal(SIGCHLD, io_chld);
}

static int start_con() {
	static const char * sh = 0;
	return -(cxt_pid || (!sh && !(sh = getenv("LF_CON"))) ||
		 (cxt_pid = launch("xterm", "!)x1", "-sb", "-e", sh, (char*)0))<0 ||
		 (msg_fd = open(tpipe_name('m'), O_WRONLY|O_APPEND))<0); }

static void kill_con() { if (cxt_pid)  kill(cxt_pid, 9); }

static void cmd_l(char *s, int n, int src) { s[n] = 0; switch(*s) {
	case 'c': if (start_con()<0) LOG("start_con failed");    return;
	case 'C': return con_stat = 1, log_sn(0, -2, 0);
	case 'f': return log_sn(0, -1, 0);
	default: LOG("unknown cmd%c 0x%x (%s)", 48+src, *s, s); return;
}}

/*** input buf ****************/
typedef void (*lfun)(char *, int, int);
typedef struct { int fd, siz, cont, arg; lfun lf; char * p; } ibuf;

#define N_IBUF 7
char bigbuf[0x7000];
static ibuf ib[7] = {{ 0, 4096, 0, 0,   &cmd_l, bigbuf},
		     {-1, 4096, 0, '.', &log_l, bigbuf+0x1000},
		     {-1, 4096, 0, 'p', &log_l, bigbuf+0x2000},
		     {-1, 4096, 0, 't', &log_l, bigbuf+0x3000},
		     {-1, 4096, 0, 'u', &log_l, bigbuf+0x4000},
		     {-1, 4096, 0, 'x', &log_l, bigbuf+0x5000},
		     {-1, 4096, 0, 'k', &log_l, bigbuf+0x6000}};

static void ib_init(int mlog_fd) {
	if ((ib[1].fd = mlog_fd)&0xffff8000) LOG("mlog_fd = %d", mlog_fd); 
	int i; for (i=2; i<N_IBUF; i++) if ((ib[i].fd = open(tpipe_name(ib[i].arg), O_RDWR))<0)
		LOG("cannot open log pipe '%c': %s", ib[i].arg, strerror(errno));
}

static void ib_read(ibuf * b) {
	if (b->fd<0) return;
	char *q, *ql, *p = b->p;
	int ag = b->arg, sz = b->siz, r = read(b->fd, p + b->cont, sz - b->cont);
	if (r<1) { LOG("inpipe(%c): %s", ag?ag:48, r?strerror(errno):"zero returned");
		   if (ag<65) { LOG("unnamed pipe fail, closing..."); b->fd = -1; return; }
		   LOG("logpipe fail, errtemp=%d", errtemp);
		   if (errtemp<99) errtemp+=10, close(b->fd),
			           b->fd = open(tpipe_name(ag), O_RDWR);
		   return; }
	for (q=p+b->cont, ql=q+r; q<ql; q++) if (*q==10) (*b->lf)(p,q-p,ag), p = q+1;
	if (p==q) { b->cont = 0; return; }
	if (4*(b->cont=q-p)>3*sz) LOG("line len>%d, thrown away...", b->cont), b->cont = 0;
	else memmove(b->p, p, b->cont);
}

/*** main *****/

static int orph = 0;
static void tick(int k) {
	static int cnt = 0; if ((cnt+=k)<10) return; else cnt = 0;
	int pp = getppid(), flg = der_alte<0 ? (pp==1) : (pp!=der_alte);
	if (orph ? ++orph : (orph = flg)) {
		LOG("ioprc orphaned (%d!=%d), %c/3", pp, der_alte, 48+orph); if (orph>2) bye(1); }}

static void hello() {
	char buf[1024]; int k = (int)readlink("/proc/self/fd/2", buf, 1023);
	if (k<0) LOG("ioprc started, /proc/self/fd/2: %s", strerror(errno));
	else buf[k] = 0, LOG("ioprc started, logfile is \"%s\"", buf);
}

#define FOR_IB for(i=0; i<N_IBUF; i++) if ((k=ib[i].fd)>=0)
int i_main(int ac, char ** av) {
	signal(SIGCHLD, io_chld);
	hello();
	ib_init(   ac<2 ? -1 : qh4r(*(int*)av[1]));
	killer_fd = (ac<3||*av[2]<48) ? 0 : qh4r(*(int*)av[2]);
	der_alte = ac<4 ? -1 : qh4r(*(int*)av[3]);
	int i, k;
	fd_set rset; struct timeval tv;
	while(1) {
		FD_ZERO(&rset);
		int maxfd = 0;
		FOR_IB { FD_SET(k, &rset); if (k>maxfd) maxfd = k; }
		tv.tv_sec=0; tv.tv_usec = 250000;
		int r = select(maxfd+1, &rset, 0, 0, &tv);
		if (r<0) { LOG("select(): %s", strerror(errno)); r = 0; }
		tick(r?r:10);
		if(!r){ if ((errtemp -= 20)<0) errtemp = 0;  continue; }
		if ((errtemp -= r)<0) errtemp = 0;
		FOR_IB if (FD_ISSET(k, &rset)) ib_read(ib+i);
	}}

/////////////////////////////////

int usage(const char *s) { return fprintf(stderr, "lf.bb: name \"%s\"%s", s,
		" unknown, valid names are lf.acv, lf.io and lf.con\n"), 1; }

int main(int ac, char ** av) {
	const char *q, *s = *av;
	for (q=s; *s; s++) if (*s=='/') q = s+1;
	if (memcmp(q, "lf.", 3)) return usage(q);
	switch(q[3]) {
		case 'c': return c_main();
		case 'a': return a_main(ac, av);
		case 'i': return i_main(ac, av);
		case 'l': return execlp("less", "less", getenv("LF_LIC"), NULL);
		default: return usage(q); }}

