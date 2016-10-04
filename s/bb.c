#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <dirent.h>
#include <unistd.h>
#include <signal.h>
#include <fcntl.h>
#include <time.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <sys/select.h>
#include <sys/inotify.h>
#include <sys/mman.h>


#define QWE_UTILC_DEF
#include "uc0.h"

/*** input buf (io, wrk) ******/
typedef void (*lfun)(char *, int, int);
typedef struct { int fd, siz, cont, arg; lfun lf; char * p; } ibuf;

static int ib_read(ibuf * b) { // -1:unixerr -2:zero_ret -3:frag_discard -4:no_fd
	if (b->fd<0) return -4;
	char *q, *ql, *p = b->p;
	int ag = b->arg, sz = b->siz, r = read(b->fd, p + b->cont, sz - b->cont); if (r<1) return -1-!r;
	for (q=p+b->cont, ql=q+r; q<ql; q++) if (*q==10) (*b->lf)(p,q-p,ag), p = q+1;
	int k = b->cont = q-p; 
	return k ? ( (4*k>3*sz) ? (b->cont=0, -3) : (memmove(b->p, p, b->cont), 0) ) : 0;
}

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

int c_main() { /* console */
	const char * cnm = tpipe_name('C');
	FILE * f = fopen(cnm, "a"); if (!f) return perror(cnm), sleep(3), 1;
	fprintf(f, "_c/proc/%d/fd/1\n", getpid()); fclose(f); while(1) sleep(60);  }

int e_main() { /* editor */
	const char * cnm = tpipe_name('C'), *jt = getenv("LF_ED_TMPF"), *ed = getenv(xapp_env[1]);
	FILE *f; int x, r = fork(); switch(r) {
		case 0: execlp(ed, ed, jt+1, NULL); perror("exec"); exit(1);
		case -1: perror("fork"), exit(1);
		default: signal(SIGHUP, SIG_IGN); waitpid(r,&x,0); 
			 exit((f=fopen(cnm, "a")) ? (fprintf(f,"-%c\n",*jt),fclose(f),0) 
			 			  : (perror(cnm),1)); }}

/**************** log i/o ****************************************************/
#define N_IBUF 9
#define LOG LOG_io
static void LOG_io(const char * fmt, ...), bye(int kf), log_l(char*,int,int), cmd_l(char*,int,int);
volatile int chld_flg = 0;
static int con_stat = 0, msg_fd = -1, errtemp = 0, killer_fd = 0, der_alte = -1, wrkdir_l = 0,
	   shutdn_flg = 8, cl_sec = 0, cl_usec = 0, selwt = 250000, gui_stat = 0, rs_ts = 0, dflg = 0;
static const char *xapp_name[N_XAPP], *wrkdir, *lf_ed;
static unsigned int ino_bv = 0u;
static int ino_epid[32], ino_oid[32], ino_fd = -1;
static char ino_inm[80], logbuf[16384], bigbuf[0x8000], xapp_nbuf[N_XAPP*256];
static int log_ix = 0, log_ox = 0, log_flg = 0;
static ibuf ib[9] = {{ 0, 3072, 0, 0,   &cmd_l, bigbuf},
		     {-1, 4096, 0, '.', &log_l, bigbuf+0x0c00},
		     {-1, 1024, 0, 'C', &cmd_l, bigbuf+0x1c00},
		     {-1, 4096, 0, 'p', &log_l, bigbuf+0x2000},
		     {-1, 4096, 0, 't', &log_l, bigbuf+0x3000},
		     {-1, 4096, 0, 'u', &log_l, bigbuf+0x4000},
		     {-1, 4096, 0, 'x', &log_l, bigbuf+0x5000},
		     {-1, 4096, 0, 'k', &log_l, bigbuf+0x6000},
		     {-1, 4096, 0, 's', &log_l, bigbuf+0x7000}};

static void ts_9(char *to) {
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

/*** shutdown --- shutdn_flg: 1-ppid 2-cmdin 4-.log 8-ulog */

static int t_usec(int set) {
	struct timeval tv; gettimeofday(&tv, NULL); if (set) return cl_sec = tv.tv_sec, cl_usec = tv.tv_usec;
	int k = tv.tv_sec - cl_sec; return (k>2100) ? INT_MAX : 1000000*k + (tv.tv_usec - cl_usec); }

static void shutdn(int ev) {
	int oldf = shutdn_flg, st = -1, us = 0;
	if (ev<0) shutdn_flg &= ev; else shutdn_flg |= ev;
	if (!(shutdn_flg&1) && (getppid()!=der_alte)) shutdn_flg |= 1;
	if (shutdn_flg & 3) {
		if (!(oldf&3)) selwt=101000, t_usec(1);
		else us = t_usec(0);
		if (!(~shutdn_flg & 14) || ((shutdn_flg&1) && us > 500000)) st = 1;
		else if (us > 1000000) st = 0;
	}
	if (shutdn_flg!=oldf || st>=0) LOG("shutdown (%d, %d)", shutdn_flg, us);
	if (st>=0) bye(st);
}

static void del_n_log(const char *s, int f, int e) {
	if ((f?rmdir:unlink)(s)<0 && errno != e) LOG("del %s: %s",s,strerror(errno));}
static void rm_r(const char * dn) {
	DIR * dir = opendir(dn); if (!dir) return LOG("rm_r: opendir: %s\n", strerror(errno));
	struct dirent * ent;  
	int l0 = strlen(dn); char buf[l0+1026]; memcpy(buf, dn, l0); buf[l0++] = '/'; buf[1024]=buf[1025]=0;
	while ((ent = readdir(dir))) {
		const char * s = ent->d_name; 
		if (*s=='.' && (!s[1] || (s[1]=='.'&&!s[2]))) continue;
		strncpy(buf+l0, s, 1024); del_n_log(buf, 0, 0);
	}
	buf[l0] = 0; del_n_log(buf,1,0);
}

static void search_and_destroy(int kfd) {
        DIR * dir = opendir("/proc");
        struct dirent * ent;
        if (!dir) { perror("/proc"); return; }
        char buf[256]; memcpy(buf, "/proc/", 6);
        char kfb[64]; int kfb_l = snprintf(kfb, 63, "/fd/%d", kfd); kfb[kfb_l++] = 0;
        const char *exe, *kf = getenv("LF_KILLER"); if (!kf) puts("???"), exit(1);
        int pid, i, k, ec, kfl = strlen(kf)+1, mypid = getpid();
        if (kfl>1023) return (void) puts("wtf");
        char buf2[1024];
        while ((ent = readdir(dir))) {
                const char * s = ent->d_name;
                for (i=pid=0; i<8 && (unsigned int)(k=s[i]-48)<10; i++) buf[6+i] = s[i], pid = 10*pid + k;
                if ( k!=-48 || pid==mypid || (memcpy(buf+6+i,kfb,kfb_l), (ec=readlink(buf,buf2,kfl+1))<0) ||
                    (buf2[ec]=0, memcmp(buf2, kf, kfl)) ) continue;
                memcpy(buf+6+i, "/exe", 5); ec = readlink(buf,buf2,1023);
                exe = (ec<0) ? "(readlink(exe) failed)" : (buf2[ec]=0, buf2);
		if (!memcmp(exe, "/usr/bin/gnome-terminal", 24)) { ec = '='; goto msg; }
                switch(kill(pid, 9)) {  case 0:     ec = '+'; break;
                                        case ESRCH: ec = '-'; break;
                                        case EPERM: ec = '!'; break;
                                        default:    ec = '?'; break; }
msg:		LOG("kill %5s %c %s", s, ec, exe);
        }
        closedir(dir);
}

static void del_workdir() {
        const char *s, *s0 = wrkdir;  int l0 = wrkdir_l;
	ino_inm[l0+4] = ino_inm[l0+44] = 0; rm_r(ino_inm+1); rm_r(ino_inm+41);
        char buf[l0+16]; memcpy(buf, s0, l0); buf[l0] = '/'; buf[l0+2] = 0;
        for (s=FIFO_LIST; *s; s++) buf[l0+1] = *s, del_n_log(buf, 0, 0);
        for (s="+,"     ; *s; s++) buf[l0+1] = *s, del_n_log(buf, 0, ENOENT);
        memcpy(buf+l0+1, "killer-file", 12); del_n_log(buf, 0, 0);
        buf[l0] = 0; del_n_log(buf, 1, 0);
}

static void close_all() {
	struct rlimit rlim; rlim.rlim_cur=9999; getrlimit(RLIMIT_NOFILE, &rlim);
	int i,n; for (i=0, n = (int)rlim.rlim_cur; i<n; i++) close(i); }

/*** output buf ***************/
static void con_end() { if (con_stat) kill(con_stat,9),con_stat=0; if (msg_fd>=0) close(msg_fd),msg_fd=-1; }

static void log_sn(const char *s, int n, int conpfx) {
	int r = 0, ix0 = log_ix, n0 = n;
	if (n<1) { switch(n) {
		case -2: if ((log_flg&1) && (r=write(msg_fd, logbuf+ix0, 16384-ix0))<0) goto err;
			 if ((r = write(msg_fd, logbuf, ix0))<0) goto err; else return;
		case -1: if ((n=ix0 - log_ox)>0 &&
			      (r = write(2, logbuf+log_ox, n), log_ox=ix0, r<1)) goto lerr;
		case 0: return;
	}}
	if ((log_ix+=n)>16383 && (log_flg|=1, log_ix&=16383)) memcpy(logbuf, s+(n=16384-ix0), log_ix);
	if (n) memcpy(logbuf + ix0, s, n);
	while ((log_ix^log_ox)&~4095) {
		int op2 = (log_ox&~4095) + 4096;
		    r   = write(2, logbuf + log_ox, op2 - log_ox);
		if (r<=0) { log_ox = log_ix; goto lerr; }
		log_ox = op2 & 16383;
	}
	if (msg_fd<0 || !con_stat || (n = n0 + conpfx)<1) return; else ix0 = log_ix - n;
	if (ix0<0 && (r=write(msg_fd, logbuf+(ix0&16383), -ix0), n+=ix0, ix0=0, r<1)) goto err;
	if ((r=write(msg_fd, logbuf+ix0, n))>0) return;
err:	return s = (r==-1) ? strerror(errno) : "zero ret", con_end(), LOG("[console write err: %s]", s);
lerr:   if (!con_stat) return;
	char errbuf[256]; n = snprintf(errbuf, 255, "log file write error: %d(%s)\n",
			r ? errno : 0, r ? strerror(errno) : "zero ret.");
	write(msg_fd, errbuf, n);
}

static void log_esc(int pfx, char *s) {
	if (pfx!='u') return LOG("invalid esc-cmd src 0x%x", pfx);
	if ((*s|1) != 49) return LOG("invalid gui esc-cmd 0x%x", *s);
	shutdn(8 ^- (gui_stat=*s&1));
}

static void log_l(char *s, int n, int pfx) {
	if (!memcmp(s, "&CMD", 4)) return log_esc(pfx, s+4);
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

static void fail(const char*s) { LOG("%s", s); log_sn(0,-1,0), exit(5); }

static void bye(int kf) {
	if (!kf) exit(1);
	con_end(); search_and_destroy(killer_fd); del_workdir();
	if (!rs_ts || time(NULL)>rs_ts+5) LOG("bye"), log_sn(0,-1,0), exit(0);
	const char *es, *fn = getenv("LF_BIN"); 
	if (!fn) LOG("restart failed (missing env)"), log_sn(0,-1,0), exit(1);
	LOG("restarting (%s)...", fn); log_sn(0, -1, 0); close_all();
	execl(fn, fn, (char*)0); es = strerror(errno);
	const char * log2 = getenv("LF_LOG"); if (!log2) exit(2);
	int fd2 = open(log2, O_WRONLY|O_APPEND); if (fd2<0) exit(3);
	write(fd2, es, strlen(es)); write(fd2, " (restart failed)\n", 18); close(fd2); exit(4);
}

/*** inotify ******************/
#define INEV_SIZ (sizeof(struct inotify_event))
static char * ino_nm(int t, int oi) { char*q = ino_inm+1+40*t; h5f(q+wrkdir_l+4,oi); return q; }
static int ino_find(int *p, int id) { BVFOR_JMC(ino_bv) if (p[j]==id) return (int)j;   return -1; }
static void ino_bye(int j) { unsigned int m = 1u<<j, k = !(m&ino_bv); 
			     if (k|(dflg&2)) LOG("BUG:ino_bye: %d"+4*!k, j);
			     ino_bv &= ~m; unlink(ino_nm(0,j)); unlink(ino_nm(1,j)); }

static void ino_ini() {
	int l = wrkdir_l;
	memcpy(ino_inm+ 1, wrkdir, l); memcpy(ino_inm+ 1+l, "/ds/00000.txt",14);
	memcpy(ino_inm+41, wrkdir, l); memcpy(ino_inm+41+l, "/d0/00000.txt",14);
	ino_inm[l+4] = ino_inm[l+44] = 0; 
	if (mkdir(ino_inm+1,0700)<0 || mkdir(ino_inm+41,0700)<0 || (ino_fd=inotify_init1(IN_CLOEXEC))<0 ||
	    inotify_add_watch(ino_fd, ino_inm+1, IN_CLOSE_WRITE)<0) 
		return ino_fd=-1, LOG("ino_ini: %s", strerror(errno));
	ino_inm[l+4] = ino_inm[l+44] = '/';
}

#define SERR(X) ((void) write(1, "_E" #X "\n", 5))
#define SERRC(X) (close(fd), (void) write(1, "_E" #X "\n", 5))
static void ino_add(int id, int edx, char * txt, int len) {
	if (dflg&2) LOG("ino_add: id=0x%x edx=%d, txt=\"%s\"", id, edx, txt);
	if (ino_fd<0) return SERR(15);
	if (ino_find(ino_oid, id)>=0) return SERR(16);
	BVFOR_JMC(~ino_bv) goto found;
	return LOG("dsc-edit table full (32)");
found:; char *nm0 = ino_nm(0, id), *nm1 = ino_nm(1, id);
	int i, fd = creat(nm1, 0600); if (fd<0) goto err2;
	for (i=0; i<len; i++) if (txt[i]==36) txt[i]=10; if (txt[len-1]!=10) txt[len++] = 10;
	txt-=7; len+=7; memcpy(txt, "#DSCR: ", 7);
	if (write(fd, txt, len)<0) goto err1; else close(fd);
	if (rename(nm1, nm0)<0 || creat(nm1, 0600)<0) goto err2;
	if (write(fd, txt, len)<0) goto err1; else close(fd);
	ino_bv |= (1u<<j); ino_oid[j] = id; 
	ino_epid[j] = edx ? (nm0[-1] = j+48, setenv("LF_ED_TMPF", nm0-1, 1),
			     launch(xapp_name[0],"!(x1","-e",lf_ed, NULL), -1)
			  : launch(xapp_name[2], "!(TT", nm0, NULL);
	return; // TODO: open ed
err1:   return close(fd), LOG("ino_add:w: %s", strerror(errno));
err2:	return LOG("ino_add:c/m: %s", strerror(errno));
}

static void ino_ev(struct inotify_event *ev) {
	const char * s = ev->name;
	if (dflg&2) LOG("ino_ev: name=\"%s\"", s);
	int i; for (i=0; i<5; i++) if ((unsigned int)(s[i]-48)>9u && (unsigned int)(s[i]-65)>5u) return;
	if (memcmp(s+5, ".txt", 5)) return;
	int id = atoi_h(s); if (dflg&2) LOG("ino_ev: close(w) for 0x%x", id);
	const char *nm0 = ino_nm(0, id), *nm1 = ino_nm(1, id);
	char buf0[2560], buf1[2560];
	int i0=8, r, r2, fd = open(nm0, O_RDONLY); if (fd<0) return SERR(19);
	if ((r =read(fd, buf0+i0, 2500))<0) return SERRC(19); else if (r==2500) return SERRC(18);
	close(fd); if ((fd=open(nm1, O_RDONLY))<0) return SERR(19);
	if ((r2=read(fd, buf1+i0, 2501))<0) return SERRC(19); else close(fd);
	if (r==r2 && !memcmp(buf0+i0, buf1+i0, r)) { if(dflg&2)LOG("ino_ev: file(0x%x)==backup",id); return;}
	fd=creat(nm1,0600); write(fd, buf0+i0, r); close(fd);
	if (buf0[i0+r-1]!=10) buf0[i0+r++] = 10;
	if (r>6 && !memcmp(buf0+i0, "#DSCR:", 6)) {
		for (i0+=6,r-=6; buf0[i0]!=10; i0++,r--);     }
	for (i=i0; i<i0+r-1; i++) if (buf0[i]==10) buf0[i]='$';
	buf0[i0-8]='H'; buf0[i0-7]='#'; h5f(buf0+i0-6, id); buf0[i0-1] = '$'; 
	write(1, buf0+i0-8, r+8);
	write(2, buf0+i0-8, r+8);
}

static void ino_read() {
        char buf[4096], *s = buf;   int r;
        if ((r = read(ino_fd, buf, 4096))<=0) { perror("read/ino"); return; }
        while(1) {
                struct inotify_event *ev = (struct inotify_event*)s;
                ino_ev(ev);
                int l = INEV_SIZ + ev->len; if ((r-=l)<=0) return; else s+=l;
        }}

/*** cmd/io_ibuf **************/
static int start_con() {
	static const char * sh = 0;
	return -((!sh && !(sh = getenv("LF_CON"))) || launch(xapp_name[0], "!(x1", "-e", sh, (char*)0)<0); }

static void con_started(char *s, int n) {
	if (memcmp(s, "_c/proc/", 8)) return LOG("_: invalid path\"%s\"", s);
	int pid = atoi(s+8); if (pid<1 || pid>65535) return LOG("/: invalid pid %d", pid);
	if (pid==con_stat) return LOG("/: already started, pid %d", pid); else con_end();
	if ((msg_fd = open(s+2, O_WRONLY|O_APPEND|O_NOCTTY))<0) return LOG("\"%s\": %s",s+2,strerror(errno));
	con_stat = pid; log_sn(0, -2, 0); --s[n-1]; s[n]=10; write(1, s, n+1);
}

static void set_xapp(int k, int l, const char *s) {
	if (dflg&4) LOG("set_xapp(%d, %d, \"%s\")", k, l, s);
	if ((unsigned int)k >= (unsigned int)N_XAPP) return LOG("set_xapp: k=%d", k);
	if ((unsigned int)l >= 255u)		     return LOG("set_xapp: l=%d", l);
	memcpy(xapp_nbuf+256*k, s, l); xapp_name[k] = xapp_nbuf+256*k;
	setenv(xapp_env[k], s, 1);
}

static void cmd_l(char *s, int n, int src) { s[n] = 0; switch(*s) {
	case 'c': if (start_con()<0) LOG("start_con failed"); return;
	case '_': return con_started(s, n);
	case '-': return ino_bye(s[1]-48);
	case 'Z': return con_end();
	case 'f': return log_sn(0, -1, 0);
	case 'r': return (void) (rs_ts = time(NULL));
	case 'x': case 'X': return set_xapp(s[1]-48, n-1, s+2);
	case 'h': return (n>7 && s[7]=='.') ? ino_add(atoi_h(s+2), s[1]&1,  s+8, n-8) : LOG("h: parse error");
	case 'd': return (void) (dflg = atoi_h(s+1));
	default: LOG("unknown cmd%c 0x%x (%s)", 48+src, *s, s); return;
}}

static void io_ib_read(ibuf *b) {
	int ec = ib_read(b); if (!ec) return;
	if (ec<-2) return LOG("ib_read: %s", "no file\0fragment discarded"+8*(ec&1));
	const char *es = (ec&1) ? strerror(errno) : "zero returned";
	int ag = b->arg, rof = (errtemp+=10)<99;
	close(b->fd); if (ag<65) return b->fd=-1, LOG("pipe(%c): %s -- bye", ag?ag:48, es), shutdn(2+2*!!ag);
	LOG("logpipe(%c): %s, eT=%d -- %s", ag, es, errtemp, rof?"reopening":"giving up"); if (!rof) return;
	b->fd = open(tpipe_name(ag), O_RDWR); LOG("logpipe(%c): reopen: %s", b->fd<0 ? strerror(errno):"OK");
}

static void io_ib_init(int mlog_fd) {
	if ((ib[1].fd = mlog_fd)&0xffff8000) LOG("mlog_fd = %d", mlog_fd); 
	int i; for (i=2; i<N_IBUF; i++) if ((ib[i].fd = open(tpipe_name(ib[i].arg), O_RDWR))<0)
		LOG("cannot open log pipe '%c': %s", ib[i].arg, strerror(errno));
}

/*** main *****/

static void i_chld(int _) { chld_flg = 1; signal(SIGCHLD, &i_chld); }

static void hello() {
	char buf[1024]; int k = (int)readlink("/proc/self/fd/2", buf, 1023);
	if (k<0) LOG("ioprc started, /proc/self/fd/2: %s", strerror(errno));
	else buf[k] = 0, LOG("ioprc started, logfile is \"%s\"", buf);
}

static void ch_bye() { 
	int j, x, pid = waitpid(-1, &x, WNOHANG);
	if (pid<=0) return LOG("wait: %s", pid<0 ? strerror(errno) : "WTF???");
 	if (dflg&1) LOG("wait: pid=%d", pid);
	if ((j=ino_find(ino_epid, pid))>=0) ino_bye(j);
}

#define FOR_IB for(i=0; i<N_IBUF; i++) if ((k=ib[i].fd)>=0)
int i_main(int ac, char ** av) {
	signal(SIGPIPE, SIG_IGN); signal(SIGHUP, SIG_IGN); signal(SIGINT, SIG_IGN); signal(SIGCHLD, i_chld);
	hello();   int i,k;
	if (!(wrkdir=getenv("LF_TMPDIR")) || (wrkdir_l=strlen(wrkdir))>20) fail("workdir");
	if (!(lf_ed=getenv("LF_ED"))) fail("ed-wrap");
	ino_ini();
	for (i=0; i<N_XAPP; i++) if (!(xapp_name[i] = getenv(xapp_env[i]))) xapp_name[i] = xapp_dflt[i];
	io_ib_init(ac<2 ? -1 : qh4r(*(int*)av[1]));
	killer_fd = (ac<3||*av[2]<48) ? 0 : qh4r(*(int*)av[2]);
	der_alte = ac<4 ? getppid() : qh4r(*(int*)av[3]);
	fd_set rset; struct timeval tv;
	while(1) {
		FD_ZERO(&rset); int maxfd = 0;
		FOR_IB 		   { FD_SET(k, &rset); if (k>maxfd) maxfd = k; }
		if ((k=ino_fd)>=0) { FD_SET(k, &rset); if (k>maxfd) maxfd = k; }
		tv.tv_sec=0; tv.tv_usec = selwt;
		int r = select(maxfd+1, &rset, 0, 0, &tv);
		if (r<0) { LOG("select(): %s", strerror(errno)); r = 0; }
		if (chld_flg) chld_flg = 0, ch_bye();
		shutdn(0);
		if(!r){ if ((errtemp -= 20)<0) errtemp = 0;  continue; }
		if ((errtemp -= r)<0) errtemp = 0;
		FOR_IB if (FD_ISSET(k, &rset)) io_ib_read(ib+i);
		if (FD_ISSET(ino_fd, &rset)) ino_read();
	}}

/**************** cmd playback (quicktest) *********************************/

static int lim_arr[96];
#define LIM(X) lim_arr[(int)(X)-32]

static int m_sleep(int x) { if (!x) return 0;
        struct timespec ts; ts.tv_sec = x/1000; ts.tv_nsec = 1000000*(x%1000); int ec;
        do ec = nanosleep(&ts, &ts); while (ec==EINTR); return ec; }

#define qh4rs(P) qh4r(*(const int*)(P))
#define IFHX(X,T) ( (unsigned int)((X)-48)<10u || (unsigned int)((X)-(T))<6u )

static void he(const char *s1, const char *s2) { fprintf(stderr, "%s: \"%s\"\n", s1, s2); }

static void qpcmd(const char *s) { switch(*s) {
	case 'L': for (++s; *s; ) {
			  while (*s==32) ++s;  if (*s==10) return;
			  int v=0, d, k = *s-32; if ((unsigned int)k>95u) return he("invalid ch",s);
			  if (*s=='*') { ++s; v = INT_MAX; }
			  else         { for (++s; (unsigned int)(d=*s-48)<10u; s++) v = 10*v + d; }
			  lim_arr[k] = v; if (!lim_arr[31]) fprintf(stderr, "lim '%c' set to %d\n", k+32, v);
		  } return;
	default: return he("unknown cmd", s);
}}

void qp1(FILE *f) {
	char buf[1024], *s, *q; buf[1023]=0;
	while (fgets(s=buf+2, 999, f)) {
		int l = 0, k = IFHX(*s, 97) ? qh4rs((s+=5)-5) : 0;
		switch(*s) {
			case '#': m_sleep(k); if (s[1]&&s[1]!=10) fprintf(stderr,"%s",s); continue;
			case '?': qpcmd(s+1); continue;
			case '^': goto sl;
			case 'm': case 'N': case 'Z': l = *s; goto lim;
			case 'Y': l = 'Z'; goto lim;
			case 'V': for (q=s+1; *q && *q!=36; q++); if (!*q) goto wtf; else q+=3;
				  switch(*(q+=2*(IFHX(q[0], 65) && IFHX(q[1], 65)))) {
					  case'c':case'e':case'k':case'g': l=*q;goto lim; default: goto wtf; }
			case 10:  continue;
			default:  goto wtf;
		}
lim:		l = LIM(l); if (k>l) k = l; goto sl;
wtf: 		fprintf(stderr,"warning: \"%s\": unknown cmd\n", s);
sl:		m_sleep(k);
		l = strlen(s); s[-2] = 'Q'; s[-1] = 'P'; if (s[l-1]!=10) s[l++] = 10;
		write(1, s-2, l+2);
	}}

static char *at_fnm, *at_fxy;
static int at_n, at_k[26], at_zf = 0;

static FILE * at_f_ij(int i, int j) {
	return at_fxy[0] = 65 + i, at_fxy[1] = j<0 ? 'z' : 48+j, fopen(at_fnm, "r"); }

static void qp1_tij(int t, int i, int j) {
	FILE * f = at_f_ij(i,j); 
	printf("_>========== %s - %s =========\n", at_fnm, f ? "ok" : strerror(errno)); fflush(stdout);
	f ? (m_sleep(t), qp1(f), fclose(f)) : perror(at_fnm);
}

static void at_dir(const char * dir) {
	int i, j, l = strlen(dir);   FILE *f;
	memcpy(at_fnm = malloc(l+6), dir, l); memcpy(at_fnm+l, "/t", 2);
	memcpy(at_fxy = at_fnm+l+2, "xy\0", 4);
	for (i=0;;i++) {
		for (j=0; (f=at_f_ij(i,j)) ;j++) fclose(f);
		if (!j) { at_n = i; break; }
		at_k[i] = j; if ((f=at_f_ij(i,-1))) fclose(f), at_zf |= 1<<i;
		fprintf(stderr," [%d: 0..%d]", i, j-1);
	}
	fprintf(stderr, " n=%d, zf=0x%x\n", at_n, at_zf);
}

static void at_random() {
	srandom(time(NULL));
	int i,k,st[32]; for (i=0; i<at_n; i++) st[i] = 1;
	while(1) {
		i = random() % at_n;
		if (!st[i]) qp1_tij(200, i, 0), st[i] = 1;
		else if (st[i]<at_k[i]) qp1_tij(200, i, st[i]), st[i]++;
		else if ((k=random()%st[i])<at_k[i]-1) qp1_tij(200, i, k+1), st[i]++;
		else qp1_tij(200, i, -1), st[i] = 0;
	}}

static void at_main(const char * opt) {
	int i,j,k, flg = (*opt=='-') ? (opt+=2, hxd2i(opt[-1])) : 0;
	at_dir(opt), write(1, "_m\n", 3);
	if (flg&1) for (i=0; i<at_n; i++) qp1_tij(200, i, 0), write(1, "_m\n", 3);
	if (flg&2) for (i=0; i<at_n; i++) for (j=1,k=at_k[i]; j<k; j++) qp1_tij(200, i, j), write(1, "_m\n", 3);
	if (flg&4) for (i=0; i<at_n; i++) at_random(), write(1, "_m\n", 3);
	if (flg&8) for (i=0; i<at_n; i++) qp1_tij(200, i, -1), write(1, "_m\n", 3);
}

static int q_main(int ac, char** av) {
	int i; for (i=0; i<96; i++) lim_arr[i] = INT_MAX;
	if (ac<2) qp1(stdin);
	FILE * f;
	for (i=1; i<ac; i++) {
		if (*av[i] == '-') { switch(av[i][1]) {
			case 'a': at_main(av[i]+2); continue;
			default: fprintf(stderr, "unknown option \"%s\"\n", av[i]); 
		}}
		if ((f=fopen(av[i],"r"))) qp1(f), fclose(f); else return perror(av[i]), 1;
	}
	return 0;
}

/* lf.tlogdump (debug tool for lflab)
   dumps tlog in a more-or-less human-readable format (time:event pairs)
   time: duration (time until next event) in microseconds unless "m" says millisec.
   event: (s) - select (with wait)         +-  - clock adjust
   	  (S) - select (without wait)      j88 - job slice (j1: longest)
   	  888 - calculate main mixer       (J) - no time for job exec    
          (p) - play (send to audio dev)
   every 4th cycle the data is send to GUI process (for CPU meter display)
   (this has no separate event, you can see that every 4th "clock adjust" takes a little longer)
*/

static void dur4(char * to, int t) {
	if (t<=9999) return (void)sprintf(to, "%4d", t);
	if (t>262141) return (void)memcpy(to, "FRVR-WTF"+4*(t&1), 4);
	if (t>99999) return (void)sprintf(to, "%dm", t/1000);
	t/=100; sprintf(to, "%dm%c", t/10, 48+(t%10));
}

static int ev4(char *to, int k) {
	if (k<128) to[0]=40, to[2]=41, to[3]=32, to[1] = k?k:63;
	else if (k<2048) return sprintf(to, "%+d", k-1090), k-1090;
	else if (k<4096) (k&=2047)<1000 ? sprintf(to, "j%d", k) 
				        : (k/=100, sprintf(to, "j%ck%c",48+k/10,48+k%10));
	else k&=4095, sprintf(to, " %d"+(k>999), k);
	return 0;
}

#define MEGA 1000000

static int  ent10(char *q, unsigned int x) { return dur4(q,x&262143), q[4]=q[9]=58, ev4(q+5, x>>18); }
static void top10(char *q, long long t) { sprintf(q, "%02d.%06d", (int)((t/MEGA)%60), (int)(t%MEGA)); }

static int t_main(int ac, char ** av) {
	long long t0 = 0LL, t1 = 0LL;
	int nrow = 30, ncol = 10, nblk = 0, atot2 = 0;
	if (ac==4) nrow=atoi(av[1]), ncol = atoi(av[2]), av += 2;
	else if (ac!=2) return fprintf(stderr,"usage: %s [nrow ncol] <binfile>\n", *av), 1;
	int fd = open(av[1], O_RDONLY); if (fd<0) return  perror(av[1]), 1;
	int blks = nrow * ncol, ibuf[blks], wid = 10*ncol, osiz0 = (nrow+1)*wid;
	char obuf[osiz0 + 80]; memset(obuf, 32, osiz0);
	while(1) {
		int i,j,k,nr = read(fd, ibuf, 4*blks), eof = 4*blks - nr;
		t1 = t0;
		if (nr<=0) return nr ? (perror("read"),1) : 0;
		if (nr&3) fprintf(stderr, "ignoring %c extra bytes\n", nr&3);
		int nrow2 = ((nr>>=2) + ncol - 1) / ncol, osiz2 = (nrow2+1)*wid;
		int adjm = 0, adjM = 0, adjtot = 0;
		for (i=0,j=0,k=0; k<nr; k++, i += (++j>=nrow2 && !(j=0))) {
			if (!j) { top10(obuf+10*i, t0); obuf[10*i+9] = (i==ncol-1)?10:'|'; }
			int a = ent10(obuf + 10*i + wid*(j+1), ibuf[k]); t0 += ibuf[k]&262143;
			if (a) { adjtot+=a; if (a<adjm) adjm=a; else if (a>adjM) adjM=a; }
		}
		for (i=0; i<osiz2; i++) if (!obuf[i]) obuf[i] = 32;
		for (i=9; i<osiz2; i+=10) obuf[i] = 32; 
		for (i=wid-1; i<osiz2; i+=wid) obuf[i] = 10;
		int sec = t0/MEGA, nx;
		nx=sprintf(obuf+osiz2, "blk%04d t%02d:%02d.%06d len:%d.%06d +-min,max,tot,atot: %d,%d,%d,%d\n",
		    	nblk++, sec/60, sec%60, (int)(t0%MEGA), (int)((t0-t1)/MEGA), (int)((t0-t1)%MEGA),
			adjm, adjM, adjtot, atot2+=adjtot);
		write(1, obuf, osiz2+nx);
		if (eof) return 0;
	}}

/**************** worker proc. (plot+tbd) **********************************/

#undef LOG
#define LOG LOG_wrk
#define LOG_E(S,X) LOG_wrk(NULL, (X), errno, (S))
#define GP_SAMPSIZ 2048
#define GP_RGB(J) ( gp_rgb + 8 * ((gp_u_perm7>>(4*(J)-4)) & 7) )
#define GP_N_CMDK (sizeof(gp_cmd_key)/sizeof(char*))

typedef int (*gp_statfun)(int, int);
typedef int (*gp_ls_fun)(char *to, const double *x, const double *y, int n);

static void LOG_wrk(const char * fmt, ...), wrk_cmd_l(char*,int,int);
static int statf_null() { return LOG("statf_null called!"), 0; }
static gp_statfun gp_cur_statfun = &statf_null;
static double * samp2gplot = NULL;
static int gp_inpipe = -1, gp_outpipe = -1, gp_pid = 0, gp_w_id = 0, gp_len = 0,
	   gp_res22 = 4, gp_u_msk = 0, gp_u_lr = 0, gp_u_perm7 = 0x6543210;
static double gp_ctr = 0.5, gp_rad = 0.5;
static char gp_w_title[64], gp_f_title[64];
static const char gp_rgb[] = "#ff0000\0#00a000\0#0000ff\0#a08000\0#009090\0#e000e0\0#000000";
static int gp_tlog_n = 0, gp_tlog_siz = 0, *gp_tlog = NULL;
static double * gp_ptf_dat = 0;
static int gp_ptf_siz, gp_ptf_bits, gp_ptf_flg;

static void wbye(int j) { LOG("bye %d", gp_pid); if (gp_pid>0) kill(gp_pid,9); exit(j); }
static inline int gp_res() { return (128+(53&-(gp_res22&1))) << (gp_res22>>1); }

static char wrk_cb_s[1024];
static ibuf wrk_cb = { 0, 1024, 0, 0, &wrk_cmd_l, wrk_cb_s };
#define GP_CKY17(J) "gbf" #J "\0" #J , "gbs" #J "\0ctrl-" #J
static const char * gp_cmd_key[] = {"gmX4\0#Left", "gmx4\0#Right", "gmc4\0#Down", "gmC4\0#Up", "gc\0c",
	"gmX@\0&#Left", "gmx@\0&#Right", "gm0\0^#Left", "gm1\0^#Right", "gmX1\0#End", "gmx1\0#PageDown",
	"gm*\0^#Down", "gm+\0^#Up", "gm<4\0#Home", "gm>4\0#PageUp", "gr-\08", "gr+\09",
	GP_CKY17(1), GP_CKY17(2), GP_CKY17(3), GP_CKY17(4), GP_CKY17(5), GP_CKY17(6), GP_CKY17(7) };

static void LOG_wrk(const char * fmt, ...) {
        char buf[1024]; int k,r; va_list ap; va_start(ap, fmt);
	if (fmt) r = vsnprintf(buf, 1023, fmt, ap); 
	else if ((r = va_arg(ap,int))>=0) return va_end(ap);
	else k = va_arg(ap,int), fmt = va_arg(ap,const char*),
	     r = snprintf(buf, 1023, "%s: ec=%d ue: %s", fmt, r, strerror(k));
        va_end(ap); if (r>=0) buf[r] = 10, write(2, buf, r+1);
}
//// w/gnuplot /////
static int gp_start() {
	static const char ini[] = "set style data line\nset autoscale fix\n"
		"set y2tics nomirror format \"%g\"\n"; // "bind c 'print \"%8K^q!\"'\n";
        if (gp_pid>0) return -1;
	if ((gp_pid = launch("gnuplot", "!><1", &gp_inpipe, &gp_outpipe, (char*)0)) < 0) return -2;
	if ((gp_inpipe|gp_outpipe)<0) return -3;
	int i; char buf[4096], *q = buf + sizeof(ini)-1; memcpy(buf, ini, sizeof(ini)-1);
	for (i=0;i<GP_N_CMDK; i++) { 
		const char *s = gp_cmd_key[i];  while (*++s); ++s;
		int j, kf = 0; for (j=0;j<3;j++) if (*s=="^&#"[j]) ++s, kf |= (1<<j);
		const char *mod = "\0(==#==)ctrl-\0__alt-\0___BUG" + 8*(kf&3);
		q +=		sprintf(q,"bind '%s%s' 'print \"%%8K^q%c\"'\n", mod, s, i+40);
		if(kf&4) q+= sprintf(q,"bind '%sKP_%s' 'print \"%%8K^q%c\"'\n", mod, s, i+40);
	}
	if (write(gp_inpipe, buf, q-buf) != q-buf) return -4;
        return 0;
}

static int gp_binplot(int n, int flg) {
	if (n<0) n = flg>>16;
	char buf[4096], *q=buf;
	int ch = ' ', yf = (flg&254) & ~gp_u_msk, y2f = yf & ((flg>>8) ^ gp_u_lr);
	if (!yf) return LOG("gp_binplot: nothing to do!"), 0;
	if (!(flg&1) && n>2) {  double *px=samp2gplot, x=px[0], stp = (px[8*n-8]-x)/(double)(n-1);
				int i; for (i=1; i<n-1; i++) *(px+=8) = (x+=stp); }
	q += sprintf(q, "set term x11 %d title '%s'\nset ytics %s%s%s", gp_w_id, gp_w_title,
			"no"+2*!y2f, "mirror\nset y2tics nomirror format '", "%g'\nplot" + 2*!y2f);
	BVFOR_JMC(yf) { q += sprintf(q, "%c'%s%s%d%s",  ch, tpipe_name(','), "' binary filetype=raw record=",
							n, " format='%double");   ch = ',';
			int i,k,rgt=!!(y2f&(1<<j)); char b4[8];
			for (i=1;i<8;i++) k=(i==(int)j), memcpy(q, "%*double%double"+8*k, 8), q+=8-k;
			b4[2-2*rgt]=0; b4[1+rgt]='|'; b4[3*rgt]=48+j; b4[4]=0;
			q += sprintf(q,"' title '%s%s%s' lc rgbcolor '%s'", b4,gp_f_title+8*j,b4+2, GP_RGB(j));
			if (rgt) memcpy(q," axes x1y2",10), q+=10;  }
	*(q++)='\n'; return write(gp_inpipe, buf, q-buf) == (q-buf) ? 0 : -1;
}

static int gp_plot(int d) { 
	if (d) gp_w_id+=d, gp_u_msk=gp_u_lr=0, gp_ctr=gp_rad=0.5;
	double dl = (double)gp_len, dc = dl*gp_ctr, dr = dl*gp_rad;
	int j0 = ivlim((int)lround(dc-dr), 0,  gp_len-1),
	    jz = ivlim((int)lround(dc+dr), j0+1, gp_len);
	return gp_binplot(-1,(*gp_cur_statfun)(j0, jz)); }

static void gp_rgb_shuffle() {
	int i, j, k, i4, j4, x, r = random(), v = gp_u_perm7;
	for (i=0; i<6; i++) k = 7-i, j=r%k, r/=k, i4=4*i, j4=4*j, x=((v>>i4)^(v>>j4))&7, v^=((x<<i4)|(x<<j4));
	gp_u_perm7 = v; }

///// w/tlog ////
static int gp_tlog_t(int ix) {  int i, i6 = ix>>6, t = gp_tlog[gp_tlog_n + i6];
				for (i=64*i6; i<ix; i++) t += gp_tlog[i]&262143;   return t; }

static int gp_statf_tlog(int j0, int jz) {
	int x, x0, k0, j, samp_cnt = 0, samp_calc = 0;
	double cpu0 = 0.0, *q = samp2gplot;
	int ttot = 0; for (j=j0; j<jz; j++) ttot += gp_tlog[j]&262143;
	int td1, t = gp_tlog_t(j0), td = 0, tdloc = 0, tdmin = ttot/(gp_res()-1), t1 = tdmin;
	for (j=j0,x0=0; j<jz; j++,x0=x) {
		int k = (x=gp_tlog[j]) >> 18, xt = x&262143;
		if (k>127) { if (k<4096) td1 = (k-1090), td += td1, tdloc += td1;  }
		else { if ((k|1)=='q' && (k0=(x0>>18)-4096)>=0) samp_cnt += k0, samp_calc += xt+(x0&262143); }
		if ((t+=xt)>t1) q[3] = !samp_cnt ? cpu0 : (cpu0=4.41*(double)samp_calc/(double)samp_cnt,
							   samp_calc=samp_cnt=0, cpu0),
				q[2] = 1e-3*(double)tdloc, tdloc = 0,
				q[0] = 1e-6*(double)t, q[1] = 1e-3*(double)td, t1 = t+tdmin, q += 8; }
	return 0x80f + ((q-samp2gplot)<<13);
}

static int gp_tlog_read(const char *fname) {
	memcpy(gp_w_title, "tlog", 5); memcpy(gp_f_title+8, "td.A\0___td.+-\0__cpu%\0__", 24);
	int fd = open(fname, O_RDONLY); if (fd<0) return -1;
	int n = lseek(fd, 0, SEEK_END); if (n<4) return -2; else gp_len=gp_tlog_n=(n>>=2);
	int r = lseek(fd, 0, SEEK_SET); if (r) return -3;
	int n6 = (n+63)>>6, siz = n + n6;
	if (siz>gp_tlog_siz) {
		if (siz>(1<<28)) return -4;
		int sz = 2*gp_tlog_siz; if (sz) free(gp_tlog); else sz = 4096;
		while (sz<n) sz+=sz;
		gp_tlog = malloc(4*(gp_tlog_siz=sz)); }
	if ((r = read(fd, gp_tlog, 4*n)) != 4*n) return close(fd), -6; else close(fd);
	int i, j, k, t = 0, *q = gp_tlog+n; 
	for (*(q++)=0, i=0; (k=i+64)<n; i=k,*(q++)=t) for (j=i; j<k; j++) t += gp_tlog[j] & 262143;
	gp_cur_statfun = &gp_statf_tlog; return 0;
}

///// w/fft ////////
inline unsigned int bit_rev(unsigned int n, int bits) {
        return ( (unsigned int)bitrev8[n>>24] + ((unsigned int)bitrev8[(n>>16)&255]<<8) +
                ((unsigned int)bitrev8[(n>>8)&255]<<16) + ((unsigned int)bitrev8[n&255]<<24) ) >> (32-bits); }

#define REV_S_BITS 5
#define REV_S_BLK (1<<REV_S_BITS)
#define REV_XC1 (t=re[y],re[y]=re[x],re[x]=t)

static void rev_small1(double *re, int bits) {
        int x, y, n = 1<<bits;  double t; for (x=0; x<n; x++) if ((y=bit_rev(x, bits))>x) REV_XC1; }

static void rev_big1(double *re, int bits) {
        int lo, hi, mi, n = 1<<bits, lhbits = (bits+1-REV_S_BITS)>>1, mbits = bits - 2*lhbits,
            x,y,nlo = 1<<lhbits, hstep = 1<<(lhbits+mbits);
        for (hi=0; hi<n; hi+=hstep) { for (lo=0; lo<nlo; lo++) {
                double t; int mi2, bs = hi+lo, rbs = bit_rev(bs, bits);
                if (bs<rbs) continue;
                if (bs==rbs) { for (mi=0; mi<hstep; mi+=nlo) if ((mi2=bit_rev(mi, bits))>mi)
                                                                  x=bs+mi, y= bs+mi2, (REV_XC1); }
                else { for (mi=0; mi<hstep; mi+=nlo) x=bs+mi, y=rbs+bit_rev(mi,bits), (REV_XC1); }}}}

static void bit_rev_blk(double *re, int bits) { (bits>REV_S_BITS?rev_big1:rev_small1)(re, bits); }

#define RE(J) ri[2*(J)]
#define IM(J) ri[2*(J)+1]

#define FFT_PASS12(P) double x00=P[j], x01=P[j+1], x02=P[j+2], x03=P[j+3], \
        x0=x00+x01, x1=x00-x01, x2=x02+x03, x3=x02-x03; \
        RE(j)   = x0+x2; IM(j)=0.0; RE(j+2) = x0-x2; IM(j+2) = 0.0;\
        RE(j+1) = RE(j+3) = x1; IM(j+1) = x3; IM(j+3) = -x3

static void fft_pass12_1(double *ri, int n) {
        int j; double *q = ri + n; for (j=0; j<n; j+=4) { FFT_PASS12(q); }}

static void fft_pass12_0(double *ri, int n) {
        int i, j, k;
        for (i=(n-1)&~4095; i>0; i-=4096) for (j=i,k=i+4096; j<k; j+=4) { FFT_PASS12(ri); }
        for (j=(n-4)&4092; j>=0; j-=4) { FFT_PASS12(ri); }}

#define ANG1_BITS 9
#define ANG1_SIZ (1<<ANG1_BITS)
static void fft2(double * ri, int n) {
        int m, m2, df = (n>(1<<19));
        double oAr[2*ANG1_SIZ], *oAi=oAr+ANG1_SIZ, circ = /*reverse ? (-2.0*M_PI) :*/ (2.0*M_PI);
        for (m=8,m2=4; m<=n; (void)(df&&(LOG("fft2: %d/%d", m,n),1)),m2=m,m*=2) {
                int mm = m2-1;
                double fi2, fi = circ / (double)m, o1i = sin(fi), o1r = cos(fi);
                if (m2<=ANG1_SIZ) {
                        int j, k, k0, l;
                        double oR = 1.0, oi = 0.0;
                        for (j=0; j<m2; j++) { if (!(j&63)) fi2 = fi*(double)(j), oR = cos(fi2), oi = sin(fi2);
                                               oAr[j] = oR; oAi[j] = oi;
                                               double tR = oR*o1r - oi*o1i; oi = oi*o1r + oR*o1i; oR = tR; }
                        for (j=0; j<n; j+=m) { for (k0=0; k0<m2; k0++) { k = k0 + j;
                                double rel = RE(l=k+m2), iml = IM(l),
                                       tr = rel*oAr[k0] - iml*oAi[k0], ti = iml*oAr[k0] + rel*oAi[k0];
                                RE(l) = RE(k)-tr; RE(k)+=tr; IM(l)=IM(k)-ti; IM(k)+=ti; }}}

                else {
                        int j, k, k0, l; double oR = 1.0, oi = 0.0;
                        for (j=0; j<n; j+=m) { for (k0=0; k0<m2; k0+=ANG1_SIZ) {
                                int k1, jk0 = j+k0;
                                for (k1 = 0; k1<ANG1_SIZ; k1++) {
                                        k = jk0 + k1; if (!(k&63)) fi2 = fi*(double)(k&mm), oR = cos(fi2), oi = sin(fi2);
                                        double rel = RE(l=k+m2), iml = IM(l),
                                               tr = rel*oR - iml*oi, ti = iml*oR + rel*oi;
                                        RE(l) = RE(k)-tr; RE(k)+=tr; IM(l)=IM(k)-ti; IM(k)+=ti;
                                        double tR = oR*o1r - oi*o1i; oi = oi*o1r + oR*o1i; oR = tR; }
                        }}}}}

static void fft(double * re_im, int bits, int flg) { // 1: high
        int n = 1<<bits, df = (bits>19), ofs = n &- (flg&1); if (df) LOG("fft: start");
        bit_rev_blk(re_im+ofs, bits); 			     if (df) LOG("fft/bitrev done");
        (ofs?fft_pass12_1:fft_pass12_0)(re_im, n);	     if (df) LOG("fft/pass12 done");
        fft2(re_im,n); }

///// w/plot_tF ////
#define GP_PT_L1(S,M) for (j=1,x=y=acc=p[0]; j<(S); j++) { acc+=(z=p[j]); if (z<x) x=z; else if (z>y) y=z; } \
		 to[0] = (M)*acc; to[1] = x; to[2] = y
static void gp_pt_mma(double *to, const double *p, int cnt, int siz, int rem) {
	double x,y,z,acc,sz_1; int i,j;
	switch(siz) { case 1: for (i=0; i<cnt; i++, to+=8) *to=p[i];   return;
		      case 2: for (i=0; i<cnt; i++, to+=8,p+=2) to[0] = .5*((x=p[0])+(y=p[1])),
				      (x<y) ? (to[1]=x,to[2]=y) : (to[1]=y,to[2]=x);
			      if (rem) to[0]=to[1]=to[2] = p[0];       return;
		      default: break; }
	for (i=0,sz_1=1.0/(double)siz; i<cnt; i++, to+=8, p+=siz) { GP_PT_L1(siz, sz_1); }
	if (rem) { GP_PT_L1(rem, 1.0/(double)rem); }}

static void gp_pt_alp(double *to, const double *p, int cnt, int siz, int rem) {
	double x, y; int i;
	if (siz==1) { for (i=0; i<cnt; i++, to+=8,p+=2) x=*p,y=p[1], to[3]=log(to[0]=sqrt(x*x+y*y)),
							to[6] = atan2(y,x);  return; }
	double z, xa, ya, za, zm, zM, sz_1 = 1.0/(double)siz;
	int j, i1 = cnt + !!rem;
	for (i=0; i<i1; i++, to+=8,p++) { 
		if (i==cnt) siz = rem, sz_1 = 1.0/(double)(rem);
		for (j=1, xa=*p,ya=*++p,zm=zM=za=xa*xa+ya*ya ; j<siz; j++) {
			x=*++p, y=*++p, xa+=x,ya+=y; za += (z = x*x+y*y); if(z>zM) zM=z; else if(z<zm) zm=z; }
		to[3] = log(to[0] = sqrt(za*sz_1)); to[4] = log(to[1] = sqrt(zm));
		to[6] = atan2(ya, xa); 		    to[5] = log(to[2] = sqrt(zM)); }}

static int gp_statf_tf(int j0, int jz) {
	if (!gp_ptf_dat) return LOG("gp_statf_tf: no data no draw"), -2;
	int nsamp = jz-j0, res = gp_res()-1, siz = 1+nsamp/res, cnt = nsamp/siz, rem = nsamp%siz,
	    gf = gp_ptf_flg, rcnt = cnt+!!rem, rv = rcnt<<16; //| (gf&8*((siz>1 ? 0x7ef0fe : 0x1280);
	LOG("gp_statf_tf: j0=%d, jz=%d, len=%d ==> rcnt=%d cnt=%d siz=%d rem=%d gsiz=%d",
			  j0,    jz,    gp_len,    rcnt,   cnt,   siz,   rem, gp_ptf_siz);
	double x, jj0 = (double)j0, jjz = (double)jz, *to = samp2gplot, *p = gp_ptf_dat+j0;
	if(gf&8) return x = 44100.0 / ((double)(1<<gp_ptf_bits)), to[0] = jj0*x, to[8*rcnt-8] = jjz*x,
			gp_pt_alp(to+1, p+j0, cnt,siz,rem), rv | ((siz>1?0xf0fe:0x9092));
	gp_pt_mma(to+1, p, cnt, siz, rem); to[0] = jj0/44100.0; to[8*rcnt-8] = jjz/44100.0;
	return (gf&1) ? (gp_pt_mma(to+4, p+gp_ptf_siz, cnt, siz, rem), rv | (siz>1?0x7e:0x12))
		      : rv | (siz>1?0xe:0x2);
}

static void gp_ptf_drop(){ if(gp_ptf_dat) munmap(gp_ptf_dat,16<<gp_ptf_bits), gp_ptf_dat=NULL, gp_ptf_bits=0;
			   if(gp_cur_statfun==&gp_statf_tf) gp_cur_statfun = &statf_null; }

static int gp_ptf_ini(const char *s) {
	int flg = hxd2i(*s), bits = s[1]-48; if (s[2]!=',') return LOG("gs: parse error"), -1;
	int len = 0; for (s+=3; *s>47; s++) len = 16*len + hxd2i(*s); 
	if (*s==',' && *++s) strncpy(gp_w_title, s, 63); else memcpy(gp_w_title, "unnamed", 8);
	LOG("gp_ptf_ini: flg=%d bits=%d len=%d", flg, bits, len);
	if (bits != gp_ptf_bits) { gp_ptf_drop(); if (!bits) return 0;
				   gp_ptf_dat = map_wdir_shm('+', 16<<(gp_ptf_bits=bits), 1); 
				   if (gp_ptf_dat==MAP_FAILED) return gp_ptf_dat=NULL, gp_ptf_bits=0, -2;  }
	int siz = 1<<bits, k = siz - len;
	if (k) { memset(gp_ptf_dat+len, 0, 8*k); if (flg&1)  memset(gp_ptf_dat+siz+len, 0, 8*k); }
	gp_ptf_siz = siz; gp_cur_statfun = &gp_statf_tf; gp_ptf_flg = flg; 
	if (flg&8) fft(gp_ptf_dat, bits, (flg>>1)&1), gp_len=siz/2,
		   memcpy(gp_f_title+8, "avg\0____min\0____max\0____lg:avg\0_lg:min\0_lg:max\0_phase\0_",56);
	else gp_len=len, memcpy(gp_f_title+8, "avg\0____min\0____max\0____avgR\0____minR\0____maxR\0___", 48); 
	return 0;
}

///// w/cmd  ///////
static int gp_move(int c, int v) {
	double g_c = gp_ctr, g_r = gp_rad, sg = (c&32) ? 1.0 : -1.0, vv = .125*(double)v,
	       x, r_min = (gp_len<127) ? 0.5 : 63.0/(double)gp_len;
	switch(c|32) {
		case ' ': goto rchk;
		case 'c': g_r *= exp(sg*M_LN2*vv); goto rchk;
		case '+': g_r = r_min;     goto cchk;
		case '*': g_c = g_r = 0.5; goto echk;
		case 'x': g_c += sg*g_r*vv; goto cchk;
		case '0': g_c = g_r; goto echk;
		case '1': g_c = 1.0-g_r; goto echk;
		case '<': x = g_c-g_r; if ((g_r*=exp(-M_LN2*vv)) < r_min) g_r = r_min; g_c = x+g_r; goto cchk;
		case '>': x = g_c+g_r; if ((g_r*=exp(-M_LN2*vv)) < r_min) g_r = r_min; g_c = x-g_r; goto cchk;
		default:  return LOG("gp_move: invalid cmd 0x%x", c), 0;
	}
rchk:	if (g_r<r_min) g_r=r_min; else if (g_r>.5) g_r = .5;
cchk:   if (g_c<g_r) g_c = g_r; else if (g_c>1.0-g_r) g_c = 1.0-g_r;
echk:	return (fabs(g_c-gp_ctr)>1e-11 || fabs(g_r-gp_rad)>1e-11) ? (gp_ctr=g_c, gp_rad=g_r, 1) : 0;
}

static int gp_cmd(const char *s) {
	LOG("gp_cmd: \"%s\"", s);
	int r; switch(*s) { 
		case 't': return ((r=gp_tlog_read(s+1))<0) ? r : gp_plot(1);
		case 's': return ((r=gp_ptf_ini  (s+1))<0) ? r : gp_plot(1);
		case 'S': return gp_ptf_drop(), 0;
		case 'c': return gp_rgb_shuffle(), gp_plot(0);
		case 'm': return gp_move(s[1], s[2]-48) ? gp_plot(0) : (LOG("gp/move: nothing happens"), 0);
		case 'r': return ((s[1]&2) ? (gp_res22<8&&++gp_res22) : (gp_res22&&gp_res22--)) ? gp_plot(0)
					   : (LOG("gp_resol: nothing happens"), 0);
		case 'b': return *((s[1]&1) ? &gp_u_lr : &gp_u_msk) ^= ((1<<(s[2]-48))&0xfe), gp_plot(0);
		default: return -9;
	}}

static void wrk_cmd_l(char *s, int n, int src) { if (n>=0) /*!!*/ s[n] = 0; switch(*s) {
	case 'g': LOG_E("gp_cmd", gp_cmd(s+1)); return;
	case 'q': bye(0);
	default: LOG("unknown command \"%s\"", s); break;
}}

static void gp_out() {
	static int state = 37;
	signed char buf[4096];
	int i, c, r = read(gp_outpipe, buf, 4096);
	if (r<=0) perror("gp_pipe/read"), exit(1);
  	LOG("gp_read: %d \"%s\"", r, r<999 ? (const char*)(buf[r-(buf[r-1]==10)]=0,buf) : "...");
	for (i=0;i<r;i++) { if	    ((c=buf[i])==state) state+=19;
			    else if (state!=132) state = (c==37) ? 56 : 37;
			    else if (state=37, (unsigned int)(c-=40) >= GP_N_CMDK) LOG("gp_key=%d,???", c);
			    else wrk_cmd_l((char*)gp_cmd_key[c],-1,0); }}
int w_main(int ac, char ** av) {
	if ((samp2gplot = (double*)map_wdir_shm(',', 64*GP_SAMPSIZ, 3))==MAP_FAILED) return perror("mmap"), 1;
	int r = gp_start(); LOG_E("gp_start", r);
	if (ac>1) ((r=gp_tlog_read(av[1]))<0) ? LOG_E(av[1],r) : gp_plot(1);
	fd_set rset; while (1) {
		FD_ZERO(&rset); FD_SET(0, &rset); FD_SET(gp_outpipe, &rset);
		r = select(gp_outpipe+1, &rset, 0, 0, NULL); if (r<0) perror("select");
		if (FD_ISSET(gp_outpipe, &rset)) gp_out();
		if (FD_ISSET(0, &rset) && (r=ib_read(&wrk_cb))<0 && (LOG_E("cmd/read",r), r>-3)) wbye(r&1); 
	}}

/*******************************/

int usage(const char *s) { return fprintf(stderr, "lf.bb: name \"%s\"%s", s,
		" unknown, valid names are lf.acv, lf.io, lf.lic, lf.con, lf.tld, lf.qplay\n"), 1; }

int main(int ac, char ** av) {
	const char *q, *s = *av;
	for (q=s; *s; s++) if (*s=='/') q = s+1;
	if (memcmp(q, "lf.", 3)) return usage(q);
	switch(q[3]) {
		case 'c': return c_main();
		case 'e': return e_main();
		case 'a': return a_main(ac, av);
		case 'i': return i_main(ac, av);
		case 'q': return q_main(ac, av);
		case 't': return t_main(ac, av);
		case 'w': return w_main(ac, av);
		case 'l': return execlp("less", "less", getenv("LF_LIC"), NULL);
		default: return usage(q); }}
