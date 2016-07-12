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

int c_main() { /* console */
	const char * cnm = tpipe_name('C');
	FILE * f = fopen(cnm, "a"); if (!f) return perror(cnm), sleep(3), 1;
	fprintf(f, "_c/proc/%d/fd/1\n", getpid()); fclose(f); while(1) sleep(60);  }

/**************** log i/o ****************************************************/

volatile int chld_flg = 0;
static int con_stat = 0, msg_fd = -1, errtemp = 0, killer_fd = 0, der_alte = -1, wrkdir_l = 0,
	   shutdn_flg = 8, cl_sec = 0, cl_usec = 0, selwt = 250000, gui_stat = 0, rs_ts = 0, dflg = 0;
static char xterm_nbuf[256];
static const char *xterm_name, *wrkdir;
static unsigned int ino_bv = 0u;
static int ino_epid[32], ino_oid[32], ino_fd = -1;
static char ino_inm[80];

static void LOG(const char * fmt, ...), bye(int kf);

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

static void del_n_log(const char *s, int f){ if ((f?rmdir:unlink)(s)<0) LOG("del %s: %s",s,strerror(errno));}
static void rm_r(const char * dn) {
	DIR * dir = opendir(dn); if (!dir) return LOG("rm_r: opendir: %s\n", strerror(errno));
	struct dirent * ent;  
	int l0 = strlen(dn); char buf[l0+1026]; memcpy(buf, dn, l0); buf[l0++] = '/'; buf[1024]=buf[1025]=0;
	while ((ent = readdir(dir))) {
		const char * s = ent->d_name; 
		if (*s=='.' && (!s[1] || (s[1]=='.'&&!s[2]))) continue;
		strncpy(buf+l0, s, 1024); del_n_log(buf, 0);
	}
	buf[l0] = 0; del_n_log(buf,1);
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
	ino_inm[l0+3] = ino_inm[l0+43] = 0; rm_r(ino_inm); rm_r(ino_inm+40);
        char buf[l0+16]; memcpy(buf, s0, l0); buf[l0] = '/'; buf[l0+2] = 0;
        for (s=FIFO_LIST; *s; s++) buf[l0+1] = *s, del_n_log(buf, 0);
        memcpy(buf+l0+1, "killer-file", 12); del_n_log(buf, 0);
        buf[l0] = 0; del_n_log(buf, 1);
}

static void close_all() {
	struct rlimit rlim; rlim.rlim_cur=9999; getrlimit(RLIMIT_NOFILE, &rlim);
	int i,n; for (i=0, n = (int)rlim.rlim_cur; i<n; i++) close(i); }

/*** output buf ***************/
static char logbuf[16384];
static int log_ix = 0, log_ox = 0, log_flg = 0;

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
static const char * ino_nm(int t, int oi) { char*q = ino_inm+40*t; h5f(q+wrkdir_l+4,oi); return q; }
static int ino_find(int *p, int id) { BVFOR_JMC(ino_bv) if (p[j]==id) return (int)j;   return -1; }
static void ino_bye(int j) { if (dflg&2) LOG("ino_bye: %d", j);
			     ino_bv &= ~(1u<<j); unlink(ino_nm(0,j)); unlink(ino_nm(1,j)); }

static void ino_ini() {
	int l = wrkdir_l;
	memcpy(ino_inm,    wrkdir, l); memcpy(ino_inm   +l, "/ds/00000.txt",14);
	memcpy(ino_inm+40, wrkdir, l); memcpy(ino_inm+40+l, "/d0/00000.txt",14);
	ino_inm[l+3] = ino_inm[l+43] = 0; 
	if (mkdir(ino_inm,0700)<0 || mkdir(ino_inm+40,0700)<0 || (ino_fd=inotify_init1(IN_CLOEXEC))<0 ||
	    inotify_add_watch(ino_fd, ino_inm, IN_CLOSE_WRITE)<0) 
		return ino_fd=-1, LOG("ino_ini: %s", strerror(errno));
	ino_inm[l+3] = ino_inm[l+43] = '/';
}

#define SERR(X) ((void) write(1, "_E" #X "\n", 5))
#define SERRC(X) (close(fd), (void) write(1, "_E" #X "\n", 5))
static void ino_add(int id, char * txt, int len) {
	if (ino_fd<0) return SERR(15);
	if (ino_find(ino_oid, id)>=0) return SERR(16);
	BVFOR_JMC(~ino_bv) goto found;
	return LOG("dsc-edit table full (32)");
found:; const char *nm0 = ino_nm(0, id), *nm1 = ino_nm(1, id);
	int i, fd = creat(nm1, 0600); if (fd<0) goto err2;
	for (i=0; i<len; i++) if (txt[i]==36) txt[i]=10; if (txt[len-1]!=10) txt[len++] = 10;
	txt-=7; len+=7; memcpy(txt, "#DSCR: ", 7);
	if (write(fd, txt, len)<0) goto err1; else close(fd);
	if (rename(nm1, nm0)<0 || creat(nm1, 0600)<0) goto err2;
	if (write(fd, txt, len)<0) goto err1; else close(fd);
	ino_bv |= (1u<<j); ino_oid[j] = id; 
	ino_epid[j] = launch("xedit", "!(TT", nm0, NULL);
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

/*** cmd **********************/
static int start_con() {
	static const char * sh = 0;
	return -((!sh && !(sh = getenv("LF_CON"))) || launch(xterm_name, "!)x1", "-e", sh, (char*)0)<0); }

static void con_started(char *s, int n) {
	if (memcmp(s, "_c/proc/", 8)) return LOG("_: invalid path\"%s\"", s);
	int pid = atoi(s+8); if (pid<1 || pid>65535) return LOG("/: invalid pid %d", pid);
	if (pid==con_stat) return LOG("/: already started, pid %d", pid); else con_end();
	if ((msg_fd = open(s+2, O_WRONLY|O_APPEND|O_NOCTTY))<0) return LOG("\"%s\": %s",s+2,strerror(errno));
	con_stat = pid; log_sn(0, -2, 0); --s[n-1]; s[n]=10; write(1, s, n+1);
}

static void cmd_l(char *s, int n, int src) { s[n] = 0; switch(*s) {
	case 'c': if (start_con()<0) LOG("start_con failed"); return;
	case '_': return con_started(s, n);
	case 'Z': return con_end();
	case 'f': return log_sn(0, -1, 0);
	case 'r': return (void) (rs_ts = time(NULL));
	case 'x': return (n>255) ? LOG("x: str too long") 
		  	         : (void) (memcpy(xterm_nbuf, s+1, n), xterm_name=xterm_nbuf);
	case 'h': return (n>6 && s[6]=='.') ? ino_add(atoi_h(s+1), s+7, n-7) : LOG("'h': parse error");
	case 'd': return (void) (dflg = atoi_h(s+1));
	default: LOG("unknown cmd%c 0x%x (%s)", 48+src, *s, s); return;
}}

/*** input buf ****************/
typedef void (*lfun)(char *, int, int);
typedef struct { int fd, siz, cont, arg; lfun lf; char * p; } ibuf;

#define N_IBUF 9
char bigbuf[0x8000];
static ibuf ib[9] = {{ 0, 3072, 0, 0,   &cmd_l, bigbuf},
		     {-1, 4096, 0, '.', &log_l, bigbuf+0x0c00},
		     {-1, 1024, 0, 'C', &cmd_l, bigbuf+0x1c00},
		     {-1, 4096, 0, 'p', &log_l, bigbuf+0x2000},
		     {-1, 4096, 0, 't', &log_l, bigbuf+0x3000},
		     {-1, 4096, 0, 'u', &log_l, bigbuf+0x4000},
		     {-1, 4096, 0, 'x', &log_l, bigbuf+0x5000},
		     {-1, 4096, 0, 'k', &log_l, bigbuf+0x6000},
		     {-1, 4096, 0, 's', &log_l, bigbuf+0x7000}};

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
		   if (ag<65) { close(b->fd); b->fd = -1; shutdn(2+2*!!ag); return; }
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
	hello(); 
	if (!(wrkdir=getenv("LF_TMPDIR")) || (wrkdir_l=strlen(wrkdir))>20) fail("workdir");
	ino_ini();
	if (!(xterm_name = getenv("LF_XTERM"))) xterm_name = "xterm";
	ib_init(   ac<2 ? -1 : qh4r(*(int*)av[1]));
	killer_fd = (ac<3||*av[2]<48) ? 0 : qh4r(*(int*)av[2]);
	der_alte = ac<4 ? getppid() : qh4r(*(int*)av[3]);
	int i, k;
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
		FOR_IB if (FD_ISSET(k, &rset)) ib_read(ib+i);
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
	else k&=4095, sprintf(to, "%d"+(k>999), k);
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

/*******************************/

int usage(const char *s) { return fprintf(stderr, "lf.bb: name \"%s\"%s", s,
		" unknown, valid names are lf.acv, lf.io, lf.lic, lf.con, lf.tld, lf.qplay\n"), 1; }

int main(int ac, char ** av) {
	const char *q, *s = *av;
	for (q=s; *s; s++) if (*s=='/') q = s+1;
	if (memcmp(q, "lf.", 3)) return usage(q);
	switch(q[3]) {
		case 'c': return c_main();
		case 'a': return a_main(ac, av);
		case 'i': return i_main(ac, av);
		case 'q': return q_main(ac, av);
		case 't': return t_main(ac, av);
		case 'l': return execlp("less", "less", getenv("LF_LIC"), NULL);
		default: return usage(q); }}

