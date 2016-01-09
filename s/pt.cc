#include "pt.h"
#include "uc0.h"
#include "util.h"
#include "glob.h"

#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/wait.h>
#include <sys/stat.h>

volatile int pt_chld_flg = 0;
int pt_cp_i2m = -1;

typedef struct { int pid, t, st; pt_wfun f; } pt_tab_t;

static int pt_cp_m2i = -1, pt_constat = 0;
static volatile pt_tab_t pt_tab[PT_COUNT];
static const char * pt_logf_name;
static int pt_acv_cur = 0;

static void log_pidstat(const char *s, int pid, int stat) {
        if (WIFEXITED(stat)) log("%s %d exited (%d)", s, pid, WEXITSTATUS(stat));
        else if (WIFSIGNALED(stat)) log("%s %d killed (%d)", s, pid, WTERMSIG(stat));
        else log("something happened to %s %d (0x%x)", s, pid, stat); }

static void pt_dlog(int flg, const char * s) {
        const char * nm = pt_logf_name;
        int fd = (flg&1) ? creat(nm, 0600) : open(nm, O_WRONLY|O_APPEND);
        if (fd<0) return; else close(1), close(2), dup2(fd, 1), dup2(fd, 2);
        if (s) write(fd, s, strlen(s)); close(fd);
}

static void pt_sigchld(int _) { while(1) {
        int i, st, pid = waitpid(-1, &st, WNOHANG); if (pid<=0) return (void) signal(SIGCHLD, &pt_sigchld);
	for (i=0; i<PT_COUNT; i++) if (pt_tab[i].pid==pid) goto found;
	log_pidstat("unreg. child", pid, st); continue;
found:  if (!i) pt_dlog(1, "\nERROR: ioprc exited, some logs will be lost...\n");
	pt_tab[i].st = st;
	log_pidstat(PT_STR+8*i, pid, st); pt_chld_flg |= (1<<i);
}}

static int io_start(int);
static int iop_dead(int pid, int stat, int td) {
	if (td<2) log("FATAL: ioprc exited again in < 2 seconds"), bye(1);
	gui_errq_add(PTE_IOCRASH);
	int pid2 = io_start(1); if (pid2<0) log("FATAL: ioprc restart failed"), bye(1);
	if (pt_constat>0) kill(pt_constat, 9); pt_constat = 0;
	return pid2;
}

static int io_start(int re) {
	static const char * exe = 0; if (!exe && !(exe=getenv("LF_IO"))) exe = "./lf.io";
	char a2[16], *aa = a2; *(int*)aa = killer_fd>0 ? qh4(killer_fd) : 33; a2[4] = 0;
	*(int*)(a2+8) = qh4(getpid()); a2[12] = 0;
	int pf1, pf2, pf3, pid = launch(exe, "!><+>", &pf1, &pf2, pt_logf_name, &pf3, a2, a2+8, (char*)0);
	if ((pid|pf1|pf2|pf3)<0) return -1;
	if (!re) pt_reg(PT_IOP, pid, &iop_dead);
	close(1), close(2), dup2(pf3, 1), dup2(pf3, 2), close(pf3);
	set_fd(&pt_cp_m2i, pf1); set_fd(&pt_cp_i2m, pf2);
	return pid;
}

static inline int iocmd_sn(const char *s, int n) {
	int r = write(pt_cp_m2i, s, n); return (r==n) ? 0 : EEE_ERRNO+(r>=0); }

static int acv_dead(int pid, int stat, int td) {
	int e, c = 0;
	if (!WIFEXITED(stat)) e = WIFSIGNALED(stat) ? PTE_ACVSIG : PTE_ACVHMM;
	else (stat=WEXITSTATUS(stat)) ? (stat==90 ? (e=PTE_ACVZERO,c=1) : (e=PTE_ACVERR)) : (e=0,c=1);
	if (e) gui_errq_add(e); if (c) gui_closewin(ACV_WIN(pt_acv_cur)); return 0;
}

//////// export ////////////////////////////////////////////

int pt_iocmd(char *s) {
	if ((pt_chld_flg & (1<<PT_IOP))) return EEE_NOPRC;
	int r, l = strlen(s), c = s[l]; s[l] = 10;
	return r = iocmd_sn(s, l+1), s[l] = c, r; }

int pt_con_op(int x) { switch(x) {
	case -2: x = pt_constat; pt_constat = 0;
		 if (x!=-4) return (x>0||x==-2) ? 0 : EEE_STATE;
	case -1: if (!pt_constat) return pt_constat = -1, iocmd_sn("c\n", 2);
		 else return pt_constat > 0 ? (kill(pt_constat, 9), pt_constat = -4, 0) : EEE_STATE;
	case -3: if (pt_constat>0) return kill(pt_constat, 9), pt_constat = -2, 0;
		 return pt_constat ? EEE_STATE : EEE_NOEFF;
	default: if (x<=0) return EEE_STATE;
		 return pt_constat = x, iocmd_sn("C\n", 2);
}}

void pt_reg(int ix, int pid, pt_wfun fun) {
	if ((unsigned int)ix >= (unsigned int)PT_COUNT) return log("BUG: pt_reg: ix=%d", ix);
        pt_tab[ix].pid = pid; pt_tab[ix].f = fun;
}

void pt_chld_act() { BVFOR_JM(pt_chld_flg) {
	volatile pt_tab_t * q = pt_tab+j;
	int t = time(0), pid = (*q->f)(q->pid, q->st, t - q->t);
	q->t = (q->pid = pid) > 0 ? t : 0;
} pt_chld_flg = 0; }

const char * pt_acv_nm(int id, int j) { 
	static char *buf2, *buf3, *buf = 0; // tmp/8.4 tmp/8.4 home/8.4
	static int tlen, hlen;
	if (!buf) {
		const char *s1 = getenv("LF_TMPROOT"), *s2 = getenv("HOME");
		tlen = s1 ? strlen(s1) : (s1="/tmp", 4);
		hlen = s2 ? strlen(s2) : (s2="/tmp", 4);
		buf = (char*)malloc(2*tlen + hlen + 45); buf2 = buf+tlen+15; buf3 = buf2+tlen+15;
		memcpy(buf , s1, tlen); buf [tlen] = '/';
		memcpy(buf2, s1, tlen); buf2[tlen] = '/';
		memcpy(buf3, s2, hlen); buf3[hlen] = '/'; tlen++, hlen++;
	}
	const char * ext = (j&1) ? "flac" : "wav";
	switch(j&6) {
		case 0: return sprintf(buf+tlen, "%08x.a20", id), buf;
		case 2: return sprintf(buf2+tlen, "%08x.%s", id, ext), buf2+tlen;
		case 6: return sprintf(buf2+tlen, "%08x.%s", id, ext), buf2;
		case 4: return sprintf(buf3+hlen, "%08x.%s", id, ext), buf3;
		default: return "WTF";
	}}

int pt_acv_op(int id, int op) {
	static const char * acv_bin = 0; if (!acv_bin && !(acv_bin=getenv("LF_ACV"))) acv_bin="lf.acv";
	int pid = 0;
	switch(op) {
		case 0:   gui_closewin(ACV_WIN(id)); return 0;
		case 0xd: pid = launch("/bin/rm", "!(kk", pt_acv_nm(id, 0), (char*)0); break;
		default:  if (op<2 || op>7) return BXE_PARSE;
			  pid = launch(acv_bin, "!(kk", "-r", pt_acv_nm(id,0), pt_acv_nm(id,op), (char*)0); break;
	}
	return  (pid<0) ? EEE_ERRNO : (pt_reg(PT_ACV, pid, &acv_dead), pt_acv_cur = id, 0);
}
			  
int pt_show_lic() { return launch("xterm", "!)((", "-e", getenv("LF_LICB"), (char*)0); }

void pt_init() {
        const char * kfn = getenv("LF_KILLER");
        int fd = -1, k = kfn ? open(kfn, O_RDONLY) : -1;
        killer_fd = k ? k : dup(0);
	if (!(pt_logf_name = getenv("LF_LOG"))) pt_logf_name = "/tmp/lf.noenv.log";
        pt_dlog(1, "(log start)\n");
	signal(SIGCHLD, &pt_sigchld);
        close(0); if ((fd=open("/dev/null", O_RDONLY))>0) dup2(fd,0), close(fd);
        for (const char* s ="ptuxmkACT"; *s; s++) if (mkfifo(tpipe_name(*s), 0600)<0)
                log("FATAL: failed mkfifo '%c': %s", *s, strerror(errno)), bye(1);
	if (io_start(0) < 0) log("FATAL: io_start failed"), bye(1);
}
