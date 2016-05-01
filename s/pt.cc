#include "pt.h"
#include "uc0.h"
#include "util.h"
#include "glob.h"
#include "cfgtab.inc"

#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/wait.h>
#include <sys/stat.h>

volatile int pt_chld_flg = 0;
int pt_cp_i2m = -1;

typedef struct { int pid, t, st; pt_wfun f; } pt_tab_t;

static int pt_cp_m2i = -1, pt_constat = 0, pt_nullfd = -1;
static volatile pt_tab_t pt_tab[PT_COUNT];
static const char * pt_logf_name;
static int pt_acv_cur = 0, pt_acv_cw = 0, *pt_con_pfd = 0;

static void log_pidstat(const char *s, int pid, int stat) {
        if (WIFEXITED(stat)) log("%s %d exited (%d)", s, pid, WEXITSTATUS(stat));
        else if (WIFSIGNALED(stat)) log("%s %d killed (%d)", s, pid, WTERMSIG(stat));
        else log("something happened to %s %d (0x%x)", s, pid, stat); }

static void pt_dlog(const char * s) {
        const char * nm = pt_logf_name;
        int fd = open(nm, O_WRONLY|O_APPEND);
	switch(fd) {
		case 1:  dup2(1, 2); break;
		case 2:  dup2(2, 1); break;
		default: if (fd<0) return;
			 dup2(fd, 1); dup2(fd, 2); close(fd); break;
	}
        if (s) write(2, s, strlen(s));
}

static void pt_sighup(int _) { pt_chld_flg |= 0x40000000; signal(SIGHUP, &pt_sighup); }
static void pt_sigint(int _) { pt_chld_flg |= 0x20000000; signal(SIGINT, &pt_sigint); }

static void pt_sigchld(int _) { while(1) {
        int i, st, pid = waitpid(-1, &st, WNOHANG); if (pid<=0) return (void) signal(SIGCHLD, &pt_sigchld);
	for (i=0; i<PT_COUNT; i++) if (pt_tab[i].pid==pid) goto found;
	log_pidstat("unreg. child", pid, st); continue;
found:  if (!i) pt_dlog("\nERROR: ioprc exited, some logs will be lost...\n");
	pt_tab[i].st = st;
	log_pidstat(PT_STR+8*i, pid, st); pt_chld_flg |= (1<<i); signal(SIGCHLD, &pt_sigchld);
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
	fflush(stdout); fflush(stderr);
	close(1), close(2), dup2(pf3, 1), dup2(pf3, 2), close(pf3);
	set_fd(&pt_cp_m2i, pf1); set_fd(&pt_cp_i2m, pf2);
	return pid;
}

static int acv_dead(int pid, int stat, int td) {
	int e, c = 0;
	if (!WIFEXITED(stat)) e = WIFSIGNALED(stat) ? PTE_ACVSIG : PTE_ACVHMM;
	else (stat=WEXITSTATUS(stat)) ? (stat==90 ? (e=PTE_ACVZERO,c=1) : (e=PTE_ACVERR)) : (e=0,c=1);
	if (e) gui_errq_add(e); if (c && pt_acv_cw) gui_closewin(ACV_WIN(pt_acv_cur)); return 0;
}

static int con_started(const char *s) {
	if (memcmp(s, "/proc/", 6) || (unsigned int)(s[6]-49)>8u) return GCE_PARSE;
	int fd, pid = atoi(s+6); if (pid&-65536) return GCE_PARSE;
	if ((fd=open(s,O_RDONLY|O_NOCTTY))<0 || dup2(fd,0)<0) return pt_constat = 0, EEE_ERRNO;
	close(fd); pt_constat = pid; *pt_con_pfd = 0; return 0;
}

//////// export ////////////////////////////////////////////

void pt_con_fd_ptr(int *p) { pt_con_pfd = p; }
int pt_iocmd_sn(const char *s, int n) { int r = write(pt_cp_m2i,s,n); return (r==n) ? 0 : EEE_ERRNO+(r>=0); }

int pt_iocmd(char *s) {
	if ((pt_chld_flg & (1<<PT_IOP))) return EEE_NOPRC;
	int r, l = strlen(s), c = s[l]; s[l] = 10;
	return r = pt_iocmd_sn(s, l+1), s[l] = c, r; }

int pt_con_op(const char *s) { 
	if (debug_flags&DFLG_PT) log("pt_con_op(%s)", s);
	if (*s!='-') return con_started(s);
	int k, l;
	switch(s[1]) {
		case '1': if (pt_constat) { if (pt_constat>0) { k = -3; goto kill; }
				  	    gui_errq_add(EEE_STATE, "con/-1"); }
			  return pt_constat = -1, pt_iocmd_sn("c\n", 2);
		case '2': case '3': if (pt_constat<1) return EEE_STATE; k = 48 - *s; goto kill;
		case '4': if (pt_constat>0) kill(pt_constat, 9), k = pt_constat = 0;
			  else if ((k=-pt_constat-2)&~1) return EEE_STATE;
			  else pt_constat = -k;
			  if (dup2(pt_nullfd, 0)<0) gui_errq_add(EEE_ERRNO, "con/-4");
			  return *pt_con_pfd = -1, pt_iocmd_sn("Z\nc\n", 2+2*k);
		default:  return EEE_RANGE;
	}
kill:	return kill(pt_constat, 9)<0 ? EEE_ERRNO : (pt_constat = k, 0);
}

void pt_reg(int ix, int pid, pt_wfun fun) {
	if ((unsigned int)ix >= (unsigned int)PT_COUNT) return log("BUG: pt_reg: ix=%d", ix);
        pt_tab[ix].pid = pid; pt_tab[ix].f = fun;
}

void pt_chld_act() {
	int k = pt_chld_flg >> 28; if (k) log("exiting on SIG%s", "INT\0HUP"+(k&4)), bye(k);
	BVFOR_JM(pt_chld_flg) {
		volatile pt_tab_t * q = pt_tab+j;
		int t = time(0), pid = (*q->f)(q->pid, q->st, t - q->t);
		q->t = (q->pid = pid) > 0 ? t : 0;
	} 
	pt_chld_flg = 0; }

int pt_acv_op(int id, int opw, const char *a1, const char *a2) {
	static const char * acv_bin = 0; if (!acv_bin && !(acv_bin=getenv("LF_BB"))) acv_bin="lf.bb";
	int sdl, tdl, pid, op = opw & 255, nwf = opw & 256;
	const char *sdir, *tdir, *src, *trg;
	if ((sdl = CFG_AO_DIR.i )) sdir = CFG_AO_DIR.s;  else sdl = tmp_dir_len, sdir = tmp_dir;
	if ((tdl = CFG_WAV_DIR.i)) tdir = CFG_WAV_DIR.s; else tdl = hsh_dir_len, tdir = hsh_dir;
	switch(op) {
		case 0xf: return nwf ? 0 : (gui_closewin(ACV_WIN(id)), 0);
		case 0xd: pid = launch("/bin/rm", "!(kk", au_file_name(sdir,sdl,id,0,0,"a20"), (char*)0);
			  pt_acv_cw = !nwf; goto ret;
		case 0xe: log("pt_acv_op: 0x0e -- why here?"); return GCE_ROUTE;
		default:  if (op>7) return GCE_PARSE; else break;
	}
	src = au_file_name(sdir, sdl, id, 0,  0,  "a20");
	trg = au_file_name(tdir, tdl, id, a1, a2, "wav\0flac" + (op&4));
	pt_acv_cw = !nwf && !(~op&3);
	switch (op&3) {
		case 0: pid = launch(acv_bin, "(kk", "lf.acv",       src, trg,         (char*)0); break;
		case 1: pid = launch(acv_bin, "(kk", "lf.acv",       src, trg, a1,     (char*)0); break;
		case 2: pid = launch(acv_bin, "(kk", "lf.acv",       src, trg, a1, a2, (char*)0); break;
		case 3: pid = launch(acv_bin, "(kk", "lf.acv", "-r", src, trg,         (char*)0); break;
	}
ret:	return (pid<0) ? EEE_ERRNO : (pt_reg(PT_ACV, pid, &acv_dead), pt_acv_cur = id, 0);
}

int pt_show_lic() { return launch(CFG_XTERM.s, "!)((", "-e", getenv("LF_LICB"), (char*)0); }

int pt_kill_pa(int flg) { return (launch("killall","!(ss","-v",(flg&2)?"-9":"-15","pulseaudio", (char*)0)<0) ?
					EEE_ERRNO : 0; }

void pt_init() {
	if (!(pt_logf_name = getenv("LF_LOG"))) pt_logf_name = "/tmp/lf.noenv.log",
		pt_dlog("no LF_LOG defined, this was a short run...\n"), exit(1);
        pt_dlog("...\n");
	signal(SIGHUP, &pt_sighup); signal(SIGINT, &pt_sigint); 
        const char * kfn = getenv("LF_KILLER"); if (!kfn) log("no killer file");
	if ((pt_nullfd = open("/dev/null", O_RDONLY)) < 0 ||
	    (pt_nullfd ? dup2(pt_nullfd, 0) : (pt_nullfd = dup(0))) < 0) perror("nullfd"), exit(1);
	if (kfn && (killer_fd = open(kfn, O_RDONLY))<0) log("\"%s\": %s", kfn, strerror(errno));
	if ((wrk_dir=getenv("LF_TMPDIR")))  wrk_dir_len=strlen(wrk_dir); else wrk_dir="/tmp", wrk_dir_len=4;
	if ((tmp_dir=getenv("LF_TMPROOT"))) tmp_dir_len=strlen(tmp_dir); else tmp_dir="/tmp", tmp_dir_len=4;
	if ((usr_dir=getenv("LF_USERDIR"))) usr_dir_len=strlen(usr_dir); else usr_dir="/tmp", usr_dir_len=4;
	if ((hsh_dir=getenv("HOME")))       hsh_dir_len=strlen(hsh_dir); else hsh_dir="/tmp", hsh_dir_len=4;
	char *s1 = (char*)malloc(usr_dir_len+10); memcpy(s1, usr_dir, usr_dir_len);
	char *s2 = (char*)malloc(usr_dir_len+13); memcpy(s2, usr_dir, usr_dir_len);
	memcpy(s1+usr_dir_len,  "/__asv.lf",    10); autosave_name   = s1;
	memcpy(s2+usr_dir_len,  "/__asv--x.lf", 13); autosave_name_x = s2;
	signal(SIGCHLD, &pt_sigchld);
	for (const char* s = FIFO_LIST; *s; s++) if (mkfifo(tpipe_name(*s), 0600)<0)
		log("FATAL: failed mkfifo '%c': %s", *s, strerror(errno)), exit(1);
	if (io_start(0) < 0) log("FATAL: io_start failed"), exit(1);
}
