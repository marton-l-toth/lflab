#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/resource.h>

#include "pt.h"
#include "uc0.h"
#include "util.h"
#include "util2.h"
#include "glob.h"

#ifdef __SSE2__
#include <xmmintrin.h>
#define FPU_INI _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); memcpy(pt_he_buf,"sse2:Y",6);
#else
#define FPU_INI memcpy(pt_he_buf,"sse2:n",6);
#endif

int cmd_report2(int c, const char *s); // cmd.cc
static int osb_hash(const char *s) { // od|sed|bc
	int r=0; while (*s) r *= (*s>99)?1000:100, r += *s, s++;   return r&0x3fffffff; }
#define QWE_DEFINE_CFGTAB
#include "cfgtab.inc"

static char pt_he_buf[40];
volatile int pt_sig_flg = 0;
int pt_nullfd = -1, pt_qenv_l[32], pt_errtemp = 0, pt_xapp_bv[N_XAPP];
const char *pt_hello = pt_he_buf, *pt_qenv_v[32];
// e: -1:no_ini(1) -2:inisv.e.(2) -4:optarg_missing 1<<30(+n):argn 1<<29(+n):ini/ln
static int ee_v[64], ee_n=0, ee_flg=0, ee_errno = 0;

typedef struct { int pid; char s[12]; } pt_ptent_t;
#define pt_cp_m2i  (pt_cp_m2iw[0])
#define pt_wrk_m2w (pt_cp_m2iw[1])
static int pt_cp_m2iw[2] = {-1, -1}, pt_constat = 0, pt_pid;
static pt_ptent_t pt_ptab[32];
static unsigned int pt_ptab_bv = 0u;
static int pt_acv_cur = 0, pt_acv_cw = 0, *pt_con_pfd = 0, *pt_io_pfd = 0;
static uid_t pt_uid, pt_gid;
static int ppath_n = 0, ppath_maxlen = 0, *ppath_len;
static const char** ppath = 0;
static struct timeval tv_zero;
static int pt_wrk_pid = 0, pt_samp_bits = 0, *pt_wrk_pfd = 0;
static double * pt_samp_buf;

static void qfail(const char * fmt, ...) {
	va_list ap; va_start(ap, fmt);
	char buf[1024]; memcpy(buf, "lflab: quick FAIL: ", 19);
	int l = 19+vsnprintf(buf+19, 1000, fmt, ap); va_end(ap); buf[l++] = 10;
	int fd = open(QENV('>'), O_WRONLY|O_APPEND|O_CREAT, 0644);
	if (fd>=0) write(fd, buf, l), close(fd);
	write(2, buf, l); exit(1);
}

static void ee_add(int ec) { if (ec<0) return (void)(ee_flg |= -ec); else ee_flg |= (ec&(7<<28));
			     if (ee_n<64) ee_v[ee_n++]=ec; }
static inline void qe_set2(int c, const char *s, int l) { QENV(c)=s; QENVL(c)=l; }
static inline void qe_set2a(int c, const char *s, int l) { 
	QENVL(c) = l++; char *q = nf_alloc(l); memcpy(q, s, l); QENV(c) = q; }
static inline void qe_len1(int c) { QENVL(c) = strlen(QENV(c)); }
static void qe_len(const char *ls) { while (*ls) qe_len1(*ls), ls++; }
static int stat2(const char *nm, int op) {
	struct stat buf; if (stat(nm, &buf)<0) return 0;
	int k = buf.st_mode; if (((op&32) ? !S_ISREG(k) : !S_ISDIR(k))) return 0;
	int flg = k & (7 | (070&-(pt_gid==buf.st_gid)) | (0700&-(pt_uid==buf.st_uid)));
	switch(op|32) { case 'r': return !!(flg&0444);
			case 'w': return !!(flg&0222);
			case 'x': return !!(flg&0111);
			default: return 0; }}

static void mk_ppath() {
        const char * p0 = getenv("PATH");
        if (!p0) { ppath_n = 2; ppath_maxlen = 8;
                   ppath = (const char**)malloc(2*sizeof(char*)); ppath_len = (int*)malloc(2*sizeof(int));
                   ppath[0]="/bin"; ppath[1]="/usr/bin"; ppath_len[0] = 4; ppath_len[1] = 8; return; }
        int i0, i, ltot = strlen(p0), np = 1;
        char * p1 = (char*)malloc(ltot+1); memcpy(p1, p0, ltot); p1[ltot]=':';
        for (i=0; i<ltot; i++) if (p1[i]==':') ++np;
        ppath = (const char**)malloc(np*sizeof(char*)); ppath_len = (int*)malloc(np*sizeof(int));
        for (i0=i=0; i<=ltot; i++) { if (p1[i]==':') {
                int l = i-i0; if (l>ppath_maxlen) ppath_maxlen = i;
                if (l>0) ppath[ppath_n] = p1+i0, ppath_len[ppath_n] = l, p1[i]=0, ++ppath_n;
                i0 = i+1; }}}

static int pfind(const char * nm) {
        if (!ppath) mk_ppath();
        int i, l = strlen(nm);
        char buf[l+ppath_maxlen+2];
        buf[ppath_maxlen]='/'; memcpy(buf+ppath_maxlen+1, nm, l+1);
        for (i=0; i<ppath_n; i++) {
                int k = ppath_len[i], j = ppath_maxlen - k;
                memcpy(buf+j, ppath[i], k); if (stat2(buf+j, 'x')) return 1; }
        return 0;
}

static const char * pfind_ls(const char ** pp) {
	while (*pp) if (pfind(*pp)) return *pp; else pp++;   return "(none)"; }

static void qe_dump() { for (int i=59; i<91; i++) if (QENV(i))
	fprintf(stderr, "QENV('%c')=\"%s\", len=%d\n", i, QENV(i), QENVL(i)); }

static int set_tmp(const char *s, int l) {
	char buf[256]; return (stat2(s,'W') && readlink(s,buf,255)<0) ? (qe_set2('t', s, l),0) : -1; }

static char wdir[24];
static void qe_ini() {
	static const char * self = "/proc/self/exe";
	static const char * rel[] = { "uh/.lflab", 
		">u/lf.log", "ju/lf.ini", "=u/lf.tlog", "au/__asv.lf", "xu/__asv--x.lf", "zw/killer-file",
		"bp/lf.bb", "qp/lf.pump", "ep/ex.lf", "cp/COPYING", "vp/lf.lic", "gp/lf.gui", "mp/lf.bin",
		"?p/help.txt", "kp/lf.con", "rp/lf.gtk.rc.def", "@p/lf.ed" };
	static const char * sete[] = { "kLF_CON", "?LF_HLP", "zLF_KILLER", "cLF_LIC", ">LF_LOG", "=LF_TLOG",
		"wLF_TMPDIR", "tLF_TMPROOT", "uLF_USERDIR", "rLF_GTKRC", "mLF_BIN", "pLF_DIR", "@LF_ED"};
	static const int nrel = sizeof(rel)/sizeof(char*), nsete = sizeof(sete)/sizeof(char*);

	QENV('>') = "/tmp/lflab.quick-fail.log";
	if (!(pt_gid=getegid(), pt_uid=geteuid())) qfail("it is a bad idea to run lflab as root");
	if ((QENV('h')=getenv("HOME"))) qe_len1('h'); else qe_set2('h', "/tmp", 4);
	if (set_tmp("/run/shm",8) && set_tmp("/dev/shm",8) && set_tmp("/tmp",4)) qfail("no tmpdir found");
	sprintf(pt_he_buf+6, ",%d", pt_pid = getpid());
	qe_set2('w', wdir, sprintf(wdir, "%s/lf.%s", QENV('t'), pt_he_buf+7));

	char pkgs[1024]; int pkl = readlink(self, pkgs, 1023);
	if ((unsigned int)pkl>999u) qfail("%s: %s", self, pkl<0 ? strerror(errno) : "linktrg too long");
	if (pkl<2 || pkgs[0]!='/') pkgs[pkl]=0, qfail("%s: \"%s\": what???", self, pkgs);
 	for (--pkl; pkgs[pkl]!='/'; --pkl);   pkgs[pkl]=0; qe_set2a('p', pkgs, pkl);

	for (int i=0; i<nrel; i++) {
		const char *s=rel[i]; char*q; int l1=QENVL(s[1]), l2=strlen(s+2), l12 = QENVL(*s) = l1+l2;
		memcpy(q=nf_alloc(l12+1), QENV(s[1]), l1); memcpy(q+l1, s+2, l2+1);     QENV( *s) = q; }
	for (int i=0; i<nsete; i++) setenv(sete[i]+1, QENV(sete[i][0]), 1);

	if (mkdir(QENV('w'), 0700)<0) qfail("mkdir(%s): %s", QENV('w'), strerror(errno));
	if (!stat2(QENV('u'), 'W') && (glob_flg|=GLF_HITHERE, (mkdir(QENV('u'),0700)<0)))
		qfail("mkdir(%s): %s", QENV('u'), strerror(errno));
	int kff=0, kf = creat(QENV('z'), 0600); if (kf>2) return (void)(killer_fd = kf);
	if (kf<0) qfail("creat(%s): %s", QENV('z'), strerror(errno));
	do { if (kff |= 1<<kf, (kf=dup(kf))<0) qfail("dup(%s): %s",QENV('z'),strerror(errno)); } while(kf<3);
	for (int i=0; i<3; i++) if (kff&(1<<i)) close(i);     killer_fd = kf;
}

void cfg_setint(cfg_ent *p, int k) { p->i = ivlim(k, p->i_m, p->i_M); }
int cfg_setstr(cfg_ent *p, const char *s) {
	int l01 = p->i_M, l0 = l01>>16, l1 = l01&65535, l = strlen(s);
        return l0<=l && l<=l1 && (memcpy(p->s, s, (p->i=l)+1), 1); }
void cfg_setx(cfg_ent *p, const char *s) {if (p->i_m==0x7fffffff) cfg_setstr(p,s); else cfg_setint(p,atoi(s));}

static void cfg0() {
	const char *s;  char *q = cfg_sbuf; 
	for (cfg_ent *p = cfg_tab+1; p->s_vd; p++) { if (p->i_m==0x7fffffff) {
		p->s=q; q+=(p->i_M&65535)+1;
		for (s = p->s_vd; *s; s++);   cfg_setstr(p, s+1);
	}}
	if (q-cfg_sbuf != sizeof(cfg_sbuf)) qfail("cfg0: diff=%d, siz=%d\n", q-cfg_sbuf, sizeof(cfg_sbuf));
}

static void il_parse(int ln, char * s) {
	int i,j,l = strlen(s); if (s[l-1]==10) s[--l] = 0;
	if (!memcmp(s,"LF_",3)) s+=3, l-=3;
	for (i=0; i<l; i++) if (s[i]=='=') goto eqfound;
	return ee_add((1<<29)+ln);
eqfound:s[i++]=0; if (s[i]=='"' && (++i<l) && s[l-1]=='"') s[--l]=0;
	if ((j=osb_find(s))<0) ee_add((1<<29)+ln); else cfg_setx(cfg_tab+j, s+i);
}

static const char * hlp_head = "usage: lflab <options> <lib>* <save_file>\n"
			       "    or lflab <options> <lib>* / <save_file>*\n"
			       "(<lib>s will be read as a collection of read-only objects)\n"
			       "options: (x:integer arg. s:string arg. [current-value|default|min..MAX])\n"
			       "-ni      skip lf.ini file (must be the first option)\n"
			       "-i     s read <s> instead of lf.ini (must be the first option)\n";
static void show_help(FILE *f) {
	fprintf(f, "%s", hlp_head);
	cfg_ent * ce; const char *he, *s;
	for (int i=1; he=cfg_help[i], (ce=cfg_tab+i)->s_vd; i++) {
		int ty = (ce->s_vd[0]=='_') ? ' ' : 'x'-5*(ce->i_m==0x7fffffff), mx=ce->i_M;
		fprintf(f, "-%s%c %s", he, ty, he+7); if (ty==' ') { fprintf(f,"\n"); continue; }
		if (ty=='s') { for (s=ce->s_vd; *s; s++);
			       fprintf(f, " [%s|%s|%d..%d]\n", ce->s, s+1, mx>>16, mx&65535); }
		else	     { fprintf(f, " [%d|%d|%d..%d]\n", ce->i, ce->i_d, ce->i_m, mx); }
}}

static void mklog() {
	int lfd = creat(QENV('>'), 0644); if (lfd<0) qfail("creat(%s): %s", QENV('>'), strerror(errno));
	struct tm * lt = localtime(&tv_zero.tv_sec);
	char lbuf[256];
	int lbl = snprintf(lbuf, 255,"-----lflab_%d.%02d----(%d)-----[%d.%02d.%02d/%02d:%02d:%02d.%06d]------",
				     v_major, v_minor, pt_pid, 1900+lt->tm_year, lt->tm_mon+1, lt->tm_mday,
				     lt->tm_hour, lt->tm_min, lt->tm_sec, (int)tv_zero.tv_usec);
	lbuf[lbl]=10; write(lfd, lbuf, lbl+1); close(lfd);
}

void cfg_export(cfg_ent *p) {
	char nm[64]; memcpy(nm  , "LF_", 3);
	int vbuf[3]; strcpy(nm+3, p->s_vd);
	if (p->i_m==0x7fffffff) return (void)setenv(nm, p->s, 1);
	if (p->i&~65535) vbuf[0]=qh4((p->i>>16)&65535), vbuf[1]=qh4(p->i&65535), vbuf[2]=0;
	else		 vbuf[0]=qh4(p->i), vbuf[1]=0;
	setenv(nm, (const char*)vbuf, 1);
}

int cfg_write(int lg) {
        backup(QENV('j'), 9);  int k = -1;  const char *v;
        FILE * f = fopen(QENV('j'), "w"); if (!f) goto err0;
        for (cfg_ent * p = cfg_tab+1; (v=p->s_vd); p++)
                if (*v!='_' && (k = (p->i_m==0x7fffffff) ? fprintf(f, "%s=\"%s\"\n", v, p->s)
							 : fprintf(f, "%s=%d\n",     v, p->i)) < 0) goto err1;
	if (fclose(f)) goto err0;
	if (lg) log("config written to \"%s\"", QENV('j'));    return 0;
err1:	fclose(f);
err0:	return lg ? (k?EEE_ERRNO:EEE_ZEROLEN) : (ee_add(-2), ee_errno = (k?errno:0), -1);
}

static const char ** cfg1(int ac, char**av) {
	int i, j, skip = 1, dfi = 1;
	const char **rv, *s, *inif = QENV('j');
	if (ac>1 && (s=av[1])[0]=='-') {
		if (s[1]=='n' && s[2]=='i' && !s[3]) skip=2, dfi=0, inif=0;
		else if (s[1]=='i' && !s[2]) (ac>2) ? (inif=av[2],dfi=0,skip=3) : (ee_add(-4),skip=2);
	}
	if (inif) { FILE * f = fopen(inif, "r");
		    if (!f) { ee_add(-1); dfi *= 2; }
		    else    { char buf[256]; for(int i=1; fgets(buf,255,f); i++) il_parse(i,buf); fclose(f); }}
	for (i=skip; i<ac; i++) {
		if (*(s=av[i])=='-') ++s; else break;
		if ((j = osb_find(s))<0) { ee_add((1<<30)+i); continue; }
		cfg_ent * ce = cfg_tab + j;
		if (ce->s_vd[0]=='_') { switch(ce->s_vd[1]) {
			case '1': ce->i = 1; continue;
			case 'H': show_help(stdout); exit(0);
			default: ee_add(-8); continue;
		}}
		if (i==ac-1) ee_add(-4); else cfg_setx(ce, av[i+1]), ++i;
	}
	debug_flags = CFG_DEBUG_FLG.i;
	for (int k=0; k<N_XAPP; k++) { cfg_ent * ent = &CFG_XTERM+k;
			               if (!*ent->s) cfg_setstr(ent, pfind_ls(xapp_ls[k]));
				       setenv(xapp_env[k], ent->s, 1); }
	if (dfi==2) cfg_write(0);
	backup(QENV('>'), CFG_LOG_BACKUP.i); mklog();
	ac -= i; rv = (const char**) malloc((ac+3)*sizeof(char*));
	if (CFG__1_NOEX.i) { if (ac) memcpy(rv, av+i, ac*sizeof(char*)); rv[ac]=0; return rv; }
	rv[0] = QENV('e'); if (ac) memcpy(rv+1, av+i, ac*sizeof(char*)), rv[ac+1]=0;
			   else    rv[ac+1] = "/"; rv[ac+2] = 0;
	return rv;
}

static int log_pidstat(const char *s, int pid, int stat) {
	int r=0, rf=0; const char * sev = "stopped/resumed/etc";
        if 	(WIFEXITED(  stat)) sev = "exited", ++rf, r = WEXITSTATUS(stat);
        else if (WIFSIGNALED(stat)) sev = "killed", ++rf, r = WTERMSIG(   stat);
	if (s) log("%s %d %s (%d)", s, pid, sev, r);
	return rf ? *sev + 256*hexc1((r>>4)&15) + 65536*hexc1(r&15) : 0;
}

static void pt_dlog(const char * s) {
        int fd = open(QENV('>'), O_WRONLY|O_APPEND);
	switch(fd) {
		case 1:  dup2(1, 2); break;
		case 2:  dup2(2, 1); break;
		default: if (fd<0) return;
			 dup2(fd, 1); dup2(fd, 2); close(fd); break;
	}
        if (s) write(2, s, strlen(s));
}

static void pt_sighup (int _) { pt_sig_flg |= 4; signal(SIGHUP,  &pt_sighup ); }
static void pt_sigint (int _) { pt_sig_flg |= 2; signal(SIGINT,  &pt_sigint ); }
static void pt_sigchld(int _) { pt_sig_flg |= 1; signal(SIGCHLD, &pt_sigchld); }

void pt_sig_act() {
	int pid, r, k, st, sf = __atomic_exchange_n(&pt_sig_flg, 0, __ATOMIC_RELAXED);
	if (sf&6) log("exiting on SIG%s", "INT\0HUP"+(sf&4)), bye(sf&6);
	while( (pid=waitpid(-1, &st, WNOHANG)) > 0 ) {
		BVFOR_JMC(pt_ptab_bv) { if (pt_ptab[j].pid==pid) goto found; }
		log_pidstat("unreg. child", pid, st); continue;
	  found:if (!(r = log_pidstat(pt_ptab[j].s, pid, st))) continue; else pt_ptab_bv &= ~(1u<<j);
		if ((k=pt_ptab[j].s[11])) cmd_report2(k, (const char*)&r);
	}
	signal(SIGCHLD, &pt_sigchld); }

static int pt_io_start(int re) {
	char a2[16], *aa = a2; *(int*)aa = killer_fd>0 ? qh4(killer_fd) : 33; a2[4] = 0;
	*(int*)(a2+8) = qh4(getpid()); a2[12] = 0;   p_close(pt_io_pfd);
	int pf1,pf3,pid = launch(QENV('b'), "><+>", &pf1,pt_io_pfd,QENV('>'),&pf3, "lf.io",a2,a2+8, (char*)0);
	if ((pid|pf1|*pt_io_pfd|pf3)<0) return -1;
	fflush(stdout); fflush(stderr);
	close(1), close(2), dup2(pf3, 1), dup2(pf3, 2), close(pf3);
	set_fd(&pt_cp_m2i, pf1, 0); 
	return pid;
}

int pt_acv_dead(const char*s) {
	int e, k, c;
	if      (*s=='k') 	 e=PTE_ACVSIG,	c=0; 
	else if (!(k=hex2(s+1))) e=0,		c=1;
	else if (k==90) 	 e=PTE_ACVZERO, c=1;
	else			 e=PTE_ACVERR,	c=0;
	if (e) gui_errq_add(e);    if (c&&pt_acv_cw) gui_closewin(ACV_WIN(pt_acv_cur));    return 0;
}

static int con_started(const char *s) {
	if (memcmp(s, "/proc/", 6) || (unsigned int)(s[6]-49)>8u) return GCE_PARSE;
	int fd, pid = atoi(s+6); if (pid&-65536) return GCE_PARSE;
	if ((fd=open(s,O_RDONLY|O_NOCTTY))<0 || dup2(fd,0)<0) return pt_constat = 0, EEE_ERRNO;
	close(fd); pt_constat = pid; *pt_con_pfd = 0; return 0;
}

//////// export ////////////////////////////////////////////

int pt_iw_cmd_sn(int j, const char *s, int n) { int fd = pt_cp_m2iw[j&1]; if (fd<0) return EEE_NOWRK;
				                int r  = write(fd,s,n); return (r==n) ? 0 : EEE_ERRNO+(r>=0); }
int pt_iocmd(char *s) {
	int r, l = strlen(s), c = s[l]; s[l] = 10;
	return r = pt_iocmd_sn(s, l+1), s[l] = c, r; }

int pt_con_op(const char *s) { 
	IFDBGX(PT) log("pt_con_op(%s)", str0(s));
	if (!s) return (pt_constat || !CFG_AUTOCON.i) ? 0 : (pt_constat = -1, pt_iocmd_sn("c\n", 2));
	if (*s!='-') return con_started(s);
	int k; switch(s[1]) {
		case '1': if (pt_constat) { if (pt_constat>0) { k = -3; goto kill; }
				  	    gui_errq_add(EEE_STATE, "con/-1"); }
			  return pt_constat = -1, pt_iocmd_sn("c\n", 2);
		case '2': case '3': if (pt_constat<1) return EEE_STATE; k = 48 - *s; goto kill;
		case '4': if (pt_constat>0) kill(pt_constat, 9), k = pt_constat = 0;
			  else if ((k=-pt_constat-2)&~1) log("con: -4 in state %d",pt_constat), k=pt_constat=0;
			  else pt_constat = -k;
			  if (dup2(pt_nullfd, 0)<0) gui_errq_add(EEE_ERRNO, "con/-4");
			  return *pt_con_pfd = -1, pt_iocmd_sn("Z\nc\n", 2+2*k);
		default:  return EEE_RANGE;
	}
kill:	return kill(pt_constat, 9)<0 ? EEE_ERRNO : (pt_constat = k, 0);
}

int pt_reg_prc(int pid, const char** av, int rprt) {
	int i = 0; BVFOR_JM(~pt_ptab_bv) { i=j; goto found; }
	log("ERROR: process table full"); return pid;
found:  const char * s = av[2*(!strcmp(*av,CFG_XTERM.s) && !memcmp(av[1],"-e",3))];
	int lmax = 12-!!rprt, l = strlen(s)+1; if (l>lmax) s+=l-lmax, l=lmax;
	memcpy(pt_ptab[i].s, s, l); pt_ptab[i].s[11] = rprt; pt_ptab_bv |= (1u<<i);
	IFDBGX(PT) log("reg. prc.: %s %d", pt_ptab[i].s, pid);
	return pt_ptab[i].pid = pid;
}

int pt_acv_op(int id, int opw, const char *a1, const char *a2) {
	int sdl, tdl, pid, op = opw & 255, nwf = opw & 256;
	const char *sdir, *tdir, *src, *trg;
	if ((sdl = CFG_AO_DIR.i )) sdir = CFG_AO_DIR.s;  else sdl = QENVL('t'), sdir = QENV('t');
	if ((tdl = CFG_WAV_DIR.i)) tdir = CFG_WAV_DIR.s; else tdl = QENVL('h'), tdir = QENV('h');
	switch(op) {
		case 0xd: if (unlink(au_file_name(sdir,sdl,id,0,0,"a20"))<0) return EEE_ERRNO; // else FT
		case 0xf: return nwf ? 0 : (gui_closewin(ACV_WIN(id)), 0);
		case 0xe: log("pt_acv_op: 0x0e -- why here?"); return GCE_ROUTE;
		default:  if (op>7) return GCE_PARSE; else break;
	}
	src = au_file_name(sdir, sdl, id, 0,  0,  "a20");
	trg = au_file_name(tdir, tdl, id, a1, a2, "wav\0flac" + (op&4));
	pt_acv_cw = !nwf && !(~op&3);
	switch (op&3) {
		case 0: pid = launch(QENV('b'), "9(kk", "lf.acv",       src, trg,         (char*)0); break;
		case 1: pid = launch(QENV('b'), "9(kk", "lf.acv",       src, trg, a1,     (char*)0); break;
		case 2: pid = launch(QENV('b'), "9(kk", "lf.acv",       src, trg, a1, a2, (char*)0); break;
		case 3: pid = launch(QENV('b'), "9(kk", "lf.acv", "-r", src, trg,         (char*)0); break;
	}
ret:	return (pid<0) ? EEE_ERRNO : (pt_acv_cur=id, 0);
}

int pt_show_lic() { return launch(CFG_XTERM.s, "!())", "-e", QENV('v'), (char*)0); }
void pt_debug()   { return qe_dump(); }
int pt_kill_pa(int flg) { return ((flg&2) ? launch("killall","!(ss","-v","-9","pulseaudio", (char*)0)
					  : launch("pulseaudio","!(ss","-k", (char*)0)) < 0 ? EEE_ERRNO : 0; }

void pt_calc_xbv() {
	for (int i=0; i<N_XAPP; i++) {
		const char **pp = xapp_ls[i];
		for (int j=0; pp[j]; j++) if(pfind(pp[j])) pt_xapp_bv[i] |= (1<<j); }}

static void pt_samp_drop() { if (pt_samp_buf) munmap(pt_samp_buf, 16<<pt_samp_bits), unlink(tpipe_name('+')),
					      pt_samp_buf = 0, pt_samp_bits = 0, pt_wrk_cmd("gx\n", 3); }

#define DEADEND(S) gui_errq_add(PTE_##S##CRASH); errtemp_cond(#S "_dead"); int pid; \
		   if ((pid = pt_##S##_start(1))<0) log("FATAL: %s restart failed", #S), bye(1); return pid
int pt_wrk_dead() { DEADEND(wrk); }
int pt_io_dead()  { pt_dlog("\nERROR: ioprc exited, some logs will be lost...\n");
		    if (pt_constat>0) kill(pt_constat, 9); pt_constat = 0;  DEADEND(io); }

double * pt_samp_shm(int bits) { log("pt_samp_shm: bits=%d", bits);
	if (bits<=pt_samp_bits) return bits ? pt_samp_buf : (pt_samp_drop(), (double*)0); else pt_samp_drop();
	return (pt_samp_buf=(double*)map_wdir_shm('+',16<<bits,3)) ? (pt_samp_bits=bits, pt_samp_buf)
		      : (log("pt_samp_shm: %s", strerror(map_errno)), pt_samp_bits = 0,  pt_samp_buf); }

int pt_wrk_start(int re) {
	p_close(pt_wrk_pfd);
	int pf1, pid = launch(QENV('b'), "><p", &pf1, pt_wrk_pfd, "lf.wrk", (char*)0);
	if ((pid|pf1)<0) return -1;
	if (re) pt_samp_drop(); 
	set_fd(&pt_wrk_m2w, pf1, 0); pt_wrk_pid = pid;
	return pid;
}

int pt_wrk_stop() { return (pt_wrk_pid<=0) ? EEE_STATE :
	(kill(pt_wrk_pid,9)<0) ? EEE_ERRNO : (pt_wrk_pid = 0); }
	
void errtemp_cond(const char *s) {
	if ((pt_errtemp+=CFG_ERRTEMP_INC.i) > CFG_ERRTEMP_INC.i*CFG_ERRTEMP_LIM.i) {
		log ("%s: giving up", s); IFDBGX(AOET) abort(); else  bye(1); }
	else IFDBGX(PT) { log("%s: errtemp=%d", s, pt_errtemp); }}


// e: -1:no_ini(1) -2:inisv.e.(2) -4:optarg_missing 1<<30(+n):argn 1<<29(+n):ini/ln
static void ee_msg() {
	if (ee_flg&1) log("ini file(%s) not found", QENV('j'));
	if (ee_flg&2) gui_errq_add(PTE_INICRE), log("unable to create \"%s\": %s",
						    QENV('j'), ee_errno?strerror(ee_errno):"disk full/bug?");
	if (ee_flg&4) log("last option needs an argument");
	if (ee_flg&((1<<30)+4)) gui_errq_add(PTE_CMDLINE);
	if (ee_flg&((1<<29))  ) gui_errq_add(PTE_INIPARSE);
	for (int k,i=0; i<ee_n; i++) k = ee_v[i]&0xffffff,
		(ee_v[i]&(1<<30)) ? log("unknown option (arg #%d)", k)
				  : log("error in \"%s\", line %d", QENV('j'), k);
}

const char ** pt_init(int ac, char ** av, int *pfd_io, int *pfd_con, int *pfd_wrk) {
	gettimeofday(&tv_zero, 0);
	struct rlimit cl; cl.rlim_cur=cl.rlim_max=1<<28; setrlimit(RLIMIT_CORE, &cl);
	FPU_INI; vstring_set(v_major, v_minor);  signal(SIGPIPE, SIG_IGN);
	pt_con_pfd = pfd_con; pt_io_pfd = pfd_io; pt_wrk_pfd = pfd_wrk;
	qe_ini(); cfg0(); 
	const char ** ret = cfg1(ac, av);
	int r = clk0.ini(17+CFG_CLK_TLBITS.i, CFG_CLK_TYPE.i, 1000000, 100); // TODO: config
	if (r) (r<0) ? (r&=1, log("FATAL: clk0: init (%s, %s)", "clk\0mmap"+4*r, r?map_errno:errno), exit(1))
		     : (log("WARNING: clk0: configured clk type nodes not work"));
	cfg_export(&CFG_CLK_TLBITS); cfg_export(&CFG_QTST_ADIF); cfg_export(&CFG_QTST_RDIF);
	signal(SIGHUP, &pt_sighup); signal(SIGINT, &pt_sigint); 
	if ((pt_nullfd = open("/dev/null", O_RDWR)) < 0 ||
	    (pt_nullfd ? dup2(pt_nullfd, 0) : (pt_nullfd = dup(0))) < 0) perror("FATAL: nullfd"), exit(1);
	signal(SIGCHLD, &pt_sigchld);
	for (const char* s = FIFO_LIST; *s; s++) if (mkfifo(tpipe_name(*s), 0600)<0)
		log("FATAL: failed mkfifo '%c'(%s): %s", *s, tpipe_name(*s), strerror(errno)), exit(1);
	if (pt_io_start(0) < 0) log("FATAL: io_start failed"), exit(1); 
	if (ee_flg) ee_msg();    return ret;
}
