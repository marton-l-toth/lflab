#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/mman.h>

#define QWE_UTILC_DEF
#include "../s/uc0.h"

static void LOG(const char * fmt, ...) {
        char buf[1024]; va_list ap; va_start(ap, fmt);
        int r = vsnprintf(buf, 1023, fmt, ap); va_end(ap); if (r<0) return;
	buf[r] = 10; write(2, buf, r+1);
}

#define LOG_E(S,X) LOG("%s: error %d, unix err: %s", (S), (X), strerror(errno));

static double * shm_gp_samp;
static int gp_inpipe = -1, gp_outpipe = -1, gp_pid = 0;
static double gp_ctr = 0.5, gp_rad = 0.5;

typedef int (*gp_ls_fun)(char *to, const double *x, const double *y, int n);

static int gp_ls_x0s(char *to, const double *x, const double *y, int n) {
	char *q = to;    int i;     double x0 = x[0], xs = x[1];
	for (i=0; i<n; i++, x0+=xs) q += sprintf(q, "%g %g\n", x0, y[i]);
	return q-to; }

static int gp_ls_xv (char *to, const double *x, const double *y, int n) {
	char *q = to;    int i;
	for (i=0; i<n; i++        ) q += sprintf(q, "%g %g\n", x[i], y[i]);
	return q-to; }

static int gp_start() {
	static const char ini[] = "set style data line\nset autoscale fix\n"
		"set y2tics nomirror format \"%g\"\n";
	static const char *dir[4] = {"Up", "Left", "Down", "Right"},
		     	  *mod[3] = {"", "ctrl-", "alt-"};
        if (gp_pid>0) return -1;
	if ((gp_pid = launch("gnuplot", "!><1", &gp_inpipe, &gp_outpipe, (char*)0)) < 0) return -2;
	if ((gp_inpipe|gp_outpipe)<0) return -3;
//	"bind Up     'print \"%8K^q\"'\n";
	char buf[4096], *q = buf + sizeof(ini)-1; memcpy(buf, ini, sizeof(ini)-1);
	int i; for (i=0; i<12; i++) 
		q += sprintf(q, "bind '%s%s' 'print \"%8K^q%c\"'\n",mod[i>>2],dir[i&3],40+i);
	if (write(gp_inpipe, buf, q-buf) != q-buf) return -4;
        return 0;
}

static unsigned int * gp_tlog;
static int gp_tlog_len, gp_tlog_siz = 0;
static int gp_resol = 512;

static int gp_tlog_read(const char *fname) {
	int fd = open(fname, O_RDONLY); if (fd<0) return -1;
	int n = lseek(fd, 0, SEEK_END); if (n<4) return -2; else n>>=2;
	int r = lseek(fd, 0, SEEK_SET); if (r) return -3;
	if ((gp_tlog_len=n)>gp_tlog_siz) {
		if (n>(1<<28)) return -4;
		int sz = 2*gp_tlog_siz; if (sz) free(gp_tlog); else sz = 4096;
		while (sz<n) sz+=sz;
		gp_tlog = malloc(4*(gp_tlog_siz=sz));
	}
	return ((r = read(fd, gp_tlog, 4*n)) == 4*n) ? (close(fd), 0) : -5;
}

static int gp_tlog_plot() {
	static const char hd[] = "set term x11\nplot '-', '-' axes x1y2 \n0 0\n";
	char buf[65536], *q=buf+sizeof(hd)-1; memcpy(buf, hd, sizeof(hd)-1);
 	int x, x0, k0, j, len = gp_tlog_len, samp_cnt = 0;
	double t = 0.0, ttot = 0.0, td = 0.0, lennn = (double)len;
	int j0 = ivlim((int)lround(lennn*(gp_ctr-gp_rad)), 0,  len-1),
	    jz = ivlim((int)lround(lennn*(gp_ctr+gp_rad)), j0+1, len);
	for (j=j0; j<jz; j++) { ttot += (double)(gp_tlog[j]&262143); }
	double tdmin = ttot / (double) gp_resol, t1 = tdmin;
	double samp_calc = 0.0;

	for (j=j0; j<jz; j++) {
		int k = (x=gp_tlog[j]) >> 18, xt = x&262143;
		if (k<4096 && k>127) { td += (double)(k-1090); }
		if ((t += (double)(xt)) > t1) 
			q+=sprintf(q,"%.15g %.15g\n",1e-6*t,1e-3*td), t1=t+tdmin;
	}
	memcpy(q, "e\n0 0\n", 6); q += 6; t = 0.0; t1 = tdmin;
	for (j=j0, x0=0; j<jz; j++,x0=x) {
		int k0, k = (x=gp_tlog[j]) >> 18, xt = x&262143;
		if ((k|1)=='q' && (k0=(x0>>18)-4096)>=0)
			samp_cnt += k0, samp_calc += (double)(xt+(x0&262143)); 
		if (((t += (double)(xt)) > t1) && samp_cnt)
			q+=sprintf(q,"%.15g %.15g\n",1e-6*t,.0441*samp_calc/(double)samp_cnt),
				samp_cnt = 0, t1=t+tdmin, samp_calc = 0.0;
	}
	q[0] = 'e'; q[1] = 10; q += 2; LOG("writing %d", q-buf);
	//write(2, buf, q-buf);
	return write(gp_inpipe, buf, q-buf) == (q-buf) ? 0 : -1;
}

static void gp_key(int c) {
	static const signed char xstp[8] = {-5,5,-30,30,-1,1,0,0};
	static const signed char zstp[8] = {-4,4,-16,16,-1,1,0,0};
	if (c>51) { return LOG("unknown key 0x%x",c); }
	if ((c-=40)<0) { return LOG("invalid key 0x%x",c+40); }
	double g_c = gp_ctr, g_r = gp_rad;
	if (c&1) g_c += .2*(double)(int)(xstp[(c>>1)&7])*g_r;
	else     g_r *= exp(M_LN2*.125*(double)(int)(zstp[(c>>1)&7]));
	if (g_r>.5) { g_r = g_c = .5; goto crok; }
	if (g_r<4.76837158203125e-07) g_r = 4.76837158203125e-07;
	if (g_c<g_r) g_c = g_r; else if (g_c>1.0-g_r) g_c = 1.0-g_r;
crok:	if (fabs(g_c-gp_ctr)>1e-11 || fabs(g_r-gp_rad)>1e-11) {
		LOG("gp_key: c=%g r=%g", gp_ctr = g_c, gp_rad = g_r);
		if (gp_tlog_len) gp_tlog_plot(); }
	else { LOG("gp_key: nothing happens"); }
}

static void gp_out() {
	static int state = 37;
	signed char buf[4096];
	int i, c, r = read(gp_outpipe, buf, 4096);
	if (r<=0) perror("gp_pipe/read"), exit(1);
//	LOG("gp_read: %d \"%s\"", r, r<999 ? (const char*)(buf[r]=0,buf) : "...");
	for (i=0; i<r; i++) {
		if ((c=buf[i])==state) { state += 19; }
		else { if (state==132) gp_key(c); state = 37; }}}

static int gp_cmd(const char *s) {
	int r; switch(*s) { 
		case 't' : return ((r=gp_tlog_read(s+1))<0) ? r : gp_tlog_plot();
		default: return -9;
	}}

static void bye(int j) { LOG("bye %d", gp_pid); if (gp_pid>0) kill(gp_pid,9); exit(j); }

int main() {
	char s[256];
	const char * fnm = tpipe_name('%');
	int r = gp_start(); if (r<0) LOG_E("gp_start", r);
	fd_set rset; struct timeval tv;
	while (1) {
		FD_ZERO(&rset); FD_SET(0, &rset); FD_SET(gp_outpipe, &rset);
		r = select(gp_outpipe+1, &rset, 0, 0, NULL);
		if (FD_ISSET(gp_outpipe, &rset)) gp_out();
		if (FD_ISSET(0, &rset)) {
			int l = read(0, s, 256); if (l<=0) perror("stdin"), bye(!!l);
			if (s[l-1]==10) s[--l] = 0;
			switch(*s) {
				case 'g': if ((r=gp_cmd(s+1))<0) LOG_E("gp_cmd",r); break;
				case 'q': bye(0); break;
				default: puts("error"); break;
		}}}}
