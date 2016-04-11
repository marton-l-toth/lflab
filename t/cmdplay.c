#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>

int lim_arr[96];
#define LIM(X) lim_arr[(int)(X)-32]

static inline int hxd2i(int c) { return 15&(c+9*(c>>6)); }

int qh4r(unsigned int x) {
        unsigned int ten = x & 0x40404040;
        x = (x&0xf0f0f0f) + 9*(ten>>6);
        x = (x&0xf000f00) | ((x&0xf000f)<<12);
        return (int)((x&0xff00) | ((x>>24)&255));
}

int m_sleep(int x) { if (!x) return;
        struct timespec ts; ts.tv_sec = x/1000; ts.tv_nsec = 1000000*(x%1000); int ec;
        do ec = nanosleep(&ts, &ts); while (ec==EINTR); return ec; }

#define qh4rs(P) qh4r(*(const int*)(P))
#define IFHX(X,T) ( (unsigned int)((X)-48)<10u || (unsigned int)((X)-(T))<6u )

void he(const char *s1, const char *s2) { fprintf(stderr, "%s: \"%s\"\n", s1, s2); }

void cmd(const char *s) { switch(*s) {
	case 'L': for (++s; *s; ) {
			  while (*s==32) ++s;  if (*s==10) return;
			  int v=0, d, k = *s-32; if ((unsigned int)k>95u) return he("invalid ch",s);
			  if (*s=='*') { ++s; v = INT_MAX; }
			  else         { for (++s; (unsigned int)(d=*s-48)<10u; s++) v = 10*v + d; }
			  lim_arr[k] = v; if (!lim_arr[31]) fprintf(stderr, "lim '%c' set to %d\n", k+32, v);
		  } return;
	default: return he("unknown cmd", s);
}}

void qwe(FILE *f) {
	char buf[1024], *s, *q; buf[1023]=0;
	while (fgets(s=buf+2, 999, f)) {
		int l = 0, k = IFHX(*s, 97) ? qh4rs((s+=5)-5) : 0;
		switch(*s) {
			case '#': m_sleep(k); continue;
			case '?': cmd(s+1); continue;
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

static void qwe_tij(int t, int i, int j) {
	FILE * f = at_f_ij(i,j);
	f ? (m_sleep(t), qwe(f), fclose(f)) : perror(at_fnm);
}

static void at_dir(const char * dir) {
	int i, j, k, l = strlen(dir);   FILE *f;
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
	int i,j,k,st[32]; for (i=0; i<at_n; i++) st[i] = 1;
	while(1) {
		i = random() % at_n;
		if (!st[i]) qwe_tij(200, i, 0), st[i] = 1;
		else if (st[i]<=at_k[i]) qwe_tij(200, i, st[i]), st[i]++;
		else if ((k=random()%st[i])<at_k[i]) qwe_tij(200, i, k), st[i]++;
		else qwe_tij(200, i, -1), st[i] = 0;
	}}

static void at_main(const char * opt) {
	int i,j,k, flg = (*opt=='-') ? (opt+=2, hxd2i(opt[-1])) : 0, len = strlen(opt);
	at_dir(opt);
	if (flg&1) for (i=0; i<at_n; i++) qwe_tij(200, i, 0);
	if (flg&2) for (i=0; i<at_n; i++) for (j=1,k=at_k[i]; j<k; j++) qwe_tij(200, i, j);
	if (flg&4) for (i=0; i<at_n; i++) at_random();
	if (flg&8) for (i=0; i<at_n; i++) qwe_tij(200, i, -1);
}

int main(int ac, char** av) {
	int i; for (i=0; i<96; i++) lim_arr[i] = INT_MAX;
	if (ac<2) qwe(stdin);
	FILE * f;
	for (i=1; i<ac; i++) {
		if (*av[i] == '-') { switch(av[i][1]) {
			case 'a': at_main(av[i]+2); continue;
			default: fprintf(stderr, "unknown option \"%s\"\n", av[i]); 
		}}
		if ((f=fopen(av[i],"r"))) qwe(f), fclose(f); else return perror(av[i]), 1;
	}
	return 0;
}
