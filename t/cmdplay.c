#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>

int lim_arr[96];
#define LIM(X) lim_arr[(int)(X)-32]


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
			  lim_arr[k] = v; fprintf(stderr, "lim for '%c' set to %d\n", k+32, v);
		  } return;
	default: return he("unknown cmd", s);
}}

void qwe(FILE *f) {
	char buf[1024], *s, *q; buf[1023]=0;
	while (fgets(s=buf+2, 999, stdin)) {
		int l = 0, k = IFHX(*s, 97) ? qh4rs((s+=5)-5) : 0;
		switch(*s) {
			case '#': m_sleep(k); continue;
			case '?': cmd(s+1); continue;
			case '^': goto sl;
			case 'm': case 'N': case 'Z': l = *s; goto lim;
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

int main(int ac, char** av) {
	int i; for (i=0; i<96; i++) lim_arr[i] = INT_MAX;
	if (ac<2) qwe(stdin);
	FILE * f;
	for (i=1; i<ac; i++) if ((f=fopen(av[i],"r"))) qwe(f), fclose(f); else return perror(av[i]), 1;
	return 0;
}
