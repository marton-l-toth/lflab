#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>

int qh4r(unsigned int x) {
        unsigned int ten = x & 0x40404040;
        x = (x&0xf0f0f0f) + 9*(ten>>6);
        x = (x&0xf000f00) | ((x&0xf000f)<<12);
        return (int)((x&0xff00) | ((x>>24)&255));
}

int u_sleep(int x) {
	long long ns = 1000ll * (long long)x;
        struct timespec ts; ts.tv_sec = ns/1000000000ll; ts.tv_nsec = ns%1000000000ll; int ec;
        do ec = nanosleep(&ts, &ts); while (ec==EINTR); return ec; }

#define qh4rs(P) qh4r(*(const int*)(P))

void qwe(FILE *f) {
	char buf[1024];
	while (fgets(buf, 1023, stdin)) {
		int k = qh4rs(buf), l = strlen(buf+5); if (k) u_sleep(1000*k);
		buf[3] = 'Q'; buf[4] = 'P'; if (buf[l+4]!=10) buf[++l+4] = 10;
		write(1, buf+3, l+2);
	}}

int main(int ac, char** av) {
	if (ac<2) qwe(stdin);
	FILE * f; int i;
	for (i=1; i<ac; i++) if ((f=fopen(av[i],"r"))) qwe(f), fclose(f); else return perror(av[i]), 1;
	return 0;
}
