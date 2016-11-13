#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdarg.h>
#define QWE_UTILC_DEF
#include "../s/uc0.h"

int main(int ac, char** av) {
	char s[256]; double x; int qw=0;
	while (fgets(s, 255, stdin)) {
		if (*s!='^' || s[1]!='<') { if (!qw) qw=1, puts("---------"); continue; }
		int i; for (qw=0,i=2; i<78; i+=11) dcd11(&x, s+i), printf(" %.12g",x); puts(""); }
	return 0;
}
