#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define FP0(S) 		fprintf(stderr, S "\n");
#define FP1(S,X) 	fprintf(stderr, S "\n", (X));
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y));
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z));
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W));

static void enc11(char* to, double x) {
        int i,j; unsigned int qw[2]; memcpy(qw, &x, 8);
        to[0] = 42 + (qw[0]>>30) + ((qw[1]>>28)&12);
        for (i=0;i<2;i++) for (j=0; j<5; j++) to[5*i+j+1] = 59 + ((qw[i]>>(6*j))&63);
}

static int dcd11(double *to, const char * s) {
        unsigned int qw[2], k = s[0] - 42u; if (k>15u) return 0;
        qw[0] = (k&3u)<<30; qw[1] = (k&12u)<<28;
        int i,j; for (i=0;i<2;i++) { for (j=0; j<5; j++) {
                if ((k=s[5*i+j+1]-59) > 63u) return 0; else qw[i] |= k << (6*j); }}
        return memcpy(to, qw, 8), 1;
}



int main(int ac, char** av) {
	char buf[1024]; double x;
	while (fgets(buf, 1023, stdin)) {
		int i,l = strlen(buf); if (buf[l-1]==10) buf[--l] = 0;
		if (buf[0]=='^'&&buf[1]=='<') for (i=2; dcd11(&x, buf+i); i+=11) enc11(buf+i, x*(.5*M_SQRT2));
		puts(buf); }
	return 0;
}
