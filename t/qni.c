#ifndef EXPR
#error this file has to be compiled with EXPR defined, see source
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

#define FP0(S) 		fprintf(stderr, S "\n")
#define FP1(S,X) 	fprintf(stderr, S "\n", (X))
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y))
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z))
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W))

static struct random_data r_dat; static char r_buf[256];


int main(int ac, char** av) {
	initstate_r(time(NULL), r_buf, 256, &r_dat);
	struct timeval tv1, tv2;
	int i, j, k, n = 1000000; double x0=0.0, x1=1.0;
	switch(ac) {
		case 4: x1 = atof(av[3]);
		case 3: x0 = atof(av[2]);
		case 2: n *= atoi(av[1]); break;
		default: FP1("usage: %s n x0 x1", *av); return 1;
	}
	gettimeofday(&tv1, NULL);
	double s, x, y, z, w, xd = (x1-x0) / (double)(n);
	FP3("x0=%g x1=%g xd=%g", x0, x1, xd);
	x = x1; s = (EXPR); x = x0; s += (EXPR); s *= .5;
	for (i=1; i<n; i++) x+=xd, s += (EXPR);
	gettimeofday(&tv2, NULL);
	FP2("t=%g sum=%.15g", 1e2*(tv2.tv_sec-tv1.tv_sec)+1e-4*(tv2.tv_usec-tv1.tv_usec), s*xd);
	return 0;
}
