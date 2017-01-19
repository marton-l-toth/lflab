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

int main(int ac, char** av) {
	struct timeval tv1, tv2;
	int i, n = 10000000; double x0=0.0, x1=1.0;
	switch(ac) {
		case 3: x1 = atof(av[2]);
		case 2: x0 = atof(av[1]);
		case 1: break;
		default: FP1("usage: %s n x0 x1", *av);
	}
	gettimeofday(&tv1, NULL);
	double s, x, y, z, w, xd = (x1-x0) / (double)(n+1);
	x = x1; s = (EXPR); x = x0; s += (EXPR); s *= .5;
	for (i=1; i<n; i++) x+=xd, s += (EXPR);
	gettimeofday(&tv2, NULL);
	FP2("t=%g sum=%.15g", 1e2*(tv2.tv_sec-tv1.tv_sec)+1e-4*(tv2.tv_usec-tv1.tv_usec), s/(double)n);
	return 0;
}
