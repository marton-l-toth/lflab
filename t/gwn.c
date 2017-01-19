#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define FP0(S) 		fprintf(stderr, S "\n")
#define FP1(S,X) 	fprintf(stderr, S "\n", (X))
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y))
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z))
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W))


static inline double f(double x) { return exp(-16.0*x*x); }



double sum(double a, int n) {
	double x = 0.0; int i;
	printf("{0.0", a);
	for (i=0; i<n; i++) x+=a/f(x), printf(",%.15g",x); puts("}");
	return x;
}

int main(int ac, char** av) {
	int i;
	if (ac!=2) return FP1("usage: %s n", *av), 1;
	int n = atoi(av[1]);
	double r, lo=.1/(double)n, hi = 0.0;
	if ((r=sum(lo,n))==1.0) return 0;
	if (r>1.0) { hi = lo; lo = 0.0; }
	else { while(1) { if ((r=sum(hi=lo+lo,n))==1.0) return 0; if (r>1.0) break; lo=hi; }}
	while(1) {
		double mid = .5*(lo+hi); if (mid>=hi || mid<=lo) return 0;
		if ((r=sum(mid,n))==1.0) return 0;
		if (r<1.0) lo=mid; else hi=mid;
	}
	return 0;
}
