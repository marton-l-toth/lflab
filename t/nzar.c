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

#include "../s/distr.h"

int ar_test(dfun1_t f, double a, int n) {  double x = 0.0; int i = 0;
					   for(;;){x += a/(*f)(x); ++i;
						   if (x>=1.0) return i<n || x>1.0;
						   if (i >= n) return -1;
					   }}
double ar_find(dfun1_t f, int n) {
	double lo0 = -1.0, hi, lo = .1/(double)n;
	int r = ar_test(f, lo, n); if (!r) return lo;
	if (r>0) { hi = lo; lo = lo0 = 0.0; }
	else 	 { while(1) { if (!(r=ar_test(f,hi=lo+lo,n))) return hi; if (r>0) break; lo=lo0=hi; }}
	while(1) {
		double mid = .5*(lo+hi); if (mid>=hi || mid<=lo) return lo0;
		if (!(r=ar_test(f, mid, n))) return mid;
		if (r<0) lo=lo0 = mid; else hi=mid;
	}}

void stat(dfun1_t f, double *q, double a, int n) {
	double v0, x = 0.0, ia = 0.0, v1 = f(0.0); int i;
	for (i=0; i<n; i++) v0=v1, v1=f(x+=a/(*f)(x)), ia += v1/v0;
	q[0] = x; q[1] = ia / (double)n;
}

void qw(int c, const char * s) {
	DISTR_NF(c); double a = ar_find(f,n), sr[2];  stat(f, sr, a, n);
	printf("{%.16g, '%c', %s}, // x1:%.16g iacc:%g\n", a, c, s, sr[0], sr[1]);
}

int main(int ac, char** av) {
	int i; for (i=1; i<ac; i++) qw(av[i][0], av[i]+1);  return 0; }
